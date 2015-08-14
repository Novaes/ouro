require("dotproduct")
local haddavx = terralib.intrinsic("llvm.x86.avx.hadd.ps.256", { vector(float,8), vector(float,8) } -> vector(float,8))
local cstdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")

-- set number to float to Single Float Point computation
local number = float
local alignment = 8

local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})

local dotune = true

-- terra naivel1conv(A : &number, B : &number, C : &number, lda : int, ldb : int, ldc : int, alpha : number)
function symmat(name,I,...)
	if not I then return symbol(name) end
	local r = {}
	for i = 0,I-1 do
		r[i] = symmat(name..tostring(i),...)
	end
	return r
end

function printMatrix(m,rows,columns)
  local matrix = m
  for i=0,rows-1 do
    for j=0,columns-1 do
      io.write(" " .. matrix[i*columns + j])
    end
    io.write("\n")
  end
  io.write("\n")
end

local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end

local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)

terra hadd(v : vector(double,8))
	var v1 = haddavx(v,v)
	var v2 = haddavx(v1,v1)
	return v2[0] + v2[4]
end

-- generate L1 convolution 
function genkernel(NB, RM, RN, V, prefetch, K, L, boundary)
	local M,N, boundaryargs
	-- if one of the parameters is lower than NB then receive the usual
	if boundary then 
		M,N = symbol(int64,"M"),symbol(int64,"N")
		boundaryargs = terralib.newlist({M,N})
	else
		boundaryargs = terralib.newlist()
		M,N = NB,NB
	end

	local cx,cy = math.floor(K/2), math.floor(L/2)
	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2*cx,RN+2*cy), symmat("b",K,L/V), symmat("c",RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y,restHadd = symbol("x"), symbol("y"), symbol("restHadd")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()


	local VP = &vector(number,V)
	local VimgP = &vector(number,V) -- scale it

	for m = 0, RM+2*cx-1 do
		for n = 0, RN+2*cy-1 do
			loadA:insert(quote
					var [a[m][n]] : vector(number,V) = unalignedload(VimgP(&A[m*lda + n*V]))
			end)
		end
	end

	for m = 0, RM-1 do
		for n = 0, RN-1 do
			loadc:insert(quote
					var [c[m][n]] = alpha * C[m*ldc + n]
			end)
			storec:insert(quote
				C[m*ldc + n] = [c[m][n]]
			end)
		end
	end

	local calcc = terralib.newlist()

	-- load full kernel
	for  k=0, K-1 do
		for l = 0, L-V do
			loadkernel:insert(quote
				var [b[k][l]] : vector(number,V) = unalignedload(VP(&B[k*ldb + l*V]))
			end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			for k=0, K-1 do
				for l = 0, L-V do
					x, y = m + (k - math.floor(K/2) ), n + (l - math.floor(L/2))
					if V > 8 then
						assert(V <= 15)
						calcc:insert(quote
							var v : vector(number,V) = @[&vector(number, V)](&[a[x+cx][y+cy]]) * [b[k][l]]
							var sum = 0
							for i=8,11 do
								sum = sum + v[i]
							end
							[c[m][n]] = [c[m][n]] + hadd(@[&vector(number,8)](&v)) + sum
						end)
					else
						calcc:insert(
							quote
								-- area regblocking not multiple of the area sizeblocking
								-- dot product
									var v : vector(number,V) = @[&vector(number, V)](&[a[x+cx][y+cy]]) * [b[k][l]]
									[c[m][n]] = [c[m][n]] + hadd(@[&vector(number,8)](&v))
								-- v = [a[x+cx][y+cy]]
								-- var sum = 0
								-- for i=0,V do
								-- 	sum = sum + v[i]
								-- end
								
								-- [c[m][n]] = [c[m][n]] + dotprod([&float](&[a[x+cx][y+cy]]), [&float](&[b[k][l]]), V)
								
						end)
					end
				end
			end
		end
	end

	return terra([A] : &number, [B] : &number, [C] : &number, [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : number, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		[loadkernel];
		for [mm] = 0, M, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=0, N, RN do 
				[loadc];
				-- llvmprefetch(A + sda*lda,0,3,1);
				[loadA];
				[calcc];
				[storec];
				A = A + RN
				C = C + RN
			end
			-- if (((N-2)/(RN)) * (RN) + 1 < N-1) then
			-- 	var offset = (((N-2)/(RN)) * (RN) + 1) + (RN)  - (N-1)
			-- 	A = A - offset
			-- 	C = C - offset
			-- end
			A = A + RM * lda - M
			C = C + RM * ldc - M
		end
	end
end

function genkernelb(NB, RM, RN, V, prefetch, K, L, boundary)
	local M,N, boundaryargs
	-- if one of the parameters is lower than NB then receive the usual
	if boundary then 
		M,N = symbol(int64,"M"),symbol(int64,"N")
		boundaryargs = terralib.newlist({M,N})
	else
		boundaryargs = terralib.newlist()
		M,N = NB,NB
	end

	local cx,cy = math.floor(K/2), math.floor(L/2)
	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2*cx,RN+2*cy), symmat("b",K,L), symmat("c",RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y = symbol("x"), symbol("y")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM+2*cx-1 do
		for n = 0, RN+2*cy-1 do
			loadA:insert(quote
					var [a[m][n]] = A[m*lda + n]
			end)
		end
	end

	for m = 0, RM-1 do
		for n = 0, RN-1 do
			loadc:insert(quote
					var [c[m][n]] = alpha * C[m*ldc + n]
			end)
			storec:insert(quote
				C[m*ldc + n] = [c[m][n]]
			end)
		end
	end

	local calcc = terralib.newlist()

	-- load full kernel
	for  k=0, K-1 do
		for l = 0, L-1 do
			loadkernel:insert(quote
				var [b[k][l]] = B[k*ldb + l]
			end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			for k=0, K-1 do
				for l = 0, L-1 do
					-- would sum mm or nn, but this position is realtive to this mini-block (rm, rn)
					x, y = m + (k - math.floor(K/2) ), n + (l - math.floor(L/2))
					--no boundary cases
					calcc:insert(
						quote
							-- area regblocking not multiple of the area sizeblocking
							[c[m][n]] = [c[m][n]] + [a[x+cx][y+cy]] * [b[k][l]]
						end
					)
				end
			end
		end
	end

	return terra([A] : &number, [B] : &number, [C] : &number, [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : number, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		[loadkernel];
		for [mm] = 0, M, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=0, N, RN do 
				[loadc];
				-- llvmprefetch(A + sda*lda,0,3,1);
				[loadA];
				[calcc];
				[storec];
				A = A + RN
				C = C + RN
			end
			-- if (((N-2)/(RN)) * (RN) + 1 < N-1) then
			-- 	var offset = (((N-2)/(RN)) * (RN) + 1) + (RN)  - (N-1)
			-- 	A = A - offset
			-- 	C = C - offset
			-- end
			A = A + RM * lda - M
			C = C + RM * ldc - M
		end
	end
end

terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

function blockedloop(M,N,blocksizes,bodyfn)
  local function generatelevel(n,ii,jj,bb0,bb1)
    if n > #blocksizes then
      return bodyfn(ii,jj)
    end
    local blocksize = blocksizes[n]
    return quote for i = ii,min(ii+bb0,M),blocksize do
                   for j = jj,min(jj+bb1,N),blocksize do
                        [ generatelevel(n+1,i,j,blocksize,blocksize) ]
           		end end end
  end
  return generatelevel(1,0,0,M,N)
end

function naive()

	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		ldc : int, kCenterX: int, kCenterY: int) 
		var e = 0
		for i=0, M do
		    for j=0, N do
		    	C[e] = 0
	        	for m=0, K do
	          		for n=0,L do
			            -- boundaries
			            var ii: int = i + (m - kCenterY)
			            var jj: int = j + (n - kCenterX)
			            if ii>=0 and ii<M and jj>=0 and jj<N then
			              C[e] = C[e] + A[ii * N + jj] * B[m * L + n]
			            end
		        	end
		      	end
		      	e = e+1
		    end
		end
	end
end

function gennaive()

	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		ldc : int, kCenterX: int, kCenterY: int) 
		var e = 0
		for i=kCenterX, M-kCenterX do
		    for j=kCenterY, N-kCenterY do
		    	C[e] = 0
	        	for m=0, K do
	          		for n=0,L do
			            -- boundaries
			            var ii: int = i + (m - kCenterY)
			            var jj: int = j + (n - kCenterX)
			            -- if ii>=0 and ii<M and jj>=0 and jj<N then
			              C[e] = C[e] + A[ii * N + jj] * B[m * L + n]
			            -- end
		        	end
		      	end
		      	e = e+1
		    end
		end
	end
end

function genconvolution(NB,NBF,RM,RN,V,K,L)
	-- need this because RM and RN are unrolled fully inside NBxNB block
	if not isinteger(NB/RN) or not isinteger(NB/RM) then
		print("3rd level blocksizes must be a multiple of the 2nd")
		return false
	end

	--5 times NB minimum by dgemm
	--local NB2 = NBF * NB
	local NB2 = NB * NBF
	local l1conv0 = genkernel(NB, RM, RN, V, false, K, L, false)
	local l1conv0b = genkernelb(NB, 1, 1, V, false, K, L, true)

	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		ldc : int, kCenterX: int, kCenterY: int) 

		M = M-kCenterX
		N = N-kCenterY
		for mm = kCenterX,M,NB2 do
			for nn = kCenterY,N,NB2 do
				for m = mm,min(mm+NB2,M),NB do
					for n = nn,min(nn+NB2,N),NB do
					 	var MM,NN = min(M-m,NB),min(N-n,NB)
	                    var isboundary = MM < NB or NN < NB
	                    var AA,CC = A + ((m-kCenterX)*lda + (n-kCenterY)),C + ((m-kCenterX)*ldc + (n-kCenterY))
	                    if isboundary then -- do not enter here YET
	                    	l1conv0b(AA,
	                         B,
	                         CC,
	                         sda,lda,ldb,ldc,0,MM,NN)
	                    else
	                        l1conv0(AA,
	                         B,
	                         CC,
	                         sda,lda,ldb,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
	                    end
					end
				end
			end
		end

         -- [ blockedloop(M,N,{NB2,NB},
         --        function(m,n) 
         --        return quote
         --            var MM,NN = min(M-m,NB),min(N-n,NB)
         --            var isboundary = MM < NB or NN < NB
         --            var AA,CC = A + (m*lda + n),C + (m*ldc + n)
         --            if isboundary then -- do not enter here YET
         --             l1conv0b(AA,
         --                 B,
         --                 CC,
         --                 sda,lda,ldb,ldc,0,MM,NN)
         --            else
         --                l1conv0(AA,
         --                 B,
         --                 CC,
         --                 sda,lda,ldb,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
         --            end
         --        end end)  
         --    ]       
	end
end

-- Different blocksizes for the same result implies in padding overheading 
-- for small blocks
local blocksizes = {8,16,24,32,40,48,56,64,80,88}
local regblocksM = {1,2,4}
local regblocksN = {1}
local vecfilters = {3}

-- initialized (defined structure of best)
local best = { gflops = 0, nb = 5, b = 24, rm = 2, rn = 2, v = 3, k = 3 }
local NB2 = {2}

if dotune then
	local tunefor = 1024 -- full size of the matrix
	--change for 10 later
	local harness = require("direct-vec-matrixtestharness")
	for _,nb in ipairs(NB2) do
		for _,b in ipairs(blocksizes) do
			for _,rm in ipairs(regblocksM) do
				for _,rn in ipairs(regblocksN) do
					for _,k in ipairs(vecfilters) do
						-- same until here
						v = k
						local my_conv = genconvolution(b,nb,rm,rn,k,k,k)
						-- local my_conv = gennaive()
						if my_conv then
							print(b,rm,rn,k)
							my_conv:compile()
							-- bellow line makes do not need boundary cases (image multiple of blocksize)
							local i = math.floor(tunefor / b) * b
							local curr_gflops = 0
							local ctyp
							local correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k, function(Me,Ne,K,L,M,N,A,B,C)
		                    	my_conv(nil,Me,Ne,K,L,1.0,A,Me,Ne,B,L,C,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
							end)
							-- if not correct then	print("<error>")  break  end
							-- print(i,unpack (exectimes),"[OK]")
							local curr_gflops = exectimes[1]
							print(curr_gflops) -- print analysis 
							if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
								best = { gflops = curr_gflops, nb = nb, b = b, rm = rm, rn = rn, v = v, k = k }
								-- terralib.tree.printraw(best)
							end
						end
					end
				end
			end
		end
	end
end

-- Naive receiving an normal sized image
local my_naivenumconv = gennaive()
terralib.saveobj("../bin/my_naivenumconv.o", {my_naivenumconv = my_naivenumconv})
-- the vectorsize is the same of vecfiltersize
local my_numconv = genconvolution(best.b,best.nb,best.rm,best.rn,best.k,best.k,best.k)
terralib.saveobj("../bin/my_numconv.o", {my_numconv = my_numconv})
