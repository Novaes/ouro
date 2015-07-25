local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
local number = double


local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})

local dotune = true

-- terra naivel1conv(A : &double, B : &double, C : &double, lda : int, ldb : int, ldc : int, alpha : double)
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

-- generate L1 convolution 
function genkernel(NB, RM, RN, V, prefetch, K, L, depth, boundary)
	local M,N, boundaryargs
	-- if one of the parameters is lower than NB then receive the usual
	if boundary then 
		M,N = symbol(int64,"M"),symbol(int64,"N")
		boundaryargs = terralib.newlist({M,N})
	else
		boundaryargs = terralib.newlist()
		M,N = NB,NB
	end

	-- assets NB/RM and NB/RN do not necessary
	-- print("parameters: "..NB .." ".. RM .." ".. RN .." ".. V .." ".. K .." ".. L)

	local VP = &vector(double,V)
	local terra vecload(data : &double, idx : int)
		var addr = &data[idx]
		return @VP(addr)
	end

	local terra vecstore(data : &double, idx : int, v : vector(double,V))
		var addr = &data[idx]
		@VP(addr) = v
	end

	local cx,cy = math.floor(K/2), math.floor(L/2)
	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2*cx,RN+2*cy), symmat("b",depth,K,L), symmat("c",depth,RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y = symbol("x"), symbol("y")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM+2*cx-1 do
		for n = 0, RN+2*cy-1 do
			loadA:insert(quote
					var [a[m][n]] = vecload(A, m*ldc + n*V)
			end)
		end
	end

	for  i=0, depth-1 do
		for m=0,RM -1 do
			for n=0,RN -1 do
				loadc:insert(quote
						-- sda*lda*i is the base when processing multiple kernel results
						var [c[i][m][n]] = alpha * vecload(C, sda*lda*i  + m*ldc + n*V)
				end)
				storec:insert(quote
					vecstore(C, sda*lda*i + m*ldc + n, [c[i][m][n]])
				end)
			end
		end
	end
	local calcc = terralib.newlist()

	-- load full kernel
	for  i=0, depth-1 do
		for  k=0, K-1 do
			for l = 0, L-1 do
				loadkernel:insert(quote
					var [b[i][k][l]] = vecload(B, ldb*ldb*i  + k*ldb + l*V)
				end)
			end
		end
	end

	-- spatial 2D convolution
	for i=0,depth-1 do -- multiple kernel
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
							[c[i][m][n]] = [c[i][m][n]] + [a[x+cx][y+cy]] * [b[i][k][l]]
							end
						)
					end
				end
			end
		end
	end

	return terra([A] : &double, [B] : &double, [C] : &double, [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		[loadkernel];
		for [mm] = 0, M, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN*V,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=0, N, RN*V do 
				[loadc];
				llvmprefetch(A + sda*lda,0,3,1);
				[loadA];
				[calcc];
				[storec];
				A = A + RN*V
				C = C + RN*V
			end
			-- if (((N-2)/(RN*V)) * (RN*V) + 1 < N-1) then
			-- 	var offset = (((N-2)/(RN*V)) * (RN*V) + 1) + (RN*V)  - (N-1)
			-- 	A = A - offset
			-- 	C = C - offset
			-- end
			A = A + RM * ldc - M
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

function gennaiveconv()
	return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
		alpha : double, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
		ldc : int, kCenterX: int, kCenterY: int) 
		var kCenterX : int = K/2
		var kCenterY : int = L/2
		for i=kCenterX, M-kCenterX do
		    for j=kCenterY, N-kCenterY do
	        	for m=0, K do
	          		for n=0,L do
			            -- boundaries
			            var ii: int = i + (m - kCenterY)
			            var jj: int = j + (n - kCenterX)
			            if ii>=0 and ii<M and jj>=0 and jj<N then
			              C[i*N + j] = C[i*N + j] + A[ii * N + jj] * B[m * L + n]
			            end
		        	end
		      	end
		    end
		end
	end
end

function genconvolution(NB,NBF,RM,RN,V,K,L,NF)
	-- register blocking does not need to be a a multiple of the blocksize anymore
	if not isinteger(NB/(RN*V)) or not isinteger(NB/RM) then
		print("3rd level blocksizes must be a multiple of the 2nd")
		return false
	end

	--5 times NB minimum by dgemm
	--local NB2 = NBF * NB
	local NB2 = NB * NBF
	local l1conv0 = genkernel(NB, RM, RN, 1, false, K, L, NF, false)
	local l1conv0b = genkernel(NB, 1, 1, 1, false, K, L, NF, true)

	return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
		alpha : double, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
		ldc : int, kCenterX: int, kCenterY: int) 

		M = M-kCenterX
		N = N-kCenterY

		for mm = kCenterX,M,NB2 do
			for nn = kCenterY,N,NB2 do
				for m = mm,min(mm+NB2,M),NB do
					for n = nn,min(nn+NB2,N),NB do
					 	var MM,NN = min(M-m,NB),min(N-n,NB)
	                    var isboundary = MM < NB or NN < NB
	                    var AA,CC = A + ((m-kCenterX)*lda + (n-kCenterY)),C + (m*ldc + n)
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
-- ending in s means SIZE
-- starting with n, means NUMBER
local blocksizes = {8}--16,24,32,40,48,56,64}
local regblocks = {1}--,2,4,8}
-- local vectors = {1,2,4,8,16}
local vectors = {1} 
local filters = {3,5}
local nfilter = {1,2,3}
-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1, k = 3, f = 3 }

if dotune then
	local tunefor = 8 -- full size of the matrix
	--change for 10 later
	local harness = require("lib/matrixtestharness")
		for _,b in ipairs(blocksizes) do
			for _,rm in ipairs(regblocks) do
				for _,rn in ipairs(regblocks) do
					for _,v in ipairs(vectors) do
						for _,f in ipairs(nfilter) do
							for _,k in ipairs(filters) do
								-- same until here
								local my_conv = genconvolution(b,1,rm,rn,v,k,k,f)
								-- local my_conv = gennaiveconv()
								if my_conv then
									print(b,rm,rn,v,k,f)
									my_conv:compile()
									-- bellow line makes do not need boundary cases (image multiple of blocksize)
									local i = math.floor(tunefor / b) * b
									local curr_gflops = 0
									local ctyp
									local correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k,f, function(M,N,K,L,A,B,C)
				                    	my_conv(nil,M,N,K,L,1.0,A,M,N,B,L,C,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2 
									end)
									if not correct then	print("<error>")  break  end
									print(i,unpack (exectimes))
									local curr_gflops = exectimes[1]
									-- print(curr_gflops) -- print analysis 
									if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
										best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v, k = k, f = f }
										terralib.tree.printraw(best)
									end
								end
							end
						end
					end
				end
			end
		end
end
local my_convolution = gennaiveconv()
-- local my_convolution = genconvolution(best.b,1,best.rm,best.rn,best.v)
terralib.saveobj("my_conv.o", {my_convolution = my_convolution})