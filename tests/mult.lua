local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
local number = double
local alignment = 8
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

local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end

local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)


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

function gencmulkernel(NB, RM, RN, V, prefetch)
	local M,N
	local A,B,C,mm,nn = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn")
	local ldc = symbol("ldc")
	local a,b,c,caddr = symmat("a",RM,RN), symmat("b",RM,RN), symmat("c",RM,RN), symmat("caddr",RM,RN)
	local load,storec,calcc = terralib.newlist(),terralib.newlist(),terralib.newlist()
	local VP = &vector(double,V)

	for m = 0, RM-1 do
		for n = 0, RN - 1 do
            load:insert(quote
            	var [caddr[m][n]] = C + m*ldc + n*V
            	var [c[m][n]] : vector(double,V) = unalignedload(VP([caddr[m][n]]))
                var [a[m][n]] : vector(double,V) = unalignedload(VP(&A[m*ldc + n*V]))
                var [b[m][n]] : vector(double,V) = unalignedload(VP(&B[m*ldc + n*V]))
            end)
            storec:insert(quote
                unalignedstore(VP([caddr[m][n]]),[c[m][n]])
            end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN - 1 do
			calcc:insert(
				quote
					[c[m][n]] = [a[m][n]] * [b[m][n]]
				end
			)
		end
	end

	-- NB must be divisible by RN*V
	return terra([A] : &double, [B] : &double, [C] : &double, [ldc] : int) 
		for [mm] = 0, NB, RM do
			for [nn]=0, NB, RN*V do
				[load];
				llvmprefetch(A + ldc*ldc,0,3,1);
				[calcc];
				[storec];
				A = A + RN*V
				B = B + RN*V
				C = C + RN*V
			end 

			A = A + RM * ldc - NB
			B = B + RM * ldc - NB
            C = C + RM * ldc - NB
		end
	end
end

function gencmul(NB,NBF,RM,RN,V)
	local NB2 = NB * NBF																																			
	local l1cmul = gencmulkernel(NB, RM, RN, V, false)

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, 
		ldc : int) 

         [ blockedloop(N,M,{NB2,NB},
                function(m,n) 
                return quote
                    var MM,NN = min(M-m,NB),min(N-n,NB)
                    var isboundary = MM < NB or NN < NB
                    var AA,BB,CC = A + (m*ldc + n), B + (m*ldc + n), C + (m*ldc + n)
                    l1cmul(AA,
                     BB,
                     CC,
                     ldc)
                end end)  
            ]       
	end
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

	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2,RN+2), symmat("b",K,L), symmat("c",RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y = symbol("x"), symbol("y")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM+1 do
		for n = 0, RN+1 do
			loadA:insert(quote
					var [a[m][n]] = vecload(A, m*ldc + n*V)
			end)
			if(m>=0 and m<RM and n>=0 and n<RN) then
				loadc:insert(quote
						var [c[m][n]] = alpha * vecload(C, (m+1)*ldc + (n+1)*V)
				end)
				storec:insert(quote
					vecstore(C, (m+1)*ldc + (n+1), [c[m][n]])
				end)
			end
		end
	end

	local calcc = terralib.newlist()

	-- load full kernel
	for  k=0, K-1 do
		for l = 0, L-1 do
			loadkernel:insert(quote
				var [b[k][l]] = vecload(B, k*ldb + l*V)
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
							if([mm] + m < NB-1 and [nn] + n < NB-1) then
								--remeber that taking the pos a[x+1][y+1], e.g. a[0][0] menas take a[-1][-1] necessary for c[0][0]
								[c[m][n]] = [c[m][n]] + [a[x+1][y+1]] * [b[k][l]]
							end
						end
					)
				end
			end
		end
	end

	return terra([A] : &double, [B] : &double, [C] : &double, [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		for [mm] = 1, NB-1, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN*V,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=1, NB-1, RN*V do
				[loadc];
				[loadkernel];
				llvmprefetch(A + sda*lda,0,3,1);
				[loadA];
				[calcc];
				[storec];
				A = A + RN*V
				C = C + RN*V
			end
			if (((NB-2)/(RN*V)) * (RN*V) + 1 < NB-1) then
				var offset = (((NB-2)/(RN*V)) * (RN*V) + 1) + (RN*V)  - (NB-1)
				A = A - offset
				C = C - offset
			end
			-- jump of two (final border one line, initial border next line)
			-- It is two because the kernel is 3, it would change for different kernel
			C = C + 2
			A = A + 2

			A = A + RM * ldc - NB
			C = C + RM * ldc - NB
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

function genconvolution(NB,NBF,RM,RN,V)
	-- register blocking does not need to be a a multiple of the blocksize anymore
	-- if not isinteger(NB/(RN*V)) or not isinteger(NB/RM) then
		-- return false
	-- end

	--5 times NB minimum by dgemm
	--local NB2 = NBF * NB


	local NB2 = NB * NBF

	local l1conv0 = genkernel(NB, RM, RN, 1, false, 3, 3, false) -- no prefetch, no boundary
	local l1conv0b = genkernel(NB, RM, RN, 1, false, 3, 3, true)

	return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
		alpha : double, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
		ldc : int, kCenterX: int, kCenterY: int) 

        var args = arrayof(int,0)
         [ blockedloop(N,M,{NB2,NB},
                function(m,n) 
                return quote
                    var MM,NN = min(M-m,NB),min(N-n,NB)
                    var isboundary = MM < NB or NN < NB
                    var AA,CC = A + (m*lda + n),C + (m*ldc + n)
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
                end end)  
            ]       
	end
end

local blocksizes = {32,48,56,64,72,80}
local regblocks = {1,2,4,8} -- blocksizes must be divisible by RN*V
local vectors = {1,2,4}

-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

if dotune then
	-- local tunefor = 1024
	local tunefor = 1024 -- full size of the matrix
	--change for 10 later
	local harness = require("mult-matrixtestharness")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocks) do
			for _,rn in ipairs(regblocks) do
				for _,v in ipairs(vectors) do
						-- same until here
					local my_cmul = gencmul(b,1,rm,rn,v)
					if my_cmul then
						print(b,rm,rn,v)
						my_cmul:compile()
						local i = math.floor(tunefor / b) * b
						local curr_gflops = 0
						local ctyp
						local correct, exectimes = harness.timefunctionsCMUL(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C)
	                    	my_cmul(nil,M,N,A,B,C,N)
						end)
						if not correct then	print("<error>")  break  end
						print(i,unpack (exectimes),"[OK]")
						local curr_gflops = exectimes[1]
						-- print(curr_gflops) -- print analysis 
						if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
							best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v }
							terralib.tree.printraw(best)
						end
					end
				end
			end
		end
	end
end

local my_convolution = genconvolution(best.b,1,best.rm,best.rn,best.v)
terralib.saveobj("my_conv.o", {my_convolution = my_convolution})