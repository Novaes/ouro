local gemm = require("gemm")
local cstdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
-- Set number to float in case of Single Float Point tests
local number = double


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

-- terra min(a : int, b : int)
-- 	return terralib.select(a < b, a, b)
-- end

-- function blockedloop(M,N,blocksizes,bodyfn)
--   local function generatelevel(n,ii,jj,bb0,bb1)
--     if n > #blocksizes then
--       return bodyfn(ii,jj)
--     end
--     local blocksize = blocksizes[n]
--     return quote for i = ii,min(ii+bb0,M),blocksize do
--                    for j = jj,min(jj+bb1,N),blocksize do
--                         [ generatelevel(n+1,i,j,blocksize,blocksize) ]
--            		end end end
--   end
--   return generatelevel(1,0,0,M,N)
-- end

function gennaiveconv()

	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		sdc : int, ldc : int, kCenterX: int, kCenterY: int, depth : int) 
		var dimKer : int = K*L
		var kCenterX : int = K/2
		var kCenterY : int = L/2
		var e = 0
		for d=0, depth do
			var baseKer = dimKer * d
			for i=kCenterX, M-kCenterX do
			    for j=kCenterY, N-kCenterY do
			    	C[e] = 0
		        	for m=0, K do
		          		for n=0,L do
				            -- boundaries
				            var ii: int = i + (m - kCenterY)
				            var jj: int = j + (n - kCenterX)
				            if ii>=0 and ii<M and jj>=0 and jj<N then
				              C[e] = C[e] + A[ii * N + jj] * B[baseKer + m * L + n]
				            end
			        	end
			      	end
			      	e = e+1
			    end
			end
		end
	end
end

function genlowimage()
   
	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, A : &number, AA : &number) 
		var kAenterX : int = K/2
		var kAenterY : int = L/2
		var e = 0
		for i=kAenterX, M-kAenterX do
		    for j=kAenterY, N-kAenterY do
	        	for m=0, K do
	          		for n=0, L do
			            var ii: int = i + (m - kAenterY)
			            var jj: int = j + (n - kAenterX)
			            if ii >= 0 and ii < M and jj >= 0 and jj < N then
			            	AA[e] = A[ii * N + jj]
			            else
			            	AA[e] = 0
			            end
			            e = e + 1
		        	end
		      	end
		    end
		end
	end
end

-- Different types according to Caffe con Troll can be implemented
function genlowkernel(loweringtype)

	return terra(gettime : {} -> number, Bs : &number, K : int, L : int, BBs : &number, Klow : int, Llow : int)
		if loweringtype == 1 then
			var base = Klow
			var count = 0
			for d=0,Llow do
				for i=0,Klow do
					BBs[i*Llow + d] = Bs[base*d + i]
				end
			end
		end
	end

end

function liftresult(loweringtype) 
	
	return terra(gettime : {} -> number, C : &number, M : int, N : int, Clow : &number, Mlow : int, Llow : int)
		if loweringtype == 1 then
			var count = 0
			for i=0,Llow do
				for j=0,Mlow do
					-- iterate over all the columns first
					C[count] = Clow[j*Llow + i]
					count = count + 1
				end
			end
		end
	end
end

terra printlifted(A : &number, M : int, N : int, depth : int)
	var dim = M*N
	for d=0,depth do
		var base = M*N
		for i=0,M do 
			for j=0,N do
				cstdio.printf("%1.f ", A[base + i*N + j])
			end
			cstdio.printf("\n")
		end
		cstdio.printf("\n")
	end
end

terra printM(A : &number, M : int, N : int)
	for i=0,M do 
		for j=0,N do
			cstdio.printf("%1.f ", A[i*N + j])
		end
		cstdio.printf("\n")
	end
	cstdio.printf("\n")
end

function genconvolution(NB,NBF,RM,RN,V,K,L)
	-- Lowering type according to CcT (http://arxiv.org/abs/1504.04343)
	local ltype = 1

	-- Lower the image is the 2D direct method without the multiplication
	-- These are not kernels, they are full operations


	local my_dgemm = generatedgemm(64,5,RM,RN,V)
	
	-- local loweredImg = genLowImage(NB, NBF, RM, RN, V, K, L)
	local loweredImg = genlowimage()
	local loweredKer = genlowkernel(ltype)
	local liftedResult = liftresult(ltype)

	return terra(gettime : {} -> number, 
		A : &number,
		M : int, N : int, K : int, L: int, 
		alpha : number, 
		B : &number, ldb : int, 
		C : &number, sdc : int, ldc : int, 
		kCenterX: int, kCenterY: int, 
		depth : int,
		AA : &number, BB : &number, CC : &number)

		-- MlowxNlow * KlowxLlow, where Nlow == Llow ---> Mlow x Llow 
		var Mlow, Nlow, Klow, Llow = sdc*ldc, K*L, K*L, depth

		-- (1) lower
		loweredKer(nil,B,K,L,BB,Klow,Llow)
		-- printM(AA,Mlow,Nlow)
		
		loweredImg(nil,M,N,K,L,A,AA)
		printM(BB,Klow,Llow)
		
		-- (2) gemm
		my_dgemm(nil,Mlow,Llow,Klow,1.0,AA,Nlow,BB,Llow,0.0,CC,Llow)

		printM(CC,Mlow,Llow)
		
		-- (3) lifting
		liftedResult(nil,C,M,N,CC,Mlow,Llow)
		-- printlifted(C,sdc,ldc,depth)
	end
end


-- function genLowImageOPT(NB,NBF,RM,RN,V,K,L)
-- 	-- register blocking does not need to be a a multiple of the blocksize anymore
-- 	if not isinteger(NB/RN) or not isinteger(NB/RM) then
-- 		print("3rd level blocksizes must be a multiple of the 2nd")
-- 		return false
-- 	end

-- 	--5 times NB minimum by dgemm
-- 	--local NB2 = NBF * NB
-- 	local NB2 = NB * NBF
-- 	local l1conv0 = genkernel(NB, RM, RN, true, K, L, false)
-- 	local l1conv0b = genkernel(NB, 1, 1, true, K, L, true)

-- 	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
-- 		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
-- 		sdc : int, ldc : int, kCenterX: int, kCenterY: int, depth : int, batches : int)


-- 		M = M-kCenterX
-- 		N = N-kCenterY

-- 		for mm = kCenterX,M,NB2 do
-- 			for nn = kCenterY,N,NB2 do
-- 				for m = mm,min(mm+NB2,M),NB do
-- 					for n = nn,min(nn+NB2,N),NB do
-- 					 	var MM,NN = min(M-m,NB),min(N-n,NB)
-- 	                    var isboundary = MM < NB or NN < NB
-- 	                    var AA,CC = A + ((m-kCenterX)*lda + (n-kCenterY)),C + ((m-kCenterX)*ldc + (n-kCenterY))
-- 	                    if isboundary then -- do not enter here YET
-- 	                    	l1conv0b(AA,
-- 	                         B,
-- 	                         CC,
-- 	                         sda,lda,ldb,sdc,ldc,0,MM,NN)
-- 	                    else
-- 	                        l1conv0(AA,
-- 	                         B,
-- 	                         CC,
-- 	                         sda,lda,ldb,sdc,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
-- 	                    end
-- 					end
-- 				end
-- 			end
-- 		end

--          -- [ blockedloop(M,N,{NB2,NB},
--          --        function(m,n) 
--          --        return quote
--          --            var MM,NN = min(M-m,NB),min(N-n,NB)
--          --            var isboundary = MM < NB or NN < NB
--          --            var AA,CC = A + (m*lda + n),C + (m*ldc + n)
--          --            if isboundary then -- do not enter here YET
--          --             l1conv0b(AA,
--          --                 B,
--          --                 CC,
--          --                 sda,lda,ldb,ldc,0,MM,NN)
--          --            else
--          --                l1conv0(AA,
--          --                 B,
--          --                 CC,
--          --                 sda,lda,ldb,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
--          --            end
--          --        end end)  
--          --    ]       
-- 	end
-- end

-- Different blocksizes for the same result implies in padding overheading 
-- ending in s means SIZE
-- starting with n, means NUMBER
local blocksizes = {8}
local regblocks = {1}--,2,4,8}

-- local vectors = {1,2,4,8,16}
local vectors = {1}
local filters = {3}--,5}
local nfilter = {2}--10,100,200,1024}--,2,3}
-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1, k = 3, f = 3 }
local NB2 = 5

if dotune then
	-- full size of the matrix
	local tunefor = 8 
	--change for 10 later
	local harness = require("lib/matrixtestharness")
		for _,b in ipairs(blocksizes) do
			for _,rm in ipairs(regblocks) do
				for _,rn in ipairs(regblocks) do
					for _,v in ipairs(vectors) do
						for _,f in ipairs(nfilter) do
							for _,k in ipairs(filters) do
								-- same until here
								local my_conv = genconvolution(b,NB2,rm,rn,v,k,k)
								-- local my_conv = gennaiveconv()
								if my_conv then
									print(b,rm,rn,v,k,f)
									my_conv:compile()
									-- bellow line makes do not need boundary cases (image multiple of blocksize)
									local i = math.floor(tunefor / b) * b
									local curr_gflops = 0
									local ctyp
									local correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k,f, 
										function(Me,Ne,K,L,M,N,A,Bs,Cs,f,AA,BB,CC)
											-- to gennaive pass the #kernels here
				                    		my_conv(nil,A,Me,Ne,K,L,1.0,Bs,L,Cs,M,N,K/2,L/2,f,AA,BB,CC) 
				                    		-- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
									end)
									if not correct then	print("<error>") break end
									print(i,unpack (exectimes),"[OK]")
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