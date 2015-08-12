local gemm = require("gemm")
local cstdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
-- Set number to float in case of Single Float Point tests
local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})

local number = double

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

terra max(a : int, b : int)
	return terralib.select(a > b, a, b)
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
function gennaivegemm()

	return terra(gettime : {} -> number, CC : &number, AA : &number, BB : &number, M : int, N : int, K : int, ix : int, iy : int)
		for m=ix, M do
			for n=iy, N do
				CC[m*N + n] = 0
				for k=0, K do
					CC[m*N + n] = CC[m*N + n] + AA[m*K + k] * BB[k*N + n]
				end
			end
		end
	end
end

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

function genlowimage(loweringtype,NB)
   
	return terra(M : int, N : int, K : int, L: int, A : &number, AA : &number) 
		if loweringtype == 1 then
			var cx : int = K/2
			var cy : int = L/2
			var base : int
			var elem : int			
		 	for mm = cx,M-cx,NB do
				for nn = cy,N-cy,NB do
					for m = mm,min(mm+NB,M-cx) do
						for n = nn,min(nn+NB,N-cy) do
							for k=0, K do
		            for l=0, L do			      	
				            var ii: int = m + (k - cy)
							      var jj: int = n + (l - cx)
								    elem = (m-cx)*(N-2*cy)+(n-cy)
							      base = K*L
							      AA[elem*base + k*K+l] = A[ii * N + jj] 
							   end
							end
						end
					end
				end
			end
		end
	end
end

-- Different types according to Caffe con Troll can be implemented
function genlowkernel(loweringtype)

	return terra( Bs : &number, K : int, L : int, BBs : &number, Klow : int, Llow : int)
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
	
	return terra(C : &number, M : int, N : int, Clow : &number, Mlow : int, Llow : int)
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

terra printMatrix(m : &number,rows : int,columns : int,depth : int)
    var dim = rows * columns
    if depth > 1 then 

    	-- print each output
        for d=0,depth-1 do
            var base = dim * d
            for i=0,rows-1 do
                for j=0,columns-1 do
                    io.write(" " .. m[base + i*columns + j])
                end
                io.write("\n")
            end
            io.write("\n")
        end

    else
        -- usual matrix print
        for i=0,rows-1 do
            for j=0,columns-1 do
              io.write(" " .. m[i*columns + j])
            end
            io.write("\n")
        end
    end
    io.write("\n")
end

terra print2D(A : &number, M : int, N : int)
	for i=0,M do 
		for j=0,N do
			cstdio.printf("%1.f ", A[i*N + j])
		end
		cstdio.printf("\n")
	end
	cstdio.printf("\n")
end

function genLoweringConv(NBL,NB,NBF,RM,RN,V,thsize)
	-- Lowering type according to CcT (http://arxiv.org/abs/1504.04343)
	local ltype = 1
	
	-- Lower the image is the 2D direct method without the multiplication
	-- These are not kernels, they are FULL OPERATIONS
	
	-- naive gemm for borders	
	-- local my_loweredimg = genLowImage(NB, NBF, RM, RN, V, K, L)
	local my_gemmopt = generatedgemm((NB*NB),5,RM,RN,V,thsize)
	local my_naivegemm = gennaivegemm()
	local my_loweredimg = genlowimage(ltype,NBL)
	local my_loweredker = genlowkernel(ltype)
	local my_liftedresult = liftresult(ltype)

	return terra(gettime : {} -> double, 
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
		-- print2D(A,M,N)
		my_loweredimg(M,N,K,L,A,AA)
		-- print2D(AA,Mlow,Klow)
		my_loweredker(B,K,L,BB,Klow,Llow)
		-- print2D(BB,Klow,Llow)
			
		-- (2) gemm
		my_gemmopt(nil,Mlow,Llow,Klow,1.0,AA ,Klow,BB ,Llow,CC,Llow)
		-- print2D(CC,Mlow,Llow)
		
		-- (3) lifting
		my_liftedresult(C,M,N,CC,Mlow,Llow)
		-- print2D(CC,Mlow,Llow)
		-- printMatrix(C,sdc,ldc,depth)
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

		-- for mm = kCenterX,M,NB2 do
		-- 	for nn = kCenterY,N,NB2 do
		-- 		for m = mm,min(mm+NB2,M),NB do
		-- 			for n = nn,min(nn+NB2,N),NB do
		-- 			 	var MM,NN = min(M-m,NB),min(N-n,NB)
	 --                    var isboundary = MM < NB or NN < NB
	 --                    var AA,CC = A + ((m-kCenterX)*lda + (n-kCenterY)),C + ((m-kCenterX)*ldc + (n-kCenterY))
	 --                    if isboundary then -- do not enter here YET
	 --                    	l1conv0b(AA,
	 --                         B,
	 --                         CC,
	 --                         sda,lda,ldb,sdc,ldc,0,MM,NN)
	 --                    else
	 --                        l1conv0(AA,
	 --                         B,
	 --                         CC,
	 --                         sda,lda,ldb,sdc,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
	 --                    end
		-- 			end
		-- 		end
		-- 	end
		-- end

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