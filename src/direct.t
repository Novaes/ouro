require("lib/dir-mthreads")
local MT = terralib.includec("pthread.h")
local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
local cmath = terralib.includec("math.h")
local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})

local number = double

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

-- generate L1 convolution 
function genDirkernel(NB, RM, RN, V, prefetch, K, L, depth, boundary)
	local M,N, boundaryargs

	if boundary then 
		M,N = symbol(int64,"M"),symbol(int64,"N")
		boundaryargs = terralib.newlist({M,N})
	else
		boundaryargs = terralib.newlist()
		M,N = NB,NB
	end

	local cx,cy = math.floor(K/2), math.floor(L/2)
	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sdc,ldc,sda,lda,ldb = symbol("sdc"),symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2*cx,RN+2*cy), symmat("b",depth,K,L), symmat("c",depth,RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y = symbol("x"), symbol("y")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	if prefetch then 
		loadA:insert(quote
			llvmprefetch(A + (RM+2*cx)*(RN+2*cy),0,3,1);
		end)
	end

	for m = 0, RM+2*cx-1 do
		for n = 0, RN+2*cy-1 do
			loadA:insert(quote
					var [a[m][n]] = A[m*lda + n]
			end)
		end
	end

	for  i=0, depth-1 do
		for m=0,RM -1 do
			for n=0,RN -1 do
				loadc:insert(quote
						-- sda*lda*i is the base when processing multiple kernel results
						var [c[i][m][n]] = alpha * C[sdc*ldc*i  + m*ldc + n]
				end)
				storec:insert(quote
					C[sdc*ldc*i + m*ldc + n] = [c[i][m][n]]
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
					-- squared kernels: ldb*ldb; so K == L
					var [b[i][k][l]] = B[ldb*ldb*i  + k*ldb + l]
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
								[c[i][m][n]] = [c[i][m][n]] + [a[x+cx][y+cy]] * [b[i][k][l]]
							end
						)
					end
				end
			end
		end
	end

	return terra([A] : &number, [B] : &number, [C] : &number, [sda] : int, [lda] : int, [ldb] : int, [sdc] : int, [ldc] : int, [alpha] : number, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		[loadkernel];
		for [mm] = 0, M, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=0, N, RN do 
				[loadc];
				[loadA];
				[calcc];
				[storec];
				A = A + RN
				C = C + RN
			end
			
			A = A + RM * lda - M
			C = C + RM * ldc - M
		end
	end
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

function maxreuse()

	return terra(gettime : {} -> number, Me : int, Ne : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		sdc : int, ldc : int, kCenterX: int, kCenterY: int, depth : int) 

		var dimKer : int = K*L
		var dimOut : int = sdc*ldc
		var kCenterX : int = K/2
		var kCenterY : int = L/2
		var e = 0
		for i=kCenterX, Me-kCenterX do
		    for j=kCenterY, Ne-kCenterY do
		    	for d=0, depth do
					var baseKer = dimKer * d
					var baseOut = dimOut * d
			    	C[baseOut + e] = 0
		        	for m=0, K do
		          		for n=0,L do
				            -- boundaries
				            var ii: int = i + (m - kCenterY)
				            var jj: int = j + (n - kCenterX)
				            if ii>=0 and ii<Me and jj>=0 and jj<Ne then
				              C[baseOut + e] = C[baseOut + e] + A[ii * Ne + jj] * B[baseKer + m * L + n]
				            end
			        	end
			      	end
		    	end
		    	e = e+1
			end
		end
	end
end

function gennaiveconv()

	return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
		alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
		ldc : int, kCenterX: int, kCenterY: int, depth : int) 
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
 
function genDirconv(NB,NBF,RM,RN,V,K,L,NF,thsize)	
	if math.floor(NB/RN) ~= NB/RN or math.floor(NB/RM) ~= NB/RM then
		print("3rd level blocksizes must be a multiple of the 2nd")
		return false
	end

	local NB2 = NB * NBF -- NBF was 5 for gemm (5 times grater than NB)
	local l1conv0 = genDirkernel(NB, RM, RN, V, true, K, L, NF, false)
    local l1conv0b = genDirkernel(NB, 1, 1, V, false, K, L, NF, true)

    return terra(gettime : {} -> number, M : int, N : int, K : int, L: int, 
        alpha : number, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number, 
        sdc : int, ldc : int, kCenterX: int, kCenterY: int) 

		var thwork : int = ldc*ldc/(NB*NB)
		
		-- thread adjusts
		var thr : int = thsize
		-- #threads is greater than the work
		while (cmath.floor(thwork / thr) == 0) do
			thr = thr - 1
		end
		var tmp : number = thwork 
		var ftaskspth : number = tmp / thr
		var taskspth : int = thwork / thr
		
		-- adjust in case of not perfect division
		-- cstdio.printf("ftaskspth: %f\n",ftaskspth)
		-- cstdio.printf("taskspth: %d\n",taskspth)
		var rr : int = 0
		if ftaskspth ~= taskspth then
			rr = thwork%thr
			if rr == 0 then
				rr = (ftaskspth-taskspth)*thr
			end
		end

		var threads : MT.pthread_t[thsize]
        var pkgs : DirL1Package[thsize]
        var added = 0
    	for i=0, thr do
        	if rr > 0 then 
        		added = 1
				rr = rr - 1
			end
			pkgs[i]:init(NB, M, N, A, sda, lda, B, ldb, C, sdc, ldc, l1conv0, l1conv0b, taskspth+added, kCenterX, kCenterY)
			added = 0
		end
		var count = 0


		-- unroll block to be done in LUA
		M = M-kCenterX
		N = N-kCenterY
		for mm = kCenterX,M,NB2 do
			for nn = kCenterY,N,NB2 do
				for m = mm,min(mm+NB2,M),NB do
					for n = nn,min(nn+NB2,N),NB do
						-- cstdio.printf("adding to thread: %d\n",count)	
				 		if pkgs[count]:addblock(m,n) then
		                	if MT.pthread_create(&threads[count], nil, dirl1MTComputation , &pkgs[count]) ~= 0 then 
	                    		-- cstdio.printf("Thread #%u creation error",threads[count])
	                    	end
	                    	-- cstdio.printf("---> thread launched: %d\n",count)
	                    	count = count + 1
	                    end
					end
				end
			end
		end
		
		for i=0,thr do
        	if MT.pthread_join(threads[i],nil) ~= 0 then
        		-- cstdio.printf("Thread #%u join error\n",i) 
        	end
        end
    end
end
