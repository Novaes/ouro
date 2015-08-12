require "lib/fftkernels"
require "lib/mthreads"
local MT = terralib.includec("pthread.h")
local cmath = terralib.includec("math.h")
local cstdio = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
local number = double
local alignment = 8
local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})
TWOPI = constant(6.28318530717959)

-- local cbit = terralib.includecstring[[
-- #include<math.h>
-- int TWOPI = 6.28318530717959;
-- void four1(double data[], int nn, int isign)
-- {
--     int n, mmax, m, j, istep, i;
--     double wtemp, wr, wpr, wpi, wi, theta;
--     double tempr, tempi;
    
--     n = nn << 1;
--     j = 1;
--     for (int i = 1; i < n; i += 2) {
-- 	if (j > i) {
-- 	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
-- 	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
-- 	}
-- 	m = n >> 1;
-- 	while (m >= 2 && j > m) {
-- 	    j -= m;
-- 	    m >>= 1;
-- 	}
-- 	j += m;
--     }
--     mmax = 2;
--     while (n > mmax) {
-- 	istep = 2*mmax;
-- 	theta = TWOPI/(isign*mmax);
-- 	wtemp = sin(0.5*theta);
-- 	wpr = -2.0*wtemp*wtemp;
-- 	wpi = sin(theta);
-- 	wr = 1.0;
-- 	wi = 0.0;
-- 	for (m = 1; m < mmax; m += 2) {
-- 	    for (i = m; i <= n; i += istep) {
-- 		j =i + mmax;
-- 		tempr = wr*data[j]   - wi*data[j+1];
-- 		tempi = wr*data[j+1] + wi*data[j];
-- 		data[j]   = data[i]   - tempr;
-- 		data[j+1] = data[i+1] - tempi;
-- 		data[i] += tempr;
-- 		data[i+1] += tempi;
-- 	    }
-- 	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
-- 	    wi = wi*wpr + wtemp*wpi + wi;
-- 	}
-- 	mmax = istep;
--     }
-- }
-- ]]

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

terra printCMatrix(str : rawstring, M : int, N : int, A : &number)
	cstdio.printf("%s\n",str)
	for i=0, M*N do
		if i ~= 0 and i % N == 0 then
			cstdio.printf("\n")
		end
		cstdio.printf(" %1.f  %1.f ",A[2*i], A[2*i + 1])
	end
	cstdio.printf("\n\n")
end

terra printMatrix(m :&number, rows : int,columns : int)
  	for i=0,rows do
  		for j=0,columns do
      		m[i*columns + j]
    	end
    	cstdio.printf("\n")
  end
  cstdio.printf("\n")
end

--[[ Use less memory for temp? ]]
function gentranskernel(M , N)
	local A, a, temp = symbol("A"), symmat("a",2*M*N), symmat("temp",M*N)
	local base = symbol("base")
	local load, calc, store = terralib.newlist(), terralib.newlist(), terralib.newlist()
	local lda = 2*N
	local tindex = 0
	local s
	for i=0,M-1 do
		local s = i+1
		for j = 2*i + 2, lda-1, 2 do
			local ii = i*lda + j
			local jj = s*lda + 2*i
			load:insert(quote -- load
				var [a[ii]] : double = A[ii + base]
				var [temp[tindex]] : double = [a[ii]]
				var [a[jj]] : double = A[jj + base]
			end)
			calc:insert(quote -- swap
				[temp[tindex]] = [a[ii]]
				[a[ii]] = [a[jj]]
				[a[jj]] = [temp[tindex]]
			end)
			store:insert(quote -- store back
				A[ii] = [a[ii]]
				A[jj] = [a[jj]]
			end)
			s = s + 1
			tindex = tindex + 1
		end
	end

	return terra([A] : &double, [base] : int)
		-- for real and complex parts, base will decide which one it is
		for i=0, 2 do
			[load];
			[calc];
			[store];
			base = base + 1
		end

	end
end

-- [[ WARNING: assume square matrix, so it works for NBxNB.
--  However, if it uses not an exactly number of NBxNB blocks it will not work (same GEMM missing case).
--  NB: l1 level block. RM and RN: register size dimension blocks ]]
function gencmulkernel(NB, RM, RN, prefetch)
	local M,N
	local A,B,C,mm,nn = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn")
	local ldc = symbol("ldc")
	local a,b,c,caddr = symmat("a",RM,RN), symmat("b",RM,RN), symmat("c",RM,RN), symmat("caddr",RM,RN)
	local loadA,loadBC,storec,calcc = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	if prefetch then 
        loadA:insert(quote
        	 llvmprefetch(A + ldc*ldc,0,3,1);
        end)
	end

	for m = 0, RM-1 do
		for n = 0, RN - 1 do
			loadA:insert(quote
                var [a[m][n]] = A[m*ldc + n]
            end)
			loadBC:insert(quote
            	var [b[m][n]] = B[m*ldc + n]
            	var [c[m][n]] = C[m*ldc + n]
            end)
            storec:insert(quote
                C[m*ldc + n] = [c[m][n]]
            end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN - 1, 2 do
			calcc:insert(
				quote
					[c[m][n]] = [a[m][n]] * [b[m][n]] - [a[m][n+1]] * [b[m][n+1]]
					[c[m][n+1]] = [a[m][n]] * [b[m][n+1]] + [b[m][n]] * [a[m][n+1]]
				    
				    -- turn around for possible turn around error on IFFT (BUG; check twidle factors)
				    if [c[m][n+1]] ~= 0 then [c[m][n+1]] = -1*[c[m][n+1]] end
				end
			)
		end
	end

	-- NB must be divisible by RN
	return terra([A] : &double, [B] : &double, [C] : &double, [ldc] : int)
		for [mm] = 0, NB, RM do
			for [nn]=0, 2*NB, RN do
				[loadBC];
				[loadA];
				[calcc];
				[storec];
				A = A + RN
				B = B + RN
				C = C + RN
			end 

			A = A + RM * ldc - NB
			B = B + RM * ldc - NB
            C = C + RM * ldc - NB
		end
	end
end

function gencmul(NB,NBF,RM,RN)
	if not isinteger(NB/RM) or not isinteger(NB/(2*RN)) then
		return false
	end

	local NB2 = NB * NBF
	local l1cmul = gencmulkernel(NB, RM, 2*RN, true)

	return terra(M : int, N : int, A : &double, B : &double, C : &double, 
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

terra transposeCMatrix(M : int, N : int, d : &number)
	var s : int = 0
	var temp : number
	-- because of imaginary part stride is 2*N
	var NN = 2*N
	-- down columns from 0 to columns - 2 along the diagonal
	for i=0,M-1 do
		s = 2*(i+1) -- next elemenet after diagonal elem
		for j=i+1,N do
			d[i*NN + s], d[j*NN + 2*i] = d[j*NN + 2*i], d[i*NN + s] -- real part
			d[i*NN + s+1], d[j*NN + (2*i)+1] = d[j*NN + (2*i)+1], d[i*NN + s+1] -- imaginary part
			s = s + 2
		end
	end
end

terra naivecmul(A : &number, B : &number, C : &number, M : int, N : int)
	for i=0, M do
		for j=0, 2*N, 2 do
			C[i*N + j] = A[i*N + j] * B[i*N + j] - A[i*N + j+1] * B[i*N + j+1] 
			C[i*N + j+1] = A[i*N + j] * B[i*N + j+1] + A[i*N + j+1] * B[i*N + j]
		end
	end
end

function genNumRecipesConv()

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, ldc : int)
		
		-- (1) 2D image convolution
		for j=0,M do
			var base = j*N*2
			four1(&A[base]-1, M, 1)
		end
		-- llvmprefetch(A + ldc*ldc,0,3,1);
		transposeCMatrix(M, N, A)
		-- l1transpose(A,0)
		for j=0,N do
			var base = j*M*2
			four1(&A[base]-1, N, 1)
		end
		transposeCMatrix(M, N, A)
		-- l1transpose(A,0)

		-- (2) 2D kernel convolution
		for j=0,M do
			var base = j*N*2
			four1(&B[base]-1, M, 1)
		end
		transposeCMatrix(M, N, B)
		-- l1transpose(B,0)
		for j=0,N do
			var base = j*M*2
			four1(&B[base]-1, N, 1)
		end
		transposeCMatrix(M, N, B)
		-- l1transpose(B,0)

		-- (3) point-wise multiplication
		printCMatrix("[A]:",M,N,A)
		printCMatrix("[B]:",M,N,B)
		naivecmul(A,B,C,M,N)
        printCMatrix("[C]:",M,N,C)

		-- (4) calculate 2D IFFT of the result
        for j=0,M do
			var base = j*N*2
			four1(&C[base]-1, M, -1)
		end
		transposeCMatrix(M, N, C)
		-- l1transpose(B,0)
		for j=0,N do
			var base = j*M*2
			four1(&C[base]-1, N, -1)
		end
		transposeCMatrix(M, N, C)

		for i=0,M*N do
			C[2*i] = C[2*i] / (M*N)
			C[2*i+1] = C[2*i+1] / (M*N)
		end
	end
end

function genfastconvMT(NB,NBF,RM,RN,V,DIMSIZE,thsize)
	local NB2 = NB * NBF

	local fullcmul = gencmul(NB,NBF,RM,RN,V)
	-- local l1transpose = gentranskernel(NB/2,NB/2) -- b/2 b/2
	local l1fft = genfftkernel(DIMSIZE,1)
	local l1invfft = genfftkernel(DIMSIZE,-1)

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, 
		ldc : int)

		-- Thread management
		var thwork = M
		var thr = thsize
		while (cmath.floor(thwork / thr) == 0) do
			thr = thr - 1
		end
		var tmp : number = thwork 
		var ftaskspth : number = tmp / thr
		var taskspth : int = thwork / thr
		
		-- adjust in case of not perfect division
		var rr : int = 0
		if ftaskspth ~= taskspth then
			rr = thwork%thr
			if rr == 0 then
				rr = (ftaskspth-taskspth)*thr
			end
		end

		var threads : MT.pthread_t[thsize]
		var pkgs : L1Package[thsize]
		var added = 0
		for i=0, thr do
			if rr > 0 then 
				added = 1
				rr = rr - 1
			end
			pkgs[i]:init(M, A, B, C, l1fft, taskspth + added, 0) -- 0 means image
			added = 0
		end

		-- (1) 2D image convolution
		var count = 0
		for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end

		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	-- cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end
		transposeCMatrix(M, N, A)


		
		count = 0
		for i=0, thr do pkgs[i]:clear()	end
		for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end
		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	-- cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end

		transposeCMatrix(M, N, A)

		-- (2) 2D kernel convolution
		-- printCMatrix("[B]:",M,N,BB)
		for i=0, thr do
			pkgs[i]:usefilter() -- it clear already
		end
		count = 0
		for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end
		
		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	-- cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end
		
		transposeCMatrix(M, N, B)
		
		for i=0, thr do	pkgs[i]:clear()	end
		count = 0
		for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end

		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	-- cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end
		transposeCMatrix(M, N, B)

		-- (3) point-wise multiplication
		-- printCMatrix("[A 2D convolved]:",M,N,A)
		-- printCMatrix("[B 2D convolved]:",M,N,B)
		fullcmul(M,N,A,B,C,N)
		-- printCMatrix("[Pointwise multiplied]:",M,N,C)

		
    	-- (4) calculate 2D IFFT of the result
    	for i=0, thr do
			pkgs[i]:useoutput() -- it clear already
		end
		count = 0

        for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end
		
		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end

		transposeCMatrix(M, N, C)
		
		for i=0, thr do	
			pkgs[i]:clear()	
		end
		count = 0

		for j=0,M do
			var base = j*N*2
			if pkgs[count]:addblock(base) then
				if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
					-- cstdio.printf("Thread #%u creation error",threads[count])
				end
				-- cstdio.printf("---> thread launched: %d\n",count)
				count = count + 1
			end
		end

		for i=0,thr do
	    	if MT.pthread_join(threads[i],nil) ~= 0 then
	        	-- cstdio.printf("Thread #%u join error\n",i) 
	    	end
		end

		transposeCMatrix(M, N, C)
		-- printCMatrix("[C]:",M,N,A)
	end

end

function genfastconv(NB,NBF,RM,RN,V,FULLIMAGESIZE)
	local NB2 = NB * NBF
	
	-- not only a cmul kernel (NBxNB), but a full point-wise multiplication
	-- only kernel imp
	-- if not isinteger(NB/RM) or not isinteger(NB/(2*RN)) then return false end
	local fullcmul = gencmul(NB,NBF,RM,RN,V)
	-- local l1transpose = gentranskernel(NB/2,NB/2) -- b/2 b/2

	local l1fft = genfftkernel(FULLIMAGESIZE,1)
	local l1invfft = genfftkernel(FULLIMAGESIZE,-1)

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, 
		ldc : int)

		-- (1) 2D image convolution
		for j=0,M do
			var base = j*N*2
			l1fft(&A[base])
			-- printCMatrix("[A]:",M,N,A)
		end
		-- llvmprefetch(A + ldc*ldc,0,3,1);
		transposeCMatrix(M, N, A)
		-- printCMatrix("[A transposed]:",M,N,A)
		for j=0,M do
			var base = j*N*2
			l1fft(&A[base])
			-- printfrintCMatrix("[A]:",M,N,A)
		end
		transposeCMatrix(M, N, A)

		-- (2) 2D kernel convolution
		-- -- printCMatrix("[B]:",M,N,BB)
		for j=0,M do
			var base = j*N*2
			l1fft(&B[base])
		end
		transposeCMatrix(M, N, B)
		-- l1transpose(B,0)
		for j=0,M do
			var base = j*N*2
			l1fft(&B[base])
		end
		transposeCMatrix(M, N, B)

		-- (3) point-wise multiplication
		-- printCMatrix("[A 2D convolved]:",M,N,A)
		-- printCMatrix("[B 2D convolved]:",M,N,B)
		fullcmul(M,N,A,B,C,N)
		-- printCMatrix("[Pointwise multiplied]:",M,N,C)

    	-- (4) calculate 2D IFFT of the result
    	-- printCMatrix("[C]:",M,N,A)
        for j=0,M do
			var base = j*N*2
			l1invfft(&C[base])
			-- printCMatrix("[C]:",M,N,C)
		end
		transposeCMatrix(M, N, C)
		-- l1transpose(BB,0)
		for j=0,M do
			var base = j*N*2
			l1invfft(&C[base])
		end
		transposeCMatrix(M, N, C)
	end
end
 
function genfftkernel(NFFT, signal)
	local A, ker, K_W = symbol("A"), symmat("ker",2,2), symbol("K_W")
	local stotal, rest, s = symbol("stotal"), symbol("rest"), symbol("s")
	local NELEMS, k, Ns = symbol("NELEMS"), symbol("k"), symbol("Ns")
	local skernel, base, ublocks =  symbol("skernel"), symbol("base"), symbol("ublocks")
	local exec, bitreversal = terralib.newlist(), terralib.newlist()	
	local ker = {}

	ker[0], ker[1] = {}, {}
	ker[0][0], ker[0][1] = 4, 2 -- kernel of 4 points; 2 stages
	ker[1][0], ker[1][1] = 2, 1 -- kernel of 2 points; 1 stage
 	stotal = math.floor(math.log(NFFT)/math.log(2)) -- log2(NFFT) to #stages
	rest = stotal
	s, k, Ns = 0,0,1 -- s: current stage
	NELEMS = NFFT
	
	repeat
		-- #stages that type of kernel absorb is lower or equal lacking #stages
		if ker[k][1] <= rest then 
			K_W = ker[k][0]
			skernel = ker[k][1]
			base  = 0
			ublocks = NELEMS/K_W
			rest = rest - skernel
			for i=0,ublocks-1 do -- loop over big ublocks
				base = i*(Ns*K_W*2) -- iterate over the big ublocks, K_W*2 each elem size
				for j=0, Ns-1 do -- inside each block
					if k == 0 then
						exec:insert(quote
							FFT_4(base + 2*j,A,Ns,NFFT,signal)
						end)
					elseif k == 1 then
						exec:insert(quote
							FFT_2(i*Ns+j,base + 2*j,A,Ns,NFFT,signal)
						end)
					end
				end
			end
			-- jump stages	(min. skernel == 1 and K_W = 2)
			s = s + skernel
			Ns = Ns * K_W
			NELEMS = NELEMS / K_W 
		else
			k = k + 1
		end
	until rest == 0
	
	bitreversal:insert(quote
		cbit.reversal(NFFT,A)
	end)

	return terra([A] : &double)
		[exec];
		[bitreversal];
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

function blockedComplexloop(M,N,blocksizes,bodyfn)
  local function generatelevel(n,ii,jj,bb0,bb1)
    if n > #blocksizes then
      return bodyfn(ii,jj)
    end
    local blocksize = blocksizes[n]
    return quote for i = ii,min(ii+bb0,M),blocksize do
                   for j = jj,min(jj+bb1,2*N),2*blocksize do
                        [ generatelevel(n+1,i,j,blocksize,blocksize) ]
           		end end end
  end
  return generatelevel(1,0,0,M,N)
end

local blocksizes = {8--[[16,24,32,40,48,56,64,1024]]}
local regblocks = {1}--{1,2,4} -- blocksizes must be divisible by RN*V
local vectors = {1--[[1,2,4,8]]}
local nthread = {2}
-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

-- kernel dependent 
local tunefor = 128 -- full size of the matrix
if dotune then
	-- local tunefor = 1024
	--change for 10 later
	local harness = require("lib/matrixtestharness")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocks) do
			for _,rn in ipairs(regblocks) do
				for _,v in ipairs(vectors) do
					for _,t in ipairs(nthread) do
						-- FOR FAST CONVOLUTION THE KERNEL IS DEPENDENT OF THE IMAGE SIZE
						local i = math.floor(tunefor / b) * b
						-- local my_fastconv = genNumRecipesConv()
						local my_fastconv = genfastconvMT(b,1,rm,rn,v,tunefor,t)
						-- local my_fastconv = genfastconv(b,1,rm,rn,v,tunefor)
						-- local my_cmul = gencmul(b,1,rm,rn,v)
						if my_fastconv then
							print(b,rm,rn,v)
							my_fastconv:compile()
							local curr_gflops = 0
							local ctyp
							local caltime, correct, exectimes = harness.timefunctionsFFT(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C) 
		                    	my_fastconv(nil,M,N,A,B,C,N)                  	
							end)
							-- local correct, exectimes = harness.timefunctionsCMUL(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C) 
		                    	-- my_cmul(nil,M,N,A,B,C,N)	                    	
							-- end)
							-- if not correct then	("<error>")  break  end
							-- print(i,unpack (exectimes))
							best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v }
							print(caltime)
							-- local curr_gflops = exectimes[1]
							-- if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
							-- 	best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v }
							-- 	terralib.tree.printraw(best)
							-- end
						end
					end
				end
			end
		end
	end
end

local my_convolution = genfastconv(best.b,1,best.rm,best.rn,best.v,tunefor)
terralib.saveobj("../bin/my_numconv.o", {my_numconv = my_convolution})
