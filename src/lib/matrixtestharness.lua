local ffi = require("ffi")
if ffi.os == "Windows" then
	return
end
local adjust = 0

-- ffi.load("/System/Library/Frameworks/Accelerate.framework/Accelerate")

function CalcTime(fn)
	local begin = terralib.currenttimeinseconds()
	local current
	local times = 0
	-- repeat
		fn()
		current = terralib.currenttimeinseconds()
		times = times + 1
	-- until (current - begin) > 0.2
	return (current - begin - adjust*times) / times 
end

-- local terra empty() end
-- adjust = CalcTime(empty) --calculate function call overhead and subtract from tests

local MTH = {}

-- find a standard convolution and make it here
-- local acc = terralib.includecstring[[
-- 	void cblas_dgemm(int, int,
--                  int, const int M, const int N,
--                  const int K, const double alpha, const double *A,
--                  const int lda, const double *B, const int ldb,
--                  const double beta, double *C, const int ldc);
-- ]]
-- local function referenceconv(M,N,K,L,A,B,C)
-- 	--acc.cblas_dgemm(101,111,111,M,N,K,1.0,A,K,B,N,0.0,C,N)
-- end

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

-- create a 1-9 with 0 padding border, used for test convolution
function fillpart(A,sM,eM,sN,eN,M,N) 
	local counter = 0
	local dim = eN-sN + 1
	local elems = dim * dim
	for m = sM, eM do 
		for n = sN, eN do
			A[m*N + n] = counter + 1
			counter = (counter + 1) % elems
		end
	end
end

function generateA(A, M, N, inX, inY)
	local nRows = M/inX - 1 -- first matrix start from 0
	local nCols = N/inY - 1
	for j=0,nRows do
		for i = 0,nCols do
			fillpart(A, 1 + j*inY, (inY-2) + j*inY, 1 + i*inX, (inX-2) + i*inX, M, N)
		end
	end
end

function naiveConvolve(out, inp, M, N, kernel, K, L)
  local kCenterX, kCenterY = math.floor(K/2), math.floor(L/2)
	for i=0, M-1 do
		for j=0, N-1 do
			out[i*N + j] = 0
		  	for m=0,K-1 do
		  		for n=0,L-1 do
			      	--boundaries
				    local ii = i + (m - kCenterY)
				    local jj = j + (n - kCenterX)
				    if ii >= 0 and ii< M and jj>=0 and jj<N then
				    	local tmp = out[i*N + j]
				    	out[i*N + j] = tmp + inp[ii*N + jj] * kernel[m*L + n] --kernel[mm+1][nn+1];
				    end
		    	end
			end
		end
  	end
end

function MTH.timefunctions(typstring,M,N,K,L,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local A = ffi.new(ctyp,M*N) 
	local B = ffi.new(ctyp,K*L)

	local CR = ffi.new(ctyp,M*N)

	--specific example A
	generateA(A, M, N, K+2, L+2) -- A has dimensions 2x the kernel and C 
	-- printMatrix(A,M,N)

	-- specific examples B
	if K == 3 then
		B[0] = 1; B[1] = 2; B[2] = 1;
		B[3] = 0; B[4] = 0; B[5] = 0;
		B[6] = -1;B[7] =-2; B[8] = -1;
	end
	if K == 5  then
		B[0] = 1; B[1] = 4; B[2] = 7; B[3] = 4; B[4] = 1;
		B[5] = 4; B[6] = 16; B[7] = 26; B[8] = 16; B[9] = 4;
		B[10] = 7; B[11] = 26; B[12] = 41; B[13] = 26; B[14] = 7;
		B[15] = 4; B[16] = 16; B[17] = 26; B[18] = 16; B[19] = 4;
		B[20] = 1; B[21] = 4; B[22] = 7; B[23] = 4; B[24] = 1;
	end
	-- printMatrix(B,K,L)
	
	printMatrix(A,M,N)
	printMatrix(B,K,L)
	
	-- randomizer A
	-- for m = 0, M-1 do
	-- 	for n = 0, N-1 do
	-- 		A[m*N + n] = math.random(0,9)
	-- 	end
	-- end

	-- randomizer B
	-- for k = 0, K-1 do
	-- 	for l = 0, L-1 do
	-- 		B[k*L + l] = math.random(0,9)
	-- 	end
	-- end

	local fns = {...}
	local Cs = {}

	-- initialize: fill the matrix C with -1
	for i,fn in ipairs(fns) do
		local C = ffi.new(ctyp,M*N)
		for j = 0, M * N - 1 do 
			C[j] = 0
		end	
		Cs[i] = C
	end

	-- compute 
	local results = {}
	for i,fn in ipairs(fns) do
		local C = Cs[i]
		local tocall = function() fn(M,N,K,L,A,B,C) end
		tocall()
		printMatrix(A,M,N)
		printMatrix(B,K,L)
		printMatrix(C,M,N)
		results[i] = M*N*K*L*2.0*1e-9 / CalcTime(tocall) -- gflop
		
		-- CORRECTNESS
		naiveConvolve(CR,A,M,N,B,K,L)
		-- ASSERT 
		
		if i ~= 1 then
			local C0 = Cs[1]
			local C1 = Cs[i]
			local c = 0
			for m = 0, M-1 do
				for n = 0, N-1 do
					if C0[c]~= C1[c] then
						return false
					end
					c = c + 1
				end
			end
		end
	end
	return true,results
end

function MTH.comparetoaccelerate(NB,myfn)
	for i = NB, math.huge, NB do
		local m,n,k = i,i,i
		io.write(("%d %d %d "):format(m,n,k))
		local s,r = MTH.timefunctions("double",m,n,k,referenceconv,myfn)
		if s then
			print(unpack(r))
		else
			print(" <error> ")
		end
		if m*n + m*k + n*k > 3*1024*1024 then
			break
		end
	end
end

return MTH