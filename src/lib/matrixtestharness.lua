local ffi = require("ffi")
if ffi.os == "Windows" then
	return
end
local adjust = 0

ffi.load("/System/Library/Frameworks/Accelerate.framework/Accelerate")

function CalcTime(fn)
	local begin = terralib.currenttimeinseconds()
	local current
	local times = 0
	repeat
		fn()
		current = terralib.currenttimeinseconds()
		times = times + 1
	until (current - begin) > 0.2
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

function MTH.timefunctions(typstring,M,N,K,L,...)

	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local A = ffi.new(ctyp,M*N) 
	local B = ffi.new(ctyp,K*L)

	--simple example
	for m = 0, M-1 do
		for n = 0, N-1 do
			if( n == 0 or m == 0 or m == M-1 or n == N-1 ) then
				A[m*N + n] = 0
			else
				A[m*N + n] = (m-1)*(N-2) + n
			end
		end
	end

	B[0] = 1; B[1] = 2; B[2] = 1;
	B[3] = 0; B[4] = 0; B[5] = 0;
	B[6] = -1;B[7] =-2; B[8] = -1;

	-- randomizing
	-- for m = 0, M-1 do
	-- 	for n = 0, N-1 do
	-- 		A[m*N + n] = math.random(0,9)
	-- 	end
	-- end

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
		 -- C[j] = -1 -- for DGEMM
		end	
		Cs[i] = C
	end

	-- local C = ffi.new(ctyp,M*N)
	-- for j = 0, M * N - 1 do 
	-- 	C[j] = 0
	-- end	


	-- compute 
	local results = {}
	for i,fn in ipairs(fns) do
		local C = Cs[i]
		local tocall = function() fn(M,N,K,L,A,B,C) end
		tocall()
		CalcTime(tocall) -- execution time
		-- printMatrix(A,M,N)
		-- printMatrix(B,K,L)
		-- printMatrix(C,M,N)
		results[i] = M*N*K*L*2.0*1e-9 / CalcTime(tocall) -- gflop
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