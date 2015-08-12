local ffi = require("ffi")
if ffi.os == "Windows" then
	return
end
local adjust = 0
local number = double
-- ffi.load("/System/Library/Frameworks/Accelerate.framework/Accelerate")

function CalcTime(fn)
	local begin = terralib.currenttimeinseconds()
	local current
	local times = 0
	--repeat
		fn()
		current = terralib.currenttimeinseconds()
		times = times + 1
	--until (current - begin) > 0.2
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

function printComplexMatrix(str, M, N, A)
	print(str)
	for i=0,M*N-1 do
		if i ~= 0 and i % N == 0 then
			io.write("\n")
		end
		io.write(A[2*i] .. " " .. A[2*i + 1] .. " ")
	end
	io.write("\n")
	io.write("\n")
end

-- [[ x: table in lua, so +1 base index ]]
function newcmatrix(M, x, NFFT, padding)
	local t = 1
	for i=0,NFFT - 1 do
	  M[2*i] = x[t]
	  M[2*i+1] = 0.0
	  t = t + 1
	end
end

function checkCmul(A,B,C,M,N)
	for i=0, M-1 do
		for j=0, N-1, 2 do
			if C[i*N + j] ~= A[i*N + j] * B[i*N + j] - A[i*N + j+1] * B[i*N + j+1] or 
				C[i*N + j+1] ~= A[i*N + j] * B[i*N + j+1] + A[i*N + j+1] * B[i*N + j]  then				
				return false
			end
		end
	end
	return true
end

function MTH.timefunctionsCMUL(typstring,M,N,K,L,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local NFFT = M*N
	local A = ffi.new(ctyp,2*NFFT) 
	local B = ffi.new(ctyp,2*NFFT)

	-- specific example A
	-- A has dimensions 2x the kernel and C 
	-- generateA(A, M, N, K+2, L+2) 
	-- printMatrix(A,M,N)

	-- specific examples B
	local sum = 1
	for i=0,NFFT do
		A[2*i] = 0--sum % 9 + 1
		A[2*i + 1] = 0
		sum = sum + 1
	end

	-- specific examples A
	sum = 1
	for i=0,NFFT-1 do
		B[2*i] = 0--sum % 9 + 1
		B[2*i + 1] = 0
		sum = sum + 1
	end

	-- printMatrix(B,K,L)

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
		local C = ffi.new(ctyp,2*NFFT)
		for j = 0, NFFT - 1 do 
			C[2*j] = 0
			C[2*j + 1] = 0
		end	
		Cs[i] = C
	end

	-- compute 
	local results = {}
	local caltime = 0
	local check = false
	for i,fn in ipairs(fns) do
		local C = Cs[i]
		local tocall = function() fn(M,N,K,L,A,B,C) end
		-- tocall()
		-- CalcTime(tocall) -- execution time
		printMatrix(A,M,2*N)
		printMatrix(B,M,2*N)
		printMatrix(C,M,2*N)
		check = checkCmul(A,B,C,M,2*N)
		caltime = CalcTime(tocall)
		results[i] = M*N*K*L*2.0*1e-9 / caltime -- gflop
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
	return caltime,check,results
end

function MTH.timefunctionsFFT(typstring,M,N,K,L,NF,...)
	-- EXAMPLE power two treatment and kernel expansion: zero-padding used
	local ctyp = typstring .. "[?] __attribute__((aligned(64)))"
	if M < 2 or N < 2 then print("Image must be >= 2x2\n") return false,{} end

    local NFFTx = math.pow(2,math.ceil(math.log(M)/math.log(2.0)))
	local NFFTy = math.pow(2,math.ceil(math.log(N)/math.log(2.0)))
	-- print("NFFTx : " .. NFFTx)
	-- print("NFFTy : " .. NFFTy)
	local NFFT = NFFTx * NFFTy

	-- A is the data, B the kernel
	local A = ffi.new(ctyp,2*NFFT)
	local B = ffi.new(ctyp,2*NFFT*NF)

	local ker = {}
	local image = {}
	for i=1,NFFT do
		image[i] = i % 9 + 1
		ker[i] = i % 9 + 1
	end

--[[
	local image = {}
	local ker = {}
	for i=1,NFFT,2 do
		image[i] = 5
		ker[i] = 5
		image[i+1] = 4
		ker[i+1] = 4
	end
]]
--[[ -- TEST  4
	local image = {5,4,5,4,
	 		   	   3,2,3,2,
	 		       5,4,5,4,
	 		   	   3,2,3,2}
	
	local ker = {5,4,5,4,
			   3,2,3,2,
			   5,4,5,4,
			   3,2,3,2}

	--expected output 
	--216 0 208 0 216 0 208 0
	--184 0 176 0 184 0 176 0
	--216 0 208 0 216 0 208 0
	--184 0 176 0 184 0 176 0
]]

	-- io.write("Image and kernel expanded to size: ".. NFFTx.." "..NFFTy.."\n")
	-- generating user A
	newcmatrix(A,image,NFFT,true)
	newcmatrix(B,ker,NFFT,true)
	--[[
	printComplexMatrix("Image",NFFTx,NFFTy,A)
	printComplexMatrix("Kernel",NFFTx,NFFTy,B)
	]]	
	local fns = {...}
	local Cs = {}
	-- initialize: fill the matrix C with -1
	for i,fn in ipairs(fns) do
		local C = ffi.new(ctyp,2*NFFT*NF)
		for j = 0, NFFT - 1, 2 do
			C[2*j] = 0
			C[2*j + 1] = 0
		end
		Cs[i] = C
	end

	local results = {}
	local caltime = -1
	for i,fn in ipairs(fns) do	
		local C = Cs[i]
		local tocall = function() fn(M,N,K,L,A,B,C,NF) end
		tocall()
		-- printComplexMatrix("Output",NFFTx,NFFTy,C)
		caltime = CalcTime(tocall)
		results[i] = 1
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
	return caltime,true,results
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
