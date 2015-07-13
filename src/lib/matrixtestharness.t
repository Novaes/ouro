IO = terralib.includec("stdio.h")
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
--  void cblas_dgemm(int, int,
--                  int, const int M, const int N,
--                  const int K, const double alpha, const double *A,
--                  const int lda, const double *B, const int ldb,
--                  const double beta, double *C, const int ldc);
-- ]]
-- local function referenceconv(M,N,K,L,A,B,C)
--  --acc.cblas_dgemm(101,111,111,M,N,K,1.0,A,K,B,N,0.0,C,N)
-- end

function createFilter(filter,rows,columns,m,pos)
	local base = pos * rows * columns
	for i=0,rows-1 do
		for j=0,rows-1 do
			m[base + i*columns + j] = filter[i*columns + j + 1]   
      -- io.write(m[i*columns + j])
  end
    -- io.write("\n")
end
  -- io.write("\n")
  return m
end

function printMatrix(m,rows,columns,depth)
	local dim = rows * columns
	if depth > 1 then 
		for d=0,depth-1 do
			local base = dim * d
			for i=0,rows-1 do
				for j=0,columns-1 do
					io.write(" " .. m[base + i*columns + j])
				end
				io.write("\n")
			end
			io.write("\n")
		end
	else
		for i=0,rows-1 do
			for j=0,columns-1 do
				io.write(" " .. m[i*columns + j])
			end
			io.write("\n")
		end
	end
end

function getpart(A,sM,eM,sN,eN,M,N,AA,u)
	local iRows, iCols = eM, eN
	local kRows,kCols = 3,3
	local kCenterX, kCenterY = math.floor(kRows/2), math.floor(kCols/2)
	for m = sM, eM do
		for n = sN, eN do
			for k=0, kRows -1 do
				for l=0, kCols-1 do
					local ii = m + (k - kCenterY)
					local jj = n + (l - kCenterX)
					AA[u] = A[ii * N + jj]
					u = u + 1
				end
			end
		end
	end
end

-- create a 1-9 with 0 padding border, used for test convolution
function fillpart(A,sM,eM,sN,eN,M,N,test) 
	local counter = 0
	local dim = eN-sN + 1
	local elems = dim * dim
	for m = sM, eM do 
		for n = sN, eN do
			if test then
				A[m*N + n] = counter + 1
				counter = (counter + 1) % elems
			else
				A[m*N + n] = math.random(0,9)
			end
		end
	end
end

function genRepMatrixA(A, M, N, inX, inY, test, AA)
    local nRows = M/inX - 1 -- first matrix start from 0
    local nCols = N/inY - 1
    for j=0,nRows do
    	for i = 0,nCols do
    		fillpart(A, 1 + j*inY, (inY-2) + j*inY, 1 + i*inX, (inX-2) + i*inX, M, N, test)
    	end
    end
end

function genFilterForLowering(nker,filter,rows,columns,m,pos)
	for i=0,rows*columns - 1 do
		m[i*nker + pos] = filter[i+1]
	end
end


function genImageForLowering(A, K, L, nker, test)
	local count = 0
	for i=0,K*L*nker - 1 do
        if test then -- specific test
        	A[i] = (count + 1) % 10
        	if A[i] == 0  then
        		count = count + 2
        		A[i] = count % 10
        	else
        		count = count + 1
        	end
        else 
        	A[i] = math.random(0,9)
        end
    end
end

function generateTestSet(A,Bs,M,N,K,L,nker,mtype)
	if mtype == "lowering" then
		local sobel, sharpen, identity = {1, 2,1, 0,0, 0,-1,-2,-1},
		{0,-1,0,-1,5,-1, 0,-1, 0},
		{0, 0,0, 0,1, 0, 0, 0, 0}
		genImageForLowering(A, K, L, nker,true)
		genFilterForLowering(nker,sobel,K,L,Bs,0)
		genFilterForLowering(nker,sharpen,K,L,Bs,1)
		genFilterForLowering(nker,identity,K,L,Bs,2)
	else -- kernel
	    genRepMatrixA(A, M, N, K+2, L+2, true) -- A has dimensions 2x the kernel and C 
	    createFilter({1,2,1,
	    	0,0,0,
	                  -1,-2,-1},K,L,Bs,0) -- edge detection (Sobel)
	    createFilter({0,-1,0,
	    	-1,5,-1,
	                  0,-1,0},K,L,Bs,1) -- sharpen
	    createFilter({0,0,0,
	    	0,1,0,
	                  0,0,0},K,L,Bs,2) -- indentity
	end
end


function generateRandomizedSet(A,Bs,M,N,K,L,depth,mtype)

	if mtype == "lowering" then
		genImageForLowering(A, K, L, nker,true)
		genFilterForLowering(nker,sobel,K,L,Bs,0) 
		genFilterForLowering(nker,sharpen,K,L,Bs,1) 
		genFilterForLowering(nker,identity,K,L,Bs,2) 
	else
        -- randomizing A
        genRepMatrixA(A, M, N, K+2, L+2, false)
        -- randomizing Bs
        for i=0,depth do
        	local b = {}
        	for j=1, K*L do
        		b[j] = math.random(0,9)
        	end
        	createFilter(b,K,L,Bs,i) 
        end
    end
end

-- inX and inY are kernel with padding
function lowerImage(AA,A,M,N,lowertype,inX,inY)
	if lowertype == "type 1" then
		local nRows = math.floor(M/inX)
		local nCols = math.floor(N/inY)
		local u = 0
		for j=0,nRows - 1 do
			for i = 0,nCols - 1 do
				getpart(A, 1 + j*inY, (inY-2) + j*inY, 1 + i*inX, (inX-2) + i*inX, M, N, AA, u)
	            u = u + (inX-2)*(inY-2) * (inX-2)*(inY-2) --inside 
	        end
	    end
	end
end
		 
function liftResult(CCs,Cs,Mlow,Llow,lowertype)
	if lowertype == "type 1" then
		local count = 0
		for i=0,Llow - 1 do
			for j=0,Mlow - 1 do

				CCs[count] = Cs[j*Llow + i]
				count = count + 1
			end
		end
	end
end

function lowerKernel(BBs,Bs,K,L,depth,lowertype)
	if lowertype == "type 1" then
		local base = K*L
		local count = 0
		for d=0,depth-1 do
			for i=0,K*L - 1 do
				BBs[i*depth + d] = Bs[base*d + i]
			end
		end
	end
end

function MTH.timefunctionsMM(typstring,M,N,depth,K,L,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
    -- necessary remotion of paddings to the lower matrix (to the lower it will be needed)
    local A = ffi.new(ctyp,M*N)
    local Bs = ffi.new(ctyp,K*L*depth)
    local nker = 3

    generateTestSet(A,Bs,M,N,K,L,depth, "direct")
    printMatrix(A,M,N,1)

    IO.printf("\n")
    printMatrix(Bs,K,L,depth)
    local blocks = M/(K+2) * M/(L+2) 
    local Mlow, Nlow = blocks * K*L, K*L
    local Klow, Llow = K*L, depth
    local AA = ffi.new(ctyp,Mlow*Nlow)
    local BBs = ffi.new(ctyp,Klow*Llow)
    lowerImage(AA,A,M,N,"type 1",K+2,L+2)
    lowerKernel(BBs,Bs,K,L,depth,"type 1")

    local fns = {...}
    local Cfns = {}
    local CCfns = {}

    -- for each function: zeroing its result (only one in use nowadays)
    for i,fn in ipairs(fns) do
    	local Cs = ffi.new(ctyp,Mlow*Llow)
    	local CCs = ffi.new(ctyp,K*L*depth)
    	for j = 0, Mlow * Llow - 1 do Cs[j] = -1 end 
    	for j = 0, M * N * depth - 1 do CCs[j] = -1 end 
    	Cfns[i] = Cs
    	CCfns[i] = CCs
    end

    local results = {}
    for i,fn in ipairs(fns) do
        -- local B = Bs[2]
        local Cs = Cfns[i]
        local CCs = CCfns[i]
        local tocall = function() fn(Mlow,Llow,Nlow,AA,BBs,Cs) end
        tocall()
        CalcTime(tocall) -- execution time
        printMatrix(AA,Mlow,Nlow,1)
        IO.printf("\n")
        printMatrix(BBs,Klow,Llow,1)
        IO.printf("\n")
        printMatrix(Cs,Mlow,Llow,1)
        IO.printf("\n")
        liftResult(CCs,Cs,Mlow,Llow,"type 1")
        printMatrix(CCs,K,L,depth)
        results[i] = CalcTime(tocall)
        results[i] = M*N*K*L*2.0*1e-9 / CalcTime(tocall) -- gflop --Mlow*Llow*Nlow
        if i ~= 1 then
        	local C0 = Cfns[1]
        	local C1 = Cfns[i]
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

function MTH.timefunctions(typstring,M,N,depth,K,L,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local A = ffi.new(ctyp,M*N)
	local Bs = ffi.new(ctyp,K*L*depth)

    -- THE TEST HAS 3 KERNELS
    -- depth = 3 
    generateTestSet(A,Bs,M,N,K,L,depth, "direct")
    -- generateRandomizedSet(A,Bs,M,N,K,L,depth)

    local fns = {...}
    local Cfns = {}

    -- initialize: fill the matrix Cs with -1
    for i,fn in ipairs(fns) do
         local Cs = ffi.new(ctyp,M*N*depth)
         for j = 0, M * N * depth - 1 do 
             Cs[j] = 0
         end 
         Cfns[i] = Cs
    end

    -- compute
    local results = {}
    for i,fn in ipairs(fns) do
        -- local B = Bs[2]
        local Cs = Cfns[i]
        local tocall = function() fn(M,N,K,L,A,Bs,Cs) end
        tocall()
        CalcTime(tocall) -- execution time
        -- printMatrix(A,M,N,1)
        -- printMatrix(Bs,K,L,depth)
        printMatrix(Cs,K,L,depth)
        results[i] = CalcTime(tocall)
        results[i] = M*N*K*L*depth*2.0*1e-9 / CalcTime(tocall) -- gflop
        if i ~= 1 then
        	local C0 = Cfns[1]
        	local C1 = Cfns[i]
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