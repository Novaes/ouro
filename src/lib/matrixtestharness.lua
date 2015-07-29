local ffi = require("ffi")
if ffi.os == "Windows" then
	return
end
local adjust = 0

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

local MTH = {}

function assertmm(CR,C,A,B,M,N,K)
	for m=0, M - 1 do
		for n=0, N - 1 do
			CR[m*N + n] = 0
			for k=0, K - 1 do
				CR[m*N + n] = CR[m*N + n] + A[m*K + k] * B[k*N + n] 
			end
			-- if CR[m*N + n] ~= C[m*N + n] then return false end
		end
	end
    return true
end

-- find a standard convolution and make it here
-- 	void cblas_dgemm(int, int,
--                  int, const int M, const int N,
--                  const int K, const double alpha, const double *A,
--                  const int lda, const double *B, const int ldb,
--                  const double beta, double *C, const int ldc);
-- ]]
-- local function referenceconv(M,N,K,L,A,B,C)
-- 	--acc.cblas_dgemm(101,111,111,M,N,K,1.0,A,K,B,N,0.0,C,N)
-- end
function asserteq(C,CR,rows,cols,depth)
	local dim = rows * cols
    for d=0,depth-1 do
        local base = dim * d
        for i=0,rows-1 do
            for j=0,cols-1 do
        		if(C[base + i*cols + j] ~= CR[base + i*cols + j]) then
					return false
				end
            end
        end
    end
    return true
end


function printMatrix(m,rows,columns,depth)
    local dim = rows * columns
    if depth > 1 then 

    	-- print each output
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

function generateTestSet(A,Me,Ne,cx,cy,Bs,K,L,NF)
	for i = cx, Me - cx - 1 do
		for j = cy, Ne - cy - 1 do
			A[i*Ne + j] = math.random(1,9)
		end
	end

	for pos=0,NF -1 do
		local filter = {}	
		for i=1, K*L do
			filter[i] = math.random(1,9)
		end		

		createFilter(filter,K,L,Bs,pos)
	end

	-- An option for a three kernels sequence
	-- createFilter({1,2,1,
	--               0,0,0,
	--               -1,-2,-1},K,L,Bs,0) -- edge detection (Sobel)

	-- createFilter({0,-1,0,
	--               -1,5,-1,
	--               0,-1,0},K,L,Bs,1) -- sharpen

	-- createFilter({0,0,0,
	--               0,1,0,
	--               0,0,0},K,L,Bs,2) -- indentity
end

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

function naiveConvolve(out, inp, Me, Ne, kernel, K, L,depth)
  local kCenterX, kCenterY = math.floor(K/2), math.floor(L/2)
  local dimKer = K*L
  local e = 0
    for d=0,depth-1 do
        local baseKer = dimKer * d
        for i= kCenterX, Me-kCenterX -1 do -- added border to compare with my result
			for j=kCenterY, Ne-kCenterY -1 do
				out[e] = 0
			  	for m=0,K-1 do
			  		for n=0,L-1 do
				      	--boundaries
					    local ii = i + (m - kCenterY)
					    local jj = j + (n - kCenterX)
					    if ii >= 0 and ii< Me and jj>=0 and jj<Ne then
					    	local tmp = out[e]
					    	out[e] = tmp + inp[ii*Ne + jj] * kernel[baseKer + m*L + n]
					    end
			    	end
				end
				e = e + 1
			end
		end
	end
end

function MTH.timefunctions(typstring,M,N,K,L,depth,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local cx,cy = math.floor(K/2), math.floor(L/2)
	local Me, Ne = M+2*cx, N+2*cy
	local A = ffi.new(ctyp,Me*Ne)
	local Bs = ffi.new(ctyp,K*L*depth)
	local CR = ffi.new(ctyp,M*N*depth)

	generateTestSet(A,Me,Ne,cx,cy,Bs,K,L,depth)

	-- Setup GEMM
	local Mlow, Nlow = M*N, K*L
	local Klow, Llow = K*L, depth
	-- result will be Mlow * Llow due to the matrix multiplication
	local AA = ffi.new(ctyp,Mlow*Nlow)
	local BB = ffi.new(ctyp,Klow*Llow)
	local CC = ffi.new(ctyp,Mlow*Llow)
	
	local fns = {...}
	local Cfns = {}

	-- initialize: fill the matrix C with -1
	for i,fn in ipairs(fns) do
		local Cs = ffi.new(ctyp,M*N*depth)
		for j = 0, M * N * depth - 1 do 
			Cs[j] = 0
		end	
		Cfns[i] = Cs
	end

	-- compute 
	local results = {}
	local checked = false
	for i,fn in ipairs(fns) do
		local Cs = Cfns[i] -- 3D

		local tocall = function() fn(Me,Ne,K,L,M,N,A,Bs,Cs,depth,AA,BB,CC) end
		-- tocall()
		results[i] = M*N*K*L*depth*2.0*1e-9 / CalcTime(tocall) -- gflop
		
		-- Print in case detailed analysis
		-- print("Image:")
		-- printMatrix(A,Me,Ne,0)
		-- print("Kernels:")
		-- printMatrix(Bs,K,L,depth)
		-- print("Outputs")
		-- printMatrix(Cs,M,N,depth)
		
		-- Check correctness to any of the function tested
		-- In this case I'm testing only the convolution

		
		

		-- CHECK CORRECTNESS
		naiveConvolve(CR,A,Me,Ne,Bs,K,L,depth)
		checked = asserteq(Cs,CR,M,N,depth)
		-- if checked == false then break end
		-- print("Naive")
		-- printMatrix(CR,M,N,depth)

		-- Print in case detailed analysis
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
	return checked,results
end

function MTH.timefunctionsGEMM(typstring,M,K,N,...)
    local ctyp = typstring.."[?] __attribute__((aligned(64)))"
    local A,B = ffi.new(ctyp,M*K), ffi.new(ctyp,K*N)

    print("M:" .. M .. " K:" .. K .. " N:" .. N)

    for m = 0, M-1 do
        for k = 0, K-1 do
            A[m*K + k] = math.random(0,9)
        end
    end
    for k = 0, K-1 do
        for n = 0, N-1 do
            B[k*N + n] = math.random(0,9)
        end
    end
    local fns = {...}
    local Cs = {}
    for i,fn in ipairs(fns) do
        local C = ffi.new(ctyp,M*N)
        for j = 0, M * N - 1 do
            C[j] = -1
        end 
        Cs[i] = C
    end

    local results = {}
    local check = false
    for i,fn in ipairs(fns) do
        local C = Cs[i]
        local tocall = function() fn(M,K,N,A,B,C) end
        tocall()
        results[i] = M*N*K*2.0*1e-9 / CalcTime(tocall)
        
        local CR = ffi.new(ctyp,M*N)
        check = assertmm(CR,C,A,B,M,N,K)
        -- printMatrix(A,M,K,0)
        -- printMatrix(B,K,N,0)
        -- printMatrix(C,M,N,0)
        -- printMatrix(CR,M,N,0)
        
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
    return check,results
end

return MTH