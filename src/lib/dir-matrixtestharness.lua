local ffi = require("ffi")
if ffi.os == "Windows" then
	return
end
local adjust = 0

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

local MTH = {}


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
	-- Open file; append mode
	for i = cx, Me - cx - 1 do
		for j = cy, Ne - cy - 1 do
			A[i*Ne + j] = math.random(1,9)
		end
	end

	-- appends a word test to the last line of the file

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

function newfiles(M,N,K,NF,pad)
	local inpfile, filterfile, outfile
	local prefix = "data/"
	local param = "k"..K.."nf"..NF.."imgx"..M.."img"..N.."p"..pad
	
	local inpname = prefix.."IMAGE-"..param..".txt"
	local filtername = prefix.."FILTERS-"..param..".txt"
	local outname = prefix.."OUTPUT-"..param..".txt"

	inpfile = io.open(inpname, "a")
	filterfile = io.open(filtername, "a")
	outfile = io.open(outname, "a")
	return inpfile, filterfile, outfile
end

function MTH.timefunctions(typstring,M,N,K,L,depth,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local cx,cy = math.floor(K/2), math.floor(L/2)
	local Me, Ne = M+2*cx, N+2*cy
	local A = ffi.new(ctyp,Me*Ne)
	local Bs = ffi.new(ctyp,K*L*depth)
	local CR = ffi.new(ctyp,M*N*depth)

	-- only in the end it will go over all inputs and outputs save in a file
	-- if you will use it often; modify it to re-use input generation
	local savefile = false
	

	generateTestSet(A,Me,Ne,cx,cy,Bs,K,L,depth)


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
	local time
	for i,fn in ipairs(fns) do
		local Cs = Cfns[i] -- 3D
		local tocall = function() fn(Me,Ne,K,L,M,N,A,Bs,Cs,depth) end

		time = CalcTime(tocall)
		results[i] = M*N*K*L*depth*2.0*1e-9 /  time-- gflop
		
		-- Print in case detailed analysis
		-- print("Image:")
		-- printMatrix(A,Me,Ne,0)
		-- print("Kernels:")
		-- printMatrix(Bs,K,L,depth)
		-- print("Outputs")
		-- printMatrix(Cs,M,N,depth)
		
		-- Check correctness to any of the function tested
		-- In this case I'm testing only the convolution

		-- print("Naive")
		-- printMatrix(CR,M,N,depth)

		-- CHECK CORRECTNESS
		naiveConvolve(CR,A,Me,Ne,Bs,K,L,depth)
		checked = asserteq(Cs,CR,M,N,depth)
		if checked == false then break end

		if savefile then 
			local inpfile, filterfile, outfile = newfiles(M,N,K,depth,cx)

			for i=0,M*N-1 do inpfile:write(A[i]) end
			for i=0,K*L*depth-1 do filterfile:write(Bs[i]) end
			for i=0,M*N*depth-1 do outfile:write(CR[i]) end

			inpfile:close() 
			filterfile:close()
			outfile:close()
		end

		
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
	return time,checked,results
end

return MTH