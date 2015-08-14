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

function asserteq(C,CR,rows,cols)
	for i = 0, rows - 1 do
		for j = 0, cols - 1 do
			if(C[i*cols + j] ~= CR[i*cols + j]) then
				return false
			end
		end
	end
	return true
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

  	local e = 0
	for i= kCenterX, M-kCenterX -1 do -- added border to compare with my result
		for j=kCenterY, N-kCenterY -1 do
			out[e] = 0
		  	for m=0,K-1 do
		  		for n=0,L-1 do
			      	--boundaries
				    local ii = i + (m - kCenterY)
				    local jj = j + (n - kCenterX)
				    -- if ii >= 0 and ii< M and jj>=0 and jj<N then
				    	local tmp = out[e]
				    	out[e] = tmp + inp[ii*N + jj] * kernel[m*L + n] --kernel[mm+1][nn+1];
				    -- end
		    	end
			end
			e = e + 1
		end
  	end
end

function MTH.timefunctions(typstring,M,N,K,L,...)
	local ctyp = typstring.."[?] __attribute__((aligned(64)))"
	local cx,cy = math.floor(K/2), math.floor(L/2)
	local A = ffi.new(ctyp,M*N)
	local Me, Ne = M+2*cx, N+2*cy
	local XX = ffi.new(ctyp,Me*Ne)
	local B = ffi.new(ctyp,K*L)
	local CR = ffi.new(ctyp,M*N)

	-- Generating image not expanded
	for i = 0, M - 1 do
		for j = 0, N - 1 do
			A[i*N + j] = math.random(1,9)
		end
	end

	-- Expanding image
	local count = 0
	for i = cx, Me - cx - 1 do 
		for j = cy, Ne - cy - 1 do
			XX[i*Ne + j] = A[count]
			count = count + 1
		end
	end


	-- specific examples B
	if K == 3 then
		B[0] = 1; B[1] = 2; B[2] = 1;
		B[3] = 0; B[4] = 0; B[5] = 0;
		B[6] = -1;B[7] =-2; B[8] = -1;
	else -- randomizer B
		for k = 0, K-1 do
			for l = 0, L-1 do
				B[k*L + l] = math.random(0,9)
			end
		end
	end
	
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
	local checked = false--false
	for i,fn in ipairs(fns) do
		local C = Cs[i]
		local tocall = function() fn(Me,Ne,K,L,M,N,XX,B,C) end
		tocall()
		results[i] = M*N*K*L*2.0*1e-9 / CalcTime(tocall) -- gflop
		
		-- Check correctness to any of the function tested
		-- In this case I'm testing only the convolution
		naiveConvolve(CR,XX,Me,Ne,B,K,L)
		checked = asserteq(C,CR,M,N)

		-- Print in case detailed analysis
		-- printMatrix(A,M,N) 
		-- printMatrix(XX,Me,Ne)
		-- printMatrix(B,K,L)
		-- printMatrix(C,M,N)
		-- printMatrix(CR,M,N)

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

return MTH
