local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")

local number = double

local function isinteger(x) return math.floor(x) == x end

local dotune = true

-- naive convolution
-- terra naivel1matmul(A : &double, B : &double, C : &double, lda : int, ldb : int, ldc : int, alpha : double)

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
function genkernel(NB, RM, RN, V, prefetch, K, L)
	-- assert: no boundary cases
	-- assert(isinteger(NB / RN)) -- number of iterations over b for matrix multiplication
	-- assert(isinteger(NB / RM)) -- number of iterations over a for matrix multiplication
	
	-- print("parameters: "..NB .. " " .. RM .." ".. RN .. " " .. V .." ".. K .." ".. L)
	
	local M,N = NB,NB -- new
	local VP = &vector(double,V)
	local terra vecload(data : &double, idx : int)
		var addr = &data[idx]
		return @VP(addr)
	end

	local terra vecstore(data : &double, idx : int, v : vector(double,V))
		var addr = &data[idx]
		@VP(addr) = v
	end

	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local lda,ldb,ldc = symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2,RN+2), symmat("b",K,L), symmat("c",RM,RN) -- changed
	local kk, ll = symbol("kk"), symbol("ll") -- new
	local x,y = symbol("x"), symbol("y") --new
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM+1 do 
		for n = 0, RN+1 do
			loadA:insert(quote
				-- mm+1 and nn+1 because I start from -1 now 
				var [a[m][n]] = vecload(A,( (mm-1) + m)*ldc + (nn-1) + n) 
				end)
			if(m>=0 and m<RM and n>=0 and n<RN) then
				loadc:insert(quote
					var [c[m][n]] = alpha * vecload(C,(mm+m)*ldc + nn + n)
					end)
				storec:insert(quote
					vecstore(C,(mm+m)*ldc + nn + n,[c[m][n]])
					end)
			end
		end
	end

	local calcc = terralib.newlist()

	-- load full kernel
	for  k=0, K-1 do
		for l = 0, L-1 do
			loadkernel:insert(quote
				var [b[k][l]] = vecload(B, k*ldb + l)
			end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			for k=0, K-1 do
				for l = 0, L-1 do
					-- would sum mm or nn, but this position is realtive to this mini-block (rm, rn)
					x, y = m + (k - math.floor(K/2) ), n + (l - math.floor(L/2))
					--if x >= -1 and x<RM+1 and y>=-1 and y<RN+1 then --RM+1 and RN+1 are the image loads, no boundaries cases
					calcc:insert(quote
						--remeber that taking the pos a[x+1][y+1], e.g. a[0][0] menas take a[-1][-1] necessary for c[0][0]
						[c[m][n]] = [c[m][n]] + [a[x+1][y+1]] * [b[k][l]]
					end)
				end
			end
		end
	end

	-- optimization point
	return terra([A] : &double, [B] : &double, [C] : &double, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double)
		-- no borders, original from 0 to NB-1
		for [mm] = 1, NB-2, RM do
			for [nn] = 1, NB-2, RN do
				[loadA];
				[loadc];
				[loadkernel];
				[calcc];
				[storec];
			end
		end
	end
end

terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

function genconvolution(NB,NBF,RM,RN,V)
	if not isinteger(NB/RN) or not isinteger(NB/RM) then -- NB/(RN*V) when vectorize
		return false
	end

	--5 times NB minimum by dgemm
	--local NB2 = NBF * NB
	local NB2 = NB * NBF

	-- EXAMPLES
	--local l1matmul = genkernel(NB, 3, 2, 1, false, 3, 3) 
	local l1matmul = genkernel(NB, 3, 3, 1, false, 3, 3)

	return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
		alpha : double, A : &double, lda : int, B : &double, ldb : int, C : &double, 
		ldc : int, kCenterX: int, kCenterY: int) 
		-- use blocking on this loop
		for mm = 0,M,NB2 do
			for nn = 0,N,NB2 do
				for m = mm,min(mm+NB2,M),NB do
					for n = nn,min(nn+NB2,N),NB do
						l1matmul(A + m*lda + n,
						         B, -- B, fixed kernel
						         C + m*ldc + n,
						         lda,ldb,ldc,0)
								-- about this last prefetch argument (instead of 0): terralib.select(k == 0,0,1) 
								-- if it is the first positon, 
								-- no because it will be in cache already, but if it is not. 
								-- Do it, because this prefetch will but the next line already 
								-- in the cache
					end
				end
			end
		end
	end
end

local blocksizes = {16,24,32,40,48,56,64,1024}
-- local blocksizes = {5}
-- local blocksizes = {1024}
local regblocks = {2,4,5,1}
-- local regblocks = {5}
-- local regblocks = {1}
-- local vectors = {1,2,4,8,16}
local vectors = {1}
-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

if dotune then
	local tunefor = 1024
	-- local tunefor = 5
	local harness = require("lib/matrixtestharness")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocks) do
			for _,rn in ipairs(regblocks) do
				for _,v in ipairs(vectors) do
					-- same until here
					local my_conv = genconvolution(b,5,rm,rn,v)
					if my_conv then
						print(b,rm,rn,v)
						my_conv:compile()
						local i = math.floor(tunefor / b) * b
						local curr_gflops = 0
						local ctyp
						local correct, exectimes = harness.timefunctions(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C)   
                        		my_conv(nil,M,N,K,L,1.0,A,N,B,L,C,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
						end)
						if not correct then	print("<error>")  break  end
						print(i,unpack (exectimes))
						local curr_gflops = exectimes[1]
						-- print(curr_gflops) -- print analysis 
						if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
							best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v }
							terralib.tree.printraw(best)
						end
					end
				end
			end
		end
	end
	terralib.tree.printraw(best)
end

local my_convolution = genconvolution(best.b,5,best.rm,best.rn,best.v)
terralib.saveobj("my_conv.o", {my_convolution = my_convolution})