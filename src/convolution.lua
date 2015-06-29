local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")
local number = double
local alignment = 8

local function isinteger(x) return math.floor(x) == x end
local llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})

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

local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end

local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)

-- generate L1 convolution 
function genkernel(NB, RM, RN, V, prefetch, K, L, boundary)
	local M,N, boundaryargs
	-- if one of the parameters is lower than NB then receive the usual
	if boundary then 
		M,N = symbol(int64,"M"),symbol(int64,"N")
		boundaryargs = terralib.newlist({M,N})
	else
		boundaryargs = terralib.newlist()
		M,N = NB,NB
	end

	-- assets NB/RM and NB/RN do not necessary
	-- print("parameters: "..NB .." ".. RM .." ".. RN .." ".. V .." ".. K .." ".. L)

	local A,B,C,mm,nn,alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
	local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2,RN+2), symmat("b",K,L/V), symmat("c",RM,RN)
	local kk, ll = symbol("kk"), symbol("ll")
	local x,y = symbol("x"), symbol("y")
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()
	local calcc = terralib.newlist()

	local VP = &vector(double,V)
	local VimgP = &vector(double,V) -- scale it

    for m = 0, RM+1 do
        for n = 0, RN-1 do -- todo remove RN+1
            loadA:insert(quote
                var [a[m][n]] : vector(double,V) = unalignedload(VimgP(&A[m*ldc + n*V])) -- it is not n*V
            end)
        end
    end

	for m = 0, RM-1 do
		for n = 0, RN-1 do
            loadc:insert(quote
                var [c[m][n]] = alpha * C[(m+1)*ldc + (n+1)]
            end)
            storec:insert(quote
                C[(m+1)*ldc + (n+1)] = [c[m][n]]
            end)
		end
	end

	-- load full kernel
	for  k=0, K-1 do
		for l = 0, L/V - 1 do
			loadkernel:insert(quote
				var [b[k][l]] : vector(double,V) = unalignedload(VP(&B[k*ldb + l*V]))
			end)
		end
	end

	-- spatial 2D convolution
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			for k=0, K-1 do
				for l = 0, L/V-1 do
					-- would sum mm or nn, but this position is realtive to this mini-block (rm, rn)
					x, y = m + (k - math.floor(K/2) ), n + (l - math.floor(L/2))
					--no boundary cases
					calcc:insert(
						quote
							-- area regblocking not multiple of the area sizeblocking
							if([mm] + m < NB-1 and [nn] + n < NB-1) then
								--remeber that taking the pos a[x+1][y+1], e.g. a[0][0] means take a[-1][-1] necessary for c[0][0]
								-- because for each block (0,0) means (1,1) for example

								-- do a function that takes this return &{vector(double,3)} -> {} and set the 
								var v : vector(double,V) = [a[x+1][y+1]] * [b[k][l]]
								var sum = 0
								for i=0,V do
									sum = sum + v[i]
								end
								[c[m][n]] = [c[m][n]] + sum

								-- IO.printf("VALUE %f ",v[0])
								-- extractsum(&([vector(double,V)]([a[x+1][y+1]]) * [b[k][l]]),V)
								-- [c[m][n]] = [c[m][n]] + extractsum([r[m][n]])
							end
						end
					)
				end
			end
		end
	end

	return terra([A] : &double, [B] : &double, [C] : &double, --[[[T] : &double]] [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double, [boundaryargs])
		-- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
		-- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
		for [mm] = 1, NB-1, RM do
			-- how it goes by blocking, it can be greater than NB-1
			-- the correct for blocking would be use min([nn]+RN*V,NB-1), 
			-- however the generation of the code could not be done first, unless many ifs would be inserted  
			for [nn]=1, NB-1, RN*V do 
				-- IO.printf("loading C...\n")
				[loadc];
				-- IO.printf("load C done\n")
                -- IO.printf("loading B...\n")
				[loadkernel];
				-- IO.printf("load B done\n")
				-- IO.printf("loading A...\n")
			    -- llvmprefetch(A + sda*lda,0,3,1);
				[loadA];
				-- IO.printf("loading A done\n")
				-- IO.printf("calculating C\n")
				[calcc];
				-- IO.printf("calculating C done\n")
				-- IO.printf("storing C back\n")
				[storec];
				A = A + RN*V
				C = C + RN*V
			end
			if ( ((NB-2)/(RN*V)) * RN*V + 1 < NB-1) then
				var offset = (((NB-2)/(RN*V)) * (RN*V) + 1) + (RN*V)  - (NB-1)
				A = A - offset
				C = C - offset
			end
			-- jump of two (final border one line, initial border next line)
			-- It is two because the kernel is 3, it would change for different kernel
			C = C + 2
			A = A + 2

			A = A + RM * ldc - NB
			C = C + RM * ldc - NB
		end
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

function genconvolution(NB,NBF,RM,RN,V,K,L)
	-- register blocking does not need to be a a multiple of the blocksize anymore
    -- needed for vector instruction
	-- if not isinteger((NB-2)/(RN*V)) or not isinteger((NB-2)/RM) then
	-- 	return false
	-- end

	--5 times NB minimum by dgemm
	--local NB2 = NBF * NB

	local NB2 = NB * NBF

	-- no prefetch, no boundary
	local l1conv0 = genkernel(NB, RM, RN, V, false, K, L, false)
	local l1conv0b = genkernel(NB, RM, RN, V, false, K, L, true)


	return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
		alpha : double, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
		ldc : int, kCenterX: int, kCenterY: int) 

         [ blockedloop(N,M,{NB2,NB},
                function(m,n) 
                return quote
                    var MM,NN = min(M-m,NB),min(N-n,NB)
                    var isboundary = MM < NB or NN < NB
                    var AA,CC = A + (m*lda + n),C + (m*ldc + n)
                    -- if isboundary then -- do not enter here YET
                    --  l1conv0b(AA,
                    --      B,
                    --      CC,
                    --      sda,lda,ldb,ldc,0,MM,NN)
                    -- else
	                    l1conv0(AA,
	                     B,
	                     CC,
	                     sda,lda,ldb,ldc,0) -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
                    -- end
                end end)  
            ]       
	end
end

-- Different blocksizes for the same result implies in padding overheading 
-- for small blocks
local blocksizes = {5,--[[10,16,24,32,40,48,56,64,1024]]}
local regblocksM = {1}
local regblocksN = {1}
-- local vectors = {1,2,4,8,16}
local vectors = {3}

-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

if dotune then
	-- local tunefor = 1024
	local tunefor = 10 -- full size of the matrix
	--change for 10 later
	local harness = require("lib/matrixtestharness")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocksM) do
			for _,rn in ipairs(regblocksN) do
				for _,v in ipairs(vectors) do
						-- same until here
					local my_conv = genconvolution(b,1,rm,rn,v,3,3)
					if my_conv then
						print(b,rm,rn,v)
						my_conv:compile()
						local i = math.floor(tunefor / b) * b
						local curr_gflops = 0
						local ctyp
						local correct, exectimes = harness.timefunctions(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C)
	                    	my_conv(nil,M,N,K,L,1.0,A,M,N,B,L,C,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
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
end

local my_convolution = genconvolution(best.b,1,best.rm,best.rn,best.v,3,3)
terralib.saveobj("my_conv.o", {my_convolution = my_convolution})