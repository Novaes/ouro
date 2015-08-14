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

local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end

local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)


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

function gencmulkernel(NB, RM, RN, prefetch)
	local M,N
	local A,B,C,mm,nn = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn")
	local ldc = symbol("ldc")
	local a,b,c,caddr = symmat("a",RM,RN), symmat("b",RM,RN), symmat("c",RM,RN), symmat("caddr",RM,RN)
	local load,storec,calcc = terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM-1 do
		for n = 0, RN - 1 do
            load:insert(quote
            	var [c[m][n]] = C[m*ldc + n]
                var [a[m][n]] = A[m*ldc + n]
                var [b[m][n]] = B[m*ldc + n]
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
				end
			)
		end
	end

	-- NB must be divisible by RN
	return terra([A] : &double, [B] : &double, [C] : &double, [ldc] : int) 
		for [mm] = 0, NB, RM do
			for [nn]=0, NB, RN do
				[load];
				-- llvmprefetch(A + ldc*ldc,0,3,1);
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
	local l1cmul = gencmulkernel(NB, RM, 2*RN, false)

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, 
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

local blocksizes = {16,32,64}
local regblocks = {1,2,4} -- blocksizes must be divisible by RN
local vectors = {1}

-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

if dotune then
	-- local tunefor = 1024
	local tunefor = 1024 -- full size of the matrix
	--change for 10 later
	local harness = require("cmult-matrixtestharness")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocks) do
			for _,rn in ipairs(regblocks) do
				for _,v in ipairs(vectors) do
						-- same until here
					local my_cmul = gencmul(b,1,rm,rn)
					if my_cmul then
						print(b,rm,rn,v)
						my_cmul:compile()
						local i = math.floor(tunefor / b) * b
						local curr_gflops = 0
						local ctyp
						local correct, exectimes = harness.timefunctionsCMUL(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C)
	                    	my_cmul(nil,M,N,A,B,C,N)
						end)
						if not correct then	print("<error>")  break  end
						print(i,unpack (exectimes),"[OK]")
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

-- local my_convolution = genconvolution(best.b,1,best.rm,best.rn,best.v)
-- terralib.saveobj("my_conv.o", {my_convolution = my_convolution})
