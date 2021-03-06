-- It is the standard GEMM from Terra, but Multi-threaded
require "lib/low-mthreads"
local stdlib = terralib.includec("stdlib.h")
local cmath = terralib.includec("math.h")
local cstdio = terralib.includec("stdio.h")
local MT = terralib.includec("pthread.h")

local number = double
local alignment = 8
local dotune = false
-- local thsize = 1
-- area tunefor/area blocks / thsize
-- local taskspth = 4

function symmat(name,I,...)
	if not I then return symbol(name) end
	local r = {}
	for i = 0,I-1 do
		r[i] = symmat(name..tostring(i),...)
	end
	return r
end

local function isinteger(x) return math.floor(x) == x end

llvmprefetch = terralib.intrinsic("llvm.prefetch",{&opaque,int,int,int} -> {})
local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end
local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)

function genkernel(NB, RM, RN, V, alpha, boundary)
	local M,N,K, boundaryargs

	-- if one of the parameters is lower than NB then receive the usual
	if boundary then 
		M,N,K = symbol(int64,"M"),symbol(int64,"N"),symbol(int64,"K")
		boundaryargs = terralib.newlist({M,N,K})
	else
		boundaryargs = terralib.newlist()
		M,N,K = NB,NB,NB
	end

	local A,B,C,mm,nn,ld = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("ld")
	local lda,ldb,ldc = symbol("lda"),symbol("ldb"),symbol("ldc")
	local a,b,c,caddr = symmat("a",RM), symmat("b",RN), symmat("c",RM,RN), symmat("caddr",RM,RN)
	local k = symbol("k")

	local loadc,storec = terralib.newlist(),terralib.newlist()
	local VT = vector(number,V)
	local VP = &VT
	
	-- it will take V elements at one c position
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			--alpha C mult, load C to fast memory (register vector), multiply it by alpha
			loadc:insert(quote
				var [caddr[m][n]] = C + m*ldc + n*V
				var [c[m][n]] = alpha * unalignedload(VP([caddr[m][n]]))
			end)
			--put back in the slow memory, from register vector to memory
			storec:insert(quote
				unalignedstore(VP([caddr[m][n]]),[c[m][n]])
			end)
		end
	end

	local calcc = terralib.newlist()

	for n = 0, RN-1 do
		calcc:insert(quote
			var [b[n]] = unalignedload(VP(&B[n*V]))
		end)
	end
	for m = 0, RM-1 do
		calcc:insert(quote
			var [a[m]] = VT(A[m*lda])
		end)
	end
	for m = 0, RM-1 do 
		for n = 0, RN-1 do
			calcc:insert(quote
				[c[m][n]] = [c[m][n]] + [a[m]] * [b[n]]
			end)
		end
	end

	local result = terra([A] : &number, [B] : &number, [C] : &number, [lda] : int64,[ldb] : int64,[ldc] : int64,[boundaryargs])
		for [mm] = 0, M, RM do
			for [nn] = 0, N,RN*V do
				[loadc];
				for [k] = 0, K do
					llvmprefetch(B + 4*ldb,0,3,1);
					[calcc];
					B = B + ldb
					A = A + 1
				end
				[storec];
				A = A - K
				B = B - ldb*K + RN*V
				C = C + RN*V
			end
			C = C + RM * ldb - N
			B = B - N
			A = A + lda*RM
		end
	end
	return result
end

local terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

function blockedloop(N,M,K,blocksizes,bodyfn)
  local function generatelevel(n,ii,jj,kk,bb0,bb1,bb2)
    if n > #blocksizes then
      return bodyfn(ii,jj,kk)
    end
    local blocksize = blocksizes[n]
    return quote for i = ii,min(ii+bb0,N),blocksize do
                   for j = jj,min(jj+bb1,M),blocksize do
                      for k = kk,min(kk+bb2,K),blocksize do
                        [ generatelevel(n+1,i,j,k,blocksize,blocksize,blocksize) ]
           end end end end 
    end
  return generatelevel(1,0,0,0,N,M,K) 
end

function generatedgemm(NB,NBF,RM,RN,V,thsize)
	--[[
	if not isinteger(NB/(RN*V)) or not isinteger(NB/RM) then
		return false
	end
]]
	local NB2 = NBF * NB
	local l1dgemm0 = genkernel(NB,RM,RN,V,0,false)
	local l1dgemm1 = genkernel(NB,RM,RN,V,1,false)
	local l1dgemm0b = genkernel(NB,1,1,1,0,true)
	local l1dgemm1b = genkernel(NB,1,1,1,1,true)

	return terra(gettime : {} -> double, M : int64, N : int64, K : int64, alpha : number, A : &number, lda : int64, B : &number, ldb : int64, 
		            C : &number, ldc : int64)
		
		-- If I do M*M/(NB*NB) it gives an overflow
		var thwork : int64 = M/NB
		var thr : int64 = thsize

		-- cstdio.printf("NB: %d\n",NB)
		-- cstdio.printf("M: %d\n",M)
		-- cstdio.printf("THWORK: %d\n",thwork)
		
		-- #threads is greater than the work
		while (cmath.floor(thwork / thr) == 0) do
			thr = thr - 1
		end
		var tmp : number = thwork 
		var ftaskspth : number = tmp / thr
		var taskspth : int = thwork / thr
		
		-- adjust in case of not perfect division
		-- cstdio.printf("ftaskspth: %f\n",ftaskspth)
		-- cstdio.printf("taskspth: %d\n",taskspth)
		var rr : int = 0
		if ftaskspth ~= taskspth then
			rr = thwork%thr
			if rr == 0 then
				rr = (ftaskspth-taskspth)*thr
			end
		end

		var threads : MT.pthread_t[thsize]
        var pkgs : L1Package[thsize]
        var added = 0
    	for i=0, thr do
        	if rr > 0 then 
        		added = 1
				rr = rr - 1
			end
			pkgs[i]:init(NB, M, N, K, A, lda, B, ldb, C, ldc, l1dgemm0b, l1dgemm0, l1dgemm1b, l1dgemm1, taskspth+added)
			added = 0
		end
		var count = 0

		for mm = 0,M,NB2 do
			for nn = 0,N,NB2 do
				for kk = 0,K, NB2 do
					for m = mm,min(mm+NB2,M),NB do
						for n = nn,min(nn+NB2,N),NB do
							for k = kk,min(kk+NB2,K),NB do
								-- cstdio.printf("adding to thread: %d\n",count)
								if pkgs[count]:addblock(m,n,k) then
				                	if MT.pthread_create(&threads[count], nil, l1MTComputation , &pkgs[count]) ~= 0 then 
			                    		-- cstdio.printf("Thread #%u creation error",threads[count])
			                    	end
			                    		-- cstdio.printf("---> thread launched: %d\n",count)
			                    	count = count + 1
			                    end
							end
						end
					end
				end
			end
		end

	    -- cstdio.printf("M %d N %d K %d\n", M,N,K)
	    --[[
		[ blockedloop(N,M,K,{NB2,NB},
				function(m,n,k) return quote
				
		end end) ]
		]]

	   	for i=0,thr do
        	if MT.pthread_join(threads[i],nil) ~= 0 then
        		-- cstdio.printf("Thread #%u join error\n",i) 
        	end
        end

	end
end

local blocksizes = {16,24,32,40,48,56,64}
local regblocks = {1,2,4}
local vectors = {1,2,4,8,16}
--local best = { gflops = 0, b = 56, rm = 4, rn = 1, v = 8 }
local best = { gflops = 0, b = 40, rm = 4, rn = 2, v = 4 }
--local best = { gflops = 0, b = 40, rm = 1, rn = 1, v = 1 }

if dotune then
	local tunefor = 1024
	local harness = require("lib/matrixtestharnessGEMM")
	for _,b in ipairs(blocksizes) do
		for _,rm in ipairs(regblocks) do
			for _,rn in ipairs(regblocks) do
				-- for_,rnn in ipairs(regblocks)
				--
				for _,v in ipairs(vectors) do
					-- same until here
					local my_dgemm = generatedgemm(b,5,rm,rn,v,1)
					if my_dgemm then
						print(b,rm,rn,v)
						my_dgemm:compile()
						local i = math.floor(tunefor / b) * b
						local avg = 0
						local ctyp
						local s, times = harness.timefunctions(tostring(number),i,i,i,function(M,K,N,A,B,C)
							--lower
							my_dgemm(nil,M,N,K,1.0,A,K,B,N,C,N)
							--lifting
						end)
						if not s then
							print("<error>")
							break
						end
						print(i,unpack(times))
						local avg = times[1]	
						if  best.gflops < avg then
							best = { gflops = avg, b = b, rm = rm, rn = rn, v = v }
							terralib.tree.printraw(best)
						end
					end
				end
			end
		end
	end
end
--[[
local my_dgemm = generatedgemm(best.b, 5, best.rm, best.rn, best.v)
if number == double then
	terralib.saveobj("my_dgemm.o", { my_dgemm = my_dgemm })
else
	terralib.saveobj("my_sgemm.o", { my_sgemm = my_dgemm })
end
]]