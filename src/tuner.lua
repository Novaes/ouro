cstdio = terralib.includec("stdio.h")
require "direct"
require "lowering"
require "fftbased"

-- General Parameters
local dotune = true
local filters = {3}
local nfilter = {20}
local nthread = {4}
local tunefor = 32
local number = double

-- FLAGS
local silent = false
local trackbest = false
local nodirect = false
local nolowering = false
local nofft = false

-- handle FLAGS 
local argc = 0
repeat 
	str = arg[argc]
	argc = argc + 1
until arg[argc] == nil


for i=0, argc-1 do
	local str = arg[i]
	if str == "--silent" then
		silent = true  
	elseif str == "--trackbest" then
		trackbest = true
	elseif str == "--nodirect" then
		nodirect = true
	elseif str == "--nolowering" then
		nolowering = true
	elseif str == "--nofft" then
		nofft = true
	elseif str == "--help" then
		io.write("execution: terra [terra-options] <source-file> [arguments]\n")
		io.write("	--silent: hide all time parameter tested results \n")
		io.write("	--trackbest: see every time you have a new best kernel, show its parameters\n")
		io.write("	--license: see license details\n")
		io.write("	--nodirect: do not compute direct\n")
		io.write("	--nolowering: do not compute lowering\n")
		io.write("	--nofft: do not compute FFT based\n")
		return
	elseif str == "--license" then
		local f = io.open("../LICENSE.txt", "rb")
    	local license = f:read("*all")
    	f:close()
		io.write(license)
		print("\n")
		return
	end
end

local method

local best = { time = 10, filter_size = 3, filters = 3, block_size = 5, rm = 5, rn = 5, vector_size = 1, threads = 3, method = "direct"}

-- Direct auto-tuner
local blocksizes = {16,32}
local regblocks = {1,2,4}
local vectors = {1}
-- todo remove  jobs from kernel init
local NBF = 5

if dotune and not nodirect then
	print("\nApplying Direct method ")
	local harness = require("lib/dir-matrixtestharness")
	for _,k in ipairs(filters) do
		for _,f in ipairs(nfilter) do
			for _,t in ipairs(nthread) do
				for _,b in ipairs(blocksizes) do
					for _,rm in ipairs(regblocks) do
						for _,rn in ipairs(regblocks) do
							for _,v in ipairs(vectors) do
									local i = math.floor(tunefor / b) * b
									local my_conv = genDirconv(b,NBF,rm,rn,v,k,k,f,t)
									-- local my_conv = gennaiveconv()
									-- local my_conv = maxreuse()
									if my_conv then
										my_conv:compile()
										-- bellow line makes do not need boundary cases (image multiple of blocksize)
										local curr_time, correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k,f, function(Me,Ne,K,L,M,N,A,B,C,f)
											-- to gennaive pass the #kernels here
					                    	my_conv(nil,Me,Ne,K,L,1.0,A,Me,Ne,B,L,C,M,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
										end)
										if not correct then	print("<error>") break end
										if not silent then print(b,rm,rn,v,k,f,curr_time, "[OK]") end
										if best.time > curr_time then --  Maximization problem (the greater gflops, the better)
											best = { time = curr_time, block_size = b, rm = rm, rn = rn, vector_size = v, filter_size = k, filters = f, threads = t, method = "direct" }
											if trackbest then 
												terralib.tree.printraw(best)
											end
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end


-- Lowering auto-tuner
local blocksizes = {16,20,32}
local regblocksM = {1}--{1,2,4}
local regblocksN = {2}--{1,2,4}
local vectors = {1,2,4,8}
local NBF = 5
local bl = 8 -- lowering blocksize

if dotune and not nolowering then
	print("\nApplying Lowering method ")
	local harness = require("lib/low-matrixtestharness")
	for _,t in ipairs(nthread) do
		for _,f in ipairs(nfilter) do
			for _,k in ipairs(filters) do
				for _,b in ipairs(blocksizes) do
					for _,rm in ipairs(regblocksM) do
						for _,rn in ipairs(regblocksN) do
							for _,v in ipairs(vectors) do				
									-- local my_conv = gennaiveconv()
								local my_conv = genLoweringConv(bl,b,NBF,rm,rn,v,t)
								-- local my_conv = generatedgemm(b,NBF,rm,rn,v)
								if my_conv then
								--	print(b,rm,rn,v,k,f)
									my_conv:compile()
									
									-- bellow line makes do not need boundary cases (image multiple of blocksize)
									local i = math.floor(tunefor / b) * b
									local ctyp
									local curr_time,correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k,f, 
										function(Me,Ne,K,L,M,N,A,Bs,Cs,f,AA,BB,CC)
											-- to gennaive pass the #kernels here
				                    		my_conv(nil,A,Me,Ne,K,L,1.0,Bs,L,Cs,M,N,K/2,L/2,f,AA,BB,CC) 
				                    		-- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
									end)
									-- test only GEMM
									-- local correct, exectimes = harness.timefunctionsGEMM(tostring(number),i,i,i,function(M,K,N,A,B,C)
									-- 	my_conv(nil,M,N,K,1.0,A,K,B,N,C,N)
									-- end)
									if not correct then	print("<error>") break end
									if not silent then print(b,rm,rn,v,k,f,curr_time, "[OK]") end
									if best.time > curr_time then --  Maximization problem (the greater gflops, the better)
										best = { time = curr_time, block_size = b, rm = rm, rn = rn, vector_size = v, filter_size = k, filters = f, threads = t, method = "lowering" }
										if trackbest then
											terralib.tree.printraw(best)
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
end


local blocksizes = {8,16,32}
local regblocks = {1,2,4} -- blocksizes must be divisible by RN*V
local vectors = {1--[[1,2,4,8]]}
local NBF = 1

-- kernel dependent 
if dotune and not nofft then
	print("\nApplying FFT based method ")
	local harness = require("lib/fft-matrixtestharness")
	for _,f in ipairs(nfilter) do
		for _,b in ipairs(blocksizes) do
			for _,rm in ipairs(regblocks) do
				for _,rn in ipairs(regblocks) do
					for _,v in ipairs(vectors) do
						for _,t in ipairs(nthread) do
							-- FOR FAST CONVOLUTION THE KERNEL IS DEPENDENT OF THE IMAGE SIZE
							local i = math.floor(tunefor / b) * b
							
							local my_fastconv = genfastconvMT(b,1,rm,rn,v,tunefor,t)
							-- local my_fastconv = genNumRecipesConv()
							-- local my_fastconv = genfastconv(b,1,rm,rn,v,tunefor)
							if my_fastconv then
								my_fastconv:compile()
								local curr_time, correct, exectimes = harness.timefunctionsFFT(tostring(number),i,i,3,3,f, function(M,N,K,L,A,B,C,NF) 
			                    	my_fastconv(nil,M,N,A,B,C,N,NF)
								end)
								-- if not correct then	("<error>")  break  end
								if not silent then print(b,rm,rn,v,t,curr_time) end
								if best.time > curr_time then --  Maximization problem (the greater gflops, the better)
									best = {time = curr_time, block_size = b, rm = rm, rn = rn, vector_size = v, filters = f, threads = t, method = "FFT based" }
									if trackbest then
										terralib.tree.printraw(best)
									end
								end
							end
						end
					end
				end
			end
		end
	end
end


print("\nBest method found: ")
terralib.tree.printraw(best)

print("\nSaving it in the bin/ folder... ")

local my_numconv

if best.method == "direct" then
	my_numconv = genDirconv(best.block_size,NBF,best.rm,best.rn,best.vector_size,best.filter_size,best.filter_size,best.filters,best.threads)
elseif best.method == "lowering" then
	my_numconv = genLoweringConv(bl,best.block_size,NBF,best.rm,best.rn,best.vector_size,best.threads)
elseif best.method == "FFT based" then
	my_numconv = genfastconvMT(best.block_size,1,best.rm,best.rn,best.vector_size,tunefor,best.threads)
end

if number == double then
	terralib.saveobj("../bin/my_dconv.o", {my_numconv = my_numconv})
else
	terralib.saveobj("../bin/my_sconv.o", {my_numconv = my_numconv})
end

print("\nSucessfully saved. ")

