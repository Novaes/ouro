require "direct"
require "lowering"

local dotune = true
local number = double

local blocksizes = {16,32,64}--{16,32,64,128}
local regblocks = {1,2,4}

-- local vectors = {1,2,4,8,16}
local vectors = {1}
local nthread = {12,24,48}
local filters = {3,5,7,9}--{25,19,25,41}
local nfilter = {1,10,20,30,40,50,70}
-- initialized (defined structure of best) 

-- todo remove  jobs from kernel init
local best = { gflops = 0, k = 3, f = 3, b = 5, rm = 5, rn = 5, v = 1, t = 3}
local NB2 = 5

if dotune then
	local tunefor = 128 -- full size of the matrix
	--change for 10 later
	local harness = require("lib/dir-matrixtestharness")
	for _,k in ipairs(filters) do
		for _,f in ipairs(nfilter) do
			for _,t in ipairs(nthread) do
				for _,b in ipairs(blocksizes) do
					for _,rm in ipairs(regblocks) do
						for _,rn in ipairs(regblocks) do
							for _,v in ipairs(vectors) do
									local i = math.floor(tunefor / b) * b
									local my_conv = genDirconv(b,NB2,rm,rn,v,k,k,f,t)
									-- local my_conv = gennaiveconv()
									-- local my_conv = maxreuse()
									if my_conv then
										print(k,f,t,b,rm,rn)
										my_conv:compile()
										-- bellow line makes do not need boundary cases (image multiple of blocksize)
										local ctyp
										local correct, exectimes = harness.timefunctions(tostring(number),i,i,k,k,f, function(Me,Ne,K,L,M,N,A,B,C,f)
											-- to gennaive pass the #kernels here
					                    	my_conv(nil,Me,Ne,K,L,1.0,A,Me,Ne,B,L,C,M,N,K/2,L/2) -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
										end)
										if not correct then	print("<error>") break end
										print(i,unpack (exectimes),"[OK]")
										local curr_gflops = exectimes[1]
										-- print(curr_gflops) -- print analysis 
										if best.gflops < curr_gflops then --  Maximization problem (the greater gflops, the better)
											best = { gflops = curr_gflops, b = b, rm = rm, rn = rn, v = v, k = k, f = f, t = t }
											--terralib.tree.printraw(best)
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







