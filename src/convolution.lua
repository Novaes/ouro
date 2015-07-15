local cstdio = terralib.includec("stdio.h")
local MT = terralib.includec("pthread.h")
local stdlib = terralib.includec("stdlib.h")
local number = double


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

-- generate L1 convolution 
function genkernel(NB, RM, RN, V, prefetch, K, L, boundary)
    local M,N--[[, boundaryargs
    -- if one of the parameters is lower than NB then receive the usual
    if boundary then 
        M,N = symbol(int64,"M"),symbol(int64,"N")
        boundaryargs = terralib.newlist({M,N})
    else
        boundaryargs = terralib.newlist()
        M,N = NB,NB
    end
]]
    -- assets NB/RM and NB/RN do not necessary
    -- print("parameters: "..NB .." ".. RM .." ".. RN .." ".. V .." ".. K .." ".. L)

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
    local sda,lda,ldb,ldc = symbol("sda"),symbol("lda"),symbol("ldb"), symbol("ldc")
    local a,b,c = symmat("a",RM+2,RN+2), symmat("b",K,L), symmat("c",RM,RN)
    local kk, ll = symbol("kk"), symbol("ll")
    local x,y = symbol("x"), symbol("y")
    local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

    for m = 0, RM+1 do
        for n = 0, RN+1 do
            loadA:insert(quote
                    var [a[m][n]] = vecload(A, m*ldc + n*V)
            end)
            if(m>=0 and m<RM and n>=0 and n<RN) then
                loadc:insert(quote
                        var [c[m][n]] = alpha * vecload(C, (m+1)*ldc + (n+1)*V)
                end)
                storec:insert(quote
                    vecstore(C, (m+1)*ldc + (n+1), [c[m][n]])
                end)
            end
        end
    end

    local calcc = terralib.newlist()

    -- load full kernel
    for  k=0, K-1 do
        for l = 0, L-1 do
            loadkernel:insert(quote
                var [b[k][l]] = vecload(B, k*ldb + l*V)
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
                    --no boundary cases
                    calcc:insert(
                        quote
                            -- area regblocking not multiple of the area sizeblocking
                            if([mm] + m < NB-1 and [nn] + n < NB-1) then 
                                --remeber that taking the pos a[x+1][y+1], e.g. a[0][0] menas take a[-1][-1] necessary for c[0][0]
                                [c[m][n]] = [c[m][n]] + [a[x+1][y+1]] * [b[k][l]]
                            end
                        end
                    )
                end
            end
        end
    end

    return terra([A] : &double, [B] : &double, [C] : &double, [sda] : int, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double--[[, [boundaryargs] ]])
        -- no borders, original from 0 to NB-1 (it is in TERRA, exclusive loop)
        -- If the kernel is different from 3x3, started indices and pointers updates will change (it can be generalized)
        for [mm] = 1, NB-1, RM do
            -- how it goes by blocking, it can be greater than NB-1
            -- the correct for blocking would be use min([nn]+RN*V,NB-1), 
            -- however the generation of the code could not be done first, unless many ifs would be inserted  
            for [nn]=1, NB-1, RN*V do 
                [loadc];
                [loadkernel];
                llvmprefetch(A + sda*lda,0,3,1);
                [loadA];
                [calcc];
                [storec];
                A = A + RN*V
                C = C + RN*V
            end
            -- adjust when regblocking is not multiple of blocksize
            if (((NB-2)/(RN*V)) * (RN*V) + 1 < NB-1) then
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

local thsize = 2
local taskspth = 4/thsize -- tunefor/thsize

struct L1Package{
    NB: int
    -- array of function pointers
    l1conv : {&double, &double, &double, int, int, int, int, double} -> {}
    M : int
    N : int
    A : &double
    sda: int
    lda : int
    B : &double
    ldb : int
    C : &double
    ldc : int
    m : int[taskspth]
    n : int[taskspth]
    curr : int
}

local pkg = L1Package

terra L1Package:init(NB: int, M : int, N : int, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
        ldc : int, l1conv : {&double, &double, &double, int, int, int, int, double} -> {})
	self.l1conv = l1conv
    self.NB = NB
    self.M = M
    self.N = N
    self.A = A
    self.sda = sda
    self.lda = lda
    self.B = B
    self.ldb = ldb
    self.C = C
    self.ldc = ldc
    self.curr = 0
end

terra L1Package:addblock(m: int, n: int)
	if self.curr >= taskspth then
		cstdio.printf("Trying to insert (%d,%d) task in a full thread\n",m,n)
	end
	self.m[self.curr] = m
	self.n[self.curr] = n
	self.curr = self.curr + 1
end

terra l1MTComputation(args: &opaque) : &opaque
    var f : &L1Package = [&L1Package](args)
    -- check received args problem
    var NB : int = (@f).NB
    var l1conv : {&double,&double,&double, int, int, int, int, double} -> {} = (@f).l1conv
    var M : int = (@f).M
    var N : int = (@f).N
    var A : &double = (@f).A
    var sda: int = (@f).sda
    var lda : int = (@f).lda
    var B : &double = (@f).B
    var ldb : int = (@f).ldb
    var C : &double = (@f).C
    var ldc : int = (@f).ldc
    var m : &int = (@f).m
    var n : &int = (@f).n
    var tasks : int = (@f).curr
    -- cstdio.printf("NB: %d NB2: %d m: %d n: %d k: %d l: %d sda: %d lda: %d ldb: %d ldc: %d kCenterX: %d kCenterY: %d\n",NB,NB2,M,N,K,L,sda,lda,ldb,ldc,kCenterX,kCenterY)
    -- cstdio.printf("M: %d\n",M)

    --compute l1sized kernel
    for i=0,tasks do
	    var MM,NN = min(M-m[i],NB),min(N-n[i],NB)
	    var isboundary = MM < NB or NN < NB
	    var AA,CC = A + (m[i]*lda + n[i]),C + (m[i]*ldc + n[i])
	                
	    l1conv(AA,
	    B,
	    CC,
	    sda,lda,ldb,ldc,0)
	end

    return nil
end

function genconvolution(NB,NBF,RM,RN,V)
    -- register blocking does not need to be a a multiple of the blocksize anymore
    -- if not isinteger(NB/(RN*V)) or not isinteger(NB/RM) then
    -- return false
    -- end

    --5 times NB minimum by dgemm
    --local NB2 = NBF * NB

    local NB2 = NB * NBF

    -- no prefetch, no boundary
    local l1conv0 = genkernel(NB, RM, RN, 1, false, 3, 3, false)
    
    -- local l1conv0b = genkernel(NB, RM, RN, 1, false, 3, 3, true)
    return terra(gettime : {} -> double, M : int, N : int, K : int, L: int, 
        alpha : double, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double, 
        ldc : int, kCenterX: int, kCenterY: int) 

        var threads : MT.pthread_t[thsize]
        var count = 0
        var pkgs : L1Package[thsize]
        for i=0, thsize do
			pkgs[i]:init(NB, M, N, A, sda, lda, B, ldb, C, ldc, l1conv0)
		end
        -- cstdio.printf("NB: %d NB2: %d m: %d n: %d k: %d l: %d sda: %d lda: %d ldb: %d ldc: %d kCenterX: 
        -- 	%d kCenterY: %d\n",pkg.NB,pkg.NB2,pkg.M,pkg.N,pkg.K,pkg.L,pkg.sda,pkg.lda,
        -- 	pkg.ldb,pkg.ldc,pkg.kCenterX,pkg.kCenterY)
        [ blockedloop(M,N,{NB2,NB},
            function(m,n)
                return quote
                	pkgs[count/thsize]:addblock(m,n)
                	cstdio.printf("adding to thread: %d\n",count/thsize)
                	if (count+1) % thsize == 0 then
                		cstdio.printf("---> thread launched: %d\n",count/thsize)
                    	if MT.pthread_create(&threads[count/thsize], nil, l1MTComputation , &pkgs[count/thsize]) ~= 0 then 
                    		cstdio.printf("Thread #%u creation error",threads[count/thsize])
                    	end
                    end
                    count = count + 1
       		 end end) ]
        for i=0,thsize do
        	if MT.pthread_join(threads[i],nil) ~= 0 then
        		cstdio.printf("Thread #%u join error",i) 
        	end
        end
    end
end

local blocksizes = {5--[[,10 ,16,24,32,40,48,56,64,1024]]}
local regblocks = {1--[[, 2,3]]}
local vectors = {1 --[[,2,4,8,16]]}

-- initialized (defined structure of best)
local best = { gflops = 0, b = 5, rm = 5, rn = 5, v = 1 }

if dotune then
    local tunefor = 10 --[[1024]] -- full size of the matrix
    --change for 10 later
    local harness = require("lib/matrixtestharness")
    for _,b in ipairs(blocksizes) do
        for _,rm in ipairs(regblocks) do
            for _,rn in ipairs(regblocks) do
                for _,v in ipairs(vectors) do
                        -- same until here
                    local my_conv = genconvolution(b,1,rm,rn,v)
                    if my_conv then
                        print(b,rm,rn,v)
                        my_conv:compile()
                        local i = math.floor(tunefor / b) * b
                        local curr_gflops = 0
                        local ctyp
                        local correct, exectimes = harness.timefunctions(tostring(number),i,i,3,3, function(M,N,K,L,A,B,C)
                            -- my_conv receives integer parameter i.e. it represents floor of K/2 and L/2
                            my_conv(nil,M,N,K,L,1.0,A,M,N,B,L,C,N,K/2,L/2) 
                        end)
                        if not correct then print("<error>")  break  end
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

local my_convolution = genconvolution(best.b,1,best.rm,best.rn,best.v)
terralib.saveobj("my_conv.o", {my_convolution = my_convolution})