local cstdlib = terralib.includec("stdlib.h")
local cstdio = terralib.includec("stdio.h")
local number = double

struct DirL1Package{
    NB: int
    l1conv : {&number, &number, &number, int, int, int, int, int, number} -> {}
    l1convb : {&number, &number, &number, int, int, int, int, int, number, int64, int64} -> {}
    M : int
    N : int
    A : &number
    sda: int
    lda : int
    B : &number
    ldb : int
    C : &number
    sdc : int
    ldc : int
    m : &int
    n : &int
    curr : int
    taskspth : int
    cx : int
    cy : int
}

terra DirL1Package:init(NB: int, M : int, N : int, A : &number, sda: int, lda : int, B : &number, ldb : int, C : &number,
        sdc : int, ldc : int, 
        l1conv : {&number, &number, &number, int, int, int, int, int, number} -> {}, 
        l1convb : {&number, &number, &number, int, int, int, int, int, number, int64, int64} -> {}, 
        taskspth : int, cx : int, cy : int)

    self.l1conv = l1conv
    self.l1convb = l1convb
    self.m = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.n = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.NB = NB
    self.M = M
    self.N = N
    self.A = A
    self.sda = sda
    self.lda = lda
    self.B = B
    self.ldb = ldb
    self.C = C
    self.sdc = sdc
    self.ldc = ldc
    self.curr = 0
    self.taskspth = taskspth
    self.cx = cx
    self.cy = cy
end

terra min(a : int, b : int)
    return terralib.select(a < b, a, b)
end

terra DirL1Package:addblock(m: int, n: int)
    -- cstdio.printf("CURRENT: %d | MAX: %d\n",self.curr,self.taskspth)
    if self.curr >= self.taskspth then
        cstdio.printf("Trying to insert (%d,%d) task in a full thread\n",m,n)
    end

    self.m[self.curr] = m
    self.n[self.curr] = n
    self.curr = self.curr + 1

    if self.curr == self.taskspth then
        return true
    end

    -- not full, over the maximum or empty
    return false 
end

terra dirl1MTComputation(args: &opaque) : &opaque
    var f : &DirL1Package = [&DirL1Package](args)
    -- check received args problem
    var NB : int = (@f).NB
    var l1conv : {&number,&number,&number, int, int, int, int, int, number} -> {} = (@f).l1conv
    var l1convb : {&number,&number,&number, int, int, int, int, int, number, int64, int64} -> {} = (@f).l1convb
    var M : int = (@f).M
    var N : int = (@f).N
    var A : &number = (@f).A
    var sda: int = (@f).sda
    var lda : int = (@f).lda
    var B : &number = (@f).B
    var ldb : int = (@f).ldb
    var C : &number = (@f).C
    var sdc : int = (@f).sdc
    var ldc : int = (@f).ldc
    var m : &int = (@f).m
    var n : &int = (@f).n
    var tasks : int = (@f).curr
    var cx : int = (@f).cx
    var cy : int = (@f).cy
    --cstdio.printf("NB: %d NB2: %d m: %d n: %d k: %d l: %d sda: %d lda: %d ldb: %d ldc: %d kCenterX: %d kCenterY: %d\n",NB,NB2,M,N,K,L,sda,lda,ldb,ldc,kCenterX,kCenterY)
    -- cstdio.printf("M: %d\n",M)

    --compute l1sized kernel
    for i=0,tasks do
        var MM,NN = min(M-m[i],NB),min(N-n[i],NB)
        var isboundary = MM < NB or NN < NB
        var AA,CC = A + (m[i]-cx)*lda + n[i]-cy,C + (m[i]-cx)*ldc + (n[i]-cy)
        if isboundary then -- do not enter here YET
            l1convb(AA,
             B,
             CC,
             sda,lda,ldb,sdc,ldc,0,MM,NN)
        else
            l1conv(AA,
                    B,
                    CC,
                    sda,lda,ldb,sdc,ldc,0) -- -- todo: analyze prefetch argument, past => terralib.select(k == 0,0,1) 
        end
        
    end

    return nil
end