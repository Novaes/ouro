cstdlib = terralib.includec("stdlib.h")
cstdio = terralib.includec("stdio.h")

struct L1Package{
    NB: int
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
    m : &int
    n : &int
    curr : int
    taskspth : int
}

terra L1Package:init(NB: int, M : int, N : int, A : &double, sda: int, lda : int, B : &double, ldb : int, C : &double,
        ldc : int, l1conv : {&double, &double, &double, int, int, int, int, double} -> {}, taskspth : int)
    self.l1conv = l1conv
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
    self.ldc = ldc
    self.curr = 0
    self.taskspth = taskspth
end

terra min(a : int, b : int)
    return terralib.select(a < b, a, b)
end

terra L1Package:addblock(m: int, n: int)
    if self.curr >= self.taskspth then
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