cstdlib = terralib.includec("stdlib.h")
cstdio = terralib.includec("stdio.h")

struct L1Package{
    NB: int
    M : int
    N : int
    K : int
    A : &double
    lda : int64
    B : &double
    ldb : int64
    C : &double
    ldc : int64
    m : &int
    n : &int
    k : &int
    l1dgemm0b : {&double,&double,&double, int64, int64, int64,  int64, int64, int64} -> {}
    l1dgemm0 : {&double,&double,&double, int64, int64, int64} -> {}
    l1dgemm1b : {&double,&double,&double, int64, int64, int64, int64, int64, int64} -> {}
    l1dgemm1 : {&double,&double,&double, int64, int64, int64} -> {}
    curr : int
    taskspth : int
}

terra L1Package:init(NB: int, M : int, N : int, K : int, A : &double, lda : int64, B : &double, ldb : int64, C : &double,
            ldc : int64,
            l1dgemm0b : {&double,&double,&double, int64, int64, int64,  int64, int64, int64} -> {}, 
            l1dgemm0 : {&double,&double,&double, int64, int64, int64} -> {},
            l1dgemm1b : {&double,&double,&double, int64, int64, int64,  int64, int64, int64} -> {},
            l1dgemm1 : {&double,&double,&double, int64, int64, int64} -> {},
            taskspth : int)

    self.l1dgemm0b = l1dgemm0b
    self.l1dgemm0  = l1dgemm0 
    self.l1dgemm1b  = l1dgemm1b
    self.l1dgemm1  = l1dgemm1 
    self.m = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.n = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.k = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.NB = NB
    self.M = M
    self.N = N
    self.K = K
    self.A = A
    self.lda = lda
    self.B = B
    self.ldb = ldb
    self.C = C
    self.ldc = ldc
    self.curr = 0
    self.taskspth = taskspth
end

terra min(a : int64, b : int64)
    return terralib.select(a < b, a, b)
end

terra L1Package:addblock(m : int, n : int, k : int)
    if self.curr >= self.taskspth then
        cstdio.printf("Trying to insert (%d, %d, %d) task in a full thread\n",m,n,k)
    end
    self.m[self.curr] = m
    self.n[self.curr] = n
    self.k[self.curr] = k
    self.curr = self.curr + 1

    if self.curr == self.taskspth then
        return true
    end

    -- not full, over the maximum or empty
    return false 
end

terra l1MTComputation(args: &opaque) : &opaque
    var f : &L1Package = [&L1Package](args)
    var NB : int = (@f).NB
    var l1dgemm0b : {&double,&double,&double, int64, int64, int64,  int64, int64, int64} -> {} = (@f).l1dgemm0b
    var l1dgemm0 : {&double,&double,&double, int64, int64, int64} -> {} = (@f).l1dgemm0
    var l1dgemm1b : {&double,&double,&double, int64, int64, int64, int64, int64, int64} -> {} = (@f).l1dgemm1b
    var l1dgemm1 : {&double,&double,&double, int64, int64, int64} -> {} = (@f).l1dgemm1
    var M : int = (@f).M
    var N : int = (@f).N
    var K : int = (@f).K
    var A : &double = (@f).A
    var lda : int64 = (@f).lda
    var B : &double = (@f).B
    var ldb : int64 = (@f).ldb
    var C : &double = (@f).C
    var ldc : int64 = (@f).ldc
    var m : &int = (@f).m
    var n : &int = (@f).n
    var k : &int = (@f).k
    var tasks : int = (@f).curr
    -- cstdio.printf("NB: %d NB2: %d m: %d n: %d k: %d l: %d sda: %d lda: %d ldb: %d ldc: %d kCenterX: %d kCenterY: %d\n",NB,NB2,M,N,K,L,sda,lda,ldb,ldc,kCenterX,kCenterY)
    -- cstdio.printf("M: %d\n",M)


    for i=0,tasks do
        var MM,NN,KK = min(M-m[i],NB),min(N-n[i],NB),min(K-k[i],NB)
        var isboundary = MM < NB or NN < NB or KK < NB
        var AA,BB,CC = A + (m[i]*lda + k[i]),B + (k[i]*ldb + n[i]),C + (m[i]*ldc + n[i])
        cstdio.printf("MM %d NN %d KK %d\n",MM,NN,KK)
        if k[i] == 0 then
            if isboundary then
                l1dgemm0b(AA,BB,CC,lda,ldb,ldc,MM,NN,KK)
            else
                l1dgemm0(AA,BB,CC,lda,ldb,ldc)
            end
        else
            if isboundary then
                l1dgemm1b(AA,BB,CC,lda,ldb,ldc,MM,NN,KK)
            else
                l1dgemm1(AA,BB,CC,lda,ldb,ldc)
            end
        end
    end

    return nil
end