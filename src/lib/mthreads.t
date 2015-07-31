local cstdlib = terralib.includec("stdlib.h")
local cstdio = terralib.includec("stdio.h")
local number = double

struct L1Package{
    NB: int
    l1dgemm0 : {&number, &number, &number, int, int, int} -> {}
    l1dgemm0b : {&number, &number, &number, int, int, int, number, number, number} -> {}
    l1dgemm1 : {&number, &number, &number, int, int, int} -> {}
    l1dgemm1b : {&number, &number, &number, int, int, int, number, number, number} -> {}
    M : int
    N : int
    K : int
    m : &int
    n : &int
    k : &int
    alpha : int
    A : &number
    lda : int
    B : &number
    ldb : int
    C : &number
    ldc : int

    curr : int
    taskspth : int
}

terra L1Package:init(NB: int, M : int, N : int, K : int, alpha : number, A : &number, lda : int, B : &number, ldb : int, C : &number,
        ldc : int, 
        l1dgemm0 : {&number, &number, &number, int, int, int} -> {},
        l1dgemm0b : {&number, &number, &number, int, int, int, number, number, number} -> {},
        l1dgemm1 : {&number, &number, &number, int, int, int} -> {},
        l1dgemm1b : {&number, &number, &number, int, int, int, number, number, number} -> {},
        taskspth : int)
    self.NB = NB
    self.l1dgemm0 = l1dgemm0
    self.l1dgemm0b = l1dgemm0b
    self.l1dgemm1 = l1dgemm1
    self.l1dgemm1b = l1dgemm1b
    self.m = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.n = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.k = [&int](cstdlib.malloc(taskspth * sizeof(int)))
    self.M = M
    self.N = N
    self.K = K
    self.alpha = alpha
    self.A = A
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

terra L1Package:addblock(m: int, n: int, k: int)
    -- cstdio.printf("CURRENT: %d | MAX: %d\n",self.curr,self.taskspth)
    if self.curr >= self.taskspth then
        cstdio.printf("Trying to insert (%d,%d) task in a full thread\n",m,n)
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
    -- check received args problem
    var NB : int = (@f).NB
    var l1dgemm0 : {&number, &number, &number, int, int, int} -> {}
    var l1dgemm0b : {&number, &number, &number, int, int, int, number, number, number} -> {} 
    var l1dgemm1 : {&number, &number, &number, int, int, int} -> {} 
    var l1dgemm1b : {&number, &number, &number, int, int, int, number, number, number} -> {} 
    var M : int = (@f).M
    var N : int = (@f).N
    var K : int = (@f).K
    var A : &number = (@f).A
    var lda : int = (@f).lda
    var B : &number = (@f).B
    var ldb : int = (@f).ldb
    var C : &number = (@f).C
    var ldc : int = (@f).ldc
    var m : &int = (@f).m
    var n : &int = (@f).n
    var k : &int = (@f).k
    var tasks : int = (@f).curr
    -- cstdio.printf("NB: %d NB2: %d m: %d n: %d k: %d l: %d lda: %d ldb: %d ldc: %d kCenterX: %d kCenterY: %d\n",NB,NB2,M,N,K,L,lda,ldb,ldc,kCenterX,kCenterY)
    -- cstdio.printf("M: %d\n",M)

    --compute l1sized kernel
    for i=0,tasks do
        var MM,NN,KK = min(M-m[i],NB),min(N-n[i],NB),min(K-k[i],NB)
        var isboundary = MM < NB or NN < NB or KK < NB
        var AA,BB,CC = A + (m[i]*lda + k[i]), B + (k[i]*ldb + n[i]), C + (m[i]*ldc + n[i])
        if k == 0 then
            if isboundary then
                --IO.printf("b0 %d %d %d\n",MM,NN,KK)
                l1dgemm0b(AA,BB,CC,lda,ldb,ldc,MM,NN,KK)

                --IO.printf("be %d %d %d\n",MM,NN,KK)
            else
                l1dgemm0(AA,BB,CC,lda,ldb,ldc)
            end
        else
            if isboundary then

                --IO.printf("b %d %d %d\n",MM,NN,KK)
                l1dgemm1b(AA,BB,CC,lda,ldb,ldc,MM,NN,KK)

                --IO.printf("be %d %d %d\n",MM,NN,KK)
            else
                l1dgemm1(AA,BB,CC,lda,ldb,ldc)
            end
        end
    end

    return nil
end