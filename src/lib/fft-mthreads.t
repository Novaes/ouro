cstdlib = terralib.includec("stdlib.h")
cstdio = terralib.includec("stdio.h")

struct FFTL1Package{
    NFFT: int
    X : &double,
    Y : &double,
    Z : &double,
    l1fft : {&double} -> {}
    bases : &int
    curr : int
    taskspth : int
    type : int
}

terra FFTL1Package:multikernels(k : int)
    self.taskspth = self.taskspth * k
    cstdlib.free(self.bases)
    self.bases = [&int](cstdlib.malloc(self.taskspth * sizeof(int)))
    self.curr = 0
    self.type = 1
end

terra FFTL1Package:usefilter()
    cstdlib.free(self.bases)
    self.bases = [&int](cstdlib.malloc(self.taskspth * sizeof(int)))
    self.curr = 0
    self.type = 1
end

terra FFTL1Package:useoutput()
    cstdlib.free(self.bases)
    self.bases = [&int](cstdlib.malloc(self.taskspth * sizeof(int)))
    self.curr = 0
    self.type = 2
end

terra FFTL1Package:clear()
    cstdlib.free(self.bases)
    self.bases = [&int](cstdlib.malloc(self.taskspth * sizeof(int)))
    self.curr = 0
end

terra FFTL1Package:init(
            NFFT : int,
            X : &double,
            Y : &double,
            Z : &double,
            l1fft : {&double} -> {}, 
            taskspth : int,
            type : int)
    
    self.NFFT = NFFT
    self.X = X
    self.Y = Y
    self.Z = Z
    self.l1fft  = l1fft
    self.taskspth = taskspth
    self.bases = [&int](cstdlib.malloc(self.taskspth * sizeof(int)))
    self.curr = 0
    self.type = type
end

terra min(a : int64, b : int64)
    return terralib.select(a < b, a, b)
end

terra FFTL1Package:addblock(b : int)
    if self.curr >= self.taskspth then
        cstdio.printf("Trying to insert (%d) task in a full thread\n",b)
    end
    self.bases[self.curr] = b
    self.curr = self.curr + 1

    if self.curr == self.taskspth then
        return true
    end

    -- not full, over the maximum or empty
    return false 
end

terra fftl1MTComputation(args: &opaque) : &opaque
    var f : &FFTL1Package = [&FFTL1Package](args)
    var NFFT : int = (@f).NFFT
    var X : &double = (@f).X
    var Y : &double = (@f).Y
    var Z : &double = (@f).Z
    var l1fft : {&double} -> {} = (@f).l1fft
    var bases : &int = (@f).bases
    var tasks : int = (@f).curr
    var type : int = (@f).type
    -- cstdio.printf("NFFT: %d tasks %d\n",NFFT,tasks)
    for i=0,tasks do
        var index = bases[i]
        if type == 0 then 
            l1fft(&X[index])
        elseif type == 1 then
            l1fft(&Y[index])
        elseif type == 2 then
            l1fft(&Z[index])
        end
    end
    return nil
end