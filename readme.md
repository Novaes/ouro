### Auto-tuner that builds optimized kernels for CNNs ###        

    The three different methods to compute the convolution: Direct Method, Lowering and Convolution by FFT 

    requirements:
        Terra (github.com/zdevito/terra)
    
    running:
        image test: use makefile (It generates imageconv.o)
        numerical tests: terra src/convolution.lua (it generates numconv.o)
        *Make sure terra is in your $PATH or you have an alias to it

    most important branches: 
        -> Direct Method (multi-threaded and auto-tuned)
        -> Lowering (multi-threaded, using optimized GEMM)
        -> FFT-based method (using kernels: FFTKERNELS, TRANSPOSE and CMULT)
        
    most important files: 
        src/tuner.lua: auto-tuner that iterates over the three methods 
        benchmarks/: benchmarking with other implementations
        src/examples/: some minimal code of implemented features

    libs:  
        image.t:image library (adapted to the project needs from github.com/jameshegarty/darkroom)
        <method>-matrixtestharness.t: time measure in the auto-tuning process
        <method>-multithreads.t: multi-thread library based on pthreads
        fft-kernels.t: n-point kernels for FFT

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: different version or with some specific feature
    
    tags:
      v1.0