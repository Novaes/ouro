### Auto-tuner that builds optimized kernels for CNNs ###        

    The three different methods to compute the convolution are each one in their branches: Direct Method (master), Lowering (wip/optlowering) and Convolution by FFT (wip/fastconvolution)

    requirements:
        Terra (github.com/zdevito/terra)
        Accelerate framework for dgemm/convolution tests
    
    running:
        image test: use makefile (It generates imageconv.o)
        numerical tests: terra src/convolution.lua (it generates numconv.o)
        *Make sure terra is in your $PATH or you have an alias to it

    most important branches: 
        -> master: by Direct method (multi-threaded and using vector instr)
        -> wip/optlowering: by Lowering (multi-threaded, using optimized GEMM)
        -> wip/fastconvolution: by FFT method (using kernels: FFTKERNELS, TRANSPOSE and CMULT)
        wip/tuneNumOfKernels: Direct method auto-tuning the number of kernels
        wip/mthreading: Direct method only multi-threading
        wip/vectinstr: Direct method only using vecinstr
        wip/benchmarks: Benchmark with MKL for Direct and FFT
        
    most important files: 
        src/convolution.lua: generates the numerical convolution over an image
        references/convolution.cpp: future benchmarking with other implementations
        src/examples/: some minimal code of implemented features

    libs:  
        image.t:image library (adapted to the project needs from github.com/jameshegarty/darkroom)
        matrixtestharness.t: time measure in the auto-tuning process
        multithreads.t: multi-thread library based on pthreads
        fftkernels.t: n-point kernels for FFT

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: different version or with some specific feature
    
    tags:
      v1.0: no pointer optimization
      v1.5: vector optimization
      v2.0: standard direct