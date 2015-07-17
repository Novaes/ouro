### Auto-tuner that builds optimized kernels for CNNs ###        
    The different methods are each one on its branch (wip/direct, wip/lowering, wip/fft)
    I'm most currently working on wip/fft

    requirements:
        Terra (github.com/zdevito/terra)
        Accelerate framework for dgemm/convolution tests
    
    running:
        image test: use makefile
        numerical tests: terra src/convolution.lua (make sure terra is in your $PATH or you have an alias to it)

    most important files: 
        src/kernel.t: gerates a convolution of size NB by NB
        src/convolve.lua: generates a convolution over an image
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
      v2.0: vector optimization