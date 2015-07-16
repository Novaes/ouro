### Auto-tuner that builds optimized kernels for CNNs ###
    
    The different methods are each one on its branch (wip/direct,wip/lowering,wip/fft)
    I'm most currently working on wip/fft

    most important files: 
        src/kernel.t: gerates a convolution of size NB by NB
        src/convolve.lua: generates a convolution over an image
        references/convolution.cpp: future benchmarking with other implementations
        src/examples/: some minimal code of implemented features

    libs:  
        image.t:image library
        matrixtestharness.t: time measure in the auto-tuning process
        multithreads.t: multi-thread library based on pthreads

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: different version or with some specific feature
    
    tags:
      v1.0: no pointer optimization
      v2.0: vector optimization