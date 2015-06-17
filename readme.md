### Auto-tuner that builds optimized kernels for CNNs ###
    
    
    most important files: 
        src/kernel.t: gerates a convolution of size NB by NB
        src/convolve.lua: generates a convolution over an image
        references/convolution.cpp: future benchmarking with other implementations


    libs:  
        image.t:image library
        matrixtestharness.t: time measure in the auto-tuning process

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: different version or with some specific feature
    
    tags:
      v1.0: no pointer optimization