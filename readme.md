### Auto-tuner that builds optimized kernels for Convolution Layers on CNNs ###        

    The three different methods to compute the convolution are each one in their branches: Direct Method (dev/direct-mt-mk), Lowering (dev/lowering-mt) and Convolution by FFT (wip/fastconvolution)

    requirements:
        Terra (github.com/zdevito/terra)

    running:
        image test: use makefile (It generates imageconv.o)
        numerical tests: terra src/convolution.lua (it generates dconv.o or sconv.o (Single/Double Precision) )
        *Make sure terra is in your $PATH or you have an alias to it

    most important branches: 
        -> dev/master: it has all three methods merged.
        -> dev/direct-mt-mk: by Direct method (multi-threaded and over multiple-kernels)
        -> dev/lowering-mt: by Lowering (multi-threaded, using optimized GEMM)
        -> wip/fastconvolution: by FFT method (using kernels: FFTKERNELS, TRANSPOSE and CMULT)
        dev/direct-vec: Direct method only using vecinstr
        dev/direct-mt: 2D (Spatial) direct convolution
        
    most important files: 
        src/tuner.lua: auto-tuner. It test the three methods. Use tuner.lua --help to see options.
        benchmarks/convolution.c: main benchmarks with other implementations
        src/imageconv.lua: generates an RGB image convolution
        src/examples/: some minimal code of implemented features

    libs:  
        image.t:image library (adapted to the project needs from github.com/jameshegarty/darkroom)
        <method>-matrixtestharness.t: time measure in the auto-tuning process. 
        <method>-mthreads: multi-thread library based on pthreads
        fft-kernels.t: n-point kernels for FFT

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: different version or with some specific feature
    
    tags:
      v1.0: no pointer optimization
      v1.5: vector optimization
      v2.0: standard direct