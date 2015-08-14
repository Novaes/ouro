### Ouro         
    It is an Auto-tuner that builds optimized kernels for ConvNets

    Given convolution parameters, It runs the three main convolution methods: Direct Method, Lowering and FFT-based. Then, it gives the best method .o and auto-tuned parameters. The project is multi-staged. It uses Lua as a high-level language and Terra as the low-level staged code. 

    requirements:
        Terra (github.com/zdevito/terra)
        Terra is a low level language created with inoperability with Lua in mind. 

    running:
        image test: use the makefile (It generates imageconv.o)
        numerical tests: terra src/tuner.lua (it generates dconv.o or sconv.o)
        You can run it as terra src/tuner.lua --help to see how it works or just see the code. 
        *Make sure terra is in your $PATH or you have an alias to it

    most important branches: 
        -> Direct Method 
        -> Lowering (multi-threaded, using optimized GEMM))
        State of the art of the method: http://arxiv.org/abs/1504.04343
        -> FFT-based method (also called: "Fast Convolution") 
        Features: Multi-threaded, using: 2-points and 4-points FFT kernels, kernel for transpose  and kernel for point-wise multiplication (blocked and auto-tuned)
        State of the art of the method: http://arxiv.org/abs/1412.7580  
    
   
    most important files: 
        src/tuner.lua: auto-tuner that iterates over the three methods 
        benchmarks/: benchmarking with other implementations
        src/examples/: some minimal code of implemented features

    libs:  
        image.t:image library (from github.com/jameshegarty/darkroom)
        <method>-matrixtestharness.t: time measure in the auto-tuning process
        <method>-multithreads.t: multi-thread library based on pthreads
        fft-kernels.t: n-point kernels for FFT

    branching tags:
      wip: works in progress
      junk: experiments
      bug: fixing a bug
      dev: intermediate developments (to be put on tests/ folder) or with some specific feature
    
    tags:
      v1.0

    Note: the FFT-based has accuracy errors and the idea is improve it and fix an upper-bound to round-off errors