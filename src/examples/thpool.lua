local IO = terralib.includec("stdio.h")
local MT = terralib.includec("pthread.h")
local POOL = terralib.includec("../include/thpool.h")
terralib.linklibrary("../lib/thpool.so")

terra task1(args: &opaque) : &opaque
    var x = MT.pthread_self()
    var y : &double = [&double](args)
    -- IO.printf("Thread #%u working on task1\n", [int64](x))
    IO.printf("Thread #%u working on task1 with value %f\n", [int64](x), @y)
end

terra main()
    var thpool : POOL.threadpool = POOL.thpool_init(4)
    var pkg : double = 3.0

    for i=0,10 do
    	-- for no arguments: pass nil instead of &pkg
        POOL.thpool_add_work(thpool, task1, &pkg) 
    end
    
    POOL.thpool_wait(thpool)

    POOL.thpool_destroy(thpool)
end

main()