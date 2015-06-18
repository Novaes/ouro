local C = terralib.includecstring [[

#include <stdint.h>
#include <sys/time.h>
#include <time.h>

const uint64_t msInNs = 1000;

inline uint64_t nowNs() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((uint64_t)tv.tv_usec) * msInNs + tv.tv_usec; 
}
]]

local S = require "std"

terra nowNs()
  return C.nowNs()
end

terra runBenchmark(fn : {&opaque} -> &opaque, arg : &opaque)
  var maxIter : int = 5
  var maxTime : int = 1000000000 -- 1 second

  -- collect execution times
  var times = [S.Vector(uint64)].salloc():init()
  var i : int = 0
  var totalTime : uint64 = 0
  while i < maxIter and totalTime < maxTime do
    var start : uint64 = nowNs()
    fn(arg)
    var elapsed : uint64 = nowNs() - start

    -- insert in order - no sort in terra :(
    var j : int = 0
    while j < i do
      if (elapsed <= times(j)) then
        break
      end
      j = j + 1
    end
    times:insert(j, elapsed)

    i = i + 1
    totalTime = totalTime + elapsed
  end

  -- discard outliers
  var size : int = times:size() - 1
  while size >= 0 do
    if times(size) <= 2 * times(0) then
      break
    end
    size = size - 1
  end
  size = size + 1

  -- computes p50 - will change to mode soon
  var p50 : uint = size * 0.5
  return times(p50)
end

 -- EXAMPLE
--[[
terra testing(p : &opaque) : &opaque
  print(p) -- example of a function
  nowNs()
  return nil
end

terra main()
  print(runBenchmark(testing, nil))
end

main()
]]