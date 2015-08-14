local dp256 = terralib.intrinsic("llvm.x86.avx.dp.ps.256", 
  { vector(float,8), vector(float,8), int32 } ->  vector(float,8) )

local stdio = terralib.includec("stdio.h")

local ffi = require("ffi")
d1 = ffi.new("float[?] __attribute__((aligned(32)))",11)
d2 = ffi.new("float[?] __attribute__((aligned(32)))",11)

for i = 0,10 do
  d1[i] = i
end

for i = 0,10 do
  d2[i] = 2
end

terra dot13(v1: &float, v2: &float, size: int)
  if (size == 1) then
    return @(v1) * @(v2)
  elseif (size == 2) then
    return dp256(@[&vector(float,8)](v1), @[&vector(float,8)](v2), 0x31)[0]
  else
    return dp256(@[&vector(float,8)](v1), @[&vector(float,8)](v2), 0x71)[0]
  end
end

terra dot4(v1: &float, v2: &float)
  return dp256(@[&vector(float,8)](v1), @[&vector(float,8)](v2), 0xF1)[0]
end

terra dot8(v1: &float, v2: &float)
  var v = dp256(@[&vector(float,8)](v1), @[&vector(float,8)](v2), 0xF1)
  return v[0] + v[4]
end

terra dotprod(v1: &float, v2: &float, size: int): float
  if (size <= 0) then
    return 0
  elseif (size <= 3) then
    return dot13(v1, v2, size)
  elseif (size < 8) then
    return dot4(v1, v2) + dotprod(v1 + 4, v2 + 4, size - 4)
  end

  return dot8(v1, v2) + dotprod(v1 + 8, v2 + 8, size - 8)
end

local x = dotprod(d1, d2, 11)
print(x)
