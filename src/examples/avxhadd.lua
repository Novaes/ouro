local IO = terralib.includec("stdio.h")
local vdots256 = terralib.intrinsic("llvm.x86.avx.dp.ps.256", { vector(float,8), vector(float,8), int32 } -> float)
local haddavx = terralib.intrinsic("llvm.x86.avx.hadd.ps.256", { vector(float,8), vector(float,8) } -> vector(float,8))

terra foo()
  var v1: vector(float, 8) = 3
  var v2: vector(float, 8) = 3
  for i=0,8 do
    haddavx()
  end
  return vdots256(v1, v2, 0xFF)
end

local x = foo()
print(x)
