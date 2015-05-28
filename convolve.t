require "lib/image"

local cstdio = terralib.includec("stdio.h")
local cstdlib = terralib.includec("stdlib.h")
local cmath = terralib.includec("math.h")

struct Filter { 
  width : int, 
  height : int, 
  stride : int, 
  divisor: int,
  floating : bool, -- is this floating point?
  weights : &double
}

terra Filter:init(
  width : int,
  height : int,
  stride : int,
  divisor: int,
  floating : bool,
  weights : &double)

  self.width = width
  self.height = height
  self.stride = stride
  self.floating = floating
  self.weights = weights
  return self
end

terra Filter:load(argc: int, argv: &rawstring)
  var width : int = cstdlib.atoi(argv[1])
  var height: int = cstdlib.atoi(argv[2])
  var divisor: int = cstdlib.atoi(argv[3])
  var weights = [&double](cstdlib.malloc(width * height * sizeof(double)))
  for i=4,argc do
    weights[i-4] = cstdlib.atoi(argv[i])
    --cstdio.printf("weights[%d]: %d\n", i-4, weights[i-4]);
  end
  return self:init(width,height,width,divisor,false,weights)
end


terra Image:slowconvolve(ker: Filter)
  var iRows = self.width
  var iCols = self.height
  var channels = self.channels
  var kRows = ker.width
  var kCols = ker.height
  var kCenterX : int, kCenterY: int = cmath.floor(kRows/2), cmath.floor(kCols/2)
  var kernel: &double = ker.weights
  var data: &uint8 = [&uint8](self.data)
  var output : &int = [&int](cstdlib.malloc(self.width * self.height * sizeof(int)))  
  
  --[[ if self.dataPtr ~= self.data then cstdio.printf("STRIDED IMAGE")  end ]]
  --[[cstdio.printf("basic data") cstdio.printf("%d %d %d %d %d %d %d",iRows, iCols, channels, kRows, kCols, kCenterX, kCenterY) ]]
  
  -- TODO test separated, using it to fit in a L1 cache during the loop
  var r = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))
  var g = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))
  var b = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))

  if(self.channels == 3) then
    var i=0
    for j=0,(self.width * self.height * channels),channels do
      r[i] = data[j]
      g[i] = data[j+1]
      b[i] = data[j+2]
      i = i+1
    end
  end

  var inputData : &uint8
  var outputData : &uint8

  --[[ -- check colors
  for i=0, self.width * self.height do
    cstdio.printf(" %d ",r[i])
  end
  ]]

  var Rout = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))
  var Gout = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))
  var Bout = [&uint8](cstdlib.malloc(self.width * self.height * sizeof(uint8)))

--[[
  cstdio.printf("weights")
  for i=0, kRows do
      for j=0, kCols do cstdio.printf("\n%d",ker.weights[kCols * i + j]) end 
      cstdio.printf("\n")
  end
]]

  for p=0,channels do

    if p == 0 then
      inputData = r
      outputData = Rout
    elseif p == 1 then
      inputData = g
      outputData = Gout
    else 
      inputData = b
      outputData = Bout
    end

    var sum: double
    var mm: int 
    var ii: int
    var jj: int

    for i=0, iRows do
      for j=0, iCols do
        sum = 0
        for m=0, kRows do
          mm = kRows - 1 - m
          for n=0,kCols do
            var nn: int = kCols - 1 - n
            -- flipped kernel
            ii = i + (m - kCenterY)
            jj = j + (n - kCenterX)
            if ii>=0 and ii<iRows and jj>=0 and jj<iCols then
              sum = sum + [double]([double](inputData[ii * iCols + jj]) * kernel[mm * kCols + nn])
            end
          end
        end
        outputData[i*iCols + j] =  sum
      end
    end
  end

  -- re-construct img data
  if (channels == 3) then
    for j=0,(self.width * self.height) do --same type to assingment
      --output[3*j], output[(3*j)+1], output[(3*j)+2] = Rout[j], Gout[j], Bout[j]
      data[3*j], data[(3*j)+1], data[(3*j)+2] = Rout[j], Gout[j], Bout[j]
    end
  end

  return output
end

terra Filter:load()
  var width : int
  var height: int
  var divisor: int
  var weights = [&double](cstdlib.malloc(width * height * sizeof(double)))
  cstdio.scanf("%d%*c",&width)
  cstdio.scanf("%d%*c",&height)
  cstdio.scanf("%d%*c",&divisor)
  for i=0,25 do
    cstdio.scanf("%d%*c",&weights[i])
  end
  return self:init(width,height,width,divisor,false,weights)
end

local terra loadAndRun(argc: int, argv: &rawstring) --gettime : {} -> double

  --loading image
  var inp: Image
  inp:load("lena.bmp")

  --loading kernel
  var kernel: Filter
  kernel:load(argc, argv)
  --kernel:load()

  inp:slowconvolve(kernel)

  var out: Image
  out = inp
  out:save("lena_out.bmp")
  out:toFloat32()
  out:save("lena_out2.bmp")
  out:free()
  return 0
end

terralib.saveobj("my_convolution.o",{ loadAndRun = loadAndRun })

--run until it does not compare with another image convolution
terralib.saveobj("my_convolution",{ main = loadAndRun })