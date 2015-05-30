require "lib/image"
require "lib/time"

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

terra memRGB(width: int, height: int)
  return ([&uint8](cstdlib.malloc(width * height * sizeof(uint8)))), 
  ([&uint8](cstdlib.malloc(width * height * sizeof(uint8)))), 
  ([&uint8](cstdlib.malloc(width * height * sizeof(uint8))))
end

terra Filter:load(argc: int, argv: &rawstring)
  var width : int = cstdlib.atoi(argv[1])
  var height: int = cstdlib.atoi(argv[2])
  var divisor: int = cstdlib.atoi(argv[3])
  var weights = [&double](cstdlib.malloc(width * height * sizeof(double)))
  for i=4,argc do
    weights[i-4] = cstdlib.atoi(argv[i])
  end
  return self:init(width,height,width,divisor,false,weights)
end

terra Filter:free()
  cstdlib.free(self.weights)
end

terra generateRGB(w: int, h: int, data: &uint8): {&uint8, &uint8, &uint8}
  var r: &uint8, g: &uint8, b: &uint8 = memRGB(w,h)
  var i=0
  for j=0,(w * h * 3),3 do
    r[i] = data[j]
  	g[i] = data[j+1]
  	b[i] = data[j+2]
  	i = i+1
  end
  return r,g,b	
end

terra free(r: &uint8, g: &uint8, b: &uint8, Rout: &uint8, Gout: &uint8, Bout: &uint8)
  cstdlib.free(r)
  cstdlib.free(g)
  cstdlib.free(b)
  cstdlib.free(Rout)
  cstdlib.free(Gout)
  cstdlib.free(Bout)
end

terra bound(d: int, max: int) : uint8
	if d < 0 then
		return 0
	elseif d > max then
		return max
	else 
		return d
	end
end

terra Image:slowconvolve(ker: Filter)
  var iRows = self.width
  var iCols = self.height
  var channels = self.channels
  var kRows = ker.width
  var kCols = ker.height
  var kCenterX : int, kCenterY: int = cmath.floor(kRows/2), cmath.floor(kCols/2)
  var kernel: &double = ker.weights
  var data: &uint8 = [&uint8](self.dataPtr)
  
  -- check arguments
  -- if self.dataPtr ~= self.data then cstdio.printf("STRIDED IMAGE")  end 
  --cstdio.printf("basic data") cstdio.printf("%d %d %d %d %d %d %d \n",iRows, iCols, channels, kRows, kCols, kCenterX, kCenterY)
  
  var r: &uint8, g: &uint8, b: &uint8 = generateRGB(self.width,self.height,data)
  var inputData : &uint8, outputData : &uint8
  var Rout: &uint8, Gout: &uint8, Bout: &uint8 = memRGB(self.width,self.height)
  
  --[[
  --check weights
  cstdio.printf("weights: \n")
  for i=0, kRows do
      for j=0, kCols do cstdio.printf(" %d ",[int](ker.weights[kCols * i + j])) end 
      cstdio.printf("\n")
  end]]


--[[for j=0,(iRows * iCols) do
      cstdio.printf(" %d ", r[j])
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
--[[  
  for j=0,(iRows * iCols) do
    cstdio.printf(" %d ", inputData[j])
  end
]]
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
        -- outputData[i*iCols + j] = bound(cmath.floor(sum/ker.divisor),255)
        outputData[i*iCols + j] = sum
      end
    end
  end

  -- re-construct img data
  if (channels == 3) then
    for j=0,(self.width * self.height) do --same type to assingment
      data[3*j], data[(3*j)+1], data[(3*j)+2] = Rout[j], Gout[j], Bout[j]
    end
  end

  free(r,b,g,Rout,Gout,Bout)
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

local terra loadAndRun(argc: int, argv: &rawstring)

  -- loading image
  var inp: Image
  inp:load("images/lena.bmp")

  -- loading kernel
  var kernel: Filter
  kernel:load(argc, argv)
  --kernel:load()
  
  -- convolve
  inp:slowconvolve(kernel)
  
  -- print(runBenchmark(testing))

  var out: Image
  out = inp
  out:save("images/lena_out.bmp")
  out:free()
  return 0
end

terralib.saveobj("my_convolution.o",{ loadAndRun = loadAndRun })
--run until it does not compare with another image convolution
terralib.saveobj("my_convolution",{ main = loadAndRun })