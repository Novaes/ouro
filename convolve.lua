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

terra Filter:print()
  cstdio.printf("weights: \n")
  for i=0, self.width do
      for j=0, self.height do 
        cstdio.printf(" %f ", self.weights[self.height * i + j]) 
      end 
      cstdio.printf("\n")
  end
end

terra Filter:flip()
    var tmp: &double  = [&double](cstdlib.malloc(self.width * self.height * sizeof(double)))
    for i=0, self.width do
      var ii: int = self.width - 1 - i
      for j=0, self.height do 
        var jj: int = self.height - 1 - j
        tmp[self.height * i + j] = self.weights[self.height * ii + jj]
      end 
    end
    self:free()
    self.weights = tmp
end

terra memRGB(width: int, height: int)
  return ([&uint8](cstdlib.calloc(width * height, sizeof(uint8)))), 
  ([&uint8](cstdlib.calloc(width * height, sizeof(uint8)))), 
  ([&uint8](cstdlib.calloc(width * height, sizeof(uint8))))
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

terra backToImage(data: &uint8, channels: int, width: int, height: int, Rout: &uint8, Gout: &uint8, Bout: &uint8)
  if (channels == 3) then
    for j=0,(width * height) do
      data[3*j], data[(3*j)+1], data[(3*j)+2] = Rout[j], Gout[j], Bout[j]
    end
  end
end

terra min(a : int, b : int)
  return terralib.select(a < b, a, b)
end

-- N to images, M to kernel
function blockedloop(N,M,blocksizes,bodyfn)
    local function generatelevel(n,ii,jj,kk,ll,bb0,bb1,bb2,bb3)
        if n > #blocksizes then
         return bodyfn(ii,jj,kk,ll)
        end
        local blocksize = blocksizes[n]
        return quote
            for i = ii,min(ii+bb0,N),blocksize do
                for j = jj,min(jj+bb1,N),blocksize do
                    for k = kk,min(kk+bb2,M),blocksize do
                        for l = ll,min(ll+bb3,M),blocksize do
                            [ generatelevel(n+1,i,j,k,l,blocksize,blocksize,blocksize,blocksize) ]
                        end
                    end
                end
            end
        end
    end
    return generatelevel(1,0,0,0,0,N,N,M,M)
end

function naiveconvolution() 
  var iRows, var iCols = self.width, self.height 
  var kRows, kCols = ker.width, ker.height
  var channels = self.channels
  var kCenterX : int, kCenterY: int = cmath.floor(kRows/2), cmath.floor(kCols/2)
  var weights: &double = ker.weights
  var data: &uint8 = [&uint8](self.dataPtr)
  var r: &uint8, g: &uint8, b: &uint8 = generateRGB(self.width,self.height,data)
  var Rout: &uint8, Gout: &uint8, Bout: &uint8 = memRGB(self.width,self.height)
  var img : &uint8, out : &uint8

  for p=0,channels do
    if p == 0 then
      img = r
      out = Rout
    elseif p == 1 then
      img = g
      out = Gout
    else 
      img = b
      out = Bout
    end
    var x: int, y: int
    var N: int, M: int = self.width, ker.width
    for i=0, iRows do
      for j=0, iCols do
        for m=0, kRows do
          for n=0,kCols do
            -- boundaries
            ii = i + (m - kCenterY)
            jj = j + (n - kCenterX)
            if ii>=0 and ii<iRows and jj>=0 and jj<iCols then
              out[i*iCols + j] = out[i*iCols + j] + [double]([double](img[ii * iCols + jj]) * weights[m * kCols + n])
            end
          end
        end
      end
    end
  end

  backToImage(data,channels,self.width,self.height,Rout,Gout,Bout)
  free(r,b,g,Rout,Gout,Bout)
end

terra Image:convolve(ker: Filter)
  var iRows = self.width
  var iCols = self.height
  var channels = self.channels
  var kRows = ker.width
  var kCols = ker.height
  var kCenterX : int, kCenterY: int = cmath.floor(kRows/2), cmath.floor(kCols/2)
  var weights: &double = ker.weights
  var data: &uint8 = [&uint8](self.dataPtr)
  -- check arguments
  -- if self.dataPtr ~= self.data then cstdio.printf("STRIDED IMAGE")  end 
  --cstdio.printf("basic data") cstdio.printf("%d %d %d %d %d %d %d \n",iRows, iCols, channels, kRows, kCols, kCenterX, kCenterY)
  var r: &uint8, g: &uint8, b: &uint8 = generateRGB(self.width,self.height,data)
  var img : &uint8, out : &uint8
  var Rout: &uint8, Gout: &uint8, Bout: &uint8 = memRGB(self.width,self.height)

  --ker:print()

  for p=0,channels do
    if p == 0 then
      img = r
      out = Rout
  --  for i=0,(kRows) do for j=0,(kCols) do cstdio.printf("%d ", out[i*kCols + j])  end cstdio.printf("\n") end 
    elseif p == 1 then
      img = g
      out = Gout
    else 
      img = b
      out = Bout
    end
    
    -- for j=0,(iRows * iCols) do cstdio.printf(" %d ", img[j])  end 
    var x: int, y: int
    var N: int, M: int = self.width, ker.width
    --todo: convolution dividing by ker.divisor in the final
    [blockedloop(N, M, {1,2,3}, function(i,j,ki,kj)
      return
        quote
          -- IO.printf("%d %d %d %d\n",i,j,ki,kj)
          x, y = i + (ki - kCenterY), j + (kj - kCenterX)
          if x >= 0 and x<iRows and y>=0 and y<iCols then
            out[i * iCols + j] = out[i * iCols + j] + img[x * iCols + y] * weights[ki * kCols + kj];
          end
        end
      end)
    ]
  end
  backToImage(data,channels,self.width,self.height,Rout,Gout,Bout)
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
  var ker: Filter
  ker:load(argc, argv)

  -- convolve
  ker:flip()
  inp:convolve(ker)

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