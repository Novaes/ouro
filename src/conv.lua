local fRows, fCols, kRows, kCols, channels, output

terra min(a : int, b : int)
  return terralib.select(a < b, a, b)
end

function blockedloop(N,blocksizes,bodyfn)
  local function generatelevel(n,ii,jj,bb)
    if n > #blocksizes then
      return bodyfn(ii,jj)
    end
    local blocksize = blocksizes[n]
    return quote
      -- receive n, size of block, body fn, -- it will use the same blocksize for x and y 
      -- TEST FOR THE BEST ONE
      for i = ii,min(ii+bb,N),blocksize do
        for j = jj,min(jj+bb,N),blocksize do
          [ generatelevel(n+1,i,j,blocksize) ]
        end
      end
    end
  end
  return generatelevel(1,0,0,N)
end

function slide(inputData, kernel, rm, rn, rmm, rnn)
  output = {{0,0,0},{0,0,0},{0,0,0}}
  local t1, t2, t3, t4 = 0, 0, 0, 0 -- index
  local kCenterX, kCenterY = math.floor(kRows/2), math.floor(kCols/2)
  for t1=0, fRows-1, t1+rm do
    for i=t1, math.min(fRows, t1+rm)-1 do
      for t2=0, fCols-1, t2+rn do
        for j=t2, math.min(fCols, t2+rn)-1 do
          for t3=0, kRows-1, t3+rmm do
            for m=t3, math.min(kRows, t3+rmm)-1 do
              for t4=0, kCols-1, t4+rnn do
                for n=t4, math.min(kCols, t4+rnn)-1 do
                  --respect constraint
                  --do what java code implementation does lower than bound return 0, greater than bound, return max
                  local ii, jj = i + (m - kCenterY), j + (n - kCenterX)
                  if ii >= 0 and ii<fRows and jj>=0 and jj<fCols then
                    output[i+1][j+1] = output[i+1][j+1] + inputData[ii+1][jj+1] * kernel[m+1][n+1];
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


function naiveConvolve(inputData, kernel)
  local kCenterX, kCenterY = math.floor(kRows/2), math.floor(kCols/2)
  output = {{0,0,0},{0,0,0},{0,0,0}}
  for i=0, fRows-1 do
    for j=0, fCols-1 do
      for m=0,kRows-1 do
        local mm = kRows - 1 - m
        for n=0,kCols-1 do
          local nn = kCols - 1 - n
          --boundaries
          local ii = i + (m - kCenterY)
          local jj = j + (n - kCenterX)
          if ii >= 0 and ii< fRows and jj>=0 and jj<fCols then
            local tmp = output[i+1][j+1]
            output[i+1][j+1] = tmp + inputData[ii+1][jj+1] * kernel[mm+1][nn+1];
          end
        end
      end
    end
  end
end

--[[
local function main()  
  --get image, decompose on RGB, do the convolution operation
  local kernel = getKernel()
  local image = getImage()
  
  --it should be a method of image that receives a kernel flipped: flip() is a kernel method to flip
  kernel = flipkernel(kernel,kRows,kCols)
  
  local blockM = {1,2}
  local blockN = {1,2}
  local blockMM = {1,2}
  local blockNN = {2,1}
  
  for _,rm in ipairs(blockM) do
    for _,rn in ipairs(blockN) do
      for _,rmm in ipairs(blockMM) do
        for _,rnn in ipairs(blockNN) do 
          for i=1,channels do
            io.write("Parameters, rm: " .. rm .. " rn: " .. rn .. " rmm: " .. rmm .. " rnn: " .. rnn .. "\n")
            slide(image, kernel,rm,rn,rmm,rnn)
            check()
          end
        end  
      end  
    end
  end 
end
]]

function equal(a,b,M,N)
  for i=0,M-1 do
    for j=0,N-1 do
      if(a[i+1][j+1] ~= b[i+1][j+1]) then
          io.write(a[i+1][j+1] .. " != " .. b[i+1][j+1] .. "\n")
          return false
      end
    end
  end
  return true
end

function printMatrix(matrix,M,N)
  for i=0,M-1 do
    for j=0,N-1 do
      io.write(" " .. matrix[i+1][j+1])
    end
    io.write("\n")
  end
  io.write("\n")
end

--flip in terra
function flipkernel(kernel,M,N)
  local newKernel = {}
  for m=0,M-1 do
    newKernel[m+1] = {}
    local mm = M - 1 - m
    for n=0,N-1 do 
      local nn = N - 1 - n
      newKernel[m+1][n+1] = kernel[mm+1][nn+1]
    end
  end
  return newKernel
end



function check()
  local test, tRows, tCols = {  {-13,-20,-17},
                                {-18,-24,-18},
                                {13,  20, 17}}, 3, 3

  for i=1,tRows do
    for j=1,tCols do
      if output[i][j] ~= test[i][j] then
        print "error"
        return false
      end
    end
  end
end

function getImage()
  fRows, fCols, channels = 3, 3, 1
  return {{1,2,3},
          {4,5,6},
          {7,8,9}}
end

function getKernel()
  kRows, kCols = 3,3
  return {{-1,-2,-1},
          {0,  0, 0},
          {1,  2, 1}}
end

function ioread()
  io.read ("*n", "*n")
  if (kRows <= 0) or bit.ban(kRows,1) ~= 1 or
    ((kCols <= 0) or ((bit.band(kCols,1)) ~= 1)) then
    print("We assume odd size filter "+
    "in order to have a center element and do the flip")
  end
end

-- multistage programming
IO = terralib.includec("stdio.h")

function sayhello()
  return (quote 
    var x:int 
    x = 2
    IO.printf("%d\n",x)
  end)
end

terra printl()
  [sayhello()]
end

printl()

terra main()
  var a : int[8][8]
  var c = 0
  var N = 8
  [blockedloop(N, {4,2,1}, function(i,j)
      return 
        quote
          --IO.printf("%d %d\n",i,j)
          a[i][j] = c
          c = c + 1
        end
    end)
  ]

  --display
  for i = 0,N do
    for j = 0,N do
      IO.printf("%d\t",a[i][j])
    end
    IO.printf("\n")
  end
end

main()