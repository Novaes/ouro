terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

function blockedloop(N,blocksizes,bodyfn)
    local function generatelevel(n,ii,jj,kk,ll,bb0,bb1,bb2,bb3)
        if n > #blocksizes then
          return bodyfn(ii,jj,kk,ll)
        end
        local blocksize = blocksizes[n]
        return quote
            for i = ii,min(ii+bb0,N),blocksize do
                for j = jj,min(jj+bb1,N),blocksize do
                    for k = kk,min(kk+bb2,N),blocksize do
                        for l = ll,min(ll+bb3,N),blocksize do
                            [ generatelevel(n+1,i,j,k,l,blocksize,blocksize,blocksize,blocksize) ]
                        end
                    end
                end
            end
        end
    end
    return generatelevel(1,0,0,0,0,N,N,N,N)
end

IO = terralib.includec("stdio.h")


terra main()
	var N = 3

--	======================================================================
	var img : int[3][3]
	for i=0,3 do
		for j=0,3 do
			img[i][j] = i*3 + j+1
		end
	end
--	======================================================================
	--flipped kernel already
	var ker : int[3][3]
	ker[0][0] = 1
	ker[0][1] = 2
	ker[0][2] = 1
	ker[1][0] = 0
	ker[1][1] = 0
	ker[1][2] = 0
	ker[2][0] = -1
	ker[2][1] = -2
	ker[2][2] = -1
--	======================================================================
	var out : int[3][3]
	for i=0,3 do
		for j=0,3 do
			out[i][j] = 0
		end
	end
--	======================================================================
	for i=0,N do
		for j=0,N do
			IO.printf("%d",img[i][j])
		end
		IO.printf("\n")
	end	
	IO.printf("\n")
	for i=0,N do
		for j=0,N do
			IO.printf("%d",ker[i][j])
		end
		IO.printf("\n")
	end
	IO.printf("\n")
	for i=0,N do
		for j=0,N do
			IO.printf("%d",out[i][j])
		end
		IO.printf("\n")
	end
--	======================================================================
	var fRows: int, fCols: int = 3, 3	
	var x: int, y: int

	-- dynamic attrs of kernel
	var kCenterX: int = 1
	var kCenterY: int = 1

	--change this N later
	[blockedloop(N, {1, 2}, function(i,j,ki,kj)
		return
			quote
				-- IO.printf("%d %d %d %d\n",i,j,ki,kj)
				x, y = i + (ki - kCenterY), j + (kj - kCenterX)
				if x >= 0 and x<fRows and y>=0 and y<fCols then
                	out[i][j] = out[i][j] + img[x][y] * ker[ki][kj];
              	end
			end
		end)
	]

	--display
	for i = 0,N do
		for j = 0,N do
			IO.printf("%d\t",out[i][j])
		end
		IO.printf("\n")
	end

end

--main()

terra conv2d()
	var N = 3
--	======================================================================
	var img : int[3][3]
	-- 1-9
	for i=0,3 do
		for j=0,3 do
			img[i][j] = i*3 + j+1
		end
	end
--	======================================================================
	--flipped kernel already
	var ker1 : int[3][1]
	ker1[0][0] = 1
	ker1[1][0] = 2
	ker1[2][0] = 1
	var ker2 : int[1][3]
	ker2[0][0] = 1
	ker2[1][0] = 2
	ker2[2][0] = 1
-- ======================================================================
	var out : int[3][3]
	for i=0,3 do
		for j=0,3 do
			out[i][j] = 0
		end
	end
	var out1 : int[1][3]
	for i=0,1 do
		for j=0,3 do
			out1[i][j] = 0
		end
	end

	for i=0,N do
		for j=0,N do
			IO.printf("%d",img[i][j])
		end
		IO.printf("\n")
	end	
	IO.printf("\n")

	for i=0,N do
		for j=0,1 do
			IO.printf("%d",ker1[i][j])
		end
		IO.printf("\n")
	end
	IO.printf("\n")
	for i=0,1 do
		for j=0,N do
			IO.printf("%d",ker2[i][j])
		end
		IO.printf("\n")
	end
	IO.printf("\n")
	for i=0,N do
		for j=0,N do
			IO.printf("%d",out[i][j])
		end
		IO.printf("\n")
	end
-- ======================================================================
	
	--borders??
		--vertical
		for i=0,N do
			for j=0,N do
				out1[0][i] = out1[0][i] + img[j][i] * ker1[j][0]
			end
		end
		var v: int = 0
		--horizontal
		for i=0,N do
			v = v + out1[0][i] * ker2[0][i]
		end
		
	-- for i=0,N do
	-- 	IO.printf("%d ",out1[i][0])
	-- end
	-- IO.printf("\n")

end
conv2d()
--[[
stdlib = terralib.includec("stdlib.h")
function Image(Spectrum)
	local struct ImageImpl {
		data : &Spectrum,
		N : int
	}
	terra ImageImpl:init(N : int)
		self.data = [&float](stdlib.malloc(N*N*sizeof(Spectrum)))
		self.N = N
	end
	ImageImpl.methods.pixel = macro(
	function(self,x,y)
		return `self.data[x*self.N + y]
	end)
	return ImageImpl
end

GreyScaleImage = Image(float)

terra laplace(input : &GreyScaleImage, output : &GreyScaleImage)
	var newN = input.N - 2 --shrink result since we do not calculate boundaries
	output:init(newN);
	[blockedloop(newN,{32,1},function(i,j)
		return quote
			output:pixel(i,j) =
              input:pixel(i+0,j+1)
            + input:pixel(i+2,j+1)
            + input:pixel(i+1,j+2)
            + input:pixel(i+1,j+0)
            - 4 * input:pixel(i+1,j+1)
		end
	end)]
end

laplace:compile()
laplace:printpretty()
]]