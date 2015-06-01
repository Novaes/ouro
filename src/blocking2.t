IO = terralib.includec("stdio.h")

terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

terra printm(m: int[3][3], M: int, N: int)
    for i=0,N do
        for j=0,N do
            IO.printf("%d",m[i][j])
        end
        IO.printf("\n")
    end 
    IO.printf("\n")
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

terra main()
	var N = 3

	--fill image
    var img : int[3][3]
	for i=0,3 do
		for j=0,3 do
			img[i][j] = i*3 + j+1
		end
	end

	--flipped kernel example
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

	var out : int[3][3]
    for i=0,3 do
		for j=0,3 do
			out[i][j] = 0
		end
	end

--	check matrices
	printm(img,N,N)
    printm(ker,N,N)
    printm(out,N,N)

--	code
	var fRows: int, fCols: int = 3, 3	
	var x: int, y: int

-- dynamic attrs of kernel
	var kCenterX: int = 1
	var kCenterY: int = 1

--change this N later
	[blockedloop(N, {1,2,3}, function(i,j,ki,kj)
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

main()