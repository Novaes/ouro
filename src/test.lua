IO = terralib.includec("stdio.h")

terra pmatrix(a: int[3][3], M: int, N: int)
	for i=0,M do
		for j=0,N do
			IO.printf("%d",a[i][j])
		end
		IO.printf("\n")
	end	
end

terra main()
	var img : int[3][3]
	for i=0,3 do
		for j=0,3 do
			img[i][j] = i*3 + j+1
		end
	end
	pmatrix(img,3,3)
end

main()
