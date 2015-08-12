cstdio = terralib.includec("stdio.h")
cstdlib = terralib.includec("stdlib.h")

--local number = double
local alignment = 8

local V = 4
local VP = &vector(double,V)

local terra vecload(data : &double, idx : int)
	var addr = &(data[idx])
	return @VP(addr)
end

local terra vecstore(data : &double, idx : int, v : vector(double,V))
	var addr = &data[idx]
	@VP(addr) = v
end

local function unalignedload(addr)
	return `terralib.attrload(addr, { align = alignment })
end

unalignedload = macro(unalignedload)

terra main()
	var A : &double = [&double](cstdlib.malloc( 8 * sizeof(double)))
	var count = 1
	for i=0,7 do
		@(A+i) = count
		count = count + 1
	end

	for i=0,7 do
		cstdio.printf("A: %f",A[i])
	end
	cstdio.printf("\n")

	var v0 : vector(double,4) = unalignedload(VP(&A[0]))
	var v1 : vector(double,4) = unalignedload(VP(&A[3]))
	for i=0,V do
		cstdio.printf("%f \n",v0[i])
	end

	for i=0,V do
		cstdio.printf("%f \n",v1[i])
	end

end

main()