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

local function unalignedstore(addr,v)
	return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore) 

terra main()
	var A : &double = [&double](cstdlib.malloc( 8 * sizeof(double)) )
	var count = 1
	for i=0,7 do
		@(A+i) = count
		count = count + 1
	end

	for i=0,7 do
		cstdio.printf("A: %f",A[i])
	end
	cstdio.printf("\n")

	--loadc 
	var v0 : vector(double,4) = unalignedload(VP(&A[0]))
	var v1 : vector(double,4) = unalignedload(VP(&A[4]))
	for i=0,V do cstdio.printf("%f \n",v0[i]) end
	for i=0,V do cstdio.printf("%f \n",v1[i]) end
	
	-- calcc
	v0 = -1 v1 = -1

	--storec
	unalignedstore(VP(&A[0]),v0)
	unalignedstore(VP(&A[4]),v1)
	for i=0,V do
		cstdio.printf("%f \n",v0[i])
	end

	for i=0,V do
		cstdio.printf("%f \n",v1[i])
	end

end

main()