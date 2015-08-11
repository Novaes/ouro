lib.includec("stdio.h")
local number = double

terra pmatrix(a: double[3][3], M: int, N: int)
  for i=0,M do
    for j=0,N do
      IO.printf("%lf ",a[i][j])
    end
    IO.printf("\n")
  end
end

local number = double
local alignment = 8

function symmat(name,I,...)
    if not I then return symbol(name) end
    local r = {}
    for i = 0,I-1 do
        r[i] = symmat(name..tostring(i),...)
    end
    return r
end

local function isinteger(x) return math.floor(x) == x end

local function unalignedload(addr)
    return `terralib.attrload(addr, { align = alignment })
end
local function unalignedstore(addr,v)
    return `terralib.attrstore(addr,v, { align = alignment })
end

unalignedload,unalignedstore = macro(unalignedload),macro(unalignedstore)

function printvec(C,N)
    for i=0,N-1 do
        print(C[i])
    end
end

-- Simple Test for vector load
function genkernel()
    --vector(TYPE, #NUMBER)
    local V = 1
    local RM = 20
    local C = symbol("C")
    local c,caddr = symmat("c",RM), symmat("caddr",RM)
    for i=0,RM-1 do
        caddr[i] = i+1
    end

    -- printvec(caddr,RM)
    local vector_type = vector(number,V)
    local vector_pointer = &vector_type
    local loadc = terralib.newlist()

    for m=0,(RM-1)*V do
        loadc:insert(quote
            var [caddr[m]] = C + m*V
            IO.printf("%d",caddr[m])
            var [c[m]] = unalignedload(vector_pointer([caddr[m]]))
        end)
    end
    return terra(C: &number) [loadc]; end
end

function genmatrix()
    local x = genkernel()
    --x(img)
end

-- Simple Test for prefetch
function prefetch()

end

-- Simple Test for varying vector instruction V


-- Simple Test for unalign load

genmatrix()

-- tests of structures on gemm.t
function l1kernel()
    local loadc = terralib.newlist()
    loadc:insert(quote
        IO.printf("::genkernel2")
    end)
    loadc:insert(quote
        IO.printf(" executed::\n")
    end)
    return terra([t]: &double) [loadc]; IO.printf("::argument passed was %d::\n",t) end
end

function genmatrix()
    local y = l1kernel()
    return y(55)
end

function testvector()
    --vector(TYPE, #NUMBER)
    local V = 5
    local vector_type = vector(double,V)
    local vector_pointer = &vector_type

    local RM,RN = 2,3
    local c = symmat("c",RM,RN)
    return (quote
        @vector_type(1) = 2
        @vector_type(6) = 3
    end)
end

--local test = require("test")
--test.eq(foo(),8)



--testvector():printpretty()
