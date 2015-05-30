terra min(a: int, b: int): int
    if a < b then return a
    else return b end
end

function blocking()
    local N = 5
    local b = 2
    local C = {}
    local x,y = 0,0
    
    for i=0,N,b do
        for j=0,N,b do
            for x=i, min(N,x+b) do
                for y=j, min(N,y+j) do
                    C[y][x] = 1
                end
            end
        end
    end

    for i = 0,N+1 do
        for j = 0,N+1 do
            print(C[i][j])
        end
    end
end
blocking()