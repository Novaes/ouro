--[[-- Changed for the convolution one, add the ll move it to before than generatedgemm

function genkernelCONV(NB, RM, RN, V, alpha, boundary) --kernel is also blocked with RM and RN
    local M,N,K,L boundaryargs

    -- if one of the parameters is lower than NB then receive the usual
    if boundary then
        M,N,K,L = symbol(int64,"M"),symbol(int64,"N"),symbol(int64,"K"),symbol(int64,"L")
        boundaryargs = terralib.newlist({M,N,K,L})
    else
        boundaryargs = terralib.newlist()
        M,N,K,L = NB,NB,NB,NB
    end

    local A,B,C,mm,nn,ld = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("ld")
    local lda,ldb,ldc = symbol("lda"),symbol("ldb"),symbol("ldc")

    --todo B can be blocked to different rmm and rnn dimensions
    local a,addr = symmat("a",RM,RN), symmat("addr",RM,RN)
    local b,baddr = symmat("b",RM,RN), symmat("baddr",RM,RN)
    local c,caddr = symmat("c",RM,RN), symmat("caddr",RM,RN)
    local k,h = symbol("k"),symbol("h")

    local loadc,storec,calcc = terralib.newlist(),terralib.newlist(),terralib.newlist()
    --todo I can put loada and loadb in calcc 
    local loada,loadb = terralib.newlist(),terralib.newlist()

    local VT = vector(number,V)
    local VP = &VT
    
    local calcc = terralib.newlist()
    --assuming V= 1
    --todo ldc == lda 
    for m = 0, RM-1 do
        for n = 0, RN-1 do
            --from memory to table variable that will later be unrolled and use registers
            loadc:insert(quote
                var [caddr[m][n]] = C + m*ldc + n*V
                var [c[m][n]] = alpha * unalignedload(VP([caddr[m][n]]))
            end)

            loada:insert(quote
                var [addr[m][n]] = A + m*lda + n*V
                var [a[m][n]] = unalignedload(VP([addr[m][n]]))
            end))

            --put back in the slow memory, from register vector to memory => from table variable to memory!
            storec:insert(quote
                unalignedstore(VP([caddr[m][n]]),[c[m][n]])
            end)
        end
    end

    -- put it is separated because later it can be extended to RMM and RNN
    -- todo beta variable multiplication factor
    -- NO V
    for k = 0, RM-1 do
        for h = 0, RN-1 do
            --V removed
            loadb:insert(quote
                var [baddr[k][h] = B + k*ldb + h
                var [b[k][h] = 1.0 * unalignedload(VP([baddr[k][h]))
            end)
        end
    end

-- 0 to rn1: var [b[n]] = unalignedload(VP(&B[n*V]))
-- 0 to rm-1: var [a[m]] = VT(A[m*lda])

    for m = 0, RM-1 do
        for n = 0, RN-1 do
            -- kernel stuff
            for k = 0, RM-1 do
                for h = 0, RN-1 do
                    calcc:insert(quote
                        [c[m][n]] = [c[m][n]] + [a[x][y]] * [b[k][h]]
                    end)
                end
            end
        end
    end

    local result = terra([A] : &number, [B] : &number, [C] : &number, kCenterX: number, kCenterY: number,  
    [lda] : int64,[ldb] : int64,[ldc] : int64,[boundaryargs])
        for [mm] = 0, M, RM do
            for [nn] = 0, N,RN do -- assume V=1 first
                [loadc]; -- taking c loading c to the variable and multiplying by alpha
                for [k] = 0, K do -- can use blocking; loop over the kernel
                    for [h] = 0,H do -- can use blocking; loop over the kernel
                        x, y = i + (k - kCenterY), j + (h - kCenterX)
                        if x >= 0 and x<iRows and y>=0 and y<iCols then
                            [loada];
                            [loadb];
                            [calcc];

                        end
                    end
                end
                [storec];
                A = A - K
                B = B - ldb*K + RN
                C = C + RN
            end
            C = C + RM * ldb - N
            B = C + RM * ldb - N
        end
    end
    return result
end

]]
