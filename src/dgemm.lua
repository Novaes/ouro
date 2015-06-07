local IO = terralib.includec("stdio.h")
local stdlib = terralib.includec("stdlib.h")

local function isinteger(x) return math.floor(x) == x end

local NB = 48
terra naivel1matmul(A : &double, B : &double, C : &double, lda : int, ldb : int, ldc : int, alpha : double)
for m = 0, NB do
	for n = 0, NB do
		C[m*ldc + n] = alpha * C[m*ldc + n]
		for k = 0, NB do
			C[m*ldc + n] = C[m*ldc + n] + A[m*lda + k] * B[k*ldb + n]
		end
	end
end
end

function symmat(name,I,...)
	if not I then return symbol(name) end
	local r = {}
	for i = 0,I-1 do
		r[i] = symmat(name..tostring(i),...)
	end
	return r
end

function genl1matmul(NB, NK, RM, RN, V, prefetch, K, L)
	
	-- assert: no boundary cases
	assert(isinteger(NB / RN)) -- number of iterations over b for matrix multiplication
	assert(isinteger(NB / RM)) -- number of iterations over a for matrix multiplication

	local M,N = NB,NB -- new

	local VP = &vector(double,V)
	local terra vecload(data : &double, idx : int)
	var addr = &data[idx]
	return @VP(addr)
end

local terra vecstore(data : &double, idx : int, v : vector(double,V))
var addr = &data[idx]
@VP(addr) = v
end

local A,B,C,mm,nn, alpha = symbol("A"),symbol("B"),symbol("C"),symbol("mn"),symbol("nn"),symbol("alpha")
local lda,ldb,ldc = symbol("lda"),symbol("ldb"), symbol("ldc")
	local a,b,c = symmat("a",RM+2,RN+2), symmat("b",K,L), symmat("c",RM,RN) -- changed
	local kk, ll = symbol("kk"), symbol("ll") -- new
	local x,y = symbol("x"), symbol("y") --new
	local loadkernel,loadA,loadc,storec = terralib.newlist(),terralib.newlist(),terralib.newlist(),terralib.newlist()

	for m = 0, RM+1 do 
		for n = 0, RN+1 do
			loadA:insert(quote
				var [a[m][n]] = vecload(A,( (mm-1) + m)*ldc + (nn-1) + n) -- mm+1 and nn+1 because I start from -1 now 
				end)
			if(m>=0 and m<RM and n>=0 and n<RN) then
				loadc:insert(quote
					var [c[m][n]] = alpha * vecload(C,(mm+m)*ldc + nn + n)
					-- var [c[m][n]] = alpha * @(C + (mm+m)*ldc + nn + n)
					end)
				storec:insert(quote
					vecstore(C,(mm+m)*ldc + nn + n,[c[m][n]])
					-- var c: &double = C + (mm+m)*ldc + nn + n
					-- c =  [c[m][n]]
					end)
			end
		end
	end



	local calcc = terralib.newlist()

	--LOAD FULL KERNEL
	for  k=0, K-1 do
		for l = 0, L-1 do
			loadkernel:insert(quote
				var [b[k][l]] = vecload(B, k*ldb + l)
				end)
		end
	end

	--OPERATION
	for m = 0, RM-1 do
		for n = 0, RN-1 do
			for k=0, K-1 do
				for l = 0, L-1 do
					-- would sum mm or nn, but this position is realtive to this mini-block (rm, rn)
					x, y = m + (k - math.floor(K/2) ), n + (l - math.floor(L/2))
					--if x >= -1 and x<RM+1 and y>=-1 and y<RN+1 then --RM+1 and RN+1 are the image loads, no boundaries cases
					calcc:insert(quote
							--remeber that taking the pos a[x+1][y+1], e.g. a[0][0] menas take a[-1][-1] necessary for c[0][0]
							[c[m][n]] = [c[m][n]] + [a[x+1][y+1]] * [b[l][k]]
							end)
					--end
				end
			end
		end
	end

	return terra([A] : &double, [B] : &double, [C] : &double, [lda] : int, [ldb] : int, [ldc] : int, [alpha] : double)
		-- no borders, original from 0 to NB-1
		for [mm] = 1, NB-2, RM do
			for [nn] = 1, NB-2, RN do
				[loadA];
				[loadc];
				[loadkernel];
				[calcc];
				[storec];
			end
		end
	end
end

local NB2 = 8 * NB
--local l1matmul = genl1matmul(NB, 4, 3, 2, 4, false)
local l1matmul = genl1matmul(NB, 4, 3, 2, 1, false, 3, 3)

terra min(a : int, b : int)
return terralib.select(a < b, a, b)
end

terra my_dgemm(gettime : {} -> double, M : int, N : int, K : int, L: int, alpha : double, A : &double, lda : int, B : &double, ldb : int,
	           beta : double, C : &double, ldc : int, kCenterX: int, kCenterY: int) --do not do borders
	--where should I check for borders (there will be two boundaries, gemm boundary, convolution boundary)
	for mm = 0,M,NB2 do
		for nn = 0,N,NB2 do--substitute here for blocking method
			for m = mm,min(mm+NB2,M),NB do
				for n = nn,min(nn+NB2,N),NB do
					l1matmul(A + m*lda + n,
					         B,-- B is always the same
					         C + m*ldc + n,
					         lda,ldb,ldc,0)
							-- terralib.select(k == 0,0,1) if it is the first posiiton, 
							-- no because it will be in cache already, but if it is not. 
							-- Do it, because this prefetch will but the next line already in the cache
						end
					end
				end
			end
		end

		terralib.saveobj("my_dgemm.o", {my_dgemm = my_dgemm})