cstdio = terralib.includec("stdio.h")
cstdlib = terralib.includec("stdlib.h")
cmath = terralib.includec("math.h")
PI = constant(3.141592653589793)
TWOPI = constant(6.28318530717959)

cbit = terralib.includecstring[[

#define SWAP(x, y) do { typeof(x) SWAP = x; x = y; y = SWAP; } while (0)

void reversalbased(int* n, int* nn, int* j, int* m, double* tempr, double* data){
	int i;
	*n = *nn << 1;
	*j = 1;
	for (i = 0; i < *n - 1; i += 2) {
		if (*j > i) {
			SWAP(data[*j - 1],data[i]);
			SWAP(data[*j],data[i+1]);
		}
		*m = *n >> 1;
		while (*m >= 2 && *j > *m) {
			*j -= *m;
			*m >>= 1;
		}
		*j += *m;
	}
}

void reversal(int N, double* data){
	int nn = N;
	int j,m,tempr,n;

	int i;
	n = nn << 1;
	j = 1;
	for (i = 0; i < n - 1; i += 2) {
		if (j > i) {
			SWAP(data[j - 1],data[i]);
			SWAP(data[j],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
}

]]

--[[

	 FFT/IFFT (Numerical Recipes in C, pages 507-508)
	 Inputs:
		&data: array of complex data points of size 2*NFFT+1
		  it MUST have length power of two
		  data[0] is unused 
		  storage: 
			  data[2*n+1] = real(x(n))
			  data[2*n+2] = img(x(n))  
		isign:  
		  if 1, computes the forward (usual) FFT
		  if -1, computes inverse FFT - in this case the output values have
			to be manually normalized by multiplying with 1/NFFT.
	 Outputs:
			data: The FFT or IFFT results are stored in data, overwriting the input. 

]]

terra four1(data : &double, nn : int, isign : int)
	var n : int32 
	var mmax : int32
	var m : int32
	var j : int32
	var istep : int32
	var i : int32
	var wtemp : double
	var wr : double
	var wpr : double
	var wpi : double
	var wi : double
	var theta : double
	var tempr : double
	var tempi : double

     -- bit reversal starting from 1
	cbit.reversalbased(&n,&nn,&j,&m,&tempr,data)
	 -- daniel-lackzson going all way down to 2-point
	 mmax = 2
	 while n > mmax do
	  istep = 2*mmax --left shift
	  theta = TWOPI/(isign*mmax) -- past used #def TWOPI = 2*M_PI
	  wtemp = cmath.sin(0.5*theta)
	  wpr = -2.0*wtemp*wtemp
	  wpi = cmath.sin(theta)
	  wr = 1.0
	  wi = 0.0
	  for m = 1,mmax,2 do
		for i = m,n+1,istep do
			j = i + mmax
			tempr = wr*data[j] - wi*data[j+1]
			tempi = wr*data[j+1] + wi*data[j]
			data[j]   = data[i]   - tempr
			data[j+1] = data[i+1] - tempi
			data[i] = data[i] + tempr
			data[i+1] = data[i+1] + tempi
		end
		wtemp = wr
		wr = wtemp*wpr - wi*wpi + wr
		wi = wi*wpr + wtemp*wpi + wi
	  end
	  mmax = istep
  end
end

terra twidle_factor(k : int, angle : double, data : &double)
	var tw : double[2], v : double[2]
	tw[0] = cmath.cos(k*angle)
	tw[1] = cmath.sin(k*angle)
	v[0] = tw[0]*data[0] - tw[1]*data[1]
	v[1] = tw[0]*data[1] + tw[1]*data[0]
	data[0] = v[0]
	data[1] = v[1]
end

terra FFT4(in0 : &double, in1 : &double, in2 : &double, in3 : &double)
	var v0 : double, v1 : double , v2 : double, v3 : double
	var in0x : double, in0y : double, in1x : double, in1y : double
	var in2x : double, in2y : double, in3x : double, in3y : double

	-- evens
	v0 = in0[0] + in2[0]
	v1 = in1[0] + in3[0]
	in0x = v0 + v1
	in2x = v0 - v1

	-- imaginary
	v0 = in0[1] + in2[1]
	v1 = in1[1] + in3[1]
	in0y = v0 + v1
	in2y = v0 - v1

	-- odds
	v2 = in0[0] - in2[0]
	v3 = in1[1] - in3[1]
	in1x = v2 + v3 -- in1.x = v2.x + v3.y
	in3x = v2 - v3 -- in3.x = v2.x - v3.y

	-- imaginary
	v2 = in0[1] - in2[1] -- v2.x
	v3 = in3[0] - in1[0]
	in1y = v2 + v3 -- in1.y = v2.y - v3.x
	in3y = v2 - v3 -- in3.y = v2.y - v3.x

	--writing it back
	in0[0] = in0x		in0[1] = in0y
	in1[0] = in1x		in1[1] = in1y
	in2[0] = in2x		in2[1] = in2y
	in3[0] = in3x		in3[1] = in3y
end

terra FFT2(in0 : &double, in1 : &double)
	var v0 : double

	v0 = in0[0]
	in0[0] = v0 + in1[0]
	in1[0] = v0 - in1[0]

	v0 = in0[1] -- same for imaginary
	in0[1] = v0 + in1[1]
	in1[1] = v0 - in1[1]
end

--[[ Ns  : number of stages, N : number of elements on this kernel , K_W : kernel size ]]
terra FFT_2(id: int, start : int, datain : &double, Ns : int, NFFT : int)
	var gId : int = id % NFFT/2
	var s : int = start
	var stride = 2*Ns -- +1 due the imaginary representation
	-- cstdio.printf("init: %d %d\n",s, Ns)

	var in0 : double[2], in1 : double[2]

	-- load reals and imaginaries
	in0[0] = datain[s] 	in0[1] = datain[s + 1]
	in1[0] = datain[(1*stride)+s]	in1[1] = datain[(1*stride)+s + 1]

	-- at first stage twiddle factors is always 0
	-- if Ns ~=1 then
		-- var angle : double = -2*PI*(gId)/(NFFT)
		-- cstdio.printf("angle: %f\n",angle)
		-- twidle_factor(1, angle, [&double](in1))
	-- end

	-- cstdio.printf("in%d: %f %f\n",0,in0[0],in0[1])
	-- cstdio.printf("in%d: %f %f\n",1,in1[0],in1[1])

	FFT2([&double](&in0), [&double](&in1))

	-- cstdio.printf("--> %d %d\n",(0*stride)+s,(1*stride)+s)	
	datain[s] = in0[0]	  datain[s + 1] = in0[1]
	datain[(1*stride)+s] = in1[0]	  datain[(1*stride)+s + 1] = in1[1]
end

-- #define Local_Stride 1
-- #define BLOCK_SIZE 16
-- added id and K_W
-- added lsize
terra FFT_64(id : int, datain : &double, dataout : &double, Ns : int, N : int, K_W : int, Local_Stride : int, BLOCK_SIZE : int, lsize : int)
	var sharedx : double[128*8]
	var sharedy : double[128*8]
	var in_stride : int = N / (128*8) 
	var gId : int = id -- get_global_id(0)
	-- get local size
	var Local_Size : int = lsize

	var in0 : double[2], in1 : double[2], in2 : double[2], in3 : double[2], in4 : double[2], in5 : double[2], in6 : double[2], in7 : double[2]

	var internal_blocksize : int = 8
	var block_id : int = gId/BLOCK_SIZE
	var internal_block_id : int = block_id%internal_blocksize
	-- var tidx : int = gId&(BLOCK_SIZE-1); -- binary and
	var tidx : int = gId 

	do -- 8-point kernel
		var new_block_id = 0 -- make map_id function
		-- var new_block_id : int = map_id(block_id, in_stride,  internal_blocksize);
		var new_gId : int = new_block_id*BLOCK_SIZE+tidx;
		in0 = datain[(0*N/8)+ new_gId]
		in1 = datain[(1*N/8)+ new_gId]
		in2 = datain[(2*N/8)+ new_gId]
		in3 = datain[(3*N/8)+ new_gId]
		in4 = datain[(4*N/8)+ new_gId]
		in5 = datain[(5*N/8)+ new_gId]
		in6 = datain[(6*N/8)+ new_gId]
		in7 = datain[(7*N/8)+ new_gId]

		var angle : double = -2*PI*(new_gId%Ns)/(Ns*8)

		twidle_factor(2, angle, [&double](in2))
		twidle_factor(3, angle, [&double](in3))
		twidle_factor(4, angle, [&double](in4))
		twidle_factor(5, angle, [&double](in5))
		twidle_factor(6, angle, [&double](in6))
		twidle_factor(7, angle, [&double](in7))

		-- FFT8(in0, in1, in2, in3, in4, in5, in6, in7)

		var Idout : int = internal_block_id*128+(tidx)*8

		sharedx[Idout+0] = in0[0]		sharedy[Idout+0 + 1] = in0[1]
		sharedx[Idout+1] = in1[0]		sharedy[Idout+1 + 1] = in1[1]
		sharedx[Idout+2] = in2[0]		sharedy[Idout+2 + 1] = in2[1]
		sharedx[Idout+3] = in3[0]		sharedy[Idout+3 + 1] = in3[1]
		sharedx[Idout+4] = in4[0]		sharedy[Idout+4 + 1] = in4[1]
		sharedx[Idout+5] = in5[0]		sharedy[Idout+5 + 1] = in5[1]
		sharedx[Idout+6] = in6[0]		sharedy[Idout+6 + 1] = in6[1]
		sharedx[Idout+7] = in7[0]		sharedy[Idout+7 + 1] = in7[1]
		-- wait all threads finish here
		-- barrier(CLK_LOCAL_MEM_FENCE)
	end

	do -- 8-point kernel
		Ns = Ns * 8
		var new_block_id = 0 -- make map_id function
		-- var new_block_id : int = map_id(block_id, Local_Stride,  internal_blocksize)
		var new_gId : int = new_block_id*BLOCK_SIZE+tidx
		var Idin : int = (internal_block_id*16+tidx)
		in0[0] = sharedx[0*Local_Size+Idin]
		in0[1] = sharedy[0*Local_Size+Idin]
		in1[0] = sharedx[1*Local_Size+Idin]
		in1[1] = sharedy[1*Local_Size+Idin]
		in2[0] = sharedx[2*Local_Size+Idin]
		in2[1] = sharedy[2*Local_Size+Idin]
		in3[0] = sharedx[3*Local_Size+Idin]
		in3[1] = sharedy[3*Local_Size+Idin]
		in4[0] = sharedx[4*Local_Size+Idin]
		in4[1] = sharedy[4*Local_Size+Idin]
		in5[0] = sharedx[5*Local_Size+Idin]
		in5[1] = sharedy[5*Local_Size+Idin]
		in6[0] = sharedx[6*Local_Size+Idin]
		in6[1] = sharedy[6*Local_Size+Idin]
		in7[0] = sharedx[7*Local_Size+Idin]
		in7[1] = sharedy[7*Local_Size+Idin]

		var angle : double = -2*PI*(new_gId%Ns)/(Ns*8)

		twidle_factor(1, angle, [&double](in1))
		twidle_factor(2, angle, [&double](in2))
		twidle_factor(3, angle, [&double](in3))
		twidle_factor(4, angle, [&double](in4))
		twidle_factor(5, angle, [&double](in5))
		twidle_factor(6, angle, [&double](in6))
		twidle_factor(7, angle, [&double](in7))

		-- FFT8(in0, in1, in2, in3, in4, in5, in6, in7)

		var Idout : int  = 0
		-- var Idout : int = get_output_Id(new_gId, Ns, 8) 

		dataout[(0*Ns)+Idout] = in0[0]	dataout[(0*Ns)+Idout + 1] = in0[0]
		dataout[(1*Ns)+Idout] = in1[1]	dataout[(1*Ns)+Idout + 1] = in1[1]
		dataout[(2*Ns)+Idout] = in2[2]	dataout[(2*Ns)+Idout + 1] = in2[2]
		dataout[(3*Ns)+Idout] = in3[3]	dataout[(3*Ns)+Idout + 1] = in3[3]
		dataout[(4*Ns)+Idout] = in4[4]	dataout[(4*Ns)+Idout + 1] = in4[4]
		dataout[(5*Ns)+Idout] = in5[5]	dataout[(5*Ns)+Idout + 1] = in5[5]
		dataout[(6*Ns)+Idout] = in6[6]	dataout[(6*Ns)+Idout + 1] = in6[6]
		dataout[(7*Ns)+Idout] = in7[7]	dataout[(7*Ns)+Idout + 1] = in7[7]
	end
end

terra FFT8(in0 : &double, in1 : &double, in2 : &double, in3 : &double, in4 : &double, 
	in5 : &double, in6 : &double, in7 : &double)
	
	var v0 : double
	-- a,b
	v0 = in0[0]
	in0[0] = v0 + in1[0]
	in1[0] = v0 - in1[0]

	v0 = in0[1] -- same for imaginary
	in0[1] = v0 + in1[1]
	in1[1] = v0 - in1[1]

end

terra FFT_8(id : int, datain : &double, dataout : &double, Ns : int, N : int, K_W : int) 
	var gId = id
	var gSize = N/8
	var stride = gSize + 1
	var in0 : double[2]
	var in1 : double[2]
	var in2 : double[2]
	var in3 : double[2]
	var in4 : double[2]
	var in5 : double[2]
	var in6 : double[2]
	var in7 : double[2]
	
	in0[0] = datain[(0*stride)+gId] 	in0[1] = datain[(0*stride)+gId + 1]
	in1[0] = datain[(1*stride)+gId]		in1[1] = datain[(1*stride)+gId + 1]
	in2[0] = datain[(2*stride)+gId] 	in2[1] = datain[(2*stride)+gId + 1]
	in3[0] = datain[(3*stride)+gId]		in3[1] = datain[(3*stride)+gId + 1]
	in4[0] = datain[(4*stride)+gId] 	in4[1] = datain[(4*stride)+gId + 1]
	in5[0] = datain[(5*stride)+gId]		in5[1] = datain[(5*stride)+gId + 1]
	in6[0] = datain[(6*stride)+gId] 	in6[1] = datain[(6*stride)+gId + 1]
	in7[0] = datain[(7*stride)+gId]		in7[1] = datain[(7*stride)+gId + 1]

	-- at first stage twiddle factors is always 0
	if Ns ~=1 then
		var angle : double = -2*PI*(gId)/(N)
		cstdio.printf("angle: %f\n",angle)
		twidle_factor(1, angle, [&double](in1))
		twidle_factor(2, angle, [&double](in2))
		twidle_factor(3, angle, [&double](in3))
		twidle_factor(4, angle, [&double](in4))
		twidle_factor(5, angle, [&double](in5))
		twidle_factor(6, angle, [&double](in6))
		twidle_factor(7, angle, [&double](in7))
	end

	-- cstdio.printf("in%d: %f %f\n",0,in0[0],in0[1])
	-- cstdio.printf("in%d: %f %f\n",1,in1[0],in1[1])

	-- FFT8([&double](&in0),[&double](in1), [&double](in2), [&double](in3), 
		-- [&double](in4), [&double](in5), [&double](in6), [&double](in7), [&double](in8))

	var Idout : int = (gId/Ns)*Ns*K_W+(gId%Ns)
	dataout[(0*(Ns+1))+Idout] = in0[0]		dataout[(0*(Ns+1))+Idout + 1] = in0[1]
	dataout[(1*(Ns+1))+Idout] = in1[0]		dataout[(1*(Ns+1))+Idout + 1] = in1[1]
	dataout[(2*(Ns+1))+Idout] = in2[0] 		dataout[(2*(Ns+1))+Idout + 1] = in2[1]
	dataout[(3*(Ns+1))+Idout] = in3[0] 		dataout[(3*(Ns+1))+Idout + 1] = in3[1]
	dataout[(4*(Ns+1))+Idout] = in4[0]		dataout[(4*(Ns+1))+Idout + 1] = in4[1]
	dataout[(5*(Ns+1))+Idout] = in5[0]		dataout[(5*(Ns+1))+Idout + 1] = in5[1]
	dataout[(6*(Ns+1))+Idout] = in6[0] 		dataout[(6*(Ns+1))+Idout + 1] = in6[1]
	dataout[(7*(Ns+1))+Idout] = in7[0] 		dataout[(7*(Ns+1))+Idout + 1] = in7[1]
end

-- terra FFT_4(id : int, datain : &double, dataout : &double, Ns : int, N : int, K_W : int)
-- 	var gId = id
-- 	var gSize = N/4
-- 	var stride = gSize + 1
-- 	var in0 : double[2], in1 : double[2], in2 : double[2], in3 : double[2]

-- 	in0[0] = datain[(0*stride)+gId] 	in0[1] = datain[(0*stride)+gId + 1]
-- 	in1[0] = datain[(1*stride)+gId]		in1[1] = datain[(1*stride)+gId + 1]
-- 	in2[0] = datain[(2*stride)+gId] 	in2[1] = datain[(2*stride)+gId + 1]
-- 	in3[0] = datain[(3*stride)+gId]		in3[1] = datain[(3*stride)+gId + 1]

-- 	-- at first stage twiddle factors is always 0
-- 	if Ns ~=1 then
-- 		var angle : double = -2*PI*(gId)/(N)
-- 		cstdio.printf("angle: %f\n",angle)
-- 		twidle_factor(1, angle, in1)
-- 		twidle_factor(2, angle, in2)
-- 		twidle_factor(3, angle, in3)
-- 	end

-- 	FFT4([&double](&in0),[&double](in1), [&double](in2), [&double](in3))

-- 	-- cstdio.printf("in0: %f %f\n", in0[0], in0[1])
-- 	-- cstdio.printf("in1: %f %f\n", in1[0], in1[1])
-- 	-- cstdio.printf("in2: %f %f\n", in2[0], in2[1])
-- 	-- cstdio.printf("in3: %f %f\n", in3[0], in3[1])

-- 	var Idout : int = (gId/Ns)*Ns*K_W+(gId%Ns)
-- 	dataout[(0*(Ns+1))+Idout] = in0[0]		dataout[(0*(Ns+1))+Idout + 1] = in0[1]
-- 	dataout[(1*(Ns+1))+Idout] = in1[0]		dataout[(1*(Ns+1))+Idout + 1] = in1[1]
-- 	dataout[(2*(Ns+1))+Idout] = in2[0] 		dataout[(2*(Ns+1))+Idout + 1] = in2[1]
-- 	dataout[(3*(Ns+1))+Idout] = in3[0] 		dataout[(3*(Ns+1))+Idout + 1] = in3[1]
-- end

terra transposeComplexMatrix(M : int, N : int, data : &double)
	var lda : int = 2*N
	var temp : double
	var s : int
	for i=0,M do
		s = i+1
		for j = 2*i + 2, lda, 2 do
			temp = data[i*lda + j]
			data[i*lda + j] = data[s*lda + 2*i]
			data[s*lda + 2*i] = temp
			s = s + 1
		end
	end
end

terra FFT_4(start : int, datain : &double, Ns : int, NFFT : int) 
	var gId : int = start/4
	var s : int = start
	var stride = 2*Ns
	var in0 : double[2], in1 : double[2], in2 : double[2], in3 : double[2]

	in0[0] = datain[(0*stride)+s] 		in0[1] = datain[(0*stride)+s + 1]
	in1[0] = datain[(1*stride)+s]		in1[1] = datain[(1*stride)+s + 1]
	in2[0] = datain[(2*stride)+s] 		in2[1] = datain[(2*stride)+s + 1]
	in3[0] = datain[(3*stride)+s]		in3[1] = datain[(3*stride)+s + 1]

	-- at first stage twiddle factors is always 0
	if Ns ~=1 then
		var angle : double = -2*PI*(s)/(NFFT)
		cstdio.printf("angle: %f\n",angle)
		twidle_factor(1, angle, [&double](in1))
		twidle_factor(2, angle, [&double](in2))
		twidle_factor(3, angle, [&double](in3))
	end

	FFT4([&double](&in0),[&double](in1), [&double](in2), [&double](in3))
	datain[(0*stride)+s] = in0[0]	  datain[(0*stride)+s + 1] = in0[1]
	datain[(1*stride)+s] = in1[0]	  datain[(1*stride)+s + 1] = in1[1]
	datain[(2*stride)+s] = in2[0]	  datain[(2*stride)+s + 1] = in2[1]
	datain[(3*stride)+s] = in3[0]	  datain[(3*stride)+s + 1] = in3[1]
end

terra printFFT(str:rawstring,NFFT: int, A: &double)
	cstdio.printf("%s\n",str)
	for i=0,NFFT do 	
		cstdio.printf("%f %f\n",A[2*i],A[2*i+1])	
	end
end

terra printFFT(str:rawstring,NFFT: int, A: &double, base : int)
	cstdio.printf("%s\n",str)
	for i=0,NFFT do 	
		cstdio.printf("%f %f\n",A[base + 2*i],A[base + 2*i+1])	
	end
end

terra printComplexMatrix(str : rawstring, M : int, N : int, A : &double)
	for i=0,M*N do
		if i ~= 0 and i % N == 0 then
			cstdio.printf("\n")
		end
		cstdio.printf("%f %f ",A[2*i],A[2*i + 1])
	end
end

function genfftlua(NFFT) -- NFFT is equal to NB*NB
	local A = symbol("A")
	local stotal, rest, s, NELEMS, k, Ns = symbol("stotal"), symbol("rest"), symbol("s"), symbol("NELEMS"), symbol("k"), symbol("Ns")
	local K_W, skernel, base, ublocks = symbol("K_W"), symbol("skernel"), symbol("base"), symbol("ublocks")
	local ker = symmat("ker",2,2)
	local exec, bitreversal = terralib.newlist(), terralib.newlist()
	
	local ker = {}
	ker[0], ker[1] = {}, {}
	ker[0][0] = 4
	ker[0][1] = 2
	ker[1][0] = 2
	ker[1][1] = 1
 	
 	stotal = math.floor(math.log(NFFT)/math.log(2)) -- log2(NFFT)
	rest = stotal
	s = 0 -- current stage
	NELEMS = NFFT
	k = 0
	Ns = 1
	
	repeat
		-- #stages that type of kernel absorb is lower or equal lacking #stages
		if ker[k][1] <= rest then 
			K_W = ker[k][0]
			skernel = ker[k][1]
			base  = 0
			ublocks = NELEMS/K_W
			rest = rest - skernel
			-- cstdio.printf("#BLOCKS: %d  K_W: %d\n",ublocks,K_W)
			-- cstdio.printf("#INNER ITER: %d\n",Ns)
			for i=0,ublocks-1 do -- loop over big ublocks
				base = i*(Ns*K_W*2) -- iterate over the big ublocks, K_W*2 each elem size
				for j=0, Ns-1 do -- inside each block
					-- cstdio.printf("BASE: %d Ns: %d\n",base + 2*j,Ns)
					if k == 0 then
						exec:insert(quote
							FFT_4(base + 2*j,A,Ns,NFFT)
						end)
					elseif k == 1 then
						exec:insert(quote
							FFT_2(i*Ns+j,base + 2*j,A,Ns,NFFT)
						end)
					end
				end
			end
			-- jump stages	(min. skernel == 1 and K_W = 2)
			s = s + skernel
			Ns = Ns * K_W
			NELEMS = NELEMS / K_W 
		else
			k = k + 1
		end
	until rest == 0
	
	-- permute
	bitreversal:insert(quote
		cbit.reversal(NFFT,A)
	end)

	return terra([A] : &double, [lda] : int)
		[exec];
		[bitreversal];
	end
end

-- [[ NFFT: FFT #points,  A: array with the NFFT points, C: output, ker: possible kernels, ftype: 1 to FFT, -1 to IFFT ]]
terra genfft(A : &double, NFFT : int, ftype : int)
	
	-- necessary total stages
	var ker : int[2][2]
	ker[0][0] = 4	ker[0][1] = 2 -- size 4
	ker[1][0] = 2 	ker[1][1] = 1 -- size 2 

	var stotal : int = [int](cmath.log2(NFFT))
	var rest : int = stotal
	var s : int = 0 -- current stage
	var NELEMS : int = NFFT
	var k : int = 0
	var Ns : int = 1
	var K_W : int, skernel : int, base : int, ublocks : int, j : int
	repeat
		-- #stages that type of kernel absorb is lower or equal lacking #stages
		if ker[k][1] <= rest then 
			K_W = ker[k][0]
			skernel = ker[k][1]
			base  = 0
			ublocks = NELEMS/K_W
			rest = rest - skernel
			-- cstdio.printf("#BLOCKS: %d  K_W: %d\n",ublocks,K_W)
			-- cstdio.printf("#INNER ITER: %d\n",Ns)
			for i=0,ublocks do -- loop over big ublocks
				base = i*(Ns*K_W*2) -- iterate over the big ublocks, K_W*2 each elem size
				for j=0, Ns do -- inside each block
					-- cstdio.printf("BASE: %d Ns: %d\n",base + 2*j,Ns)
					if k == 0 then
						FFT_4(base + 2*j,A,Ns,NFFT)
					elseif k == 1 then
						FFT_2(i*Ns+j,base + 2*j,A,Ns,NFFT)
					end
				end
			end
			-- jump stages	(min. skernel == 1 and K_W = 2)
			s = s + skernel
			Ns = Ns * K_W

			-- necessary to assembly big blocks
			-- considered big blocks change over the iteration
			NELEMS = NELEMS / K_W 
			-- printFFT("Temp. Frequency",NFFT,A)
		else
			k = k + 1
		end
	until rest == 0
	
	-- permute
	cbit.reversal(NFFT,A)
end

terra testfft4()
	var NFFT = 4 -- number of points
	var Ns = 1
	var K_W = 4
	var x = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
    x[0] = 5 x[1] = 4 x[2] = 5 x[3] = 4

  var data = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
  for i=0,NFFT do
     data[2*i] = x[i]
     data[2*i+1] = 0.0
  end

	cstdio.printf("Time: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end

	FFT_4(0,data,Ns,NFFT)

	cbit.reversal(4,data)

	cstdio.printf("Frequency: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end
end

terra testfft2()
	var NFFT = 2 -- number of points
	var Ns = 1
	var K_W = 2
	
   var data = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))

	-- -- example
	data[0] = 5.0  data[1] = 0.0 -- complex 1
	data[2] = 4.0  data[3] = 0.0 -- complex 2

	cstdio.printf("Time: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end

	FFT_2(0,data,Ns,NFFT)

	cstdio.printf("Frequency: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end
end

terra testfft2()
	var NFFT = 2 -- number of points
	var Ns = 1
	var K_W = 2
	
	var x = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
   x[0] = 5 x[1] = 4 x[2] = 5 x[3] = 4

  var data = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
  for i=0,NFFT do
     data[2*i] = x[i]
     data[2*i+1] = 0.0
  end

	-- -- example
	-- data[0] = 5.0  data[1] = 0.0 -- complex 1
	-- data[2] = 4.0  data[3] = 0.0 -- complex 2

	cstdio.printf("Time: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end

	FFT_2(0,data,Ns,NFFT)

	cstdio.printf("Frequency: \n")
	for i=0,NFFT do
		cstdio.printf("%f %f\n",data[2*i],data[2*i+1])
	end
 end

 terra test1(Z : &double, ttype : int)
	 if ttype == 0 then -- image
		Z[1] = 5  Z[3] = 4  Z[5] = 5  Z[7] = 4
		Z[9] = 3  Z[11] = 2 Z[13] = 3 Z[15] = 2
		Z[17] = 5 Z[19] = 4 Z[21] = 5 Z[23] = 4
		Z[25] = 3 Z[27] = 2 Z[29] = 3 Z[31] = 2
	 else -- kernel
		Z[1] = 0.707 Z[3] = 0.707
	 end
 end

 terra test21(Z : &double, ttype : int)
	 if ttype == 0 then -- image
		Z[1] = 1  Z[3] = 2  Z[5] = 3
		Z[7] = 4  Z[9] = 5 Z[11] = 6
		Z[13] = 7 Z[15] = 8 Z[17] = 9
	 else -- kernel
		Z[1] = -1 Z[3] = -2 Z[5] = -1 
		Z[7] = 0  Z[9] = 0 Z[11] = 0 
		Z[13] = 1 Z[15] = 2 Z[17] = 1 
	 end
 end

 terra test22(Z : &double, ttype : int)
	 if ttype == 0 then -- image
		Z[1] = 1  Z[3] = 2  Z[5] = 3  Z[7] = 0
		Z[9] = 4  Z[11] = 5 Z[13] = 6 Z[15] = 0
		Z[17] = 7 Z[19] = 8 Z[21] = 9 Z[23] = 0
		Z[25] = 0 Z[27] = 0 Z[29] = 0 Z[31] = 0
	 else -- kernel
		Z[1] = -1 Z[3] = -2 Z[5] = -1 Z[7] = 0
		Z[9] = 0  Z[11] = 0 Z[13] = 0 Z[15] = 0
		Z[17] = 1 Z[19] = 2 Z[21] = 1 Z[23] = 0
		Z[25] = 0 Z[27] = 0 Z[29] = 0 Z[31] = 0
	 end
 end

 terra test3(Z : &double)
	 Z[1] = 0.5536732287904482
	 Z[3] = 0.5896971122146935
	 Z[5] = -0.770313607771324
	 Z[7] = 0.2443430805431277
end

terra test4(Z : &double,ttype : int)
	if ttype == 0 then 
		Z[1] = 1 Z[3] = 2 Z[5] = 5 Z[7] = 0
		Z[9] = 6 Z[11] = 7 Z[13] = 11 Z[15] = 0
	else
		Z[1] = 2 Z[3] = 5 Z[5] = 9 Z[7] = 0
		Z[9] = 1 Z[11] = 2 Z[13] = 8 Z[15] = 0
	end
end

terra fft(Nx : int)
	var imsize : int 
	var NFFT : int
	var x : &double
	var X : &double

	cstdio.printf("Nx = %d\n", Nx)
	x = [&double](cstdlib.calloc(Nx, sizeof(double)))
-- for i=0,Nx do
	  x[0] =	5 x[1] = 4 x[2] = 5 x[3] = 4 -- 4 complex numbers
	  x[4] = 3 x[5] = 2 x[6] = 3 x[7] = 2 -- 4 complex numbers
	  -- x[8] =	5 x[9] = 4 x[10] = 5 x[11] = 4 -- 4 complex numbers
	  -- x[12] = 3 x[13] = 2 x[14] = 3 x[15] = 2 -- 4 complex numbers
	  -- x[16] = 5 x[17] = 4 x[18] = 5 x[19] = 4 -- 4 complex numbers
	  -- x[20] = 3 x[21] = 2 x[22] = 3 x[23] = 2 -- 4 complex numbers
	  -- x[24] = 5 x[25] = 4 x[26] = 5 x[27] = 4 -- 4 complex numbers
	  -- x[28] = 3 x[29] = 2 x[30] = 3 x[31] = 2 -- 4 complex numbers
-- end

   	 -- calculate NFFT as the next higher power of 2 >= Nx
	NFFT = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](Nx))/cmath.log(2.0))))

	cstdio.printf("NFFT = %d\n", NFFT)

	  -- allocate memory for NFFT complex numbers (note the +1) 
	X = [&double](cstdlib.calloc(2*NFFT+1, sizeof(double)))

  -- Storing x(n) in a complex array to make it work with four1. 
  -- This is needed even though x(n) is purely real in this case. 
  for i=0,Nx do
     X[2*i+1] = x[i]
     X[2*i+2] = 0.0
  end

  -- test1(X,0)
  
  -- padding the remainder of the array with zeros (0 + 0 j) 
  for i=Nx, NFFT do
     X[2*i+1] = 0.0
     X[2*i+2] = 0.0
  end

  
  cstdio.printf("\nInput complex sequence (padded to next highest power of 2):\n")
  for i=0, NFFT do
	cstdio.printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end
  -- calculate FFT 
  four1(X, NFFT, 1)
  
  cstdio.printf("\nFFT:\n")
  for i=0, NFFT do
	cstdio.printf("X[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end

  -- calculate IFFT 
  four1(X, NFFT, -1)

  -- normalize the IFFT 
  for i=0, NFFT do
	X[2*i+1] = X[2*i+1] / NFFT
	X[2*i+2] = X[2*i+2] / NFFT
  end

  cstdio.printf("\nComplex sequence reconstructed by IFFT:\n")
  for i=0, NFFT do
	cstdio.printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end

  cstdio.getchar()
end

terra fast1Dconv(imagesize : int, kersize : int)
  -- assume kernel will be always smaller than image (partical approach)
  var image : &double, imagesize : int 
  var ker: &double, kersize : int
  var X : &double, NFFT : int
  var Y : &double

  -- generating kernel
  cstdio.printf("kersize = %d\n", kersize)
  ker = ([&double](cstdlib.calloc(kersize, sizeof(double))))
  -- for i=0,kersize do ker[i] = i end


  -- generating image
  cstdio.printf("imagesize = %d\n", imagesize)
  image = [&double](cstdlib.calloc(imagesize, sizeof(double)))
  -- for i=0,imagesize do image[i] = i end

  -- calculate NFFT as the next higher power of 2 >= imagesize
  NFFT = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](2*imagesize))/cmath.log(2.0))))
  cstdio.printf("NFFT = %d\n", NFFT)
  -- allocate memory for NFFT complex numbers (note the +1)  --?? +1
  X = [&double](cstdlib.calloc(2*NFFT+1, sizeof(double)))
  Y = [&double](cstdlib.calloc(2*NFFT+1, sizeof(double)))

  -- for i=0,imagesize do  X[2*i+1] = image[i] X[2*i+2] = 0.0 end
  -- pad the remainder of the array with zeros (0 + 0 j) 
  -- for i=imagesize, NFFT do X[2*i+1] = 0.0 X[2*i+2] = 0.0 end

  -- for i=0,kersize do Y[2*i+1] = ker[i] Y[2*i+2] = 0.0 end
  -- for i=kersize, NFFT do Y[2*i+1] = 0.0 Y[2*i+2] = 0.0 end
  test1(X,0)
  test1(Y,1)

  cstdio.printf("\nInput complex sequence (padded to next highest power of 2):\n")
  cstdio.printf("Image: ")
  for i=0, NFFT do
	cstdio.printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end

  cstdio.printf("Kernel: ")
  for i=0, NFFT do
	cstdio.printf("y[%d] = (%.2f + j %.2f)\n", i, Y[2*i+1], Y[2*i+2])
  end

  -- FFT image
  four1(X, NFFT, 1)
  -- FFT kernel
  four1(Y, NFFT, 1)

  --point-wise multiply, dot product. Can be optimized with vector inst and prefetch
  for i=0,NFFT do
	var real = X[2*i+1] * Y[2*i+1] - X[2*i+2] * Y[2*i+2]
	var imag = X[2*i+1] * Y[2*i+2] + X[2*i+2] * Y[2*i+1]
	X[2*i+1] = real
	X[2*i+2] = imag
  end

  cstdio.printf("\nAfter dot product:\n")
  for i=0, NFFT do
	cstdio.printf("X[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end

  -- calculate IFFT 
  four1(X, NFFT, -1)

  -- normalize the IFFT 
  for i=0, NFFT do
	X[2*i+1] = X[2*i+1] / NFFT
	X[2*i+2] = X[2*i+2] / NFFT
  end

  cstdio.printf("\n1D convolution:\n")
  for i=0, NFFT do
	cstdio.printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i+1], X[2*i+2])
  end

  cstdio.getchar()
end

-- [[ outplace from X to X ]]
-- terra transpose(X : &double, X : &double, N : int, M : int)
-- 	var count = 1
-- 	for i=0,N*2,2 do
-- 		for j=0,M do
-- 			X[count] = X[j*(M*2) + i + 1]
-- 			X[count + 1] = X[j*(M*2) + i + 2]
-- 			count = count + 2
-- 		end
-- 	end
-- end

-- [[ kernel based fft 2D. Input: X with dim MxN. type: 1 to forward, -1 to backward ]]
terra offt2D(M : int, N : int, X : &double, type: int)
	for j=0,M do
		var base = j*N*2
		printFFT("Row",N,X,base)

		-- calculate FFT
		genfft(&X[base], N, type)

		-- printFFT("FFT",N,X,base)
	end

	 -- computing transpose
	 -- doing for squared matrices
	 cstdio.printf("\nTransposed matrix\n")
	 transposeComplexMatrix(M,N,X)

	 -- same process
	 for j=0,M do
		var base = j*N*2

		printFFT("Row",N,X,base)

		-- calculate FFT
		genfft(&X[base], N, type)
	
		-- printFFT("FFT",N,X,base)
	end

	cstdio.printf("\nTransposing it back\n") -- transposeComplexMatrix it back
	transposeComplexMatrix(M,N,X)
	
	cstdio.printf("\n")
end

terra fft2D(M : int, N : int, X : &double, ftype : int)

	for j=0,M do
		var base = j*N*2
		printFFT("Row",N,X,base+1)
		
		-- calculate FFT
		four1(&X[base], N, ftype)

		printFFT("FFT",N,X,base+1)
	end

	 -- computing transposeComplexMatrix
	 -- doing for squared matrices
	 cstdio.printf("\nTransposed matrix\n")
	 transposeComplexMatrix(M,N,X)

	 -- same process
	 for j=0,M do
		var base = j*N*2

		printFFT("Row",N,X,base+1)

		-- calculate FFT
		four1(&X[base], N, ftype)
	
		printFFT("FFT",N,X,base+1)
	end

	cstdio.printf("\nTransposing it back\n") -- transpose it back
	transposeComplexMatrix(M,N,X)
	
	cstdio.printf("\n")
end

-- [[ Inplace pointwise multiplication of X and Y. Output: X ]]
terra cmult(X : &double, Y : &double, NFFT : int, s : int)
	
  -- Can be optimized with vector inst and prefetch
  -- s: start point of X, usually 0 or 1
	for i=0,NFFT do
		var real = X[2*i+s] * Y[2*i+s] - X[2*i+s+1] * Y[2*i+s+1]
		var imag = X[2*i+s] * Y[2*i+s+1] + X[2*i+s+1] * Y[2*i+s]
		X[2*i+s] = real
		X[2*i+s+1] = imag
	end
end

function genFastConvolution(NB,NBF,RM,RN,V)
	local NB2 = NB * NBF																																			
	local l1cmul = gencmulkernel(NB, RM, RN, V, false)

	return terra(gettime : {} -> double, M : int, N : int, A : &double, B : &double, C : &double, 
		ldc : int) 

         [ blockedloop(N,M,{NB2,NB},
                function(m,n) 
                return quote
                    var MM,NN = min(M-m,NB),min(N-n,NB)
                    var isboundary = MM < NB or NN < NB
                    var AA,BB,CC = A + (m*ldc + n), B + (m*ldc + n), C + (m*ldc + n)


                    l1cmul(AA,BB,CC,ldc)


                end end)  
            ]       
	end
end

-- [[ Assuming image power of two, kernel and image blocks with same size ]]
terra genC(A : &double, B : &double, NFFTx : int, NFFTy : int)
	var NFFT = NFFTx*NFFTy

	-- image 2D fft 
	offt2D(NFFTx, NFFTy, A, 1)

	-- kernel 2D fft
	offt2D(NFFTx, NFFTy, B, 1)

	-- use the kernel
	cmult(A, B, NFFT,0)

	-- calculate IFFT 
	offt2D(NFFTx, NFFTy, A, -1)

	-- normalize IFFT
	for i=0, NFFT do
		A[2*i+1] = A[2*i+1] / NFFT
		A[2*i+2] = A[2*i+2] / NFFT
	end
	printFFT("2D convolution",NFFTx*NFFTy,A,1)
end

--[[ Fast (spatial) Convolution. Image: MxN,  Kernel: KxL]]
terra fast2Dconv(M : int, N: int, K : int, L : int)
	var x : &double, X : &double
	var y : &double, Y : &double

	-- squared in this case 
	-- exapand different dimensions
	var NFFTx = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](M))/cmath.log(2.0))))
	var NFFTy = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](N))/cmath.log(2.0))))
	var NFFT : int = NFFTx * NFFTy
	cstdio.printf("Image and kernel were expanded to %dx%d\n", NFFTx, NFFTy)

	-- +1: the vectors start from index 1
	x = ([&double](cstdlib.calloc(M*N, sizeof(double))))
	X = ([&double](cstdlib.calloc(2*NFFT+1, sizeof(double))))

	y = [&double](cstdlib.calloc(K*L, sizeof(double)))
	Y = ([&double](cstdlib.calloc(2*NFFT+1, sizeof(double))))

	test1(X,0)
	
	-- making up kernel
	-- for i=0,K*L do y[i] = i end
	-- for i=0,K*L do
	-- 	Y[2*i+1] = y[i]
	-- 	Y[2*i+2] = 0.0
	-- end

	test1(Y,1)

	-- making up image data
	-- for i=0,M*N do x[i] = i end
	-- for i=0,M*N do
	--    X[2*i+1] = x[i]
	--    X[2*i+2] = 0.0
	-- end

	printFFT("Image",NFFT,X,1)
	
	printFFT("Kernel",NFFT,Y,1)
	
	-- -- image 2D fft
	fft2D(NFFTx, NFFTy, X, 1)

	printFFT("Image convolved",NFFT,X,1)

	-- kernel 2D fft
	fft2D(NFFTx, NFFTy, Y, 1)

	printFFT("Kernel convolved",NFFT,Y,1)
	
  -- complex point-wise multiplication. Can be optimized with vector inst and prefetch
	cmult(X,Y,NFFT,1)

	printFFT("After dot product",NFFT,X,1)

	-- calculate IFFT 
	fft2D(NFFTx, NFFTy, X, -1)

	-- normalize IFF
	for i=0, NFFT do
	X[2*i+1] = X[2*i+1] / NFFT
	X[2*i+2] = X[2*i+2] / NFFT
	end

	printFFT("2D convolution",NFFT,X,1)

	cstdio.getchar()
end

function log2(number)
    return math.log(number) / math.log(2)
end

terra oFFTtest(NFFT : int)
	-- Creating data
	var x = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
	var y = [&double](cstdlib.calloc(2*NFFT+1, sizeof(double)))
	x[0] = 5  x[1] = 4 x[2] = 5 x[3] = 4
	x[4] = 3  x[5] = 2 x[6] = 3 x[7] = 2
	-- x[8] = 5  x[9] = 4 x[10] = 5 x[11] = 4
	-- x[12] = 3 x[13] = 2 x[14] = 3 x[15] = 2
	-- x[16] = 5  x[17] = 4 x[18] = 5 x[19] = 4
	-- x[20] = 3  x[21] = 2 x[22] = 3 x[23] = 2
	-- x[24] = 5  x[25] = 4 x[26] = 5 x[27] = 4
	-- x[28] = 3 x[29] = 2 x[30] = 3 x[31] = 2

	var data = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
	for i=0,NFFT do
		data[2*i] = x[i]
		data[2*i+1] = 0.0
	end

	printFFT("Time",NFFT,data)
	genfft(data,NFFT,1)
end

function testgenfftlua()
	local NB = 8
	local l1fft = genfftlua(NB)
end

--[[ Create a matrix of complexes ]]
terra newcmatrix(x : &double, NFFT: int, padding : bool) : &double
	var M : &double = [&double](cstdlib.calloc(2*NFFT, sizeof(double)))
	for i=0,NFFT do
	  M[2*i] = x[i]
	  -- M[2*i+1] = 0.0 --calloc
	end
	return M
end

terra min(a : int, b : int)
	return terralib.select(a < b, a, b)
end

terra testtranspose(M : int, N : int)
	var NFFT = M*N
	var x = [&double](cstdlib.calloc(NFFT, sizeof(double)))

	var sum = 1
	for i=0,M do
		for j=0,N do
			x[i*N + j] = sum
			sum = sum + 1
		end
	end

	var data = newcmatrix(x,NFFT,true)
	printComplexMatrix("Image",M,N,data)
	cstdio.printf("\n\n")
	var bb0 = 2
	var temp : double
	
	-- transpose
	transposeComplexMatrix(M,N,data)
	printComplexMatrix("Image",M,N,data)	
end

terra gen2Dconv(M : int, N: int, K: int, L: int)
	var x : &double, data : &double
	var y : &double, ker : &double

	-- power two treatment and kernel expansion: zero-padding used
	var NFFTx = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](M))/cmath.log(2.0))))
	var NFFTy = [int](cmath.pow(2.0, cmath.ceil(cmath.log([double](N))/cmath.log(2.0))))
	var NFFT = NFFTx * NFFTy
	cstdio.printf("Image and kernel expanded to size: %dx%d\n", NFFTx,NFFTy)

	-- generating user data
	x = [&double](cstdlib.calloc(M*N, sizeof(double)))
	x[0] = 5  x[1] = 4 
	x[2] = 5  x[3] = 4 	
	y = [&double](cstdlib.calloc(K*L, sizeof(double)))
	y[0] = 3  y[1] = 2
	y[2] = 3  y[3] = 2
	data = newcmatrix(x,NFFT,true)
	ker  = newcmatrix(y,NFFT,true)
	printFFT("Image",NFFT,data)
	printFFT("Kernel",NFFT,ker)

	-- calling fast convolsluton 2D
	genC(data, ker, M, N)
end

terra reversebits(NFFT : int)
	var x : &double, data : &double, c : &double
	x = [&double](cstdlib.calloc(NFFT, sizeof(double)))
	for i=0,NFFT do
		x[i] = i % 10
	end
	-- x[4] = 3  x[5] = 2 x[6] = 3  x[7] = 2
	data = newcmatrix(x,NFFT,true)
	c = [&double](cstdlib.calloc(2*NFFT+1, sizeof(double)))
	
	for i=0,NFFT do
		c[2*i+1] = data[2*i] 
		c[2*i+2] = data[2*i+1]
	end
	var n : int, nn = NFFT,NFFT
	var j : int, m : int, tempr : double
	cbit.reversalbased(&n,&nn,&j,&m,&tempr,c)
	cbit.reversal(NFFT,data)

	for i=0, NFFT do
		if c[2*i+1] ~= data[2*i] or c[2*i+2] ~= data[2*i + 1] then
			cstdio.printf("Error\n")
		end
	end	
end

terra testreversal()
	var NFFT : int
	var i = 2
	repeat
		NFFT = i
		reversebits(NFFT)
		i = i*2
	until i == 1024
end

	local M,N,K,L
	local NFFT
	-- fft(NFFT)
	-- fast1Dconv(M,N)
	-- M, N, K, L = 4,,1,2
	-- fast2Dconv(M,N,K,L)
	-- testfft2()
	-- testfft4()
	-- oFFTtest(NFFT)
	M, N, K, L = 2,2,2,2
	gen2Dconv(M,N,K,L)
	-- M, N = 4,4
	-- testtranspose(M,N)
	-- testreversal()

