local cstdio = terralib.includec("stdio.h")
local cstdlib = terralib.includec("stdlib.h")
local cmath = terralib.includec("math.h")
local PI = constant(3.141592653589793)
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
	  for m = 1,mmax,2 do -- 1 to #stages; 2-point butterflies
		for i = m,n+1,istep do -- 2-point butterfly
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
terra FFT_2(id: int, start : int, datain : &double, Ns : int, NFFT : int, signal : int)
	var gId : int = id % (NFFT/2) -- gsize is NFFT/2
	var s : int = start
	var stride = 2*Ns -- +1 due the imaginary representation
	-- cstdio.printf("init: %d %d\n",s, Ns)

	var in0 : double[2], in1 : double[2]

	-- load reals and imaginaries
	in0[0] = datain[s] 	in0[1] = datain[s + 1]
	in1[0] = datain[(1*stride)+s]	in1[1] = datain[(1*stride)+s + 1]

--[[
	-- at first stage twiddle factors is always 0
	 if Ns ~=1 then
		 var angle : double = signal*-2*PI*gId/NFFT
		 -- cstdio.printf("angle: %f\n",angle)
		 twidle_factor(1, angle, [&double](in1))
	 end
]]
	-- cstdio.printf("in%d: %f %f\n",0,in0[0],in0[1])
	-- cstdio.printf("in%d: %f %f\n",1,in1[0],in1[1])

	FFT2([&double](&in0), [&double](&in1))

	if signal == -1 then 
		in0[0] = in0[0] / 2     
		in0[1] = in0[1] / 2
		in1[0] = in1[0]	/ 2     
		in1[1] = in1[1] / 2
	end


	-- cstdio.printf("--> %d %d\n",(0*stride)+s,(1*stride)+s)	
	datain[s] = in0[0]	  			  datain[s + 1] = in0[1]
	datain[(1*stride)+s] = in1[0]	  datain[(1*stride)+s + 1] = in1[1]
end

terra FFT_4(start : int, datain : &double, Ns : int, NFFT : int, signal : int)
	var gId : int = start/4 -- gsize is divided by 4
	var s : int = start
	var stride = 2*Ns
	var in0 : double[2], in1 : double[2], in2 : double[2], in3 : double[2]

	-- cstdio.printf("POSITIONS: %d %d %d %d",(0*stride)+s,(1*stride)+s,(2*stride)+s,(3*stride)+s)
	in0[0] = datain[(0*stride)+s] 		in0[1] = datain[(0*stride)+s + 1]
	in1[0] = datain[(1*stride)+s]		in1[1] = datain[(1*stride)+s + 1]
	in2[0] = datain[(2*stride)+s] 		in2[1] = datain[(2*stride)+s + 1]
	in3[0] = datain[(3*stride)+s]		in3[1] = datain[(3*stride)+s + 1]

	-- at first stage twiddle factors is always 0
--[[
	if Ns ~=1 then
		var angle : double = signal*-2*PI*s / NFFT
		twidle_factor(1, angle, [&double](in1))
		twidle_factor(2, angle, [&double](in2))
		twidle_factor(3, angle, [&double](in3))
	end
]]
	FFT4([&double](&in0),[&double](in1), [&double](in2), [&double](in3))

	if signal == -1 then 
		in0[0] = in0[0]/4
		in0[1] = in0[1]/4
		in1[0] = in1[0]/4
		in1[1] = in1[1]/4
		in2[0] = in2[0]/4
		in2[1] = in2[1]/4
		in3[0] = in3[0]/4
		in3[1] = in3[1]/4
	end

	datain[(0*stride)+s] = in0[0]	  datain[(0*stride)+s + 1] = in0[1]
	datain[(1*stride)+s] = in2[0]	  datain[(1*stride)+s + 1] = in2[1]
	datain[(2*stride)+s] = in1[0]	  datain[(2*stride)+s + 1] = in1[1]
	datain[(3*stride)+s] = in3[0]	  datain[(3*stride)+s + 1] = in3[1]
end



-- terra FFT_64(id : int, datain : &double, dataout : &double, Ns : int, N : int, K_W : int, Local_Stride : int, BLOCK_SIZE : int, lsize : int)
-- 	var sharedx : double[128*8]
-- 	var sharedy : double[128*8]
-- 	var in_stride : int = N / (128*8) 
-- 	var gId : int = id -- get_global_id(0)
-- 	-- get local size
-- 	var Local_Size : int = lsize

-- 	var in0 : double[2], in1 : double[2], in2 : double[2], in3 : double[2], in4 : double[2], in5 : double[2], in6 : double[2], in7 : double[2]

-- 	var internal_blocksize : int = 8
-- 	var block_id : int = gId/BLOCK_SIZE
-- 	var internal_block_id : int = block_id%internal_blocksize
-- 	-- var tidx : int = gId&(BLOCK_SIZE-1); -- binary and
-- 	var tidx : int = gId 

-- 	do -- 8-point kernel
-- 		var new_block_id = 0 -- make map_id function
-- 		-- var new_block_id : int = map_id(block_id, in_stride,  internal_blocksize);
-- 		var new_gId : int = new_block_id*BLOCK_SIZE+tidx;
-- 		in0 = datain[(0*N/8)+ new_gId]
-- 		in1 = datain[(1*N/8)+ new_gId]
-- 		in2 = datain[(2*N/8)+ new_gId]
-- 		in3 = datain[(3*N/8)+ new_gId]
-- 		in4 = datain[(4*N/8)+ new_gId]
-- 		in5 = datain[(5*N/8)+ new_gId]
-- 		in6 = datain[(6*N/8)+ new_gId]
-- 		in7 = datain[(7*N/8)+ new_gId]

-- 		var angle : double = -2*PI*(new_gId%Ns)/(Ns*8)

-- 		twidle_factor(2, angle, [&double](in2))
-- 		twidle_factor(3, angle, [&double](in3))
-- 		twidle_factor(4, angle, [&double](in4))
-- 		twidle_factor(5, angle, [&double](in5))
-- 		twidle_factor(6, angle, [&double](in6))
-- 		twidle_factor(7, angle, [&double](in7))

-- 		-- FFT8(in0, in1, in2, in3, in4, in5, in6, in7)

-- 		var Idout : int = internal_block_id*128+(tidx)*8

-- 		sharedx[Idout+0] = in0[0]		sharedy[Idout+0 + 1] = in0[1]
-- 		sharedx[Idout+1] = in1[0]		sharedy[Idout+1 + 1] = in1[1]
-- 		sharedx[Idout+2] = in2[0]		sharedy[Idout+2 + 1] = in2[1]
-- 		sharedx[Idout+3] = in3[0]		sharedy[Idout+3 + 1] = in3[1]
-- 		sharedx[Idout+4] = in4[0]		sharedy[Idout+4 + 1] = in4[1]
-- 		sharedx[Idout+5] = in5[0]		sharedy[Idout+5 + 1] = in5[1]
-- 		sharedx[Idout+6] = in6[0]		sharedy[Idout+6 + 1] = in6[1]
-- 		sharedx[Idout+7] = in7[0]		sharedy[Idout+7 + 1] = in7[1]
-- 		-- wait all threads finish here
-- 		-- barrier(CLK_LOCAL_MEM_FENCE)
-- 	end

-- 	do -- 8-point kernel
-- 		Ns = Ns * 8
-- 		var new_block_id = 0 -- make map_id function
-- 		-- var new_block_id : int = map_id(block_id, Local_Stride,  internal_blocksize)
-- 		var new_gId : int = new_block_id*BLOCK_SIZE+tidx
-- 		var Idin : int = (internal_block_id*16+tidx)
-- 		in0[0] = sharedx[0*Local_Size+Idin]
-- 		in0[1] = sharedy[0*Local_Size+Idin]
-- 		in1[0] = sharedx[1*Local_Size+Idin]
-- 		in1[1] = sharedy[1*Local_Size+Idin]
-- 		in2[0] = sharedx[2*Local_Size+Idin]
-- 		in2[1] = sharedy[2*Local_Size+Idin]
-- 		in3[0] = sharedx[3*Local_Size+Idin]
-- 		in3[1] = sharedy[3*Local_Size+Idin]
-- 		in4[0] = sharedx[4*Local_Size+Idin]
-- 		in4[1] = sharedy[4*Local_Size+Idin]
-- 		in5[0] = sharedx[5*Local_Size+Idin]
-- 		in5[1] = sharedy[5*Local_Size+Idin]
-- 		in6[0] = sharedx[6*Local_Size+Idin]
-- 		in6[1] = sharedy[6*Local_Size+Idin]
-- 		in7[0] = sharedx[7*Local_Size+Idin]
-- 		in7[1] = sharedy[7*Local_Size+Idin]

-- 		var angle : double = -2*PI*(new_gId%Ns)/(Ns*8)

-- 		twidle_factor(1, angle, [&double](in1))
-- 		twidle_factor(2, angle, [&double](in2))
-- 		twidle_factor(3, angle, [&double](in3))
-- 		twidle_factor(4, angle, [&double](in4))
-- 		twidle_factor(5, angle, [&double](in5))
-- 		twidle_factor(6, angle, [&double](in6))
-- 		twidle_factor(7, angle, [&double](in7))

-- 		-- FFT8(in0, in1, in2, in3, in4, in5, in6, in7)

-- 		var Idout : int  = 0
-- 		-- var Idout : int = get_output_Id(new_gId, Ns, 8) 

-- 		dataout[(0*Ns)+Idout] = in0[0]	dataout[(0*Ns)+Idout + 1] = in0[0]
-- 		dataout[(1*Ns)+Idout] = in1[1]	dataout[(1*Ns)+Idout + 1] = in1[1]
-- 		dataout[(2*Ns)+Idout] = in2[2]	dataout[(2*Ns)+Idout + 1] = in2[2]
-- 		dataout[(3*Ns)+Idout] = in3[3]	dataout[(3*Ns)+Idout + 1] = in3[3]
-- 		dataout[(4*Ns)+Idout] = in4[4]	dataout[(4*Ns)+Idout + 1] = in4[4]
-- 		dataout[(5*Ns)+Idout] = in5[5]	dataout[(5*Ns)+Idout + 1] = in5[5]
-- 		dataout[(6*Ns)+Idout] = in6[6]	dataout[(6*Ns)+Idout + 1] = in6[6]
-- 		dataout[(7*Ns)+Idout] = in7[7]	dataout[(7*Ns)+Idout + 1] = in7[7]
-- 	end
-- end

-- terra FFT8(in0 : &double, in1 : &double, in2 : &double, in3 : &double, in4 : &double, 
-- 	in5 : &double, in6 : &double, in7 : &double)
	
-- 	var v0 : double
-- 	-- a,b
-- 	v0 = in0[0]
-- 	in0[0] = v0 + in1[0]
-- 	in1[0] = v0 - in1[0]

-- 	v0 = in0[1] -- same for imaginary
-- 	in0[1] = v0 + in1[1]
-- 	in1[1] = v0 - in1[1]

-- end

-- terra FFT_8(id : int, datain : &double, dataout : &double, Ns : int, N : int, K_W : int) 
-- 	var gId = id
-- 	var gSize = N/8
-- 	var stride = gSize + 1
-- 	var in0 : double[2]
-- 	var in1 : double[2]
-- 	var in2 : double[2]
-- 	var in3 : double[2]
-- 	var in4 : double[2]
-- 	var in5 : double[2]
-- 	var in6 : double[2]
-- 	var in7 : double[2]
	
-- 	in0[0] = datain[(0*stride)+gId] 	in0[1] = datain[(0*stride)+gId + 1]
-- 	in1[0] = datain[(1*stride)+gId]		in1[1] = datain[(1*stride)+gId + 1]
-- 	in2[0] = datain[(2*stride)+gId] 	in2[1] = datain[(2*stride)+gId + 1]
-- 	in3[0] = datain[(3*stride)+gId]		in3[1] = datain[(3*stride)+gId + 1]
-- 	in4[0] = datain[(4*stride)+gId] 	in4[1] = datain[(4*stride)+gId + 1]
-- 	in5[0] = datain[(5*stride)+gId]		in5[1] = datain[(5*stride)+gId + 1]
-- 	in6[0] = datain[(6*stride)+gId] 	in6[1] = datain[(6*stride)+gId + 1]
-- 	in7[0] = datain[(7*stride)+gId]		in7[1] = datain[(7*stride)+gId + 1]

-- 	-- at first stage twiddle factors is always 0
-- 	if Ns ~=1 then
-- 		var angle : double = -2*PI*(gId)/(N)
-- 		cstdio.printf("angle: %f\n",angle)
-- 		twidle_factor(1, angle, [&double](in1))
-- 		twidle_factor(2, angle, [&double](in2))
-- 		twidle_factor(3, angle, [&double](in3))
-- 		twidle_factor(4, angle, [&double](in4))
-- 		twidle_factor(5, angle, [&double](in5))
-- 		twidle_factor(6, angle, [&double](in6))
-- 		twidle_factor(7, angle, [&double](in7))
-- 	end

-- 	-- cstdio.printf("in%d: %f %f\n",0,in0[0],in0[1])
-- 	-- cstdio.printf("in%d: %f %f\n",1,in1[0],in1[1])

-- 	-- FFT8([&double](&in0),[&double](in1), [&double](in2), [&double](in3), 
-- 		-- [&double](in4), [&double](in5), [&double](in6), [&double](in7), [&double](in8))

-- 	var Idout : int = (gId/Ns)*Ns*K_W+(gId%Ns)
-- 	dataout[(0*(Ns+1))+Idout] = in0[0]		dataout[(0*(Ns+1))+Idout + 1] = in0[1]
-- 	dataout[(1*(Ns+1))+Idout] = in1[0]		dataout[(1*(Ns+1))+Idout + 1] = in1[1]
-- 	dataout[(2*(Ns+1))+Idout] = in2[0] 		dataout[(2*(Ns+1))+Idout + 1] = in2[1]
-- 	dataout[(3*(Ns+1))+Idout] = in3[0] 		dataout[(3*(Ns+1))+Idout + 1] = in3[1]
-- 	dataout[(4*(Ns+1))+Idout] = in4[0]		dataout[(4*(Ns+1))+Idout + 1] = in4[1]
-- 	dataout[(5*(Ns+1))+Idout] = in5[0]		dataout[(5*(Ns+1))+Idout + 1] = in5[1]
-- 	dataout[(6*(Ns+1))+Idout] = in6[0] 		dataout[(6*(Ns+1))+Idout + 1] = in6[1]
-- 	dataout[(7*(Ns+1))+Idout] = in7[0] 		dataout[(7*(Ns+1))+Idout + 1] = in7[1]
-- end

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
