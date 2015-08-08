#include <stdbool.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

static double CurrentTimeInSeconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

bool CalcTime(int * times, double * start) {
	if(*times == 0) {
		*start = CurrentTimeInSeconds();
	} else {
		double elapsed = CurrentTimeInSeconds() - *start;
		if(elapsed > 0.1f) {
			*start = elapsed / *times;
			return false;
		}
	}
	(*times)++;
	return true;
}

// call the convolution
//void my_convolution(double (*gettime)(),int argc, char* argv[]);

void testsize(int M, int N, int K) {
	double mytime;
	int times = 0;
	while(CalcTime(&times,&mytime)) {
		//my_convolution(CurrentTimeInSeconds,argc,argv);
	}
	double logmytime = log(M) + log(N) + log(K) + log(2) + log(1e-9) - log(mytime);
	printf("%d %d %d %f \n",M,K,N,logmytime);
}

int main() {
	int NB = 48;
	// for(int i = NB; i < 3000; i += 3*NB) {
	// 	testsize(i,i,i);
	// }
	testsize(256, 256, 5);

	return 0;
}