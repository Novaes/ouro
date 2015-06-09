#include<stdio.h>
#include<assert.h>
#include<sys/time.h>

/*
    A is a matrix of dim M x K
    B is a matrix of dim K x N
    C is a matrix of dim M x N
    lda is the stride of matrix A (normally its K, but it can be larger for alignment)
    ldb is the stride of matrix B (normally its N, ...)
    ldc is the stride of the output matrix (normally it is N, ...)
*/
void asserteq(double * A, double * B, int M, int N) {
    int i = 0;
    for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
            if(A[i] != B[i]) {
                goto printerr;
            }
            i++;
        }
    }
    return;
printerr:
    i = 0;
    for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
            if(A[i] != B[i]) {
                putchar('X');
            } else {
                putchar('*');
            }
            i++;
        }
        printf("\n");
    }
    assert(false);
}

extern "C"
void my_dgemm(double (*gettime)(), const int M, const int N, const int K,
    const int L, const double alpha, const double *A, const int lda, const double *B,
    const int ldb, const double *C, const int ldc, const int kCenterX, const int kCenterY);

static double CurrentTimeInSeconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

bool CalcTime(int * times, double * start) {
    if(*times == 0){
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

void print_matrix(double* matrix, int M, int N){
    for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
            printf("%f ", matrix[m*N + n]);
        }
        printf("\n");
    }
    printf("\n");
}

//calculates C = alpha * AB + beta * C
void testsize(int M, int N, int K, int L) {
    double * A = (double*) malloc(sizeof(double) * M * N);
    double * C = (double*) malloc(sizeof(double) * M * N);
    double * C2 = (double*) malloc(sizeof(double) * M * N);
    double * B = (double*) malloc(sizeof(double) * K * L);
    
    for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
            if( m == 0  || n == 0 ||  m == M-1  || n == N-1)
                A[m*N + n] = 0;
            else 
                A[m*N + n] = (m-1)*(N-2) + n;
        }
    }

    // flipped kernel
    B[0] = 1; B[1] = 2; B[2] = 1;
    B[3] = 0; B[4] = 0; B[5] = 0;
    B[6] = -1;B[7] = -2;B[8] = -1;

     double mytime;
     int times = 0;

    while(CalcTime(&times,&mytime))
        my_dgemm(CurrentTimeInSeconds,M,N,K,L,1.f,A,N,B,L,C,N,K/2,L/2);
    
    // test
    C2[0] = -13; C2[1] = -20; C2[2] = -17;
    C2[3] = -18; C2[4] = -24; C2[5] = -18;
    C2[6] = 13;  C2[7] = 20;  C2[8] = 17;
    //asserteq(C,C2,M,N);

    for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
            printf("%f ", C[m*N + n]);
        }
        printf("\n");
    }

    free(C);
    free(C2);
    free(A);
    free(B);
    // //double logblastime = log(M) + log(N) + log(K) + log(2) + log(1e-9) - log(blastime);
    // double logmytime = log(M) + log(N) + log(K) + log(2) + log(1e-9) - log(mytime);
    // //printf("%d %d %d %f %f %f\n",M,K,N,exp(logblastime),exp(logmytime), mytime/ blastime);
    // printf("%d %d %d %f %f \n",M,K,N,exp(logmytime), mytime);
}

int main() {
    int NB = 40;
    //  for(int i = NB; true; i += 3*NB) {
        // int m = i;
        // int n = i;
        // int k = i;
        // int k = l;
        
        // if(m*n+ m*k+n*k > 3*2048*2048) //usually max img resolution
        //  testsize(m,n,k,l);//A and B are suppose to be at max by dgemm: 1024x1024 x 3, so 2 x 1024x1024
        //  break;}
    //testsize(5000,5000,5000);
    testsize(5,5,3,3);
    return 0;
}
