#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>
// #include "mkl_vsl.h"

// int assertintel(double * A, double * B, int M, int N){
//     int i,j;
//     for (j=0; j<M; j++)
//         for (i=0; i<N; i++) {
//             double zij = A[i + N*j];
//             double eij = B[i + N*j];
//             if (fabs(zij-eij) > fabs(eij)*1e-10) {
//                 printf("ERROR: wrong results:\n");
//                 printf("    A[%2d,%2d]: %lg\n",i,j,zij);
//                 printf("    expected: %lg\n",eij);
//                 printf("EXAMPLE FAILED\n");
//                 return 1;
//             }
//         }
//     return 0;
// }

// // detailed assert
// void detassert(double * A, double * B, int M, int N) {
//     int i = 0;
//     int m,n;
//     for(m = 0; m < M; m++) {
//         for(n = 0; n < N; n++) {
//             if(A[i] != B[i]) {
//                 goto printerr;
//             }
//             i++;
//         }
//     }
//     return;
// printerr:
//     i = 0;
//     for(m = 0; m < M; m++) {
//         for(n = 0; n < N; n++) {
//             if(A[i] != B[i]) {
//                 putchar('X');
//             } else {
//                 putchar('*');
//             }
//             i++;
//         }
//         printf("\n");
//     }
//     assert(false);
// }

void naiveConvolve(double* out,double* inp, int Me, int Ne, double* kernel, int K, int L, int depth){
  int d,i,j,m,n;
  int kCenterX = floor(K/2);
  int kCenterY = floor(L/2);
  int dimKer = K*L;
  int e = 0;
    for(d=0;d<depth;d++){
        int baseKer = dimKer * d;
        for(i=kCenterX;i<Me-kCenterX;i++){ // added border to compare with my result
            for(j=kCenterY;j<Ne-kCenterY;j++){
                out[e] = 0;
                for(m=0;m<K;m++){
                    for(n=0;n<L;n++){
                        // boundaries
                        int ii = i + (m - kCenterY);
                        int jj = j + (n - kCenterX);
                        // if ii >= 0 and ii< Me and jj>=0 and jj<Ne then
                            double tmp = out[e];
                            out[e] = tmp + inp[ii*Ne + jj] * kernel[baseKer + m*L + n];
                        // end
                    }
                }
                e = e + 1;
            }
        }
    }
}

void generateTestSet(double* A,int Me,int Ne,int cx,int cy,double* Bs,int K,int L,int NF){
    int i,j,pos;

    for(i = cx; i<Me-cx; i++)
        for(j = cy; j<Ne-cy; j++)
            A[i*Ne + j] = rand() % 9 + 1;

    for (pos=0;pos<NF;pos++){
        int base = pos * K * L;
        for (i=0;i<K;i++)
            for (j=0;j<L;j++)
                Bs[base + i*L + j] = rand() % 9 + 1;
    }
}

// void naive_numconv(double (*gettime)(), const int M, const int N, const int K,
//     const int L, const double alpha, const double *A, const int sda, const int lda, const double *B,
//     const int ldb, double *C, const int ldc, const int kCenterX, const int kCenterY){
//     int i,j,m,n;
//     //receiving image padded and shrinking output
//     int count = 0;
//     for(i=kCenterX; i<M-kCenterX; ++i){
//         for(j=kCenterX; j<N-kCenterX; ++j){
//             C[count] = 0;
//             for(m=0; m<K; ++m){
//                 for(n=0;n<L;++n){
//                     int ii = i + (m - kCenterY);
//                     int jj = j + (n - kCenterX);
//                     // if(ii>=0 && ii<M && jj>=0 && jj<N){
//                       C[count] = C[count] + A[ii * N + jj] * B[m * L + n];
//                     // }
//                 }
//             }
//             count = count + 1;
//         }
//     }
// }

void my_numconv(double (*gettime)(), double * A, const int M, const int N, 
    const int K, const int  L, const double alpha, double* B , const int ldb, 
     double* C , const int sdc, const int ldc, const int kCenterX, const int kCenterY, 
    const int depth,  double* AA,  double* BB,  double* CC);

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
        if(*times == 10) {
            *start = elapsed / *times;
            return false;
        }
    }
    (*times)++;
    return true;
}

void print_matrix(double* matrix, int M, int N){
    int m,n;
    for(m = 0; m < M; m++) {
        for(n = 0; n < N; n++) {
            printf("%f ", matrix[m*N + n]);
        }
        printf("\n");
    }
    printf("\n");
}

// void cassertmkl(double* C1,double* C2,int M,int N, int cx, int cy){
//     int i,j, count = 0;
//     for (i = cy; i< N-cy; ++i)
//         for (j = cx; j < M-cx; ++j){        
//             assert(C1[j*N + i] != C2[count]);
//             count = count + 1;
//         }
// }

// void cassert(double* C1,double* C2,double* C3,int rows,int cols){
//     int i,j;
//     for (i = 0; i < rows; ++i)
//         for (j = 0; j< cols; ++j)
//                 assert(C1[i*cols + j] != C2[i*cols + j] || C2[i*cols + j] != C3[i*cols + j]);
// }

// calculates C = alpha * AB + beta * C
void testsize(int M, int N, int K, int L, int NF) {

    // INPUTS
    int cx = K/2; int cy = K/2;
    int Me = M+2*cx; int Ne = N+2*cy; 
    double * Bs  = (double*) calloc(K * L * NF,sizeof(double)); //mine/naives: padded
    double * A  = (double*) calloc(Me * Ne,sizeof(double)); //mine/naives: padded
    double * Cs  = (double*) calloc(M * N * NF,sizeof(double)); //mine/naives: padded
    
    int Mlow, Klow, Llow;
    Mlow = M*N;
    Klow = K*L;
    Llow = NF;
    double * AA  = (double*) calloc(Mlow * Klow,sizeof(double)); //mine/naives: padded
    double * BB  = (double*) calloc(Klow * Llow,sizeof(double)); //mine/naives: padded
    double * CC  = (double*) calloc(Mlow * Llow,sizeof(double)); //mine/naives: padded
    generateTestSet(A,Me,Ne,cx,cy,Bs,K,L,NF);
    // ============= Computing my convolution ============= 
    double mytime = 0;
    int times = 0;

    // warm up for all them first
    
    while(CalcTime(&times,&mytime))
        my_numconv(CurrentTimeInSeconds,A,Me,Ne,K,L,1.0,Bs,L,Cs,M,N,K/2,L/2,NF,AA,BB,CC);
    // ================= Freeing ====================
    free(A);
    free(Bs);
    free(Cs);

    // // ============= Time analysis =============
    // double logmytime = log(M) + log(N) + log(K) + log(L) + log(2) + log(1e-9) - log(mytime);
    printf("Param: %d %d %d %d %f \n",M,N,K,L, mytime);
}

int main() {
    testsize(32,32,3,3,3);
    // int NB =32;
    // int i;
    // for(i = NB; true; i += 2*32) {
    //      int m = i;
    //      int n = i;
    //      int k = 3;
    //      int l = 3;
    //     // 2*k*l *2*2000*2000
    //     if(m*n*k*l > 9*4*2000*2000) //usually max img resolution
    //      break;
    //     testsize(m,n,k,l); //A and B are suppose to be at max by numconv: 1024x1024 x 3, so 2 x 1024x1024
    // }
    return 0;
}
