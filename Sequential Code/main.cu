#include <stdio.h>
#include <math.h>
#include <time.h>
#include "fdacoefs.h"

#define FILLEN 357
#define  FS 20000
#define  N  2000000
#define  M  50
#define  D  N/M

void PolyFIR(float *Out, float *harm, float *In );
void FIR(float *out, float *h, int a);
int write(float *a , float *b);

int main(){

    register int i;
    double time_spend;
    float *harmonics = (float *) malloc(sizeof(float) * N);
    float *out = (float *) calloc(D , sizeof(float));
    float *RES = (float *) calloc(D , sizeof(float));

    //Generating Harmonics
    for(i = 0; i < N; i++)
        harmonics[i] = 2.75*sin(2*M_PI*i*50/FS + M_PI_4) + 3*sin(2*M_PI*i*400/FS + M_PI_2) + 1.25*sin(2*M_PI*i*1000/FS + 0) ;

    clock_t start_time = clock();
    PolyFIR(RES,harmonics,out);
    clock_t stop_time  = clock();

    time_spend = (double)(stop_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed Time (ms)\t %f \n" , time_spend * 1000);

    //For checking Result in Matlab
    write(harmonics , RES);
    return 0 ;
}

void PolyFIR(float *Out, float *harm, float *In ){
    register int i,j;
    for(i = 0; i < M; i++){
       FIR(In,harm,i);
       for(j = 0; j < D; j++)
        Out[j] = Out[j] + In[j];
    }

}

void FIR(float *out, float *harmonics, int start){
    register int i,j,k;
    j = 0;
    register float res = 0.0;
    register int Mul;
    register int LIM = (FILLEN - start)/M;

    for(i = 0; i < N; i += M){
        for(k = 0; k < LIM; k++){
            Mul = M*k;
            if ( (i-Mul-start)>0 )
                res = res + COE[Mul + (start-0)] * harmonics[i-Mul-start];
        }
        out[j++] = res;
        res = 0;

    }

}

int write(float *a , float *b){
    int i;
    FILE *f;
    f = fopen("data1.dat","w");
    if (f == NULL)
        return -1;
    for(i=0;i<N;i++)
        fprintf(f, "%f \n", a[i]);
    fclose(f);

    f = fopen("data2.dat","w");
    if (f == NULL)
        return -1;
    for(i=0;i<D;i++)
        fprintf(f, "%f \n", b[i]);
    fclose(f);
    return 0;
}
