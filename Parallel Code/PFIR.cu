#include <stdio.h>
#include <math.h>
#include <time.h>
#include "fdacoefs.h"

#define FILLEN 357
#define FS 20000
#define N  2000000
#define M  50
#define D  N/M
#define block_size 256

__global__ void PFIR_Kernel(float *in , float *coe ,float *out ){
	int tid;
	register int n;
	register int i,j;
	register float sum = 0;
	tid = blockIdx.x * blockDim.x + threadIdx.x;

	register int k,LIM;
	n = tid * 50;
	for (i = 0; i< 5; i+=1){
		j = 10 * i;
		LIM = FILLEN - j;
		for(k = 0; k < (LIM - 0)/M &&  n-k*M-(j+0) > 0; k++)
			sum += coe[k*M +(j+0)] *  in[n-k*M-(j+0)];		
		for(k = 0; k < (LIM - 1)/M &&  n-k*M-(j+1) > 0; k++)
			sum += coe[k*M +(j+1)] *  in[n-k*M-(j+1)];	
			
		for(k = 0; k < (LIM - 2)/M &&  n-k*M-(j+2) > 0; k++)
			sum += coe[k*M +(j+2)] *  in[n-k*M-(j+2)];	
		for(k = 0; k < (LIM - 3)/M &&  n-k*M-(j+3) > 0; k++)
			sum += coe[k*M +(j+3)] *  in[n-k*M-(j+3)];	

		for(k = 0; k < (LIM - 4)/M &&  n-k*M-(j+4) > 0; k++)
			sum += coe[k*M +(j+4)] *  in[n-k*M-(j+4)];	
				
				
		for(k = 0; k < (LIM - 5)/M &&  n-k*M-(j+5) > 0; k++)
			sum += coe[k*M +(j+5)] *  in[n-k*M-(j+5)];	

		for(k = 0; k < (LIM - 6)/M &&  n-k*M-(j+6) > 0; k++)
			sum += coe[k*M +(j+6)] *  in[n-k*M-(j+6)];	
		for(k = 0; k < (LIM - 7)/M &&  n-k*M-(j+7) > 0; k++)
			sum += coe[k*M +(j+7)] *  in[n-k*M-(j+7)];	

		for(k = 0; k < (LIM - 8)/M &&  n-k*M-(j+8) > 0; k++)
			sum += coe[k*M +(j+8)] *  in[n-k*M-(j+8)];	
		for(k = 0; k < (LIM - 9)/M &&  n-k*M-(j+9) > 0; k++)
			sum += coe[k*M +(j+9)] *  in[n-k*M-(j+9)];
		}		
	out[tid] += sum;
}
int main(){
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float time;
	register int i;
	int num_block = (D % block_size == 0) ? D/block_size : (D/block_size)+1;
	
	float *harmonics_H , *result_H , *coe_H;
	float *harmonics_D , *result_D , *coe_D;
	
	harmonics_H = (float *)malloc(sizeof(float) * N);	
	cudaMalloc((void **)&harmonics_D , sizeof(float) * N);
		
	result_H = (float *)calloc(D , sizeof(float));
	cudaMalloc((void **)&result_D   , sizeof(float) * D);
	
	coe_H  = (float *)malloc(sizeof(float)*FILLEN);
	cudaMalloc((void **)&coe_D , sizeof(float)* FILLEN);
	
	cudaMemset(result_D , 0 , sizeof(float)*D );
	
	for(i = 0; i < FILLEN; i++)
		coe_H[i] = COE[i];
		
	for(i = 0; i < N; i++)
		harmonics_H[i] = 2.75*sin(2*M_PI*i*50/FS + M_PI_4) + 3*sin(2*M_PI*i*400/FS + M_PI_2) + 1.25*sin(2*M_PI*i*1000/FS + 0);
	

	
	cudaEventRecord(start,0);
	
	cudaMemcpy(coe_D       , coe_H       ,FILLEN * sizeof(float),cudaMemcpyHostToDevice );
	cudaMemcpy(harmonics_D , harmonics_H ,sizeof(float) * N     , cudaMemcpyHostToDevice );
	PFIR_Kernel<<< num_block , block_size >>>(harmonics_D ,coe_D , result_D);
	cudaMemcpy(result_H    , result_D    ,sizeof(float)* D      , cudaMemcpyDeviceToHost );
		
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time,start,stop);


	//........................................................
	
	FILE *f = fopen("data1.dat","w");
	if(f == NULL)
		printf("ERROR <FILE> \n");
	else{
		for(i = 0;i < N; i++)
			fprintf(f, "%f \n",harmonics_H[i]);
	}
	fclose(f);
	
	f = fopen("data2.dat","w");
	if(f == NULL)
		printf("ERROR <FILE> \n");
	else{
		for(i = 0;i < D; i++)
			fprintf(f, "%f \n",result_H[i]);
	}
	fclose(f);

	printf("Elapsed time : %f ms \n" , time);
	
	cudaFree(harmonics_D);
	cudaFree(result_D);
	free(harmonics_H);
	free(result_H);

	return 0;
}
