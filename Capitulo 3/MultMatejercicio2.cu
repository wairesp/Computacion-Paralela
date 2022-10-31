%%cu
#include <cuda.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
using namespace std;

#define DIM 8
#define BlockSize 4

__global__ void multi(int *B, int *C, int *A, int width)
{
    int cvalue = 0;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;


    if (row > width || col > width) return;

    for (int e = 0; e < width; ++e){
        cvalue += B[row * width + e] * C[e * width + col];
    }
    A[row * width + col] = cvalue;

}

void matrixmulti(int B[][DIM], int C[][DIM], int A[][DIM]){
    int *dev_a, *dev_b, *dev_c;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //allocate memory on global memory of gpu
    cudaMalloc((void**)&dev_a, ((DIM)*(DIM))*sizeof(int));
    cudaMalloc((void**)&dev_b, ((DIM)*(DIM))*sizeof(int));
    cudaMalloc((void**)&dev_c, ((DIM)*(DIM))*sizeof(int));

    //Copy array B and C on device allocated memory
    cudaMemcpy(dev_a, B, ((DIM * DIM)) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, C, ((DIM * DIM)) * sizeof(int), cudaMemcpyHostToDevice);

    //two dimension threads
    dim3 dimBlock(BlockSize, BlockSize);
    dim3 dimGrid((DIM + dimBlock.x - 1) / dimBlock.x, (DIM + dimBlock.y - 1) / dimBlock.y);

    //call the kernel function multi
    cudaEventRecord(start);
    multi << < dimGrid, dimBlock >> >(dev_a, dev_b, dev_c, DIM);
    cudaEventRecord(stop);

    //retrieve array A from device memory
    cudaMemcpy(A, dev_c, ((DIM * DIM)) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaEventSynchronize(stop);

    /*for (int i = 0; i < DIM; i++){
        for (int j = 0; j < DIM; j++){
            printf("A(%d,%d) = %d \n", i, j, A[i][j]);
        }
    }*/

    //free the memory
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
}


int main(){

	srand(time(0));
	auto A = new int[DIM][DIM];
	auto B = new int[DIM][DIM];
	auto C = new int[DIM][DIM];

	//populate the arrays A and B
	for (int i = 0; i<DIM; i++){
		for (int j = 0; j < DIM; j++){
			A[i][j] = rand() % 10;
			B[i][j] = rand() % 10;
		}
	}

  printf("VALORES DE MATRIZ B \n");
	for (int y = 0; y < DIM; y++)
	{
		for (int x = 0; x < DIM; x++)
		{
			printf("%d", B[y][x]);
      printf(" ");
		}
		printf(" \n");
	}
  
  printf("UCSP\nVALORES DE MATRIZ C \n");
	for (int y = 0; y < DIM; y++)
	{
		for (int x = 0; x < DIM; x++)
		{
			printf("%d", C[y][x]);
      printf(" ");
		}
		printf(" \n");
	}

	
  matrixmulti(B,C,A);

  printf("VALORES DE MATRIZ A \n");
	for (int y = 0; y < DIM; y++)
	{
		for (int x = 0; x < DIM; x++)
		{
			printf("%d", A[y][x]);
      printf(" ");
		}
		printf(" \n");
	}

	//delete arrays
	delete[]A;
	delete[]B;
	delete[]C;
}