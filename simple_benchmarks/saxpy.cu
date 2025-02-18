#include <cuda_runtime.h>
#include <iostream>

// Simple saxpy-like kernel with a single array for measuring data movement
__global__ void saxpy(int n, float a, float * x, float y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) x[i] = a*x[i] + y;
}

int main()
{
  // Set number of threads and threadblocks
  const size_t Nthreads = 32;
  const size_t Nblocks = 100000;

  // Kernel constants
  const float x = 1.0;
  const float a = 2.0;
  const float y = 0.432;

  // Array size
  const size_t N = Nthreads*Nblocks;

  // Host FP32 array
  float * h_x;
  cudaMallocHost(&h_x, N*sizeof(float));
  for (size_t i = 0; i < N ; i++)
  {
    h_x[i] = x;
  }

  // Device array
  float * d_x;
  cudaMalloc(&d_x, N*sizeof(float));

  // Copy data to device and launch kernels
  cudaMemcpy(d_x, h_x, N*sizeof(float), cudaMemcpyHostToDevice);

  saxpy<<<Nblocks,Nthreads>>>(N, a, d_x, y);

  // Copy back, compare, and tidy up
  cudaMemcpy(h_x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
  size_t nerrors = 0;
  for (size_t i = 0; i < N ; i++)
  {
    if ( abs(h_x[i]-(a*x+y))/(a*x+y) > 1.0e-7 ) nerrors++;
  }
  std::cout << "Number of errors: " << nerrors << "\n";

  cudaFree(h_x);
  cudaFree(d_x);

  return 0;
}
