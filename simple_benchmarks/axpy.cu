#include <cuda_runtime.h>
#include <iostream>

// Set floating point type
#define FTYPE float

// Simple axpy-like kernel with a single array for measuring data movement
__global__ void axpy(int n, FTYPE a, FTYPE * x, FTYPE y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) x[i] = a*x[i] + y;
}

int main()
{
  // Set number of threads and threadblocks
  const size_t Nthreads = 32;
  const size_t Nblocks = 100000;

  std::cout << "Running with array size " << Nthreads << "x" << Nblocks << "\n";
  std::cout << "1 array element = " << sizeof(FTYPE) << " Bytes\n";

  // Kernel constants
  const FTYPE x = 1.0;
  const FTYPE a = 2.0;
  const FTYPE y = 3.0;

  // Array size
  const size_t N = Nthreads*Nblocks;

  // Host array
  FTYPE * h_x;
  cudaMallocHost(&h_x, N*sizeof(FTYPE));
  for (size_t i = 0; i < N ; i++)
  {
    h_x[i] = x;
  }

  // Device array
  FTYPE * d_x;
  cudaMalloc(&d_x, N*sizeof(FTYPE));

  // Copy data to device and launch kernels
  cudaMemcpy(d_x, h_x, N*sizeof(FTYPE), cudaMemcpyHostToDevice);

  axpy<<<Nblocks,Nthreads>>>(N, a, d_x, y);

  // Copy back, compare, and tidy up
  cudaMemcpy(h_x, d_x, N*sizeof(FTYPE), cudaMemcpyDeviceToHost);
  size_t nerrors = 0;
  for (size_t i = 0; i < N ; i++)
  {
    if ( abs((h_x[i]-(a*x+y))/(a*x+y)) > 1.0e-7 ) nerrors++;
  }
  std::cout << "Number of errors: " << nerrors << "\n";

  cudaFree(h_x);
  cudaFree(d_x);

  return 0;
}
