// cuda.cu
#include <cuda_runtime.h>
#include <iostream>
#include <stdio.h>
#include "common.h"

static void HandleError( cudaError_t err, const char *file, int line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
    exit( EXIT_FAILURE );
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

data_t d_data;

__global__ void
hough_init_kernel(data_t data)
{
  int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= size) return;

  int ix = (tid % data.bins.x);
  int iy = (tid / data.bins.x) % data.bins.y;
  int ir = (tid / (data.bins.x * data.bins.y)) % data.bins.r;
  int it = (tid / (data.bins.x * data.bins.y * data.bins.r)) % data.bins.t;

  data.map.x[tid] = data.min.x + (0.5 + ix) * (data.max.x - data.min.x) / data.bins.x;
  data.map.y[tid] = data.min.y + (0.5 + iy) * (data.max.y - data.min.y) / data.bins.y;
  data.map.r[tid] = data.min.r + (0.5 + ir) * (data.max.r - data.min.r) / data.bins.r;
  data.map.t[tid] = data.min.t + (0.5 + it) * (data.max.t - data.min.t) / data.bins.t;
  
}

__global__ void
hough_transform_kernel(data_t data)
{
  int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= size) return;

  float cx = data.map.x[tid];
  float cy = data.map.y[tid];
  float cr = data.map.r[tid];
  float ct = data.map.t[tid];

  float ts1 = 1. / (data.sigma.t * sqrtf(2. * M_PI));
  float ts2 = -0.5 / (data.sigma.t * data.sigma.t);
  
  data.hough.h[tid] = 0.;
  for (int i = 0; i < data.points.n; ++i) {
    float dx = data.points.x[i] - cx;
    float dy = data.points.y[i] - cy;
    float dt = data.points.t[i] - ct;
    float dr = hypotf(dx, dy) - cr;
    float w = 0.11398351 * expf(-0.040816327 * dr * dr  ); // sigma = 3.5 
    float wt = ts1 * expf(ts2 * dt * dt  );
    data.hough.h[tid] += (w * wt);
  }
  
}

__global__ void
find_max_kernel(data_t data)
{
  __shared__ float shm[256];
  __shared__ int shmi[256];
  
  int size = data.bins.x * data.bins.y * data.bins.r * data.bins.t;
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  if (gid >= size) return;

  int tid = threadIdx.x;
  int bid = blockIdx.x;

  shm[tid] = data.hough.h[gid];
  shmi[tid] = gid;
  __syncthreads();

  for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
    if (tid < stride) {
      if (shm[tid + stride] > shm[tid]) {
	shm[tid] = shm[tid + stride];
	shmi[tid] = shmi[tid + stride];
      }
    }
    __syncthreads();
  }

  if (tid == 0) {
    data.hough.rh[bid] = shm[0];
    data.hough.rhi[bid] = shmi[0];
  }

}

void
hough_init(data_t h_data)
{
  int size = h_data.bins.x * h_data.bins.y * h_data.bins.r * h_data.bins.t;
  int grid_size = 1 + (size - 1) / 256;
  
  /** alloc device memory for data map **/
  HANDLE_ERROR( cudaMalloc((void **)&d_data.map.x, size * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.map.y, size * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.map.r, size * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.map.t, size * sizeof(float)) );
  
  /** launch init kernel to populate data map **/
  d_data.min = h_data.min;
  d_data.max = h_data.max;
  d_data.bins = h_data.bins;
  hough_init_kernel<<<grid_size, 256>>>(d_data);

  /** copy data map from device **/
  HANDLE_ERROR( cudaMemcpy(h_data.map.x, d_data.map.x, size * sizeof(float), cudaMemcpyDeviceToHost) );
  HANDLE_ERROR( cudaMemcpy(h_data.map.y, d_data.map.y, size * sizeof(float), cudaMemcpyDeviceToHost) );
  HANDLE_ERROR( cudaMemcpy(h_data.map.r, d_data.map.r, size * sizeof(float), cudaMemcpyDeviceToHost) );
  HANDLE_ERROR( cudaMemcpy(h_data.map.t, d_data.map.t, size * sizeof(float), cudaMemcpyDeviceToHost) );
  
  /** alloc device memory for data points **/
  HANDLE_ERROR( cudaMalloc((void **)&d_data.points.x, 1024 * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.points.y, 1024 * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.points.t, 1024 * sizeof(float)) );
  
  HANDLE_ERROR( cudaMalloc((void **)&d_data.hough.h, size * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.hough.rh, grid_size * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_data.hough.rhi, grid_size * sizeof(int)) );

}

void
hough_free()
{

  /** free device memory for data map **/
  cudaFree(d_data.map.x);
  cudaFree(d_data.map.y);
  cudaFree(d_data.map.r);
  cudaFree(d_data.map.t);

  /** free device memory for data points **/
  cudaFree(d_data.points.x);
  cudaFree(d_data.points.y);
  cudaFree(d_data.points.t);
  
  /** free device memory for data hough **/
  cudaFree(d_data.hough.h);
  cudaFree(d_data.hough.rh);
  cudaFree(d_data.hough.rhi);
  
}

void
hough_transform(data_t h_data)
{
  int size = h_data.bins.x * h_data.bins.y * h_data.bins.r * h_data.bins.t;
  int grid_size = 1 + (size - 1) / 256;

  /** copy data points to device **/
  HANDLE_ERROR( cudaMemcpy(d_data.points.x, h_data.points.x, h_data.points.n * sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_data.points.y, h_data.points.y, h_data.points.n * sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_data.points.t, h_data.points.t, h_data.points.n * sizeof(float), cudaMemcpyHostToDevice) );
  
  // launch kernel
  d_data.points.n = h_data.points.n;
  d_data.sigma.t = h_data.sigma.t;
  hough_transform_kernel<<<grid_size, 256>>>(d_data);
  find_max_kernel<<<grid_size, 256>>>(d_data);
  
  // copy data from device
  HANDLE_ERROR( cudaMemcpy(h_data.hough.rh, d_data.hough.rh, grid_size * sizeof(float), cudaMemcpyDeviceToHost) );
  HANDLE_ERROR( cudaMemcpy(h_data.hough.rhi, d_data.hough.rhi, grid_size * sizeof(int), cudaMemcpyDeviceToHost) );
}

