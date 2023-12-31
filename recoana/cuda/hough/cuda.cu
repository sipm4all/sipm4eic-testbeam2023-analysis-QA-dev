// cuda.cu
#include <cuda_runtime.h>

/** 
    Hough space GPU model
    the basic computational unit is a 10x10x10 block in the (x,y,r) space
    the lattice is defined by the cell spacing (dx,dy,dr)

    the full Hough space can be mapped with Nx,Ny,Nr blocks

    for example, if we want to map a Hough space 
    x = [-64, 64]
    y = [-64, 64]
    r = [ 32, 96]
    with a cell spacing of (1,1,1) we need (Nx,Ny,Nr) = (16,16,8) blocks
    Nx = 128 / dx / 8 = 16
    Ny = 128 / dy / 8 = 16
    Nr =  64 / dr / 8 = 8

**/

float *gpu_x = nullptr;
float *gpu_y = nullptr;
float *gpu_h = nullptr;

float *gpu_xmap = nullptr;
float *gpu_ymap = nullptr;
float *gpu_rmap = nullptr;

const float x_min = -9.75;
const float x_stp = 0.5;
const float y_min = -9.75;
const float y_stp = 0.5;
const float r_min = 30.;
const float r_stp = 0.5;

__global__ void
hough_gpu_init(float *xmap, float *ymap, float *rmap, int Nx, int Ny, int Nr) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  int ix = threadIdx.x % 10;
  int iy = (threadIdx.x / 10) % 10;
  int ir = threadIdx.x / 100;
  
  int iX = blockIdx.x % Nx;
  int iY = (blockIdx.x / Nx) % Ny;
  int iR = blockIdx.x / (Nx * Ny);
  
  ix += iX * 10;
  iy += iY * 10;
  ir += iR * 10;
  
  float x = x_min + ix * x_stp;
  float y = y_min + iy * y_stp;
  float r = r_min + ir * r_stp;

  xmap[tid] = x;
  ymap[tid] = y;
  rmap[tid] = r;
  
}

__global__ void
hough_gpu_transform(float *xmap, float *ymap, float *rmap, float *x, float *y, float *h, int n) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  float cx = xmap[tid];
  float cy = ymap[tid];
  float cr = rmap[tid];

  h[tid] = 0.;
  for (int i = 0; i < n; ++i) {
    float dx = cx - x[i];
    float dy = cy - y[i];
    float dr = hypotf(dx, dy) - cr;
    float w = 0.11398351 * expf(-0.040816327 * dr * dr  );
    h[tid] += w;
  }
  
}

void
hough_init(float *cpu_xmap, float *cpu_ymap, float *cpu_rmap, int Nx, int Ny, int Nr)
{
  int Nh = 1000 * Nx * Ny * Nr;
  
  // alloc device memory
  cudaMalloc((void **)&gpu_x, 1024 * sizeof(float));
  cudaMalloc((void **)&gpu_y, 1024 * sizeof(float));
  cudaMalloc((void **)&gpu_h, Nh * sizeof(float));

  cudaMalloc((void **)&gpu_xmap, Nh * sizeof(float));
  cudaMalloc((void **)&gpu_ymap, Nh * sizeof(float));
  cudaMalloc((void **)&gpu_rmap, Nh * sizeof(float));

  // launch kernel
  dim3 block_size(1000, 1, 1);
  dim3 grid_size(Nx * Ny * Nr, 1, 1);
  hough_gpu_init<<<grid_size, block_size>>>(gpu_xmap, gpu_ymap, gpu_rmap, Nx, Ny, Nr);

  // copy data from device
  cudaMemcpy(cpu_xmap, gpu_xmap, Nh * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_ymap, gpu_ymap, Nh * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_rmap, gpu_rmap, Nh * sizeof(float), cudaMemcpyDeviceToHost);

}

void
hough_free()
{
  // free device memory
  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_h);
  
  cudaFree(gpu_xmap);
  cudaFree(gpu_ymap);
  cudaFree(gpu_rmap);
}

void
hough_transform(float *cpu_x, float *cpu_y, float *cpu_h, int cpu_n, int Nx, int Ny, int Nr)
{
  int Nh = 1000 * Nx * Ny * Nr;

  // copy data to device
  cudaMemcpy(gpu_x, cpu_x, cpu_n * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, cpu_y, cpu_n * sizeof(float), cudaMemcpyHostToDevice);
  
  // launch kernel
  dim3 block_size(1000, 1, 1);
  dim3 grid_size(Nx * Ny * Nr, 1, 1);
  hough_gpu_transform<<<grid_size, block_size>>>(gpu_xmap, gpu_ymap, gpu_rmap, gpu_x, gpu_y, gpu_h, cpu_n);

  // copy data from device
  cudaMemcpy(cpu_h, gpu_h, Nh * sizeof(float), cudaMemcpyDeviceToHost);
}

