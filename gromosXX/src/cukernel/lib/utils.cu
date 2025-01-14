/**
 * @file utils.cu
 * implementation of utils
 */

#include <cuda.h>
#include <iostream>
#include "utils.h"

int cudakernel::checkError(const char * err_msg) {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess)
    std::cout << "CUDA-ERROR " << err_msg << ": " << cudaGetErrorString(error) << std::endl;

  return (int) error;
}

