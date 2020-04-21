#ifndef __CU_DEFS_
#define __CU_DEFS_  

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cuda.h>
#include <device_launch_parameters.h>
#include "pfdtypes.h"

#define SHAREDMEM_ALLOWED 0x3000
#define MASTER_GPU 0
#define BLOCK1D 256

// macros for CUDA errors check 
inline cudaError_t CHECK(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
#endif
	return result;
}
inline void _getLstErr(const char* errorMessage, const char* file, const int line) {
#if defined(DEBUG) || defined(_DEBUG)
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) {
		fprintf(stderr, "%s(%i) : Last CUDA error : %s : (%d) %s.\n", file, line, errorMessage, static_cast<int>(err), cudaGetErrorString(err));
		cudaDeviceReset();
		exit(1);
	}
#endif
}
#define LOGERR(msg) _getLstErr(msg, __FILE__, __LINE__)
// Compute #cores per SM (only for GPU info printing) 
inline int SM2Cores(int major, int minor) {
	// Defines for GPU Architecture types (using the SM version to determine # of cores per SM)
	typedef struct {
		int SM;  // arch defined in hex
		int Cores;
	} SM;
	SM nCores[] = {
		{0x30, 192},
		{0x32, 192},
		{0x35, 192},
		{0x37, 192},
		{0x50, 128},
		{0x52, 128},
		{0x53, 128},
		{0x60,  64},
		{0x61, 128},
		{0x62, 128},
		{0x70,  64},
		{0x72,  64},
		{0x75,  64},
		{-1, -1} };

	int index = 0;

	while (nCores[index].SM != -1) {
		if (nCores[index].SM == ((major << 4) + minor)) {
			return nCores[index].Cores;
		}
		index++;
	}
	printf(
		"MapSMtoCores for SM %d.%d is undefined."
		"  Default to use %d Cores/SM\n",
		major, minor, nCores[index - 1].Cores);
	return nCores[index - 1].Cores;
}

#endif
