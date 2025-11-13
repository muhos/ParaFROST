/***********************************************************************[definitions.cuh]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

#ifndef __CU_DEFS_
#define __CU_DEFS_  

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "logging.hpp"

namespace ParaFROST {

	extern cudaDeviceProp devProp;
	extern uint32 maxGPUThreads;
	extern size_t maxGPUSharedMem;

	#if !defined(_PFROST_H_)
	#define _PFROST_H_ inline __host__
	#endif

	#if !defined(_PFROST_D_)
	#define _PFROST_D_ __device__
	#endif

	#if !defined(_PFROST_IN_D_)
	#define _PFROST_IN_D_ inline __device__
	#endif

	#if !defined(_PFROST_H_D_)
	#define _PFROST_H_D_ inline __host__ __device__
	#endif

	#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
	#define LASTERR(MESSAGE)	\
		do { \
				cudaError_t ERR = cudaGetLastError(); \
				if (cudaSuccess != ERR) { \
					LOGERRORN("%s(%i): %s due to (%d) %s", __FILE__, __LINE__, MESSAGE, static_cast<int>(ERR), cudaGetErrorString(ERR)); \
					cudaDeviceReset(); \
					exit(1); \
				} \
		} while(0)
	#else
	#define LASTERR(MESSAGE)	do { } while(0)
	#endif

	#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
		#define CHECK(FUNCCALL) \
			do { \
				const cudaError_t returncode = FUNCCALL; \
				if (returncode != cudaSuccess) { \
					LOGERRORN("CUDA runtime failure due to %s", cudaGetErrorString(returncode)); \
					cudaDeviceReset(); \
					exit(1); \
				} \
			} while(0)
	#else
		#define CHECK(FUNCCALL) FUNCCALL
	#endif

	#define SYNC(STREAM) CHECK(cudaStreamSynchronize(STREAM));

	#define SYNCALL CHECK(cudaDeviceSynchronize());

}

#endif
