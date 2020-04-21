/***********************************************************************[pfdvec.h]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Technische Universiteit Eindhoven (TU/e).

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

#ifndef __DVEC_
#define __DVEC_

#include "pfcudefs.h"

#define MAX_CAP 0xffffffff
#define MAX_BLOCK 32
#define MAX_GRID 64

template<class T>
class cuVec {
	T* _mem;
	uint32 sz, cap;
public:
	cuVec() { _mem = NULL, sz = 0, cap = 0; }
	~cuVec() { _mem = NULL, sz = 0, cap = 0; }
	void alloc(T* head, const uint32& cap) { _mem = head, this->cap = cap; }
	void copyFrom(cuVec<T>* copy);
	void copyFromIf(cuVec<T>* copy, const T& val = 0);
	void _push(const T& val) { _mem[sz++] = val; }
	void _pop(void) { sz--; }
	void _shrink(const uint32& n) { sz -= n; }
	__device__ void push(const T& val) {
		uint32 idx = atomicAdd(&sz, 1);
		assert(idx < cap);
		_mem[idx] = val;
	}
	__device__ void pop(void) { uint32 idx = atomicSub(&sz, 1); assert(idx < MAX_CAP); }
	__device__ void shrink(const uint32& n) { uint32 idx = atomicSub(&sz, n); assert(idx < MAX_CAP); }
	__host__ __device__ cuVec<T>& operator=(cuVec<T>& rhs) { assert(0); }
	__host__ __device__ const T& operator [] (const uint32& idx) const { assert(idx < sz); return _mem[idx]; }
	__host__ __device__ T&       operator [] (const uint32& idx) { assert(idx < sz); return _mem[idx]; }
	__host__ __device__ operator T* (void) { return _mem; }
	__host__ __device__ inline T* data(void) { return _mem; }
	__host__ __device__ inline T& at(const uint32& idx) { assert(idx < sz); return _mem[idx]; }
	__host__ __device__ inline bool empty(void) const { return sz == 0; }
	__host__ __device__ inline uint32 size(void) const { return sz; }
	__host__ __device__ inline uint32 capacity(void) const { return cap; }
	__host__ __device__ void resize(const uint32& n) { assert(n <= cap); sz = n; }
	__host__ __device__ void clear(void) { sz = 0; }
	__host__ __device__ void print(const bool& litType = false) { 
		printf("c | GPU Vector (size = %d)->[", sz);
		if (litType) {
			for (uint32 i = 0; i < sz; i++) {
				assert(_mem[i]);
				printf("%2d ", (_mem[i] & 1) ? -int(ABS(_mem[i])) : int(ABS(_mem[i])));
			}
		}
		else {
			for (uint32 i = 0; i < sz; i++)
				printf("%2d ", _mem[i]);
		}
		printf("]\n");
	}
};


template<class T>
__global__ void vecCopy(cuVec<T>* dest, cuVec<T>* src)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	uint32 stride = blockDim.x * gridDim.x;
	while (i < src->size()) {
		if (i == src->size() - 1) dest->resize(src->size());
		*(*dest + i) = *(*src + i);
		i += stride;
	}
}

template<class T>
__global__ void vecCopyIf(cuVec<T>* dest, cuVec<T>* src, T val)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	uint32 stride = blockDim.x * gridDim.x;
	while (i < src->size()) {
		if (src->at(i) != val) dest->push(src->at(i));
		i += stride;
	}
}


template<class T>
void cuVec<T>::copyFrom(cuVec<T>* copy) {
	assert(copy->size() <= cap);
	vecCopy<T> << <MAX_BLOCK, MAX_GRID >> > (this, copy);
	LOGERR("Device vector copy failed");
	CHECK(cudaDeviceSynchronize());
}

template <class T>
void cuVec<T>::copyFromIf(cuVec<T>* copy, const T& val) {
	assert(copy->size() <= cap);
	vecCopyIf<T> << <MAX_BLOCK, MAX_GRID >> > (this, copy, val);
	LOGERR("Device vector copy failed");
	CHECK(cudaDeviceSynchronize());
}

#endif