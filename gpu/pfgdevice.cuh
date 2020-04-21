
#ifndef __GL_DEVICE_
#define __GL_DEVICE_

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

__device__ inline int atomicAggInc(int* counter)
{
	auto g = cg::coalesced_threads();
	int res;
	if (g.thread_rank() == 0)
		res = atomicAdd(counter, g.size());
	return g.shfl(res, 0) + g.thread_rank();
}

__device__ inline unsigned int atomicAggInc(unsigned int* counter)
{
	auto g = cg::coalesced_threads();
	unsigned int res;
	if (g.thread_rank() == 0)
		res = atomicAdd(counter, g.size());
	return g.shfl(res, 0) + g.thread_rank();
}

__device__ inline unsigned long long int atomicAggInc(unsigned long long int* counter)
{
	auto g = cg::coalesced_threads();
	unsigned long long int res;
	if (g.thread_rank() == 0)
		res = atomicAdd(counter, g.size());
	return g.shfl(res, 0) + g.thread_rank();
}

#endif