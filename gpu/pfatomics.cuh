#ifndef __GL_ATOMIC_
#define __GL_ATOMIC_

#include "pfcudefs.h"

__device__ __forceinline__ uint32 lanemask() {
    uint32 lanemask;
    asm ("mov.u32 %0, %%lanemask_lt;" : "=r"(lanemask));
    return (lanemask);
}

__device__ __forceinline__ uint32 atomicAggInc(uint32* counter) {
    uint32 mask = __activemask(), total = __popc(mask), prefix = __popc(mask & lanemask());
    int lowest_lane = __ffs(mask) - 1;
    uint32 warpRes = 0;
    if (prefix == 0) warpRes = atomicAdd(counter, total);
    warpRes = __shfl_sync(mask, warpRes, lowest_lane);
    uint32 thread_offset = prefix + warpRes;
    return thread_offset;
}

__device__ __forceinline__ int atomicAggInc(int* counter) {
    uint32 mask = __activemask(), total = __popc(mask), prefix = __popc(mask & lanemask());
    int lowest_lane = __ffs(mask) - 1;
    int warpRes = 0;
    if (prefix == 0) warpRes = atomicAdd(counter, total);
    warpRes = __shfl_sync(mask, warpRes, lowest_lane);
    int thread_offset = prefix + warpRes;
    return thread_offset;
}

#endif