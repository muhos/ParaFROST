#ifndef __GRID_INFO_
#define __GRID_INFO_

#include "pfcudefs.h"

// x
_PFROST_D_ uint32 global_bx(void) { return blockDim.x * blockIdx.x; }
_PFROST_D_ uint32 global_bx_off(void) { return (blockDim.x << 1) * blockIdx.x; }
_PFROST_D_ uint32 global_tx(void) { return global_bx() + threadIdx.x; }
_PFROST_D_ uint32 global_tx_off(void) { return global_bx_off() + threadIdx.x; }
_PFROST_D_ uint32 stride_x(void) { return blockDim.x * gridDim.x; }
_PFROST_D_ uint32 stride_x_off(void) { return (blockDim.x << 1) * gridDim.x; }
// y
_PFROST_D_ uint32 global_by(void) { return blockDim.y * blockIdx.y; }
_PFROST_D_ uint32 global_ty(void) { return global_by() + threadIdx.y; }
_PFROST_D_ uint32 stride_y(void) { return blockDim.y * gridDim.y; }

#endif