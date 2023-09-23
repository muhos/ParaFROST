/***********************************************************************[count.cuh]
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

#ifndef __GPU_COUNT_
#define __GPU_COUNT_

#include "proofutils.cuh"

namespace ParaFROST {

	#define MAXREDUCEBLOCKS 256
	uint32 hostCBlocks[MAXREDUCEBLOCKS];
	uint32 hostLBlocks[MAXREDUCEBLOCKS];
	__device__ uint32 devCBlocks[MAXREDUCEBLOCKS];
	__device__ uint32 devLBlocks[MAXREDUCEBLOCKS];
	__device__ uint32 gcounter;

	// kernels
	__global__ void reset_counter() { gcounter = 0; }

	__global__ void print_counter() { printf("c gcounter = %d\n", gcounter); }

	__global__ void check_counter(const uint32 checksum) { assert(checksum == gcounter); }

	__global__ void cnt_proof_verify(const uint32* __restrict__ literals, const uint32 numLits)
	{
		grid_t tid = 0;
		while (tid < numLits) {
			addr_t lbyte = DC_PTRS->d_lbyte;
			const uint32 lit = literals[tid];
			if (lit & 0xF0000000) lbyte[lit] = 5;
			else if (lit & 0x0FE00000) lbyte[lit] = 4;
			else if (lit & 0x001FC000) lbyte[lit] = 3;
			else if (lit & 0x00003F80) lbyte[lit] = 2;
			else lbyte[lit] = 1;
			printf(" literal(%d) has %d bytes of its original\n", SIGN(lit) ? -int(ABS(lit)) : ABS(lit), lbyte[lit]);
			gcounter += lbyte[lit];
			tid++;
		}
		printf(" total = %d\n", gcounter);
	}

	__global__ void cnt_proof(const uint32* __restrict__ literals, const uint32 numLits)
	{
		uint32* sh_bytes = SharedMemory<uint32>();
		grid_t tid = global_tx_off;
		uint32 nbytes = 0;
		while (tid < numLits) {
			addr_t lbyte = DC_PTRS->d_lbyte;
			uint32 lit = literals[tid];
			countBytes(lit, lbyte[lit], nbytes);
			grid_t off = tid + blockDim.x;
			if (off < numLits) {
				lit = literals[off];
				countBytes(lit, lbyte[lit], nbytes);
			}
			tid += stride_x_off;
		}
		loadShared(sh_bytes, nbytes, numLits);
		sharedReduce(sh_bytes, nbytes);
		warpReduce(sh_bytes, nbytes);
		if (!threadIdx.x) devLBlocks[blockIdx.x] = nbytes;
	}

	__global__ void cnt_cls(const CNF* __restrict__ cnf)
	{
		uint32* sh_rCls = SharedMemory<uint32>();
		grid_t tid = global_tx_off;
		uint32 nCls = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nCls++;
			uint32 off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nCls++;
			}
			tid += stride_x_off;
		}
		loadShared(sh_rCls, nCls, cnf->size());
		sharedReduce(sh_rCls, nCls);
		warpReduce(sh_rCls, nCls);
		if (!threadIdx.x) devCBlocks[blockIdx.x] = nCls;
	}

	__global__ void cnt_lits(const CNF* __restrict__ cnf)
	{
		uint32* sh_rLits = SharedMemory<uint32>();
		grid_t tid = global_tx_off;
		uint32 nLits = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nLits += c1.size();
			uint32 off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nLits += c2.size();
			}
			tid += stride_x_off;
		}
		loadShared(sh_rLits, nLits, cnf->size());
		sharedReduce(sh_rLits, nLits);
		warpReduce(sh_rLits, nLits);
		if (!threadIdx.x) devLBlocks[blockIdx.x] = nLits;
	}

	__global__ void cnt_cls_lits(const CNF* __restrict__ cnf)
	{
		uint32* sh_rCls = SharedMemory<uint32>();
		uint32* sh_rLits = sh_rCls + blockDim.x;
		grid_t tid = global_tx_off;
		uint32 nCls = 0;
		uint32 nLits = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nCls++, nLits += c1.size();
			grid_t off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nCls++, nLits += c2.size();
			}
			tid += stride_x_off;
		}
		loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
		sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
		warpReduce(sh_rCls, nCls, sh_rLits, nLits);
		if (!threadIdx.x) {
			grid_t bx = blockIdx.x;
			devCBlocks[bx] = nCls;
			devLBlocks[bx] = nLits;
		}
	}

	__global__ void copy_if_k(uint32* __restrict__ dest, CNF* __restrict__ src)
	{
		grid_t tid = global_tx;
		while (tid < src->size()) {
			SCLAUSE& c = src->clause(tid);
			if (c.original() || c.learnt()) {
				uint32* d = dest + atomicAdd(&gcounter, c.size());
				forall_clause(c, s) { *d++ = *s; }
			}
			tid += stride_x;
		}
	}

}

#endif