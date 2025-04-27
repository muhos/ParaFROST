/***********************************************************************[count.cu]
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

#include "solver.hpp"
#include "count.cuh"
#include "shared.cuh"
#include "reduce.cuh"


namespace ParaFROST {

	uint32 hostCBlocks[MAXREDUCEBLOCKS];
	uint32 hostLBlocks[MAXREDUCEBLOCKS];

	__device__
		uint32 gcounter;
	__device__
		uint32 devCBlocks[MAXREDUCEBLOCKS];
	__device__
		uint32 devLBlocks[MAXREDUCEBLOCKS];

	__global__ void reset_counter() { gcounter = 0; }

	__global__ void print_counter() { printf("c gcounter = %d\n", gcounter); }

	__global__ void check_counter(const uint32 checksum) { assert(checksum == gcounter); }

	__global__
	void cnt_cls(const CNF* __restrict__ cnf)
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

	__global__
	void cnt_lits(const CNF* __restrict__ cnf)
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

	__global__
	void cnt_cls_lits(const CNF* __restrict__ cnf)
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

	void parcountCls(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(blockSize, sizeof(uint32));
		cnt_cls << <nBlocks, blockSize, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		inf.n_cls_after = seqreduceBlocks(hostCBlocks, nBlocks);
	}

	void parcountLits(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(blockSize, sizeof(uint32));
		cnt_lits << <nBlocks, blockSize, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		inf.n_lits_after = seqreduceBlocks(hostLBlocks, nBlocks);
	}

	void parcountAll(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(blockSize, sizeof(uint32) * 2);
		cnt_cls_lits << <nBlocks, blockSize, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		seqreduceBlocks(hostCBlocks, hostLBlocks, nBlocks);
	}

	void Solver::countCls(const bool& host)
	{
		if (host) {
			assert(!hcnf->empty());
			inf.n_cls_after = 0;
			for (uint32 i = 0; i < hcnf->size(); i++) {
				SCLAUSE& c = hcnf->clause(i);
				if (c.original() || c.learnt())
					inf.n_cls_after++;
			}
		}
		else parcountCls(cnf);
	}

	void Solver::countLits(const bool& host)
	{
		if (host) {
			assert(!hcnf->empty());
			inf.n_lits_after = 0;
			for (uint32 i = 0; i < hcnf->size(); i++) {
				SCLAUSE& c = hcnf->clause(i);
				if (c.original() || c.learnt())
					inf.n_lits_after += c.size();
			}
		}
		else parcountLits(cnf);
	}

	void Solver::countAll(const bool& host)
	{
		if (host) {
			assert(!hcnf->empty());
			inf.n_cls_after = 0, inf.n_lits_after = 0;
			for (uint32 i = 0; i < hcnf->size(); i++) {
				SCLAUSE& c = hcnf->clause(i);
				if (c.original() || c.learnt())
					inf.n_cls_after++, inf.n_lits_after += c.size();
			}
		}
		else parcountAll(cnf);
	}

	void Solver::logReductions()
	{
		int64 varsRemoved = int64(inf.n_del_vars_after) + nForced;
		int64 clsRemoved = int64(inf.nClauses) - inf.n_cls_after;
		int64 litsRemoved = int64(inf.nLiterals) - inf.n_lits_after;
		const char* header = "  %s%-10s  %-10s %-10s %-10s%s";
		LOG1(header, CREPORT, " ", "Variables", "Clauses", "Literals", CNORMAL);
		const char* rem = "  %s%-10s: %s%-9lld  %c%-8lld  %c%-8lld%s";
		const char* sur = "  %s%-10s: %s%-9d  %-9d  %-9d%s";
		LOG1(rem, CREPORT, "Removed", CREPORTVAL,
			-varsRemoved,
			clsRemoved < 0 ? '+' : '-', abs(clsRemoved),
			litsRemoved < 0 ? '+' : '-', abs(litsRemoved), CNORMAL);
		LOG1(sur, CREPORT, "Survived", CREPORTVAL,
			maxActive(),
			inf.n_cls_after,
			inf.n_lits_after, CNORMAL);
	}

}