/***********************************************************************[elimination.cu]
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
#include "grid.cuh"
#include "timer.cuh"
#include "options.cuh"
#include "bounded.cuh"
#include "subsume.cuh"
#include "blocked.cuh"
#include "redundancy.cuh"
#include "definitions.cuh"
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>
#include <cub/thread/thread_operators.cuh>

namespace ParaFROST {

	__global__ 
	void reset_id() { lastEliminatedID = -1; }

	__global__ 
	void print_id() { printf("c lastEliminatedID = %d\n", lastEliminatedID); }

	inline void vePhase1(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cuVecB* proof,
		uint32* ucnt,
		const uint32& max_cnts,
		uint32* type,
		const uint32& max_types,
		uint32* rpos,
		const uint32& max_poss,
		S_REF* rref,
		const uint32& max_refs,
		const bool& in)
	{
		LOGN2(2, "  resetting last eliminated ID to -1.. ");
		reset_id << <1, 1 >> > ();
		LOGDONE(2, 5);
		LOGN2(2, "  configuring phase 1 with ");
		grid_t nThreads = BLVE1;
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT1 * sizeof(uint32));
		OPTIMIZEBLOCKSELIM(vars->numElected, nThreads, smemSize, gopts.ve);
		LOGENDING(2, 5, "(%d/%d ths, %d bls) and %zd KB shared memory", nThreads, BLVE1, nBlocks, smemSize / KBYTE);
#if VE_DBG
		if (in) in_ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, max_cnts, type, max_types, rpos, max_poss, rref, max_refs);
		else	   ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, max_cnts, type, max_types, rpos, max_poss, rref, max_refs);
#else
		if (in) in_ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, max_cnts, type, max_types, rpos, max_poss, rref, max_refs);
		else	   ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, max_cnts, type, max_types, rpos, max_poss, rref, max_refs);
#endif
		LASTERR("BVE Phase-1 failed");
		SYNC(0);
	}

	inline void vePhase2(
		VARS* vars,
		uint32* rpos,
		S_REF* rref,
		cudaStream_t* streams,
		cuMM& cumm,
		CACHER& cacher)
	{
		const uint32 cs_offset = cumm.pinnedCNF()->size();
		const S_REF data_offset = cumm.pinnedCNF()->data().size;
		size_t tb1 = 0, tb2 = 0;
		cub::DeviceScan::ExclusiveScan(NULL, tb1, rpos, rpos, cub::Sum(), cs_offset, vars->numElected);
		cub::DeviceScan::ExclusiveScan(NULL, tb2, rref, rref, cub::Sum(), data_offset, vars->numElected);
		size_t tmpcap = tb1 + tb2;
		addr_t ts1 = NULL, ts2 = NULL;
		addr_t tmpmem = (addr_t)((tmpcap > cumm.scatterCap()) ? cacher.allocate(tmpcap) : cumm.scatter());
		ts1 = tmpmem, ts2 = ts1 + tb1;
		cub::DeviceScan::ExclusiveScan(ts1, tb1, rpos, rpos, cub::Sum(), cs_offset, vars->numElected, streams[0]);
		cub::DeviceScan::ExclusiveScan(ts2, tb2, rref, rref, cub::Sum(), data_offset, vars->numElected, streams[1]);
		LASTERR("BVE Phase-2 failed");
		SYNC(streams[0]);
		SYNC(streams[1]);
		if (tmpcap > cumm.scatterCap()) {
			assert(tmpmem != (addr_t)cumm.scatter());
			cacher.deallocate(tmpmem);
		}
	}

	inline void vePhase3(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cuVecB* proof,
		uint32* ucnt,
		uint32* type,
		uint32* rpos,
		S_REF* rref)
	{
		LOGN2(2, "  configuring phase 3 with ");
		grid_t nThreads = BLVE2;
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT2 * sizeof(uint32));
		OPTIMIZEBLOCKSELIM(vars->numElected, nThreads, smemSize, gopts.ve);
		LOGENDING(2, 5, "(%d/%d ths, %d bls) and %zd KB shared memory", nThreads, BLVE2, nBlocks, smemSize / KBYTE);
#if VE_DBG
		ve_k_2 << <1, 1, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->eligible, vars->units, vars->resolved, proof, ucnt, type, rpos, rref);
#else
		ve_k_2 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->eligible, vars->units, vars->resolved, proof, ucnt, type, rpos, rref);
#endif
		if (gopts.sync_always) {
			LASTERR("BVE Phase-3 failed");
			SYNC(0);
		}
	}

	void veAsync(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cudaStream_t* streams,
		cuVecB* proof,
		cuMM& cumm,
		CACHER& cacher,
		const cuHist& cuhist,
		const bool& in)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
		const uint32 max_refs = inf.nDualVars;
		const uint32 max_poss = inf.maxVar;
		const uint32 max_types = inf.maxVar;
		const uint32 max_cnts = inf.maxVar;
		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist;
		uint32* rpos = type + max_types;
		uint32* ucnt = cumm.resizeLits(max_cnts);
		assert(ucnt); // inf.maxVar cannot be larger than nr. of literals
		vePhase1(cnf, ot, vars, proof, ucnt, max_cnts, type, max_types, rpos, max_poss, rref, max_refs, in);
		vePhase2(vars, rpos, rref, streams, cumm, cacher);
		vePhase3(cnf, ot, vars, proof, ucnt, type, rpos, rref);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	}

	void subAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
		LOGN2(2, "  configuring SUB kernel with ");
		grid_t nThreads = BLSUB;
		OPTIMIZESHARED(nThreads, SH_MAX_SUB_IN * sizeof(uint32));
		OPTIMIZEBLOCKSELIM(vars->numElected, nThreads, smemSize, gopts.sub);
		LOGENDING(2, 5, "(%d/%d ths, %d bls) and %zd KB shared memory",
			nThreads, BLSUB, nBlocks, smemSize / KBYTE);
#if SS_DBG
		sub_k << <1, 1, smemSize >> > (cnf, ot, proof, vars->units, vars->elected, vars->eliminated);
#else
		sub_k << <nBlocks, nThreads, smemSize >> > (cnf, ot, proof, vars->units, vars->elected, vars->eliminated);
#endif
		if (gopts.profile_gpu) cutimer->stop(), cutimer->sub += cutimer->gpuTime();
		if (gopts.sync_always) {
			LASTERR("SUB Elimination failed");
			SYNCALL;
		}
	}

	void bceAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
		grid_t nThreads = BLBCE;
		OPTIMIZEBLOCKS(vars->numElected, nThreads, 0);
		bce_k << <nBlocks, nThreads >> > (cnf, ot, proof, vars->resolved, vars->elected, vars->eliminated);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->bce += cutimer->gpuTime();
		if (gopts.sync_always) {
			LASTERR("BCE Elimination failed");
			SYNCALL;
		}
	}

	void ereAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
		dim3 block2D(16, gopts.ere_min_threads);
#else 
		dim3 block2D(devProp.warpSize, gopts.ere_min_threads);
#endif

		SYNCALL;
		dim3 grid2D(1, 1, 1);
		LOGN2(2, "  configuring ERE kernel with ");
		OPTIMIZESHARED(block2D.y, SH_MAX_ERE_OUT * sizeof(uint32));
		OPTIMIZEBLOCKSERE(vars->numElected, block2D, smemSize, gopts.ere);
		grid2D.y = nBlocks;
		LOGENDING(2, 5, "((x: %d/%d, y: %d/%d) ths, (x: %d, y: %d/%d) bls) and %zd KB shared memory",
			block2D.x, devProp.warpSize, block2D.y, gopts.ere_min_threads, 
			grid2D.x, grid2D.y, hardCap, smemSize / KBYTE);
		ere_k << <grid2D, block2D, smemSize >> > (cnf, ot, proof, vars->elected, vars->eliminated);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ere += cutimer->gpuTime();
		if (gopts.sync_always) {
			LASTERR("ERE Elimination failed");
			SYNCALL;
		}
	}

	void Solver::initSharedMem()
	{
#if defined(EXTSHMEM)
		if (devProp.major < 7) LOGERROR("extending shared memory size is not supported");
		LOGN2(2, " Setting maximum shared memory to 64KB..");
		int maxbytes = 65536; // 64 KB
		cudaFuncSetAttribute(in_ve_k_1, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ve_k_1, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ve_k_2, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(sub_k, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ere_k, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		maxGPUSharedMem = maxbytes;
		LOGDONE(2, 5);
#endif
	}

	void Solver::VE()
	{
		if (opts.ve_en) {
			if (interrupted()) killSolver();
			LOG2(2, " Eliminating variables..");
			inf.numDeletedVars = 0;
			veAsync(cnf, ot, vars, streams, cuproof.gpuStream(), cumm, cacher, cuhist, stats.sigma.calls > 1);
			postVE();
			LOGREDALL(this, 2, "BVE Reductions");
		}
	}

	void Solver::postVE()
	{
		size_t bytes = 0;
		uint32* tmpmem = NULL;
		cub::DeviceSelect::If(NULL, bytes, vars->eligible, vars->electedData, tmpmem, vars->numElected, COMPACT_VARS());
		tmpmem = (uint32*)((bytes > cumm.scatterCap()) ? cacher.allocate(bytes) : cumm.scatter());
		if (!tmpmem) throw MEMOUTEXCEPTION();
		cub::DeviceSelect::If(tmpmem + 1, bytes, vars->eligible, vars->electedData, vars->electedSize, vars->numElected, COMPACT_VARS());

		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist, * rpos = type + inf.maxVar;
		veResizeCNFAsync(cnf, rref, type, rpos);

		if (bytes > cumm.scatterCap()) {
			assert(tmpmem != (uint32*)cumm.scatter());
			cacher.deallocate(tmpmem);
		}
	}

	void Solver::SUB()
	{
		if (opts.sub_en || opts.ve_plus_en) {
			if (interrupted()) killSolver();
			LOG2(2, " Eliminating (self)-subsumptions..");
			subAsync(cnf, ot, vars, cuproof.gpuStream());
			LOGREDCL(this, 2, "SUB Reductions");
		}
	}

	void Solver::BCE()
	{
		if (opts.bce_en) {
			if (interrupted()) killSolver();
			if (!vars->numElected) return;
			LOG2(2, " Eliminating blocked clauses..");
			bceAsync(cnf, ot, vars, cuproof.gpuStream());
			LOGREDCL(this, 2, "BCE Reductions");
		}
	}

	void Solver::ERE()
	{
		if (opts.ere_en) {
			if (interrupted()) killSolver();
			if (!vars->numElected) return;
			LOG2(2, " Eliminating redundances..");
			ereCls = inf.numClauses;
			cacheEliminated(streams[5]);
			ereAsync(cnf, ot, vars, cuproof.gpuStream());
			LOGREDCL(this, 2, "ERE Reductions");
			cuproof.cacheProof(0);
			cuproof.writeProof(0);
		}
	}

}