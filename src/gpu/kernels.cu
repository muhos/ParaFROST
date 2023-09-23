/***********************************************************************[kernels.cu]
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

#include "count.cuh"
#include "bounded.cuh"
#include "subsume.cuh"
#include "blocked.cuh"
#include "simplify.cuh"
#include "redundancy.cuh"
#include "occurrence.cuh"
#include <cub/device/device_scan.cuh>

using namespace cub;

namespace ParaFROST {

	//======================================================//
	//                GPU Wrappers Definitions              //
	//======================================================//
	void initDevOpts()
	{
		CHECK(cudaMemcpyToSymbol(kOpts, &gopts.hostKOpts, sizeof(KOptions), 0, cudaMemcpyHostToDevice));
	}

	void initDevVorg(const cuHist& cuhist)
	{
		DCPTR ptrs = { cuhist.d_vorg, cuhist.d_lbyte };
		CHECK(cudaMemcpyToSymbol(DC_PTRS, &ptrs, sizeof(DCPTR), 0, cudaMemcpyHostToDevice));
	}

	void initSharedMem()
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

	void printConstants()
	{
		printCMem << <1, 1 >> > ();
		LASTERR("Printing constant memory failed");
		SYNCALL;
	}

	void cuMemSetAsync(addr_t mem, const Byte& val, const size_t& size)
	{
		OPTIMIZEBLOCKS(uint32(size), BLOCK1D);
		memset_k<Byte> << <nBlocks, BLOCK1D >> > (mem, val, size);
		if (gopts.sync_always) {
			LASTERR("CUDA memory set failed");
			SYNCALL;
		}
	}

	void prepareCNFAsync(CNF* cnf, const cudaStream_t& _s)
	{
		assert(inf.nClauses);
		if (gopts.profile_gpu) cutimer->start(_s);
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		prep_cnf_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf);
		if (gopts.sync_always) {
			LASTERR("Preparing CNF failed");
			SYNCALL;
		}
		if (gopts.profile_gpu) cutimer->stop(_s), cutimer->sig += cutimer->gpuTime();
	}

	void mapFrozenAsync(VARS* vars, const uint32& size)
	{
		assert(vars->varcore == vars->eligible); // an alies of eligible
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(size, BLOCK1D);
		// 'vars->scores' is an alies for frozen vars on the GPU side
		mapfrozen_k << <nBlocks, BLOCK1D >> > (vars->scores, vars->varcore, size);
		if (gopts.sync_always) {
			LASTERR("Mapping frozen failed");
			SYNCALL;
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	}
	//=======================
	// histogram related
	//=======================
	void copyIf(uint32* dest, CNF* src)
	{
		if (gopts.profile_gpu) cutimer->start();
		reset_counter << <1, 1 >> > ();
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		copy_if_k << <nBlocks, BLOCK1D >> > (dest, src);
#if defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
		check_counter << <1, 1 >> > (inf.nLiterals);
#endif
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LASTERR("Copying literals failed");
		SYNCALL;
	}

	void copyIfAsync(uint32* dest, CNF* src)
	{
		if (gopts.profile_gpu) cutimer->start();
		reset_counter << <1, 1 >> > ();
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		copy_if_k << <nBlocks, BLOCK1D >> > (dest, src);
		if (gopts.sync_always) {
			check_counter << <1, 1 >> > (inf.nLiterals);
			LASTERR("Copying literals failed");
			SYNCALL;
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
	}

	void calcScores(VARS* vars, uint32* hist)
	{
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
		assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, inf.maxVar);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LASTERR("Assigning scores failed");
		SYNCALL;
	}

	void calcScores(VARS* vars, uint32* hist, OT* ot)
	{
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
		assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LASTERR("Assigning scores failed");
		SYNCALL;
	}
	//=======================
	// CNF measurements
	//=======================
	template <typename D, typename I>
	inline D seqreduceBlocks(const D* blocks, const I& n)
	{
		D finalcount = 0;
		for (I i = 0; i < n; ++i)
			finalcount += blocks[i];
		return finalcount;
	}

	template <typename D, typename I>
	inline void seqreduceBlocks(const D* CBlocks, const D* LBlocks, const I& n)
	{
		inf.n_cls_after = 0;
		inf.n_lits_after = 0;
		for (I i = 0; i < n; ++i) {
			inf.n_cls_after += CBlocks[i];
			inf.n_lits_after += LBlocks[i];
		}
	}

	uint32 cuPROOF::count(const uint32* literals, const uint32& numLits)
	{
		if (!proof.checkFile()) LOGERROR("host proof system is not activated");
		if (!literals) return 0;
		if (!numLits) return 0;
		enabled = true;
		LOGN2(2, " Counting proof bytes..");
		OPTIMIZEBLOCKS2(numLits, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		SYNCALL; // sync any pending kernels or transfers
		if (gopts.profile_gpu) cutimer->start();
		cnt_proof << <nBlocks, BLOCK1D, smemSize >> > (literals, numLits);
		LASTERR("Proof counting failed");
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
		const uint32 maxcap = seqreduceBlocks(hostLBlocks, nBlocks);
		assert(maxcap && maxcap < (numLits * sizeof(uint32)));
		LOGENDING(2, 5, "(%d bytes)", maxcap);
		return maxcap;
	}

	void parcountCls(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		cnt_cls << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		inf.n_cls_after = seqreduceBlocks(hostCBlocks, nBlocks);
	}

	void parcountLits(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		inf.n_lits_after = seqreduceBlocks(hostLBlocks, nBlocks);
	}

	void parcountAll(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32) * 2);
		cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		seqreduceBlocks(hostCBlocks, hostLBlocks, nBlocks);
	}
	//=======================
	// occurrence table
	//=======================
	void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
			LASTERR("Occurrence table reduction failed");
			SYNCALL;
			if (p) {
				LOGRULER('=', 30);
				LOG0("\toccurrence table");
				ot->print();
				LOGRULER('=', 30);
			}
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->rot += cutimer->gpuTime();
	}

	inline void resetOTAsync(CNF* cnf, OT* ot)
	{
		assert(cnf);
		assert(ot);
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reset_ot_k << <nBlocks, BLOCK1D >> > (ot);
		if (gopts.sync_always) {
			LASTERR("Occurrence table reset failed");
			SYNCALL;
			assert(ot->accViolation(inf.maxVar));
		}
	}

	void createOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		resetOTAsync(cnf, ot);
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		create_ot_k << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
			LASTERR("Occurrence table creation failed");
			SYNCALL;
			assert(ot->accViolation(inf.maxVar));
			if (p) {
				LOGRULER('=', 30);
				LOG0("\toccurrence table");
				ot->print();
				LOGRULER('=', 30);
			}
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->cot += cutimer->gpuTime();
	}
	//=======================
	// variable elimination
	//=======================
	inline void vePhase1(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cuVecB* proof,
		uint32* ucnt,
		uint32* type,
		uint32* rpos,
		S_REF* rref,
		const bool& in)
	{
		LOGN2(2, "  resetting last eliminated ID to -1.. ");
		reset_id << <1, 1 >> > ();
		LOGDONE(2, 5);
		LOGN2(2, "  configuring phase 1 with ");
		OPTIMIZEBLOCKSELIM(vars->numElected, BLVE1, gopts.ve);
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT1 * sizeof(uint32));
		LOGENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLVE1, nBlocks, MAXBLOCKS, smemSize / KBYTE);
#if VE_DBG
		if (in) in_ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
		else	   ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
#else
		if (in) in_ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
		else	   ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->elected, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
#endif
		LASTERR("BVE Phase-1 failed");
		SYNC(0);
	}

	inline void vePhase2(
		VARS* vars,
		uint32* rpos,
		S_REF* rref,
		cudaStream_t* streams,
		cuMM& cumm)
	{
		const uint32 cs_offset = cumm.pinnedCNF()->size();
		const S_REF data_offset = cumm.pinnedCNF()->data().size;
		size_t tb1 = 0, tb2 = 0;
		DeviceScan::ExclusiveScan(NULL, tb1, rpos, rpos, Sum(), cs_offset, vars->numElected);
		DeviceScan::ExclusiveScan(NULL, tb2, rref, rref, Sum(), data_offset, vars->numElected);
		size_t tmpcap = tb1 + tb2;
		addr_t ts1 = NULL, ts2 = NULL;
		addr_t tmpmem = (addr_t) ((tmpcap > cumm.scatterCap()) ? cacher.allocate(tmpcap) : cumm.scatter());
		ts1 = tmpmem, ts2 = ts1 + tb1;
		DeviceScan::ExclusiveScan(ts1, tb1, rpos, rpos, Sum(), cs_offset, vars->numElected, streams[0]);
		DeviceScan::ExclusiveScan(ts2, tb2, rref, rref, Sum(), data_offset, vars->numElected, streams[1]);
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
		OPTIMIZEBLOCKSELIM(vars->numElected, BLVE2, gopts.ve);
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT2 * sizeof(uint32));
		LOGENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLVE2, nBlocks, MAXBLOCKS, smemSize / KBYTE);
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
		const cuHist& cuhist,
		const bool& in)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist;
		uint32* rpos = type + inf.maxVar;
		uint32* ucnt = cumm.resizeLits(inf.maxVar);
		assert(ucnt); // inf.maxVar cannot be larger than nr. of literals
		vePhase1(cnf, ot, vars, proof, ucnt, type, rpos, rref, in);
		vePhase2(vars, rpos, rref, streams, cumm);
		vePhase3(cnf, ot, vars, proof, ucnt, type, rpos, rref);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	}

	void veResizeCNFAsync(CNF* cnf, const cuHist& cuhist)
	{
		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist, * rpos = type + inf.maxVar;
		resizeCNF_k << <1, 1 >> > (cnf, type, rpos, rref, verbose);
		if (gopts.sync_always) {
			LASTERR("Resizing CNF after BVE failed");
			SYNC(0);
		}
	}

	void subAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numElected);
		if (gopts.profile_gpu) cutimer->start();
		LOGN2(2, "  configuring SUB kernel with ");
		OPTIMIZEBLOCKSELIM(vars->numElected, BLSUB, gopts.sub);
		OPTIMIZESHARED(nThreads, SH_MAX_SUB_IN * sizeof(uint32));
		LOGENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLSUB, nBlocks, MAXBLOCKS, smemSize / KBYTE);
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
		OPTIMIZEBLOCKS(vars->numElected, BLBCE);
		bce_k << <nBlocks, BLBCE >> > (cnf, ot, proof, vars->resolved, vars->elected, vars->eliminated);
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
		dim3 block2D(16, devProp.warpSize);
#else 
		dim3 block2D(devProp.warpSize, devProp.warpSize);
#endif
        SYNCALL;
		dim3 grid2D(1, 1, 1);
		LOGN2(2, "  configuring ERE kernel with ");
		OPTIMIZEBLOCKSERE(vars->numElected, block2D, gopts.ere);
		OPTIMIZESHARED(block2D.y, SH_MAX_ERE_OUT * sizeof(uint32));
		grid2D.y = nBlocks;
		LOGENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			block2D.y, devProp.warpSize, grid2D.y, MAXBLOCKS, smemSize / KBYTE);
		ere_k << <grid2D, block2D, smemSize >> > (cnf, ot, proof, vars->elected, vars->eliminated);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ere += cutimer->gpuTime();
		if (gopts.sync_always) {
			LASTERR("ERE Elimination failed");
			SYNCALL;
		}
	}

}