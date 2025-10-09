/***********************************************************************[cnf.cu]
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
#include "options.cuh"
#include "timer.cuh"
#include "count.cuh"
#include "grid.cuh"
#include "cnf.cuh"
#include "primitives.cuh"
#include <cub/device/device_select.cuh>
#include <thrust/sort.h>

namespace ParaFROST {

	__device__ int lastEliminatedID;

	__global__
	void copy_if_k(uint32* __restrict__ dest, CNF* __restrict__ src)
	{
		for_parallel_x (tid, src->size()) {
			SCLAUSE& c = src->clause(tid);
			if (c.original() || c.learnt()) {
				uint32* d = dest + atomicAdd(&gcounter, c.size());
				forall_clause(c, s) { *d++ = *s; }
			}
		}
	}

	__global__ 
	void prep_cnf_k(CNF* cnf)
	{
		for_parallel_x (tid, cnf->size()) {
			SCLAUSE& c = cnf->clause(tid);
			devSort(c.data(), c.size());
			calcSig(c);
		}
	}

	__global__ 
	void resizeCNF_k(CNF* cnf,
		const uint32* __restrict__ type,
		const uint32* __restrict__ rpos,
		const S_REF* __restrict__ rref,
		const int verbose)
	{
		if (lastEliminatedID >= 0) {
			const uint32 lastAdded = type[lastEliminatedID];
			const uint32 lastAddedPos = rpos[lastEliminatedID];
			const S_REF  lastAddedRef = rref[lastEliminatedID];
			assert(lastAdded < NOVAR);
			assert(lastAddedPos < NOVAR);
			assert(lastAddedRef < GNOREF);
			assert(RECOVERTYPE(lastAdded) <= TYPE_MASK);
			const uint32 lastAddedCls = RECOVERADDEDCLS(lastAdded);
			const uint32 lastAddedLits = RECOVERADDEDLITS(lastAdded);
			assert(lastAddedCls && lastAddedCls <= ADDEDCLS_MAX);
			assert(lastAddedLits && lastAddedLits <= ADDEDLITS_MAX);
			const S_REF lastAddedBuckets = lastAddedLits + DC_NBUCKETS * lastAddedCls;
			const S_REF data_size = lastAddedBuckets + lastAddedRef;
			const uint32 cs_size = lastAddedCls + lastAddedPos;
			cnf->resize(data_size, cs_size);
			if (verbose > 1) printf("c   resized CNF to %d clauses and %lld data for a last ID %d\n", cs_size, data_size, lastEliminatedID);
		}
	}

	void copyIf(uint32* dest, CNF* src)
	{
		reset_counter << <1, 1 >> > ();
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.numClauses, nThreads, 0);
		copy_if_k << <nBlocks, nThreads >> > (dest, src);
#if defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
		check_counter << <1, 1 >> > (inf.numLiterals);
#endif
		LASTERR("Copying literals failed");
		SYNCALL;
	}

	void copyIfAsync(uint32* dest, CNF* src)
	{
		reset_counter << <1, 1 >> > ();
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.numClauses, nThreads, 0);
		copy_if_k << <nBlocks, nThreads >> > (dest, src);
		if (gopts.sync_always) {
			check_counter << <1, 1 >> > (inf.numLiterals);
			LASTERR("Copying literals failed");
			SYNCALL;
		}
	}

	void prepareCNFAsync(CNF* cnf, const cudaStream_t& _s)
	{
		assert(inf.numClauses);
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.numClauses, nThreads, 0);
		prep_cnf_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf);
		if (gopts.sync_always) {
			LASTERR("Preparing CNF failed");
			SYNCALL;
		}
	}

	void veResizeCNFAsync(CNF* cnf, S_REF* rref, uint32* type, uint32* rpos)
	{
		resizeCNF_k << <1, 1 >> > (cnf, type, rpos, rref, verbose);
		if (gopts.sync_always) {
			LASTERR("Resizing CNF after BVE failed");
			SYNC(0);
		}
	}

	bool Solver::reallocCNF(const bool& realloc)
	{
		if (realloc) {
			size_t maxAddedCls = opts.ve_en ? inf.numClauses : 0;
			size_t maxAddedLits = opts.ve_en ? size_t(stats.literals.original * opts.lits_mul) : 0;
			LOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
			if (!cumm.resizeCNF(cnf, inf.numClauses + maxAddedCls, inf.numLiterals + maxAddedLits)) {
				simpstate = CNFALLOC_FAIL, compacted = false;
				return false;
			}
			compacted = true;
			olist_cmp.init(cnf), compact_cmp.init(cnf);
		}
		else cumm.cacheCNFPtr(cnf), compacted = false;
		return true;
	}

	bool Solver::reallocCNF()
	{
		const int times = phase + 1;
		return reallocCNF(times > 1 && times != opts.phases && (times % opts.shrink_rate) == 0);
	}

	uint32* Solver::flattenCNF(const uint32& numLits)
	{
		assert(numLits);
		uint32* literals = cumm.resizeLits(numLits);
		if (flattened || !literals) return literals;
		LOGN2(2, " Copying survived literals..");
		if (gopts.profile_gpu) cutimer.start();
		copyIfAsync(literals, cnf);
		if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.vo += cutimer.gpuTime();
		LOGENDING(2, 5, "(%d copied)", numLits);
		flattened = true;
		return literals;
	}

	void Solver::reflectCNF(const cudaStream_t& s1, const cudaStream_t& s2)
	{
		S_REF len1 = hcnf->data().size - dataoff;
		if (!len1) return;
		uint32 len2 = hcnf->size() - csoff;
		CHECK(cudaMemcpyAsync(cumm.cnfMem() + dataoff, hcnf->data().mem + dataoff, len1 * SBUCKETSIZE, cudaMemcpyHostToDevice, s1));
		CHECK(cudaMemcpyAsync(cumm.refsMem() + csoff, hcnf->refsData() + csoff, len2 * sizeof(S_REF), cudaMemcpyHostToDevice, s2));
		dataoff = hcnf->data().size, csoff = hcnf->size();
	}

	void Solver::extractCNF(CNF* dest, BCNF& src)
	{
		for (uint32 i = 0; i < src.size(); i++) {
			CLAUSE& c = cm[src[i]];
			if (c.deleted()) continue;
			dest->newClause(c);
			inf.numClauses++, inf.numLiterals += c.size();
		}
	}

	void Solver::writeBackCNF()
	{
		int64 bliterals = maxLiterals();
		stats.literals.original = stats.literals.learnt = 0;
		for (uint32 i = 0; i < hcnf->size(); i++) {
			SCLAUSE& s = hcnf->clause(i);
			if (s.deleted()) continue; //  skip deleted clauses left by 'propFailed'
			newClause(s);
		}
		stats.clauses.original = orgs.size();
		stats.clauses.learnt = learnts.size();
		assert(maxClauses() == int64(inf.numClauses));
	}

	void Solver::cacheCNF(const cudaStream_t& s1, const cudaStream_t& s2)
	{
		// NOTE: if there are units propagated at the last phase,
		// deleted clauses will be left (not compacted), 
		// thus cnf->size() or hcnf->size() must always be used
		if (interrupted()) killSolver();
		if (simpstate == OTALLOC_FAIL) SYNCALL;
		cudaStream_t copystream;
		if (gopts.unified_access) copystream = 0, hcnf = cnf;
		else copystream = s2, cumm.mirrorCNF(hcnf);
		assert(hcnf);
		if (gopts.profile_gpu) cutimer.start(copystream);
		if (compacted) {
			countCls();
			inf.numClauses = inf.numClausesSurvived;
		}
		else {
			size_t bytes = 0;
			S_REF* tmp = NULL;
			cub::DeviceSelect::If(NULL, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), compact_cmp);
			tmp = (bytes > cumm.scatterCap()) ? (S_REF*)cacher.allocate(bytes) : cumm.scatter();
			assert(tmp);
			// *tmp will hold the new size, so tmp + 1 will be the start of the temporary array
			cub::DeviceSelect::If(tmp + 1, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), compact_cmp, s1);
			CHECK(cudaMemcpy(&inf.numClauses, tmp, sizeof(uint32), cudaMemcpyDeviceToHost));
			if (inf.numClauses) hcnf->resize(inf.numClauses);
			if (bytes > cumm.scatterCap()) {
				assert(tmp != cumm.scatter());
				cacher.deallocate(tmp);
			}
		}
		if (inf.numClauses) {
			if (!gopts.unified_access)
				CHECK(cudaMemcpyAsync(hcnf->data().mem, cumm.cnfMem(), hcnf->data().size * SBUCKETSIZE, cudaMemcpyDeviceToHost, s2));
			if (!reallocFailed() && opts.aggr_cnf_sort)
				thrust::stable_sort(thrust::cuda::par(tca).on(s1), cumm.refsMem(), cumm.refsMem() + hcnf->size(), olist_cmp);
			if (!gopts.unified_access)
				CHECK(cudaMemcpyAsync(hcnf->refsData(), cumm.refsMem(), hcnf->size() * sizeof(S_REF), cudaMemcpyDeviceToHost, s1));
		}
	}

	_PFROST_D_ S_REF* CNF::jump(S_REF& ref, const uint32& nCls, const uint32& nLits)
	{
		assert(nLits >= nCls);
		const S_REF regionSize = nLits + DC_NBUCKETS * nCls;
		ref = atomicAdd(&_data.size, regionSize);
		assert(ref < _data.cap);
		return _refs.jump(nCls);
	}

}