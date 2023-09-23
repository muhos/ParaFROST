/***********************************************************************[transfer.cu]
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

#include "solve.hpp"
#include <cub/device/device_select.cuh>
using namespace cub;
using namespace ParaFROST;

void Solver::extract(CNF* dest, BCNF& src)
{
	for (uint32 i = 0; i < src.size(); i++) {
		CLAUSE& c = cm[src[i]];
		if (c.deleted()) continue;
		dest->newClause(c);
		inf.nClauses++, inf.nLiterals += c.size();
	}
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

void Solver::newBeginning() 
{
	assert(opts.sigma_en || opts.sigma_live_en);
	assert(wt.empty());
	assert(orgs.empty());
	assert(learnts.empty());
	assert(inf.nClauses <= hcnf->size());
	if (gopts.unified_access) assert(hcnf == NULL), hcnf = cnf;
	assert(!hcnf->empty());
	cm.init(hcnf->data().size);
	if (!gopts.unified_access) {
		SYNC(streams[0]); SYNC(streams[1]); // sync CNF caching
		if (gopts.profile_gpu) cutimer->stop(streams[1]), cutimer->io += cutimer->gpuTime();
	}
	cacheResolved(streams[2]);
	writeBack();
	SYNCALL;
	if (gopts.unified_access) {
		hcnf = NULL;
		if (gopts.profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
	}
	else cumm.breakMirror(), hcnf = NULL;
}

void Solver::markEliminated(const cudaStream_t& _s)
{
	assert(vars->isEliminatedCached);

#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
	uint32 unassigned = inf.unassigned;
#endif
	
	SYNC(_s);

	forall_variables(v) {
		if (vars->cachedEliminated[v])
			markEliminated(v);
	}

#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
	assert(unassigned >= inf.unassigned);
	assert((unassigned - inf.unassigned) == vars->nMelted);
#endif
}

inline void Solver::writeBack()
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
	stats.sigma.all.literals += bliterals - maxLiterals();
	assert(maxClauses() == int64(inf.nClauses));
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
	if (gopts.profile_gpu) cutimer->start(copystream);
	if (compacted) {
		countCls();
		inf.nClauses = inf.n_cls_after;
	}
	else {
		size_t bytes = 0;
		S_REF* tmp = NULL;
		DeviceSelect::If(NULL, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), compact_cmp);
		tmp = (bytes > cumm.scatterCap()) ? (S_REF*)cacher.allocate(bytes) : cumm.scatter();
		assert(tmp);
		// *tmp will hold the new size, so tmp + 1 will be the start of the temporary array
		DeviceSelect::If(tmp + 1, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), compact_cmp, s1);
		CHECK(cudaMemcpy(&inf.nClauses, tmp, sizeof(uint32), cudaMemcpyDeviceToHost));
		if (inf.nClauses) hcnf->resize(inf.nClauses);
		if (bytes > cumm.scatterCap()) {
			assert(tmp != cumm.scatter());
			cacher.deallocate(tmp);
		}
	}
	if (inf.nClauses) {
		if (!gopts.unified_access) 
			CHECK(cudaMemcpyAsync(hcnf->data().mem, cumm.cnfMem(), hcnf->data().size * SBUCKETSIZE, cudaMemcpyDeviceToHost, s2));
		if (!reallocFailed() && opts.aggr_cnf_sort)
			thrust::stable_sort(thrust::cuda::par(tca).on(s1), cumm.refsMem(), cumm.refsMem() + hcnf->size(), olist_cmp);
		if (!gopts.unified_access) 
			CHECK(cudaMemcpyAsync(hcnf->refsData(), cumm.refsMem(), hcnf->size() * sizeof(S_REF), cudaMemcpyDeviceToHost, s1));
	}
}

void Solver::cacheUnits(const cudaStream_t& stream) 
{
	SYNC(stream);
	if ((vars->nUnits = vars->tmpUnits.size()))
		CHECK(cudaMemcpyAsync(vars->cachedUnits, cumm.unitsdPtr(), vars->nUnits * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}

void Solver::cacheEliminated(const cudaStream_t& stream)
{
	if (vars->isEliminatedCached) return;
	CHECK(cudaMemcpyAsync(vars->cachedEliminated, vars->eliminated, inf.maxVar + 1, cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
	vars->isEliminatedCached = true;
}

void Solver::cacheNumUnits(const cudaStream_t& stream)
{
	CHECK(cudaMemcpyAsync(&vars->tmpUnits, vars->units, sizeof(cuVecU), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}

void Solver::cacheResolved(const cudaStream_t& stream)
{
	cuVecU tmpObj;
	CHECK(cudaMemcpy(&tmpObj, vars->resolved, sizeof(cuVecU), cudaMemcpyDeviceToHost));
	const uint32* devStart = tmpObj.data();
	const uint32 devSize = tmpObj.size();
	if (!devSize) return;
	assert(devStart);
	const uint32 off = model.resolved.size();
	model.resolved.resize(off + devSize);
	uint32* start = model.resolved + off;
	CHECK(cudaMemcpyAsync(start, devStart, devSize * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}