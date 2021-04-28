/***********************************************************************[pftransfer.cu]
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

#include "pfsolve.h"
#include <cub/device/device_select.cuh>
using namespace cub;
using namespace pFROST;
using namespace SIGmA;


void ParaFROST::extract(CNF* dest, BCNF& src)
{
	for (uint32 i = 0; i < src.size(); i++) {
		CLAUSE& c = cm[src[i]];
		if (c.deleted()) continue;
		dest->newClause(c);
		inf.nClauses++, inf.nLiterals += c.size();
	}
}

void ParaFROST::reflectCNF(const cudaStream_t& s1, const cudaStream_t& s2) 
{
	S_REF len1 = hcnf->data().size - dataoff;
	if (!len1) return;
	uint32 len2 = hcnf->size() - csoff;
	CHECK(cudaMemcpyAsync(cumm.cnfMem() + dataoff, hcnf->data().mem + dataoff, len1 * hc_bucket, cudaMemcpyHostToDevice, s1));
	CHECK(cudaMemcpyAsync(cumm.refsMem() + csoff, hcnf->refsData() + csoff, len2 * sizeof(S_REF), cudaMemcpyHostToDevice, s2));
	dataoff = hcnf->data().size, csoff = hcnf->size();
}

void ParaFROST::newBeginning() 
{
	assert(opts.sigma_en || opts.sigma_live_en);
	assert(wt.empty());
	assert(orgs.empty());
	assert(learnts.empty());
	assert(inf.nClauses <= hcnf->size());
	if (unified_access) assert(hcnf == NULL), hcnf = cnf;
	assert(!hcnf->empty());
	cm.init(hcnf->data().size);
	if (!unified_access) {
		sync(streams[0]), sync(streams[1]); // sync CNF caching
		if (profile_gpu) cutimer->stop(streams[1]), cutimer->io += cutimer->gpuTime();
	}
	cacheResolved(streams[2]);
	writeback();
	syncAll();
	if (unified_access) {
		hcnf = NULL;
		if (profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
	}
	else cumm.breakMirror(), hcnf = NULL;
}

inline void ParaFROST::writeback()
{
	stats.literals.original = stats.literals.learnt = 0;
	for (uint32 i = 0; i < hcnf->size(); i++) {
		SCLAUSE& s = hcnf->clause(i);
		if (s.deleted()) continue; //  skip deleted clauses left by 'propFailed'
		newClause(s);
	}
	stats.clauses.original = orgs.size();
	stats.clauses.learnt = learnts.size();
	assert(maxClauses() == int64(inf.nClauses));
}

void ParaFROST::cacheCNF(const cudaStream_t& s1, const cudaStream_t& s2)
{
	// NOTE: if there are units propagated at the last phase,
	// deleted clauses will be left (not compacted), 
	// thus cnf->size() or hcnf->size() must always be used after step 1
	if (interrupted()) killSolver();
	if (sigState == OTALLOC_FAIL) syncAll();
	if (unified_access) {
		if (profile_gpu) cutimer->start();
		if (compacted) countFinal();
		else {
			// 1) compact cs w.r.t clause status on gpu
			uint32* tp = NULL;
			size_t tb = 0;
			// *tp will hold the new size, so tp + 1 will be the start of the temporary array
			DeviceSelect::If(NULL, tb, cnf->refsData(), cnf->refsData(), tp, cnf->size(), COMPACT_CMP(cnf));
			if (tb > cumm.literalsCap()) {
				tp = (uint32*)cumm.allocTemp(tb);
				if (tp == NULL) throw MEMOUTEXCEPTION();
			}
			else tp = cuhist.d_lits;
			DeviceSelect::If(tp + 1, tb, cnf->refsData(), cnf->refsData(), tp, cnf->size(), COMPACT_CMP(cnf), s1);
			countMelted();
			CHECK(cudaMemcpy(&inf.nClauses, tp, sizeof(uint32), cudaMemcpyDeviceToHost));
			if (inf.nClauses) cnf->resize(inf.nClauses);
		}
		// 2) sort cs w.r.t clause size on gpu (user-enabled)
		if (!reallocFailed() && opts.aggr_cnf_sort && inf.nClauses)
			thrust::stable_sort(thrust::cuda::par(tca).on(s1), cnf->refsData(), cnf->refsData() + cnf->size(), CNF_CMP_KEY(cnf));
	}
	else {
		cumm.mirrorCNF(hcnf);
		if (profile_gpu) cutimer->start(s2);
		if (compacted) countFinal();
		else {
			// 1) compact cs w.r.t clause status on gpu
			uint32* tp = NULL;
			size_t tb = 0;
			// *tp will hold the new size, so tp + 1 will be the start of the temporary array
			DeviceSelect::If(NULL, tb, cumm.refsMem(), cumm.refsMem(), tp, hcnf->size(), COMPACT_CMP(cnf));
			if (tb > cumm.literalsCap()) {
				tp = (uint32*)cumm.allocTemp(tb);
				if (tp == NULL) throw MEMOUTEXCEPTION();
			}
			else tp = cuhist.d_lits;
			DeviceSelect::If(tp + 1, tb, cumm.refsMem(), cumm.refsMem(), tp, hcnf->size(), COMPACT_CMP(cnf), s1);
			countMelted();
			CHECK(cudaMemcpy(&inf.nClauses, tp, sizeof(uint32), cudaMemcpyDeviceToHost));
			if (inf.nClauses) hcnf->resize(inf.nClauses);
		}
		if (inf.nClauses) {
			// 2) copy actual cnf data async.
			CHECK(cudaMemcpyAsync(hcnf->data().mem, cumm.cnfMem(), hcnf->data().size * hc_bucket, cudaMemcpyDeviceToHost, s2));
			// 3) sort cs w.r.t clause size on gpu (user-enabled)
			if (!reallocFailed() && opts.aggr_cnf_sort)
				thrust::stable_sort(thrust::cuda::par(tca).on(s1), cumm.refsMem(), cumm.refsMem() + hcnf->size(), CNF_CMP_KEY(cnf));
			// 4) copy compact cs async.
			CHECK(cudaMemcpyAsync(hcnf->refsData(), cumm.refsMem(), hcnf->size() * sizeof(S_REF), cudaMemcpyDeviceToHost, s1));
		}
	}
}