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
	// thus cnf->size() or hcnf->size() must always be used
	if (interrupted()) killSolver();
	if (sigState == OTALLOC_FAIL) syncAll();
	cudaStream_t copystream;
	if (unified_access) copystream = 0, hcnf = cnf;
	else copystream = s2, cumm.mirrorCNF(hcnf);
	assert(hcnf);
	if (profile_gpu) cutimer->start(copystream);
	if (compacted) countFinal();
	else {
		size_t bytes = 0;
		S_REF* tmp = NULL;
		// *tmp will hold the new size, so tmp + 1 will be the start of the temporary array
		DeviceSelect::If(NULL, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), COMPACT_CMP(cnf));
		tmp = cumm.scatter();
		if (bytes > cumm.scatterCap()) tmp = (S_REF*)cumm.allocTemp(bytes);
		if (!tmp) throw MEMOUTEXCEPTION();
		DeviceSelect::If(tmp + 1, bytes, cumm.refsMem(), cumm.refsMem(), tmp, hcnf->size(), COMPACT_CMP(cnf), s1);
		countMelted();
		CHECK(cudaMemcpy(&inf.nClauses, tmp, sizeof(uint32), cudaMemcpyDeviceToHost));
		if (inf.nClauses) hcnf->resize(inf.nClauses);
	}
	if (inf.nClauses) {
		if (!unified_access) 
			CHECK(cudaMemcpyAsync(hcnf->data().mem, cumm.cnfMem(), hcnf->data().size * hc_bucket, cudaMemcpyDeviceToHost, s2));
		if (!reallocFailed() && opts.aggr_cnf_sort)
			thrust::stable_sort(thrust::cuda::par(tca).on(s1), cumm.refsMem(), cumm.refsMem() + hcnf->size(), CNF_CMP_KEY(cnf));
		if (!unified_access) 
			CHECK(cudaMemcpyAsync(hcnf->refsData(), cumm.refsMem(), hcnf->size() * sizeof(S_REF), cudaMemcpyDeviceToHost, s1));
	}
}