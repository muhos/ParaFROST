/***********************************************************************[elim.cu]
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

#include "solve.h"
using namespace pFROST;

void ParaFROST::VE()
{
	if (opts.ve_en) {
		if (interrupted()) killSolver();
		PFLOG2(2, " Eliminating variables..");
		veAsync(cnf, ot, vars, streams, cuproof.gpuStream(), cumm, cuhist, stats.sigma.calls > 1);
		postVE();
		PFLREDALL(this, 2, "BVE Reductions");
	}
}

void ParaFROST::postVE()
{
	PFLOGN2(2, "  filtering out eliminated variables..");
	int n = 0, lastIdx = -1, len = vars->numPVs;
	uint32* pvs = vars->pVars->data();
	for (int i = 0; i < len; i++) {
		const uint32 x = pvs[i];
		if (ELIMINATED(x)) {
			markEliminated(RECOVERVAR(x));
			if (IS_ADDING(x) && lastIdx < i)
				lastIdx = i;
		}
		else pvs[n++] = x;
	}
	vars->pVars->resize(n);
	PFLENDING(2, 5, "(survived: %d, last index: %d)", n, lastIdx);
	if (!gopts.ve_atomic && lastIdx != -1) {
		PFLOGN2(2, "  resizing CNF to consider added resolvents..");
		assert(n < int(vars->numPVs));
		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist, * rpos = type + inf.maxVar;
		uint32	lastAdded = NOVAR, lastAddedPos = NOVAR;
		S_REF	lastAddedRef = GNOREF;
		CHECK(cudaMemcpyAsync(&lastAdded, type + lastIdx, sizeof(uint32), cudaMemcpyDeviceToHost, streams[0]));
		CHECK(cudaMemcpyAsync(&lastAddedRef, rref + lastIdx, sizeof(S_REF), cudaMemcpyDeviceToHost, streams[1]));
		CHECK(cudaMemcpyAsync(&lastAddedPos, rpos + lastIdx, sizeof(uint32), cudaMemcpyDeviceToHost, streams[2]));
		sync(streams[0]);
		assert(lastAdded < NOVAR);
		assert(RECOVERTYPE(lastAdded) < TYPE_MASK);
		const uint32 lastAddedCls = RECOVERADDEDCLS(lastAdded);
		const uint32 lastAddedLits = RECOVERADDEDLITS(lastAdded);
		assert(lastAddedCls && lastAddedCls <= ADDEDCLS_MAX);
		assert(lastAddedLits && lastAddedLits <= ADDEDLITS_MAX);
		sync(streams[1]);
		assert(lastAddedRef < GNOREF);
		const S_REF lastAddedBuckets = lastAddedLits + hc_nbuckets * lastAddedCls;
		const S_REF data_size = lastAddedBuckets + lastAddedRef;
		sync(streams[2]);
		assert(lastAddedPos < NOVAR);
		const uint32 cs_size = lastAddedCls + lastAddedPos;
		PFLENDING(2, 5, "(new clauses: %d, data: %lld)", cs_size, data_size);
		cumm.resizeCNFAsync(cnf, data_size, cs_size);
	}
	vars->numPVs = n;
}

void ParaFROST::SUB()
{
	if (opts.sub_en || opts.ve_plus_en) {
		if (interrupted()) killSolver();
		PFLOG2(2, " Eliminating (self)-subsumptions..");
		subAsync(cnf, ot, vars, cuproof.gpuStream());
		PFLREDALL(this, 2, "SUB Reductions");
	}
}

void ParaFROST::BCE()
{
	if (opts.bce_en) {
		if (interrupted()) killSolver();
		if (!vars->numPVs) return;
		PFLOG2(2, " Eliminating blocked clauses..");
		bceAsync(cnf, ot, vars, cuproof.gpuStream());
		PFLREDALL(this, 2, "BCE Reductions");
	}
}

void ParaFROST::ERE()
{
	if (opts.ere_en) {
		if (interrupted()) killSolver();
		if (!vars->numPVs) return;
		PFLOG2(2, " Eliminating redundances..");
		ereCls = inf.nClauses;
		ereAsync(cnf, ot, vars, cuproof.gpuStream());
		PFLREDALL(this, 2, "ERE Reductions");
		cuproof.cacheProof(0);
		cuproof.writeProof(0);
	}
}