/***********************************************************************[pfelim.cu]
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
#include "pfsort.h"

namespace pFROST {
	using namespace SIGmA;
	
	void ParaFROST::createOTHost(HOT& hot)
	{
		assert(hcnf != NULL);
		assert(vars->nUnits);
		assert(reallocFailed());
		for (uint32 i = 0; i < hcnf->size(); i++) {
			S_REF r = hcnf->ref(i);
			SCLAUSE& c = (*hcnf)[r];
			if (c.learnt() || c.original()) {
				assert(c.size());
				for (uint32* k = c; k != c.end(); k++) {
					assert(*k > 1);
					hot[*k].push(r);
				}
			}
		}
	}

	bool ParaFROST::propHost()
	{
		if (vars->nUnits) {
			assert(reallocFailed());
			assert(hcnf != NULL);
			assert(inf.nClauses == hcnf->size());
			HOT hot(inf.nDualVars);
			createOTHost(hot);
			// start proping on host
			nForced = sp->propagated;
			assert(vars->cachedUnits != NULL);
			uint32* t = vars->cachedUnits + vars->nUnits;
			for (uint32* u = vars->cachedUnits; u != t; u++) {
				LIT_ST val = value(*u);
				if (UNASSIGNED(val)) enqueueOrg(*u);
				else if (!val) return false; // early conflict detection
			}
			if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
			else PFLTRAIL(this, 3);
			while (sp->propagated < trail.size()) { // propagate units
				uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
				assert(assign > 1);
				HOL& hol = hot[assign], &f_hol = hot[f_assign];
				for (S_REF* i = hol; i != hol.end(); i++) (*hcnf)[*i].markDeleted();
				for (S_REF* i = f_hol; i != f_hol.end(); i++) { 
					SCLAUSE& c = (*hcnf)[*i];
					assert(c.size());
					if (c.deleted() || propClause(c, f_assign)) continue;
					assert(c.size()); // cannot be empty at this point
					if (c.size() == 1) {
						assert(*c > 1);
						if (unassigned(*c)) enqueueOrg(*c);
						else { cnfstate = UNSAT; return false; } 
					}
				}
			}
			nForced = sp->propagated - nForced;
			PFLREDALLHOST(this, 2, "BCP Reductions");
			vars->nUnits = nForced = 0;
			countCls(1), inf.nClauses = inf.n_cls_after; // only needed for writeBack assertion
		}
		return true;
	}

	bool ParaFROST::prop()
	{
		if (!enqeueCached(streams[3])) { cnfstate = UNSAT; return false; }
		while (sp->propagated < trail.size()) { // propagate units
			uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
			assert(assign > 1);
			OL& ol = (*ot)[assign], & f_ol = (*ot)[f_assign];
			for (S_REF* i = ol; i != ol.end(); i++) (*cnf)[*i].markDeleted(); // remove satisfied
			for (S_REF* i = f_ol; i != f_ol.end(); i++) { // reduce unsatisfied 
				SCLAUSE& c = (*cnf)[*i];
				assert(c.size());
				if (c.deleted() || propClause(c, f_assign)) continue; // clause satisfied
				assert(c.size()); // cannot be empty at this point
				if (c.size() == 1) {
					assert(*c > 1);
					if (unassigned(*c)) enqueueOrg(*c);
					else { cnfstate = UNSAT; return false; }  // conflict on top level
				}
			}
			(*ot)[assign].clear(true), (*ot)[f_assign].clear(true);
		}
		cleanProped();
		return true;
	}

	void ParaFROST::VE()
	{
		if (opts.ve_en) {
			if (interrupted()) killSolver();
			PFLOGN2(2, "  Eliminating variables..");
			veAsync(cnf, ot, vars, streams, cumm, cuhist, stats.sigmifications);
			postVE();
			PFLDONE(2, 5);
			PFLREDALL(this, 2, "BVE Reductions");
		}
	}

	void ParaFROST::postVE() 
	{
		int n = 0, lastIdx = -1, len = vars->numPVs;
		uint32 *pvs = vars->pVars->data();
		for (int i = 0; i < len; i++) {
			uint32 x = pvs[i];
			if (ELIMINATED(x)) {
				if (lastIdx < i) lastIdx = i;
				sp->vstate[RECOVERVAR(x)] = MELTED;
			}
			else pvs[n++] = x;
		}
		vars->pVars->resize(n);
		if (!atomic_ve && lastIdx != -1) {
			assert(n < int(vars->numPVs));
			uint32* type = cuhist.d_hist, * rpos = type + inf.maxVar;
			uint32 lastAdded = NOVAR, lastAddedPos = NOVAR;
			CHECK(cudaMemcpyAsync(&lastAdded, type + lastIdx, sizeof(uint32), cudaMemcpyDeviceToHost, streams[0]));
			CHECK(cudaMemcpyAsync(&lastAddedPos, rpos + lastIdx, sizeof(uint32), cudaMemcpyDeviceToHost, streams[1]));
			sync(streams[0]), assert(lastAdded < NOVAR), lastAdded = RECOVERADDED(lastAdded), sync(streams[1]);
			assert(lastAddedPos < NOVAR);
			uint32 cs_size = lastAdded + lastAddedPos;
			S_REF lastRef = GNOREF;
			SCLAUSE lastClause;
			CHECK(cudaMemcpy(&lastRef, cumm.csMem() + cs_size - 1, sizeof(S_REF), cudaMemcpyDeviceToHost)), assert(lastRef < GNOREF);
			CHECK(cudaMemcpy(&lastClause, cumm.cnfMem() + lastRef, sizeof(SCLAUSE), cudaMemcpyDeviceToHost)), assert(lastClause.blockSize());
			S_REF data_size = lastRef + lastClause.blockSize();
			cumm.resizeCNFAsync(cnf, data_size, cs_size);
		}
		vars->numPVs = n;
	}

	void ParaFROST::HSE()
	{
		if (opts.hse_en || opts.ve_plus_en) {
			if (interrupted()) killSolver();
			PFLOGN2(2, "  Eliminating (self)-subsumptions..");
			hseAsync(cnf, ot, vars);
			PFLDONE(2, 5);
			PFLREDALL(this, 2, "HSE Reductions");
		}
	}

	void ParaFROST::BCE()
	{
		if (opts.bce_en) {
			if (interrupted()) killSolver();
			if (!vars->numPVs) return;
			PFLOGN2(2, " Eliminating blocked clauses..");
			bceAsync(cnf, ot, vars);
			PFLDONE(2, 5);
			PFLREDALL(this, 2, "BCE Reductions");
		}
	}

	void ParaFROST::ERE()
	{
		if (opts.ere_en) {
			if (interrupted()) killSolver();
			if (!vars->numPVs) return;
			PFLOGN2(2, " Eliminating redundances..");
			ereCls = inf.nClauses;
			ereAsync(cnf, ot, vars);
			PFLDONE(2, 5);
			PFLREDALL(this, 2, "ERE Reductions");
		}
	}

	inline bool ParaFROST::propClause(SCLAUSE& c, const uint32& f_assign)
	{
		assert(c.size() > 1);
		uint32 sig = 0;
		uint32* i, * j, * end = c.end();
		for (i = c, j = i; i != end; i++) {
			uint32 lit = *i;
			if (lit != f_assign) {
				if (isTrue(lit)) return true;
				*j++ = lit;
				sig |= MAPHASH(lit);
			}
		}
		assert(j - c == c.size() - 1);
		assert(c.hasZero() < 0);
		assert(c.isSorted());
		c.set_sig(sig);
		c.pop();
		return false;
	}

	inline bool ParaFROST::enqeueCached(const cudaStream_t& stream) {
		if (vars->nUnits) {
			nForced = sp->propagated;
			sync(stream); // sync units copy
			assert(vars->cachedUnits != NULL);
			uint32* t = vars->cachedUnits + vars->nUnits;
			for (uint32* u = vars->cachedUnits; u != t; u++) {
				LIT_ST val = value(*u);
				if (UNASSIGNED(val)) enqueueOrg(*u);
				else if (!val) return false; // early conflict detection
			}
			if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
			else PFLTRAIL(this, 3);
			syncAll(); // sync ot creation
		}
		return true;
	}

	inline void	ParaFROST::cleanProped() {
		if (vars->nUnits) {
			nForced = sp->propagated - nForced;
			PFLREDALL(this, 2, "BCP Reductions");
			nForced = 0, vars->tmpObj.clear();
			assert(vars->tmpObj.data() == cumm.unitsdPtr());
			CHECK(cudaMemcpyAsync(vars->units, &vars->tmpObj, sizeof(cuVecU), cudaMemcpyHostToDevice));
		}
	}

}