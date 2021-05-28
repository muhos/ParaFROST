/***********************************************************************[lcve.cu]
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

void ParaFROST::varReorder()
{
	PFLOGN2(2, " Finding eligible variables for LCVE..");
	assert(cuhist.d_hist != NULL);
	// NOTE: OT creation will be synced in calcScores call
	if (vars->nUnits) calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
	else calcScores(vars, cuhist.d_hist);
	cuhist.cacheHist(streams[2]);
	if (profile_gpu) cutimer->start(streams[3]);
	thrust::sort(thrust::cuda::par(tca).on(streams[3]), vars->eligible, vars->eligible + inf.maxVar, GPU_LCV_CMP(vars->scores));
	PFLDONE(2, 5);
	vars->nUnits = 0;
	sync(streams[2]);
	if (profile_gpu) cutimer->stop(streams[3]), cutimer->vo += cutimer->gpuTime();
	if (verbose == 4) {
		PFLOG0(" Eligible variables:");
		for (uint32 v = 0; v < inf.maxVar; v++) {
			uint32 x = vars->eligible[v], p = V2L(x), n = NEG(p);
			PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, cuhist[p], cuhist[n], vars->scores[x]);
		}
	}
}

bool ParaFROST::LCVE()
{
	// reorder variables
	varReorder();
	// extended LCVE
	PFLOGN2(2, " Electing variables (p-mu: %d, n-mu: %d)..", opts.mu_pos << mu_inc, opts.mu_neg << mu_inc);
	vars->numPVs = 0, vars->pVars->clear();
	sp->stacktail = sp->tmp_stack;
	OT& ot = *this->ot; // cache 'ot' reference on host
	for (uint32 i = 0; i < inf.maxVar; i++) {
		const uint32 cand = vars->eligible[i];
		CHECKVAR(cand);
		if (cand > MAXVARTOELIM) continue;
		if (iassumed(cand)) continue;
		if (sp->vstate[cand].state) continue;
		if (sp->frozen[cand]) continue;
		const uint32 p = V2L(cand), n = NEG(p);
		assert(ot[p].size() >= cuhist[p]);
		assert(ot[n].size() >= cuhist[n]);
		if (!cuhist[p] && !cuhist[n]) continue;
		const uint32 pos_temp = opts.mu_pos << mu_inc, neg_temp = opts.mu_neg << mu_inc;
		if (cuhist[p] >= pos_temp && cuhist[n] >= neg_temp) break;
		assert(!sp->vstate[cand].state);
		vars->pVars->_push(cand);
		depFreeze(ot[p], cand, pos_temp, neg_temp);
		depFreeze(ot[n], cand, pos_temp, neg_temp);
	}
	vars->numPVs = vars->pVars->size();
	assert(verifyLCVE());
	clearFrozen();
	if (vars->numPVs) {
		const uint32 mcv = vars->pVars->back(), pmcv = V2L(mcv);
		PFLENDING(2, 5, "(%d elected, mcv: %d, pH: %d, nH: %d)", vars->numPVs, mcv, cuhist[pmcv], cuhist[NEG(pmcv)]);
		if (verbose == 4) { PFLOGN0(" PLCVs "); printVars(*vars->pVars, vars->numPVs, 'v'); }
	}
	if (vars->numPVs < opts.lcve_min) {
		if (!vars->numPVs) PFLDONE(2, 5);
		if (verbose > 1) PFLOGW("parallel variables not enough -> skip GPU simplifier");
		return false;
	}
	return true;
}

inline void ParaFROST::depFreeze(OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	LIT_ST* frozen = sp->frozen;
	uint32*& frozen_stack = sp->stacktail;
	CNF& cnf = *this->cnf; // cache 'cnf' reference on host
	forall_occurs(ol, i) {
		SCLAUSE& c = cnf[*i];
		if (c.deleted()) continue;
		forall_clause(c, k) {
			const uint32 lit = *k, v = ABS(lit);
			if (!frozen[v] && NEQUAL(v, cand)) {
				const uint32 p = POS(lit), n = NEG(lit);
				assert(NEG(p) == n);
				if (cuhist[p] < p_temp || cuhist[n] < n_temp) {
					frozen[v] = 1;
					assert(frozen_stack < sp->tmp_stack + inf.maxVar);
					*frozen_stack++ = v;
				}
			}
		}
	}
}

