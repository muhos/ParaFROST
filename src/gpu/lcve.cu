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

#include "solve.hpp"

using namespace ParaFROST;

void Solver::varReorder() {
    PFLOGN2(2, " Finding eligible variables for LCVE..");
    assert(cuhist.d_hist != NULL);
    // NOTE: OT creation will be synced in calcScores call
    if (vars->nUnits)
        calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
    else
        calcScores(vars, cuhist.d_hist);
    cuhist.cacheHist(streams[2]);
    if (gopts.profile_gpu) cutimer->start(streams[3]);
    cacher.insert(cumm.scatter(), cumm.scatterCap());
    thrust::sort(thrust::cuda::par(tca).on(streams[3]), vars->eligible, vars->eligible + inf.maxVar, GPU_LCV_CMP(vars->scores));
    cacher.erase(cumm.scatterCap());
    PFLDONE(2, 5);
    vars->nUnits = 0;
    sync(streams[2]);
    if (gopts.profile_gpu) cutimer->stop(streams[3]), cutimer->vo += cutimer->gpuTime();
    if (verbose == 4) {
        PFLOG0(" Eligible variables:");
        for (uint32 v = 0; v < inf.maxVar; v++) {
            uint32 x = vars->eligible[v], p = V2L(x), n = NEG(p);
            PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, cuhist[p], cuhist[n], vars->scores[x]);
        }
    }
}

bool Solver::LCVE() {
    // reorder variables
    varReorder();
    // extended LCVE
    PFLOGN2(2, " Electing variables (p-mu: %d, n-mu: %d)..", opts.mu_pos << mu_inc, opts.mu_neg << mu_inc);
    vars->numPVs = 0, vars->pVars->clear();
    sp->stacktail = sp->tmpstack;
    OT& ot = *this->ot; // cache 'ot' reference on host
    for (uint32 i = 0; i < inf.maxVar; i++) {
        const uint32 cand = vars->eligible[i];
        CHECKVAR(cand);
        if (sp->frozen[cand]) continue;
        if (sp->vstate[cand].state) continue;
        if (iassumed(cand)) continue;
        const uint32 p = V2L(cand), n = NEG(p);
        assert(ot[p].size() >= cuhist[p]);
        assert(ot[n].size() >= cuhist[n]);
        if (!cuhist[p] && !cuhist[n]) continue;
        const uint32 pMu = opts.mu_pos << mu_inc, nMu = opts.mu_neg << mu_inc;
        if (cuhist[p] > opts.lcve_max || cuhist[n] > opts.lcve_max) break;
        if (cuhist[p] >= pMu && cuhist[n] >= nMu) break;
        assert(!sp->vstate[cand].state);
        if (depFreeze(cand, ot[p]) && depFreeze(cand, ot[n]))
            vars->pVars->_push(cand);
    }
    vars->numPVs = vars->pVars->size();
    assert(verifyLCVE());
    if (vars->numPVs) {
        const uint32 mcv = vars->pVars->back(), pmcv = V2L(mcv);
        PFLENDING(2, 5, "(%d elected, mcv: %d, pH: %d, nH: %d)", vars->numPVs, mcv, cuhist[pmcv], cuhist[NEG(pmcv)]);
        if (verbose == 4) {
            PFLOGN0(" PLCVs ");
            printVars(*vars->pVars, vars->numPVs, 'v');
        }
    }
    mapFrozen();
    clearFrozen();
    if (vars->numPVs < opts.lcve_min) {
        if (opts.ve_fun_en) sync();
        if (!vars->numPVs) PFLDONE(2, 5);
        if (verbose > 1) PFLOGW("parallel variables not enough -> skip GPU simplifier");
        return false;
    }
    return true;
}

inline void Solver::mapFrozen() {
    if (!opts.ve_fun_en) return;
    const uint32* frozen = sp->tmpstack;
    const uint32* end = sp->stacktail;
    const uint32 nFrozen = uint32(end - frozen);
    if (!nFrozen) {
        vars->varcore = NULL;
        return;
    }
    assert(nFrozen <= inf.maxVar);
    CHECK(cudaMemcpy(vars->scores, frozen, nFrozen * sizeof(uint32), cudaMemcpyHostToDevice));
    mapFrozenAsync(vars, nFrozen);
}

inline bool Solver::depFreeze(const uint32& cand, OL& ol) {
    LIT_ST* frozen = sp->frozen;
    uint32*& frozentail = sp->stacktail;
    uint32* frozenstart = frozentail;
    CNF& cnf = *this->cnf; // cache 'cnf' reference on host
    forall_occurs(ol, i) {
        SCLAUSE& c = cnf[*i];
        if (c.deleted()) continue;
        if (c.size() > opts.lcve_clause_limit) {
            uint32* from = frozenstart;
            while (from != frozentail)
                sp->frozen[*from++] = 0;
            frozentail = frozenstart;
            return false;
        }
        forall_clause(c, k) {
            const uint32 lit = *k, v = ABS(lit);
            if (!frozen[v] && NEQUAL(v, cand)) {
                frozen[v] = 1;
                assert(frozentail < sp->tmpstack + inf.maxVar);
                *frozentail++ = v;
            }
        }
    }
    return true;
}
