/***********************************************************************[vmap.cpp]
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

#include "control.hpp"
#include "solve.hpp"

using namespace ParaFROST;

void Solver::map(BCNF& cnf) {
    if (cnf.empty()) return;
    assert(!vmap.empty());
    forall_cnf(cnf, i) {
        vmap.mapClause(cm[*i]);
    }
}

void Solver::map(WL& ws) {
    if (ws.empty()) return;
    forall_watches(ws, w) {
        w->imp = vmap.mapLit(w->imp);
    }
}

void Solver::map(WT& wt) {
    if (wt.empty()) return;
    assert(!vmap.empty());
    forall_variables(v) {
        const uint32 mVar = vmap.mapped(v);
        if (mVar) {
            const uint32 p = V2L(v), n = NEG(p);
            const uint32 mpos = V2L(mVar), mneg = NEG(mpos);
            if (NEQUAL(mVar, v)) { // map watch lists
                wt[mpos].copyFrom(wt[p]);
                wt[mneg].copyFrom(wt[n]);
            }
            map(wt[mpos]), map(wt[mneg]); // then map watch imps
        }
    }
    wt.resize(V2L(vmap.size()));
    wt.shrinkCap();
}

void Solver::map(const bool& sigmified) {
    assert(inf.unassigned);
    assert(conflict == NOREF);
    assert(!DL());
    assert(trail.size() == sp->propagated);
    stats.mappings++;
    int64 memBefore = sysMemUsed();
    vmap.initiate(sp);
    // map model literals
    vmap.mapOrgs(model.lits);
    // map assumptions & frozen
    assert(iconflict.empty());
    if (assumptions.size()) {
        // unfreeze unmapped
        forall_clause(assumptions, k) {
            const uint32 a = *k;
            CHECKLIT(a);
            assert(ifrozen[ABS(a)]);
            ifrozen[ABS(a)] = 0;
        }
        // map assumptions
        vmap.mapOrgs(assumptions);
        // freeze mapped
        forall_clause(assumptions, k) {
            const uint32 a = *k;
            CHECKLIT(a);
            assert(!ifrozen[ABS(a)]);
            ifrozen[ABS(a)] = 1;
        }
    }
    // map transitive start literal
    vmap.mapTransitive(last.transitive.literals);
    // map clauses and watch tables
    if (!sigmified)
        map(orgs), map(learnts), map(wt);
    else
        mapped = true, newBeginning(), mapped = false;
    // map trail, queue and heap
    vmap.mapShrinkLits(trail);
    sp->propagated = trail.size();
    const uint32 firstFrozen = vmap.firstL0();
    if (!firstFrozen) assert(trail.empty());
    vmtf.map(*vmap, firstFrozen);
    vmap.mapShrinkVars(vmtf.data());
    vmap.mapShrinkVars(bumps);
    uVec1D tmp;
    while (vsids.size()) {
        uint32 x = vsids.pop();
        if (x == firstFrozen) continue;
        uint32 mx = vmap.mapped(x);
        if (mx) tmp.push(mx);
    }
    vmap.mapShrinkVars(activity);
    vsids.rebuild(tmp);
    // map search space
    SP* newSP = new SP(vmap.size());
    vmap.mapSP(newSP);
    delete sp;
    sp = newSP;
    vmap.mapShrinkVars(vorg);
    model.init(vorg);
    if (opts.proof_en) proof.init(sp, vorg);
    // update phase-saving counters
    sp->trailpivot = 0, last.rephase.best = last.rephase.target = 0;
    for (uint32 v = 1; v <= vmap.numVars(); v++) {
        if (sp->vstate[v].state) continue;
        if (!UNASSIGNED(sp->pbest[v])) last.rephase.best++;
        if (!UNASSIGNED(sp->ptarget[v])) last.rephase.target++;
    }
    PFLOG2(2, " Variable mapping compressed %d to %d, saving %.2f KB of memory",
           inf.maxVar, vmap.numVars(), double(abs(memBefore - sysMemUsed())) / KBYTE);
    inf.maxVar = vmap.numVars();
    inf.nDualVars = V2L(inf.maxVar + 1);
    inf.maxFrozen = inf.maxMelted = inf.maxSubstituted = 0;
    vmap.destroy();
    printStats(1, 'a', CGREEN2);
}