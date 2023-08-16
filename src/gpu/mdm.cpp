/***********************************************************************[mdm.cpp]
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

inline bool Solver::verifyMDM() {
    for (uint32 i = sp->propagated; i < trail.size(); i++) {
        uint32 v = ABS(trail[i]);
        if (sp->frozen[v]) {
            PFLOG0("");
            PFLOGEN("decision(%d) is elected and frozen", v);
            printWatched(v);
            return false;
        }
    }
    return true;
}

inline bool Solver::verifySeen() {
    for (uint32 v = 0; v <= inf.maxVar; v++) {
        if (sp->seen[v]) {
            PFLOG0("");
            PFLOGEN("seen(%d) is not unseen", v);
            printWatched(v);
            return false;
        }
    }
    return true;
}

inline void Solver::clearMDM() {
    assert(verifyMDM());
    uint32 *start = trail + sp->propagated, *end = trail.end();
    while (start != end)
        sp->seen[ABS(*start++)] = 0;
    assert(verifySeen());
    assert((sp->stacktail - sp->tmpstack) <= (inf.maxVar - last.mdm.decisions));
    clearFrozen();
}

inline bool Solver::valid(WL& ws) {
    forall_watches(ws, i) {
        const WATCH w = *i;
        assert(w.imp);
        if (isTrue(w.imp)) continue; // clause satisfied
        if (w.binary())
            return false; // if 'w.imp' not satisfied then it's an implication of 'cand'
        else {
            // there cannot be falsified literal as watched,
            // so validating starts from 'c + 2'
            CLAUSE& c = cm[w.ref];
            assert(c.size() > 2);
            bool satisfied = false, unAssigned = false;
            uint32 *k = c + 2, *cend = c.end();
            while (k != cend && !satisfied && !unAssigned) {
                CHECKLIT(*k);
                const LIT_ST val = sp->value[*k];
                if (UNASSIGNED(val))
                    unAssigned = true;
                else if (val)
                    satisfied = true;
                k++;
            }
            if (!satisfied && !unAssigned) return false;
        }
    }
    return true;
}

inline bool Solver::depFreeze(WL& ws, const uint32& cand) {
    LIT_ST* frozen = sp->frozen;
    uint32*& frozen_stack = sp->stacktail;
    forall_watches(ws, i) {
        const WATCH w = *i;
        if (isTrue(w.imp)) continue;
        CLAUSE& c = cm[w.ref];
        uint32 othervar = ABS(c[0]) ^ ABS(c[1]) ^ cand;
        if (sp->seen[othervar]) return false;
        if (!frozen[othervar]) {
            frozen[othervar] = 1;
            assert(frozen_stack < sp->tmpstack + inf.maxVar);
            *frozen_stack++ = othervar;
        }
        uint32 *k = c + 2, *cend = c.end();
        while (k != cend) {
            othervar = ABS(*k++);
            if (sp->seen[othervar]) return false;
            if (!frozen[othervar]) {
                frozen[othervar] = 1;
                assert(frozen_stack < sp->tmpstack + inf.maxVar);
                *frozen_stack++ = othervar;
            }
        }
    }
    return true;
}

bool Solver::canMMD() {
    if (!opts.mdm_rounds) return false;
    if (uint64(opts.mdm_delay) > stats.conflicts) return false;
    const bool enough = varsEnough();
    const bool rounds = last.mdm.rounds;
    if (enough && !rounds && stats.conflicts >= limit.mdm) {
        last.mdm.rounds = opts.mdm_rounds;
        INCREASE_LIMIT(this, mdm, stats.mdm.calls, nlogn, true);
    }
    return enough && rounds;
}

void Solver::MDMInit() {
    if (!last.mdm.rounds) return;
    assert(inf.unassigned);
    assert(sp->propagated == trail.size());
    assert(conflict == NOREF);
    assert(UNSOLVED(cnfstate));
    if (opts.mdm_walk_en) {
        stats.mdm.walks++;
        walk();
    }
    // check if formula is solved by walk strategy
    if (!UNSOLVED(cnfstate)) return;
    stats.mdm.calls++;
    PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);
    eligible.resize(inf.maxVar);
    occurs.resize(inf.maxVar + 1);
    varOrder(); // initial variable ordering
    sp->stacktail = sp->tmpstack;
    bool skipround = false;
    if (assumptions.size()) {
        assert(sp->stacktail == sp->tmpstack);
        int level = DL();
        while (level < assumptions.size()) {
            const uint32 a = assumptions[level];
            CHECKLIT(a);
            assert(ifrozen[ABS(a)]);
            const uint32 cand = ABS(a);
            const LIT_ST val = sp->value[a];
            if (UNASSIGNED(val)) {
                level++;
                uint32 dec = a;
                if (!depFreeze(wt[dec], cand))
                    skipround = true;
                enqueueDecision(dec);
                sp->seen[cand] = 1;
            } else if (val)
                incDL(), level = DL();
            else {
                ianalyze(FLIP(a));
                cnfstate = UNSAT;
                clearMDM(), eligible.clear(true), occurs.clear(true);
                return;
            }
        }
    } else
        assert(sp->stacktail == sp->tmpstack);
    if (!skipround) {
        for (uint32 i = 0; i < inf.maxVar; i++) {
            uint32 cand = eligible[i];
            CHECKVAR(cand);
            if (sp->frozen[cand] || sp->vstate[cand].state || iassumed(cand)) continue;
            const LIT_ST pol = sp->psaved[cand];
            assert(pol >= 0);
            const uint32 dec = V2DEC(cand, pol);
            if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
                enqueueDecision(dec);
                sp->seen[cand] = 1;
            }
        }
    }
    last.mdm.decisions = trail.size() - sp->propagated;
    last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
    last.mdm.rounds--;
    stats.decisions.multiple += last.mdm.decisions;
    PFLOG2(2, " MDM %d: %d decisions are elected (%.2f%%)",
           stats.mdm.calls, last.mdm.decisions, percent(last.mdm.decisions, maxActive()));
    if (opts.mdm_vsids_pumps || opts.mdm_vmtf_pumps) pumpFrozen();
    clearMDM(), eligible.clear(true), occurs.clear(true);
    printStats(1, 'm', CMDM);
}

void Solver::MDM() {
    const bool vsidsActive = vsidsEnabled();
    if (opts.mdmvsidsonly_en && !vsidsActive) {
        last.mdm.rounds--;
        return;
    }
    assert(inf.unassigned);
    assert(sp->propagated == trail.size());
    assert(conflict == NOREF);
    assert(UNSOLVED(cnfstate));
    stats.mdm.calls++;
    PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);
    eligible.clear();
    if (vsidsActive)
        eligibleVSIDS();
    else
        eligibleVMFQ();
    assert(eligible.size() >= 1);
    const bool firstround = last.mdm.rounds == opts.mdm_rounds;
    if (opts.mdm_walk_en && !DL() && firstround) {
        stats.mdm.walks++;
        walk();
    }
    // check if formula is solved by walk strategy
    if (!UNSOLVED(cnfstate)) {
        eligible.clear();
        return;
    }
    sp->stacktail = sp->tmpstack;
    bool skipround = false;
    if (assumptions.size()) {
        assert(sp->stacktail == sp->tmpstack);
        int level = DL();
        while (level < assumptions.size()) {
            const uint32 a = assumptions[level];
            CHECKLIT(a);
            assert(ifrozen[ABS(a)]);
            const uint32 cand = ABS(a);
            const LIT_ST val = sp->value[a];
            if (UNASSIGNED(val)) {
                level++;
                uint32 dec = a;
                if (!depFreeze(wt[dec], cand))
                    skipround = true;
                enqueueDecision(dec);
                sp->seen[cand] = 1;
            } else if (val)
                incDL(), level = DL();
            else {
                ianalyze(FLIP(a));
                cnfstate = UNSAT;
                clearMDM();
                return;
            }
        }
    } else
        assert(sp->stacktail == sp->tmpstack);
    if (!skipround) {
        for (uint32 i = 0; i < eligible.size(); i++) {
            uint32 cand = eligible[i];
            CHECKVAR(cand);
            if (sp->frozen[cand] || sp->vstate[cand].state || iassumed(cand)) continue;
            uint32 dec = makeAssign(cand, useTarget());
            if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
                enqueueDecision(dec);
                sp->seen[cand] = 1;
            }
        }
    }
    last.mdm.decisions = trail.size() - sp->propagated;
    last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
    last.mdm.rounds--;
    stats.decisions.multiple += last.mdm.decisions;
    PFLOG2(2, " MDM %d: %d decisions are elected (%.2f%%)",
           stats.mdm.calls, last.mdm.decisions, percent(last.mdm.decisions, maxActive()));
    if (opts.mdm_vsids_pumps || opts.mdm_vmtf_pumps) pumpFrozen();
    clearMDM();
    printStats(firstround, 'm', CMDM);
}