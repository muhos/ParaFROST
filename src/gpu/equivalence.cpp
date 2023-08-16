/***********************************************************************[els.cpp]
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

// 'substitute' corresponding clauses having literals represented by 'smallest'
// but start with learnts first which gives priority to substituted learnts in
// the watch table when new clauses are added

inline uint32 Solver::minReachable(WL& ws, DFS* dfs, const DFS& node) {
    uint32 new_min = node.min;
    forall_watches(ws, i) {
        const WATCH w = *i;
        if (!w.binary()) continue;
        const uint32 child = w.imp;
        CHECKLIT(child);
        if (inactive(child)) continue;
        const DFS& child_dfs = dfs[child];
        if (new_min > child_dfs.min) new_min = child_dfs.min;
    }
    return new_min;
}

bool Solver::decompose() {
    if (!cnfstate) return false;
    assert(!DL());
    assert(sp->propagated == trail.size());
    assert(analyzed.empty());
    assert(minimized.empty());
    stats.decompose.calls++;
    const uint32 dfs_size = inf.nDualVars;
    DFS* dfs = pfcalloc<DFS>(dfs_size);
    uint32 *smallests = pfcalloc<uint32>(dfs_size), dfs_idx = 0;
    uVec1D& scc = analyzed;
    uVec1D& litstack = minimized;
    uint32 substituted = 0, sccs = 0, before = maxActive();
    forall_literal(root) {
        if (!cnfstate) break;
        const uint32 root_v = ABS(root);
        if (sp->vstate[root_v].state) continue;
        if (dfs[root].min == NOVAR) continue;
        assert(scc.empty());
        assert(litstack.empty());
        litstack.push(root);
        while (cnfstate && litstack.size()) {
            uint32 parent = litstack.back();
            CHECKLIT(parent);
            DFS& parent_dfs = dfs[parent];
            if (parent_dfs.min == NOVAR) {
                CHECKLIT(smallests[parent]);
                litstack.pop();
            } else { // traverse all binaries
                assert(!smallests[parent]);
                WL& ws = wt[parent];  // 'ws' holds the negations of 'parent'
                if (parent_dfs.idx) { // all children of parent visited and min reachable found
                    litstack.pop();
                    uint32 new_min = minReachable(ws, dfs, parent_dfs); // find min. reachable from the children of 'parent'
                    PFLOG2(4, " dfs search of parent(%d) with index %d reached minimum %d", l2i(parent), parent_dfs.idx, new_min);
                    if (parent_dfs.idx == new_min) { // start of SCC block
                        // find the smallest variable to represent this SCC
                        uint32 other, size = 0, smallest = parent;
                        assert(scc.size());
                        CHECKLIT(smallest);
                        const uint32 farent = FLIP(parent);
                        uint32 j = scc.size();
                        do {
                            assert(j > 0);
                            other = scc[--j];
                            CHECKLIT(other);
                            if (NEQUAL(other, farent)) {
                                if (ABS(other) < ABS(smallest))
                                    smallest = other;
                                size++;
                            } else {
                                PFLOG2(2, " Conflict as both %d and its negation in the same SCC", l2i(parent));
                                enqueueUnit(parent);
                                learnEmpty();
                            }
                        } while (cnfstate && NEQUAL(other, parent));
                        if (cnfstate) {
                            PFLOG2(4, " New SCC of size %d and smallest variable %d", size, l2i(smallest));
                            do {
                                assert(scc.size());
                                other = scc.back();
                                scc.pop();
                                CHECKLIT(other);
                                dfs[other].min = NOVAR;
                                if (iassumed(ABS(other)))
                                    smallests[other] = other;
                                else {
                                    smallests[other] = smallest;
                                    if (NEQUAL(other, smallest)) {
                                        substituted++;
                                        PFLOG2(4, "literal %d in SCC substituted by the smallest %d", l2i(other), l2i(smallest));
                                    }
                                }
                            } while (NEQUAL(other, parent));
                            sccs += size > 1;
                        }
                    } else
                        parent_dfs.min = new_min;
                } else { // traverse children of 'parent'
                    dfs_idx++;
                    assert(dfs_idx < NOVAR);
                    parent_dfs.idx = parent_dfs.min = dfs_idx;
                    scc.push(parent);
                    PFLOG2(4, " traversing all implications of parent(%d) at index %u", l2i(parent), dfs_idx);
                    forall_watches(ws, i) {
                        const WATCH w = *i;
                        if (!w.binary()) continue;
                        const uint32 child = w.imp;
                        CHECKLIT(child);
                        if (inactive(child)) continue;
                        const DFS& child_dfs = dfs[child];
                        if (child_dfs.idx) continue;
                        litstack.push(child);
                    }
                }
            }
        }
    }
    PFLOG2(2, " Decomposition %lld: %d SCCs found, %d variables substituted (%.2f%%)",
           stats.decompose.calls, sccs, substituted, percent(substituted, before));
    stats.decompose.scc += sccs;
    stats.decompose.variables += substituted;
    free(dfs), dfs = NULL;
    scc.clear();
    litstack.clear();
    bool orgsucc = false, learntsucc = false;
    if (substituted) {
        assert(reduced.empty());
        if (cnfstate) learntsucc = substitute(learnts, smallests);
        if (cnfstate) orgsucc = substitute(orgs, smallests);
        if (cnfstate && reduced.size()) {
            forall_cnf(reduced, i) {
                const C_REF r = *i;
                assert(!cm.deleted(r));
                removeClause(cm[r], r);
            }
            reduced.clear(true);
        }
    }
    if (cnfstate) {
        recycleWT(); // must be recycled before BCP
        if (sp->propagated < trail.size() && BCP()) {
            PFLOG2(2, " Propagation after substitution proved a contradiction");
            learnEmpty();
        }
        if (cnfstate) {
            assert(UNSOLVED(cnfstate));
            forall_variables(v) {
                if (sp->vstate[v].state) continue;
                const uint32 p = V2L(v);
                const uint32 other = smallests[p];
                CHECKLIT(other);
                if (NEQUAL(other, p)) {
                    assert(other < p);
                    const uint32 othervar = ABS(other);
                    assert(!MELTED(sp->vstate[othervar].state));
                    assert(!SUBSTITUTED(sp->vstate[othervar].state));
                    if (!sp->vstate[othervar].state) {
                        PFLOG2(4, " %d substituted to %d", v, othervar);
                        markSubstituted(v);
                    }
                    model.saveBinary(p, FLIP(other));
                }
            }
        }
    }
    free(smallests), smallests = NULL;
    return !cnfstate || (substituted && (orgsucc || learntsucc));
}

bool Solver::substitute(BCNF& cnf, uint32* smallests) {
    assert(UNSOLVED(cnfstate));
    assert(learntC.empty());
    bool binaries = false;
    uint32 units = 0;
    uint32 deleted = 0, replaced = 0;
    LIT_ST* values = sp->value;
    const uint32 cnfsize = cnf.size();
    for (uint32 i = 0; cnfstate && i < cnfsize; i++) {
        const C_REF ref = cnf[i];
        if (cm.deleted(ref)) continue;
        CLAUSE& c = cm[ref];
        int j, size = c.size();
        for (j = 0; j < size; j++) {
            const uint32 lit = c[j];
            if (NEQUAL(smallests[lit], lit)) break;
        }
        if (j == size) continue;
        replaced++;
        assert(learntC.empty());
        bool satisfied = false;
        for (int k = 0; !satisfied && k < size; k++) {
            const uint32 lit = c[k];
            CHECKLIT(lit);
            LIT_ST val = values[lit];
            if (UNASSIGNED(val)) {
                const uint32 other = smallests[lit];
                CHECKLIT(other);
                val = values[other];
                if (UNASSIGNED(val)) {
                    val = l2marker(other);
                    if (UNASSIGNED(val)) {
                        markLit(other);
                        learntC.push(other);
                    } else if (NEQUAL(val, SIGN(other)))
                        satisfied = true;
                } else if (val)
                    satisfied = true;
            } else if (val)
                satisfied = true;
        }
        if (satisfied) {
            PFLCLAUSE(4, c, "  satisfied after substitution");
            reduced.push(ref);
            deleted++;
        } else if (learntC.empty()) {
            PFLOG2(2, "  learnt empty clause during decomposition");
            learnEmpty();
        } else if (learntC.size() == 1) {
            PFLOG2(4, "  found unit %d after substitution", l2i(learntC[0]));
            enqueueUnit(learntC[0]);
            removeClause(c, ref);
            units++;
        } else if (NEQUAL(c[0], learntC[0]) || NEQUAL(c[1], learntC[1])) { // watches changed, 'learntC' will be added and watched
            if (opts.proof_en) proof.addClause(learntC);
            if (learntC.size() == 2) binaries = true;
            deleted++;
            uint32 last = cnf.size();
            removeClause(c, ref);
            sp->learntLBD = c.lbd();
            C_REF newref = newClause(learntC, c.learnt());
            PFLCLAUSE(4, cm[newref], "  learnt after substitution");
            assert(cnf[last] == newref);
            cnf[last] = ref;
            cnf[i] = newref;
        } else {
            if (opts.proof_en) {
                proof.addClause(learntC);
                proof.deleteClause(c);
            }
            const int csize = c.size();
            assert(csize > 2);
            int k;
            for (k = 2; k < learntC.size(); k++) c[k] = learntC[k];
            int removed = csize - k;
            if (removed) {
                PFLOG2(4, "  only shrinking clause as watches did not change");
                if (k == 2) binaries = true;
                shrinkClause(c, removed);
                if (c.original()) stats.shrunken += removed;
            } else if (keeping(c))
                markSubsume(c);
            PFLCLAUSE(4, c, "  substituted");
        }
        unmarkLearnt();
        learntC.clear();
    }
    deleted += units;
    stats.decompose.hyperunary += units;
    stats.decompose.clauses += deleted;
    PFLOG2(2, " Decomposition %lld: %d of clauses replaced %.2f%%, producing %d deleted clauses %.2f%%",
           stats.decompose.calls, replaced, percent(replaced, cnfsize), deleted, percent(deleted, replaced));
    return (units || binaries);
}

bool Solver::canELS(const bool& first) {
    if (!opts.decompose_en) return false;
    const uint64 clauses = maxClauses();
    if (first && clauses > uint64(opts.decompose_limit)) return false;
    return (3 * clauses) < (stats.searchticks + opts.decompose_min_eff);
}

void Solver::ELS(const bool& first) {
    if (!canELS(first)) return;
    int rounds = opts.decompose_min;
    if (maxClauses() > opts.decompose_limit) rounds = 1;
    bool success = true;
    for (int round = 1; success && round <= rounds; round++)
        success = decompose();
    printStats(success, 'e', CVIOLET0);
}