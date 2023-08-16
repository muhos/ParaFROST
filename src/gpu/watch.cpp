/***********************************************************************[watch.cpp]
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

void Solver::attachBins(BCNF& src, const bool& hasElim) {
    assert(!wt.empty());
    if (hasElim) {
        const VSTATE* states = sp->vstate;
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            CLAUSE& c = cm[r];
            if (c.binary()) {
                if (MELTED(states[ABS(c[0])].state)
                    || MELTED(states[ABS(c[1])].state)) {
                    removeClause(c, r);
                } else
                    attachWatch(r, c);
            }
        }
    } else {
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            const CLAUSE& c = cm[r];
            if (c.binary()) attachWatch(r, c);
        }
    }
}

void Solver::attachNonBins(BCNF& src, const bool& hasElim) {
    assert(!wt.empty());
    if (hasElim) {
        const VSTATE* states = sp->vstate;
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            CLAUSE& c = cm[r];
            if (c.binary()) continue;
            bool removed = false;
            forall_clause(c, k) {
                const uint32 lit = *k;
                if (MELTED(states[ABS(lit)].state)) {
                    removed = true;
                    break;
                }
            }
            if (removed)
                removeClause(c, r);
            else {
                sortClause(c);
                attachWatch(r, c);
            }
        }
    } else {
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            CLAUSE& c = cm[r];
            if (c.binary()) continue;
            sortClause(c);
            attachWatch(r, c);
        }
    }
}

void Solver::attachClauses(BCNF& src, const bool& hasElim) {
    assert(!wt.empty());
    if (hasElim) {
        const VSTATE* states = sp->vstate;
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            CLAUSE& c = cm[r];
            bool removed = false;
            forall_clause(c, k) {
                const uint32 lit = *k;
                if (MELTED(states[ABS(lit)].state)) {
                    removed = true;
                    break;
                }
            }
            if (removed)
                removeClause(c, r);
            else {
                if (!c.binary()) sortClause(c);
                attachWatch(r, c);
            }
        }
    } else {
        forall_cnf(src, i) {
            const C_REF r = *i;
            if (cm.deleted(r)) continue;
            assert(r < cm.size());
            CLAUSE& c = cm[r];
            if (!c.binary()) sortClause(c);
            attachWatch(r, c);
        }
    }
}

void Solver::rebuildWT(const CL_ST& code) {
    wt.resize(inf.nDualVars);
    if (PRIORALLBINS(code)) {
        if (PRIORLEARNTBINS(code)) {
            attachBins(learnts);
            attachBins(orgs);
        } else {
            attachBins(orgs);
            attachBins(learnts);
        }
        attachNonBins(orgs);
        attachNonBins(learnts);
    } else {
        attachClauses(orgs);
        attachClauses(learnts);
    }
}

void Solver::sortWT() {
    WL saved;
    forall_literal(lit) {
        assert(saved.empty());
        WL& ws = wt[lit];
        WATCH* j = ws;
        forall_watches(ws, i) {
            const WATCH w = *i;
            if (w.binary())
                *j++ = w;
            else
                saved.push(w);
        }
        ws.resize(int(j - ws));
        forall_watches(saved, i) { ws.push(*i); }
        saved.clear();
    }
    saved.clear(true);
}

void Solver::detachClauses(const bool& keepbinaries) {
    forall_literal(lit) {
        WL& ws = wt[lit];
        WATCH* j = ws;
        forall_watches(ws, i) {
            const WATCH w = *i;
            if (keepbinaries && w.binary())
                *j++ = w;
        }
        ws.resize(int(j - ws));
    }
}

void Solver::binarizeWT(const bool& keeplearnts) {
    assert(!DL());
    const LIT_ST* values = sp->value;
    forall_literal(lit) {
        const LIT_ST litval = values[lit];
        WL& ws = wt[FLIP(lit)];
        WATCH* j = ws;
        forall_watches(ws, i) {
            const WATCH w = *i;
            if (w.binary()) {
                const C_REF ref = w.ref;
                if (cm.deleted(ref)) continue;
                const uint32 imp = w.imp;
                if (litval > 0 || values[imp] > 0) {
                    if (lit < imp)
                        removeClause(cm[ref], ref);
                } else {
                    if (cm[ref].original() || keeplearnts)
                        *j++ = w;
                }
            }
        }
        ws.resize(int(j - ws));
    }
}
