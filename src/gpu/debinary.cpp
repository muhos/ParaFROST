/***********************************************************************[debins.cpp]
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

void Solver::debinary() {
    if (!opts.debinary_en) return;
    if (!cnfstate) return;
    assert(!DL());
    assert(!wt.empty());
    assert(sp->propagated == trail.size());
    stats.debinary.calls++;
    uVec1D& marked = minimized;
    int64 subsumed = 0, units = 0;
    forall_literal(lit) {
        assert(marked.empty());
        if (inactive(lit)) continue;
        uint32 unit = 0;
        WL& ws = wt[lit];
        WATCH* j = ws;
        forall_watches(ws, i) {
            const WATCH w = *j++ = *i;
            const C_REF cref = w.ref;
            if (cm.deleted(cref)) {
                j--;
                continue;
            }
            if (w.binary()) {
                const uint32 other = w.imp;
                CHECKLIT(other);
                const LIT_ST marker = l2marker(other);
                CLAUSE& c = cm[cref];
                assert(c.size() == 2);
                if (UNASSIGNED(marker)) {
                    markLit(other);
                    marked.push(other);
                } else if (NEQUAL(marker, SIGN(other))) { // found 'hyper unary'
                    unit = FLIP(lit);
                    j = ws; // the whole list is satisfied by 'unit'
                    units++;
                    break;
                } else { // found duplicate
                    PFLCLAUSE(4, c, "  found duplicated binary");
                    if (c.original()) { // find learnt duplicate if exists
                        for (WATCH* k = ws;; k++) {
                            assert(k != i);
                            if (!k->binary() || NEQUAL(k->imp, other)) continue;
                            const C_REF dref = k->ref;
                            if (cm.deleted(dref)) continue;
                            assert(cm[dref].size() == 2);
                            assert(!cm[dref].deleted());
                            removeClause(cm[dref], dref);
                            *k = w;
                            break;
                        }
                    } else
                        removeClause(c, cref);
                    subsumed++;
                    j--;
                }
            }
        }
        if (j != ws)
            ws.resize(int(j - ws));
        else
            ws.clear(true);
        forall_vector(uint32, marked, i) { unmarkLit(*i); }
        marked.clear();
        if (unit) {
            CHECKLIT(unit);
            enqueueUnit(unit);
            if (BCP()) {
                learnEmpty();
                break;
            }
        }
    }
    stats.debinary.hyperunary += units;
    stats.debinary.binaries += subsumed;
    PFLOG2(2, " Deduplicate %lld: removed %lld binaries, producing %lld hyper unaries", stats.debinary.calls, subsumed, units);
    printStats(units || subsumed, 'd', CVIOLET2);
}