/***********************************************************************[pfwatch.cpp]
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

using namespace pFROST;

void ParaFROST::attachBins(BCNF& src)
{
    assert(!wt.empty());
    forall_cnf(src, i) {
        const C_REF r = *i;
        if (cm.deleted(r)) continue;
        assert(r < cm.size());
        const CLAUSE& c = cm[r];
        if (c.binary()) attachWatch(r, c);
    }
}

void ParaFROST::attachNonBins(BCNF& src)
{
    assert(!wt.empty());
    const int level = DL();
    forall_cnf(src, i) {
        const C_REF r = *i;
        if (cm.deleted(r)) continue;
        assert(r < cm.size());
        CLAUSE& c = cm[r];
        if (c.binary()) continue;
        sortClause(c);
        attachWatch(r, c);
        if (!level) {
            if (isFalse(c[0])) {
                const uint32 p0 = sp->index[ABS(c[0])];
                if (p0 < sp->propagated) sp->propagated = p0;
            }
            if (isFalse(c[1])) {
                const uint32 p1 = sp->index[ABS(c[1])];
                if (p1 < sp->propagated) sp->propagated = p1;
            }
        }
    }
}

void ParaFROST::attachClauses(BCNF& src)
{
    assert(!wt.empty());
    const int level = DL();
    forall_cnf(src, i) {
        const C_REF r = *i;
        if (cm.deleted(r)) continue;
        assert(r < cm.size());
        CLAUSE& c = cm[r];
        if (!c.binary()) sortClause(c);
        attachWatch(r, c);
        if (!level && !c.binary()) {
            if (isFalse(c[0])) {
                const uint32 p0 = sp->index[ABS(c[0])];
                if (p0 < sp->propagated) sp->propagated = p0;
            }
            if (isFalse(c[1])) {
                const uint32 p1 = sp->index[ABS(c[1])];
                if (p1 < sp->propagated) sp->propagated = p1;
            }
        }
    }
}

void ParaFROST::rebuildWT(const CL_ST& code)
{
    wt.resize(inf.nDualVars);
    uint32 proped = sp->propagated;
    if (PRIORALLBINS(code)) {
        if (PRIORLEARNTBINS(code)) {
            attachBins(learnts);
            attachBins(orgs);
        }
        else {
            attachBins(orgs);
            attachBins(learnts);
        }
        attachNonBins(orgs);
        attachNonBins(learnts);
    }
    else {
        attachClauses(orgs);
        attachClauses(learnts);
    }
    if (proped != sp->propagated)
        PFLOG2(2, "  trail reset to position %d at literal %d in %s", sp->propagated, l2i(trail[sp->propagated]), __func__);
}

void ParaFROST::sortWT()
{
    WL saved;
    forall_literal(lit) {
        assert(saved.empty());
        WL& ws = wt[lit];
        WATCH *j = ws;
        forall_watches(ws, i) {
            const WATCH w = *i;
            if (w.binary()) *j++ = w;
            else saved.push(w);
        }
        ws.resize(int(j - ws));
        forall_watches(saved, i) { ws.push(*i); }
        saved.clear();
    }
    saved.clear(true);
}