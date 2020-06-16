/***********************************************************************[pfguess.cpp]
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

void ParaFROST::guess()
{
    // inspired by lucky phases idea in Cadical
    assert(sp->propagated == trail.size());
    assert(conflict == NOREF);
    assert(cnfstate == UNSOLVED);
    search_guess = true;
    if (sp->simplified) wtBin.recycle();
	allNegs();
	if (cnfstate == UNSOLVED) allPoss();
	if (cnfstate == UNSOLVED) forwardPoss();
	if (cnfstate == UNSOLVED) forwardNegs();
	if (cnfstate == UNSOLVED) backwardNegs();
	if (cnfstate == UNSOLVED) backwardPoss();
    search_guess = false;
    conflict = NOREF;
}

void ParaFROST::allNegs() {
    PFLOG2(2, " Guessing @%s", __func__);
    if (inf.nOrgBins) {
        for (uint32 v = 1; v <= inf.maxVar; v++) {
            WL& wBins = wtBin[neg(v2l(v))];
            for (int i = 0; i < wBins.size(); i++) {
                assert(!cm[wBins[i].ref].deleted());
                if (!sign(wBins[i].imp)) {
                    PFLOG2(2, " Found only positive binary");
                    return;
                }
            }
        }
    }
    for (uint32 i = 0; i < orgs.size(); i++) {
        CLAUSE& c = cm[orgs[i]];
        assert(!c.deleted());
        uint32* lit = c, *cend = c.end();
        LIT_ST isNeg = 0;
        while (lit != cend && !(isNeg = sign(*lit++)));
        if (isNeg) continue;
        PFLOG2(2, " Found only positive clause");
        return;
    }
    PFLOG2(2, " All clauses are negatives");
    for (uint32 v = 1; v <= inf.maxVar; v++) {
        if (!guessing(neg(v2l(v)))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    }
    stats.guess_succ = 1;
    stats.guess_who = "all negatives succeeded";
    cnfstate = SAT;
}

void ParaFROST::allPoss() {
    PFLOG2(2, " Guessing @%s", __func__);
    if (inf.nOrgBins) {
        for (uint32 v = 1; v <= inf.maxVar; v++) {
            WL& wBins = wtBin[v2l(v)];
            for (int i = 0; i < wBins.size(); i++) {
                assert(!cm[wBins[i].ref].deleted());
                if (sign(wBins[i].imp)) {
                    PFLOG2(2, " Found only negative binary");
                    return;
                }
            }
        }
    }
    for (uint32 i = 0; i < orgs.size(); i++) {
        CLAUSE& c = cm[orgs[i]];
        assert(!c.deleted());
        uint32* lit = c, * cend = c.end();
        LIT_ST isPos = 0;
        while (lit != cend && (isPos = sign(*lit++)));
        if (!isPos) continue;
        PFLOG2(2, " Found only negative clause");
        return;
    }
    PFLOG2(2, " All clauses are positives");
    for (uint32 v = 1; v <= inf.maxVar; v++) {
        if (!guessing(v2l(v))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    }
    stats.guess_succ = 1;
    stats.guess_who = "all positives succeeded";
    cnfstate = SAT;
}

void ParaFROST::forwardNegs()
{
    PFLOG2(2, " Guessing @%s", __func__);
    for (uint32 v = 1; v <= inf.maxVar; v++) {
        if (!guessing(neg(v2l(v)))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    }
    stats.guess_succ = 1;
    stats.guess_who = "forward negatives succeeded";
    cnfstate = SAT;
}

void ParaFROST::forwardPoss()
{
    PFLOG2(2, " Guessing @%s", __func__);
    for (uint32 v = 1; v <= inf.maxVar; v++) {
        if (!guessing(v2l(v))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    }
    stats.guess_succ = 1;
    stats.guess_who = "forward positives succeeded";
    cnfstate = SAT;
}

void ParaFROST::backwardNegs()
{
    PFLOG2(2, " Guessing @%s", __func__);
    uint32 v = inf.maxVar;
    do {
        if (!guessing(neg(v2l(v)))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    } while (--v > 0);
    stats.guess_succ = 1;
    stats.guess_who = "backward negatives succeeded";
    cnfstate = SAT;
}

void ParaFROST::backwardPoss()
{
    PFLOG2(2, " Guessing @%s", __func__);
    uint32 v = inf.maxVar;
    do {
        if (!guessing(v2l(v))) {
            PFLOG2(2, " Guessing @%s failed due to a conflict", __func__);
            return;
        }
    } while (--v > 0);
    stats.guess_succ = 1;
    stats.guess_who = "backward positives succeeded";
    cnfstate = SAT;
}