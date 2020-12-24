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
    C_REF* end = src.end();
    for (C_REF* i = src; i != end; i++) {
        C_REF r = *i;
        assert(r < cm.size());
        CLAUSE& c = cm[r];
        if (c.deleted() || !c.binary()) continue;
        attachWatch(r, c);
    }
}

void ParaFROST::attachNonBins(BCNF& src)
{
    assert(!wt.empty());
    C_REF* end = src.end();
    for (C_REF* i = src; i != end; i++) {
        C_REF r = *i;
        assert(r < cm.size());
        CLAUSE& c = cm[r];
        if (c.deleted() || c.binary()) continue;
        attachWatch(r, c);
    }
}

void ParaFROST::attachClauses(BCNF& src)
{
    assert(!wt.empty());
    C_REF* end = src.end();
    for (C_REF* i = src; i != end; i++) {
        C_REF r = *i;
        assert(r < cm.size());
        CLAUSE& c = cm[r];
        if (c.deleted()) continue;
        attachWatch(r, c);
    }
}

void ParaFROST::rebuildWT(const bool& binfirst)
{
	assert(!wt.empty());
    if (binfirst) {
        attachBins(orgs), attachBins(learnts);
        attachNonBins(orgs), attachNonBins(learnts);
    }
    else attachClauses(orgs), attachClauses(learnts);
}