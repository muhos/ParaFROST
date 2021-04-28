/***********************************************************************[pfaddinc.cpp]
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

uint32 ParaFROST::incVariable() {
	inf.unassigned++;
	const uint32 v = ++inf.maxVar;
	PFLOG2(4, "  adding new variable %d (%d unassigned)..", v, inf.unassigned);
	const uint32 lit = V2L(v);
	inf.nDualVars = lit + 2;
	wt.expand(lit + 2);
	ivalue.expand(lit + 2, UNDEFINED);
	bumps.expand(v + 1, 0);
	activity.expand(v + 1, 0.0);
	imarks.expand(v + 1, UNDEFINED);
	ilevel.expand(v + 1, UNDEFINED);
	ifrozen.expand(v + 1, 0);
	ivstate.expand(v + 1), ivstate[v] = VSTATE();
	vorg.expand(v + 1), vorg[v] = v;
	vmtf.init(v);
	vmtf.update(v, (bumps[v] = ++bumped));
	vsids.insert(v);
	if (sp == NULL) {
		sp = new SP();
		if (opts.proof_en)
			proof.init(sp);
	}
	sp->value = ivalue;
	sp->level = ilevel;
	sp->marks = imarks;
	sp->vstate = ivstate;
	return v;
}

bool ParaFROST::itoClause(Lits_t& c, Lits_t& org)
{
	if (org.empty()) {
		learnEmpty();
		return false;
	}
	assert(c.empty());
	bool satisfied = false;
	forall_clause(org, k) {
		const uint32 lit = *k;
		CHECKLIT(lit);
		LIT_ST marker = l2marker(lit);
		if (UNASSIGNED(marker)) {
			markLit(lit);
			LIT_ST val = sp->value[lit];
			if (UNASSIGNED(val)) c.push(lit);
			else if (val) satisfied = true;
		}
		else if (marker != SIGN(lit)) satisfied = true; // tautology
	}
	forall_clause(org, k) {
		unmarkLit(*k);
	}
	if (satisfied) {
		if (opts.proof_en) proof.deleteClause(org);
	}
	else {
		if (c.empty()) {
			learnEmpty();
			return false;
		}
		int newsize = c.size();
		if (newsize == 1) {
			const uint32 unit = *c;
			CHECKLIT(unit);
			LIT_ST val = sp->value[unit];
			if (UNASSIGNED(val)) enqueueUnit(unit), formula.units++;
			else if (!val) return false;
		}
		else if (newsize) {
			if (newsize == 2) formula.binaries++;
			else if (newsize == 3) formula.ternaries++;
			else assert(newsize > 3), formula.large++;
			if (newsize > formula.maxClauseSize)
				formula.maxClauseSize = newsize;
			const C_REF newref = newClause(c, false);
			PFLCLAUSE(4, cm[newref], "  adding new clause");
		}
		if (opts.proof_en && newsize < org.size()) {
			proof.addClause(c);
			proof.deleteClause(org);
			org.clear();
		}
	}
	c.clear(), org.clear();
	return true;
}