/***********************************************************************[pfassume.cpp]
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

bool ParaFROST::ifailed(const uint32& v)
{
	const uint32 mlit = imap(v);
	if (!mlit) return false;
	CHECKLIT(mlit);
	const int size = iconflict.size();
	assert(size);
	for (int i = 0; i < size; i++) {
		if (ABS(mlit) == v)
			return true;
	}
	return false;
}

bool ParaFROST::ieliminated(const uint32& v)
{
	const uint32 mlit = imap(v);
	if (!mlit) return true;
	CHECKLIT(mlit);
	return MELTED(sp->vstate[ABS(mlit)].state) || SUBSTITUTED(sp->vstate[ABS(mlit)].state);
}

void ParaFROST::ifreeze(const uint32& v) 
{
	const uint32 mlit = imap(v);
	CHECKLIT(mlit);
	const uint32 mvar = ABS(mlit);
	ifrozen[mvar] = 1;
	PFLOG2(3, "  freezing original variable %d (mapped to %d)..", v, mvar);
}

void ParaFROST::iunfreeze(const uint32& v)
{
	const uint32 mlit = imap(v);
	CHECKLIT(mlit);
	const uint32 mvar = ABS(mlit);
	ifrozen[mvar] = 0;
	PFLOG2(3, "  melting original variable %d (mapped to %d)..", v, mvar);
}

void ParaFROST::iassume(Lits_t& assumptions)
{
	assert(inf.maxVar);
	assert(stats.clauses.original);
	if (assumptions.empty()) return;
	PFLOGN2(2, " Adding %d assumptions..", assumptions.size());
	this->assumptions.reserve(assumptions.size());
	forall_clause(assumptions, k) {
		const uint32 a = *k, v = ABS(a);
		CHECKLIT(a);
		assert(!ieliminated(v));
		assert(!ifrozen[v]);
		ifrozen[ABS(a)] = 1;
		this->assumptions.push(a);
	}
	PFLDONE(2, 5);
}

void ParaFROST::iunassume()
{
	assert(inf.maxVar);
	assert(stats.clauses.original);
	PFLOGN2(2, " Resetting %d assumptions and solver state..", assumptions.size());
	if (assumptions.size()) {
		forall_clause(assumptions, k) {
			const uint32 a = *k;
			CHECKLIT(a);
			ifrozen[ABS(a)] = 0;
		}
		assumptions.clear(true);
	}
	cnfstate = UNSOLVED;
	PFLDONE(2, 5);
	backtrack();
}

void ParaFROST::idecide()
{
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	int level = DL();
	uint32 dec = 0;
	while (level < assumptions.size()) {
		const uint32 a = assumptions[level];
		CHECKLIT(a);
		assert(ifrozen[ABS(a)]);
		const LIT_ST val = sp->value[a];
		if (UNASSIGNED(val)) {
			dec = a;
			break;
		}
		else if (val) incDL(), level = DL();
		else {
			assert(!val);
			ianalyze(FLIP(a));
			cnfstate = UNSAT;
			return;
		}
	}
	if (!dec) {
		const uint32 cand = vsidsEnabled() ? nextVSIDS() : nextVMFQ();
		dec = makeAssign(cand, useTarget());
	}
	enqueueDecision(dec);
	stats.decisions.single++;
}

void ParaFROST::ianalyze(const uint32& failed)
{
	assert(cnfstate);
	assert(vorg);
	PFLOG2(3, " Analyzing conflict on failed assumption (%d):", l2i(failed));
	PFLTRAIL(this, 3);
	iconflict.clear();
	iconflict.push(failed);
	if (!DL()) return;
	sp->seen[ABS(failed)] = ANALYZED_M;
	assert(trail.size());
	uint32* tail = trail.end() - 1, *pivot = trail + dlevels[1];
	while (tail >= pivot) {
		const uint32 parent = *tail--;
		CHECKLIT(parent);
		const uint32 parentv = ABS(parent);
		if (sp->seen[parentv]) {
			const C_REF r = sp->source[parentv];
			if (REASON(r)) {
				CLAUSE& c = cm[r];
				PFLCLAUSE(4, c, "  analyzing %d reason", l2i(parent));
				forall_clause(c, k) {
					const uint32 other = *k;
					if (other == parent) continue;
					const uint32 v = ABS(other);
					CHECKVAR(v);
					if (sp->level[v]) 
						sp->seen[v] = ANALYZED_M;
				}
			}
			else {
				assert(sp->level[parentv]);
				iconflict.push(FLIP(parent));
			}
			sp->seen[parentv] = 0;
		}
	}
	sp->seen[ABS(failed)] = 0;
}