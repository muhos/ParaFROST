/***********************************************************************[pfclause.cpp]
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
#include "pfdimacs.h"
using namespace pFROST;

void ParaFROST::markLits(CLAUSE& c) {
	assert(c.size() > 1);
	forall_clause(c, k) {
		assert(UNASSIGNED(l2marker(*k)));
		markLit(*k);
	}
}

void ParaFROST::unmarkLits(CLAUSE& c) {
	forall_clause(c, k) {
		assert(!UNASSIGNED(l2marker(*k)));
		unmarkLit(*k);
	}
}

void ParaFROST::markSubsume(CLAUSE& c) {
	assert(keeping(c));
	forall_clause(c, k) {
		markSubsume(*k);
	}
}

bool ParaFROST::keeping(CLAUSE& c) {
	if (c.original()) return true;
	if (c.keep()) return true;
	if (c.lbd() > limit.keptlbd) return false;
	if (c.size() > limit.keptsize) return false;
	return true;
}

void ParaFROST::removeClause(CLAUSE& c, const C_REF& cref) {
	assert(cm[cref] == c);
	assert(!c.deleted());
	assert(c.size() > 1);
	const int size = c.size();
	if (c.learnt()) {
		assert(stats.clauses.learnt > 0);
		stats.clauses.learnt--;
		assert(stats.literals.learnt > 0);
		stats.literals.learnt -= size;
	}
	else {
		assert(c.original());
		assert(stats.clauses.original > 0);
		stats.clauses.original--;
		assert(stats.literals.original > 0);
		stats.literals.original -= size;
		stats.shrunken += size;
	}
	if (opts.proof_en) proof.deleteClause(c);
	c.markDeleted();
	cm.collectClause(cref, size);
}

void ParaFROST::newClause(const C_REF& cref, CLAUSE& c, const bool& learnt)
{
	assert(cm[cref] == c);
	const int size = c.size();
	assert(size > 1);
	if (learnt) {
		assert(sp->learntLBD > 0);
		int trimlbd = sp->learntLBD > size ? size : sp->learntLBD;
		c.markLearnt();
		c.set_lbd(trimlbd);
		c.set_usage(1 + (sp->learntLBD <= opts.lbd_tier2));
		if (size > 2 && trimlbd > opts.lbd_tier1) c.set_keep(0);
		learnts.push(cref);
		stats.clauses.learnt++;
		stats.literals.learnt += size;
	}
	else {
		assert(c.original());
		orgs.push(cref);
		stats.clauses.original++;
		stats.literals.original += size;
	}
	if (keeping(c)) markSubsume(c);
}

C_REF ParaFROST::newClause(const Lits_t& in_c, const bool& learnt)
{
	const C_REF r = cm.alloc(in_c);
	CLAUSE& c = cm[r];
	assert(c.keep());
	assert(!c.deleted());
	attachWatch(r, c);
	newClause(r, c, learnt);
	return r;
}

void ParaFROST::newHyper2(const bool& learnt)
{
	assert(learntC.size() == 2);
	stats.binary.resolvents++;
	const C_REF r = cm.alloc(learntC);
	CLAUSE& c = cm[r];
	const uint32 first = c[0], second = c[1];
	delayWatch(first, second, r, 2);
	delayWatch(second, first, r, 2);
	sp->learntLBD = 2;
	newClause(r, c, learnt);
	if (learnt) c.markHyper();
	learntC.clear();
}

void ParaFROST::newHyper3(const bool& learnt)
{
	const int size = learntC.size();
	assert(size > 1 && size <= 3);
	last.ternary.resolvents++;
	const C_REF r = cm.alloc(learntC);
	CLAUSE& c = cm[r];
	sp->learntLBD = size;
	newClause(r, c, learnt);
	if (learnt) c.markHyper();
	attachClause(r, c);
	PFLCLAUSE(4, c, "  added new hyper ternary resolvent");
}

inline LIT_ST ParaFROST::sortClause(CLAUSE& c, const int& start, const int& size, const bool& satonly)
{
	assert(size > 1);
	assert(start < size);
	LIT_ST* values = sp->value;
	uint32 x = c[start];
	CHECKLIT(x);
	LIT_ST xval = values[x];
	if (UNASSIGNED(xval) || (xval && satonly)) return xval;
	uint32 best = x;
	int xlevel = l2dl(x), pos = 0;
	assert(xlevel > -1);
	for (int i = start + 1; i < size; i++) {
		const uint32 y = c[i];
		CHECKLIT(y);
		const LIT_ST yval = values[y];
		if (UNASSIGNED(yval) || (yval && satonly)) {
			best = y;
			pos = i;
			xval = yval;
			break;
		}
		const int ylevel = l2dl(y);
		assert(ylevel > -1);
		bool better;
		if (!xval && yval > 0) better = true;
		else if (xval > 0 && !yval) better = false;
		else if (!xval) {
			assert(!yval);
			better = (xlevel < ylevel);
		}
		else {
			assert(xval > 0);
			assert(yval > 0);
			assert(!satonly);
			better = (xlevel > ylevel);
		}
		if (better) {
			best = y;
			pos = i;
			xval = yval;
			xlevel = ylevel;
		}
	}
	if (!pos) return xval;
	c[start] = best;
	c[pos] = x;
	return xval;
}

void ParaFROST::sortClause(CLAUSE& c)
{
	int size = c.size();
	LIT_ST val = sortClause(c, 0, size, false);
	if (size > 2) sortClause(c, 1, size, val);
}