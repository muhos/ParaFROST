/***********************************************************************[shrink.cpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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

#include "solver.hpp"
using namespace ParaFROST;

void Solver::bumpShrunken(CLAUSE& c)
{
	assert(c.learnt());
	assert(c.size() > 1);
	if (c.keep()) return;
	if (c.hyper()) return;
	int old_lbd = c.lbd();
	int new_lbd = MIN(c.size() - 1, old_lbd);
	if (new_lbd >= old_lbd) return;
	if (new_lbd <= opts.lbd_tier1) c.set_keep(1);
	else if (old_lbd > opts.lbd_tier2 && new_lbd <= opts.lbd_tier2) c.initTier2();
	c.set_lbd(new_lbd);
	LOGCLAUSE(4, c, " Bumping shrunken clause with LBD %d ", new_lbd);
}

CL_ST Solver::rootedTop(CLAUSE& c)
{
	assert(!DL());
	assert(!c.deleted());
	CL_ST st = UNDEFINED;
	const LIT_ST* values = sp->value;
	forall_clause(c, k) {
		const uint32 lit = *k;
		CHECKLIT(lit);
		const LIT_ST val = values[lit];
		if (val > 0) return 1;
		if (!val) st = 0;
	}
	return st;
}

CL_ST Solver::rooted(CLAUSE& c)
{
	assert(!c.deleted());
	CL_ST st = UNDEFINED;
	const int* levels = sp->level;
	const LIT_ST* values = sp->value;
	forall_clause(c, k) {
		const uint32 lit = *k;
		CHECKLIT(lit);
		if (!levels[ABS(lit)]) {
			const LIT_ST val = values[lit];
			assert(!UNASSIGNED(val));
			if (val > 0) return 1;
			if (!val) st = 0;
		}
	}
	return st;
}

int Solver::removeRooted(CLAUSE& c)
{
	const int* levels = sp->level;
	uint32* j = c;
	forall_clause(c, i) {
		const uint32 lit = *i;
		CHECKLIT(lit);
		if (levels[ABS(lit)]) *j++ = lit;
		else assert(!sp->value[lit]);
	}
	return int(c.end() - j);
}

void Solver::shrinkClause(CLAUSE& c, const int& remLits)
{
	assert(remLits >= 0);
	if (!remLits) return;
	c.shrink(remLits); // adjusts "pos" also
	assert(c.size() > 1);
	if (c.learnt()) {
		bumpShrunken(c); // all shrunken must get bumped to "update" new binaries
		assert(stats.literals.learnt > 0);
		stats.literals.learnt -= remLits;
	}
	else {
		assert(c.original());
		assert(stats.literals.original > 0);
		stats.literals.original -= remLits;
	}
	if (keeping(c)) markSubsume(c);
	cm.collectLiterals(remLits);
}

void Solver::shrinkClause(const C_REF& r)
{
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	assert(c.size() > 2);
	LIT_ST* values = sp->value;
	uint32* i, * cend = c.end();
	int numNonFalse = 0;
	for (i = c; numNonFalse < 2 && i != cend; i++) {
		const uint32 lit = *i;
		CHECKLIT(lit);
		if (values[lit]) {
			assert(l2dl(lit));
			numNonFalse++;
		}
	}
	if (numNonFalse < 2) return;
	if (opts.proof_en) proof.shrinkClause(c);
	shrinkClause(c, removeRooted(c));
}

bool Solver::shrink()
{
	if (sp->simplified >= inf.maxFrozen) return false;
	sp->simplified = inf.maxFrozen;
	LOGN2(2, " Shrinking all clauses..");
	assert(trail.size());
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	assert(isPropagated());
	assert(!unassigned(trail.back()));
#ifdef STATISTICS
	stats.shrink.calls++;
	int64 beforeCls = maxClauses(), beforeLits = maxLiterals();
#endif
	shrink(orgs);
	shrink(learnts);
	assert(orgs.size() == stats.clauses.original);
	assert(learnts.size() == stats.clauses.learnt);
#ifdef STATISTICS
	LOGSHRINKALL(this, 2, beforeCls, beforeLits);
#else 
	LOGDONE(2, 5);
#endif
	return true;
}

void Solver::shrinkTop(const bool& conditional)
{
	if (conditional && sp->simplified >= inf.maxFrozen) return;
	assert(IS_UNSOLVED(cnfstate));
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	assert(isPropagated());
	LOGN2(2, " Shrinking all clauses on top level..");
	if (sp->simplified < inf.maxFrozen) sp->simplified = inf.maxFrozen;
#ifdef STATISTICS
	stats.shrink.calls++;
	int64 beforeCls = maxClauses(), beforeLits = maxLiterals();
#endif
	shrinkTop(orgs), shrinkTop(learnts);
	assert(orgs.size() == stats.clauses.original);
	assert(learnts.size() == stats.clauses.learnt);
#ifdef STATISTICS
	LOGSHRINKALL(this, 2, beforeCls, beforeLits);
#else 
	LOGDONE(2, 5);
#endif
}

void Solver::shrink(BCNF& cnf)
{
	if (cnf.empty()) return;
	C_REF* j = cnf;
	forall_cnf(cnf, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		assert(!c.moved());
		CL_ST st = rooted(c);
		if (st > 0) removeClause(c, r);
		else if (!st) {
			shrinkClause(r);
			*j++ = r;
		}
		else *j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

void Solver::shrinkTop(BCNF& cnf)
{
	if (cnf.empty()) return;
	assert(!DL());
	C_REF* j = cnf;
	forall_cnf(cnf, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		assert(!c.moved());
		CL_ST st = rootedTop(c);
		if (st > 0) removeClause(c, r);
		else if (!st) {
			shrinkClause(c, removeRooted(c));
			*j++ = r;
		}
		else *j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}