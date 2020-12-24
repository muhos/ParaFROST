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
	for (uint32* k = c; k != c.end(); k++)
		markLit(*k);
}

void ParaFROST::unmarkLits(CLAUSE& c) {
	for (uint32* k = c; k != c.end(); k++)
		unmarkLit(*k);
}

void ParaFROST::markSubsume(CLAUSE& c) {
	assert(keeping(c));
	for (uint32* k = c; k != c.end(); k++)
		sp->subsume[ABS(*k)] = 1;
}

bool ParaFROST::keeping(CLAUSE& c) {
	if (c.original()) return true;
	if (c.keep()) return true;
	if (c.lbd() > lrn.keptlbd) return false;
	if (c.size() > lrn.keptsize) return false;
	return true;
}

void ParaFROST::bumpShrunken(CLAUSE& c) {
	assert(c.learnt());
	assert(c.size() > 1);
	if (c.keep()) return;
	int old_lbd = c.lbd();
	int new_lbd = std::min(c.size() - 1, old_lbd);
	if (new_lbd >= old_lbd) return;
	if (new_lbd <= opts.lbd_tier1) c.set_keep(1);
	else if (old_lbd > opts.lbd_tier2 && new_lbd <= opts.lbd_tier2) c.initTier2();
	c.set_lbd(new_lbd);
	PFLCLAUSE(4, c, " Bumping shrunken clause with LBD %d ", new_lbd);
}

CL_ST ParaFROST::rooted(CLAUSE& c) {
	assert(!c.deleted());
	CL_ST st = UNDEFINED;
	for (uint32* k = c; k != c.end(); k++) {
		uint32 lit = *k;
		if (!l2dl(lit)) {
			LIT_ST val = sp->value[lit];
			assert(!UNASSIGNED(val));
			if (val > 0) return 1;
			if (!val) st = 0;
		}
	}
	return st;
}

inline int ParaFROST::removeRooted(CLAUSE& c)
{
	uint32 *i, *j, *end = c.end();
	for (i = c, j = i; i != end; i++) {
		uint32 lit = *i;
		if (!l2dl(lit)) { assert(!value(lit)); continue; }
		*j++ = lit;
	}
	return int(end - j);
}

void ParaFROST::shrinkClause(CLAUSE& c, const int& remLits)
{
	if (!remLits) return;
	c.shrink(remLits); // adjusts "pos" also
	cm.collect(remLits);
	assert(c.size() > 1);
	if (c.original()) { assert(inf.nLiterals), inf.nLiterals -= remLits; }
	else {
		// all shrunken must get bumped to "update" new binaries
		bumpShrunken(c);
		assert(inf.nLearntLits), inf.nLearntLits -= remLits;
	}
	if (keeping(c)) markSubsume(c);
}

void ParaFROST::shrinkClause(const C_REF& r)
{
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	assert(c.size() > 2);
	uint32* i, *end = c.end();
	int numNonFalse = 0;
	for (i = c; numNonFalse < 2 && i != end; i++)
		if (value(*i)) {
			assert(l2dl(*i));
			numNonFalse++;
		}
	if (numNonFalse < 2) return;
	// TODO:: implement proof here
	shrinkClause(c, removeRooted(c));
}

void ParaFROST::removeClause(const C_REF& r) {
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	assert(c.size() > 1);
	int sz = c.size();
	if (c.original()) { lrn.elim_marked += sz, inf.nLiterals -= sz; }
	else { 
		assert(c.learnt());
		inf.nLearntLits -= sz;
	}
	if (opts.proof_en) {
		wrProof('d');
		wrProof(c, sz);
		wrProof(0);
	}
	c.markDeleted(), cm.collect(r);
}

C_REF ParaFROST::newClause(const Lits_t& in_c, const CL_ST& type)
{
	C_REF r = cm.alloc(in_c);
	int sz = in_c.size();
	CLAUSE& c = cm[r];
	assert(sz > 1);
	assert(sz == c.size());
	assert(c[0] > 1 && c[1] > 1);
	assert(c[0] <= NOVAR && c[1] <= NOVAR);
	assert(c.keep());
	c.set_status(type);
	// attach to watch table
	attachWatch(r, c);
	// attach to database
	if (ISORG(type)) {
		orgs[inf.nClauses++] = r;
		inf.nLiterals += sz;
	}
	else {
		assert(sp->learnt_lbd > 0);
		int trim_lbd = sp->learnt_lbd > sz ? sz : sp->learnt_lbd;
		c.set_lbd(trim_lbd);
		c.set_keep(trim_lbd <= opts.lbd_tier1);
		if (sp->learnt_lbd <= opts.lbd_tier1) c.initTier1(), stats.n_glues++;			// Tier1
		else if (sp->learnt_lbd <= opts.lbd_tier2) c.initTier2();						// Tier2
		else c.initTier3();																// Tier3
		learnts.push(r);
		inf.nLearntLits += sz;
	}
	if (keeping(c)) markSubsume(c);
	return r;
}

bool ParaFROST::toClause(Lits_t& c, Lits_t& org, char*& str)
{
	assert(c.empty());
	assert(org.empty());
	uint32 v = 0, s = 0;
	bool satisfied = false;
	while ((v = toInteger(str, s)) != 0) {
		if (v > inf.maxVar) PFLOGE("too many variables");
		uint32 lit = V2DEC(v, s);
		org.push(lit);
		// checking literal
		LIT_ST marker = l2marker(lit);
		if (UNASSIGNED(marker)) {
			markLit(lit);
			LIT_ST val = value(lit);
			if (UNASSIGNED(val)) c.push(lit);
			else if (val) satisfied = true;
		}
		else if (marker != SIGN(lit)) satisfied = true; // tautology
	}
	for (uint32* k = org; k != org.end(); k++)
		unmarkLit(*k);
	if (satisfied) {
		if (opts.proof_en) wrProof('d'), wrProof(org, org.size()), wrProof(0);
	}
	else {
		if (!org.size()) { if (opts.proof_en) wrProof('0'); return false; }
		int newsize = c.size();
		if (newsize == 1) {
			LIT_ST val = value(*c);
			if (UNASSIGNED(val)) enqueueOrg(*c);
			else if (!val) return false;
		}
		else if (inf.nClauses + 1 > inf.nOrgCls) PFLOGE("too many clauses");
		else if (newsize) newClause(c);
		if (opts.proof_en && newsize < org.size()) {
			wrProof('a'), wrProof(c, newsize), wrProof(0);
			wrProof('d'), wrProof(org, org.size()), wrProof(0);
			org.clear();
		}
	}
	c.clear(), org.clear();
	return true;
}