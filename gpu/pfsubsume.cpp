/***********************************************************************[pfsubsume.cpp]
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

struct SUBSUME_RANK {
	size_t operator () (const CSIZE& a) const { return a.size; }
};

inline void ParaFROST::strengthen(CLAUSE& c, const uint32& self) {
	assert(self > 1 && self <= inf.nDualVars);
	assert(c.size() > 2);
	assert(unassigned(self));
	// TODO: implement proof
	uint32* i, * j, * end = c.end();
	for (i = c, j = i; i != end; i++) {
		uint32 lit = *i;
		assert(unassigned(lit));
		if (lit == self) continue;
		*j++ = *i;
	}
	assert(j + 1 == c.end());
	shrinkClause(c, 1);
	c.initTier3();
}

inline void	ParaFROST::removeSubsumed(CLAUSE& c, const C_REF& cref, CLAUSE* s, const C_REF& sref) {
	assert(s->size() <= c.size());
	assert(c.size() > 2);
	if (c.original() && s->learnt()) {
		assert(sref < NOREF);
		assert(inf.nLearntLits);
		s->set_status(ORIGINAL);
		int ssize = s->size();
		inf.nLiterals += ssize, inf.nLearntLits -= ssize;
	}
	assert(!s->deleted());
	removeClause(cref);
}

inline bool ParaFROST::subsumeCheck(CLAUSE* subsuming, uint32& self)
{
	stats.n_subchecks++;
	assert(!self);
	uint32 prev = 0;
	uint32* end = subsuming->end();
	bool good = true;
	for (uint32* i = subsuming->data(); good && i != end; i++) {
		uint32 lit = *i;
		assert(lit > 1);
		assert(unassigned(lit));
		*i = prev;
		prev = lit;
		LIT_ST marker = l2marker(lit);
		if (UNASSIGNED(marker)) good = false;
		else if (marker == SIGN(lit)) continue;
		else if (self) good = false;
		else self = lit;
	}
	assert(prev);
	assert(!**subsuming);
	**subsuming = prev;
	return good;
}

inline CL_ST ParaFROST::subsumeClause(const C_REF& cref, CLAUSE& c, BCNF& shrunken)
{
	assert(cm[cref] == c);
	assert(!c.deleted());
	assert(c.size() > 2);
	assert(keeping(c));
	markLits(c);   
	CLAUSE* s = NULL;
	C_REF sref = NOREF;
	uint32 self = 0;
	for (uint32* k = c; k != c.end(); k++) {
		uint32 lit = *k;
		if (!sp->subsume[ABS(lit)]) continue;
		for (LIT_ST sign = 1; !s && sign >= 0; sign--) {
			assert(sign == 0 || sign == 1);
			uint32 slit = sign ? FLIP(lit) : lit;
			BOL& others = bot[slit];
			for (uint32* o = others; o != others.end(); o++) {
				self = 0;
				uint32 imp = *o;
				LIT_ST marker = l2marker(imp), impSign = SIGN(imp);
				if (UNASSIGNED(marker)) continue;
				if (marker && sign) continue; // tautology
				if (marker == !impSign) {
					if (sign) continue; // tautology
					self = imp;
				}
				else if (sign) self = slit;
				assert(subbin.original());
				assert(subbin.binary());
				subbin[0] = slit, subbin[1] = imp;
				s = &subbin, sref = NOREF; // "always original"
				break;
			}
			if (s) break;
			WOL& wol = wot[slit];
			for (C_REF* i = wol; i != wol.end(); i++) {
				CLAUSE* d = cm.clause(*i);
				if (d->deleted()) continue;
				assert(d != &c);
				assert(d->size() <= c.size());
				if (subsumeCheck(d, self)) {
					s = d, sref = *i; // can be "original or learnt" 
					break;
				}
				else self = 0;
			}
		}
		if (s) break;
	}
	unmarkLits(c);
	if (!s) return 0;
	if (self) {
		PFLCLAUSE(3, c, "  candidate ");
		strengthen(c, FLIP(self));
		PFLCLAUSE(3, (*s), "  strengthened by ");
		shrunken.push(cref);
		return -1;
	}
	else {
		PFLCLAUSE(3, c, "  candidate ");
		removeSubsumed(c, cref, s, sref);
		PFLCLAUSE(3, (*s), "  subsumed by ");
		return 1;
	}
}

void ParaFROST::schedule(BCNF& src)
{
	if (src.empty()) return;
	for (C_REF* i = src; i != src.end(); i++) {
		C_REF r = *i;
		CLAUSE& c = cm[*i];
		int size = c.size();
		if (c.deleted()) continue;
		if (size > opts.subsume_max_csize) continue;
		if (!keeping(c)) continue;
		// check marked/rooted literals
		int subsume = 0;
		bool rooted = false;
		uint32* cend = c.end();
		for (uint32* k = c; k != cend; k++) {
			uint32 lit = *k;
			assert(lit > 1);
			if (!unassigned(lit)) { rooted = true; break; }
			else if (sp->subsume[ABS(lit)]) subsume++;
			assert(value(lit) == UNDEFINED);
		}
		if (rooted) { PFLCLAUSE(4, c, "  skipping rooted clause "); continue; }
		if (subsume < 2) { PFLCLAUSE(4, c, "  skipping less than %d added literals", subsume); continue; }
		if (c.subsume()) subleftovers++;
		scheduled.push(CSIZE(r, size));
		for (uint32* k = c; k != cend; k++) subhist[*k]++;
	}
}

bool ParaFROST::subsumeAll()
{
	if (interrupted()) killSolver();
	assert(!DL());
	assert(!satisfied());
	assert(conflict == NOREF);
	assert(cnfstate != UNSAT);
	assert(wt.empty());
	int64 sub_inc = stats.n_props;
	if (sub_inc < opts.subsume_min_checks) sub_inc = opts.subsume_min_checks;
	if (sub_inc > opts.subsume_max_checks) sub_inc = opts.subsume_max_checks;
	sub_inc = std::max(sub_inc, int64(maxActive()) << 1);
	PFLOG2(2, " Subsumption trials started with limit %lld", sub_inc);
	int64 sub_limit = stats.n_subchecks + sub_inc;
	// schedule clauses
	BCNF shrunken;
	SUBSUME_OCCURS_CMP clause_cmp(subhist);
	int64 checked = 0, subsumed = 0, strengthened = 0;
	subhist.resize(inf.nDualVars, 0);
	subleftovers = 0;
	schedule(orgs);
	schedule(learnts);
	if (scheduled.empty()) goto ending;
	scheduled.shrinkCap();
	radixSort(scheduled.data(), scheduled.end(), SUBSUME_RANK());
	if (!subleftovers) {
		for (CSIZE* i = scheduled; i != scheduled.end(); i++) {
			assert(i->ref < cm.size());
			CLAUSE& c = cm[i->ref];
			if (c.size() > 2) c.markSubsume();
		}
	}
	PFLOG2(2, " Scheduled %d (%.2f %%) clauses for subsumption", scheduled.size(), 100.0 * scheduled.size() / (double)maxClauses());
	wot.resize(inf.nDualVars);
	bot.resize(inf.nDualVars);
	for (CSIZE* i = scheduled; i != scheduled.end(); i++) {
		if (interrupted()) break;
		if (stats.n_subchecks >= sub_limit) break;
		checked++;
		C_REF r = i->ref;
		CLAUSE& c = cm[r];
		assert(!c.deleted());
		PFLCLAUSE(4, c, " Subsuming ");
		if (c.size() > 2 && c.subsume()) {
			c.initSubsume();
			CL_ST st = subsumeClause(r, c, shrunken);
			if (st > 0) { subsumed++; continue; }
			if (st < 0) strengthened++;
		}
		bool subsume = true, orgbin = (c.binary() && c.original());
		uint32 minlit = 0, minhist = 0;
		int minsize = 0;
		for (uint32* k = c; k != c.end(); k++) {
			uint32 lit = *k;
			if (!sp->subsume[ABS(lit)]) subsume = false;
			const int currentsize = orgbin ? bot[lit].size() : wot[lit].size();
			if (minlit && minsize <= currentsize) continue;
			const uint32 hist = subhist[lit];
			if (minlit && minsize == currentsize && hist <= minhist) continue;
			minlit = lit, minsize = currentsize, minhist = hist;
		}
		// current scheduled clause cannot subsume more clauses
		if (!subsume) continue; 
		// attach new occurrence
		if (minsize <= opts.subsume_min_occs) {
			if (orgbin) {
				PFLOG2(4, " watching %d with %d current original binary and total %d histogram", l2i(minlit), minsize, minhist);
				assert(c.original());
				uint32 other = c[0] ^ c[1] ^ minlit;
				assert(other != minlit);
				bot[minlit].push(other);
			}
			else {
				PFLOG2(4, " watching %d with %d current and total %d histogram", l2i(minlit), minsize, minhist);
				wot[minlit].push(r);
				Sort(c.data(), c.size(), clause_cmp);
			}
		}
	}
ending:
	PFLOG2(2, " Subsumed %lld and strengthened %lld clauses", subsumed, strengthened);
	if (scheduled.size() == checked) sp->clearSubsume();
	for (C_REF* r = shrunken; r != shrunken.end(); r++) markSubsume(cm[*r]);
	shrunken.clear(true), scheduled.clear(true), subhist.clear(true);
	wot.clear(true), bot.clear(true);
	stats.n_allsubsumed += subsumed;
	stats.n_allstrengthened += strengthened;
	return (subsumed || strengthened);
}

void ParaFROST::subsume()
{
	if (orgs.empty() && learnts.empty()) return;
	stats.n_subcalls++;
	backtrack();
	if (BCP()) { cnfstate = UNSAT; return; }
	printStats(1, '-', CORANGE0);
	wt.clear(true);
	bool success = subsumeAll();
	wt.resize(inf.nDualVars);
	rebuildWT(opts.priorbins_en);
	filter(learnts, orgs, ORIGINAL);
	assert(sp->propagated == trail.size());
	int64 current_inc = int64(scale(opts.subsume_inc * (stats.n_subcalls + 1.0)));
	lrn.subsume_conf_max = nConflicts + current_inc;
	PFLOG2(2, " subsume limit increased to %lld conflicts by a weight %lld", lrn.subsume_conf_max, current_inc);
	printStats(success, 's', CORANGE1);
}