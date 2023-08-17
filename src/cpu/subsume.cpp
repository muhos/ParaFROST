/***********************************************************************[subsume.cpp]
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

#include "solve.h"
#include "histogram.h"

using namespace ParaFROST;

struct SUBSUME_RANK {
	uint32 operator () (const CSIZE& a) const { return a.size; }
};

inline void Solver::strengthen(CLAUSE& c, const uint32& self) {
	CHECKLIT(self);
	assert(c.size() > 2);
	assert(unassigned(self));
	if (opts.proof_en) proof.shrinkClause(c, self);
	uint32 *j = c;
	forall_clause(c, i) {
		const uint32 lit = *i;
		assert(unassigned(lit));
		if (NEQUAL(lit, self)) 
			*j++ = *i;
	}
	assert(j + 1 == c.end());
	shrinkClause(c, 1);
	c.initTier3();
}

inline void	Solver::removeSubsumed(CLAUSE& c, const C_REF& cref, CLAUSE* s) {
	assert(s->size() <= c.size());
	assert(c.size() > 2);
	if (c.original() && s->learnt()) {
		s->markOriginal();
		const int ssize = s->size();
		stats.clauses.original++;
		assert(stats.clauses.learnt > 0);
		stats.clauses.learnt--;
		stats.literals.original += ssize;
		assert(stats.literals.learnt > 0);
		stats.literals.learnt -= ssize;
	}
	assert(!s->deleted());
	removeClause(c, cref);
}

inline bool Solver::subsumeCheck(CLAUSE* subsuming, uint32& self)
{
	stats.subsume.checks++;
	assert(!self);
	uint32 prev = 0;
	uint32* end = subsuming->end();
	bool good = true;
	for (uint32* i = subsuming->data(); good && i != end; i++) {
		uint32 lit = *i;
		CHECKLIT(lit);
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

inline CL_ST Solver::subsumeClause(CLAUSE& c, const C_REF& cref)
{
	assert(!c.deleted());
	assert(c.size() > 2);
	assert(keeping(c));
	markLits(c);   
	CLAUSE* s = NULL;
	uint32 self = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		if (!markedSubsume(lit)) continue;
		for (LIT_ST sign = 1; !s && sign >= 0; sign--) {
			assert(sign == 0 || sign == 1);
			const uint32 slit = sign ? FLIP(lit) : lit;
			BOL& others = bot[slit];
			forall_bol(others, o) {
				self = 0;
				const uint32 imp = *o;
				const LIT_ST marker = l2marker(imp), impSign = SIGN(imp);
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
				s = &subbin; // "always original"
				break;
			}
			if (s) break;
			WOL& wol = wot[slit];
			forall_wol(wol, i) {
				CLAUSE* d = cm.clause(*i);
				if (d->deleted()) continue;
				assert(d != &c);
				assert(d->size() <= c.size());
				if (subsumeCheck(d, self)) {
					s = d; // can be "original or learnt" 
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
		PFLCLAUSE(3, (*s), "  strengthened by ");
		strengthen(c, FLIP(self));
		return -1;
	}
	else {
		PFLCLAUSE(3, c, "  candidate ");
		PFLCLAUSE(3, (*s), "  subsumed by ");
		removeSubsumed(c, cref, s);
		return 1;
	}
}

void Solver::schedule2sub(BCNF& src)
{
	if (src.empty()) return;
	const LIT_ST* values = sp->value;
	uint32* hist = vhist.data();
	forall_cnf(src, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		const int size = c.size();
		if (size > opts.subsume_max_csize) continue;
		if (!keeping(c)) continue;
		// check marked-subsume/rooted literals
		int subsume = 0;
		bool rooted = false;
		forall_clause(c, k) {
			const uint32 lit = *k;
			CHECKLIT(lit);
			if (!UNASSIGNED(values[lit])) { rooted = true; break; }
			else if (markedSubsume(lit)) subsume++;
			assert(sp->value[lit] == UNDEFINED);
		}
		if (rooted || subsume < 2) continue;
		if (c.subsume()) stats.subsume.leftovers++;
		hist_clause(c, hist);
		scheduled.push(CSIZE(r, (uint32)size));
	}
}

bool Solver::subsumeAll()
{
	if (interrupted()) killSolver();
	assert(!DL());
	assert(inf.unassigned);
	assert(conflict == NOREF);
	assert(cnfstate != UNSAT);
	assert(wt.empty());
	SET_BOUNDS(sub_limit, subsume, subsume.checks, searchprops, 0);
	// schedule clauses
	BCNF shrunken;
	HIST_LCV_CMP clause_cmp(vhist);
	int64 checked = 0, subsumed = 0, strengthened = 0;
	vhist.resize(inf.nDualVars);
	memset(vhist, 0, sizeof(uint32) * inf.nDualVars);
	stats.subsume.leftovers = 0;
	schedule2sub(orgs);
	schedule2sub(learnts);
	if (scheduled.empty()) goto ending;
	scheduled.shrinkCap();
	radixSort(scheduled.data(), scheduled.end(), SUBSUME_RANK());
	if (!stats.subsume.leftovers) {
		forall_vector(CSIZE, scheduled, i) {
			assert(i->ref < cm.size());
			CLAUSE& c = cm[i->ref];
			if (c.size() > 2) c.markSubsume();
		}
	}
	PFLOG2(2, " Scheduled %d (%.2f %%) clauses for subsumption", scheduled.size(), percent((double)scheduled.size(), (double)maxClauses()));
	wot.resize(inf.nDualVars);
	bot.resize(inf.nDualVars);
	for (CSIZE* i = scheduled; i != scheduled.end(); i++) {
		if (interrupted()) break;
		if (stats.subsume.checks >= sub_limit) break;
		checked++;
		const C_REF r = i->ref;
		CLAUSE& c = cm[r];
		assert(!c.deleted());
		PFLCLAUSE(4, c, "  subsuming ");
		if (c.size() > 2 && c.subsume()) {
			c.initSubsume();
			CL_ST st = subsumeClause(c, r);
			if (st > 0) { subsumed++; continue; }
			if (st < 0) { shrunken.push(r); strengthened++; }
		}
		bool subsume = true, orgbin = (c.binary() && c.original());
		uint32 minlit = 0, minhist = 0;
		int minsize = 0;
		forall_clause(c, k) {
			const uint32 lit = *k;
			if (!markedSubsume(lit)) subsume = false;
			const int currentsize = orgbin ? bot[lit].size() : wot[lit].size();
			if (minlit && minsize <= currentsize) continue;
			const uint32 hist = vhist[lit];
			if (minlit && minsize == currentsize && hist <= minhist) continue;
			minlit = lit, minsize = currentsize, minhist = hist;
		}
		// current scheduled clause cannot subsume more clauses
		if (!subsume) continue; 
		// attach new occurrence
		if (minsize <= opts.subsume_max_occs) {
			if (orgbin) {
				assert(c.original());
				const uint32 other = c[0] ^ c[1] ^ minlit;
				assert(other != minlit);
				bot[minlit].push(other);
			}
			else {
				wot[minlit].push(r);
				Sort(c.data(), c.size(), clause_cmp);
			}
		}
	}
ending:
	PFLOG2(2, " Subsume %lld: removed %lld and strengthened %lld clauses", stats.subsume.calls, subsumed, strengthened);
	if (scheduled.size() == checked) sp->clearSubsume();
	forall_cnf(shrunken, i) { markSubsume(cm[*i]); }
	shrunken.clear(true), scheduled.clear(true), vhist.clear(true);
	wot.clear(true), bot.clear(true);
	stats.subsume.subsumed += subsumed;
	stats.subsume.strengthened += strengthened;
	return (subsumed || strengthened);
}

void Solver::subsume()
{
	if (!stats.clauses.original && !stats.clauses.learnt) return;
	rootify();
	assert(UNSOLVED(cnfstate));
	stats.subsume.calls++;
	printStats(1, '-', CORANGE0);
	wt.clear(true);
	bool success = subsumeAll();
	rebuildWT(opts.subsume_priorbins);
	filterOrg();
	if (retrail()) PFLOG2(2, " Propagation after subsume proved a contradiction");
	INCREASE_LIMIT(subsume, stats.subsume.calls, nlognlogn, true);
	printStats(success, 'u', CORANGE1);
}

void Solver::filterOrg() {
	C_REF* j = learnts;
	forall_cnf(learnts, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		if (cm[r].learnt()) *j++ = r;
		else orgs.push(r);
	}
	assert(j >= learnts);
	learnts.resize(uint32(j - learnts));
}