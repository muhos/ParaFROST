/***********************************************************************[analyze.cpp]
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

#include "solve.hpp"
using namespace ParaFROST;

inline void	Solver::bumpReason(const uint32& lit) 
{
	CHECKLIT(lit);
	assert(isFalse(lit));
	const uint32 v = ABS(lit);
	if (!sp->level[v] || ANALYZED(sp->seen[v])) return;
	LOG2(4, "  bumping reason literal %d@%d", l2i(lit), sp->level[v]);
	assert(!sp->seen[v]);
	sp->seen[v] = ANALYZED_M;
	analyzed.push(v);
}

inline void	Solver::bumpReasons(const uint32& lit)
{
	CHECKLIT(lit);
	const uint32 v = ABS(lit);
	const C_REF r = sp->source[v];
	if (REASON(r)) {
		CLAUSE& c = cm[r];
		if (c.binary()) {
			const uint32 other = c[0] ^ c[1] ^ lit;
			bumpReason(other);
		}
		else {
			forall_clause(c, k) {
				const uint32 other = *k;
				if (NEQUAL(other, lit)) 
					bumpReason(other);
			}
		}
	}
}

inline void	Solver::bumpReasons()
{
	assert(!probed);
	forall_clause(learntC, k) {
		sp->seen[ABS(*k)] = ANALYZED_M;
	}
	forall_clause(learntC, k) {
		assert(l2dl(*k) > 0);
		bumpReasons(FLIP(*k));
	}
}

inline void	Solver::bumpVariable(const uint32& v) 
{
	CHECKVAR(v);
	assert(!sp->vstate[v].state);
	if (vsidsEnabled()) varBumpHeap(v);
	else varBumpQueue(v);
}

inline void	Solver::bumpVariables() 
{
	if (opts.bumpreason_en) bumpReasons();
	const bool vsidsEn = vsidsEnabled();
	if (!vsidsEn) rSort(analyzed, QUEUE_CMP(bumps), QUEUE_RANK(bumps));
	forall_vector(uint32, analyzed, a) {
		bumpVariable(*a);
	}
	if (vsidsEn) last.vsids.boost();
}

bool Solver::chronoAnalyze()
{
	assert(conflict != NOREF);
	const int level = DL();
	int conflictlevel = 0;
	uint32 count = 0;
	uint32 forced = 0;
	CLAUSE& c = cm[conflict];
	forall_clause(c, k) {
		const uint32 lit = *k;
		const int litlevel = l2dl(lit);
		if (litlevel > conflictlevel) {
			conflictlevel = litlevel;
			forced = lit;
			count = 1;
		}
		else if (litlevel == conflictlevel) {
			count++;
			if (conflictlevel == level && count > 1) break;
		}
	}
	assert(count);
	LOG2(3, "  found %d literals on conflict level %d", count, conflictlevel);
	if (!conflictlevel) { learnEmpty(); return false; }
	const int size = c.size();
	for (int i = 0; i < 2; i++) {
		const uint32 lit = c[i];
		uint32 maxLit = lit;
		int maxPos = i;
		int maxLevel = l2dl(maxLit);
		for (int j = i + 1; j < size; j++) {
			const uint32 other = c[j];
			int otherLevel = l2dl(other);
			if (maxLevel >= otherLevel) continue;
			maxPos = j;
			maxLit = other;
			maxLevel = otherLevel;
			if (maxLevel == conflictlevel) break;
		}
		if (maxPos == i) continue;
		if (maxPos > 1) detachWatch(FLIP(lit), conflict);
		c[maxPos] = lit;
		c[i] = maxLit;
		if (maxPos > 1) attachWatch(maxLit, c[!i], conflict, size);
	}
	if (count == 1) {
		assert(forced > 1);
		backtrack(conflictlevel - 1);
		enqueueImp(forced, conflict);
		LOGCLAUSE(3, cm[conflict], "  forced %d@%d in conflicting clause", l2i(forced), l2dl(forced));
		conflict = NOREF;
		return true;
	}
	backtrack(conflictlevel);
	return false;
}

void Solver::analyze()
{
	assert(conflict != NOREF);
	assert(learntC.empty());
	assert(analyzed.empty());
	assert(lbdlevels.empty());
	LOG2(3, " Analyzing conflict%s:", probed ? " during probing" : "");
	LOGTRAIL(this, 4);
	bool conflictchanged = true;
	while (conflictchanged) {
		stats.conflicts++;
		if (opts.chrono_en && chronoAnalyze()) return;
		if (!cnfstate) return;
		if (!DL()) { learnEmpty(); return; }
		// find first-UIP
		conflictchanged = finduip();
	}
	// update luby restart
	if (stable) lubyrest.update();
	if (!probed) {
		// minimize learnt clause 
		stats.minimize.before += learntC.size();
		if (learntC.size() > 1) minimize();
		if (learntC.size() > 1
			&& sp->learntLBD <= opts.minimize_lbd
			&& learntC.size() <= opts.minimize_min)
			minimizebin();
		stats.minimize.after += learntC.size();
		// update lbd restart mechanism
		lbdrest.update(stable, sp->learntLBD);
		// bump variable activities
		bumpVariables();
	}
	else assert(learntC.size() == 1);
	// backjump control
	C_REF added = backjump();
	// clear 
	conflict = NOREF;
	learntC.clear();
	clearLevels();
	clearAnalyzed();
	// subsume recent learnts 
	if (opts.learntsub_max && REASON(added)) subsumeLearnt(added);
	if (vsidsOnly()) printStats(stats.conflicts % opts.prograte == 0);
}

void Solver::subsumeLearnt(const C_REF& l)
{
	if (learnts.size() < 2) return;
	assert(l != NOREF);
	CLAUSE& learnt = cm[l];
	markLits(learnt);
	const int learntsize = learnt.size();
	uint64 trials = stats.subtried + opts.learntsub_max;
	C_REF* tail = learnts.end(), * head = learnts;
	uint32 nsubsumed = 0;
	while (tail != head && stats.subtried++ <= trials) {
		const C_REF t = *--tail;
		if (l == t) continue;
		if (cm.deleted(t)) continue;
		CLAUSE& c = cm[t];
		int sub = learntsize;
		forall_clause(c, k) {
			if (subsumed(*k) && !--sub) {
				nsubsumed++;
				removeClause(c, t);
				break;
			}
		}
	}
	unmarkLits(learnt);
	stats.subsume.learntfly += nsubsumed;
}