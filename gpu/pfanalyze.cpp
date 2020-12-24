/***********************************************************************[pfanalyze.cpp]
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

inline void	ParaFROST::bumpClause(CLAUSE& c) {
	assert(c.learnt());
	assert(c.size() > 2);
	CL_ST used = c.usage();
	c.initTier3();
	if (c.keep()) return;
	int old_lbd = c.lbd();
	int new_lbd = calcLBD(c);
	if (new_lbd < old_lbd) { // update old LBD
		if (new_lbd <= opts.lbd_tier1) c.set_keep(1);
		else if (old_lbd > opts.lbd_tier2 && new_lbd <= opts.lbd_tier2) c.initTier2();
		c.set_lbd(new_lbd);
		PFLCLAUSE(4, c, " Bumping clause with lbd %d ", new_lbd);
	}
	else if (used && old_lbd <= opts.lbd_tier2) c.initTier2();
}

inline void	ParaFROST::bumpVariable(const uint32& v) {
	assert(v && v <= inf.maxVar);
	if (vsidsEnabled()) varBumpHeap(v);
	else varBumpQueue(v);
}

inline bool	ParaFROST::bumpReason(const uint32& lit) {
	assert(lit);
	assert(isFalse(lit));
	uint32 v = ABS(lit);
	if (sp->seen[v] || !sp->level[v]) return false;
	PFLOG2(4, "  bumping reason literal %d@%d", l2i(lit), sp->level[v]);
	analyzed.push(v);
	sp->seen[v] = ANALYZED_M;
	return true;
}

inline void	ParaFROST::bumpReasons(const uint32& lit, const int& depth) {
	assert(lit > 1);
	assert(depth > 0);
	uint32 v = ABS(lit);
	if (!sp->level[v]) return;
	C_REF r = sp->source[v];
	if (REASON(r)) {
		CLAUSE& c = cm[r];
		if (c.binary()) {
			uint32 other = c[0] ^ c[1] ^ lit;
			if (bumpReason(other) && depth >= 2)
				bumpReasons(FLIP(other), depth - 1);
		}
		else {
			for (uint32* k = c; k != c.end(); k++) {
				uint32 other = *k;
				if (other == lit || !bumpReason(other)) continue;
				if (depth >= 2) bumpReasons(FLIP(other), depth - 1);
			}
		}
	}
}

inline void	ParaFROST::bumpReasons() {
	assert(opts.bumpreason_depth > 0);
	for (uint32* k = learntC; k != learntC.end(); k++)
		bumpReasons(FLIP(*k), opts.bumpreason_depth);
}

inline void	ParaFROST::bumpVariables() {
	if (opts.bumpreason_en) bumpReasons();
	bool vsidsEn = vsidsEnabled();
	if (!vsidsEn) rSort(analyzed, ANALYZE_CMP(bumps), ANALYZE_RANK(bumps));
	for (uint32* v = analyzed; v != analyzed.end(); v++) bumpVariable(*v);
	if (vsidsEn) decayVarAct();
}

inline int ParaFROST::calcLBD() {
	int size = learntC.size();
	if (size > 2) {
		stats.marker++;
		register int lbd = 0;
		uint32* end = learntC.end();
		for (uint32* k = learntC + 1; k != end; k++) {
			int litLevel = l2dl(*k);
			if (sp->board[litLevel] != stats.marker) { sp->board[litLevel] = stats.marker; lbd++; }
		}
		return lbd;
	}
	else if (size == 2) return 1;
	assert(size == 1);
	return 0;
}

inline int ParaFROST::calcLBD(const CLAUSE& c) {
	register int lbd = 0;
	if (c.binary()) lbd = (l2dl(c[0]) != l2dl(c[1])) + 1;
	else {
		stats.marker++;
		for (int i = 0; i < c.size(); i++) {
			int litLevel = l2dl(c[i]);
			if (sp->board[litLevel] != stats.marker) { sp->board[litLevel] = stats.marker; lbd++; }
		}
	}
	return lbd;
}

inline void	ParaFROST::clearAnalyzed() {
	for (uint32* v = analyzed; v != analyzed.end(); v++) 
		sp->seen[*v] = 0;
	analyzed.clear();
}

inline void	ParaFROST::clearMinimized() {
	for (uint32* v = minimized; v != minimized.end(); v++) 
		sp->seen[*v] = 0;
	minimized.clear();
}

bool ParaFROST::chronoAnalyze()
{
	assert(conflict != NOREF);
	int currLevel = DL();
	int confLevel = 0;
	uint32 count = 0;
	uint32 forced = 0;
	CLAUSE& c = cm[conflict];
	for (uint32* k = c; k != c.end(); k++) {
		uint32 lit = *k;
		int litLevel = l2dl(lit);
		if (litLevel > confLevel) {
			confLevel = litLevel;
			forced = lit;
			count = 1;
		}
		else if (litLevel == confLevel) {
			count++;
			if (confLevel == currLevel && count > 1) break;
		}
	}
	assert(count);
	PFLOG2(3, " Found %d literals on conflict level %d", count, confLevel);
	if (!confLevel) { cnfstate = UNSAT; return false; }
	int size = c.size();
	for (int i = 0; i < 2; i++) {
		uint32 lit = c[i], maxLit = lit;
		int maxPos = i;
		int maxLevel = l2dl(maxLit);
		for (int j = i + 1; j < size; j++) {
			uint32 other = c[j];
			int otherLevel = l2dl(other);
			if (maxLevel >= otherLevel) continue;
			maxPos = j;
			maxLit = other;
			maxLevel = otherLevel;
			if (maxLevel == confLevel) break;
		}
		if (maxPos == i) continue;
		if (maxPos > 1) detachWatch(FLIP(lit), conflict);
		c[maxPos] = lit;
		c[i] = maxLit;
		if (maxPos > 1) attachWatch(maxLit, c[!i], conflict, size);
	}
	if (count == 1) {
		assert(forced > 1);
		backtrack(confLevel - 1);
		enqueueImp(forced, conflict);
		PFLCLAUSE(3, cm[conflict], " Forced %d@%d in conflicting clause", l2i(forced), l2dl(forced));
		conflict = NOREF;
		return true;
	}
	backtrack(confLevel);
	return false;
}

inline void	ParaFROST::analyzeLit(const uint32& lit, int& track) {
	assert(lit > 1);
	uint32 v = ABS(lit);
	if (sp->seen[v]) return;
	int current_level = DL(), litLevel = sp->level[v];
	if (litLevel) {
		assert(litLevel <= current_level);
		if (litLevel == current_level) track++;
		else if (litLevel < current_level) learntC.push(lit);
		analyzed.push(v);
		sp->seen[v] = ANALYZED_M;
	}
}

inline void	ParaFROST::analyzeReason(const C_REF& r, const uint32& parent, int& track) {
	CLAUSE& c = cm[r];
	PFLCLAUSE(4, c, "  analyzing %d %s", parent ? l2i(parent) : l2i(*c), parent ? "reason" : "conflict");
	if (c.binary()) {
		if (parent) analyzeLit((c[0] ^ c[1] ^ parent), track);
		else { analyzeLit(c[0], track), analyzeLit(c[1], track); }
	}
	else {
		if (c.learnt()) bumpClause(c);
		for (uint32* k = c; k != c.end(); k++) {
			uint32 lit = *k;
			if (lit != parent)
				analyzeLit(lit, track);
		}
	}
}

inline int ParaFROST::whereToJump()
{
	int current_level = DL(), bt_level = UNDEFINED;
	if (learntC.size() == 1) bt_level = 0;
	else if (learntC.size() == 2) bt_level = l2dl(learntC[1]);
	else {
		assert(learntC.size() > 2);
		assert(current_level > 1);
		int chrono_level = current_level - 1;
		uint32* maxPos = learntC + 1, maxLit = *maxPos, * end = learntC.end();
		bt_level = l2dl(maxLit);
		for (uint32* k = learntC + 2; k != end; k++) {
			uint32 lit = *k;
			int litLevel = l2dl(lit);
			if (bt_level >= litLevel) continue;
			bt_level = litLevel;
			maxLit = lit;
			maxPos = k;
			if (litLevel == chrono_level) break;
		}
		*maxPos = learntC[1];
		learntC[1] = maxLit;
		PFLLEARNT(this, 3);
	}
	assert(current_level > bt_level);
	if (opts.chrono_en && current_level - bt_level > opts.chrono_min) {
		bt_level = current_level - 1;
		stats.cbt++;
		PFLOG2(3, " Forced chronological backtracking to level %d", bt_level);
	}
	else if (opts.chrono_en && opts.chronoreuse_en) {
		uint32 best_v = 0, best_pos = 0;
		if (vsidsEnabled()) {
			HEAP_CMP hcmp(activity);
			for (uint32 i = dlevels[bt_level + 1]; i < trail.size(); i++) {
				const uint32 v = ABS(trail[i]);
				if (best_v && !hcmp(best_v, v)) continue;
				best_v = v;
				best_pos = i;
			}
		}
		else {
			for (uint32 i = dlevels[bt_level + 1]; i < trail.size(); i++) {
				const uint32 v = ABS(trail[i]);
				if (best_v && bumps[best_v] >= bumps[v]) continue;
				best_v = v;
				best_pos = i;
			}
		}
		assert(best_v);
		PFLOG2(4, " Best variable %d at trail position %d", best_v, best_pos);
		int old_bt_level = bt_level;
		while (bt_level < current_level - 1 && dlevels[bt_level + 1] <= best_pos) bt_level++;
		if (old_bt_level == bt_level) {
			stats.ncbt++;
			PFLOG2(4, " Default non-chronological backjumping to level %d to reuse trail", bt_level);
		}
		else {
			stats.cbt++;
			PFLOG2(4, " Forced chronological backtracking to level %d to reuse trail", bt_level);
		}
	}
	else stats.ncbt++;
	assert(bt_level != UNDEFINED);
	return bt_level;
}

C_REF ParaFROST::backjump(const int& bt_level) {
	assert(trail.size());
	// cancel old assignments up to backtrack level <bt_level>
	backtrack(bt_level);
	// add learnt clause & enqueue learnt decision
	if (opts.proof_en) {
		wrProof('a');
		wrProof(learntC.data(), learntC.size());
		wrProof(0);
	}
	if (learntC.size() == 1)
		enqueue(learntC[0]), stats.n_units++;
	else {
		C_REF r = newClause(learntC, LEARNT);
		enqueue(*learntC, bt_level, r);
		return r;
	}
	return NOREF;
}

void ParaFROST::analyze()
{
	assert(conflict != NOREF);
	assert(analyzed.empty());
	PFLOG2(3, " Analyzing conflict:");
	PFLTRAIL(this, 3);
	nConflicts++;
	if (opts.chrono_en && chronoAnalyze()) return;
	if (cnfstate == UNSAT) return;
	int current_level = DL();
	if (!current_level) { cnfstate = UNSAT; return; }
	if (stats.marker < 0) sp->clearBoard(), stats.marker = 0;
	learntC.clear();
	learntC.push(0);
	sp->learnt_lbd = UNDEFINED;
	uint32 parent = 0;
	int track = 0, index = trail.size();
	C_REF r = conflict;
	while (true) {
		assert(r != NOREF);
		analyzeReason(r, parent, track);
		// find next implication clause
		parent = 0;
		while (!parent) {
			assert(index > 0);
			uint32 lit = trail[--index], v = ABS(lit);
			if (sp->seen[v] && sp->level[v] == current_level) parent = lit;
		}
		if (!--track) break;
		assert(parent > 1);
		r = l2r(parent);
	}
	assert(learntC[0] == 0);
	learntC[0] = FLIP(parent);
	PFLLEARNT(this, 3);
	// calculate lbd value & learnt stats
	assert(sp->learnt_lbd == UNDEFINED);
	sp->learnt_lbd = calcLBD();
	assert(sp->learnt_lbd >= 0);
	assert(sp->learnt_lbd < learntC.size());
	PFLOG2(4, " LBD of learnt clause = %d", sp->learnt_lbd);
	lbdrest.update(sp->learnt_lbd);
	// bump variable activities
	bumpVariables();
	// minimize learnt clause 
	stats.max_lits += learntC.size();
	if (learntC.size() > 1) minimize();
	stats.tot_lits += learntC.size();
	// backjump control
	C_REF added = backjump(whereToJump());
	// next luby sequence
	if (lrn.stable) lubyrest.update();
	// clear 
	conflict = NOREF;
	clearAnalyzed(), clearMinimized();
	// subsume recent learnts 
	if (opts.learntsub_max && REASON(added)) subsumeLearnt(added);
	printStats(vsidsOnly() && nConflicts % opts.prograte == 0);
}

void ParaFROST::subsumeLearnt(const C_REF& l)
{
	if (learnts.size() < 2) return;
	assert(l != NOREF);
	CLAUSE& learnt = cm[l];
	markLits(learnt);
	int learntSize = learnt.size();
	int64 limit = lrn.subtried + opts.learntsub_max;
	C_REF* tail = learnts.end(), * head = learnts;
	while (tail != head && lrn.subtried++ <= limit) {
		C_REF t = *--tail;
		if (l == t) continue;
		CLAUSE& c = cm[t];
		if (c.deleted()) continue;
		int sub = learntSize;
		for (uint32* k = c; k != c.end(); k++) {
			if (subsumed(*k) && !--sub) {
				PFLCLAUSE(4, c, "  found subsumed learnt");
				removeClause(t);
				stats.n_learntsubs++;
				break;
			}
		}
	}
	unmarkLits(learnt);
}