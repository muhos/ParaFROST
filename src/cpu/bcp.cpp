/***********************************************************************[bcp.cpp]
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
using namespace pFROST;

bool ParaFROST::BCP()
{
	conflict = NOREF;
	const uint32 propsbefore = sp->propagated;
	LIT_ST* values = sp->value;
	bool isConflict = false;
	while (!isConflict && sp->propagated < trail.size()) {
		const uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		const int level = l2dl(assign);
		CHECKLIT(assign);
		PFLOG2(4, "  propagating %d@%d", l2i(assign), level);
		WL& ws = wt[assign];
		uint64 ticks = cacheLines(ws.size(), sizeof(WATCH));
		WATCH* i = ws, *j = i, * wend = ws.end();
		while (i != wend) {
			const WATCH w = *j++ = *i++;
			const uint32 imp = w.imp;
			CHECKLIT(imp);
			assert(imp != f_assign);
			const LIT_ST impVal = values[imp];
			if (impVal > 0) continue; // blocking literal
			const C_REF ref = w.ref;
			//=============== binary ================//
			if (w.binary()) {
				if (cm.deleted(ref)) { j--; continue; } // use cm 'stencil' to avoid dereferencing the clause
				if (impVal) enqueue(imp, level, ref);
				else { conflict = ref; break; }
			}
			//================ large =================//
			else {
				ticks++;
				if (cm.deleted(ref)) { j--; continue; } // use cm 'stencil' to avoid dereferencing the clause
				CLAUSE& c = cm[ref];
				assert(c.size() > 2);
				assert(c[0] != c[1]);
				const uint32 other = c[0] ^ c[1] ^ f_assign; // Thanks to CaDiCaL solver
				CHECKLIT(other);
				// check if first literal is true
				const LIT_ST otherVal = values[other];
				if (otherVal > 0) 
					(j - 1)->imp = other; // satisfied, replace "w.imp" with new blocking "other"
				else {
					// === search for (un)-assigned-1 literal to watch
					uint32* cmid = c.mid(), * cend = c.end();
					uint32* k = cmid, newlit = 0;
					LIT_ST _false_ = UNDEFINED;
					while (k != cend && (_false_ = !values[newlit = *k])) k++;
					assert(_false_ != UNDEFINED);
					if (_false_) {
						k = c + 2;
						assert(c.pos() <= c.size());
						while (k != cmid && (_false_ = !values[newlit = *k])) k++;
					}
					assert(k >= c + 2 && k <= c.end());
					c.set_pos(int(k - c)); // set new position
					// ======== end of search ========
					LIT_ST val = values[newlit];
					if (val > 0) // found satisfied new literal (update "imp")
						(j - 1)->imp = newlit; 
					else if (UNASSIGNED(val)) { // found new unassigned literal to watch
						c[0] = other;
						c[1] = newlit;
						*k = f_assign;
						assert(c[0] != c[1]);
						attachWatch(newlit, other, ref, c.size());
						j--; // remove j-watch from current assignment
						ticks++;
					}
					else if (UNASSIGNED(otherVal)) { // clause is unit
						assert(!val);
						enqueueImp(other, ref);
					}
					else { // clause is conflicting
						assert(!val);
						assert(!otherVal);
						PFLCONFLICT(this, 3, other);
						conflict = ref;
						break;
					}
				}
			}
		} // end of watches loop 
		while (i != wend) *j++ = *i++;
		ws.resize(int(j - ws));
		stats.searchticks += ticks;
		isConflict = NEQUAL(conflict, NOREF);
	} // end of trail loop
	stats.searchprops += sp->propagated - propsbefore;
	if (isConflict) sp->trailpivot = dlevels.back();
	else			sp->trailpivot = sp->propagated;
	return isConflict;
}

bool ParaFROST::BCPProbe() {
	assert(UNSOLVED(cnfstate));
	assert(DL() == 1);
	conflict = NOREF;
	bool isConflict = false;
	uint32 propagatedbin = sp->propagated;
	while (!isConflict && sp->propagated < trail.size()) {
		if (propagatedbin < trail.size())
			isConflict = propbinary(trail[propagatedbin++]);
		else
			isConflict = proplarge(trail[sp->propagated++], true);
	}
	return isConflict;
}

bool ParaFROST::BCPVivify() {
	assert(UNSOLVED(cnfstate));
	conflict = NOREF;
	bool isConflict = false;
	while (!isConflict && sp->propagated < trail.size()) {
		isConflict = proplarge(trail[sp->propagated++], false);
	}
	return isConflict;
}

inline bool ParaFROST::propbinary(const uint32& assign)
{
	CHECKLIT(assign);
	assert(DL() == 1);
	const int level = l2dl(assign);
	PFLOG2(4, "  propagating %d@%d in binaries", l2i(assign), level);
	LIT_ST* values = sp->value;
	WL& ws = wt[assign];
	forall_watches(ws, i) {
		const WATCH w = *i;
		if (w.binary()) {
			const uint32 imp = w.imp;
			CHECKLIT(imp);
			assert(imp != FLIP(assign));
			const LIT_ST impVal = values[imp];
			if (impVal > 0) continue; 
			const C_REF ref = w.ref;
			if (cm.deleted(ref)) continue;
			if (impVal) enqueue(imp, level, ref);
			else { conflict = ref; return true; }
		}
	}
	return false;
}

inline bool ParaFROST::proplarge(const uint32& assign, const bool& hyper)
{
	CHECKLIT(assign);
	const bool hbr = opts.probehbr_en && hyper;
	const int level = l2dl(assign);
	const uint32 f_assign = FLIP(assign);
	PFLOG2(4, "  propagating %d@%d in large clauses", l2i(assign), level);
	LIT_ST* values = sp->value;
	WL& ws = wt[assign];
	uint64 ticks = cacheLines(ws.size(), sizeof(WATCH)) + 1;
	WATCH* i = ws, * j = i, * wend = ws.end();
	while (i != wend) {
		const WATCH w = *j++ = *i++;
		const uint32 imp = w.imp;
		CHECKLIT(imp);
		assert(imp != f_assign);
		const LIT_ST impVal = values[imp];
		if (impVal > 0) continue; // blocking literal
		const C_REF ref = w.ref;
		//=============== binary ================//
		if (w.binary()) {
			if (cm.deleted(ref)) { j--; continue; } // use cm 'stencil' to avoid dereferencing the clause
			if (impVal) enqueue(imp, level, ref);
			else { conflict = ref; break; }
		}
		//================ large =================//
		else if (NEQUAL(ref, ignore)) {
			ticks++;
			const C_REF ref = w.ref;
			if (cm.deleted(ref)) { j--; continue; }
			CLAUSE& c = cm[ref];
			assert(c.size() > 2);
			assert(c[0] != c[1]);
			const uint32 other = c[0] ^ c[1] ^ f_assign;
			CHECKLIT(other);
			const LIT_ST otherVal = values[other];
			if (otherVal > 0)
				(j - 1)->imp = other;
			else {
				// === search for (un)-assigned-1 literal to watch
				uint32* cmid = c.mid(), * cend = c.end();
				uint32* k = cmid, newlit = 0;
				LIT_ST _false_ = UNDEFINED;
				while (k != cend && (_false_ = !values[newlit = *k])) k++;
				assert(_false_ != UNDEFINED);
				if (_false_) {
					k = c + 2;
					assert(c.pos() <= c.size());
					while (k != cmid && (_false_ = !values[newlit = *k])) k++;
				}
				assert(k >= c + 2 && k <= c.end());
				c.set_pos(int(k - c));
				// ======== end of search ========
				LIT_ST val = values[newlit];
				if (val > 0)
					(j - 1)->imp = newlit;
				else if (UNASSIGNED(val)) {
					c[0] = other;
					c[1] = newlit;
					*k = f_assign;
					assert(c[0] != c[1]);
					delayWatch(newlit, other, ref, c.size()), j--;
					ticks++;
				}
				else if (UNASSIGNED(otherVal)) {
					assert(!val);
					if (hbr) {
						const uint32 dom = hyper2Resolve(c, other);
						if (dom) {
							CHECKLIT(dom);
							PFLOG2(4, "  adding hyper binary resolvent(%d %d)", l2i(dom), l2i(other));
							assert(learntC.empty());
							learntC.push(dom);
							learntC.push(other);
							if (opts.proof_en) proof.addClause(learntC);
							const int csize = c.size();
							newHyper2(); // 'c' after this line is not valid and cm[ref] should be used if needed
							delayWatch(f_assign, other, ref, csize), j--;
						}
					}
					enqueueImp(other, ref);
				}
				else {
					assert(!val);
					assert(!otherVal);
					PFLCONFLICT(this, 3, other);
					conflict = ref;
					break;
				}
			}
		}
	} // end of watches loop
	while (i != wend) *j++ = *i++;
	ws.resize(int(j - ws));
	attachDelayed();
	assert(dwatches.empty());
	stats.probeticks += ticks;
	return NEQUAL(conflict, NOREF);
}