/***********************************************************************[pfbcp.cpp]
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

bool ParaFROST::BCP()
{
	conflict = NOREF;
	uint32 propsBefore = sp->propagated;
	bool noConflict = true;
	while (noConflict && sp->propagated < trail.size()) {
		uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		int assign_dl = l2dl(assign);
		assert(assign > 1);
		PFLOG2(3, " Propagating %d@%d", l2i(assign), assign_dl);
		PFLBCPS(this, 4, assign);
		WL& ws = wt[assign];
		if (ws.size()) {
			WATCH* i = ws, *j = i, *wend = ws.end();
			while (i != wend) {
				const WATCH& w = *j++ = *i++;
				uint32 imp = w.imp;
				assert(imp > 1 && imp < NOVAR);
				assert(imp != f_assign);
				LIT_ST impVal = value(imp);
				if (impVal > 0) continue; // blocking literal
				//=============== binary ================//
				if (w.binary()) {
					C_REF r = w.ref;
					if (cm[r].deleted()) { j--; continue; } 
					if (!impVal) { conflict = r, noConflict = false; }
					else enqueue(imp, assign_dl, r);
				}
				//================ large =================//
				else {
					if (!noConflict) break; // binary conflict found
					C_REF r = w.ref;
					CLAUSE& c = cm[r];
					if (c.deleted()) { j--; continue; }
					assert(c.size() > 2);
					uint32 other = c[0] ^ c[1] ^ f_assign; // Thanks to Cadical solver
					// check if first literal is true
					LIT_ST otherVal = value(other);
					if (otherVal > 0) 
						(j - 1)->imp = other; // satisfied, replace "w.imp" with new blocking "other"
					else {
						// === search for (un)-assigned-1 literal to watch
						uint32* cmid = c.mid(), *cend = c.end();
						uint32* k = cmid, newlit = 0;
						LIT_ST _false_ = UNDEFINED; 
						while (k != cend && (_false_ = isFalse(newlit = *k))) k++;
						assert(_false_ != UNDEFINED);
						if (_false_) {
							k = c + 2;
							assert(c.pos() <= c.size());
							while (k != cmid && (_false_ = isFalse(newlit = *k))) k++;
						}
						assert(k >= c + 2 && k <= c.end());
						c.set_pos(int(k - c)); // set new position
						// ======== end of search ========
						LIT_ST val = value(newlit);
						if (val > 0) { 
							// found satisfied new literal (keep the watch & replace "imp")
							(j - 1)->imp = newlit;
						} 
						else if (UNASSIGNED(val)) {
							// found new unassigned literal to watch
							c[0] = other;
							c[1] = newlit;
							*k = f_assign;
							attachWatch(newlit, f_assign, r, c.size());
							j--; // remove j-watch from current assignment
						}
						else if (UNASSIGNED(otherVal)) {
							assert(!val);
							// clause is unit
							enqueueImp(other, r);
							if (opts.chrono_en) {
								int otherLevel = l2dl(other);
								if (otherLevel > assign_dl) {
									uint32* maxPos, *e = c.end(), maxLit = 0;
									for (maxPos = c + 2; maxPos != e; maxPos++)
										if (l2dl(maxLit = *maxPos) == otherLevel)
											break;
									assert(maxPos < e);
									assert(maxLit > 1);
									*maxPos = f_assign;
									c[0] = other;
									c[1] = maxLit;
									attachWatch(maxLit, other, r, c.size());
									j--; // remove j-watch from current assignment
								}
							}
						}
						else {
							// clause is conflicting
							assert(!val);
							assert(!otherVal);
							PFLCONFLICT(this, 3, other);
							conflict = r, noConflict = false;
							break;
						}
					}
				} 
			} // end of watches loop 
			if (j != i) {
				while (i != wend) *j++ = *i++;
				ws.resize(int(j - ws));
			}
		} // end of "ws" size check
		PFLBCPE(this, 4, assign);
	} // end of trail loop

	stats.n_props += (sp->propagated - propsBefore);
	if (noConflict) sp->trailpivot = sp->propagated;
	else sp->trailpivot = dlevels.back();
	return !noConflict;
}