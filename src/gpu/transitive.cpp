/***********************************************************************[transitive.cpp]
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

void Solver::transiting(const uint32& src, const uint64& limit, uint64& removed, uint32& units)
{
	CHECKLIT(src);
	bool failed = false;
	assert(unassigned(src));
	LOG2(4, "  performing transitive reduction on literal %d", l2i(src));
	WL& sws = wt[src];
	stats.transitiveticks += cacheLines(sws.size(), sizeof(WATCH)) + 1;
	uVec1D& marked = minimized;
	forall_watches(sws, i) {
		assert(!failed);
		const WATCH sw = *i;
		if (!sw.binary()) break;
		const C_REF cref = sw.ref;
		if (cm.deleted(cref)) continue;
		const uint32 dest = sw.imp;
		CHECKLIT(dest);
		if (!unassigned(dest)) continue;
		CLAUSE& c = cm[cref];
		LOGCLAUSE(4, c, "  finding a transitive path to %d using", l2i(dest));
		const bool learnt = c.learnt();
		assert(marked.empty());
		assert(UNASSIGNED(l2marker(src)));
		markLit(src);
		marked.push(src);
		bool transitive = false;
		uint32 propagated = 0;
		while (!transitive && !failed && propagated < marked.size()) {
			const uint32 assign = marked[propagated++];
			CHECKLIT(assign);
			assert(l2marker(assign) == SIGN(assign));
			LOG2(4, "  transitively propagating %d in:", l2i(assign));
			WL& aws = wt[assign];
			stats.transitiveticks += cacheLines(aws.size(), sizeof(WATCH)) + 1;
			forall_watches(aws, j) {
				const WATCH aw = *j;
				if (!aw.binary()) break;
				const C_REF dref = aw.ref;
				if (dref == cref) continue;
				if (cm.deleted(dref)) continue;
				const CLAUSE& d = cm[dref];
				if (!learnt && d.learnt()) continue;
				LOGCLAUSE(4, d, "  ");
				const uint32 other = aw.imp;
				CHECKLIT(other);
				if (other == dest) { 
					transitive = true; 
					break;
				}
				else {
					const LIT_ST marker = l2marker(other);
					if (UNASSIGNED(marker)) {
						LOG2(4, "  transitive assign %d", l2i(other));
						markLit(other);
						marked.push(other);
					}
					else if (NEQUAL(marker, SIGN(other))) { 
						LOG2(4, "  found both %d and %d reachable", -l2i(other), l2i(other));
						failed = true;
						break;
					}						
				}
			}
		}
		forall_vector(uint32, marked, i) { unmarkLit(*i); }
		marked.clear();
		if (transitive) {
			LOGCLAUSE(4, c, "  found transitive clause");
			removeClause(c, cref);
			removed++;
		}
		if (failed) break;
		if (stats.transitiveticks > limit) break;
	}
	if (failed) {
		units++;
		LOG2(4, "  found failed literal %d during transitive reduction", l2i(src));
		enqueueUnit(FLIP(src));
		if (BCP()) {
			LOG2(2, " Propagation within transitive reduction proved a contradiction");
			learnEmpty();
		}
	}
}

void Solver::transitive()
{
	if (!cnfstate) return;
	if (!opts.transitive_en) return;
	assert(probed);
	assert(!DL());
	assert(isPropagated());
	sortWT();
	SET_BOUNDS(this, limit, transitive, transitiveticks, searchticks, 0);
	assert(last.transitive.literals < inf.maxDualVars);
	uint32 tried = 0, units = 0;
	uint64 removed = 0;
	while (stats.transitiveticks <= limit
		&& last.transitive.literals < inf.maxDualVars
		&& cnfstate && !interrupted()) {
		const uint32 lit = last.transitive.literals++;
		CHECKLIT(lit);
		if (active(lit) && unassigned(lit)) {
			tried++;
			transiting(lit, limit, removed, units);
		}
	}
	LOG2(2, " Transitive %lld: tried %d literals, removing %lld clauses and %d units",
		stats.probe.calls, tried, removed, units);
	if (last.transitive.literals == inf.maxDualVars)
		last.transitive.literals = 2;
	stats.transitive.probed += tried;
	stats.transitive.failed += units;
	stats.transitive.removed += removed;
	printStats((units || removed), 't', CVIOLET3);
}