/***********************************************************************[autarky.cpp]
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

#include "solve.hpp"

using namespace ParaFROST;

// Autartic variables here are eliminated but 
// treated as frozen variables and get values
// like assigned variables (some exceptions
// had to be made in variable mapping though
// to make this acceptable

uint32 Solver::useAutarky(LIT_ST* autarkies)
{
	uint32 eliminated = 0;
	forall_variables(v) {
		uint32 lit = V2L(v);
		const LIT_ST val = autarkies[lit];
		if (UNASSIGNED(val)) continue;
		PFLOG2(4, "  eliminating autarkic literal %d", l2i(lit));
		markEliminated(ABS(lit));
		sp->value[lit] = val;
		sp->value[FLIP(lit)] = !val;
		eliminated++;
	}
	return eliminated;
}

uint32 Solver::autarkReasoning(LIT_ST* autarkies)
{
	assert(analyzed.empty());
	const LIT_ST* values = sp->value;
	const VSTATE* states = sp->vstate;
	// assign all variables with saved phases
	uint32 assigned = 0;
	forall_variables(v) {
		if (states[v].state) continue;
		LIT_ST pol = sp->psaved[v];
		assert(pol >= 0);
		const uint32 dec = V2DEC(v, pol);
		autarkies[dec] = 1, autarkies[FLIP(dec)] = 0;
		assigned++;
		PFLOG2(4, "   autarky assigned %d", l2i(dec));
	}
	// test autarky on large clauses and deduct unassigned
	forall_cnf(orgs, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		CLAUSE& c = cm[ref];
		if (c.binary()) continue;
		uint32 unassigned = propAutarkClause(false, ref, c, values, autarkies);
		if (unassigned) {
			assert(assigned >= unassigned);
			assigned -= unassigned;
			if (!assigned) return 0;
		}
	}
	PFLOG2(2, " Autarky %lld: initial large-clause autarky size is %d", stats.autarky.calls, assigned);
	// now propagate autarkies on binaries and leave out learnts
	forall_literal(lit) {
		if (inactive(lit)) continue;
		LIT_ST val = autarkies[lit];
		if (val > 0) continue;
		WL& ws = wt[FLIP(lit)];
		forall_watches(ws, i) {
			const WATCH w = *i;
			assert(w.binary());
			const C_REF ref = w.ref;
			assert(!cm.deleted(ref));
			if (cm[ref].learnt()) continue;
			const uint32 other = w.imp;
			assert(UNASSIGNED(values[other]));
			const LIT_ST otherval = autarkies[other];
			if (otherval > 0) continue;
			if (!val) {
				cancelAutark(true, lit, autarkies);
				val = UNDEFINED;
				assert(assigned > 0);
				if (!--assigned)
					break;
			}
			if (!otherval) {
				cancelAutark(true, other, autarkies);
				assert(assigned > 0);
				if (!--assigned)
					break;
			}
		}
		if (analyzed.size()) {
			const uint32 unassigned = propAutarky(values, autarkies);
			assert(unassigned <= assigned);
			assigned -= unassigned;
		}
		if (!assigned) return 0;
	}
	PFLOG2(2, " Autarky %lld: initial binary-clause autarky size is %d", stats.autarky.calls, assigned);
	// now propagate autarkies on large clauses
	forall_cnf(orgs, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		CLAUSE& c = cm[ref];
		if (c.binary()) continue;
		uint32 unassigned = propAutarkClause(true, ref, c, values, autarkies);
		if (unassigned) {
			unassigned += propAutarky(values, autarkies);
			assert(assigned >= unassigned);
			assigned -= unassigned;
			if (!assigned) break;
		}
		else if (!c.deleted()) {
			forall_clause(c, k) {
				const uint32 lit = *k;
				if (UNASSIGNED(values[lit]) && autarkies[lit] > 0) {
					PFLOG2(4, "   adding c(%zd) to watch list(%d)", ref, l2i(lit));
					wt[FLIP(lit)].push(WATCH(ref, c.size(), lit));
				}
			}
		}
	}
	PFLOG2(2, " Autarky %lld: final autarky size is %d", stats.autarky.calls, assigned);
	// eliminate autarkies
	if (assigned) {
		const uint32 eliminated = useAutarky(autarkies);
		assert(eliminated == assigned);
		stats.autarky.eliminated += eliminated;
		inf.maxMelted += eliminated;
		PFLOG2(2, " Autarky %lld: eliminated %d autarkic variables", stats.autarky.calls, eliminated);
		filterAutarky();
	}
	else
		detachClauses(true);
	return assigned;
}

void Solver::autarky()
{
	if (!opts.autarky_en || incremental) return;
	if (!UNSOLVED(cnfstate)) return;
	SLEEPING(sleep.autarky, opts.autarky_sleep_en);
	assert(!DL());
	stats.autarky.calls++;
	// at this step only all binaries are left
	// but learnts are ignored in reasoning 
	// which save time to reattach them again
	binarizeWT(true);
	LIT_ST* autarkies = pfmalloc<LIT_ST>(inf.nDualVars);
	memset(autarkies, UNDEFINED, inf.nDualVars);
	uint32 eliminated = autarkReasoning(autarkies);
	analyzed.clear();
	free(autarkies);
	attachNonBins(orgs, true);
	attachNonBins(learnts, true);
	if (retrail()) PFLOG2(2, " Propagation after autarky proved a contradiction");
	UPDATE_SLEEPER(this, autarky, eliminated);
	printStats(eliminated, 'k', CCYAN);
}

void Solver::filterAutarky()
{
	const VSTATE* states = sp->vstate;
	forall_literal(lit) {
		const bool litelim = MELTED(states[ABS(lit)].state);
		WL& ws = wt[FLIP(lit)];
		WATCH* j = ws;
		forall_watches(ws, i) {
			const WATCH w = *i;
			if (w.binary()) {
				const uint32 imp = w.imp;
				CHECKLIT(imp);
				const bool impelim = MELTED(states[ABS(imp)].state);
				if (litelim || impelim) {
					if (lit < imp) {
						const C_REF ref = w.ref;
						removeClause(cm[ref], ref);
					}
				}
				else
					*j++ = w;
			}
		}
		ws.resize(int(j - ws));
	}
}

uint32 Solver::propAutarky(const LIT_ST* values, LIT_ST* autarkies)
{
	assert(!DL());
	uint32 unassigned = 0;
	while (analyzed.size()) {
		const uint32 lit = analyzed.back();
		analyzed.pop();
		CHECKLIT(lit);
		assert(UNASSIGNED(autarkies[lit]));
		PFLOG2(4, "  propagating autarkic literal %d", l2i(lit));
		WL& ws = wt[FLIP(lit)];
		WATCH* j = ws;
		forall_watches(ws, i) {
			const WATCH w = *i;
			const C_REF ref = w.ref;
			CLAUSE& c = cm[ref];
			if (w.binary()) {
				*j++ = w;
				if (c.learnt()) continue;
				const uint32 imp = w.imp;
				CHECKLIT(imp);
				assert(imp != lit);
				LIT_ST impval = values[imp];
				if (impval > 0) continue;
				assert(UNASSIGNED(impval));
				if (!autarkies[imp]) {
					cancelAutark(true, imp, autarkies);
					unassigned++;
				}
			}
			else {
				assert(ABS(w.imp) == ABS(lit));
				unassigned += propAutarkClause(true, ref, c, values, autarkies);
			}
		}
		ws.resize(int(j - ws));
	}
	return unassigned;
}

inline void Solver::cancelAutark(const bool& add, const uint32& lit, LIT_ST* autarkies)
{
	const uint32 flit = FLIP(lit);
	assert(autarkies[flit] > 0);
	autarkies[lit] = autarkies[flit] = UNDEFINED;
	if (add) analyzed.push(flit);
	PFLOG2(4, "   autarky %d cancelled", l2i(flit));
}

inline uint32 Solver::propAutarkClause(const bool& add, const C_REF& ref, CLAUSE& c, const LIT_ST* values, LIT_ST* autarkies)
{
	assert(c.size() > 2);
	assert(c.original());
	assert(!c.deleted());
	PFLCLAUSE(4, c, "   autarky propagating clause");
	bool satisfied = false, falsified = false;
	forall_clause(c, k) {
		const uint32 other = *k;
		LIT_ST otherval = values[other];
		if (otherval > 0) {
			removeClause(c, ref);
			return 0;
		}
		if (otherval) { // unassigned
			otherval = autarkies[other];
			if (otherval > 0) satisfied = true;
			else if (!otherval) falsified = true;
		}
	}
	if (satisfied || !falsified) return 0;
	uint32 unassigned = 0;
	forall_clause(c, k) {
		const uint32 other = *k;
		const LIT_ST otherval = values[other];
		if (otherval) {
			assert(UNASSIGNED(otherval));
			if (!autarkies[other]) {
				cancelAutark(add, other, autarkies);
				assert(unassigned <= inf.maxVar);
				unassigned++;
			}
		}
	}
	return unassigned;
}

