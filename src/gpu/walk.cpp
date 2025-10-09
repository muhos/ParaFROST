/***********************************************************************[walk.cpp]
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

/**** cannot be touched ****/
#define WALKBASE  0.5
#define BREAKMAX  1.0
#define BREAKMIN  4.94e-324
#define EXPONENTS 1075 
/***************************/

// defined as fixed array to increase its chance
// to being kept most of the time in cache
double lookup[EXPONENTS];
double minscore = -1;

#define LOOKUPSCORE(BREAKS) \
	(BREAKS < EXPONENTS) ? lookup[BREAKS] : minscore;

WALK::WALK() :
	value(NULL)
	, limit(0)
	, initial(0)
	, current(0)
	, minimum(0)
	, flipped(0)
	, nclauses(0)
	, best(0)
{
	double breakscore = BREAKMAX, base = WALKBASE;
	for (int i = 0; breakscore; i++, breakscore *= base) {
		assert(i < EXPONENTS);
		lookup[i] = breakscore;
	}
	assert(breakscore == 0);
	minscore = lookup[EXPONENTS - 1];
}

inline void WALK::destroy()
{
	orgs.clear(true);
	unsat.clear(true);
	cinfo.clear(true);
	trail.clear();
	scores.clear();
	if (value != NULL)
		free(value), value = NULL;
	limit = 0;
	best = flipped = 0;
	minimum = current = 0;
	initial = nclauses = 0;
}

void Solver::walk()
{
	assert(!DL());
	assert(IS_UNSOLVED(cnfstate));
	assert(isPropagated());
	stats.walk.calls++;
	wt.clear(true);
	shrinkTop(true);
	bot.resize(inf.maxDualVars); // used as occurrence table for indexing 'tracker.cinfo'
	walkinit();
	walkassign();
	if (cnfstate && walkschedule()) {
		walking();
		walkstop();
	}
	else 
		tracker.destroy();
	bot.clear(true);
	rebuildWT(opts.walk_priorbins);
	if (retrail()) LOG2(2, " Propagation after walk proved a contradiction");
}

void Solver::walkinit()
{
	assert(lookup[0] == BREAKMAX);
	assert(minscore >= 0 && minscore <= BREAKMIN);
	assert(tracker.orgs.empty());
	assert(tracker.cinfo.empty());
	assert(tracker.unsat.empty());
	assert(tracker.trail.empty());
	assert(tracker.scores.empty());
	assert(!bot.empty());
	assert(stats.clauses.original < UINT32_MAX);
	const uint32 nclauses = (uint32)stats.clauses.original;
	tracker.nclauses = nclauses;
	tracker.orgs.resize(nclauses);
	tracker.cinfo.resize(nclauses);
	tracker.value = pfmalloc<LIT_ST>(inf.maxDualVars);
	assert(tracker.value != NULL);
	memset(tracker.value, UNDEFINED, inf.maxDualVars);
}

bool Solver::walkschedule()
{
	assert(tracker.nclauses == stats.clauses.original);
	assert(tracker.orgs.size() == tracker.nclauses);
	assert(tracker.cinfo.size() == tracker.nclauses);
	// schedule broken clauses
	const bool assuming = incremental && assumptions.size();
	const LIT_ST* orgvalues = sp->value;
	const LIT_ST* values = tracker.value;
	BCNF& refs = tracker.orgs;
	uVec1D& unsatclauses = tracker.unsat;
	Vec<CINFO>& cinfo = tracker.cinfo;
	uint32 scheduled = 0;
	forall_cnf(orgs, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		CLAUSE& c = cm[ref];
		forall_clause(c, k) {
			if (orgvalues[*k] > 0) { // already satisfied at root
				removeClause(c, ref);
				break;
			}
		}
		if (c.deleted()) continue;
		assert(scheduled < tracker.nclauses);
		bool notassumed = true;
		uint32 satisfied = 0;
		forall_clause(c, k) {
			const uint32 lit = *k;
			const LIT_ST val = values[lit];
			if (UNASSIGNED(val)) continue;
			bot[lit].push(scheduled);
			if (val) satisfied++;
			else if (assuming && notassumed && !iassumed(ABS(lit)))
				notassumed = false;
		}
		CINFO& info = cinfo[scheduled];
		if (!satisfied) {
			if (assuming && !notassumed) return false;
			info.unsatidx = unsatclauses.size();
			unsatclauses.push(scheduled);
		}
		info.size = satisfied;
		refs[scheduled++] = ref;
	}
	tracker.initial = unsatclauses.size();
	tracker.minimum = tracker.current = tracker.initial;
	LOG2(2, " Walk %lld: found initial %d unsatisfied large clauses (%.2f%%)",
		stats.walk.calls, tracker.initial, percent(tracker.initial, scheduled));
	return true;
}

void Solver::walkstop()
{
	assert(tracker.minimum <= tracker.initial);
	if (tracker.minimum == tracker.initial) {
		LOG2(2, " Walk %lld: no improvement as unsatisfied clauses did not change", stats.walk.calls);
		last.rephase.type = 0;
	}
	else {
		if (tracker.best && NEQUAL(tracker.best, NOVAR))
			saveTrail(tracker.value, false);
		stats.walk.improved++;
		LOG2(2, " Walk %lld: improved phases till %d unsatisfied clauses", stats.walk.calls, tracker.minimum);
		last.rephase.type = WALKPHASE;
		printStats(1, 'w', CCYAN);
	}
	tracker.destroy();
}

uint32 Solver::promoteLit()
{
	assert(assumptions.empty());
	assert(tracker.current == tracker.unsat.size());
	const uint32 unsatpos = random.irand() % tracker.current;
	const uint32 infoidx = tracker.unsat[unsatpos];
	assert(infoidx < tracker.nclauses);
	const C_REF cref = tracker.orgs[infoidx];
	CLAUSE& c = cm[cref];
	assert(tracker.scores.empty());
	LIT_ST* values = tracker.value;
	double sum = 0, score = 0;
	uint32 promoted = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		if (UNASSIGNED(values[lit])) continue;
		promoted = lit;
		const uint32 breaks = breakValue(lit);
		score = LOOKUPSCORE(breaks);
		assert(score > 0);
		tracker.scores.push(score);
		sum += score;
	}
	assert(sum);
	CHECKLIT(promoted);
	const double drand = random.drand();
	assert(drand >= 0 && drand < 1);
	const double threshold = sum * drand;
	double* scores = tracker.scores;
	sum = 0, score = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		if (UNASSIGNED(values[lit])) continue;
		score = *scores++;
		sum += score;
		if (threshold < sum) {
			promoted = lit;
			break;
		}
	}
	CHECKLIT(promoted);
	LOG2(4, "  promoted literal %d with score %.6f", l2i(promoted), score);
	tracker.scores.clear();
	return promoted;
}

uint32 Solver::ipromoteLit()
{
	assert(tracker.current == tracker.unsat.size());
	const uint32 unsatpos = random.irand() % tracker.current;
	const uint32 infoidx = tracker.unsat[unsatpos];
	assert(infoidx < tracker.nclauses);
	const C_REF cref = tracker.orgs[infoidx];
	CLAUSE& c = cm[cref];
	assert(tracker.scores.empty());
	LIT_ST* values = tracker.value;
	double sum = 0, score = 0;
	uint32 promoted = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		if (UNASSIGNED(values[lit]) || iassumed(ABS(lit))) continue;
		promoted = lit;
		const uint32 breaks = breakValue(lit);
		score = LOOKUPSCORE(breaks);
		assert(score > 0);
		tracker.scores.push(score);
		sum += score;
	}
	assert(sum);
	CHECKLIT(promoted);
	const double drand = random.drand();
	assert(drand >= 0 && drand < 1);
	const double threshold = sum * drand;
	double* scores = tracker.scores;
	sum = 0, score = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		if (UNASSIGNED(values[lit]) || iassumed(ABS(lit))) continue;
		score = *scores++;
		sum += score;
		if (threshold < sum) {
			promoted = lit;
			break;
		}
	}
	CHECKLIT(promoted);
	LOG2(4, "  promoted literal %d with score %.6f", l2i(promoted), score);
	tracker.scores.clear();
	return promoted;
}

void Solver::walkstep()
{
	tracker.flipped++;
	const uint32 lit = assumptions.empty() ? promoteLit() : ipromoteLit();
	LOG2(4, "  walking step %d by flipping promoted literal %d", tracker.flipped, l2i(lit));
	LIT_ST* values = tracker.value;
	const LIT_ST val = values[lit];
	assert(!val);
	values[lit] = !val;
	values[FLIP(lit)] = val;
	makeClauses(lit);
	breakClauses(lit);
	tracker.current = tracker.unsat.size();
	if (tracker.best < NOVAR) {
		const uint32 limit = (inf.maxVar >> 2) + 1;
		assert(limit < NOVAR);
		assert(tracker.best <= tracker.trail.size());
		if (tracker.trail.size() < limit)
			tracker.trail.push(lit);
		else if (tracker.best) {
			saveTrail(values, true);
			tracker.trail.push(lit);
			assert(tracker.trail.size() <= UINT_MAX);
		}
		else {
			tracker.trail.clear();
			tracker.best = NOVAR;
		}
	}
	else
		assert(tracker.trail.empty());
	if (tracker.current < tracker.minimum) {
		tracker.minimum = tracker.current;
		if (tracker.best == NOVAR)
			saveAll(values);
		else {
			assert(tracker.trail.size() < NOVAR);
			tracker.best = tracker.trail.size();
		}
	}
	LOG2(4, "  remained %d unsatsified clauses after flipping", tracker.current);
}

void Solver::walking()
{
	SET_BOUNDS(this, walk_limit, walk, walk.checks, searchticks, 2 * maxClauses());
	tracker.limit = walk_limit;
	while (tracker.minimum
		&& tracker.limit > stats.walk.checks
		&& !interrupted()) {
		assert(tracker.current);
		walkstep();
	}
	stats.walk.minimum = tracker.minimum;
	stats.walk.flipped += tracker.flipped;
	LOG2(2, " Walk %lld: got minimum %lld unsatisfied clauses after %d flips", stats.walk.calls, stats.walk.minimum, tracker.flipped);
}

inline void Solver::walkassign()
{
	const bool targeting = useTarget();
	const VSTATE* states = sp->vstate;
	LIT_ST* values = tracker.value;
	if (assumptions.empty()) {
		forall_variables(v) {
			if (states[v].state) continue;
			const uint32 dec = makeAssign(v, targeting);
			values[dec] = 1, values[FLIP(dec)] = 0;
		}	
	}
	else {
		const LIT_ST* orgvalues = sp->value;
		forall_clause(assumptions, k) {
			const uint32 dec = *k;
			CHECKLIT(dec);
			const LIT_ST orgval = orgvalues[dec];
			if (!orgval || !values[dec]) {
				LOG2(2, " Walk %lld: abort due to conflicting assumption %d", stats.walk.calls, l2i(dec));
				ianalyze(FLIP(dec));
				cnfstate = UNSAT;
				return;
			}
			if (!states[ABS(dec)].state && UNASSIGNED(orgval)) {
				values[dec] = 1, values[FLIP(dec)] = 0;
			}
		}
		forall_variables(v) {
			if (states[v].state || iassumed(v)) continue;
			const uint32 dec = makeAssign(v, targeting);
			values[dec] = 1, values[FLIP(dec)] = 0;
		}
	}
}

inline uint32 Solver::breakValue(const uint32& lit)
{
	CHECKLIT(lit);
	assert(!tracker.value[lit]);
	const Vec<CINFO>& cinfo = tracker.cinfo;
	BOL& list = bot[FLIP(lit)];
	uint32 breaks = 0;
	forall_bol(list, i) {
		stats.walk.checks++;
		assert(*i < tracker.nclauses);
		if (cinfo[*i].size == 1) breaks++;
	}
	return breaks;
}

inline void Solver::saveAll(const LIT_ST* values)
{
	assert(tracker.trail.empty());
	assert(tracker.best == NOVAR);
	forall_variables(v) {
		const LIT_ST val = values[V2L(v)];
		if (!UNASSIGNED(val))
			sp->psaved[v] = !val;
	}
	tracker.best = 0;
}

inline void Solver::saveTrail(const LIT_ST* values, const bool& keep)
{
	uint32* t = tracker.trail, * best = t + tracker.best;
	for (uint32* k = t; k != best; k++) {
		const uint32 lit = *k;
		LIT_ST val = values[lit];
		assert(val >= 0);
		sp->psaved[ABS(lit)] = val;
	}
	if (keep) {
		uint32* s = t, * end = tracker.trail.end();
		for (uint32* k = best; k != end; k++)
			*s++ = *k;
		assert(uint32(end - s) == tracker.best);
		assert(uint32(s - t) == (tracker.trail.size() - tracker.best));
		tracker.trail.resize(uint32(s - t));
		tracker.best = 0;
	}
}

inline bool Solver::popUnsat(const uint32& infoidx, const uint32& unsatidx, Vec<CINFO>& cinfo)
{
	assert(cinfo[infoidx].unsatidx == unsatidx);
	assert(tracker.current && tracker.current == tracker.unsat.size());
	const uint32 lastinfoidx = tracker.unsat[--tracker.current];
	tracker.unsat.pop();
	if (NEQUAL(infoidx, lastinfoidx)) {
		assert(lastinfoidx < tracker.nclauses);
		CINFO& lastinfo = cinfo[lastinfoidx];
		assert(lastinfo.unsatidx == tracker.current);
		assert(unsatidx < lastinfo.unsatidx);
		lastinfo.unsatidx = unsatidx;
		tracker.unsat[unsatidx] = lastinfoidx;
		return true;
	}
	return false;
}

inline void Solver::breakClauses(const uint32& lit)
{
	CHECKLIT(lit);
	const uint32 negated = FLIP(lit);
	assert(!tracker.value[negated]);
	LOGN2(4, "   breaking satisfied clauses with negated literal %d..", l2i(negated));
	BOL& list = bot[negated];
	uVec1D& unsatclauses = tracker.unsat;
	Vec<CINFO>& cinfo = tracker.cinfo;
	forall_bol(list, i) {
		stats.walk.checks++;
		const uint32 infoidx = *i;
		assert(infoidx < tracker.nclauses);
		CINFO& info = cinfo[infoidx];
		assert(info.size);
		if (!--info.size) {
			info.unsatidx = unsatclauses.size();
			unsatclauses.push(infoidx);
		}
	}
	LOGDONE(4, 5);
}

inline void Solver::makeClauses(const uint32& lit)
{
	CHECKLIT(lit);
	assert(tracker.value[lit] > 0);
	LOGN2(4, "   reducing unsatisfied clauses with literal %d..", l2i(lit));
	BOL& list = bot[lit];
	Vec<CINFO>& cinfo = tracker.cinfo;
	uint64& checks = stats.walk.checks;
	forall_bol(list, i) {
		checks++;
		const uint32 infoidx = *i;
		assert(infoidx < tracker.nclauses);
		CINFO& info = cinfo[infoidx];
		assert(info.size < NOVAR);
		if (!info.size++ && popUnsat(infoidx, info.unsatidx, cinfo))
			checks++;
	}
	LOGDONE(4, 5);
}