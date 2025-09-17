/***********************************************************************[probe.cpp]
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

struct PROBE_QUEUE_CMP {
	const VSTATE* states;
	const Vec<uint64>& bumped;
	PROBE_QUEUE_CMP(const VSTATE* _states, const Vec<uint64>& _bumped) :
		states(_states), bumped(_bumped) {}
	uint64 operator () (const uint32& a, const uint32& b) const {
		const uint32 av = ABS(a), bv = ABS(b);
		const bool pa = states[av].probe, pb = states[bv].probe;
		if (!pa && pb) return true;
		if (pa && !pb) return false;
		return bumped[av] < bumped[bv];
	}
};

struct PROBE_HEAP_CMP {
	const VSTATE* states;
	const Vec<double>& act;
	PROBE_HEAP_CMP(const VSTATE* _states, const Vec<double>& _act) :
		states(_states), act(_act) {}
	bool operator () (const uint32& a, const uint32& b) const {
		const uint32 av = ABS(a), bv = ABS(b);
		const bool pa = states[av].probe, pb = states[bv].probe;
		if (!pa && pb) return true;
		if (pa && !pb) return false;
		const double xact = act[av], yact = act[bv];
		if (xact < yact) return true;
		if (xact > yact) return false;
		return av < bv;
	}
};

void Solver::analyzeFailed(const uint32& failed)
{
	assert(IS_UNSOLVED(cnfstate));
	assert(DL() == 1);
	assert(conflict != NOREF);
	CHECKLIT(failed);
	analyze();
	assert(!DL());
	assert(sp->value[failed] <= 0);
	if (unassigned(failed)) {
		const uint32 unit = FLIP(failed);
		LOG2(3, "  found unassigned failed probe %d", l2i(unit));
		enqueueUnit(unit);
	}
	if (BCP()) {
		LOG2(2, "  failed probe %d proved a contradiction", l2i(failed));
		learnEmpty();
	}
}

void Solver::scheduleProbes() 
{
	assert(probes.empty());
	vhist.resize(inf.nDualVars);
	memset(vhist, 0, sizeof(uint32) * inf.nDualVars);
	histBins(orgs);
	histBins(learnts);
	VSTATE* states = sp->vstate;
	uint32 count[2] = { 0 , 0 };
	forall_variables(v) {
		if (states[v].state) continue;
		const uint32 p = V2L(v), n = NEG(p);
		const uint32 poss = vhist[p], negs = vhist[n];
		if (poss && negs) continue;
		if (!poss && !negs) continue;
		uint32 probe = negs ? p : n;
		LOG2(4, "  scheduling probe %d with binary occurs %d", l2i(probe), vhist[FLIP(probe)]);
		probes.push(probe);
		assert(states[v].probe <= 1);
		count[states[v].probe]++;
	}
	assert(probes.size() == count[0] + count[1]);
	LOG2(2, "  scheduled %d (%d prioritized) probes %.2f%%", probes.size(), count[1], percent(probes.size(), maxActive()));
}

uint32 Solver::nextProbe() 
{
	while (!probes.empty()) {
		const uint32 probe = probes.back();
		CHECKLIT(probe);
		probes.pop();
		if (inactive(probe)) continue;
		return probe;
	}
	return 0;
}

void Solver::FLE()
{
	if (!cnfstate) return;
	assert(!DL());
	assert(isPropagated());
	SLEEPING(sleep.probe, opts.probe_sleep_en);
	SET_BOUNDS(this, probe_limit, probe, probeticks, searchticks, nlogn(maxActive()));
	VSTATE* states = sp->vstate;
	ignore = NOREF;
	int64 old_hypers = stats.binary.resolvents;
	uint32 probe = 0, currprobed = 0, currfailed = 0;
	for (int round = 1; round <= opts.probe_min; round++) {
		scheduleProbes();
		if (probes.size()) {
			if (stable) Sort(probes, PROBE_HEAP_CMP(states, activity));
			else Sort(probes, PROBE_QUEUE_CMP(states, bumps));
		}
		else {
			LOG2(2, "  no candidates found to probe");
			break;
		}
		memset(vhist, 0, sizeof(uint32) * inf.nDualVars);
		stats.probe.rounds++;
		currprobed = currfailed = 0;
		while ((probe = nextProbe())
			&& stats.probeticks < probe_limit
			&& cnfstate && !interrupted())
		{
			assert(!DL());
			assert(unassigned(probe));
			assert(isPropagated());
			unmarkProbe(probe);
			if (currfailed && vhist[probe] == currfailed)
				continue;
			currprobed++;
			enqueueDecision(probe);
			const uint32 propagated = sp->propagated;
			if (BCPProbe()) {
				currfailed++;
				analyzeFailed(probe);
			}
			else {
				assert(DL() == 1);
				assert(isPropagated());
				for (uint32 i = propagated; i < trail.size(); i++)
					vhist[trail[i]] = currfailed;
				backtrack();
			}
		}
		stats.probe.probed += currprobed;
		stats.probe.failed += currfailed;
		if (!cnfstate) {
			LOG2(2, "  probing proved a contradiction");
			break;
		}
		LOG2(2, " Probe round %lld: probed %d, finding %d failed literals and %lld hyper binary resolvents",
			stats.probe.rounds, currprobed, currfailed, stats.binary.resolvents - old_hypers);
		const uint32 remained = probes.size();
		if (remained) {
			LOG2(2, "  probing hit limit at round %d with %d remaining probes", round, remained);
			bool prioritized = false;
			for (uint32 i = 0; !prioritized && i < probes.size(); i++) {
				if (markedProbe(probes[i])) prioritized = true;
			}
			if (!prioritized) {
				LOG2(2, "  prioritizing remaining %d probes at round %d", remained, round);
				while (!probes.empty()) {
					const uint32 probe = probes.back();
					CHECKLIT(probe);
					probes.pop();
					assert(!markedProbe(probe));
					markProbe(probe);
				}
			}
			break;
		}
		if (!currfailed) break;
	}
	vhist.clear(true);
	probes.clear(true);
	int64 hypers = stats.binary.resolvents - old_hypers;
	const bool success = currfailed || hypers;
	UPDATE_SLEEPER(this, probe, success);
	printStats(success, 'f', CVIOLET4);
}

void Solver::probe()
{
	rootify();
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	stats.probe.calls++;
	printStats(1, '-', CVIOLET0);
	assert(!probed);
	probed = true;
	const uint32 before = maxActive();
	ELS(true);
	ternary();  
	debinary();
	transitive();
	FLE();
	vivify();
	ELS(false);
	probed = false;
	const uint32 after = maxActive();
	const uint32 removed = before - after;
	assert(removed >= 0);
	if (removed) LOG2(2, " Probe call %lld removed %d variables %.2f%%", stats.probe.calls, removed, percent(removed, before));
	else LOG2(2, " Probe %lld: did not remove any variables", stats.probe.calls);
	INCREASE_LIMIT(this, probe, stats.probe.calls, nlogn, true);
	last.probe.reduces = stats.reduces + 1;
}