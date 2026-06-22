/***********************************************************************[elimbcp.cu]
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
#include "grid.cuh"

using namespace ParaFROST;

#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)

// Per-variable assignment: 0 = unset, 1 = true, 2 = false
#define BCP_UNSET 0u
#define BCP_TRUE  1u
#define BCP_FALSE 2u
// Control slots
#define BCP_CURR  0 // current frontier size
#define BCP_NEXT  1 // next frontier size
#define BCP_CONFL 2 // conflict flag
#define BCP_LEVEL 3 // BFS level (device-side frontier ping-pong)

_PFROST_IN_D_
uint32 bcp_lit_val(const uint32* __restrict__ state, const uint32& lit)
{
	const uint32 s = state[ABS(lit)];
	if (s == BCP_UNSET) return BCP_UNSET;
	const bool sat = SIGN(lit) ? (s == BCP_FALSE) : (s == BCP_TRUE);
	return sat ? BCP_TRUE : BCP_FALSE;
}

__global__
void bcp_seed_k(
	uint32* __restrict__ state,
	uint32* __restrict__ front,
	uint32* __restrict__ ctrl,
	Byte* __restrict__   eliminated,
	const uint32* __restrict__ units,
	const uint32 nunits)
{
	for_parallel_x(i, nunits) {
		const uint32 u = units[i];
		const uint32 v = ABS(u);
		const uint32 desired = SIGN(u) ? BCP_FALSE : BCP_TRUE;
		const uint32 old = atomicCAS(&state[v], BCP_UNSET, desired);
		if (old == BCP_UNSET) {
			MARKFORCED(eliminated[v]);
			front[atomicAggInc(&ctrl[BCP_CURR])] = FLIP(u); // falsified lit to scan
		}
		else if (old != desired) ctrl[BCP_CONFL] = 1;
	}
}

__global__
void bcp_propagate_k(
	const CNF* __restrict__ cnf,
	const OT* __restrict__  ot,
	uint32* __restrict__    state,
	Byte* __restrict__      eliminated,
	uint32* __restrict__    front_a,
	uint32* __restrict__    front_b,
	uint32* __restrict__    ctrl,
	cuVecU* __restrict__    trail)
{
	if (ctrl[BCP_CONFL]) return;
	const uint32 level = ctrl[BCP_LEVEL];
	const uint32* __restrict__ front_curr = (level & 1) ? front_b : front_a;
	uint32* __restrict__       front_next = (level & 1) ? front_a : front_b;
	const uint32 cur_size = ctrl[BCP_CURR];
	for_parallel_y(fi, cur_size) {
		const uint32 ulit = front_curr[fi]; // a falsified literal
		const OL& list = (*ot)[ulit];
		const uint32 lsz = list.size();
		for_parallel_x(j, lsz) {
			const SCLAUSE& c = (*cnf)[list[j]];
			if (!c.deleted()) {
				uint32 unit = 0; int nunset = 0; bool sat = false;
				forall_clause_const(c, k) {
					const uint32 ve = bcp_lit_val(state, *k);
					if (ve == BCP_TRUE) { sat = true; break; }
					if (ve == BCP_UNSET) { unit = *k; if (++nunset > 1) break; }
				}
				if (!sat) {
					if (!nunset) ctrl[BCP_CONFL] = 1; // all falsified -> conflict
					else if (nunset == 1) {
						const uint32 v = ABS(unit);
						const uint32 desired = SIGN(unit) ? BCP_FALSE : BCP_TRUE;
						const uint32 old = atomicCAS(&state[v], BCP_UNSET, desired);
						if (old == BCP_UNSET) {
							MARKFORCED(eliminated[v]);
							trail->insertAggr(unit);
							front_next[atomicAggInc(&ctrl[BCP_NEXT])] = FLIP(unit);
						}
						else if (old != desired) ctrl[BCP_CONFL] = 1;
					}
				}
			}
		}
	}
}

__global__
void bcp_advance_k(uint32* __restrict__ ctrl)
{
	ctrl[BCP_CURR] = ctrl[BCP_NEXT];
	ctrl[BCP_NEXT] = 0;
	ctrl[BCP_LEVEL]++;
}

__global__
void bcp_apply_k(CNF* __restrict__ cnfptr, const uint32* __restrict__ state)
{
	CNF& cnf = *cnfptr;
	for_parallel_x(i, cnf.size()) {
		const S_REF r = cnf.ref(i);
		SCLAUSE& c = cnf[r];
		if (!c.deleted() && (c.original() || c.learnt())) {
			uint32* dst = c.data();
			uint32 sig = 0; int newsz = 0; bool sat = false;
			forall_clause(c, k) {
				const uint32 lit = *k;
				const uint32 ve = bcp_lit_val(state, lit);
				if (ve == BCP_TRUE) { sat = true; break; }
				if (ve == BCP_FALSE) continue; // drop falsified literal
				*dst++ = lit, sig |= MAPHASH(lit), newsz++;
			}
			if (sat) c.markDeleted();
			else { c.set_sig(sig); c.resize(newsz); }
		}
	}
}

#endif

void Solver::createOTHost(HOT& hot)
{
	assert(hcnf != NULL);
	assert(vars->nUnits);
	assert(reallocFailed());
	for (uint32 i = 0; i < hcnf->size(); i++) {
		const S_REF r = hcnf->ref(i);
		SCLAUSE& c = (*hcnf)[r];
		if (c.learnt() || c.original()) {
			assert(c.size());
			forall_clause(c, k) {
				CHECKLIT(*k);
				hot[*k].push(r);
			}
		}
	}
}

bool Solver::propFailed()
{
	if (vars->nUnits) {
		SYNCALL; // sync 'cacheCNF'
		assert(reallocFailed());
		assert(hcnf != NULL);
		assert(inf.numClauses == hcnf->size());
		HOT hot(inf.maxDualVars);
		createOTHost(hot);
		// start proping on host
		nForced = sp->propagated;
		assert(vars->cachedUnits != NULL);
		LIT_ST* values = sp->value;
		uint32* t = vars->cachedUnits + vars->nUnits;
		for (uint32* u = vars->cachedUnits; u != t; u++) {
			const uint32 unit = *u;
			CHECKLIT(unit);
			const LIT_ST val = values[unit];
			if (UNASSIGNED(val)) enqueueDevUnit(unit);
			else if (!val) return false; // early conflict detection
		}
		if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
		else LOGN2(2, " Propagating pending forced units..");
		int64 bclauses = inf.numClauses, bliterals = inf.numLiterals;
		CNF& hcnf = *this->hcnf;
		while (sp->propagated < trail.size()) { // propagate units
			uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
			CHECKLIT(assign);
			forall_occurs(hot[assign], i) {
				hcnf[*i].markDeleted();
			}
			forall_occurs(hot[f_assign], i) {
				SCLAUSE& c = hcnf[*i];
				assert(c.size());
				if (c.deleted()) continue;
				if (propClause(values, f_assign, c))
					c.markDeleted();
				else {
					assert(c.size()); // cannot be empty at this point
					if (c.size() == 1) {
						const uint32 unit = *c;
						CHECKLIT(unit);
						if (UNASSIGNED(values[unit])) enqueueUnit(unit);
						else { learnEmpty(); return false; }
					}
				}
			}
		}
		LOGDONE(2, 5);
		nForced = sp->propagated - nForced;
		LOGREDALLHOST(this, 2, "BCP Reductions");
		countAll(1);
		inf.numClauses = inf.numClausesSurvived;
		inf.numLiterals = inf.numLiteralsSurvived;
		stats.units.forced += nForced;
		stats.sigma.all.clauses += bclauses - int64(inf.numClauses);
		stats.sigma.all.literals += bliterals - int64(inf.numLiterals);
		vars->nUnits = nForced = 0;
	}
	return true;
}

#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
bool Solver::prop()
{
	if (!vars->nUnits) return true; // nothing forced -> nothing to propagate
	LOG2(2, " Propagating %d forced units on GPU..", vars->nUnits);
	const uint32 nvars = inf.maxVar + 1;
	const size_t words = size_t(nvars) * 3 + 8; // state + front_a + front_b + ctrl
	cuPool bcp;
	if (!cumm.allocDynamic(bcp, words * sizeof(uint32), "BCP")) { simpstate = OTALLOC_FAIL; return true; }
	uint32* d_state = (uint32*)bcp.mem;
	uint32* front_a = d_state + nvars;
	uint32* front_b = front_a + nvars;
	uint32* d_ctrl  = front_b + nvars;
	CHECK(cudaMemsetAsync(d_state, 0, nvars * sizeof(uint32)));
	CHECK(cudaMemsetAsync(d_ctrl, 0, 8 * sizeof(uint32)));
	{
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(vars->nUnits, nThreads, 0);
		bcp_seed_k<<<nBlocks, nThreads>>>(d_state, front_a, d_ctrl, vars->eliminated, vars->unitsData, vars->nUnits);
		LASTERR("BCP seed failed");
	}
	// This acts like a BFS over the forced-unit closure, 
	// where the frontier is the falsified literals.
	const dim3 block2D(32, 8);          // x: scan a literal's occurrences, y: frontier literals
	const dim3 grid2D(32, 64);
	uint32 ctrl_host[3] = { 0, 0, 0 };
	for (int batch = 1; ; batch = MIN(batch << 1, 64)) {
		for (int b = 0; b < batch; b++) {
			bcp_propagate_k<<<grid2D, block2D>>>(cnf, ot, d_state, vars->eliminated, front_a, front_b, d_ctrl, vars->units);
			bcp_advance_k<<<1, 1>>>(d_ctrl);
		}
		LASTERR("BCP propagate failed");
		CHECK(cudaMemcpy(ctrl_host, d_ctrl, 3 * sizeof(uint32), cudaMemcpyDeviceToHost));
		if (ctrl_host[BCP_CONFL] || !ctrl_host[BCP_CURR]) break;
	}
	if (ctrl_host[BCP_CONFL]) { cumm.freeDynamic(bcp); learnEmpty(); return false; }
	// Clean up CNF
	{
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.numClauses, nThreads, 0);
		bcp_apply_k<<<nBlocks, nThreads>>>(cnf, d_state);
		LASTERR("BCP apply failed");
		SYNCALL;
	}
	// Enqueue units on host to update solver state
	uint32 ntrail = 0;
	CHECK(cudaMemcpy(&ntrail, vars->unitsSize, sizeof(uint32), cudaMemcpyDeviceToHost));
	assert(ntrail <= inf.maxVar);
	if (ntrail) {
		CHECK(cudaMemcpy(sp->tmpstack, vars->unitsData, ntrail * sizeof(uint32), cudaMemcpyDeviceToHost));
		for (uint32 i = 0; i < ntrail; i++) {
			const uint32 unit = sp->tmpstack[i];
			CHECKLIT(unit);
			if (i < vars->nUnits) enqueueDevUnit(unit); // BVE-origin: freeze
			else                  enqueueUnit(unit);    // derived: learn
		}
		sp->propagated = trail.size(); // device already propagated these
	}
	cumm.freeDynamic(bcp);
	stats.units.forced += ntrail;
	LOGREDALL(this, 2, "BCP Reductions");
	countAll();
	inf.numClauses = inf.numClausesSurvived, inf.numLiterals = inf.numLiteralsSurvived;
	if (inf.numLiterals) {
		flattened = false; // force flattenCNF to re-copy the post-BCP CNF
		if (!reallocOT(streams[0])) { simpstate = OTALLOC_FAIL; return true; }
	}
	createOTAsync(false, streams[0]);
	CHECK(cudaMemsetAsync(vars->unitsSize, 0, sizeof(uint32), streams[1])); // reset device units counter
	vars->nUnits = 0; // histogram already rebuilt above; varReorder uses the fresh d_hist
	SYNC(streams[0]); SYNC(streams[1]);
	return true;
}
#else
bool Solver::prop()
{
	if (!enqueueCached(streams[3])) { learnEmpty(); return false; }
	LIT_ST* values = sp->value;
	while (sp->propagated < trail.size()) {
		uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		CHECKLIT(assign);
		MARKFORCED(vars->eliminated[ABS(assign)]);
		OL& ol = (*ot)[assign];
		OL& f_ol = (*ot)[f_assign];
		CNF& cnf = *this->cnf;
		forall_occurs(f_ol, i) { // reduce unsatisfied 
			SCLAUSE& c = cnf[*i];
			assert(c.size());
			if (c.deleted()) continue;
			if (propClause(values, f_assign, c))
				c.markDeleted();
			else {
				assert(c.size()); // cannot be empty at this point
				if (c.size() == 1) {
					const uint32 unit = *c;
					CHECKLIT(unit);
					if (UNASSIGNED(values[unit])) enqueueUnit(unit);
					else { learnEmpty(); return false; }
				}
			}
		}
		forall_occurs(ol, i) { 
			SCLAUSE& c = cnf[*i];
			assert(c.size());
			c.markDeleted();
		}
		ol.clear(true), f_ol.clear(true);
	}
	cleanProped();
	return true;
}
#endif

#if defined(BCP_TEST_UNITS) && defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
void Solver::injectTestUnits(const uint32& want)
{
	if (!want) return;
	inject_units_k<<<1, 1>>>(ot, vars->units, inf.maxVar, want);
	LASTERR("inject test units failed");
	SYNCALL;
	CHECK(cudaMemcpy(&vars->nUnits, vars->unitsSize, sizeof(uint32), cudaMemcpyDeviceToHost));
	LOG2(2, " [TEST] injected %u synthetic forced units before prop", vars->nUnits);
}
#endif

inline bool Solver::propClause(const LIT_ST* values, const uint32& lit, SCLAUSE& c)
{
	assert(c.size() > 1);
	uint32 sig = 0;
	uint32* j = c;
	forall_clause(c, i) {
		const uint32 other = *i;
		if (NEQUAL(other, lit)) {
			if (values[other] > 0) return true;
			*j++ = other;
			sig |= MAPHASH(other);
		}
	}
	assert(int(j - c) == c.size() - 1);
	assert(c.hasZero() < 0);
	c.set_sig(sig);
	c.pop();
	assert(c.isSorted());
	return false;
}

inline bool Solver::enqueueCached(const cudaStream_t& stream) {
	if (vars->nUnits) {
		nForced = sp->propagated;
		SYNC(stream); // sync units copy
		assert(vars->cachedUnits != NULL);
		LIT_ST* values = sp->value;
		uint32* t = vars->cachedUnits + vars->nUnits;
		for (uint32* u = vars->cachedUnits; u != t; u++) {
			const uint32 unit = *u;
			CHECKLIT(unit);
			const LIT_ST val = values[unit];
			if (UNASSIGNED(val)) enqueueDevUnit(unit);
			else if (!val) return false; // early conflict detection
		}
		if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
		else LOG2(2, " Propagating forced units..");
		SYNCALL; // sync ot creation
	}
	return true;
}

inline void	Solver::cleanProped() {
	if (vars->nUnits) {
		nForced = sp->propagated - nForced;
		if (nForced) {
			countAll();
			inf.numClauses = inf.numClausesSurvived, inf.numLiterals = inf.numLiteralsSurvived;
		}
		LOGREDALL(this, 2, "BCP Reductions");
		nForced = 0, vars->tmpUnits.clear();
		assert(vars->tmpUnits.data() == vars->unitsData);
		if (!opts.sub_en) reduceOTAsync();
		CHECK(cudaMemcpyAsync(vars->units, &vars->tmpUnits, sizeof(cuVecU), cudaMemcpyHostToDevice));
	}
}