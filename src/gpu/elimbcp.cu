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

#define BCP_UNSET 0u
#define BCP_TRUE  1u
#define BCP_FALSE 2u

#define BCP_CURR  0 // Current frontier size
#define BCP_NEXT  1 // Next frontier size
#define BCP_CONFL 2 // Conflict flag
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
			if (c.deleted()) continue;
			uint32 unit = 0; int nunset = 0; bool sat = false;
			forall_clause_const(c, k) {
				const uint32 ve = bcp_lit_val(state, *k);
				if (ve == BCP_TRUE) { sat = true; break; }
				if (ve == BCP_UNSET) { unit = *k; if (++nunset > 1) break; }
			}
			if (sat) continue;
			if (!nunset) { ctrl[BCP_CONFL] = 1; continue; } // all falsified -> conflict
			if (nunset == 1) {
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
		if (c.deleted() || (!c.original() && !c.learnt())) continue;
		uint32* dst = c.data();
		uint32 sig = 0; int newsz = 0; bool sat = false;
		forall_clause(c, k) {
			const uint32 lit = *k;
			const uint32 ve = bcp_lit_val(state, lit);
			if (ve == BCP_TRUE) { sat = true; break; }
			if (ve == BCP_FALSE) continue; // drop falsified literal
			*dst++ = lit, sig |= MAPHASH(lit), newsz++; // keep unset literal
		}
		if (sat) c.markDeleted();
		else { c.set_sig(sig); c.resize(newsz); }
	}
}

#ifdef BCP_DEBUG

#ifndef BCP_DEBUG_MAX
#define BCP_DEBUG_MAX 64 // cap printed clauses per literal
#endif

__device__
void bcp_print_clause(const SCLAUSE& c, const S_REF r)
{
	printf("c [BCP-DBG]    ref %llu [%s sz=%d]:", (unsigned long long)r, c.deleted() ? "DEL" : "ok", c.size());
	for (int k = 0; k < c.size(); k++) {
		const uint32 l = c[k];
		printf(" %d", SIGN(l) ? -(int)ABS(l) : (int)ABS(l));
	}
	printf("\n");
}
__device__
void bcp_print_list(const CNF& cnf, const OL& list, const int signedv)
{
	printf("c [BCP-DBG]  lit %d occurs in %u clause(s)%s:\n", signedv, list.size(),
		list.size() > BCP_DEBUG_MAX ? " (truncated)" : "");
	const uint32 lim = list.size() < BCP_DEBUG_MAX ? list.size() : BCP_DEBUG_MAX;
	for (uint32 j = 0; j < lim; j++) bcp_print_clause(cnf[list[j]], list[j]);
}
__global__
void print_unit_clauses_k(
	const CNF* __restrict__ cnf, const OT* __restrict__ ot,
	const cuVecU* __restrict__ units, const uint32* __restrict__ unitsData,
	const uint32 count, const int when)
{
	if (blockIdx.x || threadIdx.x) return;
	printf("c [BCP-DBG] ===== %s propagation: %u unit(s) =====\n", when ? "AFTER" : "BEFORE", count);
	for (uint32 i = 0; i < count; i++) {
		const uint32 u = (*units)[i];
		const uint32 v = ABS(u);
		printf("c [BCP-DBG] unit %d (var %u):\n", SIGN(u) ? -(int)v : (int)v, v);
		const uint32 p = V2L(v);
		bcp_print_list(*cnf, (*ot)[p], (int)v);
		bcp_print_list(*cnf, (*ot)[NEG(p)], -(int)v);
	}
}
#endif

#endif // defined(USE_CUARENA) && defined(USE_DEVICE_CNF)

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
	if (!vars->nUnits) return true;
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
#ifdef BCP_DEBUG
	print_unit_clauses_k<<<1, 1>>>(cnf, ot, vars->units, vars->unitsData, vars->nUnits, 0);
	LASTERR("BCP debug (before) failed");
	SYNCALL;
#endif
	// Seed the initial units into the state and frontier
	{
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(vars->nUnits, nThreads, 0);
		bcp_seed_k<<<nBlocks, nThreads>>>(d_state, front_a, d_ctrl, vars->eliminated, vars->unitsData, vars->nUnits);
		LASTERR("BCP seed failed");
	}
	// Propagate the units until no more are found or a conflict is detected
	const dim3 blk2(32, 8); 
	const dim3 grid2(64, 64);
	uint32 ctrl_host[3] = { 0, 0, 0 };
	for (int batch = 1; ; batch = MIN(batch << 1, 64)) {
		for (int b = 0; b < batch; b++) {
			bcp_propagate_k<<<grid2, blk2>>>(cnf, ot, d_state, vars->eliminated, front_a, front_b, d_ctrl, vars->units);
			bcp_advance_k<<<1, 1>>>(d_ctrl);
		}
		LASTERR("BCP propagate failed");
		CHECK(cudaMemcpy(ctrl_host, d_ctrl, 3 * sizeof(uint32), cudaMemcpyDeviceToHost));
		if (ctrl_host[BCP_CONFL] || !ctrl_host[BCP_CURR]) break;
	}
	if (ctrl_host[BCP_CONFL]) { cumm.freeDynamic(bcp); learnEmpty(); return false; }
	// Delete satisfied, strengthen others
	{
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.numClauses, nThreads, 0);
		bcp_apply_k<<<nBlocks, nThreads>>>(cnf, d_state);
		LASTERR("BCP apply failed");
		SYNCALL;
	}
#ifdef BCP_DEBUG
	print_unit_clauses_k<<<1, 1>>>(cnf, ot, vars->units, vars->unitsData, vars->nUnits, 1);
	LASTERR("BCP debug (after) failed");
	SYNCALL;
#endif
	// Update solver state.
	uint32 ntrail = 0;
	CHECK(cudaMemcpy(&ntrail, vars->unitsSize, sizeof(uint32), cudaMemcpyDeviceToHost));
	assert(ntrail <= inf.maxVar);
	if (ntrail) {
		CHECK(cudaMemcpy(sp->tmpstack, vars->unitsData, ntrail * sizeof(uint32), cudaMemcpyDeviceToHost));
		for (uint32 i = 0; i < ntrail; i++) {
			const uint32 unit = sp->tmpstack[i];
			CHECKLIT(unit);
			if (i < vars->nUnits) enqueueDevUnit(unit);
			else                  enqueueUnit(unit);
		}
		sp->propagated = trail.size(); // device already propagated these
	}
	cumm.freeDynamic(bcp);
	stats.units.forced += ntrail;
	LOGREDALL(this, 2, "BCP Reductions");
	countAll();
	inf.numClauses = inf.numClausesSurvived, inf.numLiterals = inf.numLiteralsSurvived;
	createOTAsync(); // rebuild a consistent OT from the simplified CNF
	CHECK(cudaMemset(vars->unitsSize, 0, sizeof(uint32))); // reset device units header
	vars->nUnits = 0;
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