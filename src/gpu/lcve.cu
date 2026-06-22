/***********************************************************************[lcve.cu]
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
#include "key.cuh"
#include "grid.cuh"
#include "timer.cuh"
#include "table.cuh"
#include "cache.cuh"
#include "options.cuh"
#include "variables.cuh"
#include "histogram.cuh"
#include "vstate.hpp"
#include <thrust/sort.h>

using namespace ParaFROST;

#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)

__device__
bool depFreeze_d(
	const OL& ol, const CNF& cnf,
	VSTATE* __restrict__ vstate,
	uint32* __restrict__ frozenList,
	uint32& frozenTail,
	const uint32 cand,
	const uint32 pmax, const uint32 nmax,
	const int maxcsize)
{
	const uint32 savedTail = frozenTail;
	for (const S_REF* i = ol.data(), *end = ol.end(); i != end; i++) {
		const SCLAUSE& c = cnf[*i];
		if (c.deleted()) continue;
		if ((int)c.size() > maxcsize) {
			for (uint32 k = savedTail; k < frozenTail; k++)
				vstate[frozenList[k]].frozen = 0;
			frozenTail = savedTail;
			return false;
		}
		for (const uint32* k = c.data(), *cend = c.end(); k != cend; k++) {
			const uint32 v = ABS(*k);
			if (!vstate[v].frozen && v != cand) {
				vstate[v].frozen = 1;
				frozenList[frozenTail++] = v;
			}
		}
	}
	return true;
}

__global__
void lcve_k(
	const uint32* __restrict__ eligible,
	const uint32  maxVar,
	cuVecU*       elected,
	uint32*       scores,
	VSTATE* __restrict__     d_vstate,
	const uint32* __restrict__ d_hist,
	const OT* __restrict__     ot,
	const CNF* __restrict__    cnf,
	const uint32  pmax,
	const uint32  nmax,
	const uint32  maxoccurs,
	const int     maxcsize,
	const bool* __restrict__   d_assumed)
{
	uint32 frozenTail = 0;
	uint32* frozenList = scores + 1; // scores[0] = count, scores[1..] = frozen vars

	for (uint32 ei = 0; ei < maxVar; ei++) {
		const uint32 cand = eligible[ei];
		assert(cand >= 1 && cand <= maxVar);
		if (d_vstate[cand].frozen) continue;
		if (d_vstate[cand].state) continue;
		if (d_assumed && d_assumed[cand]) continue;
		const uint32 p = V2L(cand), n = NEG(p);
		const uint32 ps = d_hist[p], ns = d_hist[n];
		assert((*ot)[p].size() >= ps);
		assert((*ot)[n].size() >= ns);
		if (!ps && !ns) continue;
		if (ps > maxoccurs || ns > maxoccurs) break;
		if (ps >= pmax && ns >= nmax) break;
		if (depFreeze_d((*ot)[p], *cnf, d_vstate, frozenList, frozenTail, cand, pmax, nmax, maxcsize) &&
			depFreeze_d((*ot)[n], *cnf, d_vstate, frozenList, frozenTail, cand, pmax, nmax, maxcsize))
			elected->_push(cand);
	}
	scores[0] = frozenTail;
}

#define MIS_NONE      0u
#define MIS_UNDECIDED 1u
#define MIS_ELECTED   2u
#define MIS_FROZEN    3u

_PFROST_IN_D_
bool mis_oversize(const OL& ol, const CNF& cnf, const int& maxcsize)
{
	for (const S_REF* i = ol.data(), *end = ol.end(); i != end; i++) {
		const SCLAUSE& c = cnf[*i];
		if (!c.deleted() && (int)c.size() > maxcsize) return true;
	}
	return false;
}

_PFROST_IN_D_
bool mis_blocked(const OL& ol, const CNF& cnf, const uint32& cand, const uint32& r,
	const uint32* __restrict__ rank, const uint32* __restrict__ state)
{
	for (const S_REF* i = ol.data(), *end = ol.end(); i != end; i++) {
		const SCLAUSE& c = cnf[*i];
		if (c.deleted()) continue;
		for (const uint32* k = c.data(), *cend = c.end(); k != cend; k++) {
			const uint32 u = ABS(*k);
			if (u != cand && state[u] == MIS_UNDECIDED && rank[u] < r) return true;
		}
	}
	return false;
}

_PFROST_IN_D_
void mis_freeze(const OL& ol, const CNF& cnf, const uint32& cand,
	VSTATE* __restrict__ vstate, uint32* __restrict__ state)
{
	for (const S_REF* i = ol.data(), *end = ol.end(); i != end; i++) {
		const SCLAUSE& c = cnf[*i];
		if (c.deleted()) continue;
		for (const uint32* k = c.data(), *cend = c.end(); k != cend; k++) {
			const uint32 u = ABS(*k);
			if (u != cand) {
				vstate[u].frozen = 1;
				if (state[u] == MIS_UNDECIDED) state[u] = MIS_FROZEN;
			}
		}
	}
}

__global__
void mis_init_k(
	const uint32* __restrict__ eligible, const uint32 maxVar,
	uint32* __restrict__ rank, uint32* __restrict__ state,
	const VSTATE* __restrict__ vstate, const uint32* __restrict__ d_hist,
	const OT* __restrict__ ot, const CNF* __restrict__ cnf,
	const uint32 pmax, const uint32 nmax, const uint32 maxoccurs, const int maxcsize,
	const bool* __restrict__ d_assumed)
{
	for_parallel_x(ei, maxVar) {
		const uint32 v = eligible[ei];
		rank[v] = ei;
		uint32 st = MIS_NONE;
		if (!vstate[v].state && !(d_assumed && d_assumed[v])) {
			const uint32 p = V2L(v), n = NEG(p);
			const uint32 ps = d_hist[p], ns = d_hist[n];
			if ((ps || ns) && ps <= maxoccurs && ns <= maxoccurs && !(ps >= pmax && ns >= nmax)
				&& !mis_oversize((*ot)[p], *cnf, maxcsize)
				&& !mis_oversize((*ot)[n], *cnf, maxcsize))
				st = MIS_UNDECIDED;
		}
		state[v] = st;
	}
}

__global__
void mis_round_k(
	const uint32* __restrict__ eligible, const uint32 maxVar,
	const uint32* __restrict__ rank, const uint32* __restrict__ state,
	uint32* __restrict__ win, uint32* __restrict__ remaining,
	const OT* __restrict__ ot, const CNF* __restrict__ cnf)
{
	for_parallel_x(ei, maxVar) {
		const uint32 v = eligible[ei];
		if (state[v] == MIS_UNDECIDED) {
			*remaining = 1;
			const uint32 p = V2L(v), n = NEG(p);
			win[v] = !(mis_blocked((*ot)[p], *cnf, v, ei, rank, state) ||
			           mis_blocked((*ot)[n], *cnf, v, ei, rank, state));
		}
	}
}

__global__
void mis_freeze_k(
	const uint32* __restrict__ eligible, const uint32 maxVar,
	const uint32* __restrict__ win, uint32* __restrict__ state,
	VSTATE* __restrict__ vstate, cuVecU* __restrict__ elected,
	const OT* __restrict__ ot, const CNF* __restrict__ cnf)
{
	for_parallel_x(ei, maxVar) {
		const uint32 v = eligible[ei];
		if (state[v] == MIS_UNDECIDED && win[v]) {
			state[v] = MIS_ELECTED;
			elected->insertAggr(v);
			mis_freeze((*ot)[V2L(v)], *cnf, v, vstate, state);
			mis_freeze((*ot)[NEG(V2L(v))], *cnf, v, vstate, state);
		}
	}
}

__global__
void mis_collect_k(const uint32 maxVar, const VSTATE* __restrict__ vstate, uint32* __restrict__ scores)
{
	for_parallel_x(i, maxVar) {
		const uint32 v = i + 1;
		if (vstate[v].frozen)
			scores[atomicAggInc(&scores[0]) + 1] = v;
	}
}

#endif

__global__
void assign_scores(
	uint32* __restrict__ eligible,
	uint32* __restrict__ scores,
	const uint32* __restrict__ hist,
	uint32 size)
{
	for_parallel_x(tid, size) {
		const uint32 v = tid + 1;
		const uint32 p = V2L(v), ps = hist[p], ns = hist[NEG(p)];
		eligible[tid] = v;
		scores[v] = ps * ns;
	}
}

__global__ 
void assign_scores(
	uint32* __restrict__ eligible,
	uint32* __restrict__ scores,
	uint32* __restrict__ hist,
	const OT* __restrict__ ot,
	uint32 size)
{
	for_parallel_x(tid, size) {
		const uint32 v = tid + 1;
		const uint32 p = V2L(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
		hist[p] = ps, hist[n] = ns;
		eligible[tid] = v;
		scores[v] = ps * ns;
	}
}

__global__ 
void mapfrozen_k(const uint32* __restrict__ frozen, uint32* __restrict__ varcore, const uint32 size)
{
	for_parallel_x(tid, size) {
		assert(frozen[tid] && frozen[tid] < NOVAR);
		varcore[frozen[tid]] = tid;
	}
}

void calcScores(VARS* vars, uint32* hist)
{
	grid_t nThreads = BLOCK1D;
	OPTIMIZEBLOCKS(inf.maxVar, nThreads, 0);
	assign_scores << <nBlocks, nThreads >> > (vars->eligible, vars->scores, hist, inf.maxVar);
	LASTERR("Assigning scores failed");
	SYNCALL;
}

void calcScores(VARS* vars, uint32* hist, OT* ot)
{
	grid_t nThreads = BLOCK1D;
	OPTIMIZEBLOCKS(inf.maxVar, nThreads, 0);
	assign_scores << <nBlocks, nThreads >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
	LASTERR("Assigning scores failed");
	SYNCALL;
}

void mapFrozenAsync(VARS* vars, const uint32& size)
{
	assert(vars->varcore == vars->eligible); // an alies of eligible
	grid_t nThreads = BLOCK1D;
	OPTIMIZEBLOCKS(size, nThreads, 0);
	// 'vars->scores' is an alies for frozen vars on the GPU side
	mapfrozen_k << <nBlocks, nThreads >> > (vars->scores, vars->varcore, size);
	if (gopts.sync_always) {
		LASTERR("Mapping frozen failed");
		SYNCALL;
	}
}

void Solver::varReorder()
{
	LOGN2(2, " Finding eligible variables for LCVE..");
	assert(cuhist.d_hist != NULL);
	// NOTE: OT creation will be synced in calcScores call
	if (gopts.profile_gpu) cutimer.start();
	if (vars->nUnits) 
		calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
	else 
		calcScores(vars, cuhist.d_hist);
	if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.vo += cutimer.gpuTime();
	cuhist.cacheHist(inf.maxDualVars, streams[2]);
	if (gopts.profile_gpu) cutimer.start(streams[3]);
	cacher.insert(cumm.scatter(), cumm.scatterCap());
	thrust::sort(thrust::cuda::par(tca).on(streams[3]), vars->eligible, vars->eligible + inf.maxVar, GPU_LCV_CMP(vars->scores));
	cacher.erase(cumm.scatterCap());
	LOGDONE(2, 5);
	vars->nUnits = 0;
	SYNC(streams[2]);
	if (gopts.profile_gpu) cutimer.stop(streams[3]), stats.sigma.time.vo += cutimer.gpuTime();
	if (verbose == 4 && !cumm.isDeviceCNF()) {
		LOG0(" Eligible variables:");
		for (uint32 v = 0; v < inf.maxVar; v++) {
			uint32 x = vars->eligible[v], p = V2L(x), n = NEG(p);
			LOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, cuhist[p], cuhist[n], vars->scores[x]);
		}
	}
}

bool Solver::LCVE()
{
	// reorder variables
	varReorder();

	// extended LCVE
	LOG2(2, " Electing variables (p-mu: %d, n-mu: %d)..", opts.mu_pos << multiplier, opts.mu_neg << multiplier);
	const uint32 maxoccurs = opts.lcve_max_occurs;
	const uint32 pmax = opts.mu_pos << multiplier;
	const uint32 nmax = opts.mu_neg << multiplier;

	vars->numElected = 0;
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	cuPool assumedPool;
	bool* d_assumed = NULL;
	const bool hasAssumptions = incremental && assumptions.size();
	size_t assumedBytes = 0;
	if (hasAssumptions) {
		assumedBytes = assumed.size() * sizeof(bool);
		assert(assumed.size() > inf.maxVar);
		if (!cumm.allocDynamic(assumedPool, assumedBytes, "Assumptions")) { vars->numElected = 0; return false; }
		d_assumed = (bool*)assumedPool.mem;
	}
	cuPool mis;
	const bool useMIS = opts.lcve_fast;
	if (useMIS) {
		const uint32 nvars = inf.maxVar + 1;
		const size_t words = size_t(nvars) * 3 + 1; // rank + state + win + remaining
		if (!cumm.allocDynamic(mis, words * sizeof(uint32), "MIS")) {
			if (hasAssumptions) cumm.freeDynamic(assumedPool);
			vars->numElected = 0;
			return false;
		}
	}
	CHECK(cudaMemsetAsync(vars->electedSize, 0, sizeof(uint32), streams[0])); // GPU-side clear of elected size
	sp->stacktail = sp->tmpstack; // reset frozen list head
	// copy vstate to device (frozen bits will be set by the kernel via the frozen field)
	CHECK(cudaMemcpyAsync(cumm.deviceVstate(), sp->vstate, (inf.maxVar + 1) * sizeof(VSTATE), cudaMemcpyHostToDevice, streams[1]));
	CHECK(cudaMemsetAsync(vars->scores, 0, sizeof(uint32), streams[2])); // frozen count = scores[0]
	if (hasAssumptions)
		CHECK(cudaMemcpyAsync(d_assumed, assumed.data(), assumedBytes, cudaMemcpyHostToDevice, streams[0]));
	SYNC(streams[2]); SYNC(streams[1]); SYNC(streams[0]);
	// parallel Luby MIS election (nondeterministic, valid maximal independent set)
	if (useMIS) {
		const uint32 nvars = inf.maxVar + 1;
		uint32* rank  = (uint32*)mis.mem;
		uint32* state = rank + nvars;
		uint32* win   = state + nvars;
		uint32* d_rem = win + nvars;
		grid_t nThreads = BLOCK1D;
		OPTIMIZEBLOCKS(inf.maxVar, nThreads, 0);
		mis_init_k<<<nBlocks, nThreads>>>(vars->eligible, inf.maxVar, rank, state,
			cumm.deviceVstate(), cuhist.d_hist, ot, cnf, pmax, nmax, maxoccurs, opts.lcve_clause_max,
			d_assumed);
		LASTERR("MIS init failed");
		uint32 rem = 0;
		do {
			CHECK(cudaMemset(d_rem, 0, sizeof(uint32)));
			mis_round_k<<<nBlocks, nThreads>>>(vars->eligible, inf.maxVar, rank, state, win, d_rem, ot, cnf);
			mis_freeze_k<<<nBlocks, nThreads>>>(vars->eligible, inf.maxVar, win, state, cumm.deviceVstate(), vars->elected, ot, cnf);
			LASTERR("MIS round failed");
			CHECK(cudaMemcpy(&rem, d_rem, sizeof(uint32), cudaMemcpyDeviceToHost));
		} while (rem);
		mis_collect_k<<<nBlocks, nThreads>>>(inf.maxVar, cumm.deviceVstate(), vars->scores);
		LASTERR("MIS collect failed");
		SYNCALL;
		cumm.freeDynamic(mis);
	}
	else {
		lcve_k<<<1, 1>>>(
			vars->eligible, inf.maxVar, vars->elected, vars->scores,
			cumm.deviceVstate(), cuhist.d_hist, ot, cnf,
			pmax, nmax, maxoccurs, opts.lcve_clause_max,
			d_assumed);
		LASTERR("LCVE kernel failed");
		SYNCALL;
	}
	if (hasAssumptions) cumm.freeDynamic(assumedPool);
	CHECK(cudaMemcpy(&vars->numElected, vars->electedSize, sizeof(uint32), cudaMemcpyDeviceToHost));
#else
	vars->elected->clear();
	sp->stacktail = sp->tmpstack;
	uint32* evars = vars->eligible;
	uint32* eend = evars + inf.maxVar;
	uint32*& tail = sp->stacktail;
	LIT_ST* frozen = sp->frozen;
	const VSTATE* states = sp->vstate;
	OT& ot = *this->ot; // cache 'ot' reference on host
	while (evars != eend) {
		const uint32 cand = *evars++;
		CHECKVAR(cand);
		if (frozen[cand]) continue;
		if (states[cand].state) continue;
		if (iassumed(cand)) continue;
		const uint32 p = V2L(cand), n = NEG(p);
		const uint32 ps = cuhist[p], ns = cuhist[n];
		assert(ot[p].size() >= ps);
		assert(ot[n].size() >= ns);
		if (!ps && !ns) continue;
		if (ps > maxoccurs || ns > maxoccurs) break;
		if (ps >= pmax && ns >= nmax) break;
		if (depFreeze(ot[p], frozen, tail, cand, pmax, nmax) &&
			depFreeze(ot[n], frozen, tail, cand, pmax, nmax))
			vars->elected->_push(cand);
	}
	vars->numElected = vars->elected->size();
	assert(verifyLCVE());
#endif

	if (vars->numElected) {
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
		uint32 mcv = 0;
		CHECK(cudaMemcpy(&mcv, vars->electedData + vars->numElected - 1, sizeof(uint32), cudaMemcpyDeviceToHost));
#else
		const uint32 mcv = vars->elected->back();
#endif
		const uint32 pmcv = V2L(mcv);
		LOG2(2, " Elected %d variables, with mcv: %d, pH: %d, nH: %d.", vars->numElected, mcv, cuhist[pmcv], cuhist[NEG(pmcv)]);
	}

	if (verbose > 3 && !cumm.isDeviceCNF()) {
		LOGN0(" PLCVs ");
		printVars(*vars->elected, vars->numElected, 'v');
	}

	mapFrozen(); // async. call
	clearFrozen();

	if (vars->numElected < opts.lcve_min_vars) {
		if (gopts.hostKOpts.ve_fun_en) SYNC(0);
		if (verbose > 1) LOGWARNING("parallel variables not enough -> skip GPU simplifier");
		return false;
	}
	return true;
}

inline void	Solver::mapFrozen()
{
	if (!gopts.hostKOpts.ve_fun_en) return;
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	uint32 nFrozen = 0;
	CHECK(cudaMemcpy(&nFrozen, vars->scores, sizeof(uint32), cudaMemcpyDeviceToHost));
	if (!nFrozen) { vars->varcore = NULL; sp->stacktail = sp->tmpstack; return; }
	assert(nFrozen <= inf.maxVar);
	CHECK(cudaMemcpy(sp->tmpstack, vars->scores + 1, nFrozen * sizeof(uint32), cudaMemcpyDeviceToHost));
	sp->stacktail = sp->tmpstack + nFrozen;
	if (gopts.profile_gpu) cutimer.start();
	assert(vars->varcore == vars->eligible);
	grid_t nThreads = BLOCK1D;
	OPTIMIZEBLOCKS(nFrozen, nThreads, 0);
	mapfrozen_k<<<nBlocks, nThreads>>>(vars->scores + 1, vars->varcore, nFrozen);
	if (gopts.sync_always) {
		LASTERR("Mapping frozen failed");
		SYNCALL;
	}
	if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.ve += cutimer.gpuTime();
#else
	const uint32 *frozen = sp->tmpstack;
	const uint32 *end = sp->stacktail;
	const uint32 nFrozen = uint32(end - frozen);
	if (!nFrozen) { vars->varcore = NULL; return; }
	assert(nFrozen <= inf.maxVar);
	CHECK(cudaMemcpy(vars->scores, frozen, nFrozen * sizeof(uint32), cudaMemcpyHostToDevice));
	if (gopts.profile_gpu) cutimer.start();
	mapFrozenAsync(vars, nFrozen);
	if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.ve += cutimer.gpuTime();
#endif
}

inline bool Solver::depFreeze(OL& ol, LIT_ST* frozen, uint32*& tail, const uint32& cand, const uint32& pmax, const uint32& nmax)
{
	const int maxcsize = opts.lcve_clause_max;
	uint32* first = tail;
	CNF& cnf = *this->cnf; // cache 'cnf' reference on host
	forall_occurs(ol, i) {
		SCLAUSE& c = cnf[*i];
		if (c.deleted()) continue;
		if (c.size() > maxcsize) {
			uint32* from = first;
			while (from != tail)
				frozen[*from++] = 0;
			tail = first;
			return false;
		}
		forall_clause(c, k) {
			const uint32 v = ABS(*k);
			CHECKVAR(v);
			if (!frozen[v] && NEQUAL(v, cand)) {
				frozen[v] = 1;
				assert(tail < sp->tmpstack + inf.maxVar);
				*tail++ = v;
			}
		}
	}
	return true;
}

inline bool	Solver::verifyLCVE() {
	for (uint32 v = 0; v < vars->numElected; v++)
		if (sp->frozen[vars->elected->at(v)])
			return false;
	return true;
}
