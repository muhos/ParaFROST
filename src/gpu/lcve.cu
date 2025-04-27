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
#include <thrust/sort.h>

using namespace ParaFROST;

__global__ 
void assign_scores(
	uint32* __restrict__ eligible,
	uint32* __restrict__ scores,
	const uint32* __restrict__ hist,
	uint32 size)
{
	grid_t tid = global_tx;
	while (tid < size) {
		const uint32 v = tid + 1;
		const uint32 p = V2L(v), ps = hist[p], ns = hist[NEG(p)];
		eligible[tid] = v;
		scores[v] = ps * ns;
		tid += stride_x;
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
	grid_t tid = global_tx;
	while (tid < size) {
		const uint32 v = tid + 1;
		const uint32 p = V2L(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
		hist[p] = ps, hist[n] = ns;
		eligible[tid] = v;
		scores[v] = ps * ns;
		tid += stride_x;
	}
}

__global__ 
void mapfrozen_k(const uint32* __restrict__ frozen, uint32* __restrict__ varcore, const uint32 size)
{
	grid_t tid = global_tx;
	while (tid < size) {
		assert(frozen[tid] && frozen[tid] < NOVAR);
		varcore[frozen[tid]] = tid;
		tid += stride_x;
	}
}

void calcScores(VARS* vars, uint32* hist)
{
	if (gopts.profile_gpu) cutimer->start();
	OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
	assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, inf.maxVar);
	if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
	LASTERR("Assigning scores failed");
	SYNCALL;
}

void calcScores(VARS* vars, uint32* hist, OT* ot)
{
	if (gopts.profile_gpu) cutimer->start();
	OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
	assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
	if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
	LASTERR("Assigning scores failed");
	SYNCALL;
}

void mapFrozenAsync(VARS* vars, const uint32& size)
{
	assert(vars->varcore == vars->eligible); // an alies of eligible
	if (gopts.profile_gpu) cutimer->start();
	OPTIMIZEBLOCKS(size, BLOCK1D);
	// 'vars->scores' is an alies for frozen vars on the GPU side
	mapfrozen_k << <nBlocks, BLOCK1D >> > (vars->scores, vars->varcore, size);
	if (gopts.sync_always) {
		LASTERR("Mapping frozen failed");
		SYNCALL;
	}
	if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
}

void Solver::varReorder()
{
	LOGN2(2, " Finding eligible variables for LCVE..");
	assert(cuhist.d_hist != NULL);
	// NOTE: OT creation will be synced in calcScores call
	if (vars->nUnits) calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
	else calcScores(vars, cuhist.d_hist);
	cuhist.cacheHist(inf.nDualVars, streams[2]);
	if (gopts.profile_gpu) cutimer->start(streams[3]);
	cacher.insert(cumm.scatter(), cumm.scatterCap());
	thrust::sort(thrust::cuda::par(tca).on(streams[3]), vars->eligible, vars->eligible + inf.maxVar, GPU_LCV_CMP(vars->scores));
	cacher.erase(cumm.scatterCap());
	LOGDONE(2, 5);
	vars->nUnits = 0;
	SYNC(streams[2]);
	if (gopts.profile_gpu) cutimer->stop(streams[3]), cutimer->vo += cutimer->gpuTime();
	if (verbose == 4) {
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
	uint32* evars = vars->eligible;
	uint32* eend = evars + inf.maxVar;
	LOGN2(2, " Finding eligible variables for LCVE..");
	assert(cuhist.d_hist != NULL);
	// NOTE: OT creation will be synced in calcScores call
	if (vars->nUnits) calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
	else calcScores(vars, cuhist.d_hist);
	cuhist.cacheHist(inf.nDualVars, streams[2]);
	if (gopts.profile_gpu) cutimer->start(streams[3]);
	cacher.insert(cumm.scatter(), cumm.scatterCap());
	thrust::sort(thrust::cuda::par(tca).on(streams[3]), evars, eend, GPU_LCV_CMP(vars->scores));
	cacher.erase(cumm.scatterCap());
	LOGDONE(2, 5);
	vars->nUnits = 0;
	SYNC(streams[2]);
	if (gopts.profile_gpu) cutimer->stop(streams[3]), cutimer->vo += cutimer->gpuTime();
	if (verbose == 4) {
		LOG0(" Eligible variables:");
		for (uint32 v = 0; v < inf.maxVar; v++) {
			uint32 x = evars[v], p = V2L(x), n = NEG(p);
			LOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, cuhist[p], cuhist[n], vars->scores[x]);
		}
	}

	// extended LCVE
	LOGN2(2, " Electing variables (p-mu: %d, n-mu: %d)..", opts.mu_pos << multiplier, opts.mu_neg << multiplier);
	sp->stacktail = sp->tmpstack;
	uint32*& tail = sp->stacktail;
	LIT_ST* frozen = sp->frozen;
	const uint32 maxoccurs = opts.lcve_max_occurs;
	const uint32 pmax = opts.mu_pos << multiplier;
	const uint32 nmax = opts.mu_neg << multiplier;
	const VSTATE* states = sp->vstate;
	OT& ot = *this->ot; // cache 'ot' reference on host
	vars->numElected = 0;
	vars->elected->clear();
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

	if (vars->numElected) {
		const uint32 mcv = vars->elected->back(), pmcv = V2L(mcv);
		LOGENDING(2, 5, "(%d elected, mcv: %d, pH: %d, nH: %d)", vars->numElected, mcv, cuhist[pmcv], cuhist[NEG(pmcv)]);
	}
	else LOGDONE(2, 5);

	if (verbose > 3) {
		LOGN0(" PLCVs ");
		printVars(*vars->elected, vars->numElected, 'v');
	}

	mapFrozen(); // async. call
	clearFrozen();

	if (vars->numElected < opts.lcve_min_vars) {
		if (gopts.hostKOpts.ve_fun_en) SYNC(0); 
		if (!vars->numElected) LOGDONE(2, 5);
		if (verbose > 1) LOGWARNING("parallel variables not enough -> skip GPU simplifier");
		return false;
	}
	return true;
}

inline void	Solver::mapFrozen()
{
    if (!gopts.hostKOpts.ve_fun_en) return;
	const uint32 *frozen = sp->tmpstack;
	const uint32 *end = sp->stacktail;
	const uint32 nFrozen = uint32(end - frozen);
	if (!nFrozen) { vars->varcore = NULL; return; }
	assert(nFrozen <= inf.maxVar);
	CHECK(cudaMemcpy(vars->scores, frozen, nFrozen * sizeof(uint32), cudaMemcpyHostToDevice));
	mapFrozenAsync(vars, nFrozen);
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

