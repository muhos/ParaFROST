/***********************************************************************[simplify.cu]
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
#include "timer.cuh"
#include "options.cuh"
#include "memory.cuh"

using namespace ParaFROST;

uint32 Solver::updateNumElims() 
{
	const uint32 remainedPVs = *vars->electedSize;
	assert(remainedPVs <= vars->numElected);
	inf.currDeletedVars = vars->numElected - remainedPVs;
	return remainedPVs;
}

inline void	Solver::updateNumPVs() 
{
	uint32 remainedPVs = updateNumElims();
	vars->currMelted += inf.currDeletedVars;
	LOG2(2, "  BVE eliminated %d variables while %d survived", inf.currDeletedVars, vars->numElected);
	vars->numElected = remainedPVs;
}

inline void	Solver::cleanManaged() 
{
	if (vars != NULL) delete vars;
	cumm.freeVars(), cumm.freeCNF(), cumm.freeOT();
	vars = NULL, cnf = NULL, ot = NULL;
}


inline void	Solver::initSimplifier()
{
	cleanManaged();
	cumm.freeFixed();
	cacher.destroy();
	simpstate = AWAKEN_SUCC;
	phase = multiplier = 0;
	ereCls = nForced = 0;
	dataoff = csoff = 0;
	flattened = compacted = false;
	assert(hcnf == NULL);
	if (stats.sigma.calls > 1) {
		size_t free = 0, tot = 0;
		CHECK(cudaMemGetInfo(&free, &tot));
		cumm.reset(free);
	}
}

void Solver::simplify(const bool& keep_gpu_mem)
{
	if (alldisabled()) return;
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	assert(stats.clauses.original);
	stats.sigma.calls++;
	simplifying(keep_gpu_mem);
	INCREASE_LIMIT(this, sigma, stats.sigma.calls, nlognlogn, true);
	last.sigma.reduces = stats.reduces + 1;
	if (opts.phases > 2) {
		opts.phases--;
		LOG2(2, "  simplify phases decreased to %d", opts.phases);
	}
	if (keep_gpu_mem)
		LOG2(2, " Keeping GPU memory for future calls");
	if (!opts.solve_en) killSolver();
	timer.start();
}

void Solver::awaken()
{
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	assert(isPropagated());
	initSimplifier();
	// overapproximate the size of new formula and saved witnesses
	size_t numCls = maxClauses(), numLits = maxLiterals();
	size_t savedLits = numCls + numLits;
	if (opts.phases) {
		size_t maxAddedCls = opts.ve_en ? stats.clauses.original : 0;
		size_t maxAddedLits = opts.ve_en ? size_t(stats.literals.original * opts.lits_mul) : 0;
		LOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
		numCls += maxAddedCls, numLits += maxAddedLits;
	}
	assert(inf.maxDualVars);
	if (!cumm.allocHist(cuhist, opts.proof_en) ||
		!cumm.allocAux(numCls) ||
		!cumm.allocVars(vars, savedLits) ||
		!cumm.allocPinned(vars, cuhist) ||
		!cumm.resizeCNF(cnf, numCls, numLits)) { simpstate = CNFALLOC_FAIL; return; }
	olist_cmp.init(cnf), compact_cmp.init(cnf);
	printStats(1, '-', CGREEN0);
	inf.numClauses = inf.numLiterals = 0;
	wt.clear(true);
	if (gopts.unified_access) {
		LOGN2(2, " Extracting clauses directly to device..");
		if (gopts.profile_gpu) cutimer.start();
		extractCNF(cnf, orgs), orgs.clear(true);
		extractCNF(cnf, learnts), learnts.clear(true);
		cm.destroy();
		assert(inf.numClauses == cnf->size());
		if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.io += cutimer.gpuTime();
	}
	else {
		LOGN2(2, " Extracting clauses heterogeneously to device..");
		if (gopts.profile_gpu) cutimer.start(streams[0]);
		cumm.createMirror(hcnf, maxClauses(), maxLiterals());
		extractCNF(hcnf, orgs), reflectCNF(streams[0], streams[1]), orgs.clear(true);
		extractCNF(hcnf, learnts), reflectCNF(streams[0], streams[1]), learnts.clear(true);
		cm.destroy();
		SYNC(streams[0]); SYNC(streams[1]);
		cumm.resizeCNFAsync(cnf, hcnf->data().size, hcnf->size());
		assert(hcnf->data().size == dataoff);
		assert(inf.numClauses == hcnf->size() && hcnf->size() == csoff);
		if (gopts.profile_gpu) cutimer.stop(streams[0]), stats.sigma.time.io += cutimer.gpuTime();
	}
	LOGENDING(2, 5, "(%d clauses extracted)", inf.numClauses);
    vars->varcore = gopts.hostKOpts.ve_fun_en ? vars->eligible : NULL;
	SYNC(0); // sync device CNF resize
	initDevVorg(cuhist); // load device pointers to constant memory
	if (gopts.profile_gpu) cutimer.start(streams[0]);
	prepareCNFAsync(cnf, streams[0]);
	if (gopts.profile_gpu) cutimer.stop(streams[0]), stats.sigma.time.sig += cutimer.gpuTime();
	if (opts.proof_en) 
		cumm.cuMemSetAsync(cuhist.d_lbyte, 0, inf.maxDualVars);
	assert(vorg.size() == inf.maxVar + 1);
	cuhist.fetchVars(vorg, streams[1]);
	if (opts.proof_en) {
		const uint32 *literals = flattenCNF(inf.numLiterals);
		uint32 proofcap = cuproof.count(literals, inf.numLiterals);
		proofcap *= 1.5;
		if (!cuproof.alloc(proofcap)) { simpstate = AWAKEN_FAIL; return; }
	}
}

void Solver::simplifying(const bool& keep_gpu_mem)
{
	/********************************/
	/*         awaken sigma         */
	/********************************/
	SLEEPING(sleep.sigma, opts.sigma_sleep_en);
	rootify();
	shrinkTop(false);
	if (orgs.empty()) {
		recycleWT();
		return;
	}
	timer.stop();
	stats.time.solve += timer.cpuTime();
	timer.start();
	awaken();
	if (simpstate == CNFALLOC_FAIL || simpstate == AWAKEN_FAIL) {
		recycle();
		return;
	}
	/********************************/
	/*      reduction phases        */
	/********************************/
	assert(!phase && !multiplier);
	const int64 imelted = inf.maxMelted, iclauses = inf.numClauses, iliterals = inf.numLiterals;
	int64 cdiff = INT64_MAX, ldiff = INT64_MAX;
	int64 clsbefore = inf.numClauses, litsbefore = inf.numLiterals;
	while (inf.numClauses && inf.numLiterals && !simpstate && !interrupted()) {
		if (!reallocOT(streams[0])) break;
		SYNC(streams[0]);
		reallocCNF();
		createOTAsync();
		if (!prop()) killSolver();
		else if (!inf.numClauses) break;
		if (!LCVE()) break;
		sortOT();
		if (stop(cdiff, ldiff)) { ERE(); break; }
		SUB(), VE(), BCE();
		cuproof.cacheProof(streams[4]);
		countAll();
		updateNumPVs();
		cacheNumUnits(streams[3]);
		inf.numClauses = inf.numClausesSurvived, inf.numLiterals = inf.numLiteralsSurvived;
		cdiff = clsbefore - inf.numClauses, clsbefore = inf.numClauses;
		ldiff = litsbefore - inf.numLiterals, litsbefore = inf.numLiterals;
		cacheUnits(streams[3]);
		phase++, multiplier++;
		multiplier += phase == opts.phases;
		cuproof.writeProof(streams[4]);
		flattened = false;
	}
	/********************************/
	/*          Write Back          */
	/********************************/
	inf.maxMelted += vars->currMelted;
	bool success = (iclauses != inf.numClauses);
	stats.sigma.all.clauses += iclauses - int64(inf.numClauses);
	stats.sigma.all.literals += iliterals - int64(inf.numLiterals);
	stats.sigma.all.variables += int64(inf.maxMelted) - imelted;
	last.shrink.removed = stats.shrunken;
	if (ereCls > inf.numClauses) stats.sigma.ere.removed += ereCls - inf.numClauses;
	if (inf.maxFrozen > sp->simplified) stats.units.forced += inf.maxFrozen - sp->simplified;
	if (!inf.unassigned || !inf.numClauses) { 
		LOG2(2, "  All clauses removed, thus CNF is satisfiable");
		stats.clauses.original = 0;
		stats.clauses.learnt = 0;
		stats.literals.original = 0;
		stats.literals.learnt = 0;
		cnfstate = SAT; 
		printStats(1, 's', CGREEN); 
		return;
	}
	if (vars->nUnits && reallocFailed()) {
		// TODO:
		// If reallocation failed while unit clauses are pending,
		// one can append them to the trail and continue solving.
		LOGERROR("unit clauses are pending and reallocation of CNF/OT failed");
	}
	assert(isPropagated());
	cacheEliminated(streams[5]);             	// if ERE is enabled, this transfer would be already done
	cumm.updateMaxCap(), cacher.updateMaxCap(); // for reporting GPU memory only
	markEliminated(streams[5]);              	// must be executed before map()
	if (keep_gpu_mem) {
		SYNCALL;
		assert(!simpstate);
		if (simpstate != CNFALLOC_FAIL && !compacted) 
			reallocCNF(true);
		cacheCNF(streams[0], streams[1]);
		printStats(1, 's', CGREEN); 
		return;
	}
	cacheCNF(streams[0], streams[1]);
	if (canMap()) map(true);
	else newBeginning();
	rebuildWT(opts.sigma_priorbins);
	if (BCP()) {
		LOG2(1, " Propagation after simplify proved a contradiction");
		learnEmpty();
	}
	UPDATE_SLEEPER(this, sigma, success);
	printStats(1, 's', CGREEN);
	timer.stop(), stats.time.simp += timer.cpuTime();
}

void Solver::optSimp()
{
	LOGN2(2, " Initializing device options..");
	gopts.init(opts.proof_en && !opts.proof_nonbinary_en);
    assert(SH_MAX_BVE_OUT1 >= gopts.hostKOpts.xor_max_arity);
    assert(SH_MAX_ERE_OUT >= gopts.hostKOpts.ere_clause_max);
	initDevOpts();
	LOGDONE(2, 5);
	initSharedMem();
}

void Solver::freeSimp()
{
	SYNCALL;
	cacher.destroy();
	cumm.freeFixed();
	cumm.freePinned();
	cleanManaged();
	destroyStreams();
	cumm.breakMirror(), hcnf = NULL;
	if (opts.proof_en) 
		cuproof.destroy();
}