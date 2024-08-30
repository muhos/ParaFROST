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

using namespace ParaFROST;

inline void	Solver::cleanManaged() 
{
	if (vars != NULL) delete vars;
	cumm.freeVars(), cumm.freeCNF(), cumm.freeOT();
	vars = NULL, cnf = NULL, ot = NULL;
}

inline bool	Solver::reallocCNF()
{
	int times = phase + 1;
	if (times > 1 && times != opts.phases && (times % opts.shrink_rate) == 0) {
		size_t maxAddedCls = opts.ve_en ? inf.nClauses : 0;
		size_t maxAddedLits = opts.ve_en ? size_t(stats.literals.original * opts.lits_mul) : 0;
		LOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
		if (!cumm.resizeCNF(cnf, inf.nClauses + maxAddedCls, inf.nLiterals + maxAddedLits)) {
			simpstate = CNFALLOC_FAIL, compacted = false;
			return false;
		}
		compacted = true;
		olist_cmp.init(cnf), compact_cmp.init(cnf);
	}
	else cumm.cacheCNFPtr(cnf), compacted = false;
	return true;
}

inline bool	Solver::reallocOT(const cudaStream_t& stream)
{
	assert(inf.nLiterals);
	if (!flattenCNF(inf.nLiterals)) { simpstate = OTALLOC_FAIL; return false; }
	histSimp(inf.nLiterals);
	if (!cumm.resizeOTAsync(ot, inf.nLiterals, stream)) { simpstate = OTALLOC_FAIL; return false; }
	return true;
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

void Solver::simplify()
{
	if (alldisabled()) return;
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));
	assert(stats.clauses.original);
	stats.sigma.calls++;
	simplifying();
	INCREASE_LIMIT(this, sigma, stats.sigma.calls, nlognlogn, true);
	last.sigma.reduces = stats.reduces + 1;
	if (opts.phases > 2) {
		opts.phases--;
		LOG2(2, "  simplify phases decreased to %d", opts.phases);
	}
}

void Solver::awaken()
{
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));
	assert(sp->propagated == trail.size());
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
	assert(inf.nDualVars);
	if (!cumm.allocHist(cuhist, opts.proof_en) ||
		!cumm.allocAux(numCls) ||
		!cumm.allocVars(vars, savedLits) ||
		!cumm.allocPinned(vars, cuhist) ||
		!cumm.resizeCNF(cnf, numCls, numLits)) { simpstate = CNFALLOC_FAIL; return; }
	olist_cmp.init(cnf), compact_cmp.init(cnf);
	printStats(1, '-', CGREEN0);
	inf.nClauses = inf.nLiterals = 0;
	wt.clear(true);
	if (gopts.unified_access) {
		LOGN2(2, " Extracting clauses directly to device..");
		if (gopts.profile_gpu) cutimer->start();
		extract(cnf, orgs), orgs.clear(true);
		extract(cnf, learnts), learnts.clear(true);
		cm.destroy();
		assert(inf.nClauses == cnf->size());
		if (gopts.profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
	}
	else {
		LOGN2(2, " Extracting clauses heterogeneously to device..");
		if (gopts.profile_gpu) cutimer->start(streams[0]);
		cumm.createMirror(hcnf, maxClauses(), maxLiterals());
		extract(hcnf, orgs), reflectCNF(streams[0], streams[1]), orgs.clear(true);
		extract(hcnf, learnts), reflectCNF(streams[0], streams[1]), learnts.clear(true);
		cm.destroy();
		SYNC(streams[0]); SYNC(streams[1]);
		cumm.resizeCNFAsync(cnf, hcnf->data().size, hcnf->size());
		assert(hcnf->data().size == dataoff);
		assert(inf.nClauses == hcnf->size() && hcnf->size() == csoff);
		if (gopts.profile_gpu) cutimer->stop(streams[0]), cutimer->io += cutimer->gpuTime();
	}
	LOGENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
    vars->varcore = gopts.hostKOpts.ve_fun_en ? vars->eligible : NULL;
	SYNC(0); // sync device CNF resize
	initDevVorg(cuhist); // load device pointers to constant memory
	prepareCNFAsync(cnf, streams[0]);
	if (opts.proof_en) 
		cuMemSetAsync(cuhist.d_lbyte, 0, inf.nDualVars);
	assert(vorg.size() == inf.maxVar + 1);
	cuhist.fetchVars(vorg, streams[1]);
	if (opts.proof_en) {
		const uint32 *literals = flattenCNF(inf.nLiterals);
		uint32 proofcap = cuproof.count(literals, inf.nLiterals);
		proofcap *= 1.5;
		if (!cuproof.alloc(proofcap)) { simpstate = AWAKEN_FAIL; return; }
	}
}

void Solver::simplifying()
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
	timer.solve += timer.cpuTime();
	timer.start();
	awaken();
	if (simpstate == CNFALLOC_FAIL) {
		recycle();
		return;
	}
	/********************************/
	/*      reduction phases        */
	/********************************/
	assert(!phase && !multiplier);
	const int64 bmelted = inf.maxMelted, bclauses = inf.nClauses;
	int64 cdiff = INT64_MAX, ldiff = INT64_MAX;
	int64 clsbefore = inf.nClauses, litsbefore = inf.nLiterals;
	while (inf.nClauses && inf.nLiterals && !simpstate && !interrupted()) {
		if (!reallocOT(streams[0])) break;
		SYNC(streams[0]);
		reallocCNF();
		createOTAsync(cnf, ot, 0);
		if (!prop()) killSolver();
		if (!LCVE()) break;
		sortOT();
		if (stop(cdiff, ldiff)) { ERE(); break; }
		SUB(), VE(), BCE();
		cuproof.cacheProof(streams[4]);
		countAll();
		updateNumPVs();
		cacheNumUnits(streams[3]);
		inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
		cdiff = clsbefore - inf.nClauses, clsbefore = inf.nClauses;
		ldiff = litsbefore - inf.nLiterals, litsbefore = inf.nLiterals;
		cacheUnits(streams[3]);
		phase++, multiplier++;
		multiplier += phase == opts.phases;
		cuproof.writeProof(streams[4]);
		flattened = false;
	}
	/********************************/
	/*          Write Back          */
	/********************************/
	assert(sp->propagated == trail.size());
	cacheEliminated(streams[5]);             // if ERE is enabled, this transfer would be already done
	cumm.updateMaxCap(), cacher.updateMaxCap(); // for reporting GPU memory only
	markEliminated(streams[5]);              // must be executed before map()
	cacheCNF(streams[0], streams[1]);
	if (!propFailed()) killSolver();         // prop any pending units via host memory if 'realloc' failed
	bool success = (bclauses != inf.nClauses);
	inf.maxMelted += vars->nMelted;
	stats.sigma.all.clauses += bclauses - int64(inf.nClauses);
	stats.sigma.all.variables += int64(inf.maxMelted) - bmelted;
	last.shrink.removed = stats.shrunken;
	if (ereCls > inf.nClauses) stats.sigma.ere.removed += ereCls - inf.nClauses;
	if (inf.maxFrozen > sp->simplified) stats.units.forced += inf.maxFrozen - sp->simplified;
	if (!inf.unassigned || !inf.nClauses) { 
		LOG2(2, " All clauses removed");
		stats.clauses.original = 0;
		stats.clauses.learnt = 0;
		stats.literals.original = 0;
		stats.literals.learnt = 0;
		cnfstate = SAT; 
		printStats(1, 's', CGREEN); 
		return;
	}
	if (canMap()) map(true);
	else newBeginning();
	rebuildWT(opts.sigma_priorbins);
	if (BCP()) {
		LOG2(1, " Propagation after simplify proved a contradiction");
		learnEmpty();
	}
	UPDATE_SLEEPER(this, sigma, success);
	printStats(1, 's', CGREEN);
	timer.stop(), timer.simp += timer.cpuTime();
	if (!opts.solve_en) killSolver();
	timer.start();
}

void Solver::optSimp()
{
	LOGN2(2, " Initializing device options..");
	gopts.init(opts.proof_en && !opts.proof_nonbinary_en);
    assert(SH_MAX_BVE_OUT1 >= gopts.hostKOpts.xor_max_arity);
    assert(SH_MAX_ERE_OUT >= gopts.hostKOpts.ere_clause_max);
	cutimer = new cuTIMER;
	initDevOpts();
	LOGDONE(2, 5);
	initSharedMem();
}

void Solver::masterFree()
{
	LOG2(2, " Freeing up GPU memory..");
	SYNCALL;
	cumm.freeFixed();
	cumm.freePinned();
	cleanManaged();
	destroyStreams();
	cumm.breakMirror(), hcnf = NULL;
	if (opts.proof_en) 
		cuproof.destroy();
}

void Solver::slavesFree()
{

}