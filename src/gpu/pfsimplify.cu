/***********************************************************************[pfsimplify.cu]
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

#include "pfsolve.h"

using namespace pFROST;
using namespace SIGmA;

inline void	ParaFROST::cleanManaged()
{
	if (vars != NULL) delete vars;
	cumm.freeVars(), cumm.freeCNF(), cumm.freeOT();
	vars = NULL, cnf = NULL, ot = NULL;
}

inline bool	ParaFROST::reallocCNF()
{
	int times = phase + 1;
	if (times > 1 && times != opts.phases && (times % opts.shrink_rate) == 0) {
		size_t maxAddedCls = opts.ve_en ? inf.nClauses : 0;
		size_t maxAddedLits = opts.ve_en ? size_t(stats.literals.original * opts.lits_mul) : 0;
		PFLOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
		if (!cumm.resizeCNF(cnf, inf.nClauses + maxAddedCls, inf.nLiterals + maxAddedLits)) {
			sigState = CNFALLOC_FAIL, compacted = false;
			return false;
		}
		compacted = true;
	}
	else cumm.cacheCNFPtr(cnf), compacted = false;
	return true;
}

inline bool	ParaFROST::reallocOT(const cudaStream_t& stream)
{
	assert(inf.nLiterals);
	calcOccurs(inf.nLiterals);
	if (!cumm.resizeOTAsync(ot, inf.nLiterals, stream)) { sigState = OTALLOC_FAIL; return false; }
	return true;
}

inline void	ParaFROST::initSimplifier()
{
	cleanManaged();
	cumm.freeFixed();
	tca.destroy();
	sigState = AWAKEN_SUCC;
	dataoff = csoff = 0;
	phase = mu_inc = 0;
	ereCls = nForced = 0;
	compacted = false;
	assert(hcnf == NULL);
	if (stats.sigma.calls > 1) {
		// memory must be cleaned before this inquiry
		size_t free = 0, tot = 0;
		CHECK(cudaMemGetInfo(&free, &tot));
		cumm.reset(free);
	}
}

void ParaFROST::awaken()
{
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	assert(sp->propagated == trail.size());
	initSimplifier();
	// alloc simplifier memory 
	size_t numCls = maxClauses(), numLits = maxLiterals();
	if (opts.phases) {
		size_t maxAddedCls = opts.ve_en ? stats.clauses.original : 0;
		size_t maxAddedLits = opts.ve_en ? size_t(stats.literals.original * opts.lits_mul) : 0;
		PFLOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
		numCls += maxAddedCls, numLits += maxAddedLits;
	}
	assert(inf.nDualVars);
	if (!cumm.allocHist(cuhist) ||
		!cumm.allocAux(numCls) ||
		!cumm.allocVars(vars, numLits) ||
		!cumm.resizeCNF(cnf, numCls, numLits)) { sigState = CNFALLOC_FAIL; return; }
	printStats(1, '-', CGREEN0), inf.nClauses = inf.nLiterals = 0;
	wt.clear(true);
	if (unified_access) {
		PFLOGN2(2, " Extracting clauses directly to device..");
		if (profile_gpu) cutimer->start();
		extract(cnf, orgs), orgs.clear(true);
		extract(cnf, learnts), learnts.clear(true);
		cm.destroy();
		assert(inf.nClauses == cnf->size());
		if (profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
	}
	else {
		PFLOGN2(2, " Extracting clauses heterogeneously to device..");
		if (profile_gpu) cutimer->start(streams[0]);
		cumm.createMirror(hcnf, maxClauses(), maxLiterals());
		extract(hcnf, orgs), reflectCNF(streams[0], streams[1]), orgs.clear(true);
		extract(hcnf, learnts), reflectCNF(streams[0], streams[1]), learnts.clear(true);
		cm.destroy();
		sync(streams[0]), sync(streams[1]);
		cumm.resizeCNFAsync(cnf, hcnf->data().size, hcnf->size());
		assert(hcnf->data().size == dataoff);
		assert(inf.nClauses == hcnf->size() && hcnf->size() == csoff);
		if (profile_gpu) cutimer->stop(streams[0]), cutimer->io += cutimer->gpuTime();
	}
	PFLENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
	sync(); // sync device CNF resize
	prepareCNFAsync(cnf, streams[0]);
	assert(vorg.size() == inf.maxVar + 1);
	cuhist.fetchVars(vorg, streams[1]);
}

void ParaFROST::sigmify()
{
	if (alldisabled()) return;
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	assert(stats.clauses.original);
	stats.sigma.calls++;
	sigmifying();
	INCREASE_LIMIT(this, sigma, stats.sigma.calls, nlognlogn, true);
	last.sigma.reduces = stats.reduces + 1;
	if (opts.phases > 2) {
		opts.phases--;
		PFLOG2(2, "  sigmify phases decreased to %d", opts.phases);
	}
}

void ParaFROST::sigmifying()
{
	/********************************/
	/*		Getting ready...        */
	/********************************/
	SLEEPING(sleep.sigma, opts.sigma_sleep_en);
	rootify();
	assert(cnfstate == UNSOLVED);
	PFLOGN2(2, " Shrinking CNF before sigmification..");
	if (sp->simplified < inf.maxFrozen) sp->simplified = inf.maxFrozen;
	int64 clsbefore = maxClauses(), litsbefore = maxLiterals();
	shrinkTop(orgs), shrinkTop(learnts);
	PFLSHRINKALL(this, 2, clsbefore, litsbefore);
	assert(stats.clauses.original == orgs.size());
	assert(stats.clauses.learnt == learnts.size());
	if (orgs.empty()) {
		recycleWT();
		return;
	}
	timer.stop();
	timer.solve += timer.cpuTime();
	timer.start();
	awaken();
	if (sigState == CNFALLOC_FAIL) {
		recycle();
		return;
	}
	/********************************/
	/*      reduction phases        */
	/********************************/
	assert(!phase && !mu_inc);
	const int64 bmelted = inf.maxMelted, bclauses = inf.nClauses, bliterals = inf.nLiterals;
	int64 cdiff = INT64_MAX, ldiff = INT64_MAX;
	clsbefore = inf.nClauses, litsbefore = inf.nLiterals; 
	while (!interrupted()) {
		if (!reallocOT(streams[0])) break;
		sync(streams[0]);
		reallocCNF();
		createOTAsync(cnf, ot, 0);
		if (!prop()) killSolver();
		if (!LCVE()) break;
		sortOTAsync(cnf, ot, vars, streams);
		if (sigState == CNFALLOC_FAIL) {
			HSE();
			cacheNumUnits(streams[3]);
			cacheUnits(streams[3]);
			ERE();
			break;
		}
		if (stop(cdiff, ldiff)) { ERE(); break; }
		HSE(), VE(), BCE();
		countAll();
		cacheNumUnits(streams[3]);
		inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
		cdiff = clsbefore - inf.nClauses, clsbefore = inf.nClauses;
		ldiff = litsbefore - inf.nLiterals, litsbefore = inf.nLiterals;
		cacheUnits(streams[3]);
		phase++, mu_inc++;
		mu_inc += phase == opts.phases;
	}
	/********************************/
	/*          Write Back          */
	/********************************/
	assert(sp->propagated == trail.size());
	cumm.updateMaxCap(), tca.updateMaxCap(); // for reporting GPU memory only
	cacheCNF(streams[0], streams[1]);
	if (!propFailed()) killSolver(); // prop any pending units via host memory if 'realloc' failed
	bool success = (bliterals != inf.nLiterals);
	stats.sigma.all.variables += int64(inf.maxMelted) - bmelted;
	stats.sigma.all.clauses += bclauses - int64(inf.nClauses);
	stats.sigma.all.literals += bliterals - int64(inf.nLiterals);
	last.shrink.removed = stats.shrunken;
	if (ereCls > inf.nClauses) stats.sigma.ere.removed += ereCls - inf.nClauses;
	if (inf.maxFrozen > sp->simplified) stats.units.forced += inf.maxFrozen - sp->simplified;
	if (!inf.unassigned || !inf.nClauses) { cnfstate = SAT; printStats(1, 'p'); return; }
	if (canMap()) map(true);
	else newBeginning();
	wt.resize(inf.nDualVars);
	rebuildWT(opts.sigma_priorbins);
	if (BCP()) {
		PFLOG2(1, " Propagation after sigmify proved a contradiction");
		learnEmpty();
	}
	UPDATE_SLEEPER(this, sigma, success);
	printStats(1, 's', CGREEN);
	timer.stop(), timer.simp += timer.cpuTime();
	if (!opts.solve_en) killSolver();
	timer.start();
}

void ParaFROST::masterFree()
{
	syncAll();
	cumm.freeFixed();
	cumm.freePinned();
	tca.destroy();
	cleanManaged();
	destroyStreams();
	cumm.breakMirror(), hcnf = NULL;
}

void ParaFROST::slavesFree()
{

}

void ParaFROST::optSimp()
{
	PFLOGN2(2, " Initializing simplifier limits..");
	if (opts.phases > 1 && formula.size > 700 * MBYTE &&
		(formula.c2v > 4 || formula.maxClauseSize > 10000)) {
		if (formula.size > GBYTE) 
			opts.phases = 1, opts.ve_clause_limit = 50;
		else 
			opts.phases = 2, opts.ve_clause_limit = 20;
	}
	culimit.limits[0] = opts.hse_limit;
	culimit.limits[1] = opts.bce_limit;
	culimit.limits[2] = opts.ere_limit;
	culimit.limits[3] = opts.xor_max_arity;
	culimit.limits[4] = opts.ve_clause_limit;
	culimit.limits[5] = opts.ve_lbound_en;
	if (opts.ngpus > devCount) opts.ngpus = devCount;
	if (profile_gpu = opts.profile_simp) cutimer = new cuTIMER;
	initConstants(culimit);
	PFLDONE(2, 5);
#ifdef _DEBUG
	if (opts.ere_en) {
		PFLOGN2(2, " Disabling ERE in debug mode..");
		opts.ere_en = false;
		PFLDONE(2, 5);
	}
#endif
}