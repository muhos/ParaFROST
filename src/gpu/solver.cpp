/***********************************************************************[solver.cpp]
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
#include "control.hpp"
#include "options.cuh"
#include "timer.cuh"
#include "banner.hpp"

namespace ParaFROST {
	GOPTION			gopts;
	CNF_INFO		inf;
	cudaDeviceProp	devProp;
	uint32			maxGPUThreads = 0;
	size_t			maxGPUSharedMem = 0;
}

using namespace ParaFROST;


Solver::~Solver() 
{ 
	freeSimp();
    std::free(this->learnCallbackBuffer);
}

Solver::Solver(const std::string& path) : Solver()
{
    formula = path;
	LOGHEADER(1, 5, "Parser")
	if (!parse() || BCP()) { assert(cnfstate == UNSAT), killSolver(); }
	if (opts.parseonly_en) killSolver();
}

Solver::Solver()
  : formula()
  , timer()
  , cutimer()
  , sp(nullptr)
  , vsids(VSIDS_CMP(activity))
  , vschedule(SCORS_CMP(this))
  , bumped(0)
  , conflict(NOREF)
  , ignore(NOREF)
  , cnfstate(UNSOLVED)
  , intr(false)
  , stable(false)
  , probed(false)
  , incremental(false)
  , vars(nullptr)
  , ot(nullptr)
  , cnf(nullptr)
  , hcnf(nullptr)
  , tca(cacher)
  , cuproof(cumm, proof)
  , streams(nullptr)
  , mapped(false)
  , compacted(false)
  , flattened(false)
  , phase(0)
  , nForced(0)
  , simpstate(AWAKEN_SUCC)
  , devCount(0)
  , termCallbackState(nullptr)
  , learnCallbackState(nullptr)
  , learnCallbackBuffer(nullptr)
  , termCallback(nullptr)
  , learnCallback(nullptr)
{
    initialize(false);
}

void Solver::initialize(const bool& banner)
{
	if (banner) {
		LOGHEADER(1, 5, "Banner");
		LOGFANCYBANNER(version());
	}
	LOGHEADER(1, 5, "Build")
	getCPUInfo(stats.sysmem);
	getBuildInfo();
	initSolver();
	size_t _gfree = 0, _gpenalty = 0;
	if (!(devCount = getGPUInfo(_gfree, _gpenalty))) {
		LOGERRORN("no GPU(s) available that support CUDA");
		killSolver();
	}
	if (_gfree <= 200 * MBYTE) {
		LOGWARNING("not enough GPU memory (free = %zd MB): skip simplifier", _gfree / MBYTE);
		opts.sigma_en = opts.sigma_live_en = false;
	}
	else cumm.init(_gfree, _gpenalty);
	if (cumm.checkMemAdvice()) {
        LOG2(2, " Enabled GPU driver memory advice");
	}
	else {
        LOG2(2, " Disabled GPU driver memory advice");
	}
	if (opts.sigma_en || opts.sigma_live_en) { optSimp(), createStreams(); }
}

void Solver::allocSolver()
{
	LOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	assert(sizeof(LIT_ST) == 1);
	assert(inf.maxVar);
	assert(sp == NULL);
	const uint32 maxSize = inf.maxVar + 1;
	const C_REF initcap = inf.orgCls * sizeof(CLAUSE) + maxSize;
	sp = new SP(maxSize, opts.polarity);
	cm.init(initcap);
	wt.resize(inf.maxDualVars);
	trail.reserve(inf.maxVar);
	dlevels.reserve(inf.maxVar);
	activity.resize(maxSize, 0.0);
	bumps.resize(maxSize, 0);
	LOGDONE(2, 5);
	LOGMEMCALL(this, 2);
}

void Solver::initSolver()
{
	assert(!ORIGINAL && LEARNT && DELETED);
	assert(FROZEN_M && MELTED_M && SUBSTITUTED_M);
	assert(UNDEFINED == -1);
	assert(UNSAT == 0);
	assert(SAT == 1);
	assert(IS_UNSOLVED(cnfstate));
	forceFPU();
	opts.init();
	subbin.resize(2);
	dlevels.push(0);
	if (opts.proof_en) {
#ifdef _WIN32
		if (!opts.proof_nonbinary_en) {
			LOG2(1, " Disabling DRAT binary mode on Windows.");
			opts.proof_nonbinary_en = true;
		}
#endif
		proof.handFile(opts.proof_path, opts.proof_nonbinary_en);
	}
}

void Solver::initLimits() 
{
	LOG2(2, " Initializing solver limits..");
	formula.c2v = ratio(double(stats.clauses.original), double(inf.maxVar));
	if (opts.ternary_en && formula.ternaries < 2) {
		LOG2(2, "  Disabling hyper ternary resolution as no ternaries found");
		opts.ternary_en = false;
	}
	if (!last.vsids.inc && !last.vsids.booster) {
		last.vsids.inc = opts.var_inc;
		last.vsids.booster = 1.0 / opts.var_decay;
	}
	last.transitive.literals = 2;
	INIT_LIMIT(this, limit.reduce, opts.reduce_inc, false);
	INIT_LIMIT(this, limit.rephase, opts.rephase_inc, false);
	INIT_LIMIT(this, limit.mode.conflicts, opts.mode_inc, false);
	INIT_LIMIT(this, limit.mdm, opts.mdm_inc, true);
	INIT_LIMIT(this, limit.probe, opts.probe_inc, true);
	INIT_LIMIT(this, limit.sigma, opts.sigma_inc, true);
	INIT_LIMIT(this, limit.subsume, opts.subsume_inc, true);
	lbdrest.init(opts.lbd_rate, opts.lbd_fast, opts.lbd_slow);
	lbdrest.reset();
	stable = opts.stable_en && opts.vsidsonly_en;
	if (stable) {
		LOG2(2, "  VSIDS with initial stable mode is enabled");
		lubyrest.enable(opts.luby_inc, opts.luby_max);
	}
	else {
		LOG2(2, "  VMFQ with initial unstable mode is enabled");
		updateUnstableLimit();
	}
	if (opts.mdm_rounds) {
		LOG2(2, "  Enabling MDM with %d rounds", opts.mdm_rounds);
		last.mdm.rounds = opts.mdm_rounds;
	}
	if (opts.seed) {
		LOG2(2, "  Initial random seed is set to %d", opts.seed);
		random.init(opts.seed);
	}
	LOG2(2, " Limits initialized successfully");
}

void Solver::solve()
{
	FAULT_DETECTOR;
    LOGHEADER(1, 5, "Search");
	if (incremental) {
		LOG0("");
		LOGWARNING("isolve() should be called in incremental mode");
		return;
	}
	timer.start();
	initLimits();
	if (verbose == 1) printTable();
	if (canPreSimplify()) simplify();
	if (IS_UNSOLVED(cnfstate)) {
		MDMInit();
		while (IS_UNSOLVED(cnfstate) && !runningout()) {
			if (BCP()) analyze();
			else if (!inf.unassigned) cnfstate = SAT;
			else if (canReduce()) reduce();
			else if (canRestart()) restart();
			else if (canRephase()) rephase();
			else if (canSigmify()) simplify();
			else if (canProbe()) probe();
			else if (canMMD()) MDM();
			else decide();
		}
	}
	timer.stop(), stats.time.solve += timer.cpuTime();
	wrapup();
}

void Solver::wrapup(const CNFState& state) {
    LOGHEADER(1, 5, "Result");
	if (state != UNSOLVED) cnfstate = state;
	if (cnfstate == SAT) {
		LOGSAT("SATISFIABLE");
		assert(sp != NULL && sp->value != NULL);
		if (opts.model_en) {
			model.extend(sp->value);
			if (opts.modelprint_en) model.print();
		}
		if (opts.modelverify_en) {
			model.extend(sp->value);
			model.verify(formula.path);
		}
	}
	else if (cnfstate == UNSAT) LOGSAT("UNSATISFIABLE");
	else if (IS_UNSOLVED(cnfstate)) LOGSAT("UNKNOWN");
	if (opts.report_en) report();
}