/***********************************************************************[solve.cpp]
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

#include "solve.h" 
#include "control.h"

namespace pFROST {
	CNF_INFO inf;
	ParaFROST* pfrost = NULL;
}

using namespace pFROST;

ParaFROST::ParaFROST(const string& _path) :
	formula(_path)
	, sp(NULL)
	, vsids(VSIDS_CMP(activity))
	, vschedule(SCORS_CMP(this))
	, bumped(0)
	, conflict(NOREF)
	, ignore(NOREF)
	, cnfstate(UNSOLVED_M)
	, intr(false)
	, stable(false)
	, probed(false)
	, incremental(false)
	, mapped(false)
	, simpstate(AWAKEN_SUCC)
{
	getCPUInfo(stats.sysmem);
	getBuildInfo();
	initSolver();
	if (!parser() || BCP()) { learnEmpty(), killSolver(); }
	if (opts.parseonly_en) killSolver();
}

void ParaFROST::allocSolver()
{
	PFLOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	assert(sizeof(LIT_ST) == 1);
	assert(inf.maxVar);
	assert(sp == NULL);
	const uint32 maxSize = inf.maxVar + 1;
	const C_REF initcap = inf.nOrgCls * sizeof(CLAUSE) + maxSize;
	sp = new SP(maxSize);
	sp->initSaved(opts.polarity);
	cm.init(initcap);
	wt.resize(inf.nDualVars);
	trail.reserve(inf.maxVar);
	dlevels.reserve(inf.maxVar);
	activity.resize(maxSize, 0.0);
	bumps.resize(maxSize, 0);
	PFLDONE(2, 5);
	PFLMEMCALL(this, 2);
}

void ParaFROST::initSolver()
{
	assert(!ORIGINAL && LEARNT && DELETED);
	assert(FROZEN_M && MELTED_M && SUBSTITUTED_M);
	assert(UNDEFINED < 0);
	assert(UNSAT == 0);
	assert(SAT == 1);
	assert(UNSOLVED(cnfstate));
	forceFPU();
	opts.init();
	subbin.resize(2);
	dlevels.push(0);
	if (opts.proof_en) {
#ifdef _WIN32
		if (!opts.proof_nonbinary_en) {
			PFLOG2(1, "  Disabling DRAT binary mode on Windows");
			opts.proof_nonbinary_en = true;
		}
#endif
		proof.handFile(opts.proof_path, opts.proof_nonbinary_en);
	}
}

void ParaFROST::initLimits() 
{
	PFLOG2(2, " Initializing solver limits..");
	formula.c2v = ratio(double(stats.clauses.original), double(inf.maxVar));
	if (opts.ternary_en && formula.ternaries < 2) {
		PFLOG2(2, "  Disabling hyper ternary resolution as no ternaries found");
		opts.ternary_en = false;
	}
	if (!last.vsids.inc && !last.vsids.booster) {
		last.vsids.inc = opts.var_inc;
		last.vsids.booster = 1.0 / opts.var_decay;
	}
	last.transitive.literals = 2;
	INIT_LIMIT(limit.reduce, opts.reduce_inc, false);
	INIT_LIMIT(limit.rephase, opts.rephase_inc, false);
	INIT_LIMIT(limit.mode.conflicts, opts.mode_inc, false);
	INIT_LIMIT(limit.mdm, opts.mdm_inc, true);
	INIT_LIMIT(limit.probe, opts.probe_inc, true);
	INIT_LIMIT(limit.sigma, opts.sigma_inc, true);
	INIT_LIMIT(limit.subsume, opts.subsume_inc, true);
	lbdrest.init(opts.lbd_rate, opts.lbd_fast, opts.lbd_slow);
	lbdrest.reset();
	stable = opts.stable_en && opts.vsidsonly_en;
	if (stable) {
		PFLOG2(2, "  VSIDS with initial stable mode is enabled");
		lubyrest.enable(opts.luby_inc, opts.luby_max);
	}
	else {
		PFLOG2(2, "  VMFQ with initial unstable mode is enabled");
		updateUnstableLimit();
	}
	if (opts.mdm_rounds) {
		PFLOG2(2, "  Enabling MDM with %d rounds", opts.mdm_rounds);
		last.mdm.rounds = opts.mdm_rounds;
	}
	if (opts.seed) {
		PFLOG2(2, "  Initial random seed is set to %d", opts.seed);
		random.init(opts.seed);
	}
	PFLOG2(2, " Limits initialized successfully");
}

void ParaFROST::solve()
{
	FAULT_DETECTOR;
	timer.start();
	initLimits();
	if (verbose == 1) printTable();
	if (canPreSigmify()) sigmify();
	if (UNSOLVED(cnfstate)) {
		PFLOG2(2, "-- CDCL search started..");
		MDMInit();
		while (UNSOLVED(cnfstate) && !runningout()) {
			if (BCP()) analyze();
			else if (!inf.unassigned) cnfstate = SAT;
			else if (canReduce()) reduce();
			else if (canRestart()) restart();
			else if (canRephase()) rephase();
			else if (canSigmify()) sigmify();
			else if (canProbe()) probe();
			else if (canMMD()) MDM();
			else decide();
		}
		PFLOG2(2, "-- CDCL search completed successfully");
	}
	timer.stop(), timer.solve += timer.cpuTime();
	wrapup();
}

void ParaFROST::wrapup() 
{
	if (!quiet_en) { PFLRULER('-', RULELEN); PFLOG0(""); }
	if (cnfstate == SAT) {
		PFLOGS("SATISFIABLE");
		assert(sp != NULL && sp->value != NULL);
		if (opts.model_en) {
			model.extend(sp->value);
			if (opts.modelprint_en) 
				model.print();
		}
		if (opts.modelverify_en) {
			model.extend(sp->value);
			model.verify(formula.path);
		}
	}
	else if (cnfstate == UNSAT) PFLOGS("UNSATISFIABLE");
	else if (UNSOLVED(cnfstate)) PFLOGS("UNKNOWN");
	if (opts.report_en) report();
}