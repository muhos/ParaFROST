/***********************************************************************[pfsolve.cpp]
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

namespace pFROST {
	CNF_INFO inf;
	ParaFROST* pfrost = NULL;
}

using namespace pFROST;

ParaFROST::ParaFROST(const string& _path) :
	formula(_path)
	, vschedule(SCORS_CMP(this))
	, vsids(VSIDS_CMP(activity))
	, sp(NULL)
	, bumped(0)
	, intr(false)
	, stable(false)
	, probed(false)
	, mapped(false)
	, ignore(NOREF)
	, conflict(NOREF)
	, cnfstate(UNSOLVED)
	, sigState(AWAKEN_SUCC)
	, incremental(false)
{
	stats.sysmem = getAvailSysMem();
	getCPUInfo();
	PFLOG2(1, " Available system memory: %lld GB", stats.sysmem / GBYTE);
	initSolver();
	if (!parser() || BCP()) { learnEmpty(), killSolver(); }
	if (opts.parseonly_en) killSolver();
	if (opts.ternary_en && formula.ternaries < 2) opts.ternary_en = false;
	if (verbose == 1) printTable();
}

void ParaFROST::allocSolver()
{
	PFLOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	assert(sizeof(LIT_ST) == 1);
	assert(inf.maxVar);
	assert(sp == NULL);
	uint32 maxSize = inf.maxVar + 1;
	// search space
	sp = new SP(maxSize);
	sp->initSaved(opts.polarity);
	// input database
	cm.init(inf.nOrgCls);
	// watch table
	wt.resize(inf.nDualVars);
	// variable arrays
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
	assert(cnfstate);
	opts.init();
	subbin.resize(2);
	dlevels.push(0);
	if (opts.seed) random.init(opts.seed);
	if (opts.proof_en)
		proof.handFile(opts.proof_path, false);
}

void ParaFROST::initLimits() {
	PFLOG2(2, " Initializing solver limits..");
	formula.c2v = ratio(double(stats.clauses.original), double(inf.maxVar));
	if (!last.vsids.inc && !last.vsids.decay) {
		last.vsids.inc = opts.var_inc;
		last.vsids.decay = opts.var_decay;
	}
	last.transitive.literals = 2;
	last.mdm.rounds = opts.mdm_rounds;
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
		PFLOG2(2, "  VSIDS with initial stable mode is enabled");
		lubyrest.enable(opts.luby_inc, opts.luby_max);
	}
	else {
		PFLOG2(2, "  VMFQ with initial unstable mode is enabled");
		updateUnstableLimit();
	}
	if (opts.sigma_en && opts.phases > 1 && formula.size > 700 * MBYTE &&
		(formula.c2v > 4 || formula.maxClauseSize > 10000)) {
		if (formula.size > GBYTE)
			opts.phases = 1, opts.ve_clause_limit = 50;
		else
			opts.phases = 2, opts.ve_clause_limit = 20;
	}
	PFLOG2(2, " Limits initialized successfully");
}

void ParaFROST::solve()
{
	timer.start();
	initLimits();
	if (canPreSigmify()) sigmify();
	PFLOG2(2, "-- CDCL search started..");
	if (cnfstate == UNSOLVED) MDMInit();
	while (cnfstate == UNSOLVED && !interrupted()) {
		PFLDL(this, 3);
		if (BCP()) analyze();
		else if (!inf.unassigned) cnfstate = SAT;
		else if (canReduce()) reduce();
		else if (canRestart()) restart();
		else if (canRephase()) rephase();
		else if (canSigmify()) sigmify();
		else if (canProbe()) probe();
		else if (canMMD()) MDM();
		else decide();
		PFLTRAIL(this, 3);
	}
	timer.stop(), timer.solve += timer.cpuTime();
	PFLOG2(2, "-- CDCL search completed successfully");
	wrapup();
}

void ParaFROST::report()
{
	if (opts.report_en) {
		PFLOG0("");
		if (opts.sigma_en || opts.sigma_live_en) {
			PFLOG1("\t\t\t%sSimplifier Report%s", CREPORT, CNORMAL);
			if (!opts.profile_simp)
				PFLOG1(" %sSimplifier time        : %s%-16.3f  sec%s", CREPORT, CREPORTVAL, timer.simp, CNORMAL);
			else {
				PFLOG1(" %s - Var ordering        : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.vo, CNORMAL);
				PFLOG1(" %s - CNF compact         : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.gc, CNORMAL);
				PFLOG1(" %s - OT  creation        : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.cot, CNORMAL);
				PFLOG1(" %s - OT  sorting         : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.sot, CNORMAL);
				PFLOG1(" %s - OT  reduction       : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.rot, CNORMAL);
				PFLOG1(" %s - BVE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.ve, CNORMAL);
				PFLOG1(" %s - HSE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.hse, CNORMAL);
				PFLOG1(" %s - BCE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.bce, CNORMAL);
				PFLOG1(" %s - ERE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.ere, CNORMAL);
			}
			PFLOG1(" %sSigmifications         : %s%-10d%s", CREPORT, CREPORTVAL, stats.sigma.calls, CNORMAL);
			PFLOG1(" %s Forced units          : %s%-10d%s", CREPORT, CREPORTVAL, stats.units.forced, CNORMAL);
			PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.variables + stats.units.forced, CNORMAL);
			PFLOG1(" %s Removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.clauses, CNORMAL);
			PFLOG1(" %s Removed literals      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.literals, CNORMAL);
			PFLOG1(" %s Tried redundancies    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.tried, CNORMAL);
			PFLOG1(" %s Original redundancies : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.orgs, CNORMAL);
			PFLOG1(" %s Learnt redundancies   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.learnts, CNORMAL);
		}
		PFLOG1("\t\t\t%sSolver Report%s", CREPORT, CNORMAL);
		PFLOG1(" %sSolver time            : %s%-16.3f  sec%s", CREPORT, CREPORTVAL, timer.solve, CNORMAL);
		PFLOG1(" %sSystem memory          : %s%-16.3f  MB%s", CREPORT, CREPORTVAL, ratio(double(sysMemUsed()), double(MBYTE)), CNORMAL);
		PFLOG1(" %sFormula                : %s%-s%s", CREPORT, CREPORTVAL, formula.path.c_str(), CNORMAL);
		PFLOG1(" %s Size                  : %s%-16.3f  MB%s", CREPORT, CREPORTVAL, ratio(double(formula.size), double(MBYTE)), CNORMAL);
		PFLOG1(" %s Units                 : %s%-10d%s", CREPORT, CREPORTVAL, formula.units, CNORMAL);
		PFLOG1(" %s Binaries              : %s%-10d%s", CREPORT, CREPORTVAL, formula.binaries, CNORMAL);
		PFLOG1(" %s Ternaries             : %s%-10d%s", CREPORT, CREPORTVAL, formula.ternaries, CNORMAL);
		PFLOG1(" %s Larger                : %s%-10d%s", CREPORT, CREPORTVAL, formula.large, CNORMAL);
		PFLOG1(" %s Max clause size       : %s%-10d%s", CREPORT, CREPORTVAL, formula.maxClauseSize, CNORMAL);
		PFLOG1(" %s C2V ratio             : %s%-10.3f%s", CREPORT, CREPORTVAL, formula.c2v, CNORMAL);
		PFLOG1(" %sBacktracks             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.chrono + stats.backtrack.nonchrono, CNORMAL);
		PFLOG1(" %s Chronological         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.chrono, CNORMAL);
		PFLOG1(" %s Non-Chronological     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.nonchrono, CNORMAL);
		PFLOG1(" %s Trail reuses          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.reuses, CNORMAL);
		PFLOG1(" %sConflicts              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.conflicts, CNORMAL);
		PFLOG1(" %s Learnt literals       : %s%-16lld  ( %2.2f %% deleted )%s", CREPORT, CREPORTVAL, stats.minimize.after, percent(stats.minimize.before - stats.minimize.after, stats.minimize.before), CNORMAL);
		PFLOG1(" %s Learnt units          : %s%-10d%s", CREPORT, CREPORTVAL, stats.units.learnt, CNORMAL);
		PFLOG1(" %s Learnt OTF subsumes   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subtried, CNORMAL);
		PFLOG1(" %s Learnt OTF subsumed   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.learntfly, CNORMAL);
		PFLOG1(" %s Alluip checks         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.alluip.checks, CNORMAL);
		PFLOG1(" %s Alluip successes      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.alluip.learnts, CNORMAL);
		PFLOG1(" %sDeduplications         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.calls, CNORMAL);
		PFLOG1(" %s Hyper unaries         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.hyperunary, CNORMAL);
		PFLOG1(" %s Duplicated binaries   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.binaries, CNORMAL);
		PFLOG1(" %sDecompositions         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.calls, CNORMAL);
		PFLOG1(" %s SCCs                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.scc, CNORMAL);
		PFLOG1(" %s Hyper unaries         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.hyperunary, CNORMAL);
		PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.variables, CNORMAL);
		PFLOG1(" %s Removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.clauses, CNORMAL);
		PFLOG1(" %sHyper binary resolves  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.resolutions, CNORMAL);
		PFLOG1(" %s Added binaries        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.resolvents, CNORMAL);
		PFLOG1(" %sHyper ternary resolves : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.resolutions, CNORMAL);
		PFLOG1(" %s Added binaries        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.binaries, CNORMAL);
		PFLOG1(" %s Added ternaries       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.ternaries, CNORMAL);
		PFLOG1(" %s Subsumed ternaries    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.binaries * 2, CNORMAL);
		PFLOG1(" %sReduces                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.reduces, CNORMAL);
		PFLOG1(" %s Removed binaries      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.reduced, CNORMAL);
		PFLOG1(" %s Removed ternaries     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.reduced, CNORMAL);
		PFLOG1(" %sRestarts               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.restart.all, CNORMAL);
		PFLOG1(" %s Stable restarts       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.restart.stable, CNORMAL);
		PFLOG1(" %s Stable modes          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.stablemodes, CNORMAL);
		PFLOG1(" %s Unstable modes        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.unstablemodes, CNORMAL);
		PFLOG1(" %sRephases               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.all, CNORMAL);
		PFLOG1(" %s Original              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.org, CNORMAL);
		PFLOG1(" %s Random                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.random, CNORMAL);
		PFLOG1(" %s Invert                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.inv, CNORMAL);
		PFLOG1(" %s Flip                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.flip, CNORMAL);
		PFLOG1(" %s Best                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.best, CNORMAL);
		PFLOG1(" %sRecyclings             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.soft + stats.recycle.hard, CNORMAL);
		PFLOG1(" %s Soft                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.soft, CNORMAL);
		PFLOG1(" %s Hard                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.hard, CNORMAL);
		PFLOG1(" %sProbes calls           : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.calls, CNORMAL);
		PFLOG1(" %s Rounds                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.rounds, CNORMAL);
		PFLOG1(" %s Probed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.probed, CNORMAL);
		PFLOG1(" %s Failed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.failed, CNORMAL);
		PFLOG1(" %s Propagations          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.propagations.probe, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probeticks, CNORMAL);
		PFLOG1(" %sTransitive calls       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.calls, CNORMAL);
		PFLOG1(" %s Probed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.probed, CNORMAL);
		PFLOG1(" %s Failed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.failed, CNORMAL);
		PFLOG1(" %s Transitives           : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.removed, CNORMAL);
		PFLOG1(" %s Propagations          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.propagations.transitive, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitiveticks, CNORMAL);
		PFLOG1(" %sShrinks                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.calls, CNORMAL);
		PFLOG1(" %s removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.clauses, CNORMAL);
		PFLOG1(" %s removed literals      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.literals, CNORMAL);
		PFLOG1(" %sSubsume calls          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.calls, CNORMAL);
		PFLOG1(" %s Checks                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.checks, CNORMAL);
		PFLOG1(" %s Subsumed              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.subsumed, CNORMAL);
		PFLOG1(" %s Strengthened          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.strengthened, CNORMAL);
		PFLOG1(" %sSearch decisions       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decisions.single, CNORMAL);
		PFLOG1(" %s Propagations          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.propagations.search, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.searchticks, CNORMAL);
		PFLOG1(" %sMappings               : %s%-10d%s", CREPORT, CREPORTVAL, stats.mappings, CNORMAL);
		PFLOG1(" %sMDM calls              : %s%-10d%s", CREPORT, CREPORTVAL, stats.mdmcalls, CNORMAL);
		PFLOG1(" %s Decisions             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decisions.multiple, CNORMAL);
		PFLOG1(" %sVivification checks    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.checks, CNORMAL);
		PFLOG1(" %s Assumed               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.assumed, CNORMAL);
		PFLOG1(" %s Reused                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.reused, CNORMAL);
		PFLOG1(" %s Vivified              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.vivified, CNORMAL);
		PFLOG1(" %s Implied               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.implied, CNORMAL);
		PFLOG1(" %s Subsumed              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.subsumed, CNORMAL);
		PFLOG1(" %s Strengthened          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.strengthened, CNORMAL);
	}
}

void ParaFROST::wrapup() {
	if (!quiet_en) { PFLRULER('-', RULELEN); PFLOG0(""); }
	if (cnfstate == SAT) {
		PFLOGS("SATISFIABLE");
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
	else if (cnfstate == UNSAT) PFLOGS("UNSATISFIABLE");
	else if (cnfstate == UNSOLVED) PFLOGS("UNKNOWN");
	if (opts.report_en) report();
}