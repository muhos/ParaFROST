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
#include "pfdimacs.h"

namespace pFROST {
	CNF_INFO inf;
	ParaFROST* pfrost = NULL;
}

using namespace pFROST;
//=======================================//
//		 ParaFROST defined members       //
//=======================================//
ParaFROST::ParaFROST(const string& _path) :
	path(_path)
	, vsids(HEAP_CMP(activity))
	, starts(1)
	, nConflicts(0)
	, intr(false)
	, mapped(false)
	, conflict(NOREF)
	, cnfstate(UNSOLVED)
	, sigState(AWAKEN_SUCC)
{
	opts.init();
	stats.sysMemAvail = getAvailSysMem();
	getCPUInfo();
	PFLOG2(1, " Available system memory = %lld GB", stats.sysMemAvail / GBYTE);
	if (opts.proof_en) PFLOGE("generating proof is currently unsupported");
	if (!parser() || BCP()) { cnfstate = UNSAT, killSolver(); }
	if (opts.parse_only_en) killSolver();
	if (verbose == 1) printTable();
}

bool ParaFROST::parser() {
	struct stat st;
	stat(path.c_str(), &st);
	size_t fsz = st.st_size;
	PFLOG2(1, " Parsing CNF file \"%s%s%s\" (size: %s%.2f MB%s)",
		CREPORTVAL, path.c_str(), CNORMAL, CREPORTVAL, ratio(fsz, uint64(MBYTE)), CNORMAL);
	timer.start();
#ifdef __linux__
	int fd = open(path.c_str(), O_RDONLY, 0);
	if (fd == -1) PFLOGE("cannot open input file");
	void* buffer = mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
	char* str = (char*)buffer;
#else
	ifstream inputFile;
	inputFile.open(path, ifstream::in);
	if (!inputFile.is_open()) PFLOGE("cannot open input file");
	char* buffer = new char[fsz + 1], * str = buffer;
	inputFile.read(buffer, fsz);
	buffer[fsz] = '\0';
#endif
	Lits_t in_c, org;
	in_c.reserve(INIT_CAP);
	org.reserve(INIT_CAP);
	PFLMEMCALL(this, 2);
	char* eof = str + fsz;
	while (str < eof) {
		eatWS(str);
		if (*str == '\0' || *str == '0' || *str == '%') break;
		if (*str == 'c') eatLine(str);
		else if (*str == 'p') {
			if (!eq(str, "p cnf")) PFLOGE("header has wrong format");
			uint32 sign = 0;
			inf.orgVars = inf.maxVar = toInteger(str, sign);
			if (sign) PFLOGE("number of variables in header is negative");
			if (inf.maxVar == 0) PFLOGE("zero number of variables in header");
			if (inf.maxVar >= INT_MAX - 1) PFLOGE("number of variables not supported");
			inf.nOrgCls = toInteger(str, sign);
			if (sign) PFLOGE("number of cls in header is negative");
			if (inf.nOrgCls == 0) PFLOGE("zero number of cls in header");
			PFLOG2(1, " Found header %d %d", inf.maxVar, inf.nOrgCls);
			inf.nDualVars = V2L(inf.maxVar + 1LL);
			assert(orgs.empty());
			allocSolver();
			initSolver();
		}
		else if (!toClause(in_c, org, str)) return false;
	}
#ifdef __linux__
	if (munmap(buffer, fsz) != 0) PFLOGE("cannot clean input file %s mapping", path.c_str());
	close(fd);
#else
	delete[] buffer;
	inputFile.close();
#endif
	assert(inf.nClauses <= inf.nOrgCls);
	inf.nOrgLits = inf.nLiterals;
	if (inf.nClauses < orgs.size()) orgs.resize(inf.nClauses);
	in_c.clear(true), org.clear(true);
	timer.stop();
	timer.parse = timer.cpuTime();
	PFLOG2(1, " Read %s%d Variables%s, %s%d Clauses%s, and %s%d Literals%s in %s%.2f seconds%s",
		CREPORTVAL, inf.maxVar, CNORMAL,
		CREPORTVAL, orgs.size() + trail.size(), CNORMAL,
		CREPORTVAL, inf.nOrgLits + trail.size(), CNORMAL,
		CREPORTVAL, timer.parse, CNORMAL);
	return true;
}

void ParaFROST::allocSolver()
{
	PFLOGN2(2, " Allocating solver memory for fixed arrays..");
	assert(sizeof(LIT_ST) == 1);
	assert(inf.maxVar);
	uint32 maxSize = inf.maxVar + 1;
	// search space
	sp = new SP(maxSize);
	// input database
	cm.init(inf.nOrgCls);
	orgs.resize(inf.nOrgCls);
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
	assert(ORIGINAL && LEARNT && DELETED);
	assert(UNDEFINED < 0);
	subbin.set_status(ORIGINAL);
	subbin.resize(2);
	dlevels.push(0);
	model.init();
	lbdrest.init(opts.lbd_rate, opts.lbd_fast, opts.lbd_slow);
	if (opts.stable_en && opts.luby_inc) lubyrest.init(opts.luby_inc, opts.luby_max);
	resetSolver();
}

void ParaFROST::resetSolver() {
	PFLOG2(2, " Resetting solver..");
	lrn.lastreduce = -1;
	lrn.var_inc = opts.var_inc;
	lrn.var_decay = opts.var_decay;
	lrn.rounds = opts.mdm_rounds;
	lrn.map_conf_max = opts.map_inc;
	lrn.reduce_conf_max = opts.reduce_inc;
	lrn.rephase_conf_max = opts.rephase_inc;
	lrn.restarts_conf_max = opts.restart_inc;
	lrn.stable_conf_max = opts.stabrestart_inc;
	lrn.mdm_conf_max = int64(scale(opts.mdm_minc));
	lrn.sigma_conf_max = int64(scale(opts.sigma_inc));
	lrn.subsume_conf_max = int64(scale(opts.subsume_inc));
	lrn.stable = opts.stable_en & opts.vsidsonly_en;
	if (lrn.stable) PFLOG2(2, "  VSIDS with initial stable phasing is enabled");
	else PFLOG2(2, "  VMFQ with initial unstable phasing is enabled");
	initVars();
	initQueue();
	initHeap();
	lbdrest.reset();
	PFLOG2(2, " Solver reset successfully");
}

void ParaFROST::solve()
{
	timer.start();
	if (canPreSigmify()) sigmify();
	PFLOG2(2, "-- CDCL search started..");
	if (cnfstate == UNSOLVED) MDMInit();
	while (cnfstate == UNSOLVED && !interrupted()) {
		PFLDL(this, 3);
		if (BCP()) analyze();
		else if (satisfied()) cnfstate = SAT;
		else if (canRestart()) restart();
		else if (canRephase()) rephase();
		else if (canReduce()) reduce();
		else if (canSubsume()) subsume();
		else if (canSigmify()) sigmify();
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
		PFLOG0("\t\t\tSimplifier Report");
		if (!opts.profile_simp)
			PFLOG1(" Simplifier time        : %-10.3f  sec", timer.simp);
		else {
			PFLOG1("  - Var ordering        : %-10.2f  ms", timer.vo);
			PFLOG1("  - CNF compact         : %-10.2f  ms", timer.gc);
			PFLOG1("  - OT  creation        : %-10.2f  ms", timer.cot);
			PFLOG1("  - OT  sorting         : %-10.2f  ms", timer.sot);
			PFLOG1("  - OT  reduction       : %-10.2f  ms", timer.rot);
			PFLOG1("  - BVE                 : %-10.2f  ms", timer.ve);
			PFLOG1("  - HSE                 : %-10.2f  ms", timer.hse);
			PFLOG1("  - BCE                 : %-10.2f  ms", timer.bce);
			PFLOG1("  - ERE                 : %-10.2f  ms", timer.ere);
		}
		PFLOG0("\t\t\tSolver Report");
		PFLOG1(" Solver time            : %-10.3f  sec", timer.solve);
		PFLOG1(" System memory          : %-10.3f  MB", ((double)sysMemUsed() / MBYTE));
		PFLOG1(" Reduces                : %-10lld", stats.reduces);
		PFLOG1(" Rephases               : %-10lld", stats.n_rephs);
		PFLOG1(" Random rephases        : %-10lld", stats.n_randrephs);
		PFLOG1(" Recyclings             : %-10lld", stats.recyclings);
		PFLOG1(" Shrinkages             : %-10d", stats.shrinkages);
		PFLOG1(" Sigmifications         : %-10d", stats.sigmifications);
		PFLOG1(" Mappings               : %-10d", stats.mappings);
		PFLOG1(" Forced units           : %-10lld", stats.n_forced);
		PFLOG1(" Learnt units           : %-10lld", stats.n_units);
		PFLOG1(" Learnt glues           : %-10lld", stats.n_glues);
		PFLOG1(" Learnt subtried        : %-10lld", lrn.subtried);
		PFLOG1(" Learnt subsumed        : %-10lld", stats.n_learntsubs);
		PFLOG1(" Subsume calls          : %-10lld", stats.n_subcalls);
		PFLOG1(" Subsume checks         : %-10lld", stats.n_subchecks);
		PFLOG1(" All subsumed           : %-10lld", stats.n_allsubsumed);
		PFLOG1(" All strengthened       : %-10lld", stats.n_allstrengthened);
		PFLOG1(" Tried redundancies     : %-10lld", stats.n_triedreduns);
		PFLOG1(" Original redundancies  : %-10lld", stats.n_orgreduns);
		PFLOG1(" Learnt redundancies    : %-10lld", stats.n_lrnreduns);
		PFLOG1(" MDM calls              : %-10d", stats.mdm_calls);
		PFLOG1(" Multiple decisions     : %-10lld  (%.1f dec/sec)", stats.n_mds, stats.n_mds / timer.solve);
		PFLOG1(" Follow-Up decisions    : %-10lld  (%.1f dec/sec)", stats.n_fuds, stats.n_fuds / timer.solve);
		PFLOG1(" Propagations           : %-10lld  (%.1f prop/sec)", stats.n_props, stats.n_props / timer.solve);
		PFLOG1(" Chronological     BT   : %-10lld  (%.1f bt/sec)", stats.cbt, stats.cbt / timer.solve);
		PFLOG1(" Non-Chronological BT   : %-10lld  (%.1f bt/sec)", stats.ncbt, stats.ncbt / timer.solve);
		PFLOG1(" Trail reuses           : %-10lld  (%.1f r/sec)", stats.reuses, stats.reuses / timer.solve);
		PFLOG1(" Stable restarts        : %-10lld  (%.1f r/sec)", stats.stab_restarts, stats.stab_restarts / timer.solve);
		PFLOG1(" Restarts               : %-10d  (%.1f r/sec)", starts, starts / timer.solve);
		PFLOG1(" Conflicts              : %-10lld  (%.1f conf/sec)", nConflicts, (nConflicts / timer.solve));
		PFLOG1(" Conflict literals      : %-10lld  (%.2f %% deleted)", stats.tot_lits, (stats.max_lits - stats.tot_lits) * 100.0 / stats.tot_lits);
	}
}

void ParaFROST::wrapup() {
	if (opts.proof_en) proofFile.close();
	if (!quiet_en) { PFLRULER('-', RULELEN); PFLOG0(""); }
	if (cnfstate == SAT) {
		PFLOGS("SATISFIABLE");
		if (opts.model_en) {
			model.extend(sp->value, vorg);
			model.print();
		}
	}
	else if (cnfstate == UNSAT) PFLOGS("UNSATISFIABLE");
	else if (cnfstate == UNSOLVED) PFLOGS("UNKNOWN");
	if (opts.report_en) report();
}