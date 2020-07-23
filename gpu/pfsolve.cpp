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
#include "pfopts.h"

namespace pFROST {

	cudaDeviceProp	devProp;
	uint32			maxGPUThreads = 0;
	CNF_INFO		inf;
	ParaFROST*		pfrost = NULL;
	/***********************************/
	/*     Resources queries           */
	/***********************************/
	int64 getAvailSysMem		()
	{
#ifdef __linux__ 
		struct sysinfo memInfo;
		sysinfo(&memInfo);
		return (memInfo.freeram * (int64)memInfo.mem_unit);
#elif _WIN32
		MEMORYSTATUSEX memInfo;
		memInfo.dwLength = sizeof(MEMORYSTATUSEX);
		GlobalMemoryStatusEx(&memInfo);
		return memInfo.ullAvailPhys;
#endif
	}
	void set_timeout			(int time_limit)
	{
#ifdef __linux__ 
		if (time_limit != 0) {
			rlimit rl;
			getrlimit(RLIMIT_CPU, &rl);
			if (rl.rlim_max == RLIM_INFINITY || (rlim_t)time_limit < rl.rlim_max) {
				rl.rlim_cur = time_limit;
				if (setrlimit(RLIMIT_CPU, &rl) == -1) PFLOGW("timeout cannot be set");
			}
		}
#elif _WIN32
		PFLOGW("timeout not supported on Windows");
#endif
	}
	void handler_terminate		(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%45s", "Interrupted");
			PFLOG0("");
		}
		PFLOGS("UNKNOWN");
		if (!quiet_en) {
			PFLOG0("");
			PFLOGR('-', RULELEN);
		}
		_exit(EXIT_FAILURE);
	}
	void handler_mercy_intr		(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%45s", "Interrupted");
			PFLOG0("");
		}
		pfrost->interrupt();
	}
	void handler_mercy_timeout	(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%45s", "Timeout");
			PFLOG0("");
		}
		pfrost->interrupt();
	}
	void sig_handler			(void h_intr(int), void h_timeout(int))
	{
		signal(SIGINT, h_intr);
		signal(SIGTERM, h_intr);
#ifdef SIGXCPU
		if (h_timeout != NULL) signal(SIGXCPU, h_timeout);
#endif
	}
	//=======================================//
	//		 ParaFROST defined members       //
	//=======================================//
	ParaFROST::ParaFROST(const string& _path) :
		path					(_path)
		, timeout				(opt_timeout)
		, proof_path			(opt_proof_out)
		, seed					(opt_seed)
		, gc_perc				(opt_garbage_perc)
		, parse_only_en			(opt_parseonly_en)
		, report_en				(opt_report_en & !quiet_en)
		, prograte				(opt_progress)
		, vsids_en				(opt_vsids_en)
		, vsidsonly_en			(opt_vsidsonly_en)
		, mdmvsidsonly_en		(opt_mdmvsidsonly_en)
		, mdmfusem_en			(opt_mdmfusem_en)
		, mdmfuses_en			(opt_mdmfuses_en)
		, mcv_en				(!opt_lcv_en)
		, model_en				(opt_model_en)
		, proof_en				(opt_proof_en)
		, guess_en				(opt_guess_en)
		, reduce_en				(opt_reduce_en)
		, rephase_en			(opt_rephase_en)
		, stable_en				(opt_stable_en)
		, target_phase_en		(opt_targetphase_en)
		, sigma_en				(opt_sig_pre_en)
		, sigma_live_en			(opt_sig_live_en)
		, sigma_inc				(opt_sigma_inc)
		, sigma_min				(opt_sigma_min)
		, shrink_min			(opt_shrink_min)
		, map_min				(opt_map_min)
		, var_inc				(opt_var_inc)
		, var_decay				(opt_var_decay)
		, mdm_vsids_pumps		(opt_mdm_vsidspumps)
		, mdm_vmfq_pumps		(opt_mdm_vmfqpumps)
		, mdm_rounds			(opt_mdm_rounds)
		, mdm_freq				(opt_mdm_freq)
		, mdm_minc				(opt_mdm_minc)
		, mdm_sinc				(opt_mdm_sinc)
		, mdm_div				(opt_mdm_div)
		, stabrestart_r			(opt_stabrestart_r)
		, stabrestart_inc		(opt_stabrestart_inc)
		, luby_inc				(opt_luby_inc)
		, luby_max				(opt_luby_max)
		, restart_base			(opt_rest_base)
		, restart_inc			(opt_rest_inc)
		, polarity				(opt_polarity)
		, rephase_inc			(opt_rephase_inc)
		, lbd_tier1				(opt_lbd_tier1)
		, lbd_tier2				(opt_lbd_tier2)
		, lbd_reduce_min		(opt_lbd_min)
		, lbd_csize_min			(opt_lbd_min_size)
		, reduce_perc			(opt_reduce_perc)
		, reduce_inc			(opt_reduce_inc)
		, wt					(cm)
		, wtBin					(cm)
		, vHeap					(HEAP_CMP(varAct))
		, starts				(1)
		, nForced				(0)
		, nConflicts			(0)
		, intr					(false)
		, mapped				(false)
		, search_guess			(false)
		, conflict				(NOREF)
		, cnfstate				(UNSOLVED)
		, sigState				(AWAKEN_SUCC)
	{
		stats.sysMemAvail = getAvailSysMem();
		CHECK(cudaGetDeviceCount(&devCount));
		size_t _gfree = 0, memPenalty = 0;
		if (devCount == 0) PFLOGE("no GPU(s) available that support CUDA");
		else {
			PFLOG2(1, " Detected (%d) CUDA-enabled GPU(s)", devCount);
			CHECK(cudaGetDeviceProperties(&devProp, MASTER_GPU));
			assert(devProp.totalGlobalMem);
			if (devProp.warpSize != 32) PFLOGE("GPU warp size not supported");
#ifdef __linux__ 
			memPenalty = 200 * MBYTE;
			_gfree = devProp.totalGlobalMem - memPenalty;
#elif _WIN32
			memPenalty = 800 * MBYTE;
			_gfree = devProp.totalGlobalMem - memPenalty;
#endif
			maxGPUThreads = devProp.multiProcessorCount * devProp.maxThreadsPerMultiProcessor;
			cumem.init(_gfree);
			hcnf = NULL, cnf = NULL, ot = NULL, vars = NULL, streams = NULL;
		}
		if (!quiet_en) {
			PFLOG1(" Free system memory = %lld GB", stats.sysMemAvail / GBYTE);
			PFLOG1(" GPU: \"%s\" (compute cap: %d.%d)", devProp.name, devProp.major, devProp.minor);
			PFLOG1("  - Multiprocessors = %d MPs (%d cores/MP)", devProp.multiProcessorCount, SM2Cores(devProp.major, devProp.minor));
			PFLOG1("  - Free Global memory = %zd GB", _gfree / GBYTE);
			PFLOG1("  - Free Shared memory = %zd KB", devProp.sharedMemPerBlock / KBYTE);
			PFLOGR('-', RULELEN);
		}
		if (_gfree <= 200 * MBYTE) {
			PFLOGW("not enough GPU memory (free = %zd MB) -> skip SIGmA", _gfree / MBYTE);
			sigma_en = sigma_live_en = false;
		}
		if (sigma_en || sigma_live_en) optSimp(), createStreams();
		if (proof_en) {
			proofFile.open(proof_path, std::ofstream::binary | std::ofstream::out);
			if (!proofFile.is_open()) PFLOGE("cannot open proof file %s", proof_path.c_str());
		}
		if (parser() == UNSAT || BCP() != NOREF) { if (proof_en) wrProof('0'); cnfstate = UNSAT, killSolver(); }
		if (parse_only_en) killSolver();
		shrinkTop(), recycle();
		if (verbose == 1) printTable();
	}

	CNF_ST ParaFROST::parser() {
		struct stat st;
		stat(path.c_str(), &st);
		size_t fsz = st.st_size;
		PFLOG2(1, " Parsing CNF file \"%s\" (size: %zd KB)", path.c_str(), fsz / KBYTE);
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
		Lits_t in_c, tmpCl;
		in_c.reserve(INIT_CAP);
		PFLMEMCALL(this, 2);
		int oldSz = 0;
		char* eof = str + fsz;
		while (str < eof) {
			eatWS(str);
			if (*str == '\0' || *str == '0' || *str == '%') break;
			if (*str == 'c') eatLine(str);
			else if (*str == 'p') {
				if (!eq(str, "p cnf")) PFLOGE("header has wrong format");
				uint32 sign = 0;
				inf.maxVar = toInteger(str, sign);
				if (sign) PFLOGE("number of variables in header is negative");
				if (inf.maxVar == 0) PFLOGE("zero number of variables in header");
				if (inf.maxVar >= INT_MAX - 1) PFLOGE("number of variables not supported");
				inf.nOrgCls = toInteger(str, sign);
				if (sign) PFLOGE("number of cls in header is negative");
				if (inf.nOrgCls == 0) PFLOGE("zero number of cls in header");
				PFLOG2(1, " Found header %d %d", inf.maxVar, inf.nOrgCls);
				inf.nDualVars = v2l(inf.maxVar + 1);
				assert(orgs.empty());
				allocVars();
				initSolver();
			}
			else {
				toClause(in_c, str);
				if (proof_en) {
					oldSz = in_c.size();
					tmpCl.clear(true), tmpCl.resize(in_c.size()), tmpCl.copyFrom(in_c);
				}
				if (checkClause(in_c)) { // clause not a tautology
					if (proof_en && in_c.size() < oldSz) {
						wrProof('a'), wrProof(in_c, in_c.size()), wrProof(0);
						wrProof('d'), wrProof(tmpCl, in_c.size()), wrProof(0);
					}
					if (in_c.size() == 1) {
						LIT_ST val = value(*in_c);
						if (val == UNDEFINED) enqueue(*in_c);
						else if (!val) return UNSAT;
					}
					else if (inf.nClauses + 1 > inf.nOrgCls) PFLOGE("too many clauses");
					else newClause(in_c);
				}
			}
		}
#ifdef __linux__
		if (munmap(buffer, fsz) != 0) PFLOGE("cannot clean input file %s mapping", path.c_str());
		close(fd);
#else
		delete[] buffer;
		inputFile.close();
#endif
		assert(inf.nClauses + inf.nOrgBins <= inf.nOrgCls);
		inf.nOrgLits = inf.nLiterals + (inf.nOrgBins << 1);
		if (inf.nClauses < orgs.size()) orgs.resize(inf.nClauses);
		in_c.clear(true), tmpCl.clear(true);
		timer.stop();
		timer.parse = timer.cpuTime();
		PFLOG2(1, " Read %d Vars, %d Cls, and %d Lits in %.2f seconds", inf.maxVar, inf.nClauses + inf.nOrgBins + trail.size(), inf.nOrgLits + trail.size(), timer.parse);
		return UNSOLVED;
	}

	int64 ParaFROST::sysMemUsed()
	{
		int64 memUsed = 0;
#ifdef __linux__ 
		FILE* file = fopen("/proc/self/status", "r");
		char line[128];
		uint32 sign = 0;
		while (fgets(line, 128, file) != NULL) {
			char* str = line;
			if (eq(str, "VmRSS:")) {
				eatWS(str);
				memUsed = toInteger(str, sign);
				break;
			}
		}
		fclose(file);
		return memUsed * KBYTE;
#elif _WIN32
		PROCESS_MEMORY_COUNTERS_EX memInfo;
		GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memInfo, sizeof(PROCESS_MEMORY_COUNTERS_EX));
		memUsed = memInfo.WorkingSetSize;
#endif
		return memUsed;
		}

	void ParaFROST::killSolver()
	{
		wrapup();
		this->~ParaFROST();
		pfrost = NULL;
		exit(EXIT_SUCCESS);
	}

	void ParaFROST::allocVars()
	{
		PFLOGN2(2, " Allocating solver memory for fixed arrays..");
		assert(sizeof(C_REF) == sizeof(uint32));
		assert(sizeof(LIT_ST) == 1);
		assert(inf.maxVar);
		uint32 maxSize = inf.maxVar + 1;
		// search space
		sp = new SP(maxSize);
		// input database
		cm.init(inf.nOrgCls);
		orgs.reserve(inf.nOrgCls);
		// watch tables
		wt.resize(inf.nDualVars);
		wtBin.resize(inf.nDualVars);
		// others
		toSimp.reserve(inf.maxVar >> 2);
		varAct.resize(maxSize, 0.0);
		varBumps.resize(maxSize, 0);
		PFLDONE(2, 5);
		PFLMEMCALL(this, 2);
	}

	void ParaFROST::initSolver()
	{
		assert(ORIGINAL && LEARNT && DELETED);
		assert(UNDEFINED < 0 && ROOT_LEVEL == 0);
		srand(seed);
		model.init();
		lbdrest.init(opt_lbd_rate, opt_lbdfast, opt_lbdslow);
		if (stable_en && luby_inc) lubyrest.init(luby_inc, luby_max);
		resetSolver();
	}

	void ParaFROST::resetSolver(const bool& rst) {
		PFLOG2(2, " Resetting solver..");
		lrn.bumped = 0;
		lrn.numMDs = 0;
		lrn.nRefVars = 0;
		lrn.var_inc = var_inc;
		lrn.var_decay = var_decay;
		lrn.target = lrn.best = 0;
		lrn.lastsimplified = 0;
		lrn.rephased[0] = lrn.rephased[1] = 0;
		lrn.rounds = mdm_rounds;
		lrn.mdm_conf_max = mdm_minc;
		lrn.sigma_conf_max = sigma_inc;
		lrn.reduce_conf_max = reduce_inc;
		lrn.rephase_conf_max = rephase_inc;
		lrn.restarts_conf_max = restart_base;
		lrn.stable_conf_max = stabrestart_inc;
		lrn.stable = stable_en & (vsidsonly_en || mdm_rounds);
		if (lrn.stable) PFLOG2(2, "  VSIDS with stable phasing is enabled");
		else PFLOG2(2, "  VMFQ with initial unstable phasing is enabled");
		if (rst) {
			size_t maxSize = (size_t)inf.maxVar + 1;
			memset(sp->value, UNDEFINED, inf.nDualVars);
			memset(sp->ptarget, UNDEFINED, maxSize);
			memset(sp->pbest, UNDEFINED, maxSize);
			memset(sp->psaved, polarity, maxSize);
			for (uint32 v = 1; v <= inf.maxVar; v++)
				sp->level[v] = UNDEFINED, sp->source[v] = NOREF;
		}
		initQueue();
		initHeap();
		lbdrest.reset();
		stats.reset();
		PFLOG2(2, " Solver reset successfully");
	}

	void ParaFROST::newClause(Lits_t& in_c, const CL_ST& type)
	{
		C_REF r = cm.alloc(in_c);
		int sz = in_c.size();
		assert(sz > 1);
		assert(sz == cm[r].size());
		assert(cm[r][0] > 1 && cm[r][1] > 1);
		assert(cm[r][0] <= UINT32_MAX && cm[r][1] <= UINT32_MAX);
		if (sz == 2) {
			wtBin.newWatch(r, in_c[0], in_c[1]);
			if (type == ORIGINAL) inf.nOrgBins++;
			else inf.nLearntBins++;
		}
		else {
			assert(sz > 2);
			wt.newWatch(r, in_c[0], in_c[1]);
			if (type == ORIGINAL) {
				orgs.push(r);
				inf.nClauses++;
				inf.nLiterals += sz;
			}
			else {
				assert(sp->learnt_lbd > 0);
				CLAUSE& c = cm[r];
				if (sp->learnt_lbd > sz) c.set_lbd(sz);
				else c.set_lbd(sp->learnt_lbd);
				if (sp->learnt_lbd <= lbd_tier1) inf.nGlues++;			// Tier1
				else if (sp->learnt_lbd <= lbd_tier2) c.initTier2();	// Tier2
				else c.initTier3();										// Tier3
				learnts.push(r);
				inf.nLearntLits += sz;
			}
		}
		if (type == ORIGINAL) cm[r].set_status(ORIGINAL);
		else assert(type == LEARNT), cm[r].set_status(LEARNT), enqueue(*learntC, sp->bt_level, r);
	}

	void ParaFROST::removeClause(C_REF& r, const bool& gc) {
		CLAUSE& c = cm[r];
		assert(c.size() > 1);
		if (c.size() == 2) {
			if (c.original()) inf.nOrgBins--;
			else assert(c.learnt()), inf.nLearntBins--;
		}
		else {
			if (c.original()) inf.nLiterals -= c.size();
			else assert(c.learnt()), inf.nLearntLits -= c.size();
			if (proof_en) {
				wrProof('d');
				wrProof(c, c.size());
				wrProof(0);
			}
		}
		assert(!c.deleted());
		detachClause(r, gc), c.markDeleted(), cm.collect(r);
		depleteSource(c);
	}

	void ParaFROST::shrinkClause(C_REF& r)
	{
		CLAUSE& c = cm[r];
		assert(c.size() > 2);
		assert(unassigned(c[0]) && unassigned(c[1]));
		int sz = c.size();
		Lits_t tmpCl;
		if (proof_en) { tmpCl.resize(sz); tmpCl.copyFrom(c, sz); }
		int k = 2;
		while (k < sz) if (isFalse(c[k])) c[k] = c[--sz]; else k++;
		int remLits = c.size() - sz;
		if (remLits && sz == 2) { // clause shrinked to binary
			wt.remWatch(flip(c[0]), r);
			wt.remWatch(flip(c[1]), r);
			wtBin.newWatch(r, c[0], c[1]);
			if (c.original()) inf.nOrgBins++;
			else assert(c.learnt()), inf.nLearntBins++;
		}
		if (proof_en && remLits) {
			wrProof('a');
			wrProof(c, sz);
			wrProof(0);
			wrProof('d');
			wrProof(tmpCl, c.size());
			wrProof(0);
		}
		c.shrink(remLits); // adjusts "pos" also
		cm.collect(remLits);
		assert(c.size() > 1);
		if (c.learnt() && c.size() > 2) bumpShrinked(c);
	}

	void ParaFROST::shrinkWT()
	{
		for (uint32 i = sp->simplified; i < trail.size(); i++) {
			uint32 assign = trail[i], assign_f = flip(assign);
			wt.collect(assign), wt.collect(assign_f);
			WL& ws = wtBin[assign_f];
			for (WATCH* w = ws; w != ws.end(); w++) {
				assert(!cm[w->ref].moved());
				if (!cm[w->ref].deleted())
					removeClause(w->ref);
			}
		}
	}

	void ParaFROST::shrinkTop()
	{
		if (trail.size() == sp->simplified) return;
		PFLOGN2(2, " Shrinking CNF before solve..");
		assert(trail.size() && DL() == ROOT_LEVEL);
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		assert(sp->simplified == 0);
		assert(sp->propagated == trail.size());
		assert(!unassigned(trail.back()));
		shrinkWT();
		shrink(orgs);
		inf.nClauses = orgs.size();
		int simpVars = trail.size() - sp->simplified;
		lrn.lastsimplified += simpVars;
		PFLENDING(2, 5, "(-%d variables)", simpVars);
		sp->simplified = trail.size();
		stats.shrinkages++;
	}

	void ParaFROST::shrink()
	{
		// remove satisfied clauses & falsified literals
		PFLOGN2(2, " Shrinking CNF..");
		assert(trail.size() && DL() == ROOT_LEVEL);
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		assert(sp->propagated == trail.size());
		assert(!unassigned(trail.back()));
		shrinkWT();
		shrink(orgs);
		shrink(learnts);
		inf.nClauses = orgs.size();
		int simpVars = trail.size() - sp->simplified;
		lrn.lastsimplified += simpVars;
		PFLENDING(2, 5, "(-%d variables)", simpVars);
		sp->simplified = trail.size();
		stats.shrinkages++;
	}

	void ParaFROST::reduceLearnts() {
		assert(learnts.size());
		reduced.reserve(learnts.size());
		// keep recently used, tier1, and reason learnts otherwise reduce
		for (C_REF* r = learnts; r != learnts.end(); r++) {
			CLAUSE& c = cm[*r];
			assert(c.size() > 2);
			if (c.reason()) continue;
			if (c.usage()) c.warm();
			else if (c.lbd() > lbd_tier1) reduced.push(*r);
		}
		if (reduced.size()) {
			uint32 pivot = reduce_perc * reduced.size();
			PFLOGN2(2, " Reducing learnt database (up to %d clauses)..", pivot);
			std::stable_sort(reduced.data(), reduced.end(), LEARNT_CMP(cm));
			// remove learnts from database
			C_REF* h = reduced, * e = h + pivot;
			while (h != e) {
				assert(!cm[*h].reason());
				assert(cm[*h].lbd() > lbd_tier1);
				removeClause(*h++);
			}
			uint32 n = 0;
			for (C_REF* r = learnts; r != learnts.end(); r++)
				if (!cm[*r].deleted()) learnts[n++] = *r;
			assert(n == learnts.size() - pivot);
			learnts.resize(n);
			PFLDONE(2, 5);
		}
		reduced.clear();
	}

	void ParaFROST::reduce()
	{
		// remove satisfied clauses if possible then reduce learnts
		assert(sp->propagated == trail.size());
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		stats.reduces++;
		bool shrinked = false;
		if (canShrink()) shrink(), shrinked = true;
		reduceLearnts();
		// update counters
		double current_inc = reduce_inc * (stats.reduces + 1);
		reduceWeight(current_inc);
		lrn.reduce_conf_max = nConflicts + current_inc;
		PFLOG2(2, " Max reduce limit increased to %lld conflicts by a weight %.2f", lrn.reduce_conf_max, current_inc);
		recycle();
		if (shrinked && canMap()) map();
	}

	void ParaFROST::selfsubBin()
	{
		assert(sp->learnt_lbd != 0);
		stats.marker++;
		for (int i = 1; i < learntC.size(); i++) // mark all learnt literals except asserting
			sp->board[l2a(learntC[i])] = stats.marker;
		// check binary clauses if subsuming learnt clause
		uint32 parent = flip(learntC[0]);
		WL& wBins = wtBin[parent];
		int nLitsRem = 0;
		for (int i = 0; i < wBins.size(); i++) {
			uint32 v = l2a(wBins[i].imp);
			if (sp->board[v] == stats.marker && isTrue(wBins[i].imp)) {
				sp->board[v] = stats.marker - 1; // unmark
				nLitsRem++;
			}
		}
		int last = learntC.size() - 1;
		if (nLitsRem > 0) {
			for (int i = 1; i < learntC.size() - nLitsRem; i++) {
				if (sp->board[l2a(learntC[i])] != stats.marker) {
					swap(learntC[last], learntC[i]);
					last--, i--;
				}
			}
			learntC.shrink(nLitsRem);
			sp->learnt_lbd = calcLBD(learntC); // update learnt LBD
			PFLLEARNT(this, 3);
		}
	}

	bool ParaFROST::selfsub(const uint32& learntLit, uint32*& tail)
	{
		assert(minLevel);
		assert(l2r(learntLit) != NOREF);
		assert(l2dl(learntLit) > ROOT_LEVEL);
		uint32* tmp_tail = tail;
		toSimp.clear();
		toSimp.push(learntLit);
		for (uint32 n = 0; n < toSimp.size(); n++) {
			uint32 next = toSimp[n];
			assert(l2r(next) != NOREF);
			CLAUSE& c = cm[l2r(next)];
			if (c.size() == 2 && isFalse(*c)) {
				assert(isTrue(c[1]));
				c.swapWatched();
			}
			assert(*c == flip(next));
			for (int j = 1; j < c.size(); j++) {
				uint32 lit = c[j], v = l2a(lit);
				int litLevel = sp->level[v];
				LIT_ST& seen = sp->seen[v];
				if (litLevel == ROOT_LEVEL || seen) continue;
				if (sp->source[v] == NOREF || abortSub(litLevel)) {
					while (tail != tmp_tail) sp->seen[l2a(*--tail)] = 0;
					return false;
				}
				seen = 1;
				*tail++ = lit;
				toSimp.push(lit);
			}
		}
		return true;
	}

	void ParaFROST::selfsub()
	{
		assert(learntC.size() > 1);
		register uint32* i = sp->tmp_stack, * tail = i + learntC.size(), * j = learntC, * end = learntC.end();
		while (j != end) *i++ = *j++;
		j = learntC + 1, minLevel = 0;
		while (j != end)  minLevel |= l2hl(*j++);
		j = learntC + 1, i = j;
		for (; i != end; i++)
			if (l2r(*i) == NOREF || !selfsub(*i, tail))	*j++ = *i;
		int newSize = j - learntC;
		learntC.resize(newSize);
		PFLLEARNT(this, 3);
		// minimize further using binaries
		if (newSize > 1 && newSize <= lbd_csize_min && ((sp->learnt_lbd = calcLBD(learntC)) <= lbd_reduce_min))
			selfsubBin();
		// clear seen
		i = sp->tmp_stack;
		while (i != tail) sp->seen[l2a(*i++)] = 0;
	}

	void ParaFROST::whereToJump()
	{
		if (learntC.size() == 1) sp->bt_level = ROOT_LEVEL;
		else if (learntC.size() == 2) sp->bt_level = l2dl(learntC[1]);
		else {
			int max_k = 1;
			for (int k = 2; k < learntC.size(); k++) {
				if (l2dl(learntC[k]) > l2dl(learntC[max_k]))
					max_k = k;
			}
			uint32 max_k_lit = learntC[max_k];
			learntC[max_k] = learntC[1];
			learntC[1] = max_k_lit;
			sp->bt_level = l2dl(max_k_lit);
		}
	}

	void ParaFROST::cancelAssigns(const int& bt_level)
	{
		int level = DL();
		if (level == bt_level) return;
		savePhases(bt_level);
		uint32 from = trail_lens[bt_level], i = from;
		while (i < trail.size()) {
			uint32 lit = trail[i++];
			if (sp->level[l2a(lit)] > bt_level) cancelAssign(lit);
		}
		PFLOG2(3, " %d literals kept and %d are cancelled", from, trail.size() - from);
		trail.resize(from);
		sp->propagated = trail.size();
		assert(trail_lens[bt_level] == trail.size());
		trail_lens.shrink(level - bt_level);
		assert(DL() == bt_level);
	}

	void ParaFROST::backJump(const int& bt_level) {
		assert(trail.size() > 0);
		PFLOG2(3, " Backjumping to level %d, from literal %d at trail index %d",
			bt_level, l2i(trail[trail_lens[bt_level]]), trail_lens[bt_level]);
		// cancel old assignments up to backtrack level <bt_level>
		cancelAssigns(bt_level);
		// add learnt clause & enqueue learnt decision
		if (proof_en) {
			wrProof('a');
			wrProof(learntC.data(), learntC.size());
			wrProof(0);
		}
		if (learntC.size() == 1) enqueue(learntC[0]), stats.n_units++;
		else newClause(learntC, LEARNT);
	}

	void ParaFROST::analyze()
	{
		assert(conflict != NOREF);
		assert(analyzed.empty());
		PFLOG2(3, " Analyzing conflict:");
		PFLTRAIL(this, 3);
		nConflicts++;
		if (DL() == ROOT_LEVEL) { cnfstate = UNSAT; return; }
		learntC.clear();
		learntC.push(0);
		sp->learnt_lbd = 0;
		uint32 parent = 0;
		int track = 0, index = trail.size() - 1;
		C_REF r = conflict;
		do {
			assert(r != NOREF);
			analyzeReason(r, parent, track);
			// next implication clause
			while (!sp->seen[l2a(trail[index--])]);
			parent = trail[index + 1];
			assert(parent > 0);
			r = l2r(parent);
			sp->seen[l2a(parent)] = 0;
		} while (--track > 0);
		assert(learntC[0] == 0);
		learntC[0] = flip(parent);
		PFLLEARNT(this, 3);
		// minimize learnt clause 
		stats.max_lits += learntC.size();
		if (learntC.size() > 1) {
			selfsub();
			if (!sp->learnt_lbd) sp->learnt_lbd = calcLBD(learntC); // calculate lbd, if clause not minimized
		}
		else sp->learnt_lbd = 1;
		assert(sp->learnt_lbd != 0);
		PFLOG2(4, " lbd of learnt clause = %d", sp->learnt_lbd);
		stats.tot_lits += learntC.size();
		// adjust restart rate
		lbdrest.update(sp->learnt_lbd);
		if (lrn.stable) lubyrest.update();
		whereToJump();
		bumpVariables();
		backJump(sp->bt_level);
		stats.ncbt++;
		printStats(vsidsOnly() && nConflicts % prograte == 0);
	}

	C_REF ParaFROST::BCP()
	{
		conflict = NOREF;
		uint32 propsBefore = sp->propagated;
		while (sp->propagated < trail.size()) {
			uint32 assign = trail[sp->propagated++], assign_dl = DL();
			assert(assign > 0);
			PFLBCPS(this, 4, assign);
			// binaries
			WL& wBins = wtBin.getClean(assign);
			for (int i = 0; i < wBins.size(); i++) {
				uint32 imp = wBins[i].imp;
				LIT_ST val = value(imp);
				if (!val) { conflict = wBins[i].ref; goto bcp_exit; }
				if (val < 0) enqueue(imp, assign_dl, wBins[i].ref);
			}
			// non-binaries
			WL& ws = wt.getClean(assign);
			if (ws.size()) {
				WATCH* w_i = ws, * w_j = w_i, * w_end = ws.end();
				while (w_i != w_end) {
					uint32 wImp = w_i->imp;
					assert(wImp && wImp < UINT32_MAX);
					if (isTrue(wImp)) { *w_j++ = *w_i++; continue; }
					C_REF r = w_i->ref;
					CLAUSE& c = cm[r];
					// move assigned-0 literal to watch 1
					assert(c.size() >= 2);
					uint32 f_assign = flip(assign);
					register uint32 wlit = c[0] ^ c[1] ^ f_assign; // Thanks to Cadical solver
					c[0] = wlit, c[1] = f_assign;
					assert(c[1] == f_assign);
					w_i++;
					// check if first literal is true
					LIT_ST val0 = value(c[0]);
					WATCH w = WATCH(r, c[0]);
					if (c[0] != wImp && val0 > 0) *w_j++ = w;
					else {
						// look for (un)-assigned-1 literal to watch
						uint32* cmid = c.mid(), * cend = c.end();
						uint32* k = cmid, lit = 0;
						LIT_ST _false_ = UNDEFINED;
						while (k != cend && (_false_ = isFalse(lit = *k))) k++;
						assert(_false_ != UNDEFINED);
						if (_false_) {
							k = c + 2;
							assert(c.pos() < c.size());
							while (k != cmid && (_false_ = isFalse(lit = *k))) k++;
						}
						assert(k >= c + 2 && k <= c.end());
						c.set_pos(int(k - c)); // set new position
						if (!_false_) c[1] = lit, * k = f_assign, wt[flip(lit)].push(w); // prepare new watch
						else { // clause is unit or conflict
							*w_j++ = w;
							if (val0 < 0) {
								assert(l2dl(c[0]) == UNDEFINED);
								enqueue(c[0], assign_dl, r);
							}
							else {
								PFLCONFLICT(this, 3, c[0]);
								sp->propagated = trail.size();
								conflict = r;
								while (w_i < w_end) *w_j++ = *w_i++;
							}
						}
					}
				}
				ws.shrink(int(w_i - w_j));
			}
			PFLBCPE(this, 4, assign);
		}
	bcp_exit:
		if (!search_guess) stats.n_props += (sp->propagated - propsBefore);
		return conflict;
	}

	void ParaFROST::solve()
	{
		timer.start();
		if (canPreSigmify()) sigmify();
		if (guess_en) guess();
		if (cnfstate == UNSOLVED) MDMInit();
		PFLOG2(2, "-- CDCL search started..");
		while (cnfstate == UNSOLVED && !interrupted()) {
			PFLDL(this, 3);
			if (BCP() != NOREF) analyze();
			else if (satisfied()) cnfstate = SAT;
			else if (canRestart()) restart();
			else if (canRephase()) rephase();
			else if (canReduce()) reduce();
			else if (canSigmify()) sigmify();
			else if (canMMD()) MDM();
			else decide();
			PFLTRAIL(this, 3);
		}
		timer.stop(), timer.solve += timer.cpuTime();
		PFLOG2(2, "-- CDCL search completed successfully");
		wrapup();
	}

	void ParaFROST::restart()
	{
		assert(sp->propagated == trail.size());
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		if (lrn.stable) stats.stab_restarts++;
		starts++;
		cancelAssigns();
		sp->bt_level = ROOT_LEVEL;
		if (mdmfuses_en) MDMFuseSlave();
		lrn.restarts_conf_max = nConflicts + restart_inc;
		PFLOG2(3, " new restart limit after %lld conflicts", lrn.restarts_conf_max);
	}

	bool ParaFROST::vibrate() {
		if (!stable_en) return false;
		if (vsidsOnly()) return true;
		if (--lrn.stable_conf_max == 0) {
			lrn.stable = !lrn.stable;
			stabrestart_inc *= stabrestart_r;
			lrn.stable_conf_max = stabrestart_inc;
			PFLOG2(2, " Max stable conflicts increased to %lld", lrn.stable_conf_max);
			lbdrest.swap();
			printStats();
		}
		return lrn.stable;
	}

	void ParaFROST::recycle(CMM& new_cm)
	{
		wtBin.recycle(), wt.recycle();
		for (uint32 v = 1; v <= inf.maxVar; v++) {
			uint32 lit = what(v), flit = flip(lit);
			recycleWL(wt[lit], new_cm), recycleWL(wtBin[lit], new_cm);
			recycleWL(wt[flit], new_cm), recycleWL(wtBin[flit], new_cm);
		}
		for (uint32 i = 0; i < trail.size(); i++) {
			uint32 v = l2a(trail[i]);
			C_REF& r = sp->source[v];
			if (r != NOREF) {
				assert(!cm[r].deleted());
				if (cm[r].moved()) r = cm[r].ref();
			}
		}
		for (uint32 i = 0; i < learnts.size(); i++) cm.move(learnts[i], new_cm);
		for (uint32 i = 0; i < orgs.size(); i++) cm.move(orgs[i], new_cm);
	}

	void ParaFROST::recycle() {
		assert(sp->propagated == trail.size());
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		if (cm.garbage() > cm.size() * gc_perc) {
			PFLOGN2(2, " Recycling garbage..");
			CMM new_cm(cm.size() - cm.garbage());
			recycle(new_cm);
			PFLGCMEM(2, cm, new_cm);
			new_cm.migrate(cm);
			PFLDONE(2, 5);
			stats.recyclings++;
		}
	}

	void ParaFROST::map(WT& wt) {
		if (wt.empty()) return;
		assert(!vmap.empty());
		for (uint32 v = 1; v <= inf.maxVar; v++) {
			uint32 mVar = vmap.mapped(v);
			if (mVar) {
				uint32 p = v2l(v), n = neg(p);
				uint32 mpos = v2l(mVar), mneg = neg(mpos);
				if (mVar != v) { // map watch lists
					wt[mpos].copyFrom(wt[p]);
					wt[mneg].copyFrom(wt[n]);
				}
				map(wt[mpos]), map(wt[mneg]); // then map watch imps
			}
		}
		wt.resize(v2l(vmap.size()));
		wt.shrinkCap();
	}

	void ParaFROST::mapBins() {
		uint32 nBins = inf.nOrgBins + inf.nLearntBins;
		if (nBins == 0) return;
		BCNF bins; // for caching
		bins.reserve(nBins);
		extractBins(bins);
		mapBins(bins);
	}

	void ParaFROST::map(const bool& sigmified)
	{
		assert(!satisfied());
		assert(conflict == NOREF);
		assert(DL() == ROOT_LEVEL);
		assert(trail.size() == sp->propagated);
		stats.mappings++;
		int64 memBefore = sysMemUsed();
		// create vmap
		vmap.initiate(sp);
		// map original literals with current values
		vmap.mapOrgs(model.lits);
		// map clauses and watch tables
		if (!sigmified) {
			wtBin.recycle(), wt.recycle();
			mapBins(), map(orgs), map(learnts);
			map(wtBin), map(wt);
		}
		else mapped = true, newBeginning(), mapped = false;
		// map trail, queue and heap
		vmap.mapShrinkLits(trail);
		sp->propagated = trail.size();
		vQueue.map(*vmap, vmap.firstL0());
		vmap.mapShrinkVars(vQueue.data());
		vmap.mapShrinkVars(varBumps);
		uVec1D tmp;
		while (vHeap.size()) {
			uint32 x = vHeap.pop();
			if (x == vmap.firstL0()) continue;
			uint32 mx = vmap.mapped(x);
			if (mx) tmp.push(mx);
		}
		vmap.mapShrinkVars(varAct);
		vHeap.rebuild(tmp);
		// map phases & variable states
		SP* newSP = new SP(vmap.size());
		vmap.mapSP(newSP);
		delete sp;
		sp = newSP;
		// update counters
		lrn.best = lrn.target = 0;
		for (uint32 v = 1; v <= vmap.numVars(); v++) {
			if (sp->pbest[v] != UNDEFINED) lrn.best++;
			if (sp->ptarget[v] != UNDEFINED) lrn.target++;
		}
		PFLOG2(2, " Variable mapping compressed %d to %d, saving %.2f KB of memory",
			inf.maxVar, vmap.numVars(), double(abs(memBefore - sysMemUsed())) / KBYTE);
		inf.maxVar = vmap.numVars();
		inf.nDualVars = v2l(inf.maxVar + 1);
		inf.maxFrozen = inf.maxMelted = 0;
		vmap.destroy();
	}

	void ParaFROST::printReport()
	{
		if (report_en) {
			PFLOG0("");
			PFLOG0("\t\t\tSimplifier Report");
			PFLOG1(" Simplifier time      : %-10.3f  sec", timer.simp);
			if (profile_gpu) {
				PFLOG1("  - Var ordering      : %-10.2f  ms", cutimer->vo);
				PFLOG1("  - CNF signatures    : %-10.2f  ms", cutimer->sig);
				PFLOG1("  - CNF compact       : %-10.2f  ms", cutimer->gc);
				PFLOG1("  - CNF transfer      : %-10.2f  ms", cutimer->io);
				PFLOG1("  - OT  creation      : %-10.2f  ms", cutimer->cot);
				PFLOG1("  - OT  sorting       : %-10.2f  ms", cutimer->sot);
				PFLOG1("  - OT  reduction     : %-10.2f  ms", cutimer->rot);
				PFLOG1("  - BVE               : %-10.2f  ms", cutimer->ve);
				PFLOG1("  - HSE               : %-10.2f  ms", cutimer->hse);
				PFLOG1("  - BCE               : %-10.2f  ms", cutimer->bce);
				PFLOG1("  - HRE               : %-10.2f  ms", cutimer->hre);
			}
			PFLOG1(" Device memory        : %-10.3f  MB", ((double)cumem.maxCapacity() / MBYTE));
			PFLOG0("\t\t\tSolver Report");
			PFLOG1(" Solver time          : %-10.3f  sec", timer.solve);
			PFLOG1(" System memory        : %-10.3f  MB", ((double)sysMemUsed() / MBYTE));
			PFLOG1(" Guessed              : %-10s  (%s)", stats.guess_succ ? "yes" : "no", guess_en ? stats.guess_who : "disabled");
			PFLOG1(" Reduces              : %-10lld", stats.reduces);
			PFLOG1(" Rephases             : %-10lld", stats.n_rephs);
			PFLOG1(" Recyclings           : %-10lld", stats.recyclings);
			PFLOG1(" Shrinkages           : %-10d", stats.shrinkages);
			PFLOG1(" Sigmifications       : %-10d", stats.sigmifications);
			PFLOG1(" Mappings             : %-10d", stats.mappings);
			PFLOG1(" Forced units         : %-10lld", stats.n_forced);
			PFLOG1(" Learnt units         : %-10lld", stats.n_units);
			PFLOG1(" Learnt binaries      : %-10d", inf.nLearntBins);
			PFLOG1(" Learnt glues         : %-10d", inf.nLearntBins + inf.nGlues);
			PFLOG1(" MDM calls            : %-10d", stats.mdm_calls);
			PFLOG1(" Multiple decisions   : %-10lld  (%.1f dec/sec)", stats.n_mds, stats.n_mds / timer.solve);
			PFLOG1(" Follow-Up decisions  : %-10lld  (%.1f dec/sec)", stats.n_fuds, stats.n_fuds / timer.solve);
			PFLOG1(" Propagations         : %-10lld  (%.1f prop/sec)", stats.n_props, stats.n_props / timer.solve);
			PFLOG1(" Non-chronological BT : %-10lld  (%.1f bt/sec)", stats.ncbt, stats.ncbt / timer.solve);
			PFLOG1(" Stable restarts      : %-10lld  (%.1f r/sec)", stats.stab_restarts, stats.stab_restarts / timer.solve);
			PFLOG1(" Restarts             : %-10d  (%.1f r/sec)", starts, starts / timer.solve);
			PFLOG1(" Conflicts            : %-10lld  (%.1f conf/sec)", nConflicts, (nConflicts / timer.solve));
			PFLOG1(" Conflict literals    : %-10lld  (%.2f %% deleted)", stats.tot_lits, (stats.max_lits - stats.tot_lits) * 100.0 / stats.tot_lits);
		}
	}

	void ParaFROST::wrapup() {
		if (proof_en) proofFile.close();
		if (!quiet_en) { PFLOGR('-', RULELEN); PFLOG0(""); }
		if (cnfstate == SAT) {
			PFLOGS("SATISFIABLE");
			if (model_en) {
				model.extend(sp->value);
				model.print();
			}
		}
		else if (cnfstate == UNSAT) PFLOGS("UNSATISFIABLE");
		else if (cnfstate == UNSOLVED) PFLOGS("UNKNOWN");
		if (report_en) printReport();
	}

}