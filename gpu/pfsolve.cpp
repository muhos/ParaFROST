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

#include "pfopts.h"
#include "pfsort.h"
#include "pfsolve.h" 
#include "pfdimacs.h"

/********************************/
/*     Global Instantiation     */
/********************************/
CNF_INFO cnf_stats;
cudaDeviceProp devProp;
int maxGPUThreads = 0;
LIT_ST* assigns = NULL;
int* levels = NULL;
ParaFROST* gpfrost = NULL;
/***********************************/
/*     Resources queries           */
/***********************************/
int64 getAvailSysMem()
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

void set_timeout(int time_limit)
{
#ifdef __linux__ 
	if (time_limit != 0) {
		rlimit rl;
		getrlimit(RLIMIT_CPU, &rl);
		if (rl.rlim_max == RLIM_INFINITY || (rlim_t)time_limit < rl.rlim_max) {
			rl.rlim_cur = time_limit;
			if (setrlimit(RLIMIT_CPU, &rl) == -1) printf("WARNING - Timeout cannot be set\n");
		}
	}
#elif _WIN32
	printf("c | WARNING - timeout not supported on Windows.\n");
#endif
}

void handler_terminate(int)
{
	fflush(stdout);
	printf("c |\n");
	printf("c |%45s\n", "Interrupted");
	printf("c |\n");
	printf("s UNKNOWN\n");
	printf("c |\n");
	printf("c |--------------------------------------------------------------------------------------|\n");
	_exit(EXIT_FAILURE);
}

void handler_mercy_intr(int)
{
	fflush(stdout);
	if (!gpfrost->quiet_en) {
		printf("c |\n");
		printf("c |%45s\n", "Interrupted");
		printf("c |\n");
	}
	gpfrost->interrupt();
}

void handler_mercy_timeout(int)
{
	fflush(stdout);
	if (!gpfrost->quiet_en) {
		printf("c |\n");
		printf("c |%45s\n", "Timeout");
		printf("c |\n");
	}
	gpfrost->interrupt();
}

void sig_handler(void h_intr(int), void h_timeout(int))
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
ParaFROST::ParaFROST(const string& path) :
	timeout(opt_timeout), verbose(opt_verbose), progRate(opt_progress), initProgRate(opt_progress), seed(opt_seed),
	quiet_en(opt_quiet_en), parse_only_en(opt_par_en), rewriter_en(opt_rew_en), perf_en(opt_perf_en),
	mcv_en(!opt_lcv_en), model_en(opt_model_en), proof_en(opt_proof_en),
	fdp_en(opt_fdp_en), cbt_en(opt_cbt_en), pre_en(opt_pre_en), lpre_en(opt_lpre_en), pre_delay(opt_pre_delay),
	var_inc(opt_var_inc), var_decay(opt_var_decay), VSIDSDecayFreq(opt_var_dfreq),
	pdm_rounds(opt_pdm_rounds), pdm_freq(opt_pdm_freq), pdm_order(opt_pdm_ord),
	SH(opt_SH), blRestMin(opt_bl_rest_min),
	restart_base(opt_luby_base), restart_inc(opt_luby_inc), restPolicy(opt_restart),
	polarity(opt_pol), lbdFrozen(opt_lbd_frozen), lbdMinReduce(opt_lbd_min), lbdMinClSize(opt_lbd_min_size),
	nClsReduce(opt_init_red), incReduceSmall(opt_inc_red_sm), incReduceBig(opt_inc_red_bg),
	lbdRestBase(opt_lbd_rest_base), blockRestBase(opt_bl_rest_base), RF(opt_RF), RB(opt_RB),
	cbt_dist(opt_cbt_dist), cbt_conf_max(opt_cbt_confs),
	proof_path(opt_proof_out), gperc(opt_garbage_perc)
{
	intr = false, mapped = false;
	sysMem = NULL, var_heap = NULL, timer = NULL;
	sysMem_sz = 0ULL, sysMemCons = 0LL;
	sysMemAvail = getAvailSysMem();
	CHECK(cudaGetDeviceCount(&devCount));
	if (devCount == 0) printf("Error - no available GPU(s) that support CUDA\n"), exit(EXIT_FAILURE);
	else {
		printf("c | Detected (%d) CUDA-enabled GPU(s)\n", devCount);
		CHECK(cudaGetDeviceProperties(&devProp, MASTER_GPU));
		maxGPUThreads = devProp.multiProcessorCount * devProp.maxThreadsPerMultiProcessor;
		cnf = NULL, ot = NULL, pv = NULL, d_occurs = NULL, d_scores = NULL, raw_hist = NULL;
		gstats = NULL, streams = NULL;
	}
	if (!quiet_en) {
		printf("c | GPU: \"%s\" (compute cap: %d.%d)\n", devProp.name, devProp.major, devProp.minor);
		printf("c |  - Multiprocessors = %d MPs (%d cores/MP)\n", devProp.multiProcessorCount, SM2Cores(devProp.major, devProp.minor));
		printf("c |  - Available Global memory = %zd GB\n", devProp.totalGlobalMem / GBYTE);
		printf("c |  - Available Shared memory = %zd KB\n", devProp.sharedMemPerBlock / KBYTE);
		printf("c | Available system memory = %lld GB\n", sysMemAvail / GBYTE);
		printf("c |--------------------------------------------------------------------------------------|\n");
	}
	if (quiet_en) verbose = 0, perf_en = false;
	if (SH < 2 && restPolicy == "lbd") restPolicy = "luby";
	if (SH == 2 && restPolicy != "lbd") restPolicy = "lbd";
	if (pre_en) optSimp();
	if (proof_en) {
		proofFile.open(proof_path, std::ofstream::binary | std::ofstream::out);
		if (!proofFile.is_open()) {
			printf("Cannot open proof file %s\n", proof_path.c_str());
			exit(EXIT_FAILURE);
		}
	}
	this->path = path;
	timer = new TIMER();
	if (rewriter_en) { cnfrewriter(path); exit(EXIT_SUCCESS); }
	// parse cnf & check top state
	CNF_STATE top = parser(path);
	if (top == TERMINATE || parse_only_en) killSolver();
	else if (top == UNSAT || BCP() != NULL) { if (proof_en) write_proof('0'); killSolver(UNSAT); }
	simplify_top();
	if (verbose == 1) {
		printf("c |-------------------------------------- Progress --------------------------------------|\n");
		if (SH == 2) {
			const char* leq = u8"\u2264";
#ifdef _WIN32
			SetConsoleOutputCP(65001);
#endif
			printf("c |                    ORG                    Restarts                  Learnt           |\n");
			printf("c |     Vars      Cls     Bins(+/-)    Lits               Cls     GC(%s2)     Lits    L/C |\n", leq);
		}
		else {
			printf("c |               ORG              | Conflicts  |    Limit   |             Learnt        |\n");
			printf("c |     Vars      Cls      Lits    |            |     Cls    |     Cls      Lits     L/C |\n");
		}
		printf("c |--------------------------------------------------------------------------------------|\n");
	}
}

ParaFROST::~ParaFROST()
{
	if (verbose >= 1) cout << "c | Freeing up Host memory...";
	sysFree();
	if (verbose >= 1) cout << " done." << endl;
	if (verbose >= 1) cout << "c |--------------------------------------------------------------------------------------|" << endl;
}

CNF_STATE ParaFROST::parser(const string& path) {
	struct stat st;
	stat(path.c_str(), &st);
	size_t fsz = st.st_size;
	if (verbose >= 1) printf("c | Parsing CNF file \"%s\" (size: %zd KB)\n", path.c_str(), fsz / KBYTE);
	timer->start();
#ifdef __linux__
	int fd = open(path.c_str(), O_RDONLY, 0);
	if (fd == -1) printf("Cannot open input file\n"), exit(EXIT_FAILURE);
	void* buffer = mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
	char* str = (char*)buffer;
#else
	ifstream inputFile;
	inputFile.open(path, ifstream::in);
	if (!inputFile.is_open()) printf("Cannot open input file\n"), exit(EXIT_FAILURE);
	char* buffer = new char[fsz + 1], * str = buffer;
	inputFile.read(buffer, fsz);
	buffer[fsz] = '\0';
#endif
	uVec1D in_c, tmpCl;
	in_c.reserve(INIT_CAP);
	int oldSz = 0;
	char* eof = str + fsz;
	while (str < eof) {
		eatWS(str);
		if (*str == '\0' || *str == '0' || *str == '%') break;
		if (*str == 'c') eatLine(str);
		else if (*str == 'p') {
			if (!eq(str, "p cnf")) printf("Error - header has wrong format\n"), exit(EXIT_FAILURE);
			uint32 sign = 0;
			cnf_stats.n_org_vars = toInteger(str, sign);
			if (sign) printf("Error - number of variables in header is negative\n"), exit(EXIT_FAILURE);
			if (cnf_stats.n_org_vars == 0) printf("Error - zero number of variables in header\n"), exit(EXIT_FAILURE);
			cnf_stats.n_org_cls = toInteger(str, sign);
			if (sign) printf("Error - number of cls in header is negative\n"), exit(EXIT_FAILURE);
			if (cnf_stats.n_org_cls == 0) printf("Error - zero number of cls in header\n"), exit(EXIT_FAILURE);
			if (verbose >= 1) printf("c | Found header %d %d\n", cnf_stats.n_org_vars, cnf_stats.n_org_cls);
			assert(orgs.empty());
			orgs.resize(cnf_stats.n_org_cls);
			assert(wt.empty());
			wt.allocMem(cnf_stats.n_org_vars);
			allocSolver();
			initSolver();
		}
		else {
			toClause(in_c, str);
			if (proof_en) {
				oldSz = in_c.size();
				tmpCl.clear(true);
				tmpCl.resize(in_c.size());
				tmpCl.copyFrom(in_c);
			}
			if (checkClause(in_c)) { // clause not a tautology
				if (proof_en && in_c.size() < oldSz) {
					write_proof('a');
					write_proof(in_c, in_c.size());
					write_proof(0);
					write_proof('d');
					write_proof(tmpCl, in_c.size());
					write_proof(0);
				}
				if (in_c.size() == 1) {
					assert(*in_c > 0);
					uint32 v_idx = V2X(*in_c);
					if (sol->assign(v_idx) == UNDEFINED) { var_heap->remove(v_idx); enqueue(*in_c); }
					else if (sol->assign(v_idx) == ISNEG(*in_c)) return UNSAT;
				}
				else if (cnf_stats.global_n_cls + 1 > cnf_stats.n_org_cls) { printf("Error - too many cls\n"); return TERMINATE; }
				else {
					B_REF org = new BCLAUSE(in_c.size());
					org->copyLitsFrom(in_c);
					attachClause(org);
				}
			}
		}
	}
#ifdef __linux__
	if (munmap(buffer, fsz) != 0) {
		printf("Cannot clean input file %s mapping.\n", path.c_str());
		exit(EXIT_FAILURE);
	}
	close(fd);
#else
	delete[] buffer;
	inputFile.close();
#endif
	assert(nClauses() + nBins() <= nOrgCls());
	cnf_stats.n_org_bins = nBins();
	cnf_stats.n_org_lits = nLiterals() + ((int64)nBins() << 1);
	if ((int)nClauses() < orgs.size()) orgs.resize(nClauses());
	in_c.clear(true), tmpCl.clear(true);
	timer->stop();
	timer->par = timer->cpuTime();
	if (verbose >= 1) printf("c | Read %d Vars, %d Cls, and %lld Lits in %.3f sec.\n", nOrgVars(), nClauses() + nOrgBins() + sp->trail_size, nOrgLits() + sp->trail_size, timer->par);
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

void ParaFROST::killSolver(const CNF_STATE& status)
{
	wrapUp(status);
	if (proof_en) proofFile.close();
	this->~ParaFROST();
	gpfrost = NULL;
	exit(EXIT_SUCCESS);
}

void ParaFROST::sysFree()
{
	occurs.clear(true);
	scores.clear(true);
	learntLits.clear(true);
	learnt_cl.clear(true);
	trail_sz.clear(true);
	wt.clear(true);
	bins.clear(true);
	orgs.clear(true);
	gcr.clear(true);
	lbdQ.clear(true);
	trailQ.clear(true);
	mappedVars.clear(true);
	reverseVars.clear(true);
	if (sysMem != NULL) {
		free(sysMem);
		sysMem = NULL;
	}
	if (var_heap != NULL) {
		delete var_heap;
		var_heap = NULL;
	}
	if (timer != NULL) {
		delete timer;
		timer = NULL;
	}
}

void ParaFROST::allocSolver(const bool& re)
{
	if (re) {
		assert(sysMem != NULL);
		assert(var_heap != NULL);
		free(sysMem);
		delete var_heap;
		sysMem = NULL, var_heap = NULL;
	}
	assert(sysMem == NULL);
	assert(var_heap == NULL);
	lrn.max_cl_sz = nOrgVars();
	sysMem_sz = sizeof(SOL) + sizeof(SP) +
		nOrgVars() * (sizeof(uint32) * 3 + sizeof(LIT_ST) + sizeof(G_REF) + 4) +
		lrn.max_cl_sz * sizeof(uint32) * 2;
	sysMem = (addr_t)malloc(sysMem_sz);
	if (sysMem == NULL) throw MEMOUTEXCEPTION();
	memset(sysMem, 0, sysMem_sz);
	addr_t bottom = sysMem + sysMem_sz;
	// sol
	sol = (SOL*)sysMem;
	sysMem += sizeof(SOL);
	assert(sysMem < bottom);
	sol->allocMem(&sysMem, nOrgVars());
	assigns = sol->assigns_ptr();
	levels = sol->levels_ptr();
	assert(sysMem < bottom);
	// search space
	sp = (SP*)sysMem;
	sysMem += sizeof(SP);
	assert(sysMem < bottom);
	sp->lock = (bool*)sysMem;
	sp->frozen = sp->lock + nOrgVars();
	sp->seen = sp->frozen + nOrgVars();
	sp->pol = sp->seen + nOrgVars();
	sysMem += ((size_t)nOrgVars() << 2);
	assert(sysMem < bottom);
	sp->trail = (uint32*)sysMem;
	sysMem += (size_t)nOrgVars() * sizeof(uint32);
	assert(sysMem < bottom);
	// others
	source = (G_REF*)sysMem;
	sysMem += nOrgVars() * sizeof(G_REF);
	assert(sysMem < bottom);
	board = (int*)sysMem;
	sysMem += nOrgVars() * sizeof(int);
	assert(sysMem < bottom);
	tmp_stack = (uint32*)sysMem;
	sysMem += lrn.max_cl_sz * sizeof(uint32);
	assert(sysMem < bottom);
	simpLearnt = (uint32*)sysMem;
	sysMem += lrn.max_cl_sz * sizeof(uint32);
	assert(sysMem == bottom);
	sysMem -= sysMem_sz;
	// VSIDS heap
	var_heap = new VAR_HEAP();
	var_heap->allocMem(nOrgVars());
	// initially alloc queues & garbage collector
	if (SH == 2) { lbdQ.alloc(lbdRestBase); trailQ.alloc(blockRestBase); }
	gcr.init();
	if (verbose >= 2) printf("C | Memory consumed after allocSolver call = %lld MB\n", sysMemUsed() / MBYTE);
}

void ParaFROST::resetSolver(const bool& re) {
	nConflicts = UNKNOWN;
	ref_vars = UNKNOWN;
	cnf_stats.n_added_lits = cnf_stats.global_n_gcs = UNKNOWN;
	lrn.init();
	stats.reset();
	if (SH == 2) {
		maxConflicts = UNKNOWN;
		marker = UNKNOWN;
		reductions = 1;
		lrn.nClsReduce = nClsReduce;
		lbdSum = 0.0;
		lbdQ.reset();
		trailQ.reset();
	}
	else lrn.max_learnt_cls = nOrgCls() * lrn.size_factor_cls;
	var_heap->init(var_inc, var_decay);
	var_heap->build(nOrgVars());
	if (re) {
		for (uint32 v = 0; v < nOrgVars(); v++) { sp->pol[v] = true; sol->init(v); }
	}
	else {
		for (uint32 v = 0; v < nOrgVars(); v++) { sp->pol[v] = true; sp->lock[v] = false; source[v] = NULL; sol->init(v); }
		sp->reset_trail();
	}
}

void ParaFROST::initSolver()
{
	assert(UNKNOWN == 0);
	assert(UNDEFINED < 0 && ROOT_LEVEL == 0);
	srand(seed);
	restarts = starts = UNKNOWN;
	cl_params.init();
	resetSolver();
}

void ParaFROST::cnfrewriter(const string& path) {
	ifstream inputFile;
	inputFile.open(path + ".cnf", ifstream::in);
	if (!inputFile.is_open()) {
		cout << "c | Cannot open the input file: " << path << endl;
		exit(EXIT_FAILURE);
	}
	std::ofstream outFile;
	outFile.open(path + "_re.cnf", std::ofstream::out);
	if (!inputFile.is_open()) {
		cout << "c | Cannot open the output file: " << path << endl;
		exit(EXIT_FAILURE);
	}
	char* line = new char[MBYTE];
	if (verbose >= 1) printf("c | Parsing CNF file \"%s\"\n", path.c_str());
	outFile << "c rewritten formula by ParaFROST." << endl;
	while (inputFile.getline(line, MBYTE)) {
		char* ch = line; int len = 0;
		while (*ch) ch++;
		len = ch - line;
		if (len == 0 || (len == 1 && (*line == '0' || (*line >= '%' && *line <= '/')))) continue;
		if (*line != 'c') {
			if (*line == 'p') {
				if (!eq(line, "p cnf")) printf("Error - header has wrong format\n"), exit(EXIT_FAILURE);
				uint32 sign = 0;
				cnf_stats.n_org_vars = toInteger(line, sign);
				if (cnf_stats.n_org_vars == 0) printf("Error - zero number of variables in header\n"), exit(EXIT_FAILURE);
				cnf_stats.n_org_cls = toInteger(line, sign);
				if (cnf_stats.n_org_cls == 0) printf("Error - zero number of cls in header\n"), exit(EXIT_FAILURE);
				if (verbose >= 1) printf("c | Found header %d %d\n", cnf_stats.n_org_vars, cnf_stats.n_org_cls);
				outFile << line << endl;
			}
			else outFile << line << endl;
		}
	}
	delete[] line;
	inputFile.close();
	outFile.close();
}

void ParaFROST::attachClause(B_REF& c)
{
	CL_LEN sz = c->size();
	assert(sz > 1);
	assert(*c != NULL);
	assert(**c > 0 && *(*c + 1) > 0);
	assert(**c <= UINT32_MAX && *(*c + 1) <= UINT32_MAX);
	wt[FLIP(**c)].push(WATCH(c, *(*c + 1)));
	wt[FLIP(*(*c + 1))].push(WATCH(c, **c));
	if (sz == 2) {
		c->markOrgBin();
		bins.push(c);
		cnf_stats.global_n_bins++;
	}
	else {
		c->set_status(ORIGINAL);
		orgs[cnf_stats.global_n_cls++] = c;
		cnf_stats.global_n_lits += sz;
		if (sz > cnf_stats.max_org_cl_width) cnf_stats.max_org_cl_width = sz;
	}
}

void ParaFROST::attachClause(C_REF& c)
{
	CL_LEN sz = c->size();
	assert(sz > 2);
	assert(*c != NULL);
	assert(**c > 0 && *(*c + 1) > 0);
	assert(**c <= UINT32_MAX && *(*c + 1) <= UINT32_MAX);
	wt[FLIP(**c)].push(WATCH(c, *(*c + 1)));
	wt[FLIP(*(*c + 1))].push(WATCH(c, **c));
	c->set_status(LEARNT);
	c->markReason();
	if (SH == 1) { // csr
		double act = sz;
		if (sz > 12) act += drand();
		c->set_act(float(act));
	}
	else c->set_act(cl_params.cl_inc);
	if (SH == 2) { // lbd
		assert(sp->learnt_lbd > 0);
		if (sp->learnt_lbd <= 2) cnf_stats.global_n_gcs++; // glues
		c->melt(); // assume the clause is molten
		c->set_LBD(sp->learnt_lbd);
	}
	learnts.push(c);
	cnf_stats.n_added_lits += sz;
}

void ParaFROST::detachClause(B_REF& c, const bool& gc)
{
	if (gc) {
		wt.collect(FLIP(**c));
		wt.collect(FLIP(*(*c + 1)));
	}
	else {
		wt.remWatch(FLIP(**c), c);
		wt.remWatch(FLIP(*(*c + 1)), c);
	}
}

void ParaFROST::removeClause(B_REF& c, const bool& gc) {
	assert(c->size() >= 2);
	if (c->size() == 2) cnf_stats.global_n_bins--;
	else {
		cnf_stats.global_n_lits -= c->size();
		if (proof_en) {
			write_proof('d');
			write_proof(*c, c->size());
			write_proof(0);
		}
	}
	detachClause(c, gc);
	collectClause(c, gc);
}

void ParaFROST::removeClause(C_REF& c, const bool& gc)
{
	assert(c->status() == LEARNT);
	assert(c->size() > 2);
	cnf_stats.n_added_lits -= c->size();
	if (proof_en) {
		write_proof('d');
		write_proof(*c, c->size());
		write_proof(0);
	}
	if (gc) {
		wt.collect(FLIP(**c));
		wt.collect(FLIP(*(*c + 1)));
	}
	else {
		wt.remWatch(FLIP(**c), c);
		wt.remWatch(FLIP(*(*c + 1)), c);
	}
	collectClause(c, gc);
}

void ParaFROST::shrinkClause(G_REF gref)
{
	B_REF c = (B_REF)gref;
	assert(c->size() > 2);
	assert(assigns[V2X(**c)] == UNDEFINED && assigns[V2X(*(*c + 1))] == UNDEFINED);
	CL_LEN sz = c->size();
	if (sz == 3 && assigns[V2X((*c)[2])] == ISNEG((*c)[2])) c->pop();
	else {
		uVec1D tmpCl;
		if (proof_en) { tmpCl.resize(sz); tmpCl.copyFrom(*c, sz); }
		LIT_POS k = 2;
		while (k < sz) if (assigns[V2X((*c)[k])] == ISNEG((*c)[k])) (*c)[k] = (*c)[--sz]; else k++;
		int remLits = c->size() - sz;
		if (proof_en && remLits) {
			write_proof('a');
			write_proof(*c, sz);
			write_proof(0);
			write_proof('d');
			write_proof(tmpCl, c->size());
			write_proof(0);
		}
		c->shrink(remLits);
	}
	assert(c->size() > 1);
}

int ParaFROST::simplify(BCNF& cnf)
{
	if (cnf.size() == 0) return 0;
	int n = 0;
	for (int i = 0; i < cnf.size(); i++) {
		B_REF c = cnf[i];
		if (c->satisfied(assigns)) removeClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) cnf[n++] = c;
			else { c->markOrgBin(); c->set_status(UNKNOWN); cnf_stats.global_n_bins++; }
		}
	}
	return n;
}

int ParaFROST::simplify(LCNF& cnf)
{
	if (cnf.size() == 0) return 0;
	int n = 0;
	for (int i = 0; i < cnf.size(); i++) {
		C_REF c = cnf[i];
		if (c->satisfied(assigns)) removeClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) cnf[n++] = c;
			else { c->set_status(UNKNOWN); cnf_stats.global_n_bins++; }
		}
	}
	return n;
}

void ParaFROST::simplify()
{
	if (sp->trail_size == sp->trail_offset || lrn.simp_props > 0) return;
	if (verbose >= 2) printf("c | Simplifying CNF..");
	assert(sp->trail_head == sp->trail_size);
	assert(nClauses() > 0);
	assert(sp->trail_size > 0 && DL() == ROOT_LEVEL);
	assert(sol->assign(V2X(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	// simplify watch table
	for (int i = sp->trail_offset; i < sp->trail_size; i++) {
		uint32 assign = sp->trail[i], assign_f = FLIP(assign);
		wt.collect(assign);
		wt.collect(assign_f);
		WL& ws = wt[assign_f];
		for (int j = 0; j < ws.size(); j++) {
			B_REF c = (B_REF)ws[j].c_ref;
			if (!c->garbage() && c->size() == 2) // original or learnt (undeleted) binary
				removeClause(c);
		}
		uint32 assign_idx = V2X(assign);
		if (var_heap->has(assign_idx)) var_heap->remove(assign_idx); // remove assign from the heap
		if (pre_en) removed.push(assign); // save in case preprocessing mapped variables
	}
	// simplify input CNF
	cnf_stats.global_n_cls = simplify(orgs);
	orgs.resize(cnf_stats.global_n_cls);
	// simplify learnt CNF
	int numLearnts = simplify(learnts);
	learnts.resize(numLearnts);
	cnf_stats.global_n_del_vars += (sp->trail_size - sp->trail_offset);
	sp->trail_offset = sp->trail_size;
	lrn.simp_props = nLiterals() + nLearntLits();
	if (verbose >= 2) printf(" ==> done\n");
	if (verbose >= 2) printf("c | Recycling solver garbage..");
	recycle();
	if (verbose >= 2) printf(" ==> done\n");
}

void ParaFROST::simplify_top()
{
	if (sp->trail_size == sp->trail_offset) return;
	if (verbose >= 2) printf("c | Simplifying CNF..");
	assert(sp->trail_offset == 0);
	assert(sp->trail_head == sp->trail_size);
	assert(nClauses() > 0);
	assert(sp->trail_size > 0 && DL() == ROOT_LEVEL);
	assert(sol->assign(V2X(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	// remove root-level watch lists
	for (int i = sp->trail_offset; i < sp->trail_size; i++) {
		uint32 assign = sp->trail[i];
		wt.collect(assign);
		wt.collect(FLIP(assign));
		uint32 assign_idx = V2X(assign);
		if (var_heap->has(assign_idx)) var_heap->remove(assign_idx);
		if (pre_en) removed.push(assign);
	}
	// simplify bins
	cnf_stats.global_n_bins = 0;
	for (int i = 0; i < bins.size(); i++) {
		B_REF c = bins[i];
		if (c->satisfied(assigns)) { detachClause(c); collectClause(c); }
		else bins[cnf_stats.global_n_bins++] = c;
	}
	bins.resize(cnf_stats.global_n_bins);
	// simplify orgs
	cnf_stats.global_n_cls = 0;
	for (int i = 0; i < orgs.size(); i++) {
		B_REF c = orgs[i];
		if (c->satisfied(assigns)) removeClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) orgs[cnf_stats.global_n_cls++] = c;
			else {
				c->markOrgBin(); c->set_status(UNKNOWN);
				bins.push(c);
			}
		}
	}
	orgs.resize(cnf_stats.global_n_cls);
	cnf_stats.global_n_bins = bins.size();
	cnf_stats.global_n_del_vars += (sp->trail_size - sp->trail_offset);
	sp->trail_offset = sp->trail_size;
	assert(consistent(bins, wt));
	assert(consistent(orgs, wt));
	if (verbose >= 2) printf(" ==> done\n");
	if (verbose >= 2) printf("c | Recycling solver garbage..");
	recycle();
	if (verbose >= 2) printf(" ==> done\n");
}

void ParaFROST::enqueue(const uint32& assign, const int& pLevel, const G_REF src)
{
	assert(assign > 0);
	register uint32 assign_idx = V2X(assign);
	source[assign_idx] = src;
	sp->lock[assign_idx] = true;
	sp->trail[sp->trail_size++] = assign;
	sol->set_assign(assign_idx, !ISNEG(assign));
	sol->set_level(assign_idx, pLevel);
	if (verbose >= 3) printf("c | New %s: %d@%d.\n", src == NULL ? "decision" : "imp", (ISNEG(assign)) ? -int(ABS(assign)) : ABS(assign), pLevel);
}

void ParaFROST::binSelfsub()
{
	marker++;
	for (LIT_POS i = 1; i < learnt_cl.size(); i++) // mark all learnt literals except asserting
		board[V2X(learnt_cl[i])] = marker;
	// check binary cls if subsuming learnt clause
	uint32 parent = FLIP(learnt_cl[0]);
	WL& ws = wt[parent];
	int nLitsRem = 0;
	for (int i = 0; i < ws.size(); i++) {
		if (((B_REF)ws[i].c_ref)->size() == 2) {
			uint32 imp = ws[i].imp, imp_idx = V2X(imp);
			if (board[imp_idx] == marker && sol->assign(imp_idx) == !ISNEG(imp)) {
				board[imp_idx] = marker - 1; // unmark
				nLitsRem++;
				//printClause(((B_REF)ws[i].c_ref)->clause());
			}
		}
	}
	int last = learnt_cl.size() - 1;
	if (nLitsRem > 0) {
		//printLearnt(learnt_cl);
		for (int i = 1; i < learnt_cl.size() - nLitsRem; i++) {
			if (board[V2X(learnt_cl[i])] != marker) {
				uint32 tmpLit = learnt_cl[last];
				learnt_cl[last] = learnt_cl[i];
				learnt_cl[i] = tmpLit;
				last--; i--;
			}
		}
		//printLearnt(learnt_cl);
		learnt_cl.shrink(nLitsRem);
	}
}

bool ParaFROST::selfsub(const uint32& learntLit, uint32* tmp_stack, CL_LEN& tmp_tail, const uint32& min_level)
{
	assert(learntLit > 0);
	assert(source[V2X(learntLit)] != NULL);
	assert(sol->level(V2X(learntLit)) > ROOT_LEVEL);
	CL_LEN tmp_size = tmp_tail;
	uint32* head = simpLearnt, * sz = simpLearnt + 1;
	*head = learntLit;
	while (head != sz) {
		if (verbose >= 4) { printf("c | simpLearnt = "); printLit(*head); }
		assert(*head > 0);
		B_REF c = (B_REF)source[V2X(*head)];
		assert(c != NULL);
		assert(**c == FLIP(*head));
		for (LIT_POS l = 1; l < c->size(); l++) {
			uint32 parent = (*c)[l], parent_idx = V2X(parent);
			int parent_dl = sol->level(parent_idx);
			if (!sp->seen[parent_idx] && parent_dl > 0) {
				if (source[parent_idx] != NULL && (MAPHASH(parent_dl) & min_level) != 0) {
					sp->seen[parent_idx] = 1;
					*sz++ = parent;
					tmp_stack[tmp_tail++] = parent;
				}
				else {
					for (uint32* k = tmp_stack + tmp_size; k != tmp_stack + tmp_tail; k++)
						sp->seen[V2X(*k)] = 0;
					tmp_tail = tmp_size;
					return false;
				}
			}
		}
		head++;
	}
	return true;
}

void ParaFROST::selfsub()
{
	assert(learnt_cl.size() > 1);
	assert(lrn.max_cl_sz > learnt_cl.size());
	if (verbose >= 3) printLearnt(learnt_cl);
	register uint32* lit = learnt_cl, * tmp_lit = tmp_stack, * end = lit + learnt_cl.size();
	while (lit != end) *tmp_lit++ = *lit++;
	uint32 min_level = 0;
	for (LIT_POS k = 1; k < learnt_cl.size(); k++)  min_level |= MAPHASH(sol->level(V2X(learnt_cl[k])));
	CL_LEN new_sz = 1, tmp_sz = (CL_LEN)learnt_cl.size();
	for (LIT_POS k = 1; k < learnt_cl.size(); k++) {
		if (source[V2X(learnt_cl[k])] == NULL || !selfsub(learnt_cl[k], tmp_stack, tmp_sz, min_level))
			learnt_cl[new_sz++] = learnt_cl[k];
	}
	tmp_lit = tmp_stack; end = tmp_stack + tmp_sz;
	while (tmp_lit < end) sp->seen[V2X(*tmp_lit++)] = 0;
	learnt_cl.resize(new_sz);
}

void ParaFROST::btLevel()
{
	if (learnt_cl.size() == 1) sp->bt_level = ROOT_LEVEL;
	else if (learnt_cl.size() == 2) sp->bt_level = sol->level(V2X(learnt_cl[1]));
	else {
		LIT_POS max_k = 1;
		for (LIT_POS k = 2; k < learnt_cl.size(); k++) {
			if (sol->level(V2X(learnt_cl[k])) > sol->level(V2X(learnt_cl[max_k])))
				max_k = k;
		}
		uint32 max_k_lit = learnt_cl[max_k];
		learnt_cl[max_k] = learnt_cl[1];
		learnt_cl[1] = max_k_lit;
		sp->bt_level = sol->level(V2X(max_k_lit));
	}
}

void ParaFROST::cbtLevel(G_REF& gref)
{
	B_REF c = (B_REF)gref;
	sp->max1Found = false;
	sp->cbt_level = sol->level(V2X(**c));
	if (DL() == sp->cbt_level && DL() == sol->level(V2X(*(*c + 1)))) return;
	// find the highest level in conflicting clause beyond literal-0
	sp->max1Found = true;
	int maxIdx = 0;
	for (LIT_POS l = 1; l < c->size(); l++) {
		register int litLevel = sol->level(V2X((*c)[l]));
		if (litLevel > sp->cbt_level) {
			sp->max1Found = true;
			sp->cbt_level = litLevel;
			maxIdx = l;
		}
		else if (litLevel == sp->cbt_level && sp->max1Found == true)
			sp->max1Found = false;
	}
	if (maxIdx) { // found max. level, swap with w0-literal
		uint32 w0_lit = **c;
		**c = (*c)[maxIdx];
		(*c)[maxIdx] = w0_lit;
		if (maxIdx > 1) { // found new watch at 0-position
			assert(w0_lit == (*c)[maxIdx]);
			wt.remWatch(FLIP(w0_lit), c);
			wt[FLIP(**c)].push(WATCH(c, *(*c + 1)));
		}
	}
}

void ParaFROST::analyze(G_REF& gref)
{
	if (verbose >= 3) {
		printf("c | Analyzing conflict:\n");
		printTrail(sp->trail, sp->trail_size);
	}
	learnt_cl.clear();
	learnt_cl.push(UNKNOWN);
	int conf_level = UNDEFINED;
	int parent = 0, track = 0, index = sp->trail_size - 1;
	C_REF c = (C_REF)gref;
	if (cbt_en) conf_level = sol->level(V2X(**c));
	else conf_level = DL();
	do {
		assert(c != NULL);
		if (SH != 1 && c->status() == LEARNT) clBumpAct(c);
		if (SH == 2 && c->status() == LEARNT && c->LBD() > GLUE) { // try to update LBD of a previously learnt cls 
			uint32 currLBD = calcLBD(c);
			if (currLBD + 1 < c->LBD()) {
				if (c->LBD() <= lbdFrozen) c->freeze();
				c->set_LBD(currLBD);
			}
		}
		for (LIT_POS j = (parent == 0) ? 0 : 1; j < c->size(); j++) {
			uint32 q = (*c)[j], q_idx = V2X(q);
			int q_dl = sol->level(q_idx);
			if (!sp->seen[q_idx] && q_dl > ROOT_LEVEL) {
				var_heap->varBumpAct(q_idx);
				sp->seen[q_idx] = 1;
				if (q_dl >= conf_level) {
					track++;
					if (SH == 2 && source[q_idx] != NULL && ((C_REF)source[q_idx])->status() == LEARNT)
						learntLits.push(q);
				}
				else learnt_cl.push(q);
			}
		}
		// next implication clause
		if (cbt_en) {
			do {
				while (!sp->seen[V2X(sp->trail[index--])]);
				parent = sp->trail[index + 1];
				assert(parent > 0);
			} while (sol->level(V2X(parent)) < conf_level);
		}
		else {
			while (!sp->seen[V2X(sp->trail[index--])]);
			parent = sp->trail[index + 1];
			assert(parent > 0);
		}
		c = (C_REF)source[V2X(parent)];
		sp->seen[V2X(parent)] = 0;
		track--;
	} while (track > 0);
	assert(learnt_cl[0] == UNKNOWN);
	learnt_cl[0] = FLIP(parent);
	stats.max_lits += learnt_cl.size();
	if (learnt_cl.size() == 1) sp->bt_level = ROOT_LEVEL;
	else { // simplify learnt clause 
		selfsub();
		if (SH == 2 && learnt_cl.size() <= lbdMinClSize && calcLBD(learnt_cl) <= lbdMinReduce)
			binSelfsub();
		btLevel();
	}
	if (SH == 2) {
		sp->learnt_lbd = calcLBD(learnt_cl); // calculate LBD of minimized learnt clause
		if (verbose == 4) printf("c | LBD of learnt clause = %d\n", sp->learnt_lbd);
		if (learntLits.size() > 0) {
			for (int i = 0; i < learntLits.size(); i++) { // pump all variables having lower LBD values
				uint32 q_idx = V2X(learntLits[i]);
				c = (C_REF)source[q_idx];
				if (c->status() == LEARNT && c->LBD() < sp->learnt_lbd) // NOTE: LBD call must be protected by status check to gurantee a valid LBD value
					var_heap->varBumpAct(q_idx);
			}
			learntLits.clear();
		}
	}
	stats.tot_lits += learnt_cl.size();
	if (verbose >= 3) printLearnt(learnt_cl);
}

void ParaFROST::cancelAssigns(const int& bt_level)
{
	// Non-Chrono. BT
	if (!cbt_en && DL() > bt_level) {
		for (int i = sp->trail_size - 1; i >= trail_sz[bt_level]; i--) {
			assert(sp->trail[i] > 0);
			register uint32 assign = sp->trail[i], assign_idx = V2X(assign);
			sol->init(assign_idx);
			if (!var_heap->has(assign_idx)) var_heap->insert(assign_idx);
			if (polarity == 1) sp->pol[assign_idx] = ISNEG(assign);
			else if (polarity == -1) sp->pol[assign_idx] = !ISNEG(assign);
			if (source[assign_idx] != NULL) {
				((B_REF)source[assign_idx])->initReason();
				source[assign_idx] = NULL;
			}
			sp->lock[assign_idx] = false;
		}
		sp->trail_head = sp->trail_size = trail_sz[bt_level];
		trail_sz.shrink(trail_sz.size() - bt_level);
		assert(DL() == trail_sz.size());
	}
	// Chrono. BT
	if (cbt_en && DL() > bt_level) {
		int nLitsRem = 0;
		for (int i = sp->trail_size - 1; i >= trail_sz[bt_level]; i--) {
			assert(sp->trail[i] > 0);
			register uint32 assign = sp->trail[i], assign_idx = V2X(assign);
			if (sol->level(assign_idx) <= bt_level)
				tmp_stack[nLitsRem++] = sp->trail[i];
			else {
				sol->init(assign_idx);
				if (!var_heap->has(assign_idx)) var_heap->insert(assign_idx);
				if (polarity == 1) sp->pol[assign_idx] = ISNEG(assign);
				else if (polarity == -1) sp->pol[assign_idx] = !ISNEG(assign);
				if (source[assign_idx] != NULL) {
					((B_REF)source[assign_idx])->initReason();
					source[assign_idx] = NULL;
				}
				sp->lock[assign_idx] = false;
			}
		}
		sp->trail_head = sp->trail_size = trail_sz[bt_level];
		trail_sz.shrink(trail_sz.size() - bt_level);
		for (int i = nLitsRem - 1; i >= 0; i--) sp->trail[sp->trail_size++] = tmp_stack[i];
	}
}

void ParaFROST::backJump(const int& bt_level) {
	if (verbose >= 3) printf("c | Backjumping to level %d\n", bt_level);
	assert(sp->trail_size > 0);
	// cancel old assignments up to backtrack level <bt_level>
	cancelAssigns(bt_level);
	// add learnt clause & enqueue learnt decision
	if (proof_en) {
		write_proof('a');
		write_proof(learnt_cl.data(), learnt_cl.size());
		write_proof(0);
	}
	if (learnt_cl.size() == 1)
		enqueue(learnt_cl[0]);
	else if (learnt_cl.size() == 2) { // Note: LBD of binary learnts is always 2, no need to create C_REF clause
		cnf_stats.global_n_bins++; cnf_stats.global_n_gcs++;
		B_REF learnt_bin = new BCLAUSE(learnt_cl.size());
		learnt_bin->copyLitsFrom(learnt_cl);
		wt[FLIP(learnt_cl[0])].push(WATCH(learnt_bin, learnt_cl[1]));
		wt[FLIP(learnt_cl[1])].push(WATCH(learnt_bin, learnt_cl[0]));
		enqueue(learnt_cl[0], sp->bt_level, learnt_bin); // enqueue must always use sp->bt_level
	}
	else if (learnt_cl.size() > 2) {
		C_REF learnt = new LCLAUSE(learnt_cl.size());
		learnt->copyLitsFrom(learnt_cl);
		attachClause(learnt);
		enqueue(learnt_cl[0], sp->bt_level, learnt); // enqueue must always use sp->bt_level
	}
}

void ParaFROST::lReduce()
{
	if (verbose >= 2) printf("c | Reducing learnt database..");
	int n = 0;
	if (SH == 1) { // CSR
		Sort(learnts, LEARNT_SR());
		int pivot = learnts.size() / 2.1;
		for (int i = 0; i < learnts.size(); i++) {
			C_REF c = learnts[i];
			if (i < pivot && !c->reason()) {
				assert(c->size() > 2);
				removeClause(c);
			}
			else learnts[n++] = c;
		}
	}
	else { // activity or LBD
		Sort(learnts, LEARNT_CMP());
		if (SH == 0) {
			double act_limit = cl_params.cl_inc / learnts.size();
			int pivot = learnts.size() / 2.1;
			for (int i = 0; i < learnts.size(); i++) {
				C_REF c = learnts[i];
				if (!c->reason() && (i < pivot || c->activity() < act_limit)) {
					assert(c->size() > 2);
					removeClause(c);
				}
				else learnts[n++] = c;
			}
		}
		else {
			int pivot = learnts.size() / 2;
			if (learnts[pivot]->LBD() <= 3) lrn.nClsReduce += incReduceBig;
			for (int i = 0; i < learnts.size(); i++) {
				C_REF c = learnts[i];
				if (c->LBD() > GLUE && c->molten() && !c->reason() && i < pivot) {
					assert(c->size() > 2);
					removeClause(c);
				}
				else {
					if (!c->molten()) pivot++;
					c->melt();
					learnts[n++] = c;
				}
			}
		}
	}
	learnts.resize(n);
	if (verbose >= 2) printf(" ==> done.\n");
	if (verbose >= 2) printf("c | Recycling solver garbage..");
	recycle();
	if (verbose >= 2) printf(" ==> done\n");
}

G_REF ParaFROST::BCP()
{
	G_REF conf_ref = NULL;
	int numProps = 0;
	while (sp->trail_head < sp->trail_size) {
		numProps++;
		uint32 assign = sp->trail[sp->trail_head++], assign_idx = V2X(assign), assign_dl = sol->level(assign_idx);
		assert(assign > 0);
		if (verbose >= 4) printf("c | Propagating assign("), printLit(assign), printf("):\n");
		if (verbose >= 4) printWatched(assign_idx);
		WL& ws = wt.getClean(assign);
		if (ws.size()) {
			WATCH* w_i = ws, * w_j = w_i, * end = w_i + ws.size();
			while (w_i != end) {
				G_REF gref = w_i->c_ref;
				B_REF c = (B_REF)gref;
				uint32 wlit = w_i->imp;
				//c->print();
				assert(wlit && wlit < UINT32_MAX);
				if (sol->assign(V2X(wlit)) == !ISNEG(wlit)) { *w_j++ = *w_i++; continue; }
				// move assigned-0 literal to watch 1
				assert(c->size() >= 2);
				uint32 f_assign = FLIP(assign);
				if (**c == f_assign) c->swap_ws();
				assert(*(*c + 1) == f_assign);
				w_i++;
				// check if first literal is true
				register uint32 w0_lit_idx = V2X(**c);
				WATCH w = WATCH(gref, **c);
				if (**c != wlit && sol->assign(w0_lit_idx) == !ISNEG(**c)) *w_j++ = w;
				else {
					for (LIT_POS k = 2; k < c->size(); k++) { // look for (un)assigned-1 literal to watch
						register uint32 lit = *(*c + k);
						LIT_ST h = sol->assign(V2X(lit));
						if (h == UNDEFINED || h == !ISNEG(lit)) {
							*(*c + 1) = lit;
							*(*c + k) = f_assign;
							wt[FLIP(lit)].push(w);
							goto nextC;
						}
					}
					// clause is unit or conflict
					*w_j++ = w;
					if (sol->assign(w0_lit_idx) == UNDEFINED) {
						assert(sol->level(w0_lit_idx) == UNDEFINED);
						c->markReason();
						if (assign_dl == DL()) enqueue(**c, assign_dl, c);
						else {
							// find parent with max. level
							int maxLevel = assign_dl, maxIdx = 1;
							for (LIT_POS k = 2; k < c->size(); k++) {
								register int litLevel = sol->level(V2X(*(*c + k)));
								if (litLevel > maxLevel) { maxLevel = litLevel; maxIdx = k; }
							}
							if (maxIdx != 1) {
								c->swap(1, maxIdx);
								*w_j--; // remove (w) from assign list
								wt[FLIP((*c)[1])].push(w); // add it as new 1-literal watched 
							}
							enqueue(**c, maxLevel, c);
						}
					}
					else {
						if (verbose >= 3) printf("c | Conflict on \n"), printLit(**c);
						sp->trail_head = sp->trail_size;
						conf_ref = gref;
						while (w_i < end) *w_j++ = *w_i++;
					}
				}
			nextC:;
			}
			ws.shrink(int(w_i - w_j));
		}
		if (verbose >= 4) printWatched(assign_idx);
	}
	stats.n_props += numProps;
	lrn.simp_props -= numProps;
	return conf_ref;
}

CNF_STATE ParaFROST::decide()
{
	if (R > 0 && ((int)nOrgVars() - sp->trail_size) > ref_vars) PDM();
	if (sp->trail_size - sp->trail_head == 0) {
		uint32 dec = UNKNOWN;
		int cand = UNDEFINED;
		while (!var_heap->empty() && dec == UNKNOWN) {
			cand = var_heap->removeMin();
			assert(cand != UNDEFINED && cand < (int)nOrgVars());
			if (!sp->lock[cand]) {
				dec = (polarity != 0) ? (sp->pol[cand] ? NEG(V2D(cand + 1)) : V2D(cand + 1)) : (drand() < 0.5 ? NEG(V2D(cand + 1)) : V2D(cand + 1));
				assert(sp->trail_size < (int)nOrgVars());
			}
		}
		if (dec == UNKNOWN) return SAT;
		incDL();
		enqueue(dec, DL());
	}
	return UNSOLVED;
}

CNF_STATE ParaFROST::search()
{
	starts++;
	sp->reset_level();
	sp->learnt_lbd = UNKNOWN;
	int64 confs = UNKNOWN;
	if (lbdSum > 1e300) lbdSum *= 1e-300;
	while (!interrupted()) {
		if (verbose >= 3) printf("c | Current Decision Level = %d\n", DL());
		G_REF conf_ref = BCP();
		if (conf_ref != NULL) {
			nConflicts++; confs++;
			if (cbt_en) {
				cbtLevel(conf_ref);
				if (sp->cbt_level == ROOT_LEVEL) return UNSAT;
				if (sp->max1Found) {
					cancelAssigns(sp->cbt_level - 1);
					continue;
				}
			}
			else if (DL() == ROOT_LEVEL) return UNSAT;
			if (SH == 2) {
				maxConflicts++;
				// increase var decay over time
				if (nConflicts % VSIDSDecayFreq == 0 && var_heap->VarDecay() < 0.95)
					var_heap->incVarDecay(opt_var_decay_r);
				// block restarts
				trailQ.push(sp->trail_size);
				if (maxConflicts > blRestMin && lbdQ.full() && sp->trail_size > (RB * trailQ.average())) {
					lbdQ.reset();
					stats.nRestartStops++;
				}
			}
			// analyze conflict
			analyze(conf_ref);
			// chronological BT trigger
			if (cbt_en && nConflicts >= cbt_conf_max && (DL() - sp->bt_level) >= cbt_dist) {
				backJump(sp->cbt_level - 1);
				stats.cbt++;
			}
			else {
				backJump(sp->bt_level);
				stats.ncbt++;
			}
			// adjust restart/learnt rate
			if (SH < 2 && --lrn.adjust_cnt == 0) {
				lrn.adjust_conf *= lrn.adjust_inc;
				lrn.adjust_cnt = (int)lrn.adjust_conf;
				lrn.max_learnt_cls *= lrn.size_inc;
				printStats();
			}
			else if (SH == 2) {
				lbdQ.push(sp->learnt_lbd);
				lbdSum += sp->learnt_lbd;
				printStats(nConflicts == 1 || nConflicts % progRate == 0);
			}
			// update var/clause activities
			var_heap->VarDecayAct();
			if (SH != 1) clDecayAct();
		}
		else {
			// solving restart policy
			if (SH == 2 && lbdQ.full() && ((lbdQ.average() * RF) > (lbdSum / maxConflicts))) { // LBD restarts
				lbdQ.reset();
				cancelAssigns(ROOT_LEVEL);
				return UNSOLVED;
			}
			else if (SH < 2 && confs >= maxConflicts) { // luby/power restarts
				cancelAssigns(ROOT_LEVEL);
				return UNSOLVED;
			}
			// simplify (satisfied cls)
			if (DL() == ROOT_LEVEL) simplify();
			// learnt cls reduction
			if (SH == 2 && nConflicts >= reductions * lrn.nClsReduce) {
				assert(learnts.size() > 0);
				reductions = int(nConflicts / lrn.nClsReduce) + 1;
				lReduce();
				lrn.nClsReduce += incReduceSmall;
			}
			else if (SH < 2 && ((int64)learnts.size() - (int64)sp->trail_size) >= lrn.max_learnt_cls) lReduce();
			// choose next decision(s)
			stats.n_fuds++;
			if (decide() == SAT) return SAT;
		}
		if (verbose >= 3) printTrail(sp->trail, sp->trail_size);
	}
	return TERMINATE;
}

void ParaFROST::solve()
{
	timer->start();
	if (pre_en && !pre_delay) GPU_preprocess();
	else PDMInit();
	CNF_STATE status = UNSOLVED;
	double rest_fact = 1.0;
	while (status == UNSOLVED) {
		if (SH < 2) {
			rest_fact = restPolicy == "luby" ? luby_seq(restart_inc, restarts) : pow(restart_inc, restarts);
			maxConflicts = rest_fact * restart_base;
		}
		// special restart for preprocessing
		//if (pre_en && restarts >= 1 && restarts % pre_delay == 0) GPU_preprocess();
		status = search();
		restarts++;
		PDM_fuse();
	}
	timer->stop();
	timer->solve = timer->cpuTime();
	timer->solve -= timer->pre;
	if (proof_en) proofFile.close();
	wrapUp(status);
}

void ParaFROST::recycle() {
	if (gcr.size() > cnfSize() * gperc) {
		wt.recycle();
		gcr.recycle();
		for (uint32 v = 0; v < nOrgVars(); v++) {
			uint32 p = V2D(v + 1);
			wt[p].shrinkCap(); wt[NEG(p)].shrinkCap();
		}
		orgs.shrinkCap();
	}
}

bool ParaFROST::consistent(BCNF& cnf, WT& wt)
{
	for (int i = 0; i < cnf.size(); i++) {
		B_REF c = cnf[i];
		//c->print();
		WL& ws_0 = wt[FLIP(**c)];
		WL& ws_1 = wt[FLIP(*(*c + 1))];
		//print(ws_0);
		//print(ws_1);
		assert(ws_0.size() > 0);
		assert(ws_1.size() > 0);
		bool w0_found = false, w1_found = false;
		for (int j = 0; j < ws_0.size(); j++) {
			WATCH& w = ws_0[j];
			B_REF w_ref = (B_REF)w.c_ref;
			assert(w_ref->size() > 0);
			if (w_ref == c) {
				w0_found = true;
				break;
			}
		}
		for (int j = 0; j < ws_1.size(); j++) {
			WATCH& w = ws_1[j];
			B_REF w_ref = (B_REF)w.c_ref;
			assert(w_ref->size() > 0);
			if (w_ref == c) {
				w1_found = true;
				break;
			}
		}
		if (!w0_found || !w1_found) {
			printf(" Clause");
			c->print();
			printf(" NOT found in one of the watch lists:\n");
			wt.print(FLIP(**c));
			wt.print(FLIP(*(*c + 1)));
			return false;
		}
	}
	return true;
}

bool ParaFROST::consistent(LCNF& cnf, WT& wt)
{
	for (int i = 0; i < cnf.size(); i++) {
		C_REF c = cnf[i];
		//c->print();
		WL& ws_0 = wt[FLIP(**c)];
		WL& ws_1 = wt[FLIP(*(*c + 1))];
		//print(ws_0);
		//print(ws_1);
		assert(ws_0.size() > 0);
		assert(ws_1.size() > 0);
		bool w0_found = false, w1_found = false;
		for (int j = 0; j < ws_0.size(); j++) {
			WATCH& w = ws_0[j];
			C_REF w_ref = (C_REF)w.c_ref;
			assert(w_ref->size() > 0);
			if (w_ref == c) {
				w0_found = true;
				break;
			}
		}
		for (int j = 0; j < ws_1.size(); j++) {
			WATCH& w = ws_1[j];
			C_REF w_ref = (C_REF)w.c_ref;
			assert(w_ref->size() > 0);
			if (w_ref == c) {
				w1_found = true;
				break;
			}
		}
		if (!w0_found || !w1_found) {
			printf(" Clause");
			c->print();
			printf(" NOT found in one of the watch lists:\n");
			wt.print(FLIP(**c));
			wt.print(FLIP(*(*c + 1)));
			return false;
		}
	}
	return true;
}

void ParaFROST::printReport()
{
	if (perf_en) {
		printf("c |\nc |\t\t\tSolver Report\n");
		printf("c | Simplifier time      : %.3f sec\n", timer->pre);
		printf("c | Solver time          : %.3f sec\n", timer->solve);
		printf("C | Solver memory        : %.3f MB\n", ((double)sysMemUsed() / MBYTE));
		printf("c | PDM calls            : %-10d\n", stats.pdm_calls);
		printf("c | Restarts             : %-10d\n", starts);
		printf("c | Blocked restarts     : %-10d\n", stats.nRestartStops);
		printf("c | Chronological BT     : %-10d\n", stats.cbt);
		printf("c | Non-chronological BT : %-10d\n", stats.ncbt);
		printf("c | Parallel Decisions   : %-10lld  (%.1f dec/sec)\n", stats.n_pds, (stats.n_pds / timer->solve));
		printf("c | Follow-Up Decisions  : %-10lld  (%.1f dec/sec)\n", stats.n_fuds, (stats.n_fuds / timer->solve));
		printf("c | Propagations         : %-10lld  (%.1f prop/sec)\n", stats.n_props, (stats.n_props / timer->solve));
		printf("c | Conflicts            : %-10lld  (%.1f conf/sec)\n", nConflicts, (nConflicts / timer->solve));
		printf("c | Conflict literals    : %-10lld  (%.2f %% deleted)\n", stats.tot_lits, (stats.max_lits - stats.tot_lits) * 100 / (double)stats.tot_lits);
		printf("c |--------------------------------------------------------------------------------------|\n");
	}
}

void ParaFROST::printModel()
{
	printf("v ");
	for (int i = 0; i < removed.size(); i++) { // print any saved assignments
		sp->frozen[V2X(removed[i])] = true; // trail may not be reset
		printf("%d ", ISNEG(removed[i]) ? -int(ABS(removed[i])) : int(ABS(removed[i])));
	}
	if (mapped) { // recover mapped variables
		assert(sp->trail_size < reverseVars.size() - 1);
		for (int i = 0; i < sp->trail_size; i++) {
			int v = reverseVars[ABS(sp->trail[i])];
			if (!sp->frozen[v - 1]) printf("%d ", ISNEG(sp->trail[i]) ? -v : v);
		}
	}
	else {
		for (int i = 0; i < sp->trail_size; i++) {
			if (!sp->frozen[V2X(sp->trail[i])])
				printf("%d ", ISNEG(sp->trail[i]) ? -int(ABS(sp->trail[i])) : int(ABS(sp->trail[i])));
		}
	}
}

void ParaFROST::wrapUp(const CNF_STATE& status)
{
	// print results
	if (verbose >= 1) printf("c |--------------------------------------------------------------------------------------|\nc |\n");
	if (status == SAT) {
		if (!quiet_en) printf("c |\n");
		//printf("%s: ", path.c_str()); printf("s SATISFIABLE (time=%.3f)\n", timer->solve);
		printf("s SATISFIABLE\n");
		if (!quiet_en) printf("c |\n");
		if (model_en) {
			printModel();
			if (!quiet_en) printf("\nc |\n");
			else printf("\n");
		}
	}
	else if (status == UNSAT) {
		if (!quiet_en) printf("c |\n");
		//printf("%s: ", path.c_str()); printf("s UNSATISFIABLE (time=%.3f)\n", timer->solve);
		printf("s UNSATISFIABLE\n");
		if (!quiet_en) printf("c |\n");
	}
	else if (status == TERMINATE) {
		//printf("%s: ", path.c_str()); printf("s UNKNOWN (time=%.3f)\n", timer->solve);
		printf("s UNKNOWN\n");
		if (!quiet_en) printf("c |\n");
	}
	if (perf_en) printReport();
}