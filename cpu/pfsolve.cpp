/***********************************************************************
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
************************************************************************/
#include "pfsort.h"
#include "pfsolve.h" 
#include "pfparser.h"
#include "pfopts.h"

//=======================================//
//     container for CNF info.           //
//=======================================//
CNF_INFO cnf_stats;
/********************************/
/*     Aliases Initialization   */
/********************************/
ASSIGN_ST* assigns = NULL;
int* levels = NULL;
ParaFROST* g_pFrost = NULL;
/***********************************/
/*     Resources queries           */
/***********************************/
#ifdef __linux__ 
int64 getTotalSystemMemory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}
void set_timeout(int time_limit)
{
	if (time_limit != 0) {
		rlimit rl;
		getrlimit(RLIMIT_CPU, &rl);
		if (rl.rlim_max == RLIM_INFINITY || (rlim_t)time_limit < rl.rlim_max) {
			rl.rlim_cur = time_limit;
			if (setrlimit(RLIMIT_CPU, &rl) == -1) printf("WARNING -> Timeout cannot be set.\n");
		}
	}
}
#elif _WIN32
int64 getTotalSystemMemory()
{
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullTotalPhys;
}
void set_timeout(int time_limit)
{
	printf("C | WARNING - timeout not supported on Windows.\n");
}
#endif
// interrupt handlers
void handler_terminate(int)
{
	fflush(stdout);
	printf("c |\n");
	printf("c |%45s\n","Interrupted");
	printf("c |\n");
	printf("s UNKNOWN\n");
	printf("c |\n");
	printf("c |--------------------------------------------------------------------------------------|\n");
	_exit(EXIT_FAILURE);
}

void handler_mercy_intr(int) 
{
	fflush(stdout);
	if (!g_pFrost->quiet_en) {
		printf("c |\n");
		printf("c |%45s\n", "Interrupted");
		printf("c |\n");
	}
	g_pFrost->interrupt(); 
}

void handler_mercy_timeout(int) 
{
	fflush(stdout);
	if (!g_pFrost->quiet_en) {
		printf("c |\n");
		printf("c |%45s\n", "Timeout");
		printf("c |\n");
	}
	g_pFrost->interrupt();
}

void sig_handler(void h_intr(int), void h_timeout(int))
{
	signal(SIGINT, h_intr);
	signal(SIGTERM, h_intr);
#ifdef SIGXCPU
	if (h_timeout != NULL) signal(SIGXCPU, h_timeout);
#endif
}
//============================//
//  Specialized Comparators   //
//============================//
struct LEARNT_CMP {
	bool operator () (C_REF a, C_REF b) {
		if (a->LBD() != b->LBD()) return a->LBD() > b->LBD();
		else return a->activity() != b->activity() ? a->activity() < b->activity() : a->size() > b->size();
	}
};
struct LEARNT_SR {
	bool operator () (C_REF a, C_REF b) {
		return a->activity() > b->activity();
	}
};
//=======================================//
//		 ParaFROST defined members       //
//=======================================//
ParaFROST::ParaFROST(const string& path) :
	timeout(opt_timeout), verbose(opt_verbose), progRate(opt_progress), seed(opt_seed),
	quiet_en(opt_quiet_en), parse_only_en(opt_par_en), rewriter_en(opt_rew_en), perf_en(opt_perf_en),
	mcv_en(!opt_lcv_en), model_en(opt_model_en), proof_en(opt_proof_en),
	fdp_en(opt_fdp_en), cbt_en(opt_cbt_en), pre_en(opt_pre_en), lpre_en(opt_lpre_en), pre_delay(opt_pre_delay),
	var_inc(opt_var_inc), var_decay(opt_var_decay),
	pdm_rounds(opt_pdm_rounds), pdm_freq(opt_pdm_freq), pdm_order(opt_pdm_ord),
	restart_base(opt_luby_base), restart_inc(opt_luby_inc), restPolicy(opt_restart), SH(opt_SH),
	polarity(opt_pol), lbdFrozen(opt_lbd_frozen), lbdMinReduce(opt_lbd_min), lbdMinClSize(opt_lbd_min_size),
	nClsReduce(opt_init_red), incReduceSmall(opt_inc_red_sm), incReduceBig(opt_inc_red_bg),
	lbdRestBase(opt_lbd_rest_base), blockRestBase(opt_bl_rest_base), RF(opt_RF), RB(opt_RB),
	cbt_dist(opt_cbt_dist), cbt_conf_max(opt_cbt_confs),
	proof_path(opt_proof_out)
{
	intr = false;
	if (quiet_en) { verbose = 0; perf_en = false; }
	if (SH < 2 && restPolicy == "lbd") restPolicy = "luby";
	if (SH == 2 && restPolicy != "lbd") restPolicy = "lbd";
	if (pre_en) opt_simp();
	sysMem_sz = 0ULL;
	sysMemCons = 0.0;
	sysMemTot = (double)getTotalSystemMemory() - 1.0 * GBYTE;
	this->path = path;
	timer = new TIMER();
	if (proof_en) {
		proofFile.open(proof_path, std::ofstream::binary | std::ofstream::out);
		if (!proofFile.is_open()) {
			printf("Cannot open proof file %s\n", proof_path.c_str());
			exit(EXIT_FAILURE);
		}
	}
	if (rewriter_en) { CNF_rewriter(path); exit(EXIT_SUCCESS); }
	// parse cnf & check top state
	CNF_STATE top = CNF_parser(path);
	if (top == TERMINATE) { free_mem(); exit(EXIT_FAILURE); }
	else if (top == UNSAT || BCP() != NULL) {
		wrapUp(UNSAT);
		if (proof_en) { write_proof('0'); proofFile.close(); }
		free_mem(); exit(EXIT_SUCCESS);
	}
	if (parse_only_en) { free_mem(); exit(EXIT_SUCCESS); }
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
	free_mem();
	if (verbose >= 1) cout << " done." << endl;
	if (verbose >= 1) cout << "c |--------------------------------------------------------------------------------------|" << endl;
}

CNF_STATE ParaFROST::CNF_parser(const string& path) {
	ifstream inputFile;
	inputFile.open(path, ifstream::in);
	if (!inputFile.is_open()) {
		printf("Cannot open input file %s.\n", path.c_str());
		exit(EXIT_FAILURE);
	}
	if (verbose >= 1) printf("c | Parsing CNF file \"%s\"\n", path.c_str());
	char* tmp = new char[LIT_LEN];
	char* buffer = new char[BUFFER_SIZE];
	uVector1D header(2), tmpCl;
	uint32* in_c;
	int c_sz = 0, oldSz = 0;
	while (inputFile.getline(buffer, BUFFER_SIZE)) {
		std::streamsize len = inputFile.gcount();
		if (len == 0) { printf("Error - empty line.\n"); exit(EXIT_FAILURE); }
		if (len == 1 && *buffer == '\0') continue;
		if (len == 2 && (*buffer == '0' || (*buffer >= '%' && *buffer <= '/'))) continue;
		if (*buffer != 'c') {
			if (*buffer == 'p') {
				read_header(header, tmp, buffer);
				if (verbose >= 1) printf("c | Found header %s\n", buffer);
				cnf_stats.n_org_vars = header[0], cnf_stats.n_org_cls = header[1];
				orgs.resize(cnf_stats.n_org_cls);
				WT_alloc();
				solver_alloc();
				solver_init();
			}
			else {
				timer->start();
				char* head = buffer;
				bool c_end = false;
				while (*buffer) if (*buffer++ == ' ' && *buffer == '0') c_end = true;
				buffer = head;
				if (!c_end) { printf("Error - clause ending not recognizable.\n"); return TERMINATE; }
				c_sz = count_spaces(buffer); // predict number of lits (for faster clause storage!)
				if (c_sz == 0) { printf("Error - empty clause.\n"); return TERMINATE; }
				in_c = new uint32[c_sz];
				toLits(in_c, tmp, buffer, c_sz);
				if (proof_en) {
					oldSz = c_sz;
					tmpCl.clear(true);
					tmpCl.resize(c_sz);
					tmpCl.copyFrom(in_c, c_sz);
				}
				if (checkClause(in_c, c_sz)) { // clause not a tautology
					if (proof_en && c_sz < oldSz) {
						write_proof('a');
						write_proof(in_c, c_sz);
						write_proof(0);
						write_proof('d');
						write_proof(tmpCl, c_sz);
						write_proof(0);
					}
					if (c_sz == 1) {
						assert(*in_c > 0);
						uint32 v_idx = V2IDX(*in_c);
						if (sol->assign(v_idx) == UNDEFINED) { var_heap->remove(v_idx); enqueue(*in_c); }
						else if (sol->assign(v_idx) == ISNEG(*in_c)) return UNSAT;
					}
					else if (c_sz > MAX_CL_SZ) { 
						printf("Error - clause size \"%d\" not supported.\n", c_sz);
						return TERMINATE;
					}
					else if (cnf_stats.global_n_cls + 1 > cnf_stats.n_org_cls) {
						printf("Error - too many clauses.\n");
						return TERMINATE;
					}
					else {
						B_REF org = new BCLAUSE();
						org->copyLitsFrom(in_c, c_sz);
						attachClause(org);
					}
				}
				timer->stop();
				timer->par += timer->CPU_time();
				delete[] in_c;
			}
		}
	}
	delete[] buffer, tmp;
	assert(nClauses() + nBins() <= nOrgClauses());
	cnf_stats.n_org_bins = nBins();
	cnf_stats.n_org_lits = nLiterals() + ((int64)nBins() << 1);
	if ((int)nClauses() < orgs.size()) orgs.resize(nClauses());
	inputFile.close();
	if (verbose >= 1) printf("c | Read %d Vars, %d Cls, and %lld Lits in %.3f sec.\n", nOrgVars(), nOrgClauses(), nOrgLits() + learnt_cl.size(), timer->par);
	return UNSOLVED;
}

void ParaFROST::free_mem()
{
	occurs.clear(true);
	scores.clear(true);
	learntLits.clear(true);
	learnt_cl.clear(true);
	trail_sz.clear(true);
	wt.clear(true);
	bins.clear(true);
	orgs.clear(true);
	lbdQ.clear(true);
	trailQ.clear(true);
	scnf.clear(true);
	ot.clear(true);
	removed.clear(true);
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

void ParaFROST::WT_alloc(bool re)
{
	/* allocate host memory for the watch table */
	assert(nOrgVars() > 0);
	size_t maxVars = V2D(size_t(nOrgVars()) + 1);
	double wt_cap = (double)maxVars * sizeof(WATCH);
	assert(wt.empty());
	if (re) {
		int nOldVars = nOrgVars() + nRemVars();
		sysMemCons -= V2D(size_t(nOldVars) + 1);
		if (sysMemCons < 0) sysMemCons = 0.0;
	}
	sysMemCons += wt_cap;
	if (sysMemCons >= sysMemTot) {
		if (verbose > 1) cout << "\nc | Not enough system memory for watch table (Max: " << sysMemTot / MBYTE << ", Consumed: " << sysMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory for watch table (Max: " << sysMemTot / MBYTE << ", Consumed: " << sysMemCons / MBYTE << " MB)" << endl;
		exit(EXIT_FAILURE);
	}
	wt.resize(maxVars);
}

void ParaFROST::solver_alloc(bool re)
{
	if (verbose >= 2) cout << "c | Allocating fixed memory..";
	if (re) {
		int nOldVars = nOrgVars() + nRemVars();
		assert(sysMem != NULL);
		assert(var_heap != NULL);
		free(sysMem);
		delete var_heap;
		sysMemCons -= (sysMem_sz + sizeof(VAR_HEAP) + (double)sizeof(double) * nOldVars);
		if (sysMemCons < 0) sysMemCons = 0.0;
	}
	lrn.max_cl_sz = nOrgVars();
	sysMem_sz = sizeof(SOL) + sizeof(SP) + sizeof(PV) +
		nOrgVars() * (sizeof(uint32) * 5 + sizeof(ASSIGN_ST) + sizeof(G_REF) + 5) +
		lrn.max_cl_sz * sizeof(uint32) * 2;
	sysMemCons += sysMem_sz + sizeof(VAR_HEAP) + sizeof(double) * nOrgVars();
	if (sysMemCons >= sysMemTot) {
		cout << "| Not enough system memory for the solver (Allowed: " << sysMemTot / MBYTE << ", Consumed: " << sysMemCons / MBYTE << " MB)" << endl;
		exit(EXIT_FAILURE);
	}
	sysMem = NULL;
	sysMem = (Byte*)malloc(sysMem_sz);
	memset(sysMem, 0, sysMem_sz);
	Byte* bottom = sysMem + sysMem_sz;
	// sol
	sol = (SOL*)sysMem;
	sysMem += sizeof(SOL);
	assert(sysMem < bottom);
	sol->set_size(nOrgVars());
	sol->alloc_lits(&sysMem);
	assigns = sol->assigns_ptr();
	levels = sol->levels_ptr();
	assert(sysMem < bottom);
	// elections
	pv = (PV*)sysMem;
	sysMem += sizeof(PV);
	assert(sysMem < bottom);
	pv->PVs = (uint32*)sysMem;
	sysMem += nOrgVars() * sizeof(uint32);
	assert(sysMem < bottom);
	pv->melted = (bool*)sysMem;
	sysMem += nOrgVars();
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
	sp->free_decs = sp->trail + nOrgVars();
	sysMem += ((size_t)nOrgVars() << 1) * sizeof(uint32);
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
	var_heap = NULL;
	var_heap = new VAR_HEAP();
	var_heap->set_allocSize(nOrgVars());
	var_heap->alloc();
	if (verbose >= 2) cout << "(" << (double)sysMem_sz / MBYTE << " MB) ==> done." << endl;
}

void ParaFROST::CNF_rewriter(const string& path) {
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
	char* tmp = new char[LIT_LEN];
	char* buffer = new char[BUFFER_SIZE];
	uVector1D header(2);
	cnf_stats.global_n_cls = 0;
	register int clause_len = 0;
	cout << "c | Parsing BCNF file \"" << path << "\"" << endl;
	outFile << "c rewritten formula generated by ParaFROST." << endl;
	while (inputFile.getline(buffer, BUFFER_SIZE)) {
		std::streamsize len = inputFile.gcount();
		if (len == 0) { cout << "c | Error --> Empty line found. " << endl; exit(EXIT_FAILURE); }
		if (len == 1 && *buffer == '\0') continue;
		if (len == 2 && (*buffer == '0' || (*buffer >= '%' && *buffer <= '/'))) continue;
		if (*buffer != 'c') {
			if (*buffer == 'p') {
				read_header(header, tmp, buffer);
				cout << "c | Found header " << buffer << endl;
				cnf_stats.n_org_cls = header[1];
				outFile << buffer << endl;
			}
			else {
				char* head = buffer;
				bool clause_end = false;
				while (*buffer) if (*buffer++ == ' ' && *buffer == '0') clause_end = true;
				buffer = head;
				if (!clause_end) {
					cout << "c | Error --> clause ending not recognizable." << endl;
					cout << "c | Dump buffer --> " << buffer << endl;
					exit(EXIT_FAILURE);
				}
				clause_len = count_spaces(buffer); // predict number of lits (for faster clause storage!)
				if (clause_len == 0) {
					cout << "c | Error --> empty clause." << endl;
					exit(EXIT_FAILURE);
				}
				else if (clause_len > MAX_CL_SZ) {
					cout << "c | Error --> clause size: " << clause_len << " is unsupported." << endl;
					exit(EXIT_FAILURE);
				}
				else if (cnf_stats.global_n_cls + 1 > cnf_stats.n_org_cls) {
					cout << "c | Error --> too many clauses." << endl;
					exit(EXIT_FAILURE);
				}
				cnf_stats.global_n_cls++;
				outFile << buffer << endl;
			}
		}
	}
	delete[] buffer, tmp;
	assert(nClauses() <= nOrgClauses());
	inputFile.close();
	outFile.close();
}

void ParaFROST::attachClause(B_REF c)
{
	CL_LEN sz = c->size();
	assert(sz > 1);
	assert(*c != NULL);
	assert(**c > 0 && *(*c + 1) > 0);
	assert(**c <= UINT32_MAX && *(*c + 1) <= UINT32_MAX);
	wt[FLIP(**c)].push(WATCH(c, *(*c + 1)));
	wt[FLIP(*(*c + 1))].push(WATCH(c, **c));
	if (sz == 2) {
		c->flag_orgBin();
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

void ParaFROST::attachClause(C_REF c)
{
	CL_LEN sz = c->size();
	assert(sz > 2);
	assert(*c != NULL);
	assert(**c > 0 && *(*c + 1) > 0);
	assert(**c <= UINT32_MAX && *(*c + 1) <= UINT32_MAX);
	wt[FLIP(**c)].push(WATCH(c, *(*c + 1)));
	wt[FLIP(*(*c + 1))].push(WATCH(c, **c));
	c->set_status(LEARNT);
	c->flag_reason();
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

void ParaFROST::detachClause(B_REF c)
{
	assert(c->status() == ORIGINAL);
	assert(c->size() > 2);
	cnf_stats.global_n_lits -= c->size();
	remWatch(wt[FLIP(c->w0_lit())], c);
	remWatch(wt[FLIP(c->w1_lit())], c);
	if (proof_en) {
		write_proof('d');
		write_proof(*c, c->size());
		write_proof(0);
	}
	delete c;
}

void ParaFROST::detachClause(C_REF c)
{
	assert(c->status() == LEARNT);
	assert(c->size() > 2);
	cnf_stats.n_added_lits -= c->size();
	remWatch(wt[FLIP(c->w0_lit())], c);
	remWatch(wt[FLIP(c->w1_lit())], c);
	if (proof_en) {
		write_proof('d');
		write_proof(*c, c->size());
		write_proof(0);
	}
	delete c;
}

void ParaFROST::shrinkClause(G_REF gc)
{
	B_REF c = (B_REF)gc;
	assert(c->size() > 2);
	assert(assigns[V2IDX(c->w0_lit())] == UNDEFINED && assigns[V2IDX(c->w1_lit())] == UNDEFINED);
	uVector1D& lits = *c;
	CL_LEN sz = c->size();
	uVector1D tmpCl;
	if (proof_en) {
		tmpCl.resize(sz);
		tmpCl.copyFrom(lits);
	}
	if (sz == 3 && assigns[V2IDX(lits[2])] == ISNEG(lits[2]))
		lits.pop();
	else {
		LIT_POS k = 2;
		while (k < sz) {
			if (assigns[V2IDX(lits[k])] == ISNEG(lits[k])) {
				lits[k] = lits[sz - 1];
				sz--;
			}
			else k++;
		}
		int remLits = c->size() - sz;
		if (proof_en && remLits) {
			write_proof('a');
			write_proof(*c, sz);
			write_proof(0);
			write_proof('d');
			write_proof(tmpCl, c->size());
			write_proof(0);
		}
		lits.shrink(remLits);
	}
	assert(c->size() > 1);
}

int ParaFROST::simplify(BCNF& cnf)
{
	if (cnf.size() == 0) return 0;
	int n = 0;
	for (int i = 0; i < cnf.size(); i++) {
		B_REF c = cnf[i];
		if (c->satisfied(assigns))
			detachClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) cnf[n++] = c;
			else {
				c->flag_orgBin();
				c->set_status(UNKNOWN);
				cnf_stats.global_n_bins++;
				//printf("c | Clause became bin "); c->print();
			}
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
		if (c->satisfied(assigns))
			detachClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) cnf[n++] = c;
			else {
				c->set_status(UNKNOWN);
				cnf_stats.global_n_bins++;
				//printf("c | Clause became bin "); c->print();
			}
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
	assert(sol->assign(V2IDX(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	timer->start();
	// reduce watch table
	for (int i = sp->trail_offset; i < sp->trail_size; i++) {
		uint32 assign = sp->trail[i], assign_f = FLIP(assign);
		// remove bins from WT
		WL& ws = wt[assign_f];
		for (int j = 0; j < ws.size(); j++) {
			B_REF c = (B_REF)ws[j].c_ref;
			if (c->status() == UNKNOWN && c->size() == 2) { // original or learnt binary clause 
				//printf("c | Removing bin clause "); c->print();
				uint32 w_lit = (assign == c->w1_lit()) ? c->w0_lit() : c->w1_lit();
				remWatch(wt[FLIP(w_lit)], c);
				delete c;
				cnf_stats.global_n_bins--;
			}
		}
		// remove root-level watch lists
		wt[assign].clear(true);
		wt[assign_f].clear(true);
		uint32 assign_idx = V2IDX(assign);
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
	timer->stop();
	timer->red += timer->CPU_time();
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
	assert(sol->assign(V2IDX(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	timer->start();
	// simplify orgs
	cnf_stats.global_n_cls = 0;
	for (int i = 0; i < orgs.size(); i++) {
		B_REF c = orgs[i];
		if (c->satisfied(assigns)) detachClause(c);
		else {
			shrinkClause(c);
			if (c->size() > 2) orgs[cnf_stats.global_n_cls++] = c;
			else { c->flag_orgBin(); c->set_status(UNKNOWN); bins.push(c); }
		}
	}
	orgs.resize(cnf_stats.global_n_cls);
	// simplify bins
	cnf_stats.global_n_bins = 0;
	for (int i = 0; i < bins.size(); i++) {
		B_REF c = bins[i];
		if (c->satisfied(assigns)) {
			remWatch(wt[FLIP(c->w0_lit())], c);
			remWatch(wt[FLIP(c->w1_lit())], c);
			delete c;
		}
		else bins[cnf_stats.global_n_bins++] = c;
	}
	bins.resize(cnf_stats.global_n_bins);
	// remove root-level watch lists
	for (int i = sp->trail_offset; i < sp->trail_size; i++) {
		uint32 assign = sp->trail[i];
		wt[assign].clear(true);
		wt[FLIP(assign)].clear(true);
		uint32 assign_idx = V2IDX(assign);
		if (var_heap->has(assign_idx)) var_heap->remove(assign_idx); // remove assign from the heap
		if (pre_en) removed.push(assign); // save in case preprocessing mapped variables
	}
	cnf_stats.global_n_del_vars += (sp->trail_size - sp->trail_offset);
	sp->trail_offset = sp->trail_size;
	timer->stop();
	timer->red += timer->CPU_time();
	assert(consistent(bins, wt));
	assert(consistent(orgs, wt));
	if (verbose >= 2) printf(" ==> done\n");
}

void ParaFROST::enqueue(const uint32& assign, const int& pLevel, const G_REF src)
{
	assert(assign > 0);
	register uint32 assign_idx = V2IDX(assign);
	source[assign_idx] = src;
	sp->lock[assign_idx] = true;
	sp->trail[sp->trail_size++] = assign;
	sol->set_assign(assign_idx, !ISNEG(assign));
	sol->set_level(assign_idx, pLevel);
	if (verbose >= 3) printf("c | New %s: %d@%d.\n", src == NULL ? "decision" : "imp", (ISNEG(assign)) ? -int(ABS(assign)) : ABS(assign), pLevel);
}

void ParaFROST::simp_learnt()
{
	assert(learnt_cl.size() > 1);
	assert(lrn.max_cl_sz > learnt_cl.size());
	if (verbose >= 3) printLearnt(learnt_cl);
	register uint32* lit = learnt_cl, *tmp_lit = tmp_stack, *end = lit + learnt_cl.size();
	while (lit != end) *tmp_lit++ = *lit++;
	uint32 min_level = 0;
	for (LIT_POS k = 1; k < learnt_cl.size(); k++)  min_level |= mapHash(sol->level(V2IDX(learnt_cl[k])));
	CL_LEN new_sz = 1, tmp_sz = (CL_LEN)learnt_cl.size();
	for (LIT_POS k = 1; k < learnt_cl.size(); k++) {
		if (source[V2IDX(learnt_cl[k])] == NULL || !selfsub(learnt_cl[k], tmp_stack, tmp_sz, min_level))
			learnt_cl[new_sz++] = learnt_cl[k];
	}
	tmp_lit = tmp_stack; end = tmp_stack + tmp_sz;
	while (tmp_lit < end) sp->seen[V2IDX(*tmp_lit++)] = 0;
	learnt_cl.resize(new_sz);
}

void ParaFROST::binSelfsub() 
{
	marker++;
	for (LIT_POS i = 1; i < learnt_cl.size(); i++) // mark all learnt literals except asserting
		board[V2IDX(learnt_cl[i])] = marker;
	// check binary clauses if subsuming learnt clause
	uint32 parent = FLIP(learnt_cl[0]);
	WL& ws = wt[parent];
	int nLitsRem = 0;
	for (int i = 0; i < ws.size(); i++) {
		if (((B_REF)ws[i].c_ref)->size() == 2) {
			uint32 imp = ws[i].blocker, imp_idx = V2IDX(imp);
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
			if (board[V2IDX(learnt_cl[i])] != marker) {
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
	assert(source[V2IDX(learntLit)] != NULL);
	assert(sol->level(V2IDX(learntLit)) > ROOT_LEVEL);
	CL_LEN tmp_size = tmp_tail;
	uint32 *head = simpLearnt, *tail = simpLearnt + 1;
	*head = learntLit;
	while (head != tail) {
		if (verbose >= 4) { printf("c | simpLearnt = "); printLit(*head); }
		assert(*head > 0);
		B_REF c = (B_REF)source[V2IDX(*head)];
		assert(c != NULL);
		assert(c->w0_lit() == FLIP(*head));
		for (LIT_POS l = 1; l < c->size(); l++) {
			uint32 parent = (*c)[l], parent_idx = V2IDX(parent);
			int parent_dl = sol->level(parent_idx);
			if (!sp->seen[parent_idx] && parent_dl > 0) {
				if (source[parent_idx] != NULL && (mapHash(parent_dl) & min_level) != 0) {
					sp->seen[parent_idx] = 1;
					*tail++ = parent;
					tmp_stack[tmp_tail++] = parent;
				}
				else {
					for (uint32 *k = tmp_stack + tmp_size; k != tmp_stack + tmp_tail; k++)
						sp->seen[V2IDX(*k)] = 0;
					tmp_tail = tmp_size;
					return false;
				}
			}
		}
		head++;
	}
	return true;
}

void ParaFROST::bt_level()
{
	if (learnt_cl.size() == 1) sp->bt_level = ROOT_LEVEL;
	else if (learnt_cl.size() == 2) sp->bt_level = sol->level(V2IDX(learnt_cl[1]));
	else {
		LIT_POS max_k = 1;
		for (LIT_POS k = 2; k < learnt_cl.size(); k++) {
			if (sol->level(V2IDX(learnt_cl[k])) > sol->level(V2IDX(learnt_cl[max_k])))
				max_k = k;
		}
		uint32 max_k_lit = learnt_cl[max_k];
		learnt_cl[max_k] = learnt_cl[1];
		learnt_cl[1] = max_k_lit;
		sp->bt_level = sol->level(V2IDX(max_k_lit));
	}
}

void ParaFROST::cbt_level(C_REF c)
{
	sp->max1Found = false;
	register uint32 w0_lit = (*c)[0];
	sp->cbt_level = sol->level(V2IDX(w0_lit));
	if (DL() == sp->cbt_level && DL() == sol->level(V2IDX((*c)[1]))) return;
	// find the highest level in conflicting clause beyond literal-0
	sp->max1Found = true;
	int maxIdx = 0;
	for (LIT_POS l = 1; l < c->size(); l++) { 
		register int litLevel = sol->level(V2IDX((*c)[l]));
		if (litLevel > sp->cbt_level) {
			sp->max1Found = true;
			sp->cbt_level = litLevel;
			maxIdx = l;
		}
		else if (litLevel == sp->cbt_level && sp->max1Found == true)
			sp->max1Found = false;
	}
	if (maxIdx) { // found max. level, swap with w0-literal
		assert((*c)[0] == w0_lit);
		(*c)[0] = (*c)[maxIdx];
		(*c)[maxIdx] = w0_lit;
		if (maxIdx > 1) { // found new watch at 0-position
			assert(w0_lit == (*c)[maxIdx]);
			remWatch(wt[FLIP(w0_lit)], c);
			wt[FLIP((*c)[0])].push(WATCH(c, (*c)[1]));
		}
	}
}

void ParaFROST::analyze(C_REF c)
{
	if (verbose >= 2) printf("c | Analyzing conflict..");
	if (verbose >= 3) {
		printf("\n"); printTrail(sp->trail, sp->trail_size);
	}
	timer->start();
	learnt_cl.clear();
	learnt_cl.push(UNKNOWN);
	int conf_level = UNDEFINED;
	if (cbt_en) conf_level = sol->level(V2IDX(c->w0_lit()));
	else conf_level = DL();
	int parent = 0, track = 0, index = sp->trail_size - 1;
	do {
		assert(c != NULL);
		if (SH != 1 && c->status() == LEARNT) clBumpAct(c);
		if (SH == 2 && c->status() == LEARNT && c->LBD() > GLUE) { // try to update LBD of a previously learnt clauses 
			uint32 currLBD = calcLBD(c);
			if (currLBD + 1 < c->LBD()) {
				if (c->LBD() <= lbdFrozen) c->freeze();
				c->set_LBD(currLBD);
			}
		}
		for (LIT_POS j = (parent == 0) ? 0 : 1; j < c->size(); j++) {
			uint32 q = (*c)[j], q_idx = V2IDX(q);
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
				while (!sp->seen[V2IDX(sp->trail[index--])]);
				parent = sp->trail[index + 1];
				assert(parent > 0);
			} while (sol->level(V2IDX(parent)) < conf_level);
		}
		else {
			while (!sp->seen[V2IDX(sp->trail[index--])]);
			parent = sp->trail[index + 1];
			assert(parent > 0);
		}
		c = (C_REF)source[V2IDX(parent)];
		sp->seen[V2IDX(parent)] = 0;
		track--;
	} while (track > 0);
	assert(learnt_cl[0] == UNKNOWN);
	learnt_cl[0] = FLIP(parent);
	stats.max_lits += learnt_cl.size();
	if (learnt_cl.size() == 1) sp->bt_level = ROOT_LEVEL;
	else { // simplify learnt clause 
		simp_learnt(); 
		if (SH == 2 && learnt_cl.size() <= lbdMinClSize && calcLBD(learnt_cl) <= lbdMinReduce)
			binSelfsub();
		bt_level();
	}
	if (SH == 2) {
		sp->learnt_lbd = calcLBD(learnt_cl); // calculate LBD of minimized learnt clause
		if (verbose == 4) printf("c | LBD of learnt clause = %d\n", sp->learnt_lbd);
		if (learntLits.size() > 0) {
			for (int i = 0; i < learntLits.size(); i++) { // pump all variables having lower LBD values
				uint32 q_idx = V2IDX(learntLits[i]);
				c = (C_REF)source[q_idx];
				if (c->status() == LEARNT && c->LBD() < sp->learnt_lbd) // NOTE: LBD call must be protected by status check to gurantee a valid LBD value
					var_heap->varBumpAct(q_idx);
			}
			learntLits.clear();
		}
	}
	stats.tot_lits += learnt_cl.size();
	timer->stop();
	timer->bj += timer->CPU_time();
	if (verbose >= 3) printLearnt(learnt_cl);
}

void ParaFROST::cancel_assigns(const int& bt_level)
{
	// Non-Chrono. BT
	if (!cbt_en && DL() > bt_level) {
		for (int i = sp->trail_size - 1; i >= trail_sz[bt_level]; i--) {
			assert(sp->trail[i] > 0);
			register uint32 assign = sp->trail[i], assign_idx = V2IDX(assign);
			sol->init(assign_idx);
			if (!var_heap->has(assign_idx)) var_heap->insert(assign_idx);
			if (polarity == 1) sp->pol[assign_idx] = ISNEG(assign);
			else if (polarity == -1) sp->pol[assign_idx] = !ISNEG(assign);
			if (source[assign_idx] != NULL) {
				((B_REF)source[assign_idx])->init_reason();
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
			register uint32 assign = sp->trail[i], assign_idx = V2IDX(assign);
			if (sol->level(assign_idx) <= bt_level) 
				tmp_stack[nLitsRem++] = sp->trail[i];
			else {
				sol->init(assign_idx);
				if (!var_heap->has(assign_idx)) var_heap->insert(assign_idx);
				if (polarity == 1) sp->pol[assign_idx] = ISNEG(assign);
				else if (polarity == -1) sp->pol[assign_idx] = !ISNEG(assign);
				if (source[assign_idx] != NULL) {
					((B_REF)source[assign_idx])->init_reason();
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
	if (verbose >= 2) printf("c | Backjumping to level %d\n", bt_level);
	timer->start();
	assert(sp->trail_size > 0);
	// cancel old assignments up to backtrack level <bt_level>
	cancel_assigns(bt_level);
	// add learnt clause & enqueue learnt decision
	if (proof_en) {
		write_proof('a');
		write_proof(learnt_cl.d_ptr(), learnt_cl.size());
		write_proof(0);
	}
	if (learnt_cl.size() == 1) 
		enqueue(learnt_cl[0]);
	else if (learnt_cl.size() == 2) { // Note: LBD of binary learnts is always 2, no need to create C_REF clause
		cnf_stats.global_n_bins++; cnf_stats.global_n_gcs++;
		B_REF learnt_bin = new BCLAUSE();
		learnt_bin->copyLitsFrom(learnt_cl, learnt_cl.size());
		wt[FLIP(learnt_cl[0])].push(WATCH(learnt_bin, learnt_cl[1]));
		wt[FLIP(learnt_cl[1])].push(WATCH(learnt_bin, learnt_cl[0]));
		enqueue(learnt_cl[0], sp->bt_level, learnt_bin); // enqueue must always use sp->bt_level
	}
	else if (learnt_cl.size() > 2) {
		C_REF learnt = new LCLAUSE();
		learnt->copyLitsFrom(learnt_cl, learnt_cl.size());
		attachClause(learnt);
		enqueue(learnt_cl[0], sp->bt_level, learnt); // enqueue must always use sp->bt_level
	}
	timer->stop();
	timer->bj += timer->CPU_time();
}

void ParaFROST::remWatch(WL& ws, const G_REF gc)
{
	if (ws.size() == 0) return;     
	if (ws.size() > 1) {
		int c_idx = 0;
		while (ws[c_idx].c_ref != gc) c_idx++;
		assert(c_idx < ws.size());
		while (c_idx < ws.size() - 1) {
			ws[c_idx] = ws[c_idx + 1];
			c_idx++;
		}
	}
	ws.pop();
}

void ParaFROST::lReduce()
{
	if (verbose >= 2) printf("c | Reducing learnt database..");
	timer->start();
	/*printf("c | before reduce:\n");
	printCNF(learnts, 0);*/
	int n = 0;
	if (SH == 1) { // CSR
		Sort(learnts, LEARNT_SR());
		int pivot = learnts.size() / 2.1;
		for (int i = 0; i < learnts.size(); i++) {
			C_REF c = learnts[i];
			if (i < pivot && !c->isReason()) {
				assert(c->size() > 2);
				detachClause(c);
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
				if (!c->isReason() && (i < pivot || c->activity() < act_limit)) {
					assert(c->size() > 2);
					detachClause(c);
				}
				else learnts[n++] = c;
			}
		}
		else {
			int pivot = learnts.size() / 2;
			if (learnts[pivot]->LBD() <= 3) lrn.nClsReduce += incReduceBig;
			for (int i = 0; i < learnts.size(); i++) {
				C_REF c = learnts[i];
				if (c->LBD() > GLUE && c->molten() && !c->isReason() && i < pivot) {
					assert(c->size() > 2);
					detachClause(c);
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
	timer->stop();
	timer->red += timer->CPU_time();
	assert(consistent(learnts, wt));
	if (verbose >= 2) printf(" ==> done.\n");
	/*printf("c | after reduce:\n");
	printCNF(learnts, 0);*/
}

C_REF ParaFROST::BCP()
{
	timer->start();
	C_REF conf_ref = NULL;
	int numProps = 0;
	while (sp->trail_head < sp->trail_size) { // propagate assignments
		numProps++;
		uint32 assign = sp->trail[sp->trail_head++], assign_idx = V2IDX(assign), assign_dl = sol->level(assign_idx);
		assert(assign > 0);
		if (verbose >= 4) printf("c | Propagating assign(%d@%d):\n", (ISNEG(assign)) ? -int(ABS(assign)) : ABS(assign), assign_dl);
		if (verbose >= 4) printWatched(assign_idx);
		WL& ws = wt[assign];
		if (ws.size()) {
			WATCH* w_i = ws, * w_j = w_i, * end = w_i + ws.size();
			while (w_i != end) {
				uint32 blocker = w_i->blocker;
				assert(blocker && blocker < UINT32_MAX);
				if (sol->assign(V2IDX(blocker)) == !ISNEG(blocker)) {
					*w_j++ = *w_i++;
					continue;
				}
				// move assigned-0 literal to watch 1
				C_REF c = (C_REF)w_i->c_ref;
				//c->print();
				assert(c->size() >= 2);
				uint32 f_assign = FLIP(assign);
				if (c->w0_lit() == f_assign) c->swap_ws(); 
				assert(c->w1_lit() == f_assign);
				w_i++;
				// check if first literal is true
				uint32 w0_lit = c->w0_lit(), w0_lit_idx = V2IDX(w0_lit);
				WATCH w = WATCH(c, w0_lit);
				if (w0_lit != blocker && sol->assign(w0_lit_idx) == !ISNEG(w0_lit))
					*w_j++ = w;
				else {
					for (LIT_POS k = 2; k < c->size(); k++) { // look for (un)assigned-1 literal to watch
						register uint32 lit = (*c)[k];
						ASSIGN_ST h = sol->assign(V2IDX(lit));
						if (h == UNDEFINED || h == !ISNEG(lit)) { 
							c->set_w1(lit);
							(*c)[k] = f_assign;
							assert(c->w1_lit() == lit && (*c)[k] == f_assign);
							wt[FLIP(lit)].push(w);
							goto nextC;
						}
					}
					// clause is unit or conflict
					*w_j++ = w;
					if (sol->assign(w0_lit_idx) == UNDEFINED) {
						assert(sol->level(w0_lit_idx) == UNDEFINED);
						c->flag_reason();
						if (assign_dl == DL())
							enqueue(w0_lit, assign_dl, c);
						else {
							// find parent with max. level
							int maxLevel = assign_dl, maxIdx = 1;
							for (LIT_POS k = 2; k < c->size(); k++) { 
								register int litLevel = sol->level(V2IDX((*c)[k]));
								if (litLevel > maxLevel) { maxLevel = litLevel; maxIdx = k; }
							}
							if (maxIdx != 1) {
								c->swap(1, maxIdx);
								*w_j--; // remove (w) from assign list
								wt[FLIP((*c)[1])].push(w); // add it as new 1-literal watched 
							}
							enqueue(w0_lit, maxLevel, c);
							//printClause(*c); printf("c | max level in c = %d\n", maxLevel);
						}
					}
					else {
						if (verbose >= 3) printf("c | Conflict on (%d@%d)\n", w0_lit_idx + 1, sol->level(w0_lit_idx));
						sp->trail_head = sp->trail_size;
						conf_ref = c;
						while (w_i < end) *w_j++ = *w_i++;
					}
				}
			nextC:;
			}
			ws.shrink(w_i - w_j);
		}
		if (verbose >= 4) printWatched(assign_idx);
	}
	stats.n_props += numProps;
	lrn.simp_props -= numProps;
	timer->stop();
	timer->bcp += timer->CPU_time();
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
	while (!interrupted()) {
		if (verbose >= 2) printf("c | Current Decision Level = %d.\n", DL());
		C_REF conf_ref = BCP();
		if (conf_ref != NULL) {
			nConflicts++; confs++; 
			if (cbt_en) {
				cbt_level(conf_ref);
				if (sp->cbt_level == ROOT_LEVEL) return UNSAT;
				if (sp->max1Found) {
					cancel_assigns(sp->cbt_level - 1);
					continue;
				}
			}
			else if (DL() == ROOT_LEVEL) return UNSAT;
			if (SH == 2) {
				maxConflicts++;
				// increase var decay over time
				if (nConflicts % blockRestBase == 0 && var_heap->VarDecay() < 0.95) 
					var_heap->incVarDecay(opt_var_decay_r);
				// block restarts
				trailQ.push(sp->trail_size);
				if (maxConflicts > BL_RESTART_MIN && lbdQ.isFull() && sp->trail_size > (RB * trailQ.average())) {
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
			if (SH == 2 && lbdQ.isFull() && ((lbdQ.average() * RF) > (lbdSum / maxConflicts))) { // LBD restarts
				lbdQ.reset();
				cancel_assigns(ROOT_LEVEL);
				return UNSOLVED;
			}
			else if (SH < 2 && confs >= maxConflicts) { // luby/power restarts
				cancel_assigns(ROOT_LEVEL);
				return UNSOLVED;
			}
			// simplify (satisfied clauses)
			if (DL() == ROOT_LEVEL) simplify();
			// learnt clauses reduction
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
	if (pre_en && !pre_delay) {	preprocess(); pre_en = false; }
	PDM_init();
	CNF_STATE status = UNSOLVED;
	double rest_fact = 1.0;
	while (status == UNSOLVED) {
		if (SH < 2) {
			rest_fact = restPolicy == "luby" ? luby_seq(restart_inc, restarts) : pow(restart_inc, restarts);
			maxConflicts = rest_fact * restart_base;
		}
		// special restart for preprocessing
		if (pre_en && restarts >= 1 && restarts % pre_delay == 0) {
			preprocess(); PDM_init();
			pre_en = false; progRate = opt_progress;
		}
		status = search();
		restarts++;
		PDM_fuse();
	}
	if (proof_en) proofFile.close();
	wrapUp(status);
}

void ParaFROST::solver_init()
{
	assert(UNKNOWN == 0);
	assert(UNDEFINED < 0 && ROOT_LEVEL == 0);
	srand(seed);
	restarts = starts = UNKNOWN;
	if (SH == 2) { lbdQ.alloc(lbdRestBase); trailQ.alloc(blockRestBase); }
	cl_params.init();
	init();
}

bool ParaFROST::consistent(BCNF& cnf, WT& wt)
{
	for (int i = 0; i < cnf.size(); i++) {
		B_REF c = cnf[i];
		//c->print();
		WL& ws_0 = wt[FLIP(c->w0_lit())];
		WL& ws_1 = wt[FLIP(c->w1_lit())];
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
			printWL(ws_0);
			printWL(ws_1);
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
		WL& ws_0 = wt[FLIP(c->w0_lit())];
		WL& ws_1 = wt[FLIP(c->w1_lit())];
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
			printWL(ws_0);
			printWL(ws_1);
			return false;
		}
	}
	return true;
}

void ParaFROST::print_reports()
{
	if (perf_en) {
		printf("c |\nc |                      Solver Report                        \n");
		if (simp_perf_en) {
			printf("c | VO + OT              : %.3f sec\n", (double)timer->vo + timer->ot);
			printf("c | BVE                  : %.3f sec\n", timer->bve);
			printf("c | SUB                  : %.3f sec\n", timer->hse);
			printf("c | BCE                  : %.3f sec\n", timer->bce);
			printf("c | HRE                  : %.3f sec\n", timer->hre);
		}
		printf("c | PDM                  : %.3f sec\n", timer->pdm);
		printf("c | BCP                  : %.3f sec\n", timer->bcp);
		printf("c | Backjump             : %.3f sec\n", timer->bj);
		printf("c | CNF Reduce           : %.3f sec\n", timer->red);
		printf("c | -----------------------------------\n");
		double simp_t = 0;
		if (simp_perf_en) {
			simp_t = (double)timer->vo + timer->ot + timer->bve + timer->hse + timer->bce + timer->hre;
			printf("c | Simpifier time       : %.3f sec\n", simp_t);
		}
		double sol_t = (double)timer->pdm + timer->red + timer->bcp + timer->bj;
		printf("c | Solving time         : %.3f sec\n", sol_t + simp_t);
		printf("c | PDM calls            : %-10d\n", stats.pdm_calls);
		printf("c | Restarts             : %-10d\n", starts);
		printf("c | Blocked restarts     : %-10d\n", stats.nRestartStops);
		printf("c | Chronological BT     : %-10d\n", stats.cbt);
		printf("c | Non-chronological BT : %-10d\n", stats.ncbt);
		printf("c | Parallel Decisions   : %-10lld  (%.1f dec/sec)\n", stats.n_pds, (stats.n_pds / sol_t));
		printf("c | Follow-Up Decisions  : %-10lld  (%.1f dec/sec)\n", stats.n_fuds, (stats.n_fuds / sol_t));
		printf("c | Propagations         : %-10lld  (%.1f prop/sec)\n", stats.n_props, (stats.n_props / sol_t));
		printf("c | Conflicts            : %-10lld  (%.1f conf/sec)\n", nConflicts, (nConflicts / sol_t));
		printf("c | Conflict literals    : %-10lld  (%.2f %% deleted)\n", stats.tot_lits, (stats.max_lits - stats.tot_lits) * 100 / (double)stats.tot_lits);
		printf("c |--------------------------------------------------------------------------------------|\n");
	}
}

void ParaFROST::print_model()
{
	printf("v ");
	for (int i = 0; i < removed.size(); i++) { // print any saved assignments
		sp->frozen[V2IDX(removed[i])] = true; // trail may not be reset
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
			if (!sp->frozen[V2IDX(sp->trail[i])]) 
				printf("%d ", ISNEG(sp->trail[i]) ? -int(ABS(sp->trail[i])) : int(ABS(sp->trail[i])));
		}
	}
}

void ParaFROST::wrapUp(const CNF_STATE &status)
{
	// print results
	if (verbose >= 1) printf("c |--------------------------------------------------------------------------------------|\nc |\n");
	double simp_t = (double)timer->vo + timer->ot + timer->bve + timer->hse + timer->bce + timer->hre;
	double sol_t = (double)timer->pdm + timer->red + timer->bcp + timer->bj;
	if (status == SAT) {
		if (!quiet_en) printf("c |\n");
		printf("%s: ", path.c_str()); printf("s SATISFIABLE (time=%.3f)\n", simp_t + sol_t);
		if (!quiet_en) printf("c |\n");
		if (model_en) {
			print_model();
			if (!quiet_en) printf("\nc |\n");
		}
	}
	else if (status == UNSAT) {
		if (!quiet_en) printf("c |\n");
		printf("%s: ", path.c_str()); printf("s UNSATISFIABLE (time=%.3f)\n", simp_t + sol_t);
		if (!quiet_en) printf("c |\n");
	}
	else if (status == TERMINATE) {
		printf("%s: ", path.c_str()); printf("s UNKNOWN (time=%.3f)\n", simp_t + sol_t);
		if (!quiet_en) printf("c |\n");
	}
	if (perf_en) print_reports();
}