/***********************************************************************[pfsolve.h]
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

#ifndef __SOLVE_
#define __SOLVE_

#include "pfsimp.h"
#include "pfvec.h"
#include "pfheap.h"
#include "pfqueue.h"
#include "pfdefs.h"


class SOL
{
	LIT_ST* _assigns;
	int* _level;
	uint32 numVars;
public:
	SOL() { _assigns = NULL, _level = NULL, numVars = 0; }
	~SOL() { _assigns = NULL, _level = NULL, numVars = 0; }
	void allocMem(addr_t* _mem, const uint32& sz) {
		numVars = sz;
		assert(numVars > 0);
		_assigns = (LIT_ST*)(*_mem);
		*_mem += numVars * sizeof(LIT_ST);
		_level = (int*)(*_mem);
		*_mem += numVars * sizeof(int);
	}
	inline uint32 size() { return numVars; }
	inline uint32 size() const { return numVars; }
	inline LIT_ST* assigns_ptr() { return _assigns; }
	inline int* levels_ptr() { return _level; }
	inline void init(const uint32& idx) { _assigns[idx] = UNDEFINED, _level[idx] = UNDEFINED; }
	inline void init_assign(const uint32& idx) { _assigns[idx] = UNDEFINED; }
	inline void init_level(const uint32& idx) { _level[idx] = UNDEFINED; }
	inline void set_assign(const uint32& idx, const LIT_ST& asg) { _assigns[idx] = asg; }
	inline LIT_ST assign(const uint32& idx) { return _assigns[idx]; }
	inline int level(const uint32& idx) { return _level[idx]; }
	inline void set_level(const uint32& idx, const int& dl) { _level[idx] = dl; }

	inline void copyAssignsTo(LIT_ST* dest) {
		LIT_ST* a = _assigns, * dest_a = dest, * a_e = a + numVars;
		while (a != a_e) *dest_a++ = *a++;
	}

	inline void print_assign(const uint32& idx) {
		printf("c | var(%d@%d) = %d\n", idx + 1, _level[idx], (int)_assigns[idx]);
	}

	inline void print_sol(uint32* vars, const uint32& numVars) {
		printf("c | Assigned vars: \n");
		for (uint32 i = 0; i < numVars; i++) {
			uint32 x = vars[i];
			if (_assigns[x] != UNDEFINED)
				printf("c | var(%d@%d) = %d\n", x + 1, _level[x], (int)_assigns[x]);
		}
		printf("c |\n");
	}
};
extern LIT_ST* assigns;
extern int* levels;

class BCLAUSE {
	uint32* _lits;
	CL_LEN _sz;
	int8_t _st;
public:
	BCLAUSE() { _lits = NULL, _sz = UNKNOWN, _st = UNKNOWN; }
	explicit BCLAUSE(const CL_LEN& newSz) {
		_lits = NULL, _st = UNKNOWN;
		_sz = newSz;
		_lits = new uint32[_sz];
	}
	virtual ~BCLAUSE() { clear(true); }
	uint32& operator [] (LIT_POS i) { assert(i < _sz); return _lits[i]; }
	uint32  operator [] (LIT_POS i) const { assert(i < _sz); return _lits[i]; }
	uint32 lit(uint32 i) { return _lits[i]; }
	operator uint32* (void) { assert(_sz != 0); return _lits; }
	inline uint32* d_ptr() { return _lits; }
	inline CL_LEN size() const { return _sz; }
	inline void pop() { _sz--; }
	inline void shrink(const CL_LEN& n) { _sz -= n; }
	inline void resize(const CL_LEN& newSz) { _sz = newSz; }
	inline void set_status(const int8_t& status) { assert(status <= ST_MASK); _st = (_st & ST_RST) | status; }
	inline void markDeleted() { assert(DELETED == ST_MASK); _st |= DELETED; }
	inline int8_t status() const { return (_st & ST_MASK); }
	inline void initGarbage() { _st &= GAR_RST; }
	inline void markGarbage() { _st |= GAR_MASK; }
	inline bool garbage() const { return (_st & GAR_MASK); }
	inline void initReason() { _st &= IMP_RST; }
	inline void markReason(void) { _st |= IMP_MASK; }
	inline bool reason() const { return (_st & IMP_MASK); }
	inline void initOrgBin(void) { _st &= BIN_RST; }
	inline void markOrgBin(void) { _st |= BIN_MASK; } // set as original binary
	inline bool orgBin(void) const { return (_st & BIN_MASK); }
	inline void melt(void) { _st |= DEL_MASK; } // set for deletion
	inline void freeze() { _st &= DEL_RST; } // prevent from deletion
	inline bool molten() const { return (_st & DEL_MASK); } // is removable
	inline bool isWatchedVar(const uint32& v_idx) const { assert(_sz > 1); return (V2X(*_lits) == v_idx || V2X(*(_lits + 1)) == v_idx); }
	inline void copyLitsFrom(uVec1D& src) {
		assert(_sz <= src.size());
		uint32* d = _lits, * end = d + _sz, * s = src;
		while (d != end) *d++ = *s++;
		assert(int(d - _lits) == _sz);
	}
	inline bool satisfied(const LIT_ST* assigns) const {
		if (status() == DELETED) return true;
		uint32* h = _lits;
		uint32* e = h + _sz;
		while (h != e) {
			assert(*h > 0);
			if (assigns[V2X(*h)] == !ISNEG(*h)) return true;
			h++;
		}
		return false;
	}
	inline void swap_ws() { // swap watched literals
		uint32 lit0 = *_lits;
		*_lits = _lits[1];
		_lits[1] = lit0;
	}
	inline void swap(const LIT_POS& p1, const LIT_POS& p2) { // swap two literals
		uint32 tmp = _lits[p1];
		_lits[p1] = _lits[p2];
		_lits[p2] = tmp;
	}
	inline void clear(bool _free = false) {
		if (_free && _lits != NULL) { delete _lits; _lits = NULL; }
		_sz = UNKNOWN; _st = UNKNOWN;
	}
	void print() const {
		printf("(");
		for (LIT_POS l = 0; l < _sz; l++) {
			int lit = int(ABS(_lits[l]));
			lit = (ISNEG(_lits[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (this->status() == ORIGINAL) printf(") O, r=%d\n", reason());
		else if (this->status() == LEARNT) printf(") A, r=%d\n", reason());
		else if (this->status() == DELETED) printf(") X, r=%d\n", reason());
		else printf(") U, r=%d\n", reason());
	}
};
typedef BCLAUSE* B_REF; // ref to BCLAUSE
typedef Vec<B_REF> BCNF;

class LCLAUSE : public BCLAUSE {
	float _act;
	uint32 _lbd;
public:
	LCLAUSE() : BCLAUSE(), _lbd(UNKNOWN), _act(0.0) {}
	explicit LCLAUSE(const CL_LEN& sz) : BCLAUSE(sz), _lbd(UNKNOWN), _act(0.0) {}
	~LCLAUSE() { _act = 0.0; _lbd = UNKNOWN; }

	inline void set_LBD(const uint32& lbd) { _lbd = lbd; }
	inline uint32 LBD() const { return _lbd; }
	inline float activity() const { return _act; }
	inline void inc_act(const float& act) { _act += act; }
	inline void sca_act(const float& act) { _act *= act; }
	inline void set_act(const float& act) { _act = act; }
	void print() {
		printf("(");
		for (LIT_POS l = 0; l < size(); l++) {
			int lit = int(ABS((*this)[l]));
			lit = (ISNEG((*this)[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (this->status() == ORIGINAL) printf(") O:%d, g=%d, a=%.2f\n", reason(), LBD(), activity());
		else if (this->status() == LEARNT) printf(") A:%d, g=%d, a=%.2f\n", reason(), LBD(), activity());
		else if (this->status() == DELETED) printf(") X:%d, g=%d, a=%.2f\n", reason(), LBD(), activity());
		else printf(") U:%d, g=%d, a=%.2f\n", reason(), LBD(), activity());
	}
};
typedef LCLAUSE* C_REF; // ref to LCLAUSE
typedef Vec<C_REF> LCNF;
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

struct WATCH {
	G_REF c_ref;
	uint32  imp;
	WATCH() { c_ref = NULL; imp = 0; }
	WATCH(G_REF ref, uint32 lit) { c_ref = ref; imp = lit; }
	~WATCH() { c_ref = NULL; imp = 0; }
};
typedef Vec<WATCH> WL;
class WT {
	Vec<WL, int64> wt;
	Vec<uint32> garbage;
	Vec<bool> collected;
	uint32 nVars;
public:
	WT() { nVars = 0; }
	WL& operator[](const uint32& lit) { assert(lit > 1); return wt[lit]; }
	const WL& operator[](const uint32& lit) const { assert(lit > 1); return wt[lit]; }
	inline void allocMem(const uint32& sz) { assert(sz); nVars = sz; wt.resize((int64)V2D(nVars + 1LL)); collected.resize((int64)V2D(nVars + 1LL), 0); }
	inline int64 size() const { return wt.size(); }
	inline bool empty() const { return wt.size() == 0; }
	inline void destroy(const uint32& lit) { wt[lit].clear(true); }
	inline WL& getClean(const uint32& lit) { assert(lit > 1); if (collected[lit]) recycle(lit); return wt[lit]; }
	inline void collect(const uint32& lit) { assert(lit > 1); if (!collected[lit]) { collected[lit] = 1; garbage.push(lit); } }
	void recycle(const uint32& lit) {
		WL& ws = wt[lit];
		if (ws.empty()) { ws.clear(true); collected[lit] = 0; return; }
		WATCH* w_i, * w_j, * end = ws + ws.size();
		for (w_i = ws, w_j = ws; w_i != end; w_i++)
			if (!(((B_REF)w_i->c_ref)->garbage())) *w_j++ = *w_i;
		ws.shrink(int(w_i - w_j));
		if (ws.empty()) ws.clear(true);
		collected[lit] = 0;
	}
	void recycle() {
		for (int i = 0; i < garbage.size(); i++) if (collected[garbage[i]]) recycle(garbage[i]);
		garbage.clear();
	}
	void remWatch(const uint32& lit, const G_REF gref)
	{
		WL& ws = wt[lit];
		if (ws.size() == 0) return;
		if (ws.size() > 1) {
			int c_idx = 0;
			while (ws[c_idx].c_ref != gref) c_idx++;
			assert(c_idx < ws.size());
			while (c_idx < ws.size() - 1) { ws[c_idx] = ws[c_idx + 1]; c_idx++; }
		}
		ws.pop();
	}
	void print(const uint32& lit) {
		WL& ws = wt[lit];
		for (int i = 0; i < ws.size(); i++) {
			printf("c | ");
			if (((B_REF)ws[i].c_ref)->status() == LEARNT) {
				C_REF cl = (C_REF)ws[i].c_ref;
				cl->print();
			}
			else
				((B_REF)ws[i].c_ref)->print();
		}
	}
	void print(void) {
		printf("c | Negative occurs:\n");
		for (uint32 v = 0; v < nVars; v++) {
			uint32 p = V2D(v + 1);
			if (wt[p].size() != 0) printf("c | Var(%d):\n", -int(v + 1)), print(p);
		}
		printf("c | Positive occurs:\n");
		for (uint32 v = 0; v < nVars; v++) {
			uint32 n = NEG(V2D(v + 1));
			if (wt[n].size() != 0) printf("c | Var(%d):\n", v + 1), print(n);
		}
	}
	void clear(bool _free = true) {
		wt.clear(_free);
		collected.clear(_free);
		garbage.clear(_free);
		nVars = 0;
	}
};

class GC {
	Vec<G_REF> garbage;
public:
	GC() {}
	inline void init(const int& cap = KBYTE) { garbage.reserve(cap, 0); }
	inline void clear(bool _free = false) { garbage.clear(_free); }
	inline int size() const { return garbage.size(); }
	inline bool empty() const { return garbage.size() == 0; }
	inline void collect(G_REF c) { garbage.push(c); }
	inline void recycle() {
		for (int i = 0; i < garbage.size(); i++) {
			assert(garbage[i] != NULL);
			if (((B_REF)garbage[i])->status() != LEARNT) delete (B_REF)garbage[i];
			else delete (C_REF)garbage[i];
			garbage[i] = NULL;
		}
		// start clean
		garbage.clear(true);
		init();
	}
	void print() {
		if (empty()) { printf("c | Garbage is empty\n"); return; }
		printf("c |\t\tgarbage (size = %d):\n", size());
		printf("c | %11s\t %10s\n", "address", "clause");
		for (int i = 0; i < size(); i++) {
			printf("c | %11llx->", (int64)garbage[i]);
			if (((B_REF)garbage[i])->status() != LEARNT) ((B_REF)garbage[i])->print();
			else ((C_REF)garbage[i])->print();
		}
		printf("c |--------------------------------------------------------------------------------------|\n");
	}
};

struct SP {
	uint32* trail;
	bool* lock, * seen, * frozen, * pol;
	int bt_level, cbt_level;
	int trail_head, trail_size, trail_offset;
	uint32 learnt_lbd;
	bool max1Found;
	void reset_trail() {
		trail_head = trail_size = UNKNOWN;
		trail_offset = UNKNOWN;
	}
	void reset_level() {
		bt_level = ROOT_LEVEL;
		cbt_level = UNDEFINED;
		max1Found = false;
	}
};
struct STATS {
	int64 n_units;
	int64 n_props;
	int64 n_fuds;
	int64 n_pds;
	int64 tot_lits, max_lits;
	int pdm_calls;
	int nRestartStops;
	int ncbt, cbt;
	void reset() {
		pdm_calls = UNKNOWN;
		n_pds = UNKNOWN;
		n_props = UNKNOWN;
		n_units = UNKNOWN;
		n_fuds = UNKNOWN;
		max_lits = UNKNOWN;
		tot_lits = UNKNOWN;
		nRestartStops = UNKNOWN;
		ncbt = cbt = UNKNOWN;
	}
};
struct CL_PARAM {
	float cl_inc, cl_decay;
	void init() {
		cl_inc = 1.0;
		cl_decay = float(0.999);
	}
};
struct LEARN {
	int64 simp_props;
	int64 max_learnt_cls;
	int64 nClsReduce;
	int max_cl_sz, adjust_cnt;
	double size_factor_cls;
	double adjust_start_conf, adjust_conf;
	double size_inc, adjust_inc;
	void init(void) {
		size_inc = 1.1;
		adjust_inc = 1.5;
		adjust_start_conf = 100;
		size_factor_cls = 1.0 / 3.1;
		simp_props = UNKNOWN;
		adjust_conf = adjust_start_conf;
		adjust_cnt = (uint32)adjust_conf;
	}
};

class cuMM {
	addr_t gMemPV, gMemSol, gMemVOrd, gMemStats; // fixed pools
	addr_t gMemCNF, gMemOT;  // dynamic pools
	size_t gMemPV_sz, gMemSol_sz, gMemVOrd_sz, gMemStats_sz, gMemCNF_sz, gMemOT_sz;
	size_t cap;
public:
	cuMM() { 
		gMemPV = NULL, gMemSol = NULL, gMemVOrd = NULL, gMemStats = NULL, gMemCNF = NULL, gMemOT = NULL;
		gMemPV_sz = 0ULL, gMemSol_sz = 0ULL, gMemVOrd_sz = 0ULL, gMemStats_sz = 0ULL;
		gMemCNF_sz = 0ULL, gMemOT_sz = 0ULL, cap = 0ULL;
	}
	~cuMM() { }
	bool allocPV(PV*, const uint32&);
	bool allocVO(OCCUR*&, SCORE*&, const uint32&);
	bool allocStats(GSTATS*&, const uint32&);
	bool resizeCNF(CNF*&, const uint32&, const uint64&, const int& phase = 0);
	OT* resizeOT(uint32*, const uint32&, const int64&);
	_PFROST_H_D_ bool empty(void) const { return cap == 0; }
	_PFROST_H_D_ size_t capacity(void) const { return cap; }
};

class ParaFROST {
protected:
	string path;
	TIMER* timer;
	addr_t sysMem; 
	size_t sysMem_sz;
	int64 sysMemAvail, sysMemCons;
	Vec<OCCUR> occurs;
	Vec<SCORE> scores;
	BCNF orgs, bins;
	LCNF learnts;
	WT wt;
	GC gcr;
	SP* sp; 
	SOL* sol;
	int* board;
	uint32* tmp_stack, * simpLearnt;
	uVec1D learnt_cl, learntLits;
	Vec1D trail_sz;
	G_REF* source;
	CL_PARAM cl_params;
	VAR_HEAP* var_heap;
	LEARN lrn;
	STATS stats;
	int64 maxConflicts, nConflicts;
	int starts, restarts, R, ref_vars, reductions, marker;
	BQUEUE<uint32> lbdQ, trailQ;
	double lbdSum;
	std::ofstream proofFile;
	bool intr, units_f;
public:
	/**********************/
	/*   Solver methods   */
	/**********************/
	inline void interrupt(void) { intr = true; }
	inline bool interrupted(void) { return intr; }
	inline void write_proof(const Byte& byte) { proofFile << byte; }
	inline double drand(void) const { return ((double)rand() / (double)RAND_MAX); };
	inline void incDL(void) { trail_sz.push(sp->trail_size); }
	inline int DL(void) const { return trail_sz.size(); }
	inline void clDecayAct() { cl_params.cl_inc *= (1 / cl_params.cl_decay); }
	inline int64 cnfSize() { return (int64)nBins() + nClauses() + learnts.size(); }
	inline void collectClause(B_REF& c, const bool& gc = true) {
		if (gc) { assert(!c->garbage()); c->markGarbage(); gcr.collect(c); }
		else { assert(c != NULL); delete c; c = NULL; }
	}
	inline void collectClause(C_REF& c, const bool& gc = true) {
		if (gc) { assert(!c->garbage()); c->markGarbage(); gcr.collect(c); }
		else { assert(c != NULL); delete c; c = NULL; }
	}
	inline uint32 calcLBD(uVec1D& c) {
		marker++;
		register uint32 lbd = 0;
		for (LIT_POS i = 0; i < c.size(); i++) {
			int litLevel = sol->level(V2X(c[i]));
			if (board[litLevel] != marker) { board[litLevel] = marker; lbd++; }
		}
		return lbd;
	}
	inline uint32 calcLBD(C_REF& c) {
		marker++;
		register uint32 lbd = 0;
		for (LIT_POS i = 0; i < c->size(); i++) {
			int litLevel = sol->level(V2X((*c)[i]));
			if (board[litLevel] != marker) { board[litLevel] = marker; lbd++; }
		}
		return lbd;
	}
	inline void PDM_fuse() {
		if ((!pdm_freq && SH == 2 && (starts % ((maxConflicts / 1000) + 1)) == 0) ||
			(!pdm_freq && SH < 2 && starts % maxConflicts == 0) ||
			(pdm_freq && starts % pdm_freq == 0)) {
			ref_vars = UNKNOWN;
			R = pdm_rounds;
		}
	}
	inline void clHist(const G_REF& gref) {
		B_REF c = (B_REF)gref;
		for (LIT_POS i = 0; i < c->size(); i++) {
			assert((*c)[i] > 0);
			if (ISNEG((*c)[i])) occurs[V2X((*c)[i])].ns++;
			else occurs[V2X((*c)[i])].ps++;
		}
	}
	inline void clBumpAct(C_REF& c) {
		c->inc_act(cl_params.cl_inc);
		if (c->activity() > (float)1e30) {
			for (int i = 0; i < learnts.size(); i++) learnts[i]->sca_act((float)1e-30);
			cl_params.cl_inc *= (float)1e-30;
		}
	}
	inline double luby_seq(double y, int x) {
		int size, seq;
		for (size = 1, seq = 0; size < x + 1; seq++, size = (size << 1) + 1);
		while (size - 1 != x) {
			size = (size - 1) >> 1;
			seq--;
			x = x % size;
		}
		return pow(y, seq);
	}
	inline void printWatched(const uint32& v_idx) {
		printf("c |\t\tWL of v(%d)\n", v_idx + 1);
		uint32 p = V2D(v_idx + 1);
		wt.print(p), wt.print(NEG(p));
		printf("c |---------------------------------\n");
	}
	inline void printStats() {
		if (verbose == 1) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			printf("c | %9d %9d %10lld | %10lld %10lld %10lld %10lld %7.0f |\n",
				remainedVars, nClauses(), nLiterals(),
				nConflicts, lrn.max_learnt_cls,
				(int64)learnts.size(), nLearntLits(), (float)nLearntLits() / learnts.size());
		}
	}
	inline void printStats(bool p) {
		if (verbose == 1 && p) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			int l2c = learnts.size() == 0 ? 0 : int(nLearntLits() / learnts.size());
			printf("c | %9d %9d %9d %10lld %6d %9d %8d %10lld %6d |\n",
				remainedVars, nClauses(), nBins(), nLiterals(),
				starts,
				learnts.size(), nGlues(), nLearntLits(), l2c);
			progRate += (progRate >> 2);
		}
	}
	inline void write_proof(uint32* lits, const CL_LEN& len) {
		assert(len > 0);
		uint32* lit = lits, * end = lits + len;
		while (lit != end) {
			register uint32 b = 0;
			if (mapped) b = revLit(*lit++);
			else b = *lit++;
			while (b > 127) { write_proof(Byte(128 | (b & 127))); b >>= 7; }
			write_proof(Byte(b));
		}
	}
	//=============================================
	ParaFROST(const string&);
	~ParaFROST(void);
	int64 sysMemUsed(void);
	void sysFree(void);
	void killSolver(const CNF_STATE& status = TERMINATE);
	void resetSolver(const bool& re = true);
	void allocSolver(const bool& re = false);
	void initSolver(void);
	void recycle(void);
	void cnfrewriter(const string&);
	void varOrder(void);
	void hist(BCNF&, const bool& rst = false);
	void hist(LCNF&, const bool& rst = false);
	void attachClause(B_REF&);
	void attachClause(C_REF&);
	void shrinkClause(G_REF);
	void removeClause(B_REF&, const bool& gc = true);
	void removeClause(C_REF&, const bool& gc = true);
	void detachClause(B_REF&, const bool& gc = true);
	int simplify(BCNF&);
	int simplify(LCNF&);
	void simplify(void);
	void simplify_top(void);
	void solve(void);
	CNF_STATE parser(const string&);
	CNF_STATE search(void);
	CNF_STATE decide(void);
	void pumpFrozen(void);
	void PDMInit(void);
	void PDM(void);
	bool depFreeze_init(WL&, const uint32&);
	bool depFreeze(WL&, const uint32&);
	void eligibility(void);
	G_REF BCP(void);
	void analyze(G_REF&);
	void cbtLevel(G_REF&);
	void enqueue(const uint32&, const int& pLevel = ROOT_LEVEL, const G_REF = NULL);
	void btLevel(void);
	void binSelfsub(void);
	void selfsub(void);
	bool selfsub(const uint32&, uint32*, CL_LEN&, const uint32&);
	void cancelAssigns(const int&);
	void backJump(const int&);
	void lReduce(void);
	bool consistent(BCNF&, WT&);
	bool consistent(LCNF&, WT&);
	void printReport(void);
	void printModel(void);
	void wrapUp(const CNF_STATE&);
	/*****************************/
	/*  SIGmA Simplifier methods */
	/*****************************/
protected:
	cuMM cuMem;
	CNF* cnf;
	OT* ot;
	PV* pv;
	OCCUR* d_occurs;
	SCORE* d_scores;
	uTHRVector histogram, rawLits;
	uint32* raw_hist;
	GSTATS* gstats;
	cudaStream_t* streams;
	bool mapped;
	uVec1D removed;
	Vec1D mappedVars, reverseVars;
	std::ofstream outputFile;
	int devCount;
public:
	inline void createStreams(void) { 
		assert(streams == NULL);
		streams = new cudaStream_t[nstreams];
		for (int i = 0; i < nstreams; i++) cudaStreamCreate(streams + i); 
	}
	inline void destroyStreams(void) { 
		for (int i = 0; i < nstreams; i++) cudaStreamDestroy(streams[i]); 
		if (streams != NULL) delete[] streams;
	}
	inline uint32 mapLit(const uint32& orgLit) {
		static int maxVar = 1;
		assert(!mappedVars.empty());
		register int orgVar = ABS(orgLit);
		assert(!pv->sol->value[orgVar - 1]);
		assert(orgVar > 0 && orgVar < mappedVars.size());
		if (mappedVars[orgVar] == 0) { reverseVars[maxVar] = orgVar; mappedVars[orgVar] = maxVar++; }
		assert(maxVar - 1 <= (int)nOrgVars());
		return (V2D(mappedVars[orgVar]) | ISNEG(orgLit));
	}
	inline uint32 revLit(const uint32& mapLit) { return (V2D(reverseVars[ABS(mapLit)]) | ISNEG(mapLit)); }
	inline bool verifyLCVE(void) { 
		for (uint32 v = 0; v < pv->numPVs; v++) if (sp->frozen[pv->pVars->at(v)]) return false;
		return true;
	}
	inline void logReductions(void) {
		printf("c |  Parallel variables = %d\n", pv->numPVs);
		printf("c |  Variables after = %d (-%d)\n", nOrgVars() - cnf_stats.n_del_vars, cnf_stats.n_del_vars);
		printf("c |  Clauses after = %d (-%d)\n", cnf_stats.n_cls_after, nClauses() - cnf_stats.n_cls_after);
		if (cnf_stats.n_lits_after > nLiterals()) printf("c |  Literals after = %lld (+%lld)\n", cnf_stats.n_lits_after, cnf_stats.n_lits_after - nLiterals());
		else printf("c |  Literals after = %lld (-%lld)\n", cnf_stats.n_lits_after, nLiterals() - cnf_stats.n_lits_after);
	}
	inline void printPStats(void) {
		if (verbose == 1 && SH == 2) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			printf("c |p %8d %9d %9s %10lld %6d %9s %8s %10s %6s |\n",
				remainedVars, nClauses(), "---", nLiterals(), starts, "---", "---", "---", "---");
		}
		if (verbose == 1 && SH < 2) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			printf("c |p %8d %9d %10lld | %10lld %10s %10s %10s %7s |\n",
				remainedVars, nClauses(), nLiterals(), nConflicts, "---", "---", "---", "---");
		}
	}
	//========================================
	void free_master(void);
	void free_slaves(void);
	void optSimp(void);
	void extractBins(void);
	bool awaken(void);
	void _simplify(void);
	void depFreeze(const OL&, const uint32&, const uint32&, const uint32&);
	bool LCVE(void);
	void GPU_hist(void);
	void GPU_occurs(const int64&, const bool& init = false);
	void GPU_VO(void);
	void GPU_CNF_fill(void);
	CNF_STATE GPU_VE(void);
	void GPU_SUB(void);
	void GPU_HRE(void);
	void GPU_BCE(void);
	void GPU_preprocess(void);
	// options
	bool quiet_en, parse_only_en, rewriter_en, perf_en;
	bool mcv_en, model_en, proof_en, fdp_en, cbt_en, csr_en;
	int progRate, initProgRate, verbose, seed, timeout;
	int pdm_rounds, pdm_freq, pdm_order;
	int SH, polarity, restart_base, blRestMin, lbdRestBase, blockRestBase, VSIDSDecayFreq;
	int nClsReduce, lbdFrozen, lbdMinReduce, lbdMinClSize, incReduceSmall, incReduceBig;
	int cbt_dist, cbt_conf_max;
	double var_inc, var_decay;
	double RF, RB, restart_inc, gperc;
	string restPolicy, proof_path;
	bool p_cnf_en, p_ot_en;
	bool pre_en, lpre_en, ve_en, ve_plus_en, res_en, subst_en, sub_en, bce_en, hre_en, all_en;
	int pre_delay, mu_pos, mu_neg, phases, cnf_free_freq, ngpus, nstreams;
};
/******************/
/*     Aliases    */
/******************/
extern ParaFROST* gpfrost; // global solver
extern uint32 verb;
//====================================================//
//                  HELPER Functions                  //
//====================================================//
inline void printLit(const uint32& l) { printf("%d@%d", (ISNEG(l)) ? -int(ABS(l)) : ABS(l), levels[ABS(l) - 1]); }
inline void printVars(const uint32* arr, const int& size)
{
	printf("c | [");
	for (int i = 0; i < size; i++) {
		printf("%4d ", arr[i] + 1);
	}
	printf("]\n");
}
inline void printAssigns(const uint32* arr, const int& size)
{
	printf("c | [");
	for (int i = 0; i < size; i++) {
		printf("%4d ", ISNEG(arr[i]) ? -int(ABS(arr[i])) : ABS(arr[i]));
	}
	printf("]\n");
}
inline void printTrail(const uint32* trail, const int& size)
{
	for (int i = 0; i < size; i++)
		printf("c | trail[%d] = %d@%d\n", i, ISNEG(trail[i]) ? -int(ABS(trail[i])) : ABS(trail[i]), levels[V2X(trail[i])]);
}
inline void printClause(const uVec1D& cl)
{
	printf("c | (");
	for (int i = 0; i < cl.size(); i++) {
		printf("%4d ", ISNEG(cl[i]) ? -int(ABS(cl[i])) : ABS(cl[i]));
	}
	printf(")\n");
}
inline void printLearnt(const uVec1D& cl)
{
	printf("c | Learnt(");
	for (LIT_POS n = 0; n < cl.size(); n++) {
		printf(" %d@%d ", (ISNEG(cl[n])) ? -int(ABS(cl[n])) : ABS(cl[n]), levels[V2X(cl[n])]);
	}
	printf(")\n");
}
template<class T>
inline void printCNF(const Vec<T>& cnf, const int& offset = 0) {
	printf("c |               Host CNF(sz = %d)\n", cnf.size());
	for (int c = offset; c < cnf.size(); c++) {
		if (cnf[c]->size() > 0) {
			printf("c | C(%d)->", c);
			cnf[c]->print();
		}
	}
	printf("c |=================================================|\n");
}

#endif // !__CUSAT_SIMP_