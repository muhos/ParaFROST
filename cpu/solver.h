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

#ifndef __SOLVE_
#define __SOLVE_

#include "defs.h"
#include "args.h"
#include "Heap.h"
#include "BQueue.h"

//====================================================//
//       Solver Data structures primitives            //
//====================================================//

/*****************************************************/
/*  Name:     SOL                                    */
/*  Usage:    store assignment info.                 */
/*  Scope:    host and device                        */
/*  Memory:   managed memory (initially on host)     */
/*  Dependency:  none                                */
/*****************************************************/
class SOL
{
	ASSIGN_ST* _assigns;
	int* _level;
	uint32 numVars;
public:
	SOL()
	{
		_assigns = NULL;
		_level = NULL;
		numVars = 0;
	}

	~SOL() {
		_assigns = NULL;
		_level = NULL;
		numVars = 0;
	}

	void alloc_lits(Byte** _mem) {
		assert(numVars > 0);
		_assigns = (ASSIGN_ST*)(*_mem);
		*_mem += numVars * sizeof(ASSIGN_ST);
		_level = (int*)(*_mem);
		*_mem += numVars * sizeof(int);
	}

	inline void set_size(const uint32& numVars) {
		this->numVars = numVars;
	}

	inline uint32 size() {
		return numVars;
	}

	inline ASSIGN_ST* assigns_ptr() {
		return _assigns;
	}

	inline void copyAssignsTo(ASSIGN_ST* dest) {
		ASSIGN_ST* a = _assigns, *dest_a = dest, * a_e = a + numVars;
		while (a != a_e) *dest_a++ = *a++;
	}

	inline int* levels_ptr() {
		return _level;
	}

	inline void init(const uint32& idx) {
		*(_assigns + idx) = UNDEFINED;
		*(_level + idx) = UNDEFINED;
	}

	inline void init_assign(const uint32& idx) {
		*(_assigns + idx) = UNDEFINED;
	}

	inline void init_level(const uint32& idx) {
		*(_level + idx) = UNDEFINED;
	}

	inline void set_assign(const uint32& idx, const ASSIGN_ST& asg) {
		*(_assigns + idx) = asg;
	}

	inline ASSIGN_ST assign(const uint32& idx) {
		return *(_assigns + idx);
	}

	inline int level(const uint32& idx) {
		return *(_level + idx);
	}

	inline void set_level(const uint32& idx, const int& dl) {
		*(_level + idx) = dl;
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
extern ASSIGN_ST* assigns;
extern int* levels;
/*****************************************************/
/*  Name:     BCLAUSE                                */
/*  Usage:    abstract clause for solving on host    */
/*  Dependency:  none                                */
/*****************************************************/
class BCLAUSE {
	uVector1D _lits;
	int8_t _st;
public:
	BCLAUSE() { _st = UNKNOWN; }
	BCLAUSE(const LIT_POS& nLits) { 
		_lits.incMem(nLits);
		_st = UNKNOWN;
	}
	~BCLAUSE() {
		_lits.clear(true);
		_st = UNKNOWN;
	}
	uint32& operator [] (LIT_POS i) { assert(i < _lits.size()); return _lits[i]; }
	uint32  operator [] (LIT_POS i) const { assert(i < _lits.size()); return _lits[i]; }
	operator uint32* (void) { assert(_lits.size() != 0); return _lits; }
	operator uVector1D& (void) { assert(_lits.size() != 0); return _lits; }

	inline uint32* d_ptr() { return _lits; }

	inline void copyLitsFrom(uint32* src, const CL_LEN& sz) {
		_lits.incMem(sz);
		register uint32* h = _lits, * e = h + sz, * s = src;
		while (h != e) *h++ = *s++;
	}

	inline void push(const uint32& lit) { _lits.push(lit); }

	inline void pop() { _lits.pop(); }

	inline void resize(const CL_LEN& newSz) { _lits.resize(newSz); }

	inline void clear(bool dealloc = false) { _lits.clear(dealloc); _st = UNKNOWN; }

	inline bool satisfied(const ASSIGN_ST* assigns) {
		if (status() == DELETED) return true;
		uint32* h = _lits;
		uint32* e = h + _lits.size();
		while (h != e) {
			assert(*h > 0);
			if (assigns[V2IDX(*h)] == !ISNEG(*h)) return true;
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

	inline void reset() {
		_lits.clear();
		_st = UNKNOWN;
	}

	inline void set_w1(const uint32& lit) {
		assert(_lits.size() > 1);
		_lits[1] = lit;
	}

	inline uint32 w1_lit() {
		assert(_lits.size() > 1);
		return _lits[1];
	}

	inline uint32 w0_lit() {
		assert(_lits.size() > 0);
		return *_lits;
	}

	inline void print() {
		printf("(");
		for (LIT_POS l = 0; l < _lits.size(); l++) {
			int lit = int(ABS(_lits[l]));
			lit = (ISNEG(_lits[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (this->status() == ORIGINAL) printf(") O, r=%d\n", isReason());
		else if (this->status() == LEARNT) printf(") A, r=%d\n", isReason());
		else if (this->status() == DELETED) printf(") X, r=%d\n", isReason());
		else printf(") U, r=%d\n", isReason());
	}

	inline void set_status(const int8_t& status) { assert(status <= ST_MASK); _st = (_st & ST_RST) | status; }

	inline int8_t status() { return (_st & ST_MASK); }

	inline int8_t status() const { return (_st & ST_MASK); }

	inline void init_reason() { _st &= IMP_RST; }

	inline void flag_reason(void) { _st |= IMP_MASK; }

	inline bool isReason() { return (_st & IMP_MASK); }

	inline void reset_orgBin(void) { _st &= BIN_RST; }

	inline void flag_orgBin(void) { _st |= BIN_MASK; } // set as original binary

	inline bool orgBin(void) { return (_st & BIN_MASK); }

	inline void melt(void) { _st |= DEL_MASK; } // set for deletion

	inline void freeze() { _st &= DEL_RST; } // prevent from deletion

	inline bool molten() { return (_st & DEL_MASK); } // is removable

	inline CL_LEN size() const { return (CL_LEN)_lits.size(); }

	inline bool isWatched(const uint32& lit) {
		assert(_lits.size() > 1);
		return (*_lits == lit || *(_lits + 1) == lit);
	}

	inline bool isWatchedVar(const uint32& v_idx) {
		return (V2IDX(*_lits) == v_idx || V2IDX(*(_lits + 1)) == v_idx);
	}
};
typedef BCLAUSE* B_REF; // ref to BCLAUSE
/*****************************************************/
/*  Name:     LCLAUSE                                */
/*  Usage:    Extended BCLAUSE with clause heuristics*/
/*  Dependency:  BCLAUSE                             */
/*****************************************************/
class LCLAUSE: public BCLAUSE {
	float _act;
	uint32 _lbd;
public:
	LCLAUSE() {
		BCLAUSE();
		_act = 0.0;
		_lbd = UNKNOWN;
	}
	~LCLAUSE() {
		clear(true);
		_act = 0.0;
		_lbd = 0;
	}

	inline void set_LBD(const uint32& lbd) { _lbd = lbd; }

	inline uint32 LBD() { return _lbd; }

	inline void init_act() {
		_act = 0.0;
	}

	inline void inc_act(const float& act) {
		_act += act;
	}

	inline void sca_act(const float& act) {
		_act *= act;
	}

	inline void set_act(const float& act) {
		_act = act;
	}

	inline float activity() {
		return _act;
	}

	inline void print() {
		printf("(");
		for (LIT_POS l = 0; l < size(); l++) {
			int lit = int(ABS((*this)[l]));
			lit = (ISNEG((*this)[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (this->status() == ORIGINAL) printf(") O:%d, g=%d, a=%.2f\n", isReason(), LBD(), activity());
		else if (this->status() == LEARNT) printf(") A:%d, g=%d, a=%.2f\n", isReason(), LBD(), activity());
		else if (this->status() == DELETED) printf(") X:%d, g=%d, a=%.2f\n", isReason(), LBD(), activity());
		else printf(") U:%d, g=%d, a=%.2f\n", isReason(), LBD(), activity());
	}
};
typedef LCLAUSE* C_REF; // ref to LCLAUSE
/*****************************************************/
/*  Name:     SCLAUSE                                */
/*  Usage:    Extended BCLAUSE with signature        */
/*  Dependency:  BCLAUSE                             */
/*****************************************************/
class SCLAUSE : public BCLAUSE {
	uint32 _sig;
public:
	SCLAUSE() {
		BCLAUSE();
		_sig = UNKNOWN;
	}
	~SCLAUSE() {
		clear(true);
		_sig = UNKNOWN;
	}

	inline uint32 lit(const LIT_POS& idx) { assert(idx < size()); return (*this)[idx]; }

	inline bool has(const uint32& lit) { // binary search
		if (size() == 2) {
			if (this->lit(0) == lit || this->lit(1) == lit) return true;
			else return false;
		}
		else {
			assert(this->isSorted());
			LIT_POS low = 0, high = size() - 1, mid;
			uint32 first = this->lit(low), last = this->lit(high);
			while (first <= lit && last >= lit) {
				mid = (low + high) >> 1;
				uint32 m = this->lit(mid);
				if (m < lit) first = this->lit(low = mid + 1);
				else if (m > lit) last = this->lit(high = mid - 1);
				else return true; // found
			}
			if (this->lit(low) == lit) return true; // found
			else return false; // Not found
		}
	}

	inline bool isSorted() {
		for (LIT_POS i = 0; i < size(); i++) {
			if (i > 0 && this->lit(i) < this->lit(i - 1)) return false;
		}
		return true;
	}

	inline LIT_POS hasZero() {
		for (LIT_POS l = 0; l < size(); l++) {
			if (this->lit(l) == 0) return l;
		}
		return -1;
	}

	inline void set_sig(const uint32& sig) { _sig = sig; }

	inline void calcSig(const uint32& init_sig = 0) {
		_sig = init_sig;
		for (LIT_POS l = 0; l < this->size(); l++)
			_sig |= (1UL << HASH((*this)[l]));
	}

	inline uint32 sig() { return _sig; }

	inline void print() {
		printf("(");
		for (LIT_POS l = 0; l < size(); l++) {
			int lit = int(ABS((*this)[l]));
			lit = (ISNEG((*this)[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (this->status() == ORIGINAL) printf(") O, s=0x%X\n", sig());
		else if (this->status() == LEARNT) printf(") A, s=0x%X\n", sig());
		else if (this->status() == DELETED) printf(") X, s=0x%X\n", sig());
		else printf(") U, s=0x%X\n", sig());
	}
};
typedef SCLAUSE* S_REF; // ref to SCLAUSE
/*****************************************************/
/*  Name:     WATCH                                  */
/*  Usage:    stores literal & pointer to clauses    */
/*  Dependency:  none                                */
/*****************************************************/
struct WATCH {
	G_REF c_ref;
	uint32  blocker;
	WATCH() { c_ref = NULL; blocker = 0; }
	WATCH(G_REF ref, uint32 lit) { c_ref = ref; blocker = lit; }
	~WATCH(){ c_ref = NULL; blocker = 0; }
};
/*****************************************************/
// define types for new structures
typedef Vec<B_REF> BCNF;
typedef Vec<C_REF> LCNF;
typedef Vec<WATCH> WL;
typedef Vec<WL> WT;
typedef Vec<S_REF> SCNF;
typedef Vec<S_REF> OL;
typedef Vec<OL> OT;
typedef struct
{
	uint32* PVs;
	bool* melted;
	int nPVs, mu_inc;
}PV;
struct SP {
	uint32* trail, * free_decs;
	bool* lock, * seen, * frozen, * pol;
	int bt_level, cbt_level;
	int trail_head, trail_size, trail_offset, numFree;
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
		max_learnt_cls = nClauses() * size_factor_cls;
		adjust_conf = adjust_start_conf;
		adjust_cnt = (uint32)adjust_conf;
	}
};

//====================================================//
//                  HELPER Functions                  //
//====================================================//
inline void printLit(const uint32& lit)
{
	printf("%d@", (ISNEG(lit)) ? -int(ABS(lit)) : ABS(lit));
}
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
		printf("c | trail[%d] = %d@%d\n", i, ISNEG(trail[i]) ? -int(ABS(trail[i])) : ABS(trail[i]), levels[V2IDX(trail[i])]);
}
inline void printClause(const uVector1D& cl)
{
	printf("c | (");
	for (int i = 0; i < cl.size(); i++) {
		printf("%4d ", ISNEG(cl[i]) ? -int(ABS(cl[i])) : ABS(cl[i]));
	}
	printf(")\n");
}
inline void printLearnt(const uVector1D& cl)
{
	printf("c | Learnt(");
	for (LIT_POS n = 0; n < cl.size(); n++) {
		printf(" %d@%d ", (ISNEG(cl[n])) ? -int(ABS(cl[n])) : ABS(cl[n]), levels[V2IDX(cl[n])]);
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
inline void printWL(const WL& list) {
	for (int i = 0; i < list.size(); i++) {
		printf("c | ");
		if (((B_REF)list[i].c_ref)->status() == LEARNT) {
			C_REF cl = (C_REF)list[i].c_ref;
			cl->print();
		}
		else
			((B_REF)list[i].c_ref)->print();
	}
}
inline void printWT(const WT& wt) {
	printf("c | Negative occurs:\n");
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 p = V2D(v + 1);
		if (wt[p].size() != 0) {
			printf("c | Var(%d):\n", -int(v + 1));
			printWL(wt[p]);
		}
	}
	printf("c | Positive occurs:\n");
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 n = NEG(V2D(v + 1));
		if (wt[n].size() != 0) {
			printf("c | Var(%d):\n", v + 1);
			printWL(wt[n]);
		}
	}
}
inline void printOL(const OL& list) {
	for (int i = 0; i < list.size(); i++) {
		printf("c | ");
		list[i]->print();
	}
}
inline void printOT(const OT& ot) {
	printf("c | Positive occurs:\n");
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 p = V2D(v + 1);
		if (ot[p].size() != 0) {
			printf("c | Var(%d):\n", v + 1);
			printOL(ot[p]);
		}
	}
	printf("c | Negative occurs:\n");
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 n = NEG(V2D(v + 1));
		if (ot[n].size() != 0) {
			printf("c | Var(%d):\n", -int(v + 1));
			printOL(ot[n]);
		}
	}
}

/*****************************************************/
/*  Name:     ParaFROST                              */
/*  Usage:    global handler for solver/simplifier   */
/*  Scope:    host only                              */
/*  Memory:   system memory                          */
/*****************************************************/
class ParaFROST {
protected:
	string path;
	TIMER* timer;
	Vec<OCCUR> occurs; Vec<SCORE> scores;
	WT wt; BCNF orgs, bins; LCNF learnts;
	VAR_HEAP* var_heap; PV* pv; SP* sp; SOL* sol;
	G_REF* source; 
	int* board;	uint32* tmp_stack, *simpLearnt;
	uVector1D learnt_cl, learntLits; vector1D trail_sz;
	LEARN lrn; CL_PARAM cl_params; STATS stats;
	Byte* sysMem; size_t sysMem_sz;
	double sysMemTot, sysMemCons;
	int64 maxConflicts, nConflicts;
	int starts, restarts, R, ref_vars, reductions, marker;
	BQUEUE<uint32> lbdQ, trailQ; float lbdSum; 
	std::ofstream proofFile;
	bool intr, units_f;
public:
	ParaFROST(const string&);
	~ParaFROST(void);
	//============== inline methods ===============
	inline void interrupt(void) { intr = true; }
	inline bool interrupted(void) { return intr; }
	inline void write_proof(const Byte& byte) { proofFile << byte; }
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
	inline double drand(void) const { return ((double)rand() / (double)RAND_MAX); };
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
	inline void incDL(void) { trail_sz.push(sp->trail_size); }
	inline int DL(void) const { return trail_sz.size(); }
	inline void clHist(const G_REF gc) {
		B_REF c = (B_REF)gc;
		for (LIT_POS i = 0; i < c->size(); i++) {
			assert((*c)[i] > 0);
			if (ISNEG((*c)[i])) occurs[V2IDX((*c)[i])].ns++;
			else occurs[V2IDX((*c)[i])].ps++;
		}
	}
	inline void clBumpAct(C_REF c) {
		c->inc_act(cl_params.cl_inc);
		if (c->activity() > (float)1e20) {
			for (int i = 0; i < learnts.size(); i++) learnts[i]->sca_act((float)1e-20);
			cl_params.cl_inc *= (float)1e-20;
		}
	}
	inline void clDecayAct() { cl_params.cl_inc *= (1 / cl_params.cl_decay); }
	inline uint32 calcLBD(uVector1D& c) {
		marker++;
		register uint32 lbd = 0;
		for (LIT_POS i = 0; i < c.size(); i++) {
			int litLevel = sol->level(V2IDX(c[i]));
			if (board[litLevel] != marker) { board[litLevel] = marker; lbd++; }
		}
		return lbd;
	}
	inline uint32 calcLBD(C_REF c) {
		marker++;
		register uint32 lbd = 0;
		for (LIT_POS i = 0; i < c->size(); i++) {
			int litLevel = sol->level(V2IDX((*c)[i]));
			if (board[litLevel] != marker) { board[litLevel] = marker; lbd++; }
		}
		return lbd;
	}
	inline void init(bool re = true) {
		if (SH == 2) {
			maxConflicts = UNKNOWN;
			marker = UNKNOWN;
			reductions = 1;
			lrn.nClsReduce = nClsReduce;
			lbdQ.reset();
			trailQ.reset();
		}
		nConflicts = UNKNOWN;
		ref_vars = UNKNOWN;
		cnf_stats.n_added_lits = cnf_stats.global_n_gcs = UNKNOWN;
		lrn.init();
		stats.reset();
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
	inline void printWatched(const uint32& v_idx) {
		uint32 p = V2D(v_idx + 1);
		WL& pos_list = wt[NEG(p)], &neg_list = wt[p];
		printf("c |\t\tWL of v(%d)\n", v_idx + 1);
		printWL(pos_list); printWL(neg_list);
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
			printf("c | %9d %9d %8d %10lld %6d %9d %8d %10lld %7.0f |\n",
				remainedVars, nClauses(), nBins(), nLiterals(),
				starts,
				learnts.size(), nGlues(), nLearntLits(), learnts.size() == 0 ? 0 : nLearntLits() / (float)learnts.size());
			progRate += progRate / 4;
		}
	}
	//=============================================
	void free_mem(void);
	void CNF_rewriter(const string&);
	void var_order(void);
	void hist(BCNF&, bool rst = false);
	void hist(LCNF&, bool rst = false);
	void attachClause(B_REF);
	void attachClause(C_REF);
	void WT_alloc(bool re = false);
	void solver_alloc(bool re = false);
	void solver_init(void);
	int simplify(BCNF&);
	int simplify(LCNF&);
	void simplify(void);
	void simplify_top(void);
	void solve(void);
	CNF_STATE CNF_parser(const string&);
	CNF_STATE search(void);
	CNF_STATE decide(void);
	void pumpFrozen(void);
	void PDM_init(void);
	void PDM(void);
	void PDM_fuse(void);
	bool depFreeze_init(WL&, const uint32&);
	bool depFreeze(WL&, const uint32&);
	void eligibility(void);
	C_REF BCP(void);
	void enqueue(const uint32&, const int& pLevel = ROOT_LEVEL, const G_REF = NULL);
	void bt_level(void);
	void binSelfsub(void);
	void simp_learnt(void);
	bool selfsub(const uint32&, uint32*, CL_LEN&, const uint32&);
	void cbt_level(C_REF);
	void analyze(C_REF);
	void cancel_assigns(const int&);
	void backJump(const int&);
	void remWatch(WL&, const G_REF);
	void shrinkClause(G_REF);
	void detachClause(B_REF);
	void detachClause(C_REF);
	void lReduce(void);
	bool consistent(BCNF&, WT&);
	bool consistent(LCNF&, WT&);
	void print_reports(void);
	void print_model(void);
	void wrapUp(const CNF_STATE&);
	/* flags */
	bool quiet_en, parse_only_en, rewriter_en, perf_en;
	bool mcv_en, model_en, proof_en, fdp_en, cbt_en;
	int progRate, verbose, seed, timeout;
	int pdm_rounds, pdm_freq, pdm_order;
	int SH, polarity, restart_base, lbdRestBase, blockRestBase;
	int nClsReduce, lbdFrozen, lbdMinReduce, lbdMinClSize, incReduceSmall, incReduceBig;
	int cbt_dist, cbt_conf_max;
	double var_inc, var_decay;
	double RF, RB, restart_inc;
	string restPolicy, proof_path;
	/********************************************/
	/*                Simplifier                */
	/********************************************/
protected:
	uVector1D removed; vector1D mappedVars, reverseVars;
	SCNF scnf; OT ot;
	bool mapped;
public:
	/* inline helpers */
	inline uint32 mapLit(const uint32& orgLit) {
		static int maxVar = 1;
		assert(!mappedVars.empty());
		register int orgVar = ABS(orgLit);
		assert(!pv->melted[orgVar - 1]);
		assert(orgVar > 0 && orgVar < mappedVars.size());
		if (mappedVars[orgVar] == 0) { reverseVars[maxVar] = orgVar; mappedVars[orgVar] = maxVar++; }
		assert(maxVar - 1 <= nOrgVars());
		return (V2D(mappedVars[orgVar]) | ISNEG(orgLit));
	}
	inline uint32 revLit(const uint32& mappedLit) {
		return (V2D(reverseVars[ABS(mappedLit)]) | ISNEG(mappedLit));
	}
	inline void mapClause(SCLAUSE& s, BCLAUSE& d) {
		assert(d.size());
		for (LIT_POS l = 0; l < s.size(); l++) d[l] = mapLit(s[l]);
	}
	inline void clHist(const S_REF c) {
		assert(c->size());
		assert(c->status() == ORIGINAL || c->status() == LEARNT);
		for (LIT_POS i = 0; i < c->size(); i++) {
			assert(c->lit(i));
			if (ISNEG(c->lit(i))) occurs[V2IDX(c->lit(i))].ns++;
			else occurs[V2IDX(c->lit(i))].ps++;
		}
	}
	inline void cnt_cls()
	{
		cnf_stats.n_cls_after = 0;
		cnf_stats.n_lits_after = 0;
		for (int i = 0; i < scnf.size(); i++)
			if (scnf[i]->status() != DELETED) {
				cnf_stats.n_cls_after++;
				cnf_stats.n_lits_after += scnf[i]->size();
			}
	}
	inline void cnt_lits()
	{
		cnf_stats.n_lits_after = 0;
		for (int i = 0; i < scnf.size(); i++)
			if (scnf[i]->status() != DELETED)
				cnf_stats.n_lits_after += scnf[i]->size();
	}
	inline void eval_reds()
	{
		cnf_stats.n_del_vars = 0;
		cnf_stats.n_cls_after = 0;
		cnf_stats.n_lits_after = 0;
		for (int i = 0; i < scnf.size(); i++) {
			if (scnf[i]->status() != DELETED) {
				cnf_stats.n_cls_after++;
				cnf_stats.n_lits_after += scnf[i]->size();
				for (LIT_POS k = 0; k < scnf[i]->size(); k++)
					sp->frozen[V2IDX(scnf[i]->lit(k))] = true;
			}
		}
		bool* f_end = sp->frozen + nOrgVars();
		for (bool* f = sp->frozen; f != f_end; f++) {
			if (*f == 0) cnf_stats.n_del_vars++;
			else *f = 0;
		}
	}
	inline void logReductions()
	{
		printf("c |  Parallel variables = %d\n", pv->nPVs);
		printf("c |  Variables after = %d (-%d)\n", nOrgVars() - cnf_stats.n_del_vars, cnf_stats.n_del_vars);
		printf("c |  Clauses after = %d (-%d)\n", cnf_stats.n_cls_after, nClauses() - cnf_stats.n_cls_after);
		if (cnf_stats.n_lits_after > nLiterals()) printf("c |  Literals after = %lld (+%lld)\n", cnf_stats.n_lits_after, cnf_stats.n_lits_after - nLiterals());
		else printf("c |  Literals after = %lld (-%lld)\n", cnf_stats.n_lits_after, nLiterals() - cnf_stats.n_lits_after);
	}
	inline void print_pstats()
	{
		if (verbose == 1 && SH == 2) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			printf("c |p %8d %9d %8s %10lld %6d %9s %8s %10s %7s |\n",
				remainedVars, nClauses(), "---", nLiterals(), starts, "---", "---", "---", "---");
		}
		if (verbose == 1 && SH < 2) {
			int remainedVars = (int)nOrgVars() - (int)nRemVars();
			printf("c |p %8d %9d %10lld | %10lld %10s %10s %10s %7s |\n",
				remainedVars, nClauses(), nLiterals(), nConflicts, "---", "---", "---", "---");
		}
	}
	/******************/
	void opt_simp();
	void hist(const SCNF&, bool rst = false);
	void var_reorder(void);
	void extractBins(void);
	bool awaken(void);
	void cleanSlate(void);
	void create_ot(bool rst = true);
	CNF_STATE prop(void);
	void strengthen(S_REF, const uint32&);
	bool propClause(S_REF c, const uint32&);
	void attachClause(S_REF, const bool& added = true);
	void reattachClause(B_REF);
	void shrinkCNF(const bool countVars = false);
	void _simplify(void);
	void preprocess(void);
	void depFreeze(const OL&, const uint32&, const int&, const int&);
	void LCVE(void);
	void bve(void);
	void VE(void);
	void HSE(void);
	void HRE(void);
	void BCE(void);
	void SUB(void);
	bool consistent(const SCNF&, const OT&);
	/* flags */
	bool pre_en, lpre_en, simp_perf_en, ve_en, ve_plus_en, sub_en, bce_en, hre_en, all_en;
	int pre_delay, mu_pos, mu_neg, phases, cnf_free_freq;
};
extern ParaFROST* g_pFrost;

#endif 