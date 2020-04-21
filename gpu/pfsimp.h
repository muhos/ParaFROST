#ifndef __SIGMA_SIMP_
#define __SIGMA_SIMP_

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/counting_iterator.h>
#include "pfcudefs.h"
#include "pfdefs.h"
#include "pfdvec.h"
typedef thrust::device_vector<uint32> uTHRVector;
extern cudaDeviceProp devProp;
extern int maxGPUThreads;

#define CE_POS_LMT 256
#define CE_NEG_LMT 256
#define MIN_PVARS 2
#define LIT_REM_THR 10
#define MAX_GL_RES_LEN 60
#define MAX_SH_RES_LEN 120
#define FAN_LMT 64
#define SHARED_CL_LEN 180
#define SH_MAX_HRE_IN 100
#define SH_MAX_HRE_OUT 250
#define CL_MAX_LEN_BCE 180 
#define BLUB 256
#define BLSIMP 64

class SCLAUSE {
	uint32* _lits;
	CL_LEN _sz;
	uint32 _sig;
	int8_t _st;
public:
	__host__ __device__ SCLAUSE() { _lits = NULL, _sz = 0, _st = 0, _sig = 0; }
	__host__ __device__ ~SCLAUSE() { clear(true); }
	__host__ __device__ uint32& operator [] (const LIT_POS& i) { assert(i < _sz); return _lits[i]; }
	__host__ __device__ uint32  operator [] (const LIT_POS& i) const { assert(i < _sz); return _lits[i]; }
	__host__ __device__ operator uint32* (void) { assert(_sz != 0); return _lits; }
	__host__ __device__ inline uint32* data(const uint32& idx) { return (_lits + idx); }
	__host__ __device__ inline CL_LEN size() const { return _sz; }
	__host__ __device__ inline void resize(const LIT_POS& size) { _sz = size; }
	__host__ __device__ inline void push(const uint32& lit) { _lits[_sz++] = lit; }
	__host__ __device__ inline void set_ptr(uint32* head) { _lits = head; }
	__host__ __device__ inline void set_sig(const uint32& sig) { _sig = sig; }
	__host__ __device__ inline void set_status(const int8_t& status) { assert(status <= ST_MASK); _st = (_st & ST_RST) | status; }
	__host__ __device__ inline void pop() { _sz--; }
	__host__ __device__ inline uint32 sig() const { return _sig; }
	__host__ __device__ inline int8_t status() const { return (_st & ST_MASK); }
	__host__ __device__ inline void markDeleted() { assert(DELETED == ST_MASK); _st |= DELETED; }
	__host__ __device__ inline bool molten() const { return (_st & DEL_MASK); } // is removable
	__host__ __device__ inline void melt(void) { _st |= DEL_MASK; } // set for deletion
	__host__ __device__ inline void freeze() { _st &= DEL_RST; } // prevent from deletion
	__host__ __device__ inline void copyFrom(SCLAUSE& src) {
		_sig = src.sig();
		_sz = src.size();
		uint32* d = _lits, * e = d + _sz, * s = src;
		while (d != e) *d++ = *s++;
	}
	__host__ __device__ inline void copyFrom(uint32* src, const CL_LEN& size) {
		_sz = size;
		uint32* d = _lits, * e = d + _sz, * s = src;
		while (d != e) *d++ = *s++;
	}
	__host__ __device__ inline void filter(void) {
		_sig = 0;
		CL_LEN newSz = 1;
		for (LIT_POS k = 1; k < _sz; k++) {
			uint32 next = _lits[k];
			if (_lits[k - 1] != next) {
				_lits[newSz++] = next;
				_sig |= MAPHASH(next);
			}
		}
		_sz = newSz;
	}
	__host__ __device__ inline void copyShared(uint32* src, const CL_LEN& size) {
		_sig = 0, _sz = 1;
		*_lits = *src;
		for (LIT_POS k = 1; k < size; k++) {
			uint32 next = src[k];
			if (src[k - 1] != next) {
				_lits[_sz++] = next;
				_sig |= MAPHASH(next);
			}
		}
	}
	__host__ __device__ inline void shareTo(uint32* dest) {
		assert(_sz > 1);
		uint32* d = dest, *s = _lits, *e = s + _sz;
		while (s != e) *d++ = *s++;
		assert(d - dest <= SHARED_CL_LEN && d - dest == _sz);
	}
	__host__ __device__ inline bool has(const uint32& lit) const { // binary search
		if (_sz == 2) {
			if (_lits[0] == lit || _lits[1] == lit) return true;
			else return false;
		}
		else {
			assert(this->isSorted());
			LIT_POS low = 0, high = _sz - 1, mid;
			uint32 first = _lits[low], last = _lits[high];
			while (first <= lit && last >= lit) {
				mid = (low + high) >> 1;
				uint32 m = _lits[mid];
				if (m < lit) first = _lits[low = mid + 1];
				else if (m > lit) last = _lits[high = mid - 1];
				else return true; // found
			}
			if (_lits[low] == lit) return true; // found
			else return false; // Not found
		}
	}
	__host__ __device__ inline void calcSig(const uint32& init_sig = 0) {
		_sig = init_sig;
		for (LIT_POS l = 0; l < _sz; l++) _sig |= MAPHASH(_lits[l]);
	}
	__host__ __device__ inline void clear(const bool& _free = false) {
		if (_free) _lits = NULL;
		_sz = 0;
	}
	__host__ __device__ inline bool isSorted() const {
		for (LIT_POS i = 0; i < _sz; i++) {
			if (i > 0 && _lits[i] < _lits[i - 1]) return false;
		}
		return true;
	}
	__host__ __device__ inline LIT_POS hasZero() const {
		for (LIT_POS i = 0; i < _sz; i++) {
			if (_lits[i] == 0) return i;
		}
		return -1;
	}
	__host__ __device__ inline void print() {
		printf("(");
		for (LIT_POS l = 0; l < _sz; l++) {
			int lit = int(ABS(_lits[l]));
			lit = (ISNEG(_lits[l])) ? -lit : lit;
			printf("%4d ", lit);
		}
		if (status() == ORIGINAL) printf(") O, s=%#0lx\n", _sig);
		else if (status() == LEARNT) printf(") A, s=%#0lx\n", _sig);
		else if (status() == DELETED) printf(") X, s=%#0lx\n", _sig);
		else printf(") U, s=%d\n", _sig);
	}
};
typedef SCLAUSE* S_REF;

class CNF {
	S_REF cls;
	uint32* lits;
	uint32 nCls_cap, n_cls;
	uint64 nLits_cap, n_lits;
public:
	__host__ __device__ CNF() { cls = NULL, lits = NULL, nCls_cap = 0, n_cls = 0, nLits_cap = 0, n_lits = 0; }
	__host__ __device__ ~CNF() { for (uint32 i = 0; i < n_cls; i++) cls[i].~SCLAUSE(); cls = NULL, lits = NULL; }
	__host__ __device__ inline void allocMem(const uint32& clsCap, const uint64& litsCap) {
		assert(clsCap > 0);
		assert(litsCap > 0);
		nCls_cap = clsCap, nLits_cap = litsCap;
		cls = (S_REF)(this + 1);
		lits = (uint32*)(cls + nCls_cap);
	}
	__host__ __device__ SCLAUSE& operator [] (const uint32& i) { return cls[i]; }
	__host__ __device__ SCLAUSE operator [] (const uint32& i) const { return cls[i]; }
	__host__ __device__ operator S_REF (void) { return cls; }
	__host__ __device__ void attach(uint32* clLits, const CL_LEN& clSize) {
		SCLAUSE* c = cls + n_cls;
		c->set_ptr(lits + n_lits);
		c->copyFrom(clLits, clSize);
		c->set_status(ORIGINAL);
		n_cls++, n_lits += clSize;
	}
	__device__ inline uint32 jumpCls(const uint32& offset);
	__device__ inline uint64 jumpLits(const uint64& offset);
	__device__ inline void resize(const uint32& nCls) { n_cls = nCls; }
	__device__ inline void resizeData(const uint64& nLits) { n_lits = nLits; }
	__host__ __device__ inline uint32 size() const { return n_cls; }
	__host__ __device__ inline uint32 capacity() const { return nCls_cap; }
	__host__ __device__ inline uint64 numLits() const { return n_lits; }
	__host__ __device__ inline uint64 litsCapacity() const { return nLits_cap; }
	__host__ __device__ inline uint32* data(const uint64& idx) { return (lits + idx); }
	__host__ __device__ inline void copyFrom(CNF* src) {
		uint32 sz = src->size(), nCls = 0;
		uint64 nLits = 0;
		for (uint32 i = 0; i < sz; i++) {
			SCLAUSE& s = (*src)[i];
			if (s.status() == LEARNT || s.status() == ORIGINAL) {
				S_REF d = cls + nCls;
				d->set_ptr(lits + nLits);
				d->set_status(ORIGINAL);
				d->copyFrom(s);
				nCls++, nLits += s.size();
			}
		}
		n_cls = nCls, n_lits = nLits;
	}
	__host__ __device__ inline void print_remained() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size() && cls[c].status() != DELETED) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	__host__ __device__ inline void print() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size()) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	__host__ __device__ inline void print_sig() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size()) printf("c | C(%d)->sig(%d)\n", c, cls[c].sig());
		}
	}
	__host__ __device__ inline void print_deleted() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].status() == DELETED) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	__host__ __device__ inline void dump() {
		printf("c | [");
		for (uint64 i = 0; i < n_lits; i++) {
			int lit = int(ABS(lits[i]));
			lit = (ISNEG(lits[i])) ? -lit : lit;
			printf("%6d", lit);
		}
		printf("]\n");
	}
};

class OL {
	uint32* ol;
	uint32 sz, cap;
public:
	__host__ __device__ OL() : ol(NULL), sz(0), cap(0) {}
	__host__ __device__ ~OL() { clear(true); };
	__host__ __device__ inline void allocList(uint32* head, const uint32& cap) { ol = head, this->cap = cap; }
	__device__ inline void push(const uint32& cl_idx);
	__host__ __device__ uint32& operator [] (const uint32& i) { assert(i < sz); return ol[i]; }
	__host__ __device__ uint32 operator [] (const uint32& i) const { assert(i < sz); return ol[i]; }
	__host__ __device__ operator uint32* (void) { return ol; }
	__host__ __device__ inline void resize(const uint32& sz) { this->sz = sz; }
	__host__ __device__ inline void shrink(const uint32& n) { assert(n <= sz); sz -= n; }
	__host__ __device__ inline bool empty() const { return sz == 0; }
	__host__ __device__ inline uint32 size() const { return sz; }
	__host__ __device__ inline uint32 capacity() const { return cap; }
	__host__ __device__ inline void clear(const bool& _free = false) {
		if (_free) ol = NULL, cap = 0;
		sz = 0;
	}
	__host__ __device__ inline void remove(const uint32& i) { sz--; if (sz != i) ol[i] = ol[sz]; }
	__host__ __device__ inline void print() {
		printf("[");
		for (uint32 i = 0; i < sz; i++) {
			printf(" %3d", ol[i]);
			if (i < sz - 1) printf(",");
		}
		printf("]\n");
	}
	__host__ __device__ inline void print(CNF& cnf) {
		for (uint32 i = 0; i < sz; i++)
			printf("c | "), cnf[ol[i]].print();
	}
	__device__ inline void copyFrom(OL* src) {
		sz = src->size();
		assert(sz <= cap);
		uint32* d = ol, * e = d + sz, * s = *src;
		while (d != e)
			*d++ = *s++;
	}
};

class OT {
	OL* lists;
	uint32* occurs;
	int64 maxLists, maxEntries;
public:
	__host__ __device__ OT() {
		lists = NULL, occurs = NULL;
		maxLists = 0, maxEntries = 0;
	}
	__host__ __device__ ~OT() {
		for (int64 i = 0; i < maxLists; i++) lists[i].~OL();
		lists = NULL, occurs = NULL;
		maxLists = 0, maxEntries = 0;
	}
	__host__ __device__ void allocMem(addr_t* _mem, const int64& maxLists, const int64& maxEntries) {
		assert(maxLists > 0);
		assert(maxEntries > 0);
		this->maxLists = maxLists;
		this->maxEntries = maxEntries;
		lists = (OL*)(*_mem);
		*_mem += maxLists * sizeof(OL);
		occurs = (uint32*)(*_mem);
		*_mem += maxEntries * sizeof(uint32);
	}
	__host__ __device__ void allocMem(const int64& nlists, const int64& nEntries) {
		assert(nlists > 0);
		assert(nEntries > 0);
		maxLists = nlists;
		maxEntries = nEntries;
		lists = (OL*)(this + 1);
		occurs = (uint32*)(lists + maxLists);
	}
	__host__ __device__ OL& operator [] (const int64& i) { assert(i < maxLists); return lists[i]; }
	__host__ __device__ OL operator [] (const int64& i) const { assert(i < maxLists); return lists[i]; }
	__host__ __device__ operator OL* (void) { return lists; }
	__host__ __device__ inline int64 size() const { return maxLists; }
	__host__ __device__ inline int64 capacity() const { return maxEntries; }
	__host__ __device__ inline uint32* data(const int64& i) { return occurs + i; }
	__host__ __device__ inline bool accViolation() {
		for (int64 v = 2; v < maxLists; v++) {
			if (lists[v].size() > lists[v].capacity()) {
				int64 sign_v = ISNEG(v) ? -int64(ABS(v)) : ABS(v);
				printf("\nc | Allocated list size exceeded for variable (%lld)\n", sign_v);
				printf("c | Allocated size = %d, sz = %d\n", lists[v].capacity(), lists[v].size());
				return false;
			}
		}
		return true;
	}
	__host__ __device__ inline void print() {
		for (int64 v = 2; v < size(); v++) {
			int64 sign_v = ISNEG(v) ? -int64(ABS(v)) : ABS(v);
			if (lists[v].size() != 0) {
				printf("c | Var(%lld)->", sign_v);
				lists[v].print();
			}
		}
	}
};

template<typename T>
class SharedMemory
{
public:
	__device__ inline operator T* () {
		extern __shared__ int __smem[];
		return (T*)__smem;
	}

	__device__ inline operator const T* () const {
		extern __shared__ int __smem[];
		return (T*)__smem;
	}
};

struct GSTATS {
	Byte* seen;
	uint64 numLits;	
	uint32 numDelVars, numClauses;
};
struct GSOL {
	cuVec<uint32>* assigns;
	LIT_ST* value;
	uint32 head;
};
struct PV {
	GSOL* sol;
	cuVec<uint32>* pVars;
	uint32 numPVs, mu_inc;
	PV() : sol(NULL), pVars(NULL), numPVs(0), mu_inc(0) {}
	~PV() { sol = NULL, pVars = NULL; }
};
//=================================================================================//
//                           GPU Wrappers Declaration                              //
//=================================================================================//
void mem_set(addr_t, const Byte&, const size_t&);
void mem_set(LIT_ST*, const LIT_ST&, const size_t&);
void copy(uint32*, CNF*, const int64&);
void copyIf(uint32*, CNF*, GSTATS*);
void calc_vscores(OCCUR*, SCORE*, uint32*);
void calc_added(CNF*, OT*, PV*, GSTATS*);
void calc_sig(CNF*, const uint32&, const uint32&);
void create_ot(CNF*, OT*, const bool&);
void ve(CNF*, OT*, PV*);
void hse(CNF*, OT*, PV*);
void bce(CNF*, OT*, PV*);
void hre(CNF*, OT*, PV*);
void evalReds(CNF*, GSTATS*);
void countCls(CNF*, GSTATS*);
void countLits(CNF*, GSTATS*);

#endif // !__SIGMA_SIMP_
