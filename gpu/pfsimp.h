#ifndef __SIGMA_SIMP_
#define __SIGMA_SIMP_

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/counting_iterator.h>
#include "pfcudefs.h"
#include "pfdefs.h"
typedef thrust::device_vector<uint32> uTHRVector;
extern cudaDeviceProp devProp;
extern uint32 maxGPUThreads;

#define CE_POS_LMT 256
#define CE_NEG_LMT 256
#define MIN_PVARS 2
#define LIT_REM_THR 10
#define MAX_GL_RES_LEN 60
#define MAX_SH_RES_LEN 120
#define FAN_LMT 64
#define SH_MAX_HSE_IN 180
#define SH_MAX_BCE_IN 180 
#define SH_MAX_HRE_IN 100
#define SH_MAX_HRE_OUT 250
#define BLUB 256
#define BLSIMP 64

template<typename T>
class SharedMemory
{
public:
	_PFROST_D_ operator T* () {
		extern __shared__ int _smem[];
		return (T*)_smem;
	}
	_PFROST_D_ operator const T* () const {
		extern __shared__ int _smem[];
		return (T*)_smem;
	}
};

template<typename T>
class cuVec {
	T* _mem;
	uint32 sz, cap;
public:
	_PFROST_H_D_ cuVec() { _mem = NULL, sz = 0, cap = 0; }
	_PFROST_H_D_ ~cuVec() { clear(true); }
	_PFROST_D_ void shrink(const uint32&);
	_PFROST_D_ void insert(const T&);
	_PFROST_D_ void push(const T&);
	_PFROST_D_ void pop(void);
	_PFROST_D_ void copyFrom(cuVec<T>& src) {
		sz = src.size();
		assert(sz <= cap);
		T* d = _mem, * e = d + sz, * s = src;
		while (d != e) *d++ = *s++;
	}
	_PFROST_H_D_ void _pop(void) { sz--; }
	_PFROST_H_D_ void _shrink(const uint32& n) { sz -= n; }
	_PFROST_H_D_ void _push(const T& val) { _mem[sz++] = val; }
	_PFROST_H_D_ void alloc(T* head, const uint32& cap) { _mem = head, this->cap = cap; }
	_PFROST_H_D_ cuVec<T>& operator=(cuVec<T>& rhs) { return *this; }
	_PFROST_H_D_ const T& operator [] (const uint32& idx) const { assert(idx < sz); return _mem[idx]; }
	_PFROST_H_D_ T& operator [] (const uint32& idx) { assert(idx < sz); return _mem[idx]; }
	_PFROST_H_D_ operator T* (void) { return _mem; }
	_PFROST_H_D_ T* data(void) { return _mem; }
	_PFROST_H_D_ T& at(const uint32& idx) { assert(idx < sz); return _mem[idx]; }
	_PFROST_H_D_ bool empty(void) const { return sz == 0; }
	_PFROST_H_D_ uint32 size(void) const { return sz; }
	_PFROST_H_D_ uint32 capacity(void) const { return cap; }
	_PFROST_H_D_ void resize(const uint32& n) { assert(n <= cap); sz = n; }
	_PFROST_H_D_ void clear(const bool& _free = false) {
		if (_free) _mem = NULL, cap = 0;
		sz = 0;
	}
	_PFROST_H_D_ void print(const bool& litType = false) {
		printf("c | GPU Vector(size = %d)->[", sz);
		if (litType) {
			for (uint32 i = 0; i < sz; i++) {
				assert(_mem[i]);
				printf("%2d  ", (_mem[i] & 1) ? -int(ABS(_mem[i])) : int(ABS(_mem[i])));
				if (i && i < sz - 1 && i % 10 == 0) printf("\nc |\t\t\t");
			}
		}
		else {
			for (uint32 i = 0; i < sz; i++) {
				printf("%2d  ", _mem[i]);
				if (i && i < sz - 1 && i % 10 == 0) printf("\nc |\t\t\t");
			}
		}
		printf("]\n");
	}
};
typedef cuVec<uint32> cuVecU;
typedef cuVecU OL;

class SCLAUSE {
	uint32* _lits;
	CL_LEN _sz;
	uint32 _sig;
	int8_t _st;
public:
	_PFROST_H_D_ SCLAUSE() { _lits = NULL, _sz = 0, _st = 0, _sig = 0; }
	_PFROST_H_D_ ~SCLAUSE() { clear(true); }
	_PFROST_H_D_ uint32& operator [] (const LIT_POS& i) { assert(i < _sz); return _lits[i]; }
	_PFROST_H_D_ uint32  operator [] (const LIT_POS& i) const { assert(i < _sz); return _lits[i]; }
	_PFROST_H_D_ operator uint32* (void) { assert(_sz != 0); return _lits; }
	_PFROST_H_D_ uint32* data(const uint32& idx) { return (_lits + idx); }
	_PFROST_H_D_ CL_LEN size() const { return _sz; }
	_PFROST_H_D_ void resize(const LIT_POS& size) { _sz = size; }
	_PFROST_H_D_ void push(const uint32& lit) { _lits[_sz++] = lit; }
	_PFROST_H_D_ void set_ptr(uint32* head) { _lits = head; }
	_PFROST_H_D_ void set_sig(const uint32& sig) { _sig = sig; }
	_PFROST_H_D_ void set_status(const int8_t& status) { assert(status <= ST_MASK); _st = (_st & ST_RST) | status; }
	_PFROST_H_D_ void pop() { _sz--; }
	_PFROST_H_D_ uint32 sig() const { return _sig; }
	_PFROST_H_D_ int8_t status() const { return (_st & ST_MASK); }
	_PFROST_H_D_ void markDeleted() { assert(DELETED == ST_MASK); _st |= DELETED; }
	_PFROST_H_D_ bool molten() const { return (_st & DEL_MASK); } // is removable
	_PFROST_H_D_ void melt(void) { _st |= DEL_MASK; } // set for deletion
	_PFROST_H_D_ void freeze() { _st &= DEL_RST; } // prevent from deletion
	_PFROST_H_D_ void copyFrom(SCLAUSE& src) {
		_sig = src.sig();
		_sz = src.size();
		uint32* d = _lits, * e = d + _sz, * s = src;
		while (d != e) *d++ = *s++;
	}
	_PFROST_H_D_ void copyFrom(uint32* src, const CL_LEN& size) {
		_sz = size;
		uint32* d = _lits, * e = d + _sz, * s = src;
		while (d != e) *d++ = *s++;
	}
	_PFROST_H_D_ void filter(void) {
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
	_PFROST_H_D_ void copyShared(uint32* src, const CL_LEN& size) {
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
	_PFROST_H_D_ void shareTo(uint32* dest) {
		assert(_sz > 1);
		uint32* d = dest, *s = _lits, *e = s + _sz;
		while (s != e) *d++ = *s++;
		assert(d - dest == _sz);
	}
	_PFROST_H_D_ bool has(const uint32& lit) const { // binary search
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
	_PFROST_H_D_ void calcSig(const uint32& init_sig = 0) {
		_sig = init_sig;
		for (LIT_POS l = 0; l < _sz; l++) _sig |= MAPHASH(_lits[l]);
	}
	_PFROST_H_D_ void clear(const bool& _free = false) {
		if (_free) _lits = NULL;
		_sz = 0;
	}
	_PFROST_H_D_ bool isSorted() const {
		for (LIT_POS i = 0; i < _sz; i++) {
			if (i > 0 && _lits[i] < _lits[i - 1]) return false;
		}
		return true;
	}
	_PFROST_H_D_ LIT_POS hasZero() const {
		for (LIT_POS i = 0; i < _sz; i++) {
			if (_lits[i] == 0) return i;
		}
		return -1;
	}
	_PFROST_H_D_ void print() {
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
	CNF() { cls = NULL, lits = NULL, nCls_cap = 0, n_cls = 0, nLits_cap = 0, n_lits = 0; }
	~CNF();
	_PFROST_D_ void allocMem(const uint32& clsCap, const uint64& litsCap) {
		assert(clsCap > 0);
		assert(litsCap > 0);
		nCls_cap = clsCap, nLits_cap = litsCap;
		cls = (S_REF)(this + 1);
		lits = (uint32*)(cls + nCls_cap);
	}
	_PFROST_D_ uint32 jumpCls(const uint32& offset);
	_PFROST_D_ uint64 jumpLits(const uint64& offset);
	_PFROST_D_ void resize(const uint32& nCls) { n_cls = nCls; }
	_PFROST_D_ void resizeData(const uint64& nLits) { n_lits = nLits; }
	_PFROST_H_D_ SCLAUSE& operator [] (const uint32& i) { return cls[i]; }
	_PFROST_H_D_ SCLAUSE operator [] (const uint32& i) const { return cls[i]; }
	_PFROST_H_D_ operator S_REF (void) { return cls; }
	_PFROST_H_D_ void attach(uint32* clLits, const CL_LEN& clSize) {
		SCLAUSE* c = cls + n_cls;
		c->set_ptr(lits + n_lits);
		c->copyFrom(clLits, clSize);
		c->set_status(ORIGINAL);
		n_cls++, n_lits += clSize;
	}
	_PFROST_H_D_ uint32 size() const { return n_cls; }
	_PFROST_H_D_ uint32 capacity() const { return nCls_cap; }
	_PFROST_H_D_ uint64 numLits() const { return n_lits; }
	_PFROST_H_D_ uint64 litsCapacity() const { return nLits_cap; }
	_PFROST_H_D_ uint32* data(const uint64& idx) { return (lits + idx); }
	_PFROST_H_D_ void copyFrom(CNF* src) {
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
	_PFROST_H_D_ void print_remained() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size() && cls[c].status() != DELETED) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	_PFROST_H_D_ void print() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size()) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	_PFROST_H_D_ void print_sig() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].size()) printf("c | C(%d)->sig(%d)\n", c, cls[c].sig());
		}
	}
	_PFROST_H_D_ void print_deleted() {
		for (uint32 c = 0; c < size(); c++) {
			if (cls[c].status() == DELETED) {
				printf("c | C(%d)->", c);
				cls[c].print();
			}
		}
	}
	_PFROST_H_D_ void dump() {
		printf("c | [");
		for (uint64 i = 0; i < n_lits; i++) {
			int lit = int(ABS(lits[i]));
			lit = (ISNEG(lits[i])) ? -lit : lit;
			printf("%6d", lit);
		}
		printf("]\n");
	}
};

class OT {
	OL* lists;
	uint32* occurs;
	int64 maxLists, maxEntries;
public:
	OT() { lists = NULL, occurs = NULL, maxLists = 0, maxEntries = 0; }
	~OT();
	_PFROST_D_ void allocMem(const int64& nlists, const int64& nEntries) {
		assert(nlists > 0);
		assert(nEntries > 0);
		maxLists = nlists;
		maxEntries = nEntries;
		lists = (OL*)(this + 1);
		occurs = (uint32*)(lists + maxLists);
	}
	_PFROST_H_D_ OL& operator [] (const int64& i) { assert(i < maxLists); return lists[i]; }
	_PFROST_H_D_ OL operator [] (const int64& i) const { assert(i < maxLists); return lists[i]; }
	_PFROST_H_D_ operator OL* (void) { return lists; }
	_PFROST_H_D_ int64 size() const { return maxLists; }
	_PFROST_H_D_ int64 capacity() const { return maxEntries; }
	_PFROST_H_D_ uint32* data(const int64& i) { return occurs + i; }
	_PFROST_H_D_ bool accViolation() {
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
	_PFROST_H_D_ void print() {
		for (int64 v = 2; v < size(); v++) {
			int64 sign_v = ISNEG(v) ? -int64(ABS(v)) : ABS(v);
			if (lists[v].size() != 0) {
				printf("c | Var(%lld)->", sign_v);
				lists[v].print();
			}
		}
	}
};

struct GSTATS {
	Byte* seen;
	uint64 numLits;	
	uint32 numDelVars, numClauses;
	GSTATS() : seen(NULL), numDelVars(0), numClauses(0), numLits(0) {}
	~GSTATS() { seen = NULL; }
};
struct GSOL {
	cuVecU* assigns;
	LIT_ST* value;
	uint32 head;
	GSOL() : assigns(NULL), value(NULL), head(0) {}
	~GSOL() { assigns->~cuVec(), assigns = NULL, value = NULL; }
};
struct PV {
	GSOL* sol;
	cuVecU* pVars;
	uint32 numPVs, mu_inc;
	PV() : sol(NULL), pVars(NULL), numPVs(0), mu_inc(0) {}
	~PV() { sol->~GSOL(), pVars->~cuVec(), sol = NULL, pVars = NULL; }
};
//=================================================================================//
//                           GPU Wrappers Declaration                              //
//=================================================================================//
void clearClauses(CNF*);
void clearLists(OT*);
void mem_set(addr_t, const Byte&, const size_t&);
void mem_set(LIT_ST*, const LIT_ST&, const size_t&);
void copy(uint32*, CNF*, const int64&);
void copyIf(uint32*, CNF*, GSTATS*);
void calc_vscores(OCCUR*, SCORE*, uint32*);
void calc_added(CNF*, OT*, PV*, GSTATS*);
void calc_sig(CNF*, const uint32&, const uint32&);
void create_ot(CNF*, OT*, const bool&);
void reduce_ot(CNF*, OT*, const bool&);
CNF_STATE ve(CNF*, OT*, PV*);
CNF_STATE hse(CNF*, OT*, PV*);
void bce(CNF*, OT*, PV*);
void hre(CNF*, OT*, PV*);
void evalReds(CNF*, GSTATS*);
void countCls(CNF*, GSTATS*);
void countLits(CNF*, GSTATS*);

#endif 
