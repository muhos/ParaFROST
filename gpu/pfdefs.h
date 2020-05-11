#ifndef __GL_DEFS_
#define __GL_DEFS_
//=======================================//
//            C++ directives             //
//=======================================//
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <locale>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <climits>
#include <csignal>
#include "pfdtypes.h"
#include "pflogging.h"

using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::sort;
using std::fstream;
using std::ifstream;
// Platform-related directives
#ifdef __linux__ 
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#elif _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif

//=======================================//
//            Global Parameters          //
//=======================================//
#define MBYTE 0x00100000
#define KBYTE 0x00000400
#define GBYTE 0x40000000
#define INIT_CAP 32
#define DEFINED 2
#define UNSOLVED -1
#define UNDEFINED -1
#define TERMINATE -2
#define ROOT_LEVEL 0
#define UNKNOWN 0
#define UNSAT 0
#define SAT 1
#define GLUE 2
//======== DANGER ZONE =========
#define NEG_SIGN   0x00000001
#define HASH_MASK  0x0000001F
#define ABS(x)     ((x) >> 1)
#define V2D(x)     ((x) << 1)
#define V2X(x)     (ABS(x) - 1)
#define ISNEG(x)   (x & NEG_SIGN)
#define NEG(x)     (x | NEG_SIGN)
#define FLIP(x)    (x ^ NEG_SIGN)
#define HASH(x)    (x & HASH_MASK)
#define MAPHASH(x) (1UL << HASH(x))
#define POS(x)     (x & 0xFFFFFFFE)
#define MIN(x, y)  (x < y ? x : y)
#define ORIGINAL (int8_t)0x01
#define LEARNT   (int8_t)0x02
#define DELETED  (int8_t)0x03
#define ST_MASK  (int8_t)0x03  // 0000-0011
#define IMP_MASK (int8_t)0x04  // 0000-0100
#define DEL_MASK (int8_t)0x08  // 0000-1000
#define BIN_MASK (int8_t)0x10  // 0001-0000
#define GAR_MASK (int8_t)0x20  // 0010-0000
#define ST_RST   (int8_t)0xFC  // xxxx-xx00
#define IMP_RST  (int8_t)0xFB  // xxxx-x0xx
#define DEL_RST  (int8_t)0xF7  // xxxx-0xxx
#define BIN_RST  (int8_t)0xEF  // xxx0-xxxx
#define GAR_RST  (int8_t)0xDF  // xx0x-xxxx
//==============================
// interrupt handlers
void set_timeout(int);
void handler_terminate(int);
void handler_mercy_intr(int);
void handler_mercy_timeout(int);
void sig_handler(void h_intr(int), void h_timeout(int) = NULL);

struct OCCUR { uint32 ps, ns; };
struct SCORE { uint32 v, sc; };
struct CNF_INFO {
	uint32 n_org_vars, n_org_cls, n_org_bins, n_del_vars, n_cls_after, max_added_cls;
	uint32 global_n_del_vars, global_n_bins, global_n_gcs, global_n_cls;
	int64 n_org_lits, global_n_lits, n_added_lits, n_dual_vars, max_added_lits, n_lits_after;
	CNF_INFO() {
		max_added_cls = 0, max_added_lits = 0;
		n_del_vars = 0, n_cls_after = 0, n_lits_after = 0, n_added_lits = 0;
		n_org_vars = 0, n_org_cls = 0, n_org_bins = 0, n_org_lits = 0, n_dual_vars = 0;
		global_n_del_vars = 0, global_n_bins = 0, global_n_cls = 0, global_n_gcs = 0, global_n_lits = 0;
	}
};
extern CNF_INFO cnf_stats;

class TIMER {
private:
	clock_t _start, _stop;
	float _cpuTime;

public:
	float par, solve, pre;

	TIMER() {
		_start = 0, _stop = 0, _cpuTime = 0;
		par = 0, solve = 0, pre = 0;
	}
	~TIMER() { _cpuTime = 0; }

	void start() { _start = clock(); }
	void stop() { _stop = clock(); }
	float cpuTime() { return _cpuTime = ((float)abs(_stop - _start)) / CLOCKS_PER_SEC; }
};

inline uint32 nOrgVars() { return cnf_stats.n_org_vars; }
inline int64 nDualVars() { return cnf_stats.n_dual_vars; }
inline uint32 nOrgCls() { return cnf_stats.n_org_cls; }
inline uint32 nOrgBins() { return cnf_stats.n_org_bins; }
inline int64 nOrgLits() { return cnf_stats.n_org_lits; }
inline uint32 maxAddedCls() { return cnf_stats.max_added_cls; }
inline uint32 nClauses() { return cnf_stats.global_n_cls; }
inline uint32 nBins() { return cnf_stats.global_n_bins; }
inline uint32 nGlues() { return cnf_stats.global_n_gcs; }
inline uint32 nVarsDeleted() { return cnf_stats.global_n_del_vars; }
inline uint32 nVarsRemained() { return cnf_stats.n_org_vars - cnf_stats.global_n_del_vars; }
inline int64 maxAddedLits() { return cnf_stats.max_added_lits; }
inline int64 nLiterals() { return cnf_stats.global_n_lits; }
inline int64 nLearntLits() { return cnf_stats.n_added_lits; }
template<class T> inline bool eq(T& in, arg_t ref) {
	while (*ref) { if (*ref != *in) return false; ref++; in++; }
	return true;
}

#endif // !__GL_DEFS_

