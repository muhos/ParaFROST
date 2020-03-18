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
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <random>
#include <sched.h> 
#include <csignal>
#include "Vec.h"
using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::fstream;
using std::ifstream;
//=======================================//
//      ParaFROST Parameters & Macros    //
//=======================================//
#define MBYTE 0x00100000
#define KBYTE 0x00000400
#define GBYTE 0x40000000
#define LIT_LEN 32
#define BUFFER_SIZE (MBYTE << 3)
#define MAX_CLAUSE_LEN 0xffffui16
#define NEG_SIGN 0x00000001
#define POS_SIGN 0x00000000
#define HASH_MASK 0x3F
#define TAUTOLOGY 0
#define CMP(a,b) (b == 0 ? a : b)
#define ABS(x) ((x) >> 1)
#define V2D(x) ((x) << 1)
#define V2IDX(x) (ABS(x) - 1)
#define ISNEG(x) (x & NEG_SIGN)
#define NEG(x) (x | NEG_SIGN)
#define FLIP(x) (x ^ NEG_SIGN)
#define HASH(x) (x & HASH_MASK)
#define UNKNOWN 0
#define ORIGINAL 1
#define LEARNT 2
#define GLUE 2
#define DELETED 3
#define ST_MASK (int8_t)0x03  // 0000-0011
#define IMP_MASK (int8_t)0x04 // 0000-0100
#define DEL_MASK (int8_t)0x08 // 0000-1000
#define BIN_MASK (int8_t)0x10 // 0001-0000
#define ST_RST (int8_t)0xFC   // xxxx-xx00
#define IMP_RST (int8_t)0xFB  // xxxx-x0xx
#define DEL_RST (int8_t)0xF7  // xxxx-0xxx
#define BIN_RST (int8_t)0xEF  // xxx0-xxxx
#define MAX(x,y) ((x > y) ? x : y)
#define MAP_LIT(x) (1ULL << HASH(x))
#define POS(x) (x & 0xFFFFFFFE)
#define ROOT_LEVEL 0
#define UNDEFINED -1
#define CUT_OFF -2
#define UNSOLVED -1
#define UNSAT 0
#define SAT 1
#define BL_RESTART_MIN 10000

// data types
typedef unsigned char Byte;
typedef unsigned int uint32;
typedef signed long long int int64;
typedef unsigned long long int uint64;
typedef int LIT_POS;
typedef LIT_POS CL_LEN;
typedef signed char CNF_STATE;
typedef signed char ASSIGN_ST;
typedef void* G_REF;
typedef Vec<uint32> uVector1D;
typedef Vec<uVector1D> uVector2D;
typedef Vec<int> vector1D;

void set_timeout(int);
void handler_exit(int);
void sig_handler(void handler(int));
/********************************************/
/* Platform-related directives */
#ifdef __linux__ 
#include <sys/resource.h>
#include <unistd.h>
#elif _WIN32
#include <windows.h>
#endif

//====================================================//
//       Global Data structures primitives            //
//====================================================//
/*****************************************************/
/*  Name:     CNF_INFO                               */
/*  Usage:    collecting information about BCNF       */
/*  Scope:    host only                              */
/*  Memory:   system memory                          */
/*  Dependency:  none                                */
/*****************************************************/
struct CNF_INFO {
	int max_org_cl_width;
	uint32 max_org_vars;
	uint32 max_org_cls;
	uint32 n_org_vars;
	uint32 n_org_cls;
	uint32 n_org_bins;
	uint32 n_del_vars;
	uint32 n_cls_after;
	uint32 global_n_del_vars;
	uint32 global_n_bins;
	uint32 global_n_cls;
	uint32 global_n_gcs;
	int64 global_n_lits;
	int64 n_org_lits;
	int64 n_added_lits;
	int64 n_lits_after;
};
extern CNF_INFO cnf_stats;
inline uint32 nOrgVars() {
	return cnf_stats.n_org_vars;
}
inline uint32 nOrgClauses() {
	return cnf_stats.n_org_cls;
}
inline uint32 nOrgBins() {
	return cnf_stats.n_org_bins;
}
inline int64 nOrgLits() {
	return cnf_stats.n_org_lits;
}
inline uint32 nClauses() {
	return cnf_stats.global_n_cls;
}
inline uint32 nBins() {
	return cnf_stats.global_n_bins;
}
inline uint32 nGlues() {
	return cnf_stats.global_n_gcs;
}
inline uint32 nRemVars() {
	return cnf_stats.global_n_del_vars;
}
inline int64 nLiterals() {
	return cnf_stats.global_n_lits;
}
inline int64 nLearntLits() {
	return cnf_stats.n_added_lits;
}

struct OCCUR {
	uint32 ps, ns;
	OCCUR() : ps(0), ns(0) {}
	void reset(void) { ps = 0; ns = 0; }
};

struct SCORE {
	uint32 v, sc;
	SCORE() : v(0), sc(0) {}
};

class TIMER {
private:
	clock_t _start, _stop;
	float cpuTime;

public:
	float par, vo, pdm;
	float bcp, bj, red;
	float ot, lcve, bve, hse, bce, hre, simp;

	TIMER() {
		_start = 0, _stop = 0, cpuTime = 0;
		par = 0, vo = 0, pdm = 0;
		bcp = 0, bj = 0, red = 0;
		lcve = 0, ot = 0, bve = 0, hse = 0, bce = 0, hre = 0, simp = 0;
	}

	~TIMER() {
		cpuTime = 0;
	}

	void start() {
		_start = clock();
	}

	void stop() {
		_stop = clock();
	}

	float CPU_time() {
		return cpuTime = ((float)abs(_stop - _start)) / CLOCKS_PER_SEC;
	}
};

#endif // __GL_DEFS_

