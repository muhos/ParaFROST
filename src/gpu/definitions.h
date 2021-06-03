/***********************************************************************[definitions.h]
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

#ifndef __GL_DEFS_
#define __GL_DEFS_
//=======================================//
//            C++ directives             //
//=======================================//
#include <iostream>
#include <algorithm>
#include <cstring>
#include <locale>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <csignal>
#include "logging.h"
#include "datatypes.h"
#include "const.h"

#ifdef __linux__ 
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <cpuid.h>
#elif _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#include <intrin.h>
#include <Winnt.h>
#include <io.h>
#endif
#undef ERROR
#undef hyper 
#undef SET_BOUNDS
using std::swap;
using std::string;
using std::ifstream;

namespace pFROST {

    #define RESETSTRUCT(MEMPTR) \
	{ \
		assert(MEMPTR); \
		memset(MEMPTR, 0, sizeof(*MEMPTR)); \
	}

	struct OCCUR { uint32 ps, ns; };
	struct CNF_INFO {
		uint32 orgVars, maxVar, maxFrozen, maxMelted, maxSubstituted, nDualVars, unassigned, n_del_vars_after;
		uint32 nOrgCls, n_cls_after, n_lits_after, nClauses, nLiterals;
		CNF_INFO() { RESETSTRUCT(this); }
	};
	extern CNF_INFO inf;

	class TIMER {
	private:
		clock_t _start, _stop;
		clock_t _start_p, _stop_p;
		float _cpuTime;
	public:
		float parse, solve, simp;
		float vo, ve, sub, bce, ere, cot, rot, sot, gc, io;
		TIMER			() { RESETSTRUCT(this); }
		void start		() { _start = clock(); }
		void stop		() { _stop = clock(); }
		float cpuTime	() { return _cpuTime = ((float)abs(_stop - _start)) / CLOCKS_PER_SEC; }
		void pstart		() { _start_p = clock(); }
		void pstop		() { _stop_p = clock(); }
		float pcpuTime	() { return _cpuTime = (((float)abs(_stop_p - _start_p)) / CLOCKS_PER_SEC) * float(1000.0); }
	};
	//====================================================//
	//                 iterators & checkers               //
	//====================================================//
	template <class T>
	inline bool  _checkvar(const T VAR) {
		const bool invariant = VAR <= 0 || VAR > inf.maxVar;
		if (invariant)
			PFLOGEN("invariant \"VAR > 0 && VAR <= inf.maxVar\" failed on variable (%lld), bound = %lld",
				int64(VAR), int64(inf.maxVar));
		return !invariant;
	}
	template <class T>
	inline bool  _checklit(const T LIT) {
		const bool invariant = LIT <= 1 || LIT >= inf.nDualVars;
		if (invariant)
			PFLOGEN("invariant \"LIT > 1 && LIT < inf.nDualVars\" failed on literal (%lld), bound = %lld",
				int64(LIT), int64(inf.nDualVars));
		return !invariant;
	}
	#define CHECKVAR(VAR) assert(_checkvar(VAR))

	#define CHECKLIT(LIT) assert(_checklit(LIT))

	#define forall_variables(VAR) for (uint32 VAR = 1; VAR <= inf.maxVar; VAR++)

	#define forall_literal(LIT) for (uint32 LIT = 2; LIT < inf.nDualVars; LIT++)

	#define forall_vector(TYPE, VEC, PTR) \
		for (TYPE* PTR = VEC, *VEND = VEC.end(); PTR != VEND; PTR++)

	#define forall_clause(C, PTR) \
		for (uint32* PTR = C, *CEND = C.end(); PTR != CEND; PTR++)

	#define forall_cnf(INCNF, PTR) \
		for (C_REF* PTR = INCNF, *CNFEND = INCNF.end(); PTR != CNFEND; PTR++)

	#define forall_occurs(LIST, PTR) \
		for (S_REF* PTR = LIST, *LISTEND = LIST.end(); PTR != LISTEND; PTR++)
	//====================================================//
	//                 Global Inline helpers              //
	//====================================================//
	template<class T>
	inline bool		eq				(T& in, arg_t ref) {
		while (*ref) { if (*ref != *in) return false; ref++; in++; }
		return true;
	}
	template<class T>
	inline bool		eqn				(T in, arg_t ref, const bool& lower = false) {
		if (lower) {
			while (*ref) { 
				if (tolower(*ref) != tolower(*in))
					return false; 
				ref++; in++;
			}
		}
		else {
			while (*ref) { if (*ref != *in) return false; ref++; in++; }
		}
		return true;
	}
	inline uint64	linear			(uint64 n) { return n; }
	inline uint64	quadratic		(uint64 n) { return n * n; }
	inline uint64	logn			(uint64 n) { return (uint64)log10(n + 10); }
	inline uint64	nlogn			(uint64 n) { return n * logn(n); }
	inline uint64	nbylogn			(uint64 n) { return n / logn(n); }
	inline uint64	lognlogn		(uint64 n) {
		double val = log10(n + 10);
		return uint64(val * val);
	}
	inline uint64	nlognlogn		(uint64 n) { return n * lognlogn(n); }
	inline double	ratio			(const double& x, const double& y) { return y ? x / y : 0; }
	inline uint64	ratio			(const uint64& x, const uint64& y) { return y ? x / y : 0; }
	inline double	percent			(const double& x, const double& y) { return ratio(100 * x, y); }
	inline int		l2i				(const uint32& lit) { CHECKLIT(lit); return SIGN(lit) ? -int(ABS(lit)) : int(ABS(lit)); }
	inline uint32	maxInactive		() { return inf.maxMelted + inf.maxFrozen + inf.maxSubstituted; }
	inline uint32	maxActive		() { assert(inf.maxVar >= maxInactive()); return inf.maxVar - maxInactive(); }
	
}

#endif // __GL_DEFS_

