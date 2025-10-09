/***********************************************************************[definitions.hpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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
#include <cstdint>
#include <climits>
#include <cstdlib>
#include <csignal>
#include "logging.hpp"
#include "datatypes.hpp"
#include "constants.hpp"

#if defined(__linux__)
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <cpuid.h>
#elif defined(__CYGWIN__)
#include </usr/include/sys/resource.h>
#include </usr/include/sys/mman.h>
#include </usr/include/sys/sysinfo.h>
#include </usr/include/sys/unistd.h>
#elif defined(_WIN32)
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
using std::string;
using std::ifstream;

namespace ParaFROST {

	struct OCCUR { uint32 ps, ns; };
	struct CNF_INFO {
		uint32 orgVars, orgCls, maxVar, maxDualVars; // relative to input CNF.
		uint32 maxFrozen, maxMelted, maxSubstituted, unassigned; // during solving.
		// relative to simplifications.
		uint32 currDeletedVars, prevDeletedVars; 
		uint32 numClausesSurvived, numClauses;
		uint32 numLiteralsSurvived, numLiterals;
		CNF_INFO() { RESETSTRUCT(this); }
	};
	extern CNF_INFO inf;

	class TIMER {
	private:
		clock_t _start, _stop;
		float _cpuTime;
	public:
		TIMER			() : _start(0), _stop(0), _cpuTime(0.0) {}
		void start		() { _start = clock(); }
		void stop		() { _stop = clock(); }
		float cpuTime	() { return _cpuTime = ((float)abs(_stop - _start)) / CLOCKS_PER_SEC; }
	};
	//====================================================//
	//                 iterators & checkers               //
	//====================================================//
	template <class T>
	inline bool  _checkvar(const T VAR) {
		const bool invariant = VAR <= 0 || VAR > inf.maxVar;
		if (invariant)
			LOGERRORN("invariant \"VAR > 0 && VAR <= inf.maxVar\" failed on variable (%lld), bound = %lld",
				int64(VAR), int64(inf.maxVar));
		return !invariant;
	}
	template <class T>
	inline bool  _checklit(const T LIT) {
		const bool invariant = LIT <= 1 || LIT >= inf.maxDualVars;
		if (invariant)
			LOGERRORN("invariant \"LIT > 1 && LIT < inf.maxDualVars\" failed on literal (%lld), bound = %lld",
				int64(LIT), int64(inf.maxDualVars));
		return !invariant;
	}
	#define CHECKVAR(VAR) assert(_checkvar(VAR))

	#define CHECKLIT(LIT) assert(_checklit(LIT))

	#define forall_variables(VAR) for (uint32 VAR = 1; VAR <= inf.maxVar; VAR++)

	#define forall_literal(LIT) for (uint32 LIT = 2; LIT < inf.maxDualVars; LIT++)

	#define forall_vector(TYPE, VEC, PTR) \
		for (TYPE* PTR = VEC, *VEND = VEC.end(); PTR != VEND; PTR++)

	#define forall_clause(C, PTR) \
		for (uint32* PTR = C, *CEND = C.end(); PTR != CEND; PTR++)

	#define forall_clause_const(C, PTR) \
		for (const uint32* PTR = C.data(), *CEND = C.end(); PTR != CEND; PTR++)

	#define forall_cnf(INCNF, PTR) \
		for (C_REF* PTR = INCNF, *CNFEND = INCNF.end(); PTR != CNFEND; PTR++)

	#define forall_occurs(LIST, PTR) \
		for (S_REF* PTR = LIST, *LISTEND = LIST.end(); PTR != LISTEND; PTR++)

	#define forall_occurs_const(LIST, PTR) \
		for (const S_REF* PTR = LIST.data(), *LISTEND = LIST.end(); PTR != LISTEND; PTR++)
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
	inline bool		hasstr			(const char* in, const char* ref)
	{
		size_t count = 0;
		const size_t reflen = strlen(ref);
		while (*in) {
			if (ref[count] != *in)
				count = 0;
			else
				count++;
			in++;
			if (count == reflen)
				return true;
		}
		return false;
	}
	inline double	ratio			(const double& x, const double& y) { return y ? x / y : 0; }
	inline uint64	ratio			(const uint64& x, const uint64& y) { return y ? x / y : 0; }
	inline double	percent			(const double& x, const double& y) { return ratio(100 * x, y); }
	inline int		l2i				(const uint32& lit) { return SIGN(lit) ? -int(ABS(lit)) : int(ABS(lit)); }
	inline uint32	maxInactive		() { return inf.maxMelted + inf.maxFrozen + inf.maxSubstituted; }
	inline uint32	maxActive		() { assert(inf.maxVar >= maxInactive()); return inf.maxVar - maxInactive(); }
	
}

#endif // __GL_DEFS_

