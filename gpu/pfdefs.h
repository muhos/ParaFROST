/***********************************************************************[pfdefs.h]
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
#include <random>
#include "pflogging.h"
#include "pfdtypes.h"
#include "pfconst.h"

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
using std::swap;
using std::string;
using std::ostream;
using std::fstream;
using std::ifstream;

namespace pFROST {

	//===================================================//
	//       Global Data structures primitives           //
	//===================================================//
	struct OCCUR { uint32 ps, ns; };
	struct CNF_INFO {
		uint32 orgVars, maxVar, maxFrozen, maxMelted, nDualVars, n_del_vars_after;
		uint32 nOrgCls, nOrgLits, n_cls_after, n_lits_after;
		uint32 nClauses, nLiterals, nLearntLits;
		CNF_INFO() { memset(this, 0, sizeof(*this)); }
	};
	extern CNF_INFO inf;

	class TIMER {
	private:
		clock_t _start, _stop;
		clock_t _start_p, _stop_p;
		float _cpuTime;
	public:
		float parse, solve, simp;
		float vo, ve, hse, bce, ere, cot, rot, sot, gc, io;
		TIMER			() { memset(this, 0, sizeof(*this)); }
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
	#define CHECKVAR(VAR) assert(VAR > 0 && VAR <= inf.maxVar)

	#define CHECKLIT(LIT) assert(LIT > 1 && LIT < inf.nDualVars)

	#define forall_variables(VAR) for (uint32 VAR = 1; VAR <= inf.maxVar; VAR++)

	#define forall_literals(LIT) for (uint32 LIT = 2; LIT < inf.nDualVars; LIT++)
	//====================================================//
	//                 Global Inline helpers              //
	//====================================================//
	template<class T>
	inline bool		eq				(T& in, arg_t ref) {
		while (*ref) { if (*ref != *in) return false; ref++; in++; }
		return true;
	}
	template<class T>
	inline bool		eqn(T in, arg_t ref, const bool& lower = false) {
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
	inline double	ratio			(const double& x, const double& y) { return y ? x / y : 0; }
	inline int		l2i				(const uint32& lit) { assert(lit > 1); return SIGN(lit) ? -int(ABS(lit)) : int(ABS(lit)); }
	inline uint32	maxInactive		() { return inf.maxMelted + inf.maxFrozen; }
	inline uint32	maxActive		() { assert(inf.maxVar >= maxInactive()); return inf.maxVar - maxInactive(); }
	inline uint32	maxLiterals		() { return inf.nLiterals + inf.nLearntLits; }
}

#endif // __GL_DEFS_

