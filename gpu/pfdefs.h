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
#include <cstdlib>
#include <fstream>
#include <climits>
#include <csignal>
#include "pfdtypes.h"
#include "pflogging.h"
#include "pfconst.h"
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

using std::swap;
using std::string;
using std::ostream;
using std::fstream;
using std::ifstream;

namespace pFROST {

	// interrupt handlers
	void set_timeout(int);
	void handler_terminate(int);
	void handler_mercy_intr(int);
	void handler_mercy_timeout(int);
	void sig_handler(void h_intr(int), void h_timeout(int) = NULL);
	//===================================================//
	//       Global Data structures primitives           //
	//===================================================//
	struct OCCUR { uint32 ps, ns; };
	struct CNF_INFO {
		uint32 maxVar, maxFrozen, maxMelted, nDualVars;
		uint32 maxAddedCls, maxAddedLits;
		uint32 nOrgCls, nOrgBins, nOrgLits, n_del_vars_after, n_cls_after, n_lits_after;
		uint32 nClauses, nGlues, nLiterals, nLearntBins, nLearntLits;
		CNF_INFO() {
			nOrgCls = 0, nOrgBins = 0, nOrgLits = 0;
			maxVar = 0, maxFrozen = 0, maxMelted = 0, nDualVars = 0;
			maxAddedCls = 0, maxAddedLits = 0;
			nLearntBins = 0, nClauses = 0, nGlues = 0, nLiterals = 0;
			n_del_vars_after = 0, n_cls_after = 0, n_lits_after = 0, nLearntLits = 0;
		}
	};
	extern CNF_INFO inf;

	class TIMER {
	private:
		clock_t _start, _stop;
		float _cpuTime;
	public:
		float parse, solve, simp;
		TIMER() {
			_start = 0, _stop = 0, _cpuTime = 0;
			parse = 0, solve = 0, simp = 0;
		}
		void start() { _start = clock(); }
		void stop() { _stop = clock(); }
		float cpuTime() { return _cpuTime = ((float)abs(_stop - _start)) / CLOCKS_PER_SEC; }
	};
	//====================================================//
	//                 Global Inline helpers              //
	//====================================================//
	template<class T>
	inline bool		eq				(T& in, arg_t ref) {
		while (*ref) { if (*ref != *in) return false; ref++; in++; }
		return true;
	}
	inline double	ratio			(const double& x, const double& y) { return y ? x / y : 0; }
	inline LIT_ST	flip			(const LIT_ST& sign) { return FLIP(sign); }
	inline LIT_ST	sign			(const uint32& lit) { assert(lit > 1); return LIT_ST(ISNEG(lit)); }
	inline uint32	flip			(const uint32& lit) { assert(lit > 1); return FLIP(lit); }
	inline uint32	neg				(const uint32& lit) { assert(lit > 1); return NEG(lit); }
	inline uint32	l2a				(const uint32& lit) { assert(lit > 1); return ABS(lit); }
	inline uint32	l2x				(const uint32& lit) { assert(lit > 1); return V2X(lit); }
	inline int		l2i				(const uint32& lit) { assert(lit > 1); return sign(lit) ? -int(l2a(lit)) : int(l2a(lit)); }
	inline uint32	v2l				(const uint32& v) { assert(v); return V2D(v); }
	inline uint32	v2dec			(const uint32& v, const LIT_ST phase) { assert(v); return (v2l(v) | phase); }
	inline void		printLit		(const uint32& lit) { fprintf(stdout, "%6d", l2i(lit)); }
	inline void		printVars		(const uint32* arr, const uint32& size, const LIT_ST& type = 'x') {
		fprintf(stdout, "(size = %d)->[", size);
		for (uint32 i = 0; i < size; i++) {
			if (type == 'l') printLit(arr[i]);
			else if (type == 'v') fprintf(stdout, "%5d", arr[i]);
			else fprintf(stdout, "%5d", arr[i] + 1);
			fprintf(stdout, "  ");
			if (i && i < size - 1 && i % 8 == 0) { putc('\n', stdout); PFLOGN0("\t\t\t"); }
		}
		putc(']', stdout), putc('\n', stdout);
	}
	inline uint32	maxActive		() { return inf.maxVar - (inf.maxMelted + inf.maxFrozen); }
	inline uint32	maxOrgLits		() { return inf.nLiterals + (inf.nOrgBins << 1); }
	inline uint32	maxLearntLits	() { return inf.nLearntLits + (inf.nLearntBins << 1); }
	inline uint32	maxLiterals		() { return maxOrgLits() + maxLearntLits(); }
}

#endif // !__GL_DEFS_

