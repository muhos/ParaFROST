/***********************************************************************[pfspace.h]
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

#ifndef __SPACE_
#define __SPACE_

#include "pfralloc.h"
#include "pfdefs.h"

namespace pFROST {
	/*****************************************************/
	/*  Usage:    Information of search space            */
	/*  Dependency: none                                 */
	/*****************************************************/
	class SP {
		addr_t		_mem;
		size_t		_sz, _cap;
		template <class T>
		size_t		calcBytes(const uint32& sz, const uint32& nVecs) {
			assert(sz); return size_t(sz) * nVecs * sizeof(T);
		}
	public:
		int* level, * board, learnt_lbd;
		LIT_ST* locked, * seen, * frozen, * vstate, * marks, * subsume;
		LIT_ST* value, * psaved, * ptarget, * pbest;
		uint32* tmp_stack, propagated, trailpivot, simplified;
		C_REF* source;
		SP(const uint32& size) :
			_mem(NULL)
			, _sz(size)
			, _cap(0)
			, trailpivot(0)
			, propagated(0)
			, simplified(0)
			, learnt_lbd(0)
		{
			size_t vec1Bytes = calcBytes<LIT_ST>(size, 11);
			size_t vec4Bytes = calcBytes<uint32>(size, 3);
			size_t vec8Bytes = calcBytes<C_REF>(size, 1);
			_cap = vec1Bytes + vec4Bytes + vec8Bytes;
			assert(_cap);
			pfalloc(_mem, _cap);
			assert(_mem != NULL);
			memset(_mem, 0, _cap);
			// 1-byte arrays
			value = (LIT_ST*)_mem;
			locked = value + _sz + _sz;
			frozen = locked + _sz;
			seen = frozen + _sz;
			psaved = seen + _sz;
			ptarget = psaved + _sz;
			pbest = ptarget + _sz;
			vstate = pbest + _sz;
			marks = vstate + _sz;
			subsume = marks + _sz;
			// 4-byte arrays
			level = (int*)(_mem + vec1Bytes);
			board = level + _sz;
			tmp_stack = (uint32*)(board + _sz);
			// 8-byte arrays
			source = (C_REF*)(_mem + vec4Bytes + vec1Bytes);
			assert(_mem + _cap == addr_t(source) + vec8Bytes);
		}
		void		printStates(const uint32& size) {
			PFLOGN1(" States->[");
			for (uint32 v = 1; v <= size; v++) {
				fprintf(stdout, "%5d:%d ", v, vstate[v]);
				if (v > 1 && v < size - 1 && v % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		printValues(const uint32& size) {
			PFLOGN1(" Values->[");
			for (uint32 lit = 2; lit <= V2L(size); lit++) {
				fprintf(stdout, "%5d:%d ", l2i(lit), value[lit]);
				if (lit > 2 && lit < size - 1 && lit % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		printLevels(const uint32& size) {
			PFLOGN1(" Levels->[");
			for (uint32 v = 1; v <= size; v++) {
				fprintf(stdout, "%5d@%d ", v, level[v]);
				if (v > 1 && v < size - 1 && v % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		lockMelted(const uint32& size) {
			for (uint32 v = 1; v <= size; v++)
				if (vstate[v] == MELTED)
					locked[v] = 1;
		}
		void		clearBoard() { memset(board, 0, _sz); }
		void		clearSubsume() { memset(subsume, 0, _sz); }
		void		destroy() { if (_mem != NULL) std::free(_mem); }
		~SP() { destroy(); }
	};
}

#endif