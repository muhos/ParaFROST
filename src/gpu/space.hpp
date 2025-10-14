/***********************************************************************[space.hpp]
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

#ifndef __SPACE_
#define __SPACE_

#include "malloc.hpp"
#include "definitions.hpp"
#include "vstate.hpp"

namespace ParaFROST {
	/*****************************************************/
	/*  Usage:    Information of search space            */
	/*  Dependency: none                                 */
	/*****************************************************/
	class SP {
		addr_t		_mem;
		size_t		_sz, _cap;
		template <class T>
		inline size_t calcBytes(const uint32& sz, const uint32& nVecs) const {
			assert(sz); 
			return sz * nVecs * sizeof(T);
		}

		#define forall_space(X) for (uint32 X = 1; X < _sz; X++)
		#define	breakline(X) if (X > 1 && X < _sz - 2 && X % 10 == 0) { \
						PUTCH('\n'); LOGN0("\t\t"); }
	public:
		// arrays
		int* level;
		uint32* tmpstack, *stacktail;
		uint64* board;
		C_REF* source;
		VSTATE* vstate;
		LIT_ST* seen, * frozen, * marks;
		LIT_ST* value, * psaved, * ptarget, * pbest;
		// scalers
		int learntLBD;
		int reasonsize, resolventsize;
		int conflictdepth, conflictsize;
		uint32 trailpivot;
		uint32 simplified;
		uint32 propagated;
		bool selfallocated;
		//================
				 SP	() { RESETSTRUCT(this); }
		explicit SP	(const size_t& size, const LIT_ST& pol) 
		{
			RESETSTRUCT(this);
			assert(sizeof(C_REF) == sizeof(uint64));
			assert(sizeof(VSTATE) == sizeof(Byte));
			const size_t vec8Bytes = calcBytes<C_REF>(size, 2);
			const size_t vec4Bytes = calcBytes<uint32>(size, 2);
			const size_t vec1Bytes = calcBytes<LIT_ST>(size, 9);
			_sz = size;
			_cap = vec1Bytes + vec4Bytes + vec8Bytes;
			assert(_cap);
			pfralloc(_mem, align_up(_cap, 64));
			assert(_mem != NULL);
			memset(_mem, 0, _cap);
			// 8-byte arrays
			source = (C_REF*)_mem;
			board = (uint64*)(source + _sz);
			// 4-byte arrays
			level = (int*)(_mem + vec8Bytes);
			tmpstack = (uint32*)(level + _sz);
			// 1-byte arrays
			value = (LIT_ST*)(_mem + vec8Bytes + vec4Bytes);
			frozen = value + _sz + _sz;
			seen = frozen + _sz;
			psaved = seen + _sz;
			ptarget = psaved + _sz;
			pbest = ptarget + _sz;
			marks = pbest + _sz;
			vstate = (VSTATE*)(marks + _sz);
			assert(_mem + _cap == addr_t(vstate) + _sz);
			// initialize with custom values
			memset(value, UNDEFINED, _sz + _sz);
			memset(marks, UNDEFINED, _sz);
			memset(ptarget, UNDEFINED, _sz);
			memset(pbest, UNDEFINED, _sz);
			memset(psaved, pol, _sz);
			forall_space(v) {
				level[v] = UNDEFINED;
				source[v] = NOREF;
			}
			selfallocated = true;
		}
		size_t	size		() const { return _sz; }
		size_t	capacity	() const { return _cap; }
		// Assume memory is already allocated and initialized properly.
		void	copyFrom	(SP* src)
		{
			assert(src->size() <= _sz);
			selfallocated = src->selfallocated;
			propagated = src->propagated;
			trailpivot = src->trailpivot;
			simplified = src->simplified;
			for (size_t v = 1; v < src->size(); v++) {
				const uint32 p = V2L(v), n = NEG(p);
				value[p] = src->value[p];
				value[n] = src->value[n];
				level[v] = src->level[v];
				source[v] = src->source[v];
				vstate[v] = src->vstate[v];
				psaved[v] = src->psaved[v];
			}
		}
		void 	printPhases	() {
			LOGN1(" Phases->[");
			forall_space(v) {
				PRINT("%5d:%d ", v, psaved[v]);
				breakline(v);
			}
			putc(']', stdout), PUTCH('\n');
		}
		void	printStates	() {
			LOGN1(" States->[");
			forall_space(v) {
				PRINT("%5d:%d ", v, vstate[v].state);
				breakline(v);
			}
			putc(']', stdout), PUTCH('\n');
		}
		void	printValues	() {
			LOGN1(" Values->[");
			forall_space(v) {
				uint32 lit = V2L(v);
				PRINT("%5d:%d ", l2i(lit), value[lit]);
				breakline(v);
			}
			putc(']', stdout), PUTCH('\n');
		}
		void	printLevels	() {
			LOGN1(" Levels->[");
			forall_space(v) {
				PRINT("%5d@%d ", v, level[v]);
				breakline(v);
			}
			putc(']', stdout), PUTCH('\n');
		}
		void	clearSubsume() { forall_space(v) vstate[v].subsume = 0; }
		void	destroy		() { if (_mem != NULL) std::free(_mem); _mem = NULL; }
				~SP			() { destroy(); }
	};
}

#endif