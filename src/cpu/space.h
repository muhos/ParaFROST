/***********************************************************************[space.h]
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

#include "malloc.h"
#include "definitions.h"

namespace pFROST {
	/*****************************************************/
	/*  Usage:    Information of search space            */
	/*  Dependency: none                                 */
	/*****************************************************/
	struct VSTATE {
		Byte state : 3, dlcount : 2, subsume : 1, probe : 1;
		VSTATE() : state(0), dlcount(0), subsume(0), probe(0) {}
	};
	class SP {
		addr_t		_mem;
		size_t		_sz, _cap;
		template <class T>
		inline size_t calcBytes(const uint32& sz, const uint32& nVecs) const {
			assert(sz); 
			return sz * nVecs * sizeof(T);
		}

		#define forall_space(x) for (uint32 x = 1; x < _sz; x++)
		#define	breakline(x) if (x > 1 && x < _sz - 2 && x % 10 == 0) { \
						putc('\n', stdout); PFLOGN0("\t\t"); }
	public:
		// arrays
		int* level;
		uint32* index;
		uint32* tmp_stack;
		uint64* board;
		C_REF* source;
		VSTATE* vstate;
		LIT_ST* seen, * frozen, * marks;
		LIT_ST* value, * psaved, * ptarget, * pbest;
		// scalers
		int learntLBD;
		uint32 stacktail;
		uint32 trailpivot;
		uint32 simplified;
		uint32 propagated;
		//================
		SP() { memset(this, 0, sizeof * this); }
		SP(const uint32& size) 
		{
			memset(this, 0, sizeof * this);
			assert(sizeof(C_REF) == sizeof(uint64));
			assert(sizeof(VSTATE) == sizeof(Byte));
			const size_t vec8Bytes = calcBytes<C_REF>(size, 2);
			const size_t vec4Bytes = calcBytes<uint32>(size, 3);
			const size_t vec1Bytes = calcBytes<LIT_ST>(size, 9);
			_sz = size;
			_cap = vec1Bytes + vec4Bytes + vec8Bytes;
			assert(_cap);
			pfralloc(_mem, _cap);
			assert(_mem != NULL);
			memset(_mem, 0, _cap);
			// 8-byte arrays
			source = (C_REF*)_mem;
			board = (uint64*)(source + _sz);
			// 4-byte arrays
			level = (int*)(_mem + vec8Bytes);
			index = (uint32*)(level + _sz);
			tmp_stack = index + _sz;
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
			forall_space(v) {
				level[v] = UNDEFINED;
				source[v] = NOREF;
			}
		}
		size_t	size		() const { return _sz; }
		size_t	capacity	() const { return _cap; }
		void	initSaved	(const LIT_ST& pol) {
			memset(psaved, pol, _sz);
		}
		void	copyFrom	(SP* src)
		{
			propagated = src->propagated;
			trailpivot = src->trailpivot;
			simplified = src->simplified;
			forall_space(v) {
				const uint32 p = V2L(v), n = NEG(p);
				value[p] = src->value[p];
				value[n] = src->value[n];
				level[v] = src->level[v];
				vstate[v] = src->vstate[v];
			}
		}
		void	printStates	() {
			PFLOGN1(" States->[");
			forall_space(v) {
				fprintf(stdout, "%5d:%d ", v, vstate[v].state);
				breakline(v);
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void	printValues	() {
			PFLOGN1(" Values->[");
			forall_space(v) {
				uint32 lit = V2L(v);
				fprintf(stdout, "%5d:%d ", l2i(lit), value[lit]);
				breakline(v);
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void	printLevels	() {
			PFLOGN1(" Levels->[");
			forall_space(v) {
				fprintf(stdout, "%5d@%d ", v, level[v]);
				breakline(v);
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void	clearSubsume() { forall_space(v) vstate[v].subsume = 0; }
		void	destroy		() { if (_mem != NULL) std::free(_mem); }
				~SP			() { destroy(); }
	};
}

#endif