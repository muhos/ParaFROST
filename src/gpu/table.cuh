/***********************************************************************[table.cuh]
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

#ifndef __GPU_TABLE_
#define __GPU_TABLE_

#include "vector.hpp"
#include "vector.cuh"
#include "definitions.cuh"

namespace ParaFROST {

	typedef cuVec<S_REF> OL;
	typedef Vec<S_REF> HOL;
	typedef Vec<HOL> HOT;

	class OT {

		OL* _lists;
		S_REF* _occurs;
		uint32 maxLists;

	public:

		~OT() { _lists = NULL, _occurs = NULL; }
		_PFROST_H_D_			OT				() : _lists(NULL), _occurs(NULL), maxLists(0) {}
		_PFROST_H_D_			OT				(const uint32& nlists) : maxLists(nlists) {
			assert(nlists);
			_lists = (OL*)(this + 1);
			_occurs = (S_REF*)(_lists + maxLists);
		}
		_PFROST_H_D_	S_REF*	data			(const uint32& i = 0) { return _occurs + i; }
		_PFROST_H_D_	OL&		operator []		(const uint32& i) { assert(i < maxLists); return _lists[i]; }
		_PFROST_H_D_	OL		operator []		(const uint32& i) const { assert(i < maxLists); return _lists[i]; }
		_PFROST_H_D_			operator OL*	() { return _lists; }
		_PFROST_H_D_	uint32	size			() const { return maxLists; }
		_PFROST_H_D_	void	print			() {
			for (uint32 v = 2; v < size(); v++) {
				int sign_v = SIGN(v) ? -int(ABS(v)) : ABS(v);
				printf("c  list[ %d ][cap = %d]", sign_v, _lists[v].capacity()), _lists[v].print();
			}
		}
		_PFROST_H_		bool	accViolation	(const uint32& size) {
			for (uint32 v = 1; v <= size; v++) {
				uint32 p = V2L(v), n = NEG(p);
				if (_lists[p].size() > _lists[p].capacity()) {
					LOGERRORN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
						v, _lists[v].capacity(), _lists[v].size());
					return false;
				}
				if (_lists[n].size() > _lists[n].capacity()) {
					LOGERRORN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
						v, _lists[v].capacity(), _lists[v].size());
					return false;
				}
			}
			return true;
		}
	};

	#define forall_occurs(LIST, PTR) \
		for (S_REF* PTR = LIST, *LISTEND = LIST.end(); PTR != LISTEND; PTR++)

}

#endif
