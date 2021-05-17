/***********************************************************************[pfsolvertypes.h]
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

#ifndef __SOLVER_TYPES_
#define __SOLVER_TYPES_

#include "pfspace.h"
#include "pfwatch.h"
#include "pfsclause.h"

namespace pFROST {
	
	typedef Vec<C_REF> BCNF;
	typedef Vec<WATCH, int> WL;
	typedef Vec<WL> WT;
	typedef Vec<uint32, int> BOL;
	typedef Vec<C_REF, int> WOL;
	typedef Vec<S_REF, int> OL;
	typedef Vec<OL> OT;
	typedef Vec<S_REF, size_t> SCNF;
	
	struct CSIZE {
		C_REF ref;
		size_t size;
		CSIZE() {}
		CSIZE(const C_REF& _r, const int& _s) : ref(_r), size(_s) {}
	};

	struct DFS {
		uint32 idx, min;
		DFS() : idx(0), min(0) { }
	};

	#define forall_occurs(LIST, PTR) \
		for (S_REF* PTR = LIST, *END = LIST.end(); PTR != END; PTR++)
}

#endif