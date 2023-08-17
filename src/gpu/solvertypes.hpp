/***********************************************************************[solvertypes.hpp]
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

#include "space.hpp"
#include "watch.hpp"

namespace ParaFROST {
	
	typedef Vec<C_REF> BCNF;
	typedef Vec<WATCH, int> WL;
	typedef Vec<WL> WT;
	typedef Vec<uint32, int> BOL;
	typedef Vec<C_REF, int> WOL;
	
	struct CSIZE {
		C_REF ref;
		uint32 size;
		CSIZE() {}
		CSIZE(const C_REF& _r, const uint32& _s) : ref(_r), size(_s) {}
	};

	struct DFS {
		uint32 idx, min;
		DFS() : idx(0), min(0) { }
	};

	#define forall_bol(BLIST, PTR) \
		for (uint32* PTR = BLIST, *END = BLIST.end(); PTR != END; PTR++)

	#define forall_wol(WLIST, PTR) \
		for (C_REF* PTR = WLIST, *END = WLIST.end(); PTR != END; PTR++)
}

#endif