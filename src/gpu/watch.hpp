/***********************************************************************[watch.hpp]
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

#ifndef __WATCH_
#define __WATCH_

#include "clause.hpp"
#include "vector.hpp"
#include "definitions.hpp"

namespace ParaFROST {

	struct WATCH {
		C_REF	ref;
		uint32	imp;
		int		size;

		inline		WATCH	() : ref(NOREF), imp(0), size(0) { }
		inline		WATCH	(const C_REF& cref, const int& _size, const uint32& _imp) :
			ref(cref), imp(_imp), size(_size) {}
		inline bool binary	() const { return size == 2; }
	};

	struct DWATCH {
		C_REF	ref;
		uint32	lit, imp;
		int		size;

		inline		DWATCH	() : ref(NOREF), lit(0), imp(0), size(0) {}
		inline		DWATCH	(const uint32 _lit, const uint32 _imp, const C_REF _ref, const int& _size) :
			ref(_ref), lit(_lit), imp(_imp), size(_size) {}
	};

	#define PRIORALLBINS(CODE) (CODE & 1)

	#define PRIORLEARNTBINS(CODE) (CODE & 2)

	#define forall_watches(WS, PTR) \
		for (WATCH* PTR = WS, *WSEND = WS.end(); PTR != WSEND; PTR++)

	#define forall_dwatches(DWS, PTR) \
		for (DWATCH* PTR = DWS, *DWSEND = DWS.end(); PTR != DWSEND; PTR++)
}

#endif