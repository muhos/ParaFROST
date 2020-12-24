/***********************************************************************[pfwatch.h]
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

#include "pfclause.h"
#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {

	struct WATCH {
		C_REF	ref;
		uint32	imp;
		int		size;

		inline		WATCH	() { ref = NOREF, size = 0, imp = 0; }
		inline		WATCH	(const C_REF& cref, const int& sz, const uint32& lit) { ref = cref, size = sz, imp = lit; }
		inline bool binary	() const { return size == 2; }
	};
}

#endif