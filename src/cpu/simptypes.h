/***********************************************************************[simptypes.h]
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

#ifndef __SIMP_TYPES_
#define __SIMP_TYPES_

#include "datatypes.h"
#include "sclause.h"
#include "vector.h"

namespace ParaFROST {

typedef Vec<S_REF, int> OL;
typedef Vec<OL> OT;
typedef Vec<S_REF, size_t> SCNF;

#define forall_occurs(LIST, PTR) \
  for (S_REF* PTR = LIST, *END = LIST.end(); PTR != END; PTR++)
} // namespace ParaFROST

#endif