/***********************************************************************[mdmassign.hpp]
Copyright(c) 2022, Muhammad Osama - Anton Wijs,
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

#ifndef __MDM_ASSIGN_
#define __MDM_ASSIGN_

#include "solvertypes.h"

namespace ParaFROST {

#define mdm_prefetch(VALUES, STATES, FROZEN, TAIL) \
  const LIT_ST* VALUES = sp->value;                \
  const VSTATE* STATES = sp->vstate;               \
  LIT_ST* FROZEN = sp->frozen;                     \
  sp->stacktail = sp->tmpstack;                    \
  uint32*& TAIL = sp->stacktail;

#define mdm_assign(CAND, DEC)                                           \
  assert(CAND == ABS(DEC));                                             \
  WL& ws = wt[DEC];                                                     \
  if (valid(values, ws) && depFreeze(CAND, values, frozen, tail, ws)) { \
    enqueueDecision(DEC);                                               \
    sp->seen[CAND] = 1;                                                 \
  }

} // namespace ParaFROST

#endif