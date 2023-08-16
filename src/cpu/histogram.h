/***********************************************************************[histogram.hpp]
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

#include "definitions.h"

namespace ParaFROST {

#define hist_score(V, OCCURS) (OCCURS + V)->ps*(OCCURS + V)->ns

#define hist_clause(C, HIST) \
  {                          \
    forall_clause(C, k) {    \
      HIST[*k]++;            \
    }                        \
  }

#define write_scores(VARS, SCORES, OCCURS) \
  uint32* SCORES = sp->tmpstack;           \
  OCCUR* OCCURS = occurs.data();           \
  uint32* VARS = eligible.data();          \
  {                                        \
    forall_variables(v) {                  \
      VARS[v - 1] = v;                     \
      SCORES[v] = hist_score(v, OCCURS);   \
    }                                      \
  }

#define count_occurs(C, OCCURS)    \
  {                                \
    forall_clause(C, k) {          \
      const uint32 lit = *k;       \
      CHECKLIT(lit);               \
      if (SIGN(lit))               \
        (OCCURS + ABS(lit))->ns++; \
      else                         \
        (OCCURS + ABS(lit))->ps++; \
    }                              \
  }

} // namespace ParaFROST