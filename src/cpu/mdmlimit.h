/***********************************************************************[mdmlimit.hpp]
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

#ifndef __MDM_LIMIT_
#define __MDM_LIMIT_

namespace SeqFROST {

#define mdm_update                                                                         \
  {                                                                                        \
    last.mdm.decisions = trail.size() - sp->propagated;                                    \
    last.mdm.unassigned = inf.maxVar - last.mdm.decisions;                                 \
    last.mdm.rounds--;                                                                     \
    stats.decisions.multiple += last.mdm.decisions;                                        \
    PFLOG2(2, " MDM %d: %d decisions are elected (%.2f%%)",                                \
           stats.mdm.calls, last.mdm.decisions, percent(last.mdm.decisions, maxActive())); \
  }

} // namespace SeqFROST

#endif