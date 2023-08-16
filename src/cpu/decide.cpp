/***********************************************************************[decide.cpp]
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

#include "solve.h"
using namespace ParaFROST;

uint32 Solver::nextVSIDS() {
  assert(inf.unassigned);
  VSTATE* states = sp->vstate;
  uint32 cand = 0;
  while (!vsids.empty()) {
    cand = vsids.top();
    CHECKVAR(cand);
    if (!states[cand].state && UNASSIGNED(sp->value[V2L(cand)])) break;
    vsids.pop();
  }
  assert(cand);
  PFLOG2(4, " Next heap choice %d, activity %e", cand, activity[cand]);
  return cand;
}

uint32 Solver::nextVMFQ() {
  assert(inf.unassigned);
  VSTATE* states = sp->vstate;
  uint32 free = vmtf.free();
  assert(free);
  if (states[free].state || !UNASSIGNED(sp->value[V2L(free)])) {
    do { free = vmtf.previous(free); } while (states[free].state || !UNASSIGNED(sp->value[V2L(free)]));
    vmtf.update(free, bumps[free]);
  }
  PFLOG2(4, " Next queue choice %d, bumped %lld", free, bumps[free]);
  return free;
}

uint32 Solver::makeAssign(const uint32& v, const bool& tphase) {
  CHECKVAR(v);
  LIT_ST pol = UNDEFINED;
  if (tphase) pol = sp->ptarget[v];
  if (UNASSIGNED(pol)) pol = sp->psaved[v];
  assert(pol >= 0);
  return V2DEC(v, pol);
}

void Solver::decide() {
  assert(inf.unassigned);
  assert(sp->propagated == trail.size());
  assert(conflict == NOREF);
  assert(UNSOLVED(cnfstate));
  uint32 cand = vsidsEnabled() ? nextVSIDS() : nextVMFQ();
  uint32 dec = makeAssign(cand, useTarget());
  enqueueDecision(dec);
  stats.decisions.single++;
}