/***********************************************************************[redundancy.h]
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

#ifndef __ERE_
#define __ERE_

#include "simplify.h"
using namespace ParaFROST;

inline void find_best(const uint32 x, const SCLAUSE& c, const OT& ot, int& minsize, uint32& sig, uint32& best) {
  assert(c.size() > 1);
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i];
    if (NEQUAL(ABS(lit), x)) {
      int lsize = ot[lit].size();
      if (lsize < minsize) minsize = lsize, best = lit;
      sig |= MAPHASH(lit);
    }
  }
}

inline bool merge_ere(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, const SCLAUSE& target) {
  CHECKVAR(x);
  assert(!c1.deleted());
  assert(!c2.deleted());
  assert(!target.deleted());
  const int n1 = c1.size(), n2 = c2.size(), n = target.size();
  int it1 = 0, it2 = 0, it = 0;
  uint32 lit1, lit2, v1, v2;
  while (it < n && it1 < n1 && it2 < n2) {
    const uint32 lit = target[it];
    if (ABS(lit) == x) // target is one of the antecedents
      return false;
    lit1 = c1[it1], lit2 = c2[it2];
    v1 = ABS(lit1), v2 = ABS(lit2);
    if (v1 == x) {
      it1++;
    } else if (v2 == x) {
      it2++;
    } else if (v1 < v2) {
      if (NEQUAL(lit1, lit))
        return false;
      it1++, it++;
    } else if (v2 < v1) {
      if (NEQUAL(lit2, lit))
        return false;
      it2++, it++;
    } else { // repeated literal
      assert(lit1 == lit2);
      if (NEQUAL(lit1, lit))
        return false;
      it1++, it2++, it++;
    }
  }
  while (it < n && it1 < n1) {
    lit1 = c1[it1];
    if (NEQUAL(ABS(lit1), x) && NEQUAL(lit1, target[it++]))
      return false;
    it1++;
  }
  while (it < n && it2 < n2) {
    lit2 = c2[it2];
    if (NEQUAL(ABS(lit2), x) && NEQUAL(lit2, target[it++]))
      return false;
    it2++;
  }
  return true;
}

inline void forward_equ(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, const OT& ot, const int& maxsize, ERESTATS& erestats) {
  CHECKVAR(x);
#ifdef STATISTICS
  erestats.tried++;
#endif
  const int n1 = c1.size(), n2 = c2.size();
  int len = n1 + n2 - 2;
  int it1 = 0, it2 = 0;
  int lsize, minsize = INT_MAX;
  uint32 msig = 0, best = 0;
  while (it1 < n1 && it2 < n2) {
    const uint32 lit1 = c1[it1], lit2 = c2[it2];
    const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
    if (v1 == x)
      it1++;
    else if (v2 == x)
      it2++;
    else if ((lit1 ^ lit2) == NEG_SIGN)
      return;
    else if (v1 < v2) {
      it1++, lsize = ot[lit1].size();
      if (lsize < minsize) minsize = lsize, best = lit1;
      msig |= MAPHASH(lit1);
    } else if (v2 < v1) {
      it2++, lsize = ot[lit2].size();
      if (lsize < minsize) minsize = lsize, best = lit2;
      msig |= MAPHASH(lit2);
    } else { // repeated literal
      assert(len > 0);
      assert(lit1 == lit2);
      it1++, it2++, len--, lsize = ot[lit1].size();
      if (lsize < minsize) minsize = lsize, best = lit1;
      msig |= MAPHASH(lit1);
    }
  }
  if (len > 1 && (!maxsize || len <= maxsize)) {
    const CL_ST type = (c1.learnt() || c2.learnt()) ? LEARNT : ORIGINAL;
    CHECKLIT(best);
    const OL& minlist = ot[best];
    assert(minsize == minlist.size());
    for (int i = 0; i < minsize; i++) {
      S_REF c = minlist[i];
      if (len == c->size()
          && (c->learnt() || (c->status() == type))
          && sub(msig, c->sig())
          && merge_ere(x, c1, c2, *c)) {
#ifdef STATISTICS
        if (c->learnt())
          erestats.learnts++;
        else
          erestats.orgs++;
#endif
        solver->removeClause(c);
        break;
      }
    }
  }
}

#endif
