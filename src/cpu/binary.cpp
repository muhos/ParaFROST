/***********************************************************************[binary.cpp]
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

inline uint32 Solver::analyzeReason(const C_REF& ref, const uint32& lit) {
  CHECKLIT(lit);
  assert(REASON(ref));
  int* levels = sp->level;
  const uint32 fit = FLIP(lit);
  CLAUSE& reason = cm[ref];
  PFLCLAUSE(4, reason, "   checking %d reason", l2i(fit));
  uint32 dom = 0;
  forall_clause(reason, k) {
    const uint32 other = *k, v = ABS(other);
    if (NEQUAL(other, fit) && levels[v]) {
      assert(isFalse(other));
      if (dom) return 0;
      dom = other;
    }
  }
  return dom;
}

uint32 Solver::hyper2Resolve(CLAUSE& c, const uint32& lit) {
  assert(DL() == 1);
  int* levels = sp->level;
  int nonRoots = 0;
  uint32 child = 0;
  forall_clause(c, k) {
    const uint32 other = *k, v = ABS(other);
    if (NEQUAL(other, lit) && levels[v]) {
      assert(isFalse(other));
      if (!nonRoots++) child = other;
    }
  }
  assert(nonRoots);
  if (nonRoots < 2) return 0;
  CHECKLIT(child);
  assert(analyzed.empty());
  PFLCLAUSE(4, c, "  Finding first dominator for %d via", l2i(child));
  stats.binary.resolutions++;
  uint32 dom = child, vom = ABS(dom), prev = 0;
  sp->seen[vom] = ANALYZED_M;
  analyzed.push(dom);
  C_REF r = sp->source[vom];
  while (REASON(r)) {
    prev = dom;
    if (!(dom = analyzeReason(r, prev))) break;
    assert(dom != prev);
    CHECKLIT(dom);
    assert(isFalse(dom));
    vom = ABS(dom);
    assert(!sp->seen[vom]);
    sp->seen[vom] = ANALYZED_M;
    analyzed.push(dom);
    r = sp->source[vom];
  }
  PFLOG2(4, "   found dominator %d of child %d", dom ? l2i(dom) : l2i(prev), l2i(child));
  const uint32 depth = analyzed.size();
  uint32 reset = 0;
  forall_clause(c, k) {
    const uint32 q = *k;
    CHECKLIT(q);
    if (q == lit || q == child || !levels[ABS(q)]) continue;
    assert(isFalse(q));
    PFLOG2(4, "  Finding next dominator for %d:", l2i(q));
    dom = q, vom = ABS(dom);
    r = sp->source[vom];
    while (!sp->seen[vom] && REASON(r)) {
      prev = dom;
      if (!(dom = analyzeReason(r, prev))) break;
      assert(dom != prev);
      CHECKLIT(dom);
      assert(isFalse(dom));
      vom = ABS(dom);
      r = sp->source[vom];
    }
    PFLOG2(4, "   found dominator %d of child %d", dom ? l2i(dom) : l2i(prev), l2i(q));
    while (reset < depth) {
      const uint32 a = analyzed[reset];
      CHECKLIT(a);
      if (a == dom) break;
      const uint32 av = ABS(a);
      assert(sp->seen[av]);
      sp->seen[av] = 0;
      reset++;
    }
    if (reset == depth) break;
  }
  dom = 0;
  while (reset < depth) {
    const uint32 a = analyzed[reset];
    CHECKLIT(a);
    if (!dom) dom = a;
    const uint32 av = ABS(a);
    assert(sp->seen[av]);
    sp->seen[av] = 0;
    reset++;
  }
  analyzed.clear();
  return dom;
}
