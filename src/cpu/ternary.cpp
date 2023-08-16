/***********************************************************************[ternary.cpp]
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

inline bool SCORS_CMP::operator()(const uint32& a, const uint32& b) const {
  const double as = solver->livescore(a), bs = solver->livescore(b);
  if (as < bs) return true;
  if (as > bs) return false;
  return a > b;
}

inline bool Solver::findBinary(uint32 first, uint32 second) {
  assert(wot.size());
  assert(active(first));
  assert(active(second));
  CHECKLIT(first);
  CHECKLIT(second);
  if (wot[first].size() > wot[second].size()) swap(first, second);
  CHECKLIT(first);
  WOL& list = wot[first];
  stats.ternary.checks += cacheLines(list.size(), sizeof(WATCH)) + 1;
  forall_wol(list, i) {
    const CLAUSE& c = cm[*i];
    if (!c.binary()) continue;
    const uint32 other = c[0] ^ c[1] ^ first;
    if (other == second) return true;
  }
  return false;
}

inline bool Solver::findTernary(uint32 first, uint32 second, uint32 third) {
  assert(wot.size());
  assert(active(first));
  assert(active(second));
  CHECKLIT(first);
  CHECKLIT(second);
  CHECKLIT(third);
  if (wot[second].size() > wot[third].size()) swap(second, third);
  if (wot[first].size() > wot[second].size()) swap(first, second);
  CHECKLIT(first);
  WOL& list = wot[first];
  stats.ternary.checks += cacheLines(list.size(), sizeof(WATCH)) + 1;
  forall_wol(list, i) {
    const CLAUSE& c = cm[*i];
    if (c.binary()) {
      const uint32 other = c[0] ^ c[1] ^ first;
      assert(other != first);
      CHECKLIT(other);
      if (other == second) return true;
      if (other == third) return true;
    } else {
      assert(c.size() == 3);
      stats.ternary.checks++;
      const uint32 x = c[0], y = c[1], z = c[2];
      if (x == first) {
        if (y == second && z == third) return true;
        if (z == second && y == third) return true;
      } else if (y == first) {
        if (x == second && z == third) return true;
        if (z == second && x == third) return true;
      } else if (z == first) {
        if (x == second && y == third) return true;
        if (y == second && x == third) return true;
      }
    }
  }
  return findBinary(second, third);
}

bool Solver::hyper3Resolve(CLAUSE& pos, CLAUSE& neg, const uint32& p) {
  CHECKLIT(p);
  assert(pos.size() == 3), assert(neg.size() == 3);
  assert(learntC.empty());
  PFLCLAUSE(4, pos, "  hyper ternary resolving %6d ", l2i(p));
  PFLCLAUSE(4, neg, "  with\t\t\t    ");
  stats.ternary.resolutions++;
  forall_clause(pos, k) {
    const uint32 lit = *k;
    assert(unassigned(lit));
    if (NEQUAL(lit, p)) learntC.push(lit);
  }
  assert(learntC.size() == 2);
  const uint32 n = NEG(p);
  const uint32 first = learntC[0], second = learntC[1];
  uint32 third = 0;
  forall_clause(neg, k) {
    const uint32 lit = *k;
    assert(unassigned(lit));
    if (NEQUAL(lit, n)) {
      const uint32 mask1 = NEQUAL(lit, first);
      if (mask1 == NEG_SIGN) return false;
      const uint32 mask2 = NEQUAL(lit, second);
      if (mask2 == NEG_SIGN) return false;
      if (mask1 && mask2) { // unique
        third = lit;
        learntC.push(lit);
      }
    }
  }
  const int size = learntC.size();
  if (size > 3) return false;
  if (size == 3 && findTernary(first, second, third)) return false;
  if (size == 2 && findBinary(first, second)) return false;
  return true;
}

void Solver::ternaryResolve(const uint32& p, const uint64& limit) {
  CHECKLIT(p);
  WOL &poss = wot[p], &negs = wot[NEG(p)];
  for (int i = 0; i < poss.size(); i++) {
    if (stats.ternary.checks > limit) break;
    const C_REF pref = poss[i];
    if (cm.deleted(pref)) continue;
    CLAUSE* pos = cm.clause(pref);
    if (pos->binary()) continue;
    stats.ternary.checks++;
    for (int j = 0; j < negs.size(); j++) {
      const C_REF nref = negs[j];
      assert(pref != nref);
      if (cm.deleted(nref)) continue;
      CLAUSE& neg = cm[nref];
      if (neg.binary()) continue;
      stats.ternary.checks++;
      if (hyper3Resolve(*pos, neg, p)) {
        if (opts.proof_en) proof.addClause(learntC);
        const int size = learntC.size();
        bool learnt = false;
        if (size == 3) {
          learnt = true;
          stats.ternary.ternaries++;
        } else {
          assert(size == 2);
          learnt = pos->learnt() && neg.learnt();
          PFLOG2(4, "  hyper ternary resolvent subsumes resolved clauses");
          removeClause(*pos, pref);
          removeClause(neg, nref);
          stats.ternary.binaries++;
        }
        newHyper3(learnt);
        pos = cm.clause(pref); // update if cm memory is reallocated
      }
      learntC.clear();
      assert(cm.clause(pref) == pos);
      if (pos->deleted()) break;
    }
  }
}

void Solver::scheduleTernary(LIT_ST* use) {
  assert(vschedule.empty());
  VSTATE* states = sp->vstate;
  forall_variables(v) {
    if (states[v].state) continue;
    const uint32 p = V2L(v), n = NEG(p);
    if (use[p] && use[n]) vschedule.insert(v);
  }
  PFLOG2(2, " Ternary %lld: scheduled %d variables %.2f%%",
         stats.ternary.calls, vschedule.size(), percent(vschedule.size(), maxActive()));
}

void Solver::attachTernary(BCNF& cnf, LIT_ST* use) {
  LIT_ST* value = sp->value;
  forall_cnf(cnf, i) {
    const C_REF ref = *i;
    if (cm.deleted(ref)) continue;
    CLAUSE& c = cm[ref];
    const int size = c.size();
    if (size <= 3 && UNASSIGNED(value[c[0]]) && UNASSIGNED(value[c[1]])) {
      if (size == 2)
        attachBinary(ref, c);
      else if (UNASSIGNED(value[c[2]])) {
        assert(size > 1);
        attachClause(ref, c);
        use[c[0]] = use[c[1]] = use[c[2]] = 1;
      }
    }
  }
}

void Solver::ternarying(const uint64& resolvents_limit, const uint64& checks_limit) {
  uint32 scheduled = vschedule.size();
  while (!vschedule.empty()) {
    if (stats.ternary.checks > checks_limit) break;
    if (last.ternary.resolvents > resolvents_limit) break;
    const uint32 cand = vschedule.pop();
    CHECKVAR(cand);
    assert(!sp->vstate[cand].state);
    ternaryResolve(V2L(cand), checks_limit);
  }
  uint32 processed = scheduled - vschedule.size();
  PFLOG2(2, " Ternary %lld: processed %d candidates %.2f%%",
         stats.ternary.calls, processed, percent(processed, vschedule.size()));
}

void Solver::ternary() {
  if (!cnfstate) return;
  if (interrupted()) return;
  if (!canTernary()) return;
  assert(!DL());
  assert(!last.ternary.resolvents);
  assert(sp->propagated == trail.size());
  SLEEPING(sleep.ternary, opts.ternary_sleep_en);
  stats.ternary.calls++;
  wt.clear(true);
  wot.resize(inf.nDualVars);
  LIT_ST* use = pfcalloc<LIT_ST>(inf.nDualVars);
  attachTernary(orgs, use);
  attachTernary(learnts, use);
  const int64 numClauses = maxClauses();
  const int64 resolvents_limit = numClauses * opts.ternary_perc;
  scheduleTernary(use);
  uint32 scheduled = vschedule.size();
  if (scheduled) {
    SET_BOUNDS(checks_limit, ternary, ternary.checks, searchticks, 2 * numClauses + nlogn(scheduled));
    ternarying(resolvents_limit, checks_limit);
  }
  free(use);
  wot.clear(true);
  vschedule.destroy();
  rebuildWT(opts.ternary_priorbins);
  if (retrail()) PFLOG2(2, " Propagation after ternary proved a contradiction");
  const int64 subsumed = numClauses + last.ternary.resolvents - maxClauses();
  PFLOG2(2, " Ternary %lld: added %lld resolvents %.2f%% and subsumed %lld clauses %.2f%%",
         stats.ternary.calls, last.ternary.resolvents, percent((double)last.ternary.resolvents, (double)numClauses),
         subsumed, percent(subsumed, numClauses));
  UPDATE_SLEEPER(ternary, last.ternary.resolvents);
  printStats(last.ternary.resolvents, 'h', CVIOLET1);
  last.ternary.resolvents = 0;
}

bool Solver::canTernary() {
  return opts.ternary_en && (uint64(maxClauses() << 1) < (stats.searchticks + opts.ternary_min_eff));
}