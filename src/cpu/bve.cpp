/***********************************************************************[bve.cpp]
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

#include "and.h"
#include "equivalence.h"
#include "ifthenelse.h"
#include "redundancy.h"
#include "subsume.h"
#include "xor.h"

using namespace ParaFROST;

void Solver::bve() {
  assert(hc_isize == sizeof(uint32));
  assert(hc_scsize == sizeof(SCLAUSE));
  if (interrupted()) killSolver();
  if (opts.profile_simp) timer.pstart();
  Lits_t out_c;
  out_c.reserve(opts.ve_clause_limit);
#ifdef STATISTICS
  BVESTATS& bvestats = stats.sigma.bve;
#endif
  for (uint32 i = 0; i < PVs.size(); i++) {
    uint32 v = PVs[i];
    assert(v);
    assert(!sp->vstate[v].state);
    const uint32 p = V2L(v), n = NEG(p);
    OL &poss = ot[p], &negs = ot[n];
    int pOrgs = 0, nOrgs = 0;
    countOrgs(poss, pOrgs), countOrgs(negs, nOrgs);
    // pure-literal
    if (!pOrgs || !nOrgs) {
      toblivion(p, pOrgs, nOrgs, poss, negs, model);
#ifdef STATISTICS
      bvestats.pures++;
#endif
      v = 0;
    }
    // Equiv/NOT-gate Reasoning
    else if (uint32 def = find_BN_gate(p, poss, negs)) {
#ifdef STATISTICS
      bvestats.inverters++;
#endif
      save_BN_gate(p, pOrgs, nOrgs, poss, negs, model);
      if (substitute_single(p, def, ot)) {
        PFLOG2(2, "  BVE proved a contradiction");
        learnEmpty();
        killSolver();
      }
      v = 0;
    } else {
      assert(pOrgs && nOrgs);
      out_c.clear();
      const int nOrgCls = pOrgs + nOrgs;
      int nAddedCls = 0;
      Byte type = 0;
      if (nOrgCls > 2) {
        // AND/OR-gate Reasoning
        if (find_AO_gate(n, nOrgCls, ot, out_c, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.andors++;
#endif
        } else if (!nAddedCls && find_AO_gate(p, nOrgCls, ot, out_c, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.andors++;
#endif
        }
      }
      if (!type && nOrgCls > 3) {
        // ITE-gate Reasoning
        if (find_ITE_gate(p, nOrgCls, ot, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.ites++;
#endif
        } else if (!nAddedCls && find_ITE_gate(n, nOrgCls, ot, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.ites++;
#endif
        }
        // XOR-gate Reasoning
        else if (find_XOR_gate(p, nOrgCls, ot, out_c, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.xors++;
#endif
        } else if (!nAddedCls && find_XOR_gate(n, nOrgCls, ot, out_c, nAddedCls)) {
          type = SUBSTITUTION;
#ifdef STATISTICS
          bvestats.xors++;
#endif
        }
      }
      // n-by-m resolution
      if (!type && !nAddedCls && !countResolvents(v, nOrgCls, poss, negs, nAddedCls)) {
        type = RESOLUTION;
#ifdef STATISTICS
        bvestats.resolutions++;
#endif
      }
      //=======================
      // resolve or substitute
      //=======================
      if (type & SUBSTITUTION) {
        if (nAddedCls) xsubstitute(v, out_c);
        toblivion(p, pOrgs, nOrgs, poss, negs, model);
        v = 0;
      } else if (type & RESOLUTION) {
        if (nAddedCls) xresolve(v, out_c);
        toblivion(p, pOrgs, nOrgs, poss, negs, model);
        v = 0;
      }
    }
    if (!v) {
      markEliminated(PVs[i]);
      PVs[i] = 0;
    }
  }
  if (opts.profile_simp) timer.pstop(), timer.ve += timer.pcpuTime();
}

inline void Solver::xsubstitute(const uint32& x, Lits_t& out_c) {
  CHECKVAR(x);
  PFLOG2(4, " Substituting(%d):", x);
  PFLOCCURS(solver, 4, x);
  uint32 dx = V2L(x), fx = NEG(dx);
  if (ot[dx].size() > ot[fx].size()) swap(dx, fx);
  OL &me = ot[dx], &other = ot[fx];
  forall_occurs(me, i) {
    SCLAUSE& ci = **i;
    if (ci.original()) {
      const bool a = ci.molten();
      forall_occurs(other, j) {
        SCLAUSE& cj = **j;
        if (cj.original()) {
          const bool b = cj.molten();
          if (NEQUAL(a, b) && merge(x, ci, cj, out_c))
            newResolvent(out_c);
        }
      }
    }
  }
}

inline void Solver::xresolve(const uint32& x, Lits_t& out_c) {
  CHECKVAR(x);
  PFLOG2(4, " Resolving(%d):", x);
  PFLOCCURS(solver, 4, x);
  uint32 dx = V2L(x), fx = NEG(dx);
  if (ot[dx].size() > ot[fx].size()) swap(dx, fx);
  OL &me = ot[dx], &other = ot[fx];
  forall_occurs(me, i) {
    SCLAUSE& ci = **i;
    if (ci.original()) {
      forall_occurs(other, j) {
        SCLAUSE& cj = **j;
        if (cj.original() && merge(x, ci, cj, out_c))
          newResolvent(out_c);
      }
    }
  }
}

inline void Solver::newResolvent(const Lits_t& resolvent) {
  const int size = resolvent.size();
  assert(size);
  const size_t bytes = hc_scsize + (size - 1) * hc_isize;
  S_REF added = (S_REF) new Byte[bytes];
  added->init(resolvent);
  assert(added->size() == size);
  assert(added->hasZero() < 0);
  assert(added->original());
  assert(added->isSorted());
  if (size == 1) {
    const uint32 unit = **added;
    const LIT_ST val = sp->value[unit];
    if (UNASSIGNED(val)) {
      enqueueUnit(unit);
    } else if (!val) {
      PFLOG2(2, "  BVE proved a contradiction");
      learnEmpty();
      killSolver();
    }
  } else {
    if (opts.proof_en)
      proof.addResolvent(*added);
    added->calcSig();
    added->markAdded();
    scnf.push(added);
    PFLCLAUSE(4, (*added), " Resolvent");
  }
}