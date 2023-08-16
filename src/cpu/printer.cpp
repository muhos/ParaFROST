/***********************************************************************[printer.cpp]
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

void Solver::printTable() {
  const char* header = " Progress ";
  size_t len = strlen(header);
  if (RULELEN < len) PFLOGE("ruler length is smaller than the table title");
  size_t gap = (RULELEN - strlen(header)) / 2;
  PUTCH('c');
  REPCH('-', gap);
  PRINT("%s%s%s", CHEADER, header, CNORMAL);
  REPCH('-', gap);
  PUTCH('\n');
  string h = "";
  h = "              ORG                Conflicts   Restarts            Learnt           V(%)  C(%)";
  if (RULELEN < h.size()) PFLOGE("ruler length is smaller than the table header");
  PUTCH('c');
  PRINT("%s%s%s", CHEADER, h.c_str(), CNORMAL);
  REPCH(' ', RULELEN - h.size());
  PUTCH('\n');
  h = "      V         C         L                                C          L     L/C";
  if (RULELEN < h.size())
    PFLOGE("ruler length is smaller than the table header");
  PUTCH('c');
  PRINT("%s%s%s", CHEADER, h.c_str(), CNORMAL);
  REPCH(' ', RULELEN - h.size());
  PUTCH('\n');
  solLine = "  %9d %9d %10d %10d %8d %9d %10d %6d  %3d%s %3d%s", solLineLen = 91;
  if (RULELEN < solLineLen)
    PFLOGE("ruler length is smaller than the progress line");
  PFLRULER('-', RULELEN);
}

void Solver::printStats(const bool& _p, const Byte& _t, const char* _c) {
  if (verbose == 1 && _p) {
    const int l2c = (int)ratio(stats.literals.learnt, stats.clauses.learnt);
    const int vr = (int)percent(maxActive(), inf.orgVars);
    const int cr = (int)percent((double)stats.clauses.original, inf.nOrgCls);
    solLine[0] = _t;
    PFLOGN0("");
    SETCOLOR(_c, stdout);
    PRINT(solLine.c_str(),
          maxActive(), stats.clauses.original, stats.literals.original,
          stats.conflicts, stats.restart.all, stats.clauses.learnt, stats.literals.learnt,
          l2c, vr, "%", cr, "%");
    SETCOLOR(CNORMAL, stdout);
    REPCH(' ', RULELEN - solLineLen);
    PUTCH('\n');
  }
}

void Solver::printVars(const uint32* arr, const uint32& size, const LIT_ST& type) {
  PRINT("(size = %d)->[", size);
  for (uint32 i = 0; i < size; i++) {
    if (type == 'l') {
      PRINT("%4d", l2i(arr[i]));
    } else if (type == 'v') {
      PRINT("%4d", arr[i]);
    } else {
      PRINT("%4d", arr[i] + 1);
    }
    PRINT("  ");
    if (i && i < size - 1 && i % 10 == 0) {
      PUTCH('\n');
      PFLOGN0("\t\t\t");
    }
  }
  PUTCH(']');
  PUTCH('\n');
}

void Solver::printClause(const Lits_t& c) {
  PUTCH('(');
  for (int i = 0; i < c.size(); i++)
    PRINT("%4d", l2i(c[i]));
  PUTCH(')');
  PUTCH('\n');
}

void Solver::printTrail(const uint32& off) {
  if (trail.empty()) return;
  PFLOGN1(" Trail (size = %d)->[", trail.size());
  for (uint32 i = off; i < trail.size(); i++) {
    PRINT("%4d@%-4d", l2i(trail[i]), l2dl(trail[i]));
    if (i && i < trail.size() - 1 && i % 8 == 0) {
      PUTCH('\n');
      PFLOGN0("\t\t\t");
    }
  }
  PUTCH(']');
  PUTCH('\n');
}

void Solver::printCNF(const BCNF& cnf, const int& off) {
  PFLOG1("\tHost CNF(size = %d)", cnf.size());
  for (uint32 i = off; i < cnf.size(); i++) {
    if (cm[cnf[i]].size()) {
      PFLCLAUSE(1, cm[cnf[i]], " C(%d)->", i);
    }
  }
}

void Solver::printCNF(const SCNF& cnf, const int& off, const char* t) {
  PFLOG1("\tHost CNF(size = %zd)", cnf.size());
  if (!strcmp(t, "added")) {
    for (uint32 c = off; c < cnf.size(); c++) {
      if (cnf[c]->size() && cnf[c]->added()) {
        PFLOGN1(" C(%d)->", c);
        cnf[c]->print();
      }
    }
  } else if (!strcmp(t, "deleted")) {
    for (uint32 c = off; c < cnf.size(); c++) {
      if (cnf[c]->size() && cnf[c]->deleted()) {
        PFLOGN1(" C(%d)->", c);
        cnf[c]->print();
      }
    }
  } else {
    for (uint32 c = off; c < cnf.size(); c++) {
      if (cnf[c]->size()) {
        PFLOGN1(" C(%d)->", c);
        cnf[c]->print();
      }
    }
  }
}

void Solver::printOL(const OL& list) {
  for (int i = 0; i < list.size(); i++) {
    PFLOGN0(" ");
    list[i]->print();
  }
}

void Solver::printOL(const uint32& lit) {
  CHECKLIT(lit);
  if (ot[lit].empty()) return;
  PFLOG1(" List(%d):", l2i(lit));
  printOL(ot[lit]);
}

void Solver::printOccurs(const uint32& v) {
  CHECKVAR(v);
  uint32 p = V2L(v);
  printOL(p), printOL(NEG(p));
}

void Solver::printWL(const uint32& lit, const bool& bin) {
  CHECKLIT(lit);
  const WL& ws = wt[lit];
  if (ws.size()) PFLOG1("  list(%d):", -l2i(lit));
  for (int i = 0; i < ws.size(); i++) {
    if (!ws[i].binary() && bin) continue;
    PFLCLAUSE(1, cm[ws[i].ref], "  %sW(r: %-4zd, sz: %-4d, i: %-4d)->%s",
              CLOGGING, ws[i].ref, ws[i].size, l2i(ws[i].imp), CNORMAL);
  }
}

void Solver::printWL(const WL& ws, const bool& bin) {
  for (int i = 0; i < ws.size(); i++) {
    if (!ws[i].binary() && bin) continue;
    PFLCLAUSE(1, cm[ws[i].ref], "  %sW(r: %-4zd, sz: %-4d, i: %-4d)->%s",
              CLOGGING, ws[i].ref, ws[i].size, l2i(ws[i].imp), CNORMAL);
  }
}

void Solver::printWatched(const uint32& v) {
  CHECKVAR(v);
  uint32 p = V2L(v);
  printWL(p), printWL(NEG(p));
}

void Solver::printBinaries(const uint32& v) {
  CHECKVAR(v);
  uint32 p = V2L(v);
  printWL(p, 1), printWL(NEG(p), 1);
}

void Solver::printWT() {
  PFLOG0(" Watches:");
  for (uint32 lit = 2; lit < wt.size(); lit++)
    printWL(lit);
}

void Solver::printOT() {
  PFLOG0(" Positive occurs:");
  forall_variables(v) {
    printOL(V2L(v));
  }
  PFLOG0(" Negative occurs:");
  forall_variables(v) {
    printOL(NEG(V2L(v)));
  }
}

void Solver::printHeap() {
  PFLOG1(" Heap (size = %d):", vsids.size());
  for (uint32 i = 0; i < vsids.size(); i++)
    PFLOG1(" h(%d)->(v: %d, a: %g)", i, vsids[i], activity[vsids[i]]);
}

void Solver::printSource() {
  for (uint32 i = 0; i < trail.size(); i++) {
    assert(trail[i] > 1);
    uint32 v = ABS(trail[i]);
    C_REF r = sp->source[v];
    if (REASON(r))
      PFLCLAUSE(1, cm[r], " Source(v:%d, r:%zd)->", v, r);
  }
}

void Solver::printLearnt() {
  PFLOGN1(" %sLearnt(", CCONFLICT);
  for (int i = 0; i < learntC.size(); i++)
    PRINT("%d@%-6d", l2i(learntC[i]), l2dl(learntC[i]));
  PRINT(")%s\n", CNORMAL);
}

void Solver::printSortedStack(const int& tail) {
  if (!tail) return;
  if (vhist.empty()) return;
  PFLOGN1("  %ssorted(", CMAGENTA);
  for (int i = 0; i < tail; i++) {
    const uint32 lit = sp->tmpstack[i];
    PRINT("%d:%-4d", l2i(lit), vhist[lit]);
  }
  PRINT(")%s\n", CNORMAL);
}