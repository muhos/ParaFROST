/***********************************************************************[printer.cpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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

#include "solver.hpp"

using namespace ParaFROST;

void Solver::printTable()
{
	string h = "";
	h = "              ORG                Conflicts     Restarts            Learnt         V(%)  C(%)";
	if (RULELEN < h.size()) LOGERROR("ruler length is smaller than the table header");
	PUTCH('c');
	PRINT("%s%s%s", CHEADER, h.c_str(), CNORMAL);
	REPCH(' ', RULELEN - h.size());
	PUTCH('\n');
	h = "      V         C         L                                C          L       L/C";
	if (RULELEN < h.size())
		LOGERROR("ruler length is smaller than the table header");
	PUTCH('c');
	PRINT("%s%s%s", CHEADER, h.c_str(), CNORMAL);
	REPCH(' ', RULELEN - h.size());
	PUTCH('\n');
    LOGRULER('-', RULELEN);
	solLine = " %9d %9d %10d %10d %8d %9d %10d %6d  %3d%s  %3d%s", solLineLen = 91;
	if (RULELEN < solLineLen)
		LOGERROR("ruler length is smaller than the progress line");
}

void Solver::printStats(const bool& _p, const Byte& _t, const char* _c)
{
	if (verbose == 1 && _p) {
		const int l2c = (int)ratio(stats.literals.learnt, stats.clauses.learnt);
		const int vr = (int)percent(maxActive(), inf.orgVars);
		const int cr = (int)percent((double)stats.clauses.original, inf.orgCls);
		solLine[0] = _t;
		LOGN0("");
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

void Solver::printVars(const uint32* arr, const uint32& size, const LIT_ST& type)
{
	PRINT("(size = %d)->[", size);
	for (uint32 i = 0; i < size; i++) {
		if (type == 'l') {
			PRINT("%4d", l2i(arr[i]));
		}
		else if (type == 'v') {
			PRINT("%4d", arr[i]);
		}
		else {
			PRINT("%4d", arr[i] + 1);
		}
		PRINT("  ");
		if (i && i < size - 1 && i % 10 == 0) { PUTCH('\n'); LOGN0("\t\t\t"); }
	}
	PUTCH(']');
	PUTCH('\n');
}

void Solver::printClause(const Lits_t& c)
{
	PUTCH('(');
	for (int i = 0; i < c.size(); i++)
		PRINT("%4d", l2i(c[i]));
	PUTCH(')'); PUTCH('\n');
}

void Solver::printTrail(const uint32& off)
{
	if (trail.empty()) return;
	LOGN1(" Trail (size = %d)->[", trail.size());
	for (uint32 i = off; i < trail.size(); i++) {
		PRINT("%4d@%-4d", l2i(trail[i]), l2dl(trail[i]));
		if (i && i < trail.size() - 1 && i % 8 == 0) { PUTCH('\n'); LOGN0("\t\t\t"); }
	}
	PUTCH(']'); PUTCH('\n');
}

void Solver::printCNF(const BCNF& cnf, const int& off)
{
	LOG1("\tHost CNF(size = %d)", cnf.size());
	for (uint32 i = off; i < cnf.size(); i++) {
		if (cm[cnf[i]].size()) {
			LOGCLAUSE(1, cm[cnf[i]], " C(%d)->", i);
		}
	}
}

void Solver::printOL(const OL& list, const bool& host) {
	CNF* cs = host ? hcnf : cnf;
	for (uint32 i = 0; i < list.size(); i++) {
		LOGN1(" C(%lld)->", uint64(list[i]));
		(*cs)[list[i]].print();
	}
}

void Solver::printOL(const uint32& lit, const bool& host) {
	if ((*ot)[lit].size() == 0) return;
	LOG1(" List(%d):", l2i(lit));
	printOL((*ot)[lit], host);
}

void Solver::printWL(const uint32& lit, const bool& bin)
{
	CHECKLIT(lit);
	const WL& ws = wt[lit];
	if (ws.size()) LOG1("  list(%d):", -l2i(lit));
	for (int i = 0; i < ws.size(); i++) {
		if (!ws[i].binary() && bin) continue;
		LOGCLAUSE(1, cm[ws[i].ref], "  %sW(r: %-4zd, sz: %-4d, i: %-4d)->%s",
			CLOGGING, ws[i].ref, ws[i].size, l2i(ws[i].imp), CNORMAL);
	}
}

void Solver::printWL(const WL& ws, const bool& bin)
{
	for (int i = 0; i < ws.size(); i++) {
		if (!ws[i].binary() && bin) continue;
		LOGCLAUSE(1, cm[ws[i].ref], "  %sW(r: %-4zd, sz: %-4d, i: %-4d)->%s",
			CLOGGING, ws[i].ref, ws[i].size, l2i(ws[i].imp), CNORMAL);
	}
}

void Solver::printWatched(const uint32& v)
{
	CHECKVAR(v);
	uint32 p = V2L(v);
	printWL(p), printWL(NEG(p));
}

void Solver::printBinaries(const uint32& v)
{
	CHECKVAR(v);
	uint32 p = V2L(v);
	printWL(p, 1), printWL(NEG(p), 1);
}

void Solver::printWT()
{
	LOG0(" Watches:");
	for (uint32 lit = 2; lit < wt.size(); lit++)
		printWL(lit);
}

void Solver::printHeap()
{
	LOG1(" Heap (size = %d):", vsids.size());
	for (uint32 i = 0; i < vsids.size(); i++)
		LOG1(" h(%d)->(v: %d, a: %g)", i, vsids[i], activity[vsids[i]]);
}

void Solver::printSource()
{
	for (uint32 i = 0; i < trail.size(); i++) {
		assert(trail[i] > 1);
		uint32 v = ABS(trail[i]);
		C_REF r = sp->source[v];
		if (REASON(r))
			LOGCLAUSE(1, cm[r], " Source(v:%d, r:%zd)->", v, r);
	}
}

void Solver::printLearnt()
{
	LOGN1(" %sLearnt(", CCONFLICT);
	for (int i = 0; i < learntC.size(); i++)
		PRINT("%d@%-6d", l2i(learntC[i]), l2dl(learntC[i]));
	PRINT(")%s\n", CNORMAL);
}

void Solver::printSortedStack(const int& tail)
{
	if (!tail) return;
	if (vhist.empty()) return;
	LOGN1("  %ssorted(", CMAGENTA);
	for (int i = 0; i < tail; i++) {
		const uint32 lit = sp->tmpstack[i];
		PRINT("%d:%-4d", l2i(lit), vhist[lit]);
	}
	PRINT(")%s\n", CNORMAL);
}