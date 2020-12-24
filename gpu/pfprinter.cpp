/***********************************************************************[pfprinter.cpp]
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

#include "pfsolve.h"
using namespace pFROST;

void ParaFROST::printTable()
{
	const char* header = " Progress ";
	size_t len = strlen(header);
	if (RULELEN < len) PFLOGE("ruler length is smaller than the table title");
	size_t gap = (RULELEN - strlen(header)) / 2;
	PFLOGN0(""); REPCH('-', gap);
	fprintf(stdout, "%s%s%s", CHEADER, header, CNORMAL);
	REPCH('-', gap); putc('|', stdout), putc('\n', stdout);
	string h = "";
	h = "              ORG                Conflicts   Restarts            Learnt           V(%)";
	if (RULELEN < h.size()) PFLOGE("ruler length is smaller than the table header");
	PFLOGN0("");
	fprintf(stdout, "%s%s%s", CHEADER, h.c_str(), CNORMAL);
	REPCH(' ', RULELEN - h.size()), putc('|', stdout), putc('\n', stdout);
	h = "      V         C         L                                C          L     L/C";
	if (RULELEN < h.size()) PFLOGE("ruler length is smaller than the table header");
	PFLOGN0("");
	fprintf(stdout, "%s%s%s", CHEADER, h.c_str(), CNORMAL);
	REPCH(' ', RULELEN - h.size()), putc('|', stdout), putc('\n', stdout);
	solLine = "  %9d %9d %10d %10d %8d %9d %10d %6d  %3d%s", solLineLen = 86;
	if (RULELEN < solLineLen) PFLOGE("ruler length is smaller than the progress line");
	PFLRULER('-', RULELEN);
}

void ParaFROST::printStats(const bool& _p, const Byte& _t, const char* _c)
{
	if (verbose == 1 && _p) {
		int l2c = (int)ratio(inf.nLearntLits, learnts.size());
		int vr = int(100.0 * double(maxActive()) / double(inf.orgVars));
		solLine[0] = _t;
		PFLOGN0("");
		SETCOLOR(_c, stdout);
		fprintf(stdout, solLine.c_str(),
			maxActive(), orgs.size(), inf.nLiterals,
			nConflicts, starts - 1, learnts.size(), inf.nLearntLits,
			l2c, vr, "%");
		SETCOLOR(CNORMAL, stdout);
		REPCH(' ', RULELEN - solLineLen), putc('|', stdout), putc('\n', stdout);
	}
}

void ParaFROST::printVars(const uint32* arr, const uint32& size, const LIT_ST& type)
{
	printf("(size = %d)->[", size);
	for (uint32 i = 0; i < size; i++) {
		if (type == 'l') printf("%4d", l2i(arr[i]));
		else if (type == 'v') printf("%4d", arr[i]);
		else printf("%4d", arr[i] + 1);
		printf("  ");
		if (i && i < size - 1 && i % 10 == 0) { putc('\n', stdout); PFLOGN0("\t\t\t"); }
	}
	putc(']', stdout), putc('\n', stdout);
}

void ParaFROST::printClause(const Lits_t& c)
{
	putc('(', stdout);
	for (int i = 0; i < c.size(); i++)
		printf("%4d", l2i(c[i]));
	putc(')', stdout), putc('\n', stdout);
}

void ParaFROST::printTrail(const uint32& off)
{
	if (trail.empty()) return;
	PFLOGN1(" Trail (size = %d)->[", trail.size());
	for (uint32 i = off; i < trail.size(); i++) {
		printf("%4d@%-4d", l2i(trail[i]), l2dl(trail[i]));
		if (i && i < trail.size() - 1 && i % 8 == 0) { putc('\n', stdout); PFLOGN0("\t\t\t"); }
	}
	putc(']', stdout), putc('\n', stdout);
}

void ParaFROST::printCNF(const BCNF& cnf, const int& off)
{
	PFLOG1("\tHost CNF(size = %d)", cnf.size());
	for (uint32 i = off; i < cnf.size(); i++) {
		if (cm[cnf[i]].size()) {
			PFLCLAUSE(1, cm[cnf[i]], " C(%d)->", i);
		}
	}
}

void ParaFROST::printOL(const OL& list, const bool& host) {
	CNF* cs = host ? hcnf : cnf;
	for (uint32 i = 0; i < list.size(); i++) {
		PFLOGN1(" C(%lld)->", uint64(list[i]));
		(*cs)[list[i]].print();
	}
}

void ParaFROST::printOL(const uint32& lit, const bool& host) {
	if ((*ot)[lit].size() == 0) return;
	PFLOG1(" List(%d):", l2i(lit));
	printOL((*ot)[lit], host);
}

void ParaFROST::printWL(const uint32& lit, const bool& bin)
{
	CHECKLIT(lit);
	const WL& ws = wt[lit];
	if (ws.size()) PFLOG1("  list(%d):", -l2i(lit));
	for (int i = 0; i < ws.size(); i++) {
		if (!ws[i].binary() && bin) continue;
		PFLCLAUSE(1, cm[ws[i].ref], "  %sW(r: %-4zd, sz: %-4d, i: %-4d)->%s",
			CLOGGING, ws[i].ref, ws[i].size, l2i(ws[i].imp), CNORMAL);
	}
}

void ParaFROST::printWatched(const uint32& v)
{
	CHECKVAR(v);
	uint32 p = V2L(v);
	printWL(p), printWL(NEG(p));
}

void ParaFROST::printBinaries(const uint32& v)
{
	CHECKVAR(v);
	uint32 p = V2L(v);
	printWL(p, 1), printWL(NEG(p), 1);
}

void ParaFROST::printWT()
{
	PFLOG0(" Watches:");
	for (uint32 lit = 2; lit < wt.size(); lit++)
		printWL(lit);
}

void ParaFROST::printOT()
{
	PFLOG0(" Positive occurs:");
	forall_variables(v) {
		printOL(V2L(v));
	}
	PFLOG0(" Negative occurs:");
	forall_variables(v) {
		printOL(NEG(V2L(v)));
	}
}

void ParaFROST::printHeap()
{
	PFLOG1(" Heap (size = %d):", vsids.size());
	for (uint32 i = 0; i < vsids.size(); i++)
		PFLOG1(" h(%d)->(v: %d, a: %g)", i, vsids[i], activity[vsids[i]]);
}

void ParaFROST::printSource()
{
	for (uint32 i = 0; i < trail.size(); i++) {
		assert(trail[i] > 1);
		uint32 v = ABS(trail[i]);
		C_REF r = sp->source[v];
		if (REASON(r))
			PFLCLAUSE(1, cm[r], " Source(v:%d, r:%zd)->", v, r);
	}
}

void ParaFROST::printLearnt()
{
	PFLOGN1(" %sLearnt(", CCONFLICT);
	for (int i = 0; i < learntC.size(); i++)
		printf("%d@%-6d", l2i(learntC[i]), l2dl(learntC[i]));
	fprintf(stdout, ")%s\n", CNORMAL);
}