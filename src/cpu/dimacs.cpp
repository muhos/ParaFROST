/***********************************************************************[dimacs.cpp]
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
#include "dimacs.h"
#include "control.h"

using namespace pFROST;

bool ParaFROST::parser() 
{
	FAULT_DETECTOR;
	struct stat st;
	if (!canAccess(formula.path.c_str(), st)) PFLOGE("cannot access the input file");
	const uint64 fsz = formula.size = st.st_size;
	PFLOG2(1, " Parsing CNF file \"%s%s%s\" (size: %s%lld MB%s)",
		CREPORTVAL, formula.path.c_str(), CNORMAL, CREPORTVAL, ratio(fsz, uint64(MBYTE)), CNORMAL);
	timer.start();
#if defined(__linux__) || defined(__CYGWIN__)
	int fd = open(formula.path.c_str(), O_RDONLY, 0);
	if (fd == -1) PFLOGE("cannot open input file");
	void* buffer = mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
	char* str = (char*)buffer;
#else
	ifstream inputFile;
	inputFile.open(formula.path, ifstream::in);
	if (!inputFile.is_open()) PFLOGE("cannot open input file");
	char* buffer = pfcalloc<char>(fsz + 1), * str = buffer;
	inputFile.read(buffer, fsz);
	buffer[fsz] = '\0';
#endif
	Lits_t in_c, org;
	in_c.reserve(INIT_CAP);
	org.reserve(INIT_CAP);
	PFLMEMCALL(this, 2);
	char* eof = str + fsz;
	while (str < eof) {
		eatWS(str);
		if (*str == '\0' || *str == '0' || *str == '%') break;
		if (*str == 'c') eatLine(str);
		else if (*str == 'p') {
			if (!eq(str, "p cnf")) PFLOGE("header has wrong format");
			uint32 sign = 0;
			inf.orgVars = toInteger(str, sign);
			if (!opts.parseincr_en) {
				inf.unassigned = inf.maxVar = inf.orgVars;
				inf.nDualVars = V2L(inf.orgVars + 1);
			}
			if (sign) PFLOGE("number of variables in header is negative");
			if (inf.orgVars == 0) PFLOGE("zero number of variables in header");
			if (inf.orgVars >= INT_MAX - 1) PFLOGE("number of variables not supported");
			inf.nOrgCls = toInteger(str, sign);
			if (sign) PFLOGE("number of clauses in header is negative");
			if (inf.nOrgCls == 0) PFLOGE("zero number of clauses in header");
			PFLOG2(1, " Found header %s%d %d%s", CREPORTVAL, inf.orgVars, inf.nOrgCls, CNORMAL);
			assert(orgs.empty());
			if (!opts.parseincr_en) {
				allocSolver();
				initQueue();
				initHeap();
				initVars();
				assert(vorg.size() == inf.maxVar + 1);
				model.init(vorg);
				if (opts.proof_en) 
					proof.init(sp, vorg);
			}
		}
		else if (!toClause(in_c, org, str)) return false;
		else if (opts.parseincr_en) {
			incremental = true;
			uint32 v = 0, s = 0;
			while ((v = toInteger(str, s)) != 0) {
				while (v > inf.maxVar) iadd();
				org.push(V2DEC(v, s));
			}
			if (!itoClause(in_c, org)) return false;
		}
	}
#if defined(__linux__) || defined(__CYGWIN__)
	if (munmap(buffer, fsz) != 0) PFLOGE("cannot clean input file %s mapping", formula.path.c_str());
	close(fd);
#else
	free(buffer);
	inputFile.close();
#endif
	assert(stats.clauses.original == orgs.size());
	assert(orgs.size() <= inf.nOrgCls);
	orgs.shrinkCap();
	in_c.clear(true), org.clear(true);
	timer.stop();
	timer.parse = timer.cpuTime();
	PFLOG2(1, " Read %s%d Variables%s, %s%d Clauses%s, and %s%lld Literals%s in %s%.2f seconds%s",
		CREPORTVAL, inf.maxVar, CNORMAL,
		CREPORTVAL, orgs.size() + trail.size(), CNORMAL,
		CREPORTVAL, stats.literals.original + trail.size(), CNORMAL,
		CREPORTVAL, timer.parse, CNORMAL);
	PFLOG2(1, "  found %s%d units%s, %s%d binaries%s, %s%d ternaries%s, %s%d larger%s", 
		CREPORTVAL, formula.units, CNORMAL, 
		CREPORTVAL, formula.binaries, CNORMAL, 
		CREPORTVAL, formula.ternaries, CNORMAL, 
		CREPORTVAL, formula.large, CNORMAL);
	PFLOG2(1, "  maximum clause size: %s%d%s", CREPORTVAL, formula.maxClauseSize, CNORMAL);
	return true;
}

bool ParaFROST::toClause(Lits_t& c, Lits_t& org, char*& str)
{
	assert(c.empty());
	assert(org.empty());
	uint32 v = 0, s = 0;
	bool satisfied = false;
	while ((v = toInteger(str, s)) != 0) {
		if (v > inf.maxVar) PFLOGE("too many variables");
		uint32 lit = V2DEC(v, s);
		CHECKLIT(lit);
		org.push(lit);
		// checking literal
		LIT_ST marker = l2marker(lit);
		if (UNASSIGNED(marker)) {
			markLit(lit);
			LIT_ST val = sp->value[lit];
			if (UNASSIGNED(val)) c.push(lit);
			else if (val) satisfied = true;
		}
		else if (marker != SIGN(lit)) satisfied = true; // tautology
	}
	forall_clause(org, k) {
		unmarkLit(*k);
	}
	if (satisfied) {
		if (opts.proof_en) proof.deleteClause(org);
	}
	else {
		if (org.empty()) { 
			if (opts.proof_en) proof.addEmpty();
			return false; 
		}
		int newsize = c.size();
		if (newsize == 1) {
			const uint32 unit = *c;
			CHECKLIT(unit);
			LIT_ST val = sp->value[unit];
			if (UNASSIGNED(val)) enqueueUnit(unit), formula.units++;
			else if (!val) return false;
		}
		else if (orgs.size() + 1 > inf.nOrgCls) PFLOGE("too many clauses");
		else if (newsize) {
			if (newsize == 2) formula.binaries++;
			else if (newsize == 3) formula.ternaries++;
			else assert(newsize > 3), formula.large++;
			if (newsize > formula.maxClauseSize)
				formula.maxClauseSize = newsize;
			newClause(c, false);
		}
		if (opts.proof_en && newsize < org.size()) {
			proof.addClause(c);
			proof.deleteClause(org);
			org.clear();
		}
	}
	c.clear(), org.clear();
	return true;
}