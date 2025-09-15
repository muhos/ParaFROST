/***********************************************************************[dimacs.cpp]
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
#include "dimacs.hpp"
#include "control.hpp"

using namespace ParaFROST;

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#if defined(__linux__) || defined(__CYGWIN__)
#define EOF_MSG "(exit: cntrl + D)"
#elif defined(_WIN32)
#define EOF_MSG "(exit: cntrl + Z)"
#endif

bool Solver::parser() 
{
	FAULT_DETECTOR;
	struct stat st;
	bool readfile = canAccess(formula.path.c_str(), st);
	const uint64 fsz = formula.size = st.st_size;
	if (readfile)
		LOG2(1, " Parsing CNF file \"%s%s%s\" (size: %s%lld MB%s)",
			CREPORTVAL, formula.path.c_str(), CNORMAL, CREPORTVAL, ratio(fsz, uint64(MBYTE)), CNORMAL);
	else {
		formula.size = st.st_size = 0;
		formula.path.assign("stdin");
		LOG2(1, " Reading DIMACS from \"%sstdin%s\" " EOF_MSG "...", CREPORTVAL, CNORMAL);
	}
	timer.start();
	Lits_t in_c, org;
	if (readfile) {
		char* str = NULL;
#if defined(__linux__) || defined(__CYGWIN__)
		int fd = open(formula.path.c_str(), O_RDONLY, 0);
		if (fd == -1) LOGERROR("cannot open input file");
		void* buffer = mmap(NULL, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
		str = (char*) buffer;
#else
		ifstream inputFile;
		inputFile.open(formula.path, ifstream::in);
		if (!inputFile.is_open()) LOGERROR("cannot open input file");
		char* buffer = pfcalloc<char>(fsz + 1);
		str = buffer;
		inputFile.read(buffer, fsz);
		buffer[fsz] = '\0';
#endif
		in_c.reserve(INIT_CAP);
		org.reserve(INIT_CAP);
		// Initially 'inf' is reset but if parser() is called again
		// by an incremental solving procedure, then 'inf' has to be
		// manually reset
		RESETSTRUCT(&inf);
		char* eof = str + fsz;
		while (str < eof) {
			eatWS(str);
			if (*str == '\0' || *str == '%') break;
			if (*str == 'c') eatLine(str);
			else if (*str == 'p') {
				if (!eq(str, "p cnf")) LOGERROR("header has wrong format");
				uint32 sign = 0;
				inf.orgVars = toInteger(str, sign);
				if (!opts.parseincr_en) {
					inf.unassigned = inf.maxVar = inf.orgVars;
					inf.nDualVars = V2L(inf.orgVars + 1);
				}
				if (sign) LOGERROR("number of variables in header is negative");
				if (inf.orgVars == 0) LOGERROR("zero number of variables in header");
				if (inf.orgVars >= INT_MAX - 1) LOGERROR("number of variables not supported");
				inf.orgCls = toInteger(str, sign);
				if (sign) LOGERROR("number of clauses in header is negative");
				if (inf.orgCls == 0) LOGERROR("zero number of clauses in header");
				LOG2(1, " Found header %s%d %d%s", CREPORTVAL, inf.orgVars, inf.orgCls, CNORMAL);
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
			else if (opts.parseincr_en) {
				incremental = true;
				uint32 v = 0, s = 0;
				while ((v = toInteger(str, s)) != 0) {
					while (v > inf.maxVar) iadd();
					org.push(V2DEC(v, s));
				}
				if (!itoClause(in_c, org)) return false;
			}
			else if (!toClause(in_c, org, str)) return false;
		}
#if defined(__linux__) || defined(__CYGWIN__)
		if (munmap(buffer, fsz) != 0) LOGERROR("cannot clean input file %s mapping", formula.path.c_str());
		close(fd);
#else
		delete[] buffer;
		inputFile.close();
#endif
	}
	else {
		in_c.reserve(INIT_CAP);
		org.reserve(INIT_CAP);
		int ch;
		formula.eatComment(ch);
		// read CNF header
		if (ch != 'p') 
			LOGERROR("expected 'c' or 'p'");
		ch = formula.get();
		if (!isSpace(ch))
			LOGERROR("expected space after 'p'");
		else 
			formula.eatWS(ch);
		if (ch == 'c') {
			if ((ch = formula.get()) != 'n') 
				LOGERROR("expected 'n' after 'p c'");
			if ((ch = formula.get()) != 'f') 
				LOGERROR("expected 'f' after 'p cn'");
			ch =  formula.get();
			if (!isSpace(ch)) 
				LOGERROR("expected space after 'p cnf'");
			formula.eatWS(ch);
			if (isSpace(ch)) LOGERROR("expected digit after 'p cnf '");
			uint32 sign = 0;
			inf.orgVars = formula.toInteger(ch, sign);
			if (!opts.parseincr_en) {
				inf.unassigned = inf.maxVar = inf.orgVars;
				inf.nDualVars = V2L(inf.orgVars + 1);
			}
			if (sign) LOGERROR("number of variables in header is negative");
			if (inf.orgVars == 0) LOGERROR("zero number of variables in header");
			if (inf.orgVars >= INT_MAX - 1) LOGERROR("number of variables not supported");
			if (!isSpace(ch)) 
				LOGERROR("expected space after 'p cnf %d'", inf.orgVars);
			formula.eatWS(ch);
			if (!isDigit(ch)) 
				LOGERROR("expected digit after 'p cnf %d '", inf.orgVars);
			inf.orgCls = formula.toInteger(ch, sign);
			if (sign) LOGERROR("number of clauses in header is negative");
			if (inf.orgCls == 0) LOGERROR("zero number of clauses in header");
			while (ch != '\n') {
				if (!isSpace(ch))
					LOGERROR("expected newline after 'p cnf %d %d'", inf.orgVars, inf.orgCls);
				ch = formula.get();
			}
			LOG2(1, " Found header %s%d %d%s", CREPORTVAL, inf.orgVars, inf.orgCls, CNORMAL);
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
		else
			LOGERROR("expected 'c' after 'p '");
		// now read clauses
		while ((ch = formula.get()) != EOF) {
			if (isSpace(ch)) continue;
			else if (opts.parseincr_en) {
				incremental = true;
				uint32 v = 0, s = 0;
				while ((v = formula.toInteger(ch, s)) != 0) {
					while (v > inf.maxVar) iadd();
					org.push(V2DEC(v, s));
				}
				if (!itoClause(in_c, org)) return false;
			}
			else if (!toClause(in_c, org, ch)) return false;
		}
	}
	assert(stats.clauses.original == orgs.size());
	assert(orgs.size() <= inf.orgCls);
	orgs.shrinkCap();
	in_c.clear(true), org.clear(true);
	timer.stop();
	timer.parse = timer.cpuTime();
	LOG2(1, " Read %s%d Variables%s, %s%d Clauses%s, and %s%lld Literals%s in %s%.2f seconds%s",
		CREPORTVAL, inf.maxVar, CNORMAL,
		CREPORTVAL, orgs.size() + trail.size(), CNORMAL,
		CREPORTVAL, stats.literals.original + trail.size(), CNORMAL,
		CREPORTVAL, timer.parse, CNORMAL);
	LOG2(1, "  found %s%d units%s, %s%d binaries%s, %s%d ternaries%s, %s%d larger%s", 
		CREPORTVAL, formula.units, CNORMAL, 
		CREPORTVAL, formula.binaries, CNORMAL, 
		CREPORTVAL, formula.ternaries, CNORMAL, 
		CREPORTVAL, formula.large, CNORMAL);
	LOG2(1, "  maximum clause size: %s%d%s", CREPORTVAL, formula.maxClauseSize, CNORMAL);
	return true;
}

bool Solver::toClause(Lits_t& c, Lits_t& org, char*& str)
{
	assert(c.empty());
	assert(org.empty());
	uint32 v = 0, s = 0;
	bool satisfied = false;

	while ((v = toInteger(str, s)) != 0) {
		if (v > inf.maxVar) LOGERROR("too many variables");
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
		int newsize = c.size();
		if (!newsize) {
			learnEmpty();
			LOG2(1, " Found empty clause");
			return false; 
		}
		else if (newsize == 1) {
			const uint32 unit = *c;
			CHECKLIT(unit);
			LIT_ST val = sp->value[unit];
			if (UNASSIGNED(val)) enqueueUnit(unit), formula.units++;
			else if (!val) {
				learnEmpty();
				LOG2(1, "  Found conflicting unit");
				return false;
			}
		}
		else if (orgs.size() + 1 > inf.orgCls) LOGERROR("too many clauses");
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

bool Solver::toClause(Lits_t& c, Lits_t& org, int& ch)
{
	assert(c.empty());
	assert(org.empty());
	uint32 v = 0, s = 0;
	bool satisfied = false;
	while ((v = formula.toInteger(ch, s)) != 0) {
		formula.eatWS(ch);
		if (v > inf.maxVar) LOGERROR("too many variables");
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
		int newsize = c.size();
		if (!newsize) {
			learnEmpty();
			LOG2(1, "  Found empty clause");
			return false;
		}
		else if (newsize == 1) {
			const uint32 unit = *c;
			CHECKLIT(unit);
			LIT_ST val = sp->value[unit];
			if (UNASSIGNED(val)) enqueueUnit(unit), formula.units++;
			else if (!val) {
				learnEmpty();
				LOG2(1, "  Found conflicting unit");
				return false;
			}
		}
		else if (orgs.size() + 1 > inf.orgCls) LOGERROR("too many clauses");
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

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic pop
#endif