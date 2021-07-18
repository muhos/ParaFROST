/***********************************************************************[logging.h]
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

#ifndef __LOGGING_
#define __LOGGING_

#include <cstdio>
#include <cstring>
#include "const.h"
#include "color.h"

#define RULELEN     92
#define PREFIX      "c "
#define UNDERLINE	"\u001b[4m"

#if defined(__linux__) || defined(__CYGWIN__)
#pragma GCC system_header
#endif

#define PUTCH(CH) putc(CH, stdout);

#define PRINT(FORMAT, ...) fprintf(stdout, FORMAT, ## __VA_ARGS__);

inline void REPCH(const char& ch, const size_t& size, const size_t& off = 0) {
    for (size_t i = off; i < size; i++) PUTCH(ch);
}

#define PFLOGS(RESULT) \
    do { \
        if (!quiet_en) PFLOG0(""); \
        PRINT("s %s\n", RESULT); \
    } while (0)

#define PFLOGE(FORMAT, ...) \
  do { \
     SETCOLOR(CERROR, stderr); \
     fprintf(stderr, "ERROR - "); \
     fprintf(stderr, FORMAT, ## __VA_ARGS__); \
     putc('\n', stderr); \
     SETCOLOR(CNORMAL, stderr); \
     exit(1); \
  } while (0)

#define PFLOGEN(FORMAT, ...) \
  do { \
     SETCOLOR(CERROR, stderr); \
     fprintf(stderr, "ERROR - "); \
     fprintf(stderr, FORMAT, ## __VA_ARGS__); \
     putc('\n', stderr); \
     SETCOLOR(CNORMAL, stderr); \
  } while (0)

#define PFLOGW(FORMAT, ...) \
  do { \
     SETCOLOR(CWARNING, stderr); \
     fprintf(stderr, "WARNING - ");\
     fprintf(stderr, FORMAT, ## __VA_ARGS__);\
     PUTCH('\n'); \
     SETCOLOR(CNORMAL, stderr); \
  } while (0)

#define PFLRULER(CH, TIMES) \
  do { \
     PUTCH(PREFIX[0]); \
     REPCH(CH, TIMES);      \
     PUTCH('\n'); \
  } while (0)

#define PFNAME(NAME, VER) \
  PFLRULER('-', RULELEN); \
  do { \
     const char* suffix = "SAT Solver (version "; \
     size_t len = strlen(NAME) + strlen(suffix) + strlen(VER) + 1; \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the title"); \
     size_t gap = (RULELEN - len - 3) / 2; \
     PRINT(PREFIX); PUTCH(' '); \
     REPCH(' ', gap); \
     PRINT("%s%s%s %s%s%s)%s", UNDERLINE, CSOLVER, NAME, CSOLVER, suffix, VER, CNORMAL); \
     REPCH(' ', RULELEN + 1, len + gap + 3), PUTCH('\n'); \
  } while (0); \

#define PFAUTHORS(AUTHORS) \
  do { \
     const char *prefix = "Authors: "; \
     size_t len = strlen(prefix) + strlen(AUTHORS); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the authors"); \
     size_t gap = RULELEN - len - 1; \
     PRINT(PREFIX); PUTCH(' '); \
     PRINT("%s%s%s%s", CAUTHOR, prefix, AUTHORS, CNORMAL); \
     REPCH(' ', gap), PUTCH('\n'); \
  } while (0); \

#define PFRIGHTS(RIGHTS) \
  do { \
     const char *suffix = ", all rights reserved"; \
     size_t len = strlen(RIGHTS) + strlen(suffix); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the rights"); \
     size_t gap = RULELEN - len - 1; \
     PRINT(PREFIX); PUTCH(' ');\
     PRINT("%s%s%s%s", CRIGHTS, RIGHTS, CNORMAL, suffix); \
     REPCH(' ', gap), PUTCH('\n'); \
  } while (0); \

#define PFLOG0(MESSAGE) do { PRINT(PREFIX); PRINT("%s\n", MESSAGE); } while (0)

#define PFLOGN0(MESSAGE) do { PRINT(PREFIX); PRINT("%s", MESSAGE); } while (0)

#define PFLOG1(FORMAT, ...) \
    do { \
        PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); PUTCH('\n'); \
    } while (0)

#define PFLOGN1(FORMAT, ...) \
    do { \
        PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); \
    } while (0)

#define PFLOG2(VERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); PUTCH('\n'); } \
    } while (0)

#define PFLOGN2(VERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); } \
    } while(0)

#define PFPRINT(VERBOSITY, MAXVERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) { PRINT(FORMAT, ## __VA_ARGS__); } \
    } while (0)

#define PFLDONE(VERBOSITY, MAXVERBOSITY) if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) PRINT("done.\n");

#define PFLENDING(VERBOSITY, MAXVERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) { \
            PRINT(FORMAT, ## __VA_ARGS__); \
            PRINT(" done.\n"); \
        } \
    } while(0)

#define PFLMEMCALL(SOLVER, VERBOSITY) PFLOG2(VERBOSITY, " Memory used in %s call = %lld MB", __func__, sysMemUsed() / MBYTE);

#define PFLGCMEM(VERBOSITY, oldB, newB) \
    if (verbose >= VERBOSITY) { \
        double diff = abs(double(oldB.size() - newB.size())) * oldB.bucket(); \
        PRINT("(%.2f KB saved) ", diff / KBYTE); \
    }

#define PFLREDALL(SOLVER, VERBOSITY, MESSAGE) \
    if (verbose >= VERBOSITY) { \
        SOLVER->evalReds(); \
        PFLOG1("\t\t %s%s%s", CLBLUE, MESSAGE, CNORMAL); \
        SOLVER->logReductions(); }

#define PFLREDCL(SOLVER, VERBOSITY, MESSAGE) \
    if (verbose >= VERBOSITY) { \
        inf.n_del_vars_after = 0; \
        SOLVER->countAll(); \
        PFLOG1("\t\t %s%s%s", CLBLUE, MESSAGE, CNORMAL); \
        SOLVER->logReductions(); }

#define PFORGINF(SOLVER, CLS, LITS) \
    int64 CLS = SOLVER->stats.clauses.original; \
    int64 LITS = SOLVER->stats.literals.original; \

#define PFLEARNTINF(SOLVER, CLS, LITS) \
    int64 CLS = SOLVER->stats.clauses.learnt; \
    int64 LITS = SOLVER->stats.literals.learnt; \

#define PFLSHRINKALL(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - maxClauses(), RLITS = BLITS - maxLiterals(); \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        PFLENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#define PFLSHRINKORG(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - SOLVER->stats.clauses.original, RLITS = BLITS - SOLVER->stats.literals.original; \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        PFLENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#define PFLSHRINKLEARNT(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - SOLVER->stats.clauses.learnt, RLITS = BLITS - SOLVER->stats.literals.learnt; \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        PFLENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#ifdef LOGGING

#define PFLBCPS(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		PFLOG1("\t Before BCP(%d)", l2i(LIT)); \
		SOLVER->printWatched(ABS(LIT)); }

#define PFLBCP(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		PFLOG1("\t BCP(%d)", l2i(LIT)); \
		SOLVER->printOL(LIT); \
        SOLVER->printOL(FLIP(LIT)); }

#define PFLTRAIL(SOLVER, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printTrail();

#define PFLLEARNT(SOLVER, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printLearnt();

#define PFLSORTED(SOLVER, SIZE, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printSortedStack(SIZE);

#define PFLCLAUSE(VERBOSITY, CLAUSE, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { \
            PRINT(PREFIX);\
            PRINT(FORMAT, ## __VA_ARGS__); \
            SETCOLOR(CLOGGING, stdout);\
            CLAUSE.print(); \
            SETCOLOR(CNORMAL, stdout);\
        } \
    } while (0)

#define PFLDL(SOLVER, VERBOSITY) PFLOG2(VERBOSITY, " Current decision level: %d", SOLVER->DL());

#define PFLBCPE(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		PFLOG1("\t After BCP(%d)", l2i(LIT)); \
		SOLVER->printWatched(ABS(LIT)); \
        PFLRULER('-', 30); }

#define PFLOCCURS(SOLVER, VERBOSITY, VAR) \
     if (verbose >= VERBOSITY) { \
		PFLOG1("\t Full Occurrence LIST(%d)", VAR); \
		SOLVER->printOccurs(VAR); \
        PFLRULER('-', 30); }

#define PFLNEWLIT(SOLVER, VERBOSITY, SRC, LIT) \
    do { \
        PFLOG2(VERBOSITY, "   %sNew %s( %d@%d )%s", CREPORTVAL, SRC == NOREF ? !SOLVER->DL() ? "forced unit" : "decision" : "unit", l2i(LIT), l2dl(LIT), CNORMAL); \
    } while (0)

#define PFLCONFLICT(SOLVER, VERBOSITY, LIT) \
    do { \
        PFLOG2(VERBOSITY, " %sConflict detected in literal( %d@%d )%s", CCONFLICT, l2i(LIT), l2dl(LIT), CNORMAL); \
    } while (0)

#else // NO LOGGING

#define PFLBCPS(SOLVER, VERBOSITY, LIT) do { } while (0)

#define PFLBCP(SOLVER, VERBOSITY, LIT) do { } while (0)

#define PFLTRAIL(SOLVER, VERBOSITY) do { } while (0)

#define PFLLEARNT(SOLVER, VERBOSITY) do { } while (0)

#define PFLSORTED(SOLVER, SIZE, VERBOSITY) do { } while (0)

#define PFLCLAUSE(VERBOSITY, CLAUSE, FORMAT, ...) do { } while (0)

#define PFLDL(SOLVER, VERBOSITY) do { } while (0)

#define PFLBCPE(SOLVER, VERBOSITY, LIT) do { } while (0)

#define PFLOCCURS(SOLVER, VERBOSITY, VAR) do { } while (0)

#define PFLNEWLIT(SOLVER, VERBOSITY, SRC, LIT) do { } while (0)

#define PFLCONFLICT(SOLVER, VERBOSITY, LIT) do { } while (0)
    
#endif // NO LOGGING

#endif