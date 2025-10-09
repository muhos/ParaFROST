/***********************************************************************[logging.hpp]
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

#ifndef __LOGGING_
#define __LOGGING_

#include <cstdio>
#include <cstring>
#include <string>
#include "constants.hpp"
#include "color.hpp"

#define STARTLEN    10
#define RULELEN     92
#define PREFIX      "c "
#define UNDERLINE	"\u001b[4m"

#if defined(__linux__) || defined(__CYGWIN__)
#pragma GCC system_header
#endif

#define PUTCH(CH, ...) putc(CH, stdout);

#define PRINT(FORMAT, ...) fprintf(stdout, FORMAT, ## __VA_ARGS__);

#define LOGRULER(CH, TIMES) \
  do { \
     PUTCH(PREFIX[0]); \
     REPCH(CH, TIMES);      \
     PUTCH('\n'); \
  } while (0)

#define LOGERROR(FORMAT, ...) \
  do { \
     SETCOLOR(CERROR, stderr); \
     fprintf(stderr, "ERROR - "); \
     fprintf(stderr, FORMAT, ## __VA_ARGS__); \
     putc('\n', stderr); \
     SETCOLOR(CNORMAL, stderr); \
     exit(1); \
  } while (0)

#define LOGERRORN(FORMAT, ...) \
  do { \
     SETCOLOR(CERROR, stderr); \
     fprintf(stderr, "ERROR - "); \
     fprintf(stderr, FORMAT, ## __VA_ARGS__); \
     putc('\n', stderr); \
     SETCOLOR(CNORMAL, stderr); \
  } while (0)

#define LOGWARNING(FORMAT, ...) \
  do { \
     SETCOLOR(CWARNING, stderr); \
     fprintf(stderr, "WARNING - ");\
     fprintf(stderr, FORMAT, ## __VA_ARGS__);\
     putc('\n', stderr); \
     SETCOLOR(CNORMAL, stderr); \
  } while (0)

inline void REPCH(const char& ch, const size_t& size, const size_t& off = 0) {
    for (size_t i = off; i < size; i++) PUTCH(ch);
}

#define LOGSAT(RESULT) \
    do { \
        if (!quiet_en) LOG0(""); \
        SETCOLOR(CWHITE, stdout); \
        PRINT("s %s\n", RESULT); \
        SETCOLOR(CNORMAL, stdout); \
    } while (0)

#define LOGHEADER(VERBOSITY, MAXVERBOSITY, HEAD) \
  do { \
    if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) { \
	  size_t len = strlen(HEAD) + 4; \
	  if (RULELEN < len) LOGERROR("ruler length is smaller than header line (%zd)", len); \
      SETCOLOR(CNORMAL, stdout); \
	  PUTCH('c'); \
	  REPCH('-', STARTLEN); \
	  PRINT("[ %s ]", HEAD); \
	  REPCH('-', (RULELEN - len - STARTLEN)); \
	  PUTCH('\n'); \
    } \
  } while (0); \

#define LOG0(MESSAGE) do { PRINT(PREFIX); PRINT("%s\n", MESSAGE); } while (0)

#define LOGN0(MESSAGE) do { PRINT(PREFIX); PRINT("%s", MESSAGE); } while (0)

#define LOG1(FORMAT, ...) \
    do { \
        PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); PUTCH('\n'); \
    } while (0)

#define LOGN1(FORMAT, ...) \
    do { \
        PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); \
    } while (0)

#define LOG2(VERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); PUTCH('\n'); } \
    } while (0)

#define LOGN2(VERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { PRINT(PREFIX); PRINT(FORMAT, ## __VA_ARGS__); } \
    } while(0)

#define PRINT2(VERBOSITY, MAXVERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) { PRINT(FORMAT, ## __VA_ARGS__); } \
    } while (0)

#define LOGDONE(VERBOSITY, MAXVERBOSITY) if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) PRINT("done.\n");

#define LOGENDING(VERBOSITY, MAXVERBOSITY, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY && verbose < MAXVERBOSITY) { \
            PRINT(FORMAT, ## __VA_ARGS__); \
            PRINT(" done.\n"); \
        } \
    } while(0)

#define LOGMEMCALL(SOLVER, VERBOSITY) LOG2(VERBOSITY, " Memory used in %s call: %lld MB", __func__, sysMemUsed() / MBYTE);

#define LOGGCMEM(VERBOSITY, oldB, newB) \
    if (verbose >= VERBOSITY) { \
        double diff = abs(double(oldB.size() - newB.size())) * oldB.bucket(); \
        PRINT("(%.2f KB saved) ", diff / KBYTE); \
    }

#define LOGREDALL(SOLVER, VERBOSITY, MESSAGE) \
    if (verbose >= VERBOSITY) { \
        SOLVER->countAll(); \
        SOLVER->updateNumElims(); \
        inf.prevDeletedVars += inf.currDeletedVars; \
        LOG1("\t\t %s%s%s", CLBLUE, MESSAGE, CNORMAL); \
        SOLVER->logReductions(); }

#define LOGREDALLHOST(SOLVER, VERBOSITY, MESSAGE) \
    if (verbose >= VERBOSITY) { \
        SOLVER->countAll(1); \
        SOLVER->updateNumElims(); \
        LOG1("\t\t %s", MESSAGE); \
        SOLVER->logReductions(); }

#define LOGREDCL(SOLVER, VERBOSITY, MESSAGE) \
    if (verbose >= VERBOSITY) { \
        SOLVER->countAll(); \
        inf.currDeletedVars = 0; \
        LOG1("\t\t %s%s%s", CLBLUE, MESSAGE, CNORMAL); \
        SOLVER->logReductions(); }

#define CNFORGINF(SOLVER, CLS, LITS) \
    int64 CLS = SOLVER->stats.clauses.original; \
    int64 LITS = SOLVER->stats.literals.original; \

#define CNFLEARNTINF(SOLVER, CLS, LITS) \
    int64 CLS = SOLVER->stats.clauses.learnt; \
    int64 LITS = SOLVER->stats.literals.learnt; \

#define LOGSHRINKALL(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - maxClauses(), RLITS = BLITS - maxLiterals(); \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        LOGENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#define LOGSHRINKORG(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - SOLVER->stats.clauses.original, RLITS = BLITS - SOLVER->stats.literals.original; \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        LOGENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#define LOGSHRINKLEARNT(SOLVER, VERBOSITY, BCLS, BLITS) \
    do { \
        int64 RCLS = BCLS - SOLVER->stats.clauses.learnt, RLITS = BLITS - SOLVER->stats.literals.learnt; \
        SOLVER->stats.shrink.clauses += RCLS, SOLVER->stats.shrink.literals += RLITS; \
        LOGENDING(VERBOSITY, 5, "(-%lld clauses, -%lld literals)", RCLS, RLITS); \
    } while (0)

#ifdef LOGGING

#define LOGBCPS(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		LOG1("\t Before BCP(%d)", l2i(LIT)); \
		SOLVER->printWatched(ABS(LIT)); }

#define LOGBCP(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		LOG1("\t BCP(%d)", l2i(LIT)); \
		SOLVER->printOL(LIT); \
        SOLVER->printOL(FLIP(LIT)); }

#define LOGTRAIL(SOLVER, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printTrail();

#define LOGLEARNT(SOLVER, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printLearnt();

#define LOGSORTED(SOLVER, SIZE, VERBOSITY) if (verbose >= VERBOSITY) SOLVER->printSortedStack(SIZE);

#define LOGCLAUSE(VERBOSITY, CLAUSE, FORMAT, ...) \
    do { \
        if (verbose >= VERBOSITY) { \
            PRINT(PREFIX);\
            PRINT(FORMAT, ## __VA_ARGS__); \
            SETCOLOR(CLOGGING, stdout);\
            CLAUSE.print(); \
            SETCOLOR(CNORMAL, stdout);\
        } \
    } while (0)

#define LOGDL(SOLVER, VERBOSITY) LOG2(VERBOSITY, " Current decision level: %d", SOLVER->DL());

#define LOGBCPE(SOLVER, VERBOSITY, LIT) \
     if (verbose >= VERBOSITY) { \
		LOG1("\t After BCP(%d)", l2i(LIT)); \
		SOLVER->printWatched(ABS(LIT)); \
        LOGRULER('-', 30); }

#define LOGOCCURS(SOLVER, VERBOSITY, VAR) \
     if (verbose >= VERBOSITY) { \
		LOG1("\t Full Occurrence LIST(%d)", VAR); \
		SOLVER->printOccurs(VAR); \
        LOGRULER('-', 30); }

#define LOGNEWLIT(SOLVER, VERBOSITY, SRC, LIT) \
    do { \
        LOG2(VERBOSITY, "   %sNew %s( %d@%d )%s", CREPORTVAL, SRC == NOREF ? !SOLVER->DL() ? "forced unit" : "decision" : "unit", l2i(LIT), l2dl(LIT), CNORMAL); \
    } while (0)

#define LOGCONFLICT(SOLVER, VERBOSITY, LIT) \
    do { \
        LOG2(VERBOSITY, " %sConflict detected in literal( %d@%d )%s", CCONFLICT, l2i(LIT), l2dl(LIT), CNORMAL); \
    } while (0)

#else // NO LOGGING

#define LOGBCPS(SOLVER, VERBOSITY, LIT) do { } while (0)

#define LOGBCP(SOLVER, VERBOSITY, LIT) do { } while (0)

#define LOGTRAIL(SOLVER, VERBOSITY) do { } while (0)

#define LOGLEARNT(SOLVER, VERBOSITY) do { } while (0)

#define LOGSORTED(SOLVER, SIZE, VERBOSITY) do { } while (0)

#define LOGCLAUSE(VERBOSITY, CLAUSE, FORMAT, ...) do { } while (0)

#define LOGDL(SOLVER, VERBOSITY) do { } while (0)

#define LOGBCPE(SOLVER, VERBOSITY, LIT) do { } while (0)

#define LOGOCCURS(SOLVER, VERBOSITY, VAR) do { } while (0)

#define LOGNEWLIT(SOLVER, VERBOSITY, SRC, LIT) do { } while (0)

#define LOGCONFLICT(SOLVER, VERBOSITY, LIT) do { } while (0)

#endif // NO LOGGING

#endif