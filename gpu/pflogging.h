/***********************************************************************[pflogging.h]
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
#include "pfconst.h"

#define RULELEN 90

inline void REPCH(const char& ch, const size_t& size, const size_t& off = 0) {
    for (size_t i = off; i < size; i++) putc(ch, stdout);
}

#define PFLDONE(lvl, max) if (verbose >= lvl && verbose < max) fprintf(stdout, "done.\n");

#define PFLENDING(lvl, max, format, ...) \
    do { \
        if (verbose >= lvl && verbose < max) { \
            fprintf(stdout, format, ## __VA_ARGS__); \
            fprintf(stdout, " done.\n"); \
        } \
    }while(0)

#define PFLMEMCALL(s, lvl) PFLOG2(lvl, " Memory used in %s call = %lld MB", __func__, s->sysMemUsed() / MBYTE);

#define PFLGCMEM(lvl, oldB, newB) \
    if (verbose >= lvl) { \
        double diff = abs(double(oldB.size() - newB.size())) * oldB.bucket(); \
        fprintf(stdout, "(%.2f KB saved) ", diff / KBYTE); \
    }

#define PFLOG0(msg) fprintf(stdout, "c |%s\n", msg);

#define PFLOG1(format, ...) \
    do { \
        fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__), putc('\n', stdout); \
    } while (0)

#define PFLOG2(lvl, format, ...) \
    do { \
        if (verbose >= lvl) { fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__), putc('\n', stdout); } \
    } while (0)

#define PFLOGN0(msg) fprintf(stdout, "c |%s", msg);

#define PFLOGN1(format, ...) \
    do { \
        fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__); \
    } while (0)

#define PFLOGN2(lvl, format, ...) \
    do { \
        if (verbose >= lvl) { fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__); } \
    } while(0)

#define PFLOGS(sat) \
    do { \
        if (!quiet_en) PFLOG0(""); \
        fprintf(stdout, "s %s\n", sat); \
    } while (0)

#define PFLOGE(format, ...) \
  do { \
     fprintf(stderr, "ERROR - "), fprintf(stderr, format, ## __VA_ARGS__), putc('\n', stdout); \
     exit(1); \
  } while (0)

#define PFLOGEN(format, ...) \
  do { \
     fprintf(stderr, "ERROR - "), fprintf(stderr, format, ## __VA_ARGS__), putc('\n', stdout); \
  } while (0)

#define PFLOGW(format, ...) \
  do { \
     fprintf(stderr, "c | WARNING - "), fprintf(stderr, format, ## __VA_ARGS__), putc('\n', stdout); \
  } while (0)

#define PFLOGR(ch, x) \
  do { \
     fprintf(stdout, "c |"); \
     REPCH(ch, x);      \
     fprintf(stdout, "|\n"); \
  } while (0)

#define PFNAME(name) \
  PFLOGR('-', RULELEN); \
  do { \
     const char* suffix = "SAT Solver"; \
     size_t len = strlen(name) + strlen(suffix); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the title"); \
     size_t gap = (RULELEN - len) / 2; \
     PFLOGN0(" "); \
     REPCH(' ', gap); \
     fprintf(stdout, "%s %s", name, suffix); \
     REPCH(' ', RULELEN + 1, len + gap + 3), fprintf(stdout, "|\n"); \
  } while (0); \

#define PFAUTHORS(authors) \
  do { \
     const char *prefix = "Authors: "; \
     size_t len = strlen(prefix) + strlen(authors); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the authors"); \
     size_t gap = RULELEN - len - 1; \
     PFLOGN1(" %s%s", prefix, authors); \
     REPCH(' ', gap), fprintf(stdout, "|\n"); \
  } while (0); \

#define PFRIGHTS(rights) \
  do { \
     const char *suffix = ", all rights reserved."; \
     size_t len = strlen(rights) + strlen(suffix); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the rights"); \
     size_t gap = RULELEN - len - 1; \
     PFLOGN1(" %s%s", rights, suffix); \
     REPCH(' ', gap), fprintf(stdout, "|\n"); \
  } while (0); \

#define PFLBCPS(s, lvl, lit) \
     if (verbose >= lvl) { \
		PFLOG1("\t Before BCP(%d)", l2i(lit)); \
		s->printWatched(ABS(lit)); }

#define PFLBCP(s, lvl, lit) \
     if (verbose >= lvl) { \
		PFLOG1("\t BCP(%d)", l2i(lit)); \
		s->printOL(lit); \
        s->printOL(FLIP(lit)); }

#define PFLTRAIL(s, lvl) if (verbose >= lvl) s->printTrail();

#define PFLLEARNT(s, lvl) if (verbose >= lvl) s->printLearnt();

#define PFLCLAUSE(lvl, c, format, ...) \
    do { \
        if (verbose >= lvl) { \
            fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__); \
            c.print(); \
        } \
    } while (0)

#define PFLDL(s, lvl) PFLOG2(lvl, " Current decision level: %d", s->DL());

#define PFLMH(ch) \
    do { \
        fprintf(stdout, "%c ", ch);\
    } while (0); \

#define PFLMLIT(v, val) fprintf(stdout, "%c%d ", val ? ' ' : '-', v); 

#define PFLBCPE(s, lvl, lit) \
     if (verbose >= lvl) { \
		PFLOG1("\t After BCP(%d)", l2i(lit)); \
		s->printWatched(ABS(lit)); \
        PFLOGR('-', 30); }

#define PFLNEWLIT(s, lvl, src, lit) \
    do { \
        PFLOG2(lvl, " New %s( %d@%d )", src == NOREF ? s->DL() == ROOT_LEVEL ? "forced unit" : "decision"  : "unit", l2i(lit), l2dl(lit)); \
    }while(0)

#define PFLCONFLICT(s, lvl, lit) PFLOG2(lvl, " Conflict detected in literal( %d@%d )", l2i(lit), l2dl(lit));

#define PFLREDALL(s, lvl, msg) \
    if (verbose >= lvl) { \
        s->evalReds(); \
        PFLOG1("\t\t %s",msg); \
        s->logReductions();\
    }

#define PFLREDCL(s, lvl, msg) \
    if (verbose >= lvl) { \
        inf.n_del_vars_after = 0; \
        s->countCls(); \
        PFLOG1("\t\t %s", msg); \
        s->logReductions(); \
    }

#endif