#ifndef __LOGGING_
#define __LOGGING_

#include <cstdio>
#include <cstring>

#define RULELEN 90

inline void REPCH(const char& ch, const size_t& size, const size_t& off = 0) {
    for (size_t i = off; i < size; i++) putc(ch, stdout);
}

#define PFLOG(format, ...) \
  do { \
     fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__), putc('\n', stdout); \
  } while (0)

#define PFLOGN(format, ...) \
  do { \
     fprintf(stdout, "c |"), fprintf(stdout, format, ## __VA_ARGS__); \
  } while (0)

#define PFLOGS(sat) fprintf(stdout, "s %s\n", sat);

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
     const char* title = "SAT Solver"; \
     size_t len = strlen(title) + strlen(name); \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the title"); \
     size_t gap = (RULELEN - len) / 2; \
     PFLOGN(" "); \
     REPCH(' ', gap); \
     fprintf(stdout, "%s %s", name, title); \
     REPCH(' ', RULELEN + 1, len + gap + 3), fprintf(stdout, "|\n"); \
     const char* rights = "Technische Universiteit Eindhoven (TU/e), all rights reserved."; \
     len = strlen(rights) - 1; \
     if (RULELEN < len) PFLOGE("ruler length is smaller than the rights"); \
     gap = RULELEN - strlen(rights) - 1; \
     PFLOGN(" %s", rights); \
     REPCH(' ', gap), fprintf(stdout, "|\n"); \
  } while (0); \
  PFLOGR('-', RULELEN);

#endif