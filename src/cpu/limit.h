/***********************************************************************[limit.h]
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

#ifndef __LIMIT_
#define __LIMIT_

#include "datatypes.h"
#include "scale.h"

namespace ParaFROST {

struct LAST {
  struct {
    uint64 conflicts, count;
    uint32 best, target;
    LIT_ST type;
  } rephase;
  struct {
    uint32 decisions, unassigned;
    int rounds;
  } mdm;
  struct {
    uint64 reduces;
  } sigma, probe;
  struct {
    uint64 resolvents;
  } ternary;
  struct {
    uint32 literals;
  } transitive;
  struct {
    double inc, booster;
    inline void boost() {
      assert(booster >= 1);
      inc *= booster;
    }
    inline void scale(const double& val) {
      assert(val > 0 && val < 1);
      inc *= val;
      assert(inc <= (1.0 / val));
    }
  } vsids;
  struct {
    int64 removed;
  } shrink;

  LAST() { RESETSTRUCT(this); }
};
struct LIMIT {
  uint64 mdm;
  uint64 sigma;
  uint64 probe;
  uint64 reduce;
  uint64 rephase;
  uint64 subsume;
  struct {
    uint64 ticks, conflicts;
  } mode;
  struct {
    uint64 conflicts;
  } restart;
  int keptsize, keptlbd;

  LIMIT() { RESETSTRUCT(this); }
};

struct MONITOR {
  uint32 now, all;
};

struct SLEEP {
  MONITOR sigma, probe, ternary, autarky, debinary;
  SLEEP() { RESETSTRUCT(this); }
};

#define INIT_LIMIT(RESULT, OPTION_INC, SCALE_INCREASE) \
  RESULT = (SCALE_INCREASE) ? relscale(stats.clauses.original, OPTION_INC) : OPTION_INC;

#define INCREASE_LIMIT(OPTION, N, SCALING_FUNC, SCALE_INCREASE)                                        \
  do {                                                                                                 \
    uint64 INC = opts.OPTION##_inc;                                                                    \
    INC *= SCALING_FUNC(N) + 1;                                                                        \
    const uint64 SCALED = (SCALE_INCREASE) ? relscale(stats.clauses.original, INC) : INC;              \
    limit.OPTION = stats.conflicts + SCALED;                                                           \
    PFLOG2(2, "  %s limit increased to %lld conflicts by a weight %lld", __func__, limit.OPTION, INC); \
  } while (0)

#define SET_BOUNDS(RESULT, OPTION, START, REFERENCE, SCALE)                                             \
  uint64 RESULT = stats.START;                                                                          \
  do {                                                                                                  \
    const uint64 STATREF = stats.REFERENCE;                                                             \
    const uint64 MINIMUM = opts.OPTION##_min_eff;                                                       \
    const uint64 MAXIMUM = MINIMUM * (opts.OPTION##_max_eff);                                           \
    const double RELEFF = (double)(opts.OPTION##_rel_eff) * 1e-3;                                       \
    uint64 INCREASE = uint64(STATREF * RELEFF) + SCALE;                                                 \
    if (INCREASE < MINIMUM) INCREASE = MINIMUM;                                                         \
    if (INCREASE > MAXIMUM) INCREASE = MAXIMUM;                                                         \
    RESULT += INCREASE;                                                                                 \
    PFLOG2(2, "  %s efficiency bounds increased to %lld by a weight %lld", __func__, RESULT, INCREASE); \
  } while (0)

#define SLEEPING(SMONITOR, ENABLED)                                         \
  do {                                                                      \
    if (!ENABLED) break;                                                    \
    MONITOR& monitor = SMONITOR;                                            \
    assert(monitor.all <= monitor.now);                                     \
    if (!monitor.all) break;                                                \
    PFLOG2(2, "  %s still sleeping for %d time(s)", __func__, monitor.all); \
    monitor.all--;                                                          \
    return;                                                                 \
  } while (0)

#define UPDATE_SLEEPER(SMONITOR, SUCCESS)                                            \
  do {                                                                               \
    if (!opts.SMONITOR##_sleep_en) break;                                            \
    MONITOR& monitor = sleep.SMONITOR;                                               \
    const uint32 PERIOD = opts.nap;                                                  \
    assert(monitor.all <= monitor.now);                                              \
    if (SUCCESS) {                                                                   \
      if (monitor.now) {                                                             \
        PFLOG2(2, "  %s is waking up after %d passed", __func__, monitor.now);       \
        monitor.now = monitor.all = 0;                                               \
      } else                                                                         \
        assert(!monitor.all);                                                        \
    } else {                                                                         \
      if (monitor.now < PERIOD) {                                                    \
        monitor.now++;                                                               \
        PFLOG2(2, "  %s sleeping period is increased to %d", __func__, monitor.now); \
      } else                                                                         \
        monitor.all = monitor.now;                                                   \
    }                                                                                \
    assert(monitor.all <= monitor.now);                                              \
  } while (0)
} // namespace ParaFROST

#endif