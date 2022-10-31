/***********************************************************************[statistics.h]
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

#ifndef __STATS_
#define __STATS_

#include "datatypes.h"
#include <cassert>
#include <cstring>

namespace ParaFROST {

	struct SIGMASTATS {
		uint32 calls;
		struct { uint64 removed; } ere;
		struct { uint64 variables, clauses; int64 literals; } all;
	};

	struct STATS {
		uint64 sysmem;
		uint64 reuses;
		uint64 reduces;
		uint64 subtried;
		int64 shrunken;
		struct { uint64 ticks; } mode;
		struct { uint64 hard, soft; } recycle;
		struct { uint64 all, stable; } restart;
		struct { uint64 before, after; } minimize;
		struct { uint32 learnt, forced; } units;
		struct { uint64 original, learnt; } clauses, literals;
		struct { uint64 chrono, nonchrono; } backtrack;
		struct { uint64 calls, eliminated; } autarky;
		struct { uint64 probed, failed, removed; } transitive;
		struct { uint64 calls, clauses, literals; } shrink;
		struct { uint32 calls, vmtf, vsids, walks; } mdm;
		struct { uint64 single, multiple, massumed; } decisions;
		struct { uint64 calls, binaries, hyperunary; } debinary;
		struct { uint64 calls, rounds, failed, probed; } probe;
		struct { uint64 resolutions, resolvents, reduced; } binary;
		struct { uint64 all, random, best, inv, org, flip; } rephase;
		struct { uint64 calls, checks, minimum, flipped, improved; } walk;
		struct { uint64 calls, scc, variables, hyperunary, clauses; } decompose;
		struct {
			uint64 checks, calls;
			uint64 resolutions, binaries, ternaries, reduced;
		} ternary;
		struct {
			uint64 checks, leftovers, calls;
			uint64 subsumed, strengthened, learntfly;
			uint64 subsumedfly, strengthenedfly;
		} subsume;
		struct {
			uint64 assumed, reused;
			uint64 checks, vivified, implied, subsumed, strengthened;
		} vivify;
		SIGMASTATS sigma;
		uint64 searchprops, searchticks, probeticks, transitiveticks;
		uint64 marker;
		uint64 conflicts;
		uint64 stablemodes;
		uint64 unstablemodes;
		uint32 mappings;

		STATS() { memset(this, 0, sizeof(*this)); }
	};
}

#endif