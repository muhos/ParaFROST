/***********************************************************************[segsort.cu]
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

#if defined(_WIN32) && !defined(__CUDACC__)
#define __CUDACC__
#endif
#include <moderngpu/kernel_segsort.hxx>
#include "modernalloc.cuh"
#include "options.cuh"
#include "timer.cuh"

using namespace ParaFROST;

__global__ void print_segments(const S_REF* segs, const uint32* hist, const uint32 n) {
	for (int i = 0; i < n; i++) {
		printf("c   Segment %d starts at %lld with count %d\n", i, segs[i], hist[i]);
	}
}

void Solver::sortOT() {
	assert(cumm.occurs());
	if (gopts.profile_gpu) cutimer.start();
	static MCA mca(cacher, 0, 0);
	const int offset = 3; // first three elements in occurs = zero
	segmented_sort(cumm.occurs(), inf.numLiterals, cuhist.d_segs + offset, inf.maxDualVars - offset, olist_cmp, mca);
	if (gopts.sync_always) {
		LASTERR("Sorting OT failed");
		SYNCALL;
	}
	if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.sot += cutimer.gpuTime();
}