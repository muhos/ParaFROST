/***********************************************************************[segsort.cu]
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

#if defined(_WIN32) && !defined(__CUDACC__)
#define __CUDACC__
#endif
#include <moderngpu/kernel_segsort.hxx>
#include "modernalloc.cuh"

using namespace mgpu;
using namespace pFROST;

MCA context(0, 0); // for segmented sort

void ParaFROST::sortOT() {
	assert(cumm.occurs());
	if (gopts.profile_gpu) cutimer->start();
	const int offset = 3; // first three elements in occurs = zero
	segmented_sort(cumm.occurs(), inf.nLiterals, cuhist.d_segs + offset, inf.nDualVars - offset, olist_cmp, context);
	if (gopts.sync_always) {
		LOGERR("Sorting OT failed");
		syncAll();
	}
	if (gopts.profile_gpu) cutimer->stop(), cutimer->sot += cutimer->gpuTime();
}