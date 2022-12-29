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

#include "solve.hpp"

#if defined(_WIN32) && !defined(__CUDACC__)
#define __CUDACC__
#endif
#include <moderngpu/kernel_segsort.hxx>
#include "modernalloc.cuh"

using namespace mgpu;
using namespace ParaFROST;

__global__ void ptx_dummy_k() { }

MCA context(0);

void MCA::init() {
	cudaFuncAttributes attr;
	cudaError_t result = cudaFuncGetAttributes(&attr, ptx_dummy_k);
	if (cudaSuccess != result) throw mgpu::cuda_exception_t(result);
	_ptx_version = attr.ptxVersion;
	_props = devProp;
}

void Solver::segsortOTAsync() {
	assert(cumm.occurs());
	if (gopts.profile_gpu) cutimer->start();
	if (!context.ptx_version()) context.init();
	const int offset = 3; // first three elements in occurs = zero
	segmented_sort(cumm.occurs(), inf.nLiterals, cuhist.d_segs + offset, inf.nDualVars - offset, olist_cmp, context);
	if (gopts.sync_always) {
		LOGERR("Sorting OT failed");
		syncAll();
	}
	if (gopts.profile_gpu) cutimer->stop(), cutimer->sot += cutimer->gpuTime();
}