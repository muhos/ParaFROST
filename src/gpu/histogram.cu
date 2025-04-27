/***********************************************************************[histogram.cu]
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
#include "options.cuh"
#include "timer.cuh"
#include "count.cuh"
#include "grid.cuh"
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/device_ptr.h>


namespace ParaFROST {

	__constant__ DCPTR DC_PTRS[1];

	__global__
		void printCMem()
	{
		printf("c   Variable mapping array %p\n", DC_PTRS->d_vorg);
		printf("c   Literal  bytes array %p\n", DC_PTRS->d_lbyte);
	}

	void printConstants()
	{
		printCMem << <1, 1 >> > ();
		LASTERR("Printing constant memory failed");
		SYNCALL;
	}

	void initDevVorg(const cuHist& cuhist)
	{
		DCPTR ptrs = { cuhist.d_vorg, cuhist.d_lbyte };
		CHECK(cudaMemcpyToSymbol(DC_PTRS, &ptrs, sizeof(DCPTR), 0, cudaMemcpyHostToDevice));
	}

	void Solver::histSimp(const uint32& numLits)
	{
		LOGN2(2, " Computing histogram on %d elements..", numLits);
		assert(numLits);
		cuLits& culits = cumm.literals();
		assert(culits.size >= numLits);
		thrust::device_ptr<uint32> thrust_lits = thrust::device_ptr<uint32>(culits.mem);
		thrust::device_ptr<uint32> thrust_hist = thrust::device_ptr<uint32>(cuhist.d_hist);
		SYNC(0); // sync 'flattenCNF'
		if (gopts.profile_gpu) cutimer->start();
		cacher.insert(cumm.scatter(), cumm.scatterCap());
		thrust::sort(thrust::cuda::par(tca), thrust_lits, thrust_lits + numLits);
		thrust::counting_iterator<size_t> search_begin(0);
		thrust::upper_bound(thrust::cuda::par(tca), thrust_lits, thrust_lits + numLits, search_begin, search_begin + inf.nDualVars, thrust_hist);
		thrust::adjacent_difference(thrust::cuda::par(tca), thrust_hist, thrust_hist + inf.nDualVars, thrust_hist);
		cacher.erase(cumm.scatterCap());
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LOGDONE(2, 5);
	}

}