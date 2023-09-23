/***********************************************************************[histogram.cuh]
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

#ifndef __GPU_HISTOGRAM_
#define __GPU_HISTOGRAM_

#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include "definitions.cuh"
#include "vector.hpp"

namespace ParaFROST {

	typedef thrust::device_ptr<uint32> t_iptr;

	struct cuHist {

		S_REF   * d_segs;
		uint32  * d_hist, * h_hist;
		uint32	* d_vorg;
		Byte	* d_lbyte;		// only used if proof is enabled
		t_iptr	thrust_hist;

		cuHist() : 
			d_segs(NULL)
			, d_hist(NULL)
			, h_hist(NULL)
			, d_vorg(NULL) 
			, d_lbyte(NULL)
		{ }

		~cuHist() 
		{
			h_hist = NULL, d_hist = NULL;
			d_segs = NULL, d_vorg = NULL;
			d_lbyte = NULL; 
		}

		inline uint32	operator[]	(const uint32& i) const {
			assert(h_hist); 
			return h_hist[i];
		}
		inline uint32&	operator[]	(const uint32& i) { 
			assert(h_hist);
			return h_hist[i];
		}
		inline void		cacheHist	(const uint32& size, const cudaStream_t& _s = 0) {
			CHECK(cudaMemcpyAsync(h_hist, d_hist, size * sizeof(uint32), cudaMemcpyDeviceToHost, _s));
		}
		inline void		fetchVars	(uVec1D& vorg, const cudaStream_t& _s = 0) {
			CHECK(cudaMemcpyAsync(d_vorg, vorg.data(), vorg.size() * sizeof(uint32), cudaMemcpyHostToDevice, _s));
		}
	};

}

#endif
