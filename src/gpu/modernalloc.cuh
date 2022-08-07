/***********************************************************************[modernalloc.cuh]
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

#ifndef __MODERNGPU_
#define __MODERNGPU_  

#include <exception>
#include <cuda_runtime.h>
#include "moderngpu/context.hxx"
#include "logging.h"
#include "cache.cuh"

constexpr size_t MAXMEMBLOCK = 10 * MBYTE;

using namespace mgpu;

namespace pFROST {

	/*****************************************************/
	/*  Usage:    Modern-GPU cached memory allocator     */
	/*  Dependency: context_t                            */
	/*****************************************************/

	class MCA : public context_t {
	protected:
		cudaDeviceProp _props;
		int _ptx_version;
		cudaStream_t _stream;

		template<int dummy_arg = 0>
		void init() {
			cudaFuncAttributes attr;
			cudaError_t result = cudaFuncGetAttributes(&attr, dummy_k<0>);
			if (cudaSuccess != result) throw cuda_exception_t(result);
			_ptx_version = attr.ptxVersion;

			int ord;
			cudaGetDevice(&ord);
			cudaGetDeviceProperties(&_props, ord);
		}

	public:
		MCA(bool print_prop = true, cudaStream_t stream_ = 0) :
			context_t(), _stream(stream_) {

			init();
			if (print_prop) {
				printf("%s\n", device_prop_string(_props).c_str());
			}
		}

		virtual const cudaDeviceProp& props	() const { return _props; }
		virtual int				ptx_version	() const { return _ptx_version; }
		virtual cudaStream_t	stream		() { return _stream; }
		virtual void*			alloc		(size_t size, memory_space_t space) {
			return cacher.allocate(size);
		}
		virtual void			free		(void* p, memory_space_t space) {
			cacher.deallocate(p, MAXMEMBLOCK);
		}
		virtual void		synchronize		() {
			cudaError_t result = _stream ?
				cudaStreamSynchronize(_stream) :
				cudaDeviceSynchronize();
			if (cudaSuccess != result) throw cuda_exception_t(result);
		}
		virtual cudaEvent_t event			() { return NULL; }
		virtual void		timer_begin		() {}
		virtual double		timer_end		() { return 0; }

	};

}

#endif