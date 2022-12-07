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

namespace ParaFROST {

	/*****************************************************/
	/*  Usage:    Modern-GPU cached memory allocator     */
	/*  Dependency: context_t                            */
	/*****************************************************/

	class MCA : public mgpu::context_t {
	protected:
		cudaDeviceProp _props;
		int _ptx_version;
		cudaStream_t _stream;

	public:

		MCA(cudaStream_t stream_ = 0) :
			context_t()
			, _ptx_version(0)
			, _stream(stream_) { }

		void init();
		virtual const cudaDeviceProp& props		() const { return _props; }
		virtual int				ptx_version		() const { return _ptx_version; }
		virtual cudaStream_t	stream			() { return _stream; }
		virtual void*			alloc			(size_t size, mgpu::memory_space_t space) {
			return cacher.allocate(size);
		}
		virtual void			free			(void* p, mgpu::memory_space_t space) {
			cacher.deallocate(p, MAXMEMBLOCK);
		}
		virtual void			synchronize		() {
			cudaError_t result = _stream ?
				cudaStreamSynchronize(_stream) :
				cudaDeviceSynchronize();
			if (cudaSuccess != result) throw mgpu::cuda_exception_t(result);
		}
		virtual cudaEvent_t		event			() { return NULL; }
		virtual void			timer_begin		() {}
		virtual double			timer_end		() { return 0; }

	};

}

#endif