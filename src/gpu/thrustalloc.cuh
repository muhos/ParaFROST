/***********************************************************************[thrustalloc.cuh]
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

#ifndef __THRUST_MEMORY_
#define __THRUST_MEMORY_

#include <cuda_runtime.h>
#include <thrust/system/cuda/vector.h>
#include <thrust/device_ptr.h>
#include "logging.hpp"
#include "cache.cuh"

namespace ParaFROST {

	/*****************************************************/
	/*  Usage:    Thrust cached memory allocator         */
	/*  Dependency: None                                 */
	/*****************************************************/

	class TCA {

	public:
		typedef char value_type;

		char* allocate(size_t size) {
			return (char*)cacher.allocate(size);
		}

		void deallocate(char* p, size_t) {
			cacher.deallocate(p);
		}
	};

}

#endif
