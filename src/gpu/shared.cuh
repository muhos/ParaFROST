/***********************************************************************[shared.cuh]
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

#ifndef __GPU_SHARED_
#define __GPU_SHARED_

#include "definitions.cuh"

namespace ParaFROST {

	template<class T>
	class SharedMemory
	{
	public:
		_PFROST_D_ operator T* () {
			extern __shared__ int _smem[];
			return (T*)_smem;
		}
		_PFROST_D_ operator const T* () const {
			extern __shared__ int _smem[];
			return (T*)_smem;
		}
	};

}

#endif
