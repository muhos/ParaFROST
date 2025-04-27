/***********************************************************************[vector.cu]
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

#include <cuda_runtime.h>
#include "vector.cuh"

namespace ParaFROST {

	template<>
	_PFROST_D_ void cuVec<S_REF>::insert(const S_REF& val) {
		const uint32 idx = atomicInc(&sz, cap);
		assert(checkAtomicBound(idx, cap));
		_mem[idx] = val;
	}

	
	template<>
	_PFROST_D_ void cuVec<uint32>::insert(const uint32& val) {
		const uint32 idx = atomicInc(&sz, cap);
		assert(checkAtomicBound(idx, cap));
		_mem[idx] = val;
	}

	template<>
	_PFROST_D_ S_REF* cuVec<S_REF>::jump(const uint32& n)
	{
		const uint32 idx = atomicAdd(&sz, n);
		assert(checkAtomicBound(idx, cap));
		return _mem + idx;
	}

	template<>
	_PFROST_D_ uint32* cuVec<uint32>::jump(const uint32& n)
	{
		const uint32 idx = atomicAdd(&sz, n);
		assert(checkAtomicBound(idx, cap));
		return _mem + idx;
	}

	template<>
	_PFROST_D_ Byte* cuVec<Byte>::jump(const uint32& n)
	{
		const uint32 idx = atomicAdd(&sz, n);
		assert(checkAtomicBound(idx, cap));
		return _mem + idx;
	}

}