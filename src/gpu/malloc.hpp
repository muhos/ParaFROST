/***********************************************************************[malloc.hpp]
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

#ifndef __RALLOC_
#define __RALLOC_

#include "datatypes.hpp"
#include "logging.hpp"
#include <cstdlib>
#include <cstring>

namespace ParaFROST {

	class MEMOUTEXCEPTION {};

	inline std::size_t align_up(const size_t& n, const size_t& a) {
		return (n + (a - 1)) & ~(a - 1);
	}

	template <class T>
	T* pfmalloc(const size_t& numElements) {
		if (!numElements) LOGERROR("catched zero-memory size at %s", __func__);
		T* _mem = (T*)std::malloc(numElements * sizeof(T));
		if (_mem == NULL) throw MEMOUTEXCEPTION();
		return _mem;
	}

	template <class T>
	T* pfcalloc(const size_t& numElements) {
		if (!numElements) LOGERROR("catched zero-memory size at %s", __func__);
		T* _mem = (T*)std::calloc(numElements, sizeof(T));
		if (_mem == NULL) throw MEMOUTEXCEPTION();
		return _mem;
	}

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

	template <class T>
	void pfralloc(T*& mem, const size_t& bytes) {
		if (!bytes) LOGERROR("catched zero-memory size at %s", __func__);
		T* _mem = (T*)std::realloc(mem, bytes);
		if (_mem == NULL) throw MEMOUTEXCEPTION();
		mem = _mem;
	}

	template <class T>
	void pfshrinkAlloc(T*& mem, const size_t& bytes) {
		if (!bytes) LOGERROR("catched zero-memory size at %s", __func__);
		T* _mem = NULL;
		_mem = (T*)std::realloc(_mem, bytes);
		if (_mem == NULL) throw MEMOUTEXCEPTION();
		std::memcpy(_mem, mem, bytes);
		std::free(mem);
		mem = _mem;
	}

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic pop
#endif
}

#endif