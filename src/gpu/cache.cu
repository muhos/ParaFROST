/***********************************************************************[cache.cu]
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

#include <cuda_runtime.h>
#include <cassert>
#include "cache.cuh"
#include "logging.h"

using namespace pFROST;

void CACHER::destroy() {
	for (free_cache_t::iterator i = free_cache.begin(); i != free_cache.end(); i++) {
		if (cudaFree(i->second) != cudaSuccess) {
			PFLOGEN("cannot deallocate cached free memory block %p", i->second);
			throw CACHEMEMOUT();
		}
		i->second = NULL;
	}
	for (alloc_cache_t::iterator i = alloc_cache.begin(); i != alloc_cache.end(); i++) {
		if (cudaFree(i->first) != cudaSuccess) {
			PFLOGEN("cannot deallocate cached memory block %p", i->first);
			throw CACHEMEMOUT();
		}
	}
	free_cache.clear();
	alloc_cache.clear();
	used = 0;
}

void* CACHER::allocate(size_t size) {
	void* p = NULL;
	free_cache_t::iterator free_block = free_cache.lower_bound(size);
	// found free block
	if (free_block != free_cache.end()) {
		p = free_block->second;
		size = free_block->first;
		free_cache.erase(free_block);
	}
	// no free blocks, allocate new one
	else {
		if (cudaMalloc((void**)&p, size) != cudaSuccess) {
			PFLOGEN("cannot allocate new memory block via cache allocator");
			throw CACHEMEMOUT();
		}
		used += size;
	}
	assert(p);
	alloc_cache.insert(std::make_pair(p, size)); // cache new block
	return p;
}

void CACHER::deallocate(void* p) {
	alloc_cache_t::iterator allocated_block = alloc_cache.find(p);
	if (allocated_block == alloc_cache.end()) {
		PFLOGEN("memory block %p is not allocated via cache allocator", p);
		throw CACHEMEMOUT();
	}
	const size_t size = allocated_block->second;
	alloc_cache.erase(allocated_block);
	free_cache.insert(std::make_pair(size, p)); // cache free block
}

void CACHER::deallocate(void* p, const size_t bound) {
	alloc_cache_t::iterator allocated_block = alloc_cache.find(p);
	if (allocated_block == alloc_cache.end()) {
		PFLOGEN("memory block %p is not allocated via cache allocator", p);
		throw CACHEMEMOUT();
	}
	const size_t size = allocated_block->second;
	alloc_cache.erase(allocated_block);
	if (size < bound) free_cache.insert(std::make_pair(size, p)); // cache free block
	else { // deallocate free block
		if (cudaFree(p) != cudaSuccess) throw CACHEMEMOUT();
		used -= size;
	}
}