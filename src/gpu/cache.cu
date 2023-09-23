/***********************************************************************[cache.cu]
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
#include <cassert>
#include "cache.cuh"
#include "logging.hpp"

using namespace ParaFROST;

void CACHER::destroy() {
	for (free_cache_t::iterator i = free_cache.begin(); i != free_cache.end(); i++) {
		CACHEMEMCHECK(cudaFree(i->second));
		i->second = NULL;
	}
	for (alloc_cache_t::iterator i = alloc_cache.begin(); i != alloc_cache.end(); i++) {
		CACHEMEMCHECK(cudaFree(i->first));
	}
	free_cache.clear();
	alloc_cache.clear();
	used = 0;
}

void* CACHER::allocate(size_t size) {
	if (!size)
		return NULL;
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
		CACHEMEMCHECK(cudaMalloc((void**)&p, size));
		used += size;
	}
	assert(p);
	alloc_cache.insert(std::make_pair(p, size)); // cache new block
	return p;
}

void CACHER::deallocate(void* p) {
	alloc_cache_t::iterator allocated_block = alloc_cache.find(p);
	if (allocated_block == alloc_cache.end()) {
		LOGERROR("memory block %p is not allocated via cache allocator", p);
	}
	const size_t size = allocated_block->second;
	alloc_cache.erase(allocated_block);
	free_cache.insert(std::make_pair(size, p)); // cache free block
}

void CACHER::deallocate(void* p, const size_t bound) {
	alloc_cache_t::iterator allocated_block = alloc_cache.find(p);
	if (allocated_block == alloc_cache.end()) {
		LOGERROR("memory block %p is not allocated via cache allocator", p);
	}
	const size_t size = allocated_block->second;
	alloc_cache.erase(allocated_block);
	if (size < bound) free_cache.insert(std::make_pair(size, p)); // cache free block
	else { // deallocate free block
		CACHEMEMCHECK(cudaFree(p));
		used -= size;
	}
}