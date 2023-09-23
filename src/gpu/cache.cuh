/***********************************************************************[cache.cuh]
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

#ifndef __GPU_CACHE_
#define __GPU_CACHE_  

#include <map>
#include "logging.hpp"
#include "constants.hpp"

namespace ParaFROST {

	// memory cache types
	typedef std::map		< void*, size_t > alloc_cache_t;
	typedef std::multimap	< size_t, void*  > free_cache_t;

	constexpr size_t MAXMEMBLOCK = 10 * MBYTE;

#define CACHEMEMCHECK(X) \
	if (X != cudaSuccess) {	\
		cudaDeviceReset(); \
		LOGERROR("cannot (de)allocate new memory block via the cache allocator"); \
	}

	// memory cache container
	class CACHER {
		size_t			used, maxcap;
		free_cache_t    free_cache;
		alloc_cache_t   alloc_cache;

	public:

		CACHER() : used(0), maxcap(0) { }
		~CACHER() { destroy(); }
		inline size_t	maxCapacity		() const { return maxcap; }
		inline void		updateMaxCap	() { if (maxcap < used) maxcap = used; }
		inline void		insert			(void* p, const size_t size) {
			assert(p);
			assert(size);
			free_cache.insert(std::make_pair(size, p));
		}
		inline void		erase			(const size_t size) {
			assert(size);
			free_cache_t::iterator free_block = free_cache.lower_bound(size);
			if (free_block != free_cache.end())
				free_cache.erase(free_block);
		}
		void			destroy			();
		void*			allocate		(size_t size);
		void			deallocate		(void* p);
		void			deallocate		(void* p, const size_t bound);

	};

	extern CACHER cacher;

}

#endif