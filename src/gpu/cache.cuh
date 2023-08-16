/***********************************************************************[cache.cuh]
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

#ifndef __CU_CACHE_
#define __CU_CACHE_

#include <map>
#include "logging.hpp"

namespace ParaFROST {

// memory cache types
typedef std::map<void*, size_t> alloc_cache_t;
typedef std::multimap<size_t, void*> free_cache_t;

// memory exception types
class CACHEMEMOUT {};

// memory cache container
class CACHER {
    size_t used, maxcap;
    free_cache_t free_cache;
    alloc_cache_t alloc_cache;

  public:
    CACHER() : used(0), maxcap(0) {}

    inline size_t maxCapacity() const { return maxcap; }
    inline void updateMaxCap() {
        if (maxcap < used) maxcap = used;
    }
    inline void insert(void* p, const size_t size) {
        assert(p);
        assert(size);
        free_cache.insert(std::make_pair(size, p));
    }
    inline void erase(const size_t size) {
        assert(size);
        free_cache_t::iterator free_block = free_cache.lower_bound(size);
        if (free_block != free_cache.end())
            free_cache.erase(free_block);
        else {
            PFLOGEN("could not find the cached free memory block to erase");
            throw CACHEMEMOUT();
        }
    }
    void destroy();
    void* allocate(size_t size);
    void deallocate(void* p);
    void deallocate(void* p, const size_t bound);
};

extern CACHER cacher;

} // namespace ParaFROST

#endif