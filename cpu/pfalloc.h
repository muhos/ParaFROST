/***********************************************************************[pfalloc.h]
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

#ifndef __ALLOC_
#define __ALLOC_

#include "pfralloc.h"
#include "pfdtypes.h"
#include "pflogging.h"
#include <cstdint>
#include <limits>
#include <cassert>

namespace pFROST {
    /******************************************************/
    /*  Usage: global memory manager with garbage monitor */
    /*  Dependency:  none                                 */
    /******************************************************/
    template<class T, class S = uint32>
    class SMM
    {
        T* _mem;
        S sz, cap, maxCap;
        S _junk;
        size_t _bucket;
        bool check(const S& d) const {
            if (d >= sz) {
                PFLOGEN("memory index (%d) violates memory boundary (%d)", d, sz);
                return false;
            }
            return true;
        }
        bool check(const G_REF& d) const {
            if (d < _mem) {
                PFLOGEN("memory access violation at location: %p, base address: %p", d, _mem);
                return false;
            }
            if (d > (sz - 1) + _mem) {
                PFLOGEN("memory access violation at location: %p, end address: %p", d, ((sz - 1) + _mem));
                return false;
            }
            return true;
        }
        bool checkSize(const S& newSz) const {
            if (sz != 0 && newSz <= sz) {
                PFLOGEN("size overflow during memory allocation: (new = %d, old = %d)", newSz, sz);
                return false;
            }
            return true;
        }
    public:
                        ~SMM        () { dealloc(); }
                        SMM         () {
            maxCap = std::numeric_limits<S>::max();
            assert(maxCap > INT8_MAX);
            _mem = NULL, _bucket = sizeof(T), sz = 0LL, cap = 0LL, _junk = 0LL;
        }
        explicit        SMM         (const S& _cap) {
            maxCap = std::numeric_limits<S>::max();
            assert(maxCap > INT8_MAX);
            _mem = NULL, _bucket = sizeof(T), sz = 0LL, cap = 0LL, _junk = 0LL, init(_cap);
        }
        inline void     dealloc     () { if (_mem != NULL) std::free(_mem), _mem = NULL; sz = cap = 0, _junk = 0; }
        inline size_t   bucket      () const { assert(_bucket); return _bucket; }
        inline S        size        () const { return sz; }
        inline S        garbage     () const { return _junk; }
        inline T&       operator[]  (const S& idx) { assert(check(idx)); return _mem[idx]; }
        inline const T& operator[]  (const S& idx) const { assert(check(idx)); return _mem[idx]; }
        inline T*       address     (const S& idx) { assert(check(idx)); return _mem + idx; }
        inline const T* address     (const S& idx) const { assert(check(idx)); return _mem + idx; }
        inline void     collect     (const S& size) { _junk += size; }
        inline void     init        (const S& init_cap) {
            if (init_cap == 0) return;
            assert(_bucket);
            if (init_cap > maxCap) {
                printf("Error - initial size exceeds maximum memory size: (max = %d, size = %d)\n", maxCap, init_cap);
                throw MEMOUTEXCEPTION();
            }
            cap = init_cap;
            pfalloc(_mem, _bucket * cap);
        }
        inline void     reserve     (const S& min_cap) {
            if (cap >= min_cap) return;
            if (cap > (maxCap - cap)) cap = min_cap;
            else { cap <<= 1; if (cap < min_cap) cap = min_cap; }
            assert(_bucket);
            pfalloc(_mem, _bucket * cap);
        }
        inline S        alloc       (const S& size) {
            assert(size > 0);
            assert(checkSize(sz + size));
            reserve(sz + size);
            S oldSz = sz;
            sz += size;
            assert(sz > 0);
            return oldSz;
        }
        inline void     migrate     (SMM& newBlock) {
            if (newBlock._mem != NULL) std::free(newBlock._mem);
            newBlock._mem = _mem, newBlock.sz = sz, newBlock.cap = cap, newBlock._junk = _junk;
            _mem = NULL, sz = 0, cap = 0, _junk = 0;
        }
    };

}

#endif