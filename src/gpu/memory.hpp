/***********************************************************************[memory.hpp]
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

#ifndef __MEMORY_
#define __MEMORY_

#include "clause.hpp"
#include "vector.hpp"
#include "logging.hpp"
#include <cstdint>
#include <limits>
#include <cassert>

namespace ParaFROST {
    /******************************************************/
    /*  Usage: global memory manager with garbage monitor */
    /*         and boundary checker                       */
    /*  Dependency:  none                                 */
    /******************************************************/
    template<class T, class S = size_t>
    class SMM
    {
        T* _mem;
        S sz, cap, maxCap;
        S _junk;
        size_t _bucket;
    protected:
        bool check(const S& d) const {
            if (d >= sz) {
                LOGERRORN("memory index (%zd) violates memory boundary (%zd)", d, sz);
                return false;
            }
            return true;
        }
        bool check(const G_REF& d) const {
            if (d < _mem) {
                LOGERRORN("memory access violation at location: %p, base address: %p", d, _mem);
                return false;
            }
            if (d > (sz - 1) + _mem) {
                LOGERRORN("memory access violation at location: %p, end address: %p", d, ((sz - 1) + _mem));
                return false;
            }
            return true;
        }
        bool checkSize(const S& newSz) const {
            if (sz != 0 && newSz <= sz) {
                LOGERRORN("size overflow during memory allocation: (new = %zd, old = %zd)", newSz, sz);
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
            if (!init_cap) return;
            assert(_bucket);
            if (init_cap > maxCap) {
                LOGERRORN("initial size exceeds maximum memory size: (max = %zd, size = %zd)", maxCap, init_cap);
                throw MEMOUTEXCEPTION();
            }
            cap = init_cap;
            pfralloc(_mem, _bucket * cap);
        }
        inline void     reserve     (const S& min_cap) {
            if (cap >= min_cap) return;
            cap = (cap > (maxCap - cap)) ? min_cap : (cap << 1);
            if (cap < min_cap) cap = min_cap;
            assert(_bucket);
            pfralloc(_mem, _bucket * cap);
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
        inline void     migrateTo   (SMM& newBlock) {
            if (newBlock._mem != NULL) std::free(newBlock._mem);
            newBlock._mem = _mem, newBlock.sz = sz, newBlock.cap = cap, newBlock._junk = _junk;
            _mem = NULL, sz = 0, cap = 0, _junk = 0;
        }
    };


    /*****************************************************/
    /*  Usage: memory manager for CNF clauses            */
    /*  Dependency:  CLAUSE, SMM                         */
    /*****************************************************/
    typedef SMM<Byte, C_REF> CTYPE;
    class CMM : public CTYPE
    {
        Vec<bool, C_REF> stencil;
    public:
        CMM() { 
            assert(CTYPE::bucket() == 1);
            assert(hc_isize == sizeof(uint32));
            assert(hc_csize == sizeof(CLAUSE)); 
        }
        explicit				CMM             (const C_REF& init_cap) : CTYPE(init_cap), stencil(init_cap, 0) { assert(CTYPE::bucket() == 1); }
        inline void				init            (const C_REF& init_cap) { CTYPE::init(init_cap), stencil.resize(init_cap, 0); }
        inline		 CLAUSE&    operator[]		(const C_REF& r) { return (CLAUSE&)CTYPE::operator[](r); }
        inline const CLAUSE&    operator[]		(const C_REF& r) const { return (CLAUSE&)CTYPE::operator[](r); }
        inline		 CLAUSE*    clause          (const C_REF& r) { return (CLAUSE*)address(r); }
        inline const CLAUSE*    clause          (const C_REF& r) const { return (CLAUSE*)address(r); }
        inline bool				deleted         (const C_REF& r) const { assert(check(r)); return stencil[r]; }
        inline void				collectClause   (const C_REF& r, const int& size) { CTYPE::collect(bytes(size)); assert(check(r)); stencil[r] = true; }
        inline void				collectLiterals (const int& size) { CTYPE::collect(size * hc_isize); }
        inline void				migrateTo       (CMM& dest) {
            CTYPE::migrateTo(dest);
            stencil.migrateTo(dest.stencil);
        }
        template <class SRC>
        inline C_REF			alloc           (const SRC& src) {
            assert(src.size() > 1);
            size_t cBytes = bytes(src.size());
            C_REF r = CTYPE::alloc(cBytes);
            new (clause(r)) CLAUSE(src);
            assert(clause(r)->capacity() == cBytes);
            assert(src.size() == clause(r)->size());
            stencil.expand(r + 1, 0);
            return r;
        }
        inline C_REF			alloc           (const int& size) {
            assert(size > 1);
            size_t cBytes = bytes(size);
            C_REF r = CTYPE::alloc(cBytes);
            new (clause(r)) CLAUSE(size);
            assert(clause(r)->capacity() == cBytes);
            assert(size == clause(r)->size());
            stencil.expand(r + 1, 0);
            return r;
        }
        inline size_t			bytes           (const int& size) {
            assert(size > 1);
            return (hc_csize + (size_t(size) - 2) * hc_isize);
        }
        inline void				destroy         () { dealloc(), stencil.clear(true); }
    };

}

#endif