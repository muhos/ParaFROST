/***********************************************************************[random.hpp]
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

#ifndef __RANDOM_
#define __RANDOM_

#include <cassert>
#include "datatypes.hpp"

namespace ParaFROST {

    class RANDOM {
        uint32 _seed_;

    public:
                      RANDOM          () : _seed_(1) {}
                      RANDOM          (const uint32& seed) : _seed_(seed) { assert(_seed_); }
        inline void   init            (const uint32& seed) { _seed_ = seed; assert(_seed_); }
        inline uint32 seed            () const { return _seed_; }
        inline uint32 irand           () {
            _seed_ ^= _seed_ << 13;
            _seed_ ^= _seed_ >> 17;
            _seed_ ^= _seed_ << 5;
            return _seed_;
        }
        inline double drand           () {
            return irand() * 2.328306e-10;
        }
        inline bool   brand           () {
            const double fraction = drand();
            assert(fraction >= 0 && fraction < 1);
            return uint32(2 * fraction);
        }
    };

}

#endif
