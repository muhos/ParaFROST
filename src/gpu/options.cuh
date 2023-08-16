/***********************************************************************[options.cuh]
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

#ifndef __GPU_OPTIONS_
#define __GPU_OPTIONS_

#include "datatypes.hpp"

namespace ParaFROST {

struct GOPTION {
    GOPTION();
    void init();

    bool unified_access;
    bool profile_gpu;
    bool sync_always;
    bool ve_atomic;
    bool gc_gpu;

    uint32 ve_min_threads;
    uint32 sub_min_threads;
    uint32 ere_min_threads;

    double ve_min_blocks;
    double sub_min_blocks;
    double ere_min_blocks;
};

extern GOPTION gopts;

} // namespace ParaFROST

#endif