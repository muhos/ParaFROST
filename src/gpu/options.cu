/***********************************************************************[options.cu]
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

#include "input.hpp"
#include "options.cuh"

using namespace ParaFROST;

BOOL_OPT opt_ve_atomic_en("veatomic", "enable atomic version of BVE (produces out-of-order resolvents)", false);
BOOL_OPT opt_gc_gpu_en("gcgpu", "enable device garbage collection in parallel (stream compaction approach)", true);
BOOL_OPT opt_sync_always_en("syncalways", "synchronize all device calls with the host (useful for debugging)", false);
BOOL_OPT opt_profile_gpu_en("profilegpu", "profile simplifications", false);
BOOL_OPT opt_unified_access_en("unifiedaccess", "enable unified access in transferring data to/from device", false);

INT_OPT opt_ve_min_threads("veminthreads",
                           "set the minimum x-threads allowed per block in BVE (the tuner may optimize it to a higher value)", 4, INT32R(2, 1024));
INT_OPT opt_sub_min_threads("subminthreads",
                            "set the minimum x-threads allowed per block in SUB (the tuner may optimize it to a higher value)", 4, INT32R(2, 1024));
INT_OPT opt_ere_min_threads("ereminthreads",
                            "set the minimum y-threads allowed per block in ERE (the tuner may optimize it to a higher value)", 4, INT32R(2, 32));
DOUBLE_OPT opt_ve_min_blocks("veminblocks",
                             "set the minimum percentage of supported blocks to trigger the embedded tuner in BVE", 0.5, FP64R(0, 1));
DOUBLE_OPT opt_sub_min_blocks("subminblocks",
                              "set the minimum percentage of supported blocks to trigger the embedded tuner in SUB", 0.5, FP64R(0, 1));
DOUBLE_OPT opt_ere_min_blocks("ereminblocks",
                              "set the minimum percentage of supported blocks to trigger the embedded tuner in ERE", 0.5, FP64R(0, 1));

GOPTION::GOPTION() {
    RESETSTRUCT(this);
}

void GOPTION::init() {
    gc_gpu = opt_gc_gpu_en;
    ve_atomic = opt_ve_atomic_en;
    ve_min_threads = opt_ve_min_threads;
    ve_min_blocks = opt_ve_min_blocks;
    sub_min_threads = opt_sub_min_threads;
    sub_min_blocks = opt_sub_min_blocks;
    ere_min_threads = opt_ere_min_threads;
    ere_min_blocks = opt_ere_min_blocks;
    sync_always = opt_sync_always_en;
    profile_gpu = opt_profile_gpu_en;
    unified_access = opt_unified_access_en;
}