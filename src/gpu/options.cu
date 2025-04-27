/***********************************************************************[options.cu]
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

#include "options.cuh"
#include "input.hpp"

namespace ParaFROST {
	__constant__ KOptions kOpts[1];

	void initDevOpts()
	{
		CHECK(cudaMemcpyToSymbol(kOpts, &gopts.hostKOpts, sizeof(KOptions), 0, cudaMemcpyHostToDevice));
	}
}

using namespace ParaFROST;

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

// kernel options
BOOL_OPT opt_ve_fun_en("vefunction", "enable function table reasoning", true);
BOOL_OPT opt_ve_lbound_en("velitsbound", "skip variables resulting in more literals than original", false);

INT_OPT opt_bce_max_occurs("bcemaxoccurs", "maximum occurrence list size to scan in BCE", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_sub_max_occurs("submaxoccurs", "maximum occurrence list size to scan in SUB", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_sub_clause_max("subclausemax", "maximum clause size to check in SUB", 100, INT32R(0, INT32_MAX));
INT_OPT opt_ere_extend("ereextend", "extend ERE with clause strengthening (0: no extend, 1: originals, 2: all)", 1, INT32R(0, 3));
INT_OPT opt_ere_max_occurs("eremaxoccurs", "maximum occurrence list size to scan in ERE", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_ere_clause_max("ereclausemax", "maximum resolvent size for forward check in ERE", 250, INT32R(2, INT32_MAX));
INT_OPT opt_ve_clause_max("resolventmax", "maximum resolvent size in BVE (0: no limit)", 100, INT32R(0, INT32_MAX));
INT_OPT opt_xor_max_arity("xormaxarity", "maximum XOR fanin size", 10, INT32R(2, 20));

GOPTION::GOPTION() {
	RESETSTRUCT(this);
}

void GOPTION::init(const bool& opt_proof_en) 
{
	ve_min_threads		= opt_ve_min_threads;
	ve_min_blocks		= opt_ve_min_blocks;
	sub_min_threads		= opt_sub_min_threads;
	sub_min_blocks		= opt_sub_min_blocks;
	ere_min_threads		= opt_ere_min_threads;
	ere_min_blocks		= opt_ere_min_blocks;
	sync_always			= opt_sync_always_en;
	profile_gpu			= opt_profile_gpu_en;
	unified_access		= opt_unified_access_en;

	hostKOpts.proof_en			= opt_proof_en;
	hostKOpts.ve_fun_en			= opt_ve_fun_en;
	hostKOpts.ve_lbound_en		= opt_ve_lbound_en;
	hostKOpts.ve_clause_max		= opt_ve_clause_max;
	hostKOpts.xor_max_arity		= opt_xor_max_arity;
	hostKOpts.ere_extend_en		= opt_ere_extend > 0;
	hostKOpts.ere_extend_all_en = opt_ere_extend == 2;
	hostKOpts.ere_clause_max	= MIN(opt_ere_clause_max, SH_MAX_ERE_OUT);
	hostKOpts.ere_max_occurs	= opt_ere_max_occurs;
	hostKOpts.sub_max_occurs	= opt_sub_max_occurs;
	hostKOpts.sub_clause_max	= opt_sub_clause_max;
	hostKOpts.bce_max_occurs	= opt_bce_max_occurs;
}