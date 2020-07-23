/***********************************************************************[pfsimpopts.h]
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

#ifndef __SIGMA_OPTS_
#define __SIGMA_OPTS_

#include "pfargs.h"

namespace pFROST {

#define MAX_GPU_COUNT (31 + 1) // index starts from GPU 1 (first slave)
#define MAX_STREAMS 32

	BOOL_OPT opt_ve_en("ve", "enable bounded variable elimination (BVE)", true);
	BOOL_OPT opt_ve_plus_en("ve+", "enable HSE before BVE", true);
	BOOL_OPT opt_sub_en("sub", "enable hybrid subsumption elimination (HSE) with high bounds in LCVE", false);
	BOOL_OPT opt_bce_en("bce", "enable blocked clause elimination", false);
	BOOL_OPT opt_hre_en("hre", "enable hidden redundancy elimination", false);
	BOOL_OPT opt_all_en("all", "enable all simplifications", false);
	BOOL_OPT opt_unified_access_en("unified-access", "enable unified access in transferring CNF to/from device", false);
	BOOL_OPT opt_profile_gpu_en("profile-gpu", "profile device kernels", false);
	BOOL_OPT opt_sync_always_en("sync-always", "synchronize all device calls with the host (useful for debugging)", false);
	BOOL_OPT opt_atomic_ve_en("atomic-ve", "enable atomic version of BVE (produces out-of-order resolvents)", false);
	BOOL_OPT opt_sort_cnf_en("sort-simp", "sort simplified formula before writing to host", true);
	BOOL_OPT opt_gc_par_en("gc-par", "enable device garbage collection in parallel (stream compaction approach)", true);
	BOOL_OPT opt_solve_en("solve", "proceed with solving after simplifications", true);
	INT_OPT opt_lcve_min("lcve-min", "minimum parallel variables to simplify", 2, INT32R(1, INT32_MAX));
	INT_OPT opt_ve_phase_min("ve-phase-min", "minimum removed literals to stop stage-1 reductions ", 500, INT32R(1, INT32_MAX));
	INT_OPT opt_mu_pos("mu-pos", "set the positive freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
	INT_OPT opt_mu_neg("mu-neg", "set the negative freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
	INT_OPT opt_xor_max_arity("xor-max-arity", "maximum XOR fanin size", 64, INT32R(2, INT32_MAX));
	INT_OPT opt_hse_max_occurs("hse-max", "maximum occurrence list size to scan in HSE", 30000, INT32R(100, INT32_MAX));
	INT_OPT opt_bce_max_occurs("bce-max", "maximum occurrence list size to scan in BCE", 30000, INT32R(100, INT32_MAX));
	INT_OPT opt_hre_max_occurs("hre-max", "maximum occurrence list size to scan in HRE", 30000, INT32R(100, INT32_MAX));
	INT_OPT opt_phases("phases", "set the number of phases in stage-1 reductions", 5, INT32R(0, INT32_MAX));
	INT_OPT opt_cnf_free("gc-freq", "set the frequency of device garbage collection", 2, INT32R(0, 3));
	INT_OPT opt_gpus("ngpus", "number of GPU(s) to be activated", 1, INT32R(1, MAX_GPU_COUNT));
	INT_OPT opt_streams("nstreams", "number of GPU streams to be created", 6, INT32R(1, MAX_STREAMS));

}

#endif