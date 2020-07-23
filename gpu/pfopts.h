/***********************************************************************[pfopts.h]
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

#ifndef __SOLVE_OPTS_
#define __SOLVE_OPTS_

#include "pfargs.h"

namespace pFROST {

	BOOL_OPT opt_sig_pre_en("sigma", "enable preprocessing using SIGmA", false);
	BOOL_OPT opt_sig_live_en("sigma-live", "enable live SIGmA (inprocessing)", false);
	BOOL_OPT opt_report_en("report", "allow performance report on stdout", true);
	BOOL_OPT opt_parseonly_en("parseonly", "parse only the input formula", false);
	BOOL_OPT opt_lcv_en("lcv", "use least-constrained variables to make multiple decisions", false);
	BOOL_OPT opt_vsids_en("vsids", "enable VSIDS (VMFQ otherwise)", true);
	BOOL_OPT opt_vsidsonly_en("vsidsonly", "enable VSIDS only (VMFQ disabled)", false);
	BOOL_OPT opt_mdmvsidsonly_en("eligible-vsidsonly", "enable VSIDS only in MDM (VMFQ disabled)", false);
	BOOL_OPT opt_mdmfusem_en("mdm-fuse-master", "enable MDM fusing in master mode (relies on conflicts only)", false);
	BOOL_OPT opt_mdmfuses_en("mdm-fuse-slave", "enable MDM fusing in slave mode (slave to restarts and conflicts)", true);
	BOOL_OPT opt_guess_en("guess", "enable variable phase guessing", true);
	BOOL_OPT opt_polarity("polarity", "initial variable polarity", true);
	BOOL_OPT opt_rephase_en("rephase", "enable variable rephasing", true);
	BOOL_OPT opt_reduce_en("reduce", "enable learnt database reduction", true);
	BOOL_OPT opt_stable_en("stable", "enable variable phases stabilization based on restarts", true);
	BOOL_OPT opt_targetphase_en("target-phase", "use target phase in decision making", true);
	BOOL_OPT opt_model_en("model", "print model on stdout", false);
	BOOL_OPT opt_proof_en("proof", "generate proof in binary DRAT format", false);
	INT_OPT opt_progress("progress-rate", "progress rate to print search statistics", 15000, INT32R(1, INT32_MAX));
	INT_OPT opt_mdm_vsidspumps("mdm-vsids-pumps", "set the number of follow-up decision pumps using VSIDS activity", 1, INT32R(0, INT32_MAX));
	INT_OPT opt_mdm_vmfqpumps("mdm-vmfq-pumps", "set the number of follow-up decision pumps using VMFQ activity", 1, INT32R(0, INT32_MAX));
	INT_OPT opt_mdm_rounds("mdm", "set the number of mdm rounds in a single search", 3, INT32R(0, INT32_MAX));
	INT_OPT opt_mdm_freq("mdm-freq", "MDM frequency based on conflicts", 100000, INT32R(1, INT32_MAX));
	INT_OPT opt_mdm_div("mdm-div", "MDM frequency divider", 1, INT32R(0, INT32_MAX));
	INT_OPT opt_mdm_sinc("mdm-sinc", "MDM divider increment in slave mode", 5, INT32R(0, INT32_MAX));
	INT_OPT opt_mdm_minc("mdm-minc", "MDM conflicts limit increment in master mode", 50000, INT32R(0, INT32_MAX));
	INT_OPT opt_timeout("timeout", "set the timeout in seconds", 0, INT32R(0, INT32_MAX));
	INT_OPT opt_seed("seed", "seed value for random generation", 9999991, INT32R(1, INT32_MAX));
	INT_OPT opt_map_min("map-min", "minimum variables to map", 2000, INT32R(0, INT32_MAX));
	INT_OPT opt_shrink_min("shrink-min", "minimum variables to shrink", 50, INT32R(0, INT32_MAX));
	INT_OPT opt_sigma_min("sigma-min", "minimum root variables to awaken live sigma", 100, INT32R(1, INT32_MAX));
	INT_OPT opt_sigma_inc("sigma-inc", "live sigma increment value based on conflicts", 4000, INT32R(1, INT32_MAX));
	INT_OPT opt_reduce_inc("reduce-inc", "increment value of learnt clauses deletion", 300, INT32R(10, INT32_MAX));
	INT_OPT opt_lbd_tier1("lbd-tier1", "lbd value of tier 1 learnts", 2, INT32R(1, INT32_MAX));
	INT_OPT opt_lbd_tier2("lbd-tier2", "lbd value of tier 2 learnts", 6, INT32R(3, INT32_MAX));
	INT_OPT opt_lbd_min_size("lbd-min-size", "minimum clause size to minimize a learnt", 30, INT32R(3, INT32_MAX));
	INT_OPT opt_lbd_min("lbd-min", "minimum lbd value to minimize a learnt", 6, INT32R(3, INT32_MAX));
	INT_OPT opt_rephase_inc("rephase-inc", "rephasing increment value", 1000, INT32R(100, INT32_MAX));
	INT_OPT opt_luby_inc("luby-inc", "luby increment value", 1024, INT32R(1, INT32_MAX));
	INT_OPT opt_luby_max("luby-max", "luby sequence maximum value", 1048576, INT32R(1, INT32_MAX));
	INT_OPT opt_rest_inc("pow-restart-inc", "power restart increment value", 2, INT32R(0, INT32_MAX));
	INT_OPT opt_rest_base("pow-restart-base", "power restart initial value", 2, INT32R(0, INT32_MAX));
	INT_OPT opt_stabrestart_r("stable-restart-rate", "stable restart increase rate", 2, INT32R(1, INT32_MAX));
	INT_OPT opt_stabrestart_inc("stable-restart-inc", "stable restart increment value", 1000, INT32R(1, INT32_MAX));
	INT_OPT opt_lbdfast("lbd-fast", "initial lbd fast window", 50, INT32R(1, 100));
	INT_OPT opt_lbdslow("lbd-slow", "initial lbd slow window", 100000, INT32R(1000, INT32_MAX));
	DOUBLE_OPT opt_reduce_perc("reduce-perc", "percentage of learnt clauses to reduce", 0.75, FP64R(0.1, 1));
	DOUBLE_OPT opt_lbd_rate("lbd-rate", "slow rate in firing lbd restarts", 1.1, FP64R(1, 10));
	DOUBLE_OPT opt_var_inc("var-inc", "VSIDS increment value", 1.0, FP64R(1, 10));
	DOUBLE_OPT opt_var_decay("var-decay", "VSIDS decay value", 0.95, FP64R(0, 1));
	DOUBLE_OPT opt_garbage_perc("gc-perc", "collect garbage if its percentage exceeds this value", 0.4, FP64R(0, 1));
	STRING_OPT opt_proof_out("proof-file", "output file to write binary proof", "proof.out");

}

#endif

