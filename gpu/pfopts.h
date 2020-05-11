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

#ifndef __PFRST_OPTS_
#define __PFRST_OPTS_

#include "pfargs.h"

BOOL_OPT opt_quiet_en("q", "enable quiet mode, same as verbose=0", false);
BOOL_OPT opt_pre_en("pre", "enable preprocessing", false);
BOOL_OPT opt_lpre_en("lpre", "enable preprocessing on learnt cls if exist", false);
BOOL_OPT opt_perf_en("perf", "allow performance report on stdout", true);
BOOL_OPT opt_par_en("par", "parse only the input formula", false);
BOOL_OPT opt_lcv_en("lcv", "use least-constrained variables to make parallel decisions", false);
BOOL_OPT opt_fdp_en("fdp", "enable follow-up decision prioritization", true);
BOOL_OPT opt_cbt_en("cbt", "enable chronological backtracking", false);
BOOL_OPT opt_model_en("m", "allow removed print on stdout", true);
BOOL_OPT opt_proof_en("p", "enable proof generation in binary DRAT format", true);
INT_OPT opt_verbose("verbose", "set the verbosity", 1, INT32R(0, 4));
INT_OPT opt_progress("progress-rate", "progress rate to print search statistics", 10000, INT32R(1, INT32_MAX));
INT_OPT opt_pdm_rounds("pdm", "set the number <n> of pdm rounds in a single search", 3, INT32R(0, INT32_MAX));
INT_OPT opt_pdm_freq("pdm-freq", "PDM execution frequency per restarts (0 for conflicts slave)", 0, INT32R(0, INT32_MAX));
INT_OPT opt_pdm_ord("pdm-ord", "set the pdm ordering scheme to either 1 (hist + act), 0 (high act) or -1 (low act)", 1, INT32R(-1, 1));
INT_OPT opt_pol("pol", "set polarity saving to either 1 (same), 0 (random) or -1 (revert)", 1, INT32R(-1, 1));
INT_OPT opt_timeout("timeout", "set the timeout in seconds", 0, INT32R(0, INT32_MAX));
INT_OPT opt_luby_base("luby-base", "set the restart base", 100, INT32R(32, INT32_MAX));
INT_OPT opt_SH("sh", "define the search heuristics, where:\nc |\
	0 -> clause activity heuristic only,\nc |\
	1 -> clause-sized random deletion heuristic,\nc |\
	2 -> LBD heuristic.", 2, INT32R(0, 2));
INT_OPT opt_seed("seed", "seed value for random generation", 9453921, INT32R(1, INT32_MAX));
INT_OPT opt_pre_delay("pre-delay", "delay for applying preprocessing (in restarts)", 50, INT32R(0, INT32_MAX));
INT_OPT opt_init_red("init-reduce", "initial number of conflicts for learnts reduction", 2000, INT32R(0, INT32_MAX));
INT_OPT opt_inc_red_sm("inc-small", "small step for learnt cls deletion", 300, INT32R(0, INT32_MAX));
INT_OPT opt_inc_red_bg("inc-big", "large step for learnt cls deletion", 1000, INT32R(0, INT32_MAX));
INT_OPT opt_lbd_frozen("lbd-frozen", "freeze cls if their LBD decreases than this value", 30, INT32R(0, INT32_MAX));
INT_OPT opt_lbd_min_size("lbd-min-size", "minimum size required to minimize clause", 30, INT32R(3, INT32_MAX));
INT_OPT opt_lbd_min("lbd-min", "minimum LBD required to minimize clause", 6, INT32R(3, INT32_MAX));
INT_OPT opt_lbd_rest_base("lbd-restart", "base of restarts (initial LBD queue size)", 50, INT32R(10, INT32_MAX));
INT_OPT opt_bl_rest_min("bl-rest-min", "minimum conflicts to block restarts", 10000, INT32R(1000, INT32_MAX));
INT_OPT opt_bl_rest_base("bl-rest-base", "base of block restarts (initial trail queue size)", 5000, INT32R(10, INT32_MAX));
INT_OPT opt_var_dfreq("var-dfreq", "minimum conflicts to increase VSIDS decay", 5000, INT32R(1000, INT32_MAX));
INT_OPT opt_cbt_dist("cbt-dist", "maximum distance (level difference) to activate CBT", 500, INT32R(GLUE, INT32_MAX));
INT_OPT opt_cbt_confs("cbt-conf", "maximum number of conflicts to activate CBT", 5000, INT32R(0, INT32_MAX));
DOUBLE_OPT opt_RF("rf-rate", "restart rate for trigger", 0.8, FP64R(0, 1));
DOUBLE_OPT opt_RB("rb-rate", "restart rate for defuse", 1.4, FP64R(1, 5));
DOUBLE_OPT opt_luby_inc("luby-inc", "luby restart increment value", 2.0);
DOUBLE_OPT opt_var_inc("var-inc", "VSIDS increment value", 1.0, FP64R(1, 10));
DOUBLE_OPT opt_var_decay("var-decay", "VSIDS decay value", 0.7, FP64R(0, 1));
DOUBLE_OPT opt_var_decay_r("var-decay-r", "VSIDS decay rate value", 0.006, FP64R(0, 1));
DOUBLE_OPT opt_garbage_perc("gc-perc", "Collect garbage if its percentage exceeds this value", 0.25, FP64R(0, 1));
STRING_OPT opt_restart("restart", "enables <policy> restart, where:\nc |\
	pow  -> geometric restarts,\nc |\
	luby -> luby restarts,\nc |\
	lbd  -> dynamic lbd restarts.", "lbd");
STRING_OPT opt_proof_out("proof-file", "output file to write binary proof", "proof.out");

#endif

