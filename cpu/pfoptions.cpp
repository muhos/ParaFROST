/***********************************************************************[pfoptions.cpp]
Copyright(c) 2020, Muhammad Osama  Anton Wijs,
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

#include "pfoptions.h"

using namespace pFROST;

// simplifier options
BOOL_OPT opt_ve_en("ve", "enable bounded variable elimination (BVE)", true);
BOOL_OPT opt_ve_plus_en("ve+", "enable HSE + BVE", true);
BOOL_OPT opt_hse_en("hse", "enable hybrid subsumption elimination", true);
BOOL_OPT opt_bce_en("bce", "enable blocked clause elimination", false);
BOOL_OPT opt_ere_en("ere", "enable eager redundancy elimination", true);
BOOL_OPT opt_all_en("all", "enable all simplifications", false);
BOOL_OPT opt_profile_simp_en("profilesimp", "profile simplifications", false);
BOOL_OPT opt_aggr_cnf_sort("aggresivesort", "sort simplified formula with aggresive key before writing to host", false);
BOOL_OPT opt_solve_en("solve", "proceed with solving after simplifications", true);
INT_OPT opt_lcve_min("lcvemin", "minimum parallel variables to simplify", 2, INT32R(1, INT32_MAX));
INT_OPT opt_ve_phase_min("vephasemin", "minimum removed literals to stop stage1 reductions ", 500, INT32R(1, INT32_MAX));
INT_OPT opt_mu_pos("mupos", "set the positive freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_mu_neg("muneg", "set the negative freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_xor_max_arity("xormaxarity", "maximum XOR fanin size", 20, INT32R(2, INT32_MAX));
INT_OPT opt_hse_max_occurs("hsemax", "maximum occurrence list size to scan in HSE", 30000, INT32R(100, INT32_MAX));
INT_OPT opt_bce_max_occurs("bcemax", "maximum occurrence list size to scan in BCE", 30000, INT32R(100, INT32_MAX));
INT_OPT opt_ere_max_occurs("eremax", "maximum occurrence list size to scan in ERE", 30000, INT32R(100, INT32_MAX));
INT_OPT opt_phases("phases", "set the number of phases in stage1 reductions", 5, INT32R(0, INT32_MAX));
INT_OPT opt_cnf_free("gcfreq", "set the frequency of CNF memory shrinkage in SIGmA", 2, INT32R(0, 5));

// solver options
BOOL_OPT opt_sig_pre_en("sigma", "enable preprocessing using SIGmA", false);
BOOL_OPT opt_sig_live_en("sigmalive", "enable live SIGmA (inprocessing)", false);
BOOL_OPT opt_subsume_en("subsume", "enable forward subsumption elimination", true);
BOOL_OPT opt_report_en("report", "allow performance report on stdout", true);
BOOL_OPT opt_parseonly_en("parseonly", "parse only the input formula", false);
BOOL_OPT opt_vsids_en("vsids", "enable VSIDS (VMFQ otherwise)", true);
BOOL_OPT opt_vsidsonly_en("vsidsonly", "enable VSIDS only (VMFQ disabled)", false);
BOOL_OPT opt_mdm_lcv_en("mdmlcv", "use leastconstrained variables to make multiple decisions", false);
BOOL_OPT opt_mdmvsidsonly_en("mdmvsidsonly", "enable VSIDS only in MDM (VMFQ disabled)", false);
BOOL_OPT opt_mdmfusem_en("mdmfusemaster", "enable MDM fusing in master mode (relies on conflicts only)", true);
BOOL_OPT opt_mdmfuses_en("mdmfuseslave", "enable MDM fusing in slave mode (slave to restarts and conflicts)", false);
BOOL_OPT opt_polarity("polarity", "initial variable polarity", true);
BOOL_OPT opt_rephase_en("rephase", "enable variable rephasing", true);
BOOL_OPT opt_reduce_en("reduce", "enable learnt database reduction", true);
BOOL_OPT opt_reusetrail_en("reusetrail", "reuse part of the trail upon restarts", true);
BOOL_OPT opt_stable_en("stable", "enable variable phases stabilization based on restarts", true);
BOOL_OPT opt_targetphase_en("targetphase", "use target phase in decision making", true);
BOOL_OPT opt_bumpreason_en("bumpreason", "bump reason literals via learnt clause", true);
BOOL_OPT opt_chrono_en("chrono", "enable chronological backtracking", true);
BOOL_OPT opt_chronoreuse_en("chronoreusetrail", "enable reuse trail when chronological backtracking", false);
BOOL_OPT opt_model_en("model", "print model on stdout", false);
BOOL_OPT opt_proof_en("proof", "generate proof in binary DRAT format", false);
BOOL_OPT opt_priorbins_en("priorbins", "prioritize binaries in watch table", true);
INT_OPT opt_progress("progressrate", "progress rate to print search statistics", 15000, INT32R(1, INT32_MAX));
INT_OPT opt_seed("seed", "seed value for random generation", 0, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_vsidspumps("mdmvsidspumps", "set the number of followup decision pumps using VSIDS activity", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_vmfqpumps("mdmvmfqpumps", "set the number of followup decision pumps using VMFQ activity", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_rounds("mdm", "set the number of mdm rounds in a single search", 0, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_freq("mdmfreq", "MDM frequency based on conflicts", 100000, INT32R(1, INT32_MAX));
INT_OPT opt_mdm_div("mdmdiv", "MDM frequency divider", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_sinc("mdmsinc", "MDM divider increment in slave mode", 5, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_minc("mdmminc", "MDM conflicts limit increment in master mode", 20000, INT32R(0, INT32_MAX));
INT_OPT opt_map_min("mapmin", "minimum variables to map", 100, INT32R(0, INT32_MAX));
INT_OPT opt_map_inc("mapinc", "mapping increment value based on conflicts", 2000, INT32R(0, INT32_MAX));
INT_OPT opt_sigma_inc("sigmainc", "live sigma increment value based on conflicts", 2000, INT32R(1, INT32_MAX));
INT_OPT opt_sigma_min("sigmamin", "minimum root variables shrunken to awaken SIGmA", 4000, INT32R(1, INT32_MAX));
INT_OPT opt_chrono_min("chronomin", "minimum distance to trigger chronological backtracking", 100, INT32R(0, INT32_MAX));
INT_OPT opt_subsume_inc("subsumeinc", "forward subsumption increment value based on conflicts", 10000, INT32R(100, INT32_MAX));
INT_OPT opt_subsume_min_occs("subsumeminoccurs", "minimum occurrences to subsume or strengthen", 30000, INT32R(10, INT32_MAX));
INT_OPT opt_subsume_max_csize("subsumemaxcsize", "maximum subsuming clause size", 1000, INT32R(2, INT32_MAX));
INT_OPT opt_subsume_max_checks("subsumemaxchecks", "maximum number of clauses to scan in subsume", 100000000, INT32R(100, INT32_MAX));
INT_OPT opt_subsume_min_checks("subsumeminchecks", "minimum number of clauses to scan in subsume", 1000000, INT32R(100, INT32_MAX));
INT_OPT opt_reduce_inc("reduceinc", "increment value of clauses reduction based on conflicts", 300, INT32R(10, INT32_MAX));
INT_OPT opt_lbd_tier1("lbdtier1", "lbd value of tier 1 learnts", 2, INT32R(1, INT32_MAX));
INT_OPT opt_lbd_tier2("lbdtier2", "lbd value of tier 2 learnts", 6, INT32R(3, INT32_MAX));
INT_OPT opt_lbd_fast("lbdfast", "initial lbd fast window", 33, INT32R(1, 100));
INT_OPT opt_lbd_slow("lbdslow", "initial lbd slow window", 100000, INT32R(1000, INT32_MAX));
INT_OPT opt_luby_inc("lubyinc", "luby increment value based on conflicts", 1 << 10, INT32R(1, INT32_MAX));
INT_OPT opt_luby_max("lubymax", "luby sequence maximum value", 1 << 20, INT32R(1, INT32_MAX));
INT_OPT opt_learntsub_max("subsumelearntmax", "maximum learnt clauses to subsume", 20, INT32R(0, INT32_MAX));
INT_OPT opt_minimizebin_max("minimizebinmax", "maximum learnt clause size to minimize using binaries", 0, INT32R(0, INT32_MAX));
INT_OPT opt_minimize_depth("minimizedepth", "minimization depth to explore", 1000, INT32R(1, INT32_MAX));
INT_OPT opt_bumpreason_depth("bumpreasondepth", "bumping depth to explore", 1, INT32R(1, 5));
INT_OPT opt_rephase_inc("rephaseinc", "rephasing increment value based on conflicts", 1000, INT32R(100, INT32_MAX));
INT_OPT opt_powrestart_inc("powerrestartinc", "power restart increment value based on conflicts", 2, INT32R(0, INT32_MAX));
INT_OPT opt_stabrestart_inc("stablerestartinc", "stable restart increment value based on conflicts", 1000, INT32R(1, INT32_MAX));
DOUBLE_OPT opt_stabrestart_rate("stablerestartrate", "stable restart increase rate", 2.0, FP64R(1, 5));
DOUBLE_OPT opt_lbd_rate("lbdrate", "slow rate in firing lbd restarts", 1.1, FP64R(1, 10));
DOUBLE_OPT opt_map_perc("mapperc", "minimum percentage of variables to map", 0.1, FP64R(0, 1));
DOUBLE_OPT opt_reduce_perc("reduceperc", "percentage of learnt clauses to reduce", 0.75, FP64R(0.1, 1));
DOUBLE_OPT opt_var_inc("varinc", "VSIDS increment value", 1.0, FP64R(1, 10));
DOUBLE_OPT opt_var_decay("vardecay", "VSIDS decay value", 0.95, FP64R(0, 1));
DOUBLE_OPT opt_garbage_perc("gcperc", "collect garbage if its percentage exceeds this value", 0.25, FP64R(0, 1));
STRING_OPT opt_proof_out("prooffile", "output file to write binary proof", "proof.out");

void OPTION::init() {
	parse_only_en = opt_parseonly_en;
	priorbins_en = opt_priorbins_en;
	proof_path = opt_proof_out;
	proof_en = opt_proof_en;
	prograte = opt_progress;
	polarity = opt_polarity;
	vsids_en = opt_vsids_en;
	vsidsonly_en = opt_vsidsonly_en;
	var_inc = opt_var_inc;
	var_decay = opt_var_decay;
	model_en = opt_model_en;
	minimizebin_max = opt_minimizebin_max;
	minimize_depth = opt_minimize_depth;
	mdmvsidsonly_en = opt_mdmvsidsonly_en;
	mdmfusem_en = opt_mdmfusem_en;
	mdmfuses_en = opt_mdmfuses_en;
	mdm_mcv_en = !opt_mdm_lcv_en;
	mdm_vsids_pumps = opt_mdm_vsidspumps;
	mdm_vmfq_pumps = opt_mdm_vmfqpumps;
	mdm_rounds = opt_mdm_rounds;
	mdm_freq = opt_mdm_freq;
	mdm_minc = opt_mdm_minc;
	mdm_sinc = opt_mdm_sinc;
	mdm_div = opt_mdm_div;
	map_perc = opt_map_perc;
	map_min = opt_map_min;
	map_inc = opt_map_inc;
	chrono_en = opt_chrono_en;
	chronoreuse_en = opt_chronoreuse_en;
	chrono_min = opt_chrono_min;
	report_en = opt_report_en & !quiet_en;
	reduce_en = opt_reduce_en;
	reduce_perc = opt_reduce_perc;
	reduce_inc = opt_reduce_inc;
	rephase_en = opt_rephase_en;
	rephase_inc = opt_rephase_inc;
	reusetrail_en = opt_reusetrail_en;
	restart_inc = opt_powrestart_inc;
	stabrestart_rate = opt_stabrestart_rate;
	stabrestart_inc = opt_stabrestart_inc;
	stable_en = opt_stable_en;
	bumpreason_en = opt_bumpreason_en;
	bumpreason_depth = opt_bumpreason_depth;
	sigma_en = opt_sig_pre_en;
	sigma_live_en = opt_sig_live_en;
	sigma_inc = opt_sigma_inc;
	sigma_min = opt_sigma_min;
	subsume_en = opt_subsume_en;
	subsume_inc = opt_subsume_inc;
	subsume_min_occs = opt_subsume_min_occs;
	subsume_min_checks = opt_subsume_min_checks;
	subsume_max_checks = opt_subsume_max_checks;
	subsume_max_csize = opt_subsume_max_csize;
	seed = opt_seed;
	lbd_tier1 = opt_lbd_tier1;
	lbd_tier2 = opt_lbd_tier2;
	lbd_fast = opt_lbd_fast;
	lbd_slow = opt_lbd_slow;
	lbd_rate = opt_lbd_rate;
	luby_inc = opt_luby_inc;
	luby_max = opt_luby_max;
	learntsub_max = opt_learntsub_max;
	target_phase_en = opt_targetphase_en;
	gc_perc = opt_garbage_perc;
	// initialize simplifier options
	if (sigma_en || sigma_live_en) {
		solve_en = opt_solve_en;
		aggr_cnf_sort = opt_aggr_cnf_sort;
		profile_simp = opt_profile_simp_en;
		ve_plus_en = opt_ve_plus_en;
		bce_en = opt_bce_en;
		ere_en = opt_ere_en;
		all_en = opt_all_en;
		phases = opt_phases;
		mu_pos = opt_mu_pos;
		mu_neg = opt_mu_neg;
		lcve_min = opt_lcve_min;
		lits_min = opt_ve_phase_min;
		shrink_rate = opt_cnf_free;
		hse_limit = opt_hse_max_occurs;
		bce_limit = opt_bce_max_occurs;
		ere_limit = opt_ere_max_occurs;
		xor_max_arity = opt_xor_max_arity;
		ve_en = opt_ve_en || ve_plus_en;
		hse_en = opt_hse_en || ve_plus_en;
		if (all_en) ve_en = 1, ve_plus_en = 1, bce_en = 1, ere_en = 1;
		if (!phases && (ve_en | hse_en | bce_en)) phases = 1; // at least 1 phase needed
		if (phases && !(ve_en | hse_en | bce_en)) phases = 0;
		if (phases > 1 && !ve_en) phases = 1;
	}
}