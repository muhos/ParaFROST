/***********************************************************************[pfoptions.cpp]
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

#include "pfoptions.h"
#include "pfconst.cuh"

using namespace pFROST;

// simplifier options
BOOL_OPT opt_ve_en("ve", "enable bounded variable elimination (BVE)", true);
BOOL_OPT opt_ve_lbound_en("velbound", "skip variables resulting in more literals than original", false);
BOOL_OPT opt_ve_plus_en("ve+", "enable HSE + BVE", true);
BOOL_OPT opt_hse_en("hse", "enable hybrid subsumption elimination", true);
BOOL_OPT opt_bce_en("bce", "enable blocked clause elimination", false);
BOOL_OPT opt_ere_en("ere", "enable eager redundancy elimination", true);
BOOL_OPT opt_all_en("all", "enable all simplifications", false);
BOOL_OPT opt_gc_gpu_en("gcgpu", "enable device garbage collection in parallel (stream compaction approach)", true);
BOOL_OPT opt_profile_simp_en("profilesimp", "profile simplifications", false);
BOOL_OPT opt_aggr_cnf_sort("aggresivesort", "sort simplified formula with aggresive key before writing to host", false);
BOOL_OPT opt_unified_access_en("unifiedaccess", "enable unified access in transferring CNF to/from device", false);
BOOL_OPT opt_sync_always_en("syncalways", "synchronize all device calls with the host (useful for debugging)", false);
BOOL_OPT opt_atomic_ve_en("atomicve", "enable atomic version of BVE (produces out-of-order resolvents)", false);
BOOL_OPT opt_solve_en("solve", "proceed with solving after simplifications", true);
INT_OPT opt_lcve_min("lcvemin", "minimum parallel variables to simplify", 2, INT32R(1, INT32_MAX));
INT_OPT opt_ve_clause_max("veclausemax", "maximum resolvent size (0: no limit)", 100, INT32R(0, INT32_MAX));
INT_OPT opt_ve_phase_min("vephasemin", "minimum removed literals to stop reductions", 500, INT32R(1, INT32_MAX));
INT_OPT opt_mu_pos("mupos", "set the positive freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_mu_neg("muneg", "set the negative freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_xor_max_arity("xormaxarity", "maximum XOR fanin size", 10, INT32R(2, INT32_MAX));
INT_OPT opt_hse_max_occurs("hsemax", "maximum occurrence list size to scan in HSE", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_bce_max_occurs("bcemax", "maximum occurrence list size to scan in BCE", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_ere_max_occurs("eremax", "maximum occurrence list size to scan in ERE", 3e3, INT32R(100, INT32_MAX));
INT_OPT opt_phases("phases", "set the number of phases in to run reductions", 5, INT32R(0, INT32_MAX));
INT_OPT opt_cnf_free("gcfreq", "set the frequency of CNF memory shrinkage on GPU", 2, INT32R(0, 5));
INT_OPT opt_gpus("ngpus", "number of GPU(s) to be activated", 1, INT32R(1, 32));
INT_OPT opt_streams("nstreams", "number of GPU streams to be created", 6, INT32R(1, 32));
DOUBLE_OPT opt_lits_mul("literalsmul", "initial literals multiplier", 1.2, FP64R(0, 4));

// solver options
BOOL_OPT opt_boundsearch_en("boundsearch", "activate search bounds on decisions and/or conflicts", false);
BOOL_OPT opt_chrono_en("chrono", "enable chronological backtracking", true);
BOOL_OPT opt_chronoreuse_en("chronoreusetrail", "enable reuse trail when chronological backtracking", false);
BOOL_OPT opt_bumpreason_en("bumpreason", "bump reason literals via learnt clause", true);
BOOL_OPT opt_debinary_en("debinary", "remove duplicated binaries", true);
BOOL_OPT opt_decompose_en("decompose", "decompose binary implication gragh into SCCs", true);
BOOL_OPT opt_targetphase_en("targetphase", "use target phase in decision making", true);
BOOL_OPT opt_ternary_en("ternary", "enable hyper ternary resolution", true);
BOOL_OPT opt_ternary_sleep_en("ternarysleep", "allow hyper ternary resolution to sleep", true);
BOOL_OPT opt_transitive_en("transitive", "enable transitive reduction on binary implication graph", true);
BOOL_OPT opt_parseonly_en("parseonly", "parse only the input formula", false);
BOOL_OPT opt_polarity("polarity", "initial variable polarity", true);
BOOL_OPT opt_proof_en("proof", "generate proof in binary DRAT format", false);
BOOL_OPT opt_probe_en("probe", "enable failed literal probing", true);
BOOL_OPT opt_probe_sleep_en("probesleep", "allow failed literal probing to sleep", true);
BOOL_OPT opt_probehbr_en("probehyper", "learn hyper binary clauses", true);
BOOL_OPT opt_model_en("model", "extend model with eliminated variables", true);
BOOL_OPT opt_modelprint_en("modelprint", "print model on stdout", false);
BOOL_OPT opt_modelverify_en("modelverify", "verify model on input formula", false);
BOOL_OPT opt_mdmlcv_en("mdmlcv", "use least-constrained variables to make multiple decisions", false);
BOOL_OPT opt_mdmvsidsonly_en("mdmvsidsonly", "enable VSIDS only in MDM (VMFQ disabled)", false);
BOOL_OPT opt_mdmfusem_en("mdmfusemaster", "enable MDM fusing in master mode (relies on conflicts only)", true);
BOOL_OPT opt_mdmfuses_en("mdmfuseslave", "enable MDM fusing in slave mode (slave to restarts and conflicts)", false);
BOOL_OPT opt_mdmassume_en("mdmassume", "choose multiple decisions based on given assumptions (incremental mode)", false);
BOOL_OPT opt_report_en("report", "allow performance report on stdout", true);
BOOL_OPT opt_rephase_en("rephase", "enable variable rephasing", true);
BOOL_OPT opt_reduce_en("reduce", "enable learnt database reduction", true);
BOOL_OPT opt_reusetrail_en("reusetrail", "reuse part of the trail upon restarts", true);
BOOL_OPT opt_sigpre_en("sigma", "enable preprocessing using SIGmA", true);
BOOL_OPT opt_siglive_en("sigmalive", "enable live SIGmA (inprocessing)", true);
BOOL_OPT opt_sigsleep_en("sigmasleep", "allow SIGmA to sleep", true);
BOOL_OPT opt_subsume_en("subsume", "enable forward subsumption elimination", true);
BOOL_OPT opt_stable_en("stable", "enable variable phases stabilization based on restarts", true);
BOOL_OPT opt_vsids_en("vsids", "enable VSIDS (VMFQ otherwise)", true);
BOOL_OPT opt_vsidsonly_en("vsidsonly", "enable VSIDS only (VMFQ disabled)", false);
BOOL_OPT opt_vivify_en("vivify", "enable vivification", true);

INT_OPT opt_chrono_min("chronomin", "minimum distance to trigger chronological backtracking", 100, INT32R(0, INT32_MAX));
INT_OPT opt_decompose_min("decomposemin", "minimum rounds to decompose", 2, INT32R(1, 10));
INT_OPT opt_decompose_limit("decomposelimit", "decompose round limit", 1e7, INT32R(0, 10));
INT_OPT opt_decompose_min_eff("decomposemineff", "decompose minimum efficiency", 1e7, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_vsidspumps("mdmvsidspumps", "set the number of follow-up decision pumps using VSIDS activity", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_vmfqpumps("mdmvmfqpumps", "set the number of follow-up decision pumps using VMFQ activity", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_rounds("mdmrounds", "set the number of mdm rounds in a single search", 0, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_freq("mdmfreq", "MDM frequency based on conflicts", 1e5, INT32R(1, INT32_MAX));
INT_OPT opt_mdm_div("mdmdiv", "MDM frequency divider", 1, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_sinc("mdmsinc", "MDM divider increment in slave mode", 5, INT32R(0, INT32_MAX));
INT_OPT opt_mdm_minc("mdminc", "MDM conflicts limit increment in master mode", 2e3, INT32R(0, INT32_MAX));
INT_OPT opt_minimize_depth("minimizedepth", "minimization depth to explore", 1e3, INT32R(1, INT32_MAX));
INT_OPT opt_mode_inc("modeinc", "mode increment value based on conflicts", 1e3, INT32R(1, INT32_MAX));
INT_OPT opt_nap("nap", "maximum naping period", 2, INT32R(0, 10));
INT_OPT opt_ternary_priorbins("ternarypriorbins", "prioritize binaries in watch table after hyper ternary resolution (1: enable, 2: prioritize learnts)", 0, INT32R(0, 2));
INT_OPT opt_ternary_max_eff("ternarymaxeff", "maximum hyper ternary resolution efficiency", 1e2, INT32R(0, INT32_MAX));
INT_OPT opt_ternary_min_eff("ternarymineff", "minimum hyper ternary resolution efficiency", 1e6, INT32R(0, INT32_MAX));
INT_OPT opt_ternary_rel_eff("ternaryreleff", "relative hyper ternary resolution efficiency per mille", 40, INT32R(0, INT32_MAX));
INT_OPT opt_transitive_max_eff("transitivemaxeff", "maximum transitive efficiency", 1e2, INT32R(0, INT32_MAX));
INT_OPT opt_transitive_min_eff("transitivemineff", "minimum transitive efficiency", 1e6, INT32R(0, INT32_MAX));
INT_OPT opt_transitive_rel_eff("transitivereleff", "relative transitive efficiency per mille", 20, INT32R(0, INT32_MAX));
INT_OPT opt_restart_inc("restartinc", "restart increment value based on conflicts", 1, INT32R(1, INT32_MAX));
INT_OPT opt_reduce_inc("reduceinc", "increment value of clauses reduction based on conflicts", 300, INT32R(10, INT32_MAX));
INT_OPT opt_rephase_inc("rephaseinc", "rephasing increment value based on conflicts", 1e3, INT32R(100, INT32_MAX));
INT_OPT opt_progress("progressrate", "progress rate to print search statistics", 1e4, INT32R(1, INT32_MAX));
INT_OPT opt_probe_inc("probeinc", "probe increment value based on conflicts", 100, INT32R(1, INT32_MAX));
INT_OPT opt_probe_min("probemin", "minimum rounds to probe", 2, INT32R(1, 10));
INT_OPT opt_probe_max_eff("probemaxeff", "maximum probe efficiency", 1e2, INT32R(0, INT32_MAX));
INT_OPT opt_probe_min_eff("probemineff", "minimum probe efficiency", 5e5, INT32R(0, INT32_MAX));
INT_OPT opt_probe_rel_eff("probereleff", "relative probe efficiency per mille", 2, INT32R(0, INT32_MAX));
INT_OPT opt_seed("seed", "seed value for random generation", 0, INT32R(0, INT32_MAX));
INT_OPT opt_sigma_inc("sigmainc", "live sigma increment value based on conflicts", 500, INT32R(1, INT32_MAX));
INT_OPT opt_sigma_min("sigmamin", "minimum root variables shrunken to awaken SIGmA", 4e3, INT32R(1, INT32_MAX));
INT_OPT opt_sigma_priorbins("sigmapriorbins", "prioritize binaries in watch table after sigmification (1: enable, 2: prioritize learnts)", 1, INT32R(0, 2));
INT_OPT opt_subsume_priorbins("subsumepriorbins", "prioritize binaries in watch table after subsume (1: enable, 2: prioritize learnts)", 1, INT32R(0, 2));
INT_OPT opt_subsume_inc("subsumeinc", "forward subsumption increment value based on conflicts", 2e3, INT32R(100, INT32_MAX));
INT_OPT opt_subsume_min_occs("subsumeminoccurs", "minimum occurrences to subsume or strengthen", 3e3, INT32R(10, INT32_MAX));
INT_OPT opt_subsume_max_csize("subsumemaxcsize", "maximum subsuming clause size", 1e3, INT32R(2, INT32_MAX));
INT_OPT opt_subsume_max_eff("subsumemaxeff", "maximum number of clauses to scan in subsume", 1e2, INT32R(0, INT32_MAX));
INT_OPT opt_subsume_min_eff("subsumemineff", "minimum number of clauses to scan in subsume", 1e6, INT32R(0, INT32_MAX));
INT_OPT opt_subsume_rel_eff("subsumereleff", "relative subsume efficiency per mille", 1e4, INT32R(0, INT32_MAX));
INT_OPT opt_vivify_priorbins("vivifypriorbins", "prioritize binaries in watch table before vivification (1: enable, 2: prioritize learnts)", 0, INT32R(0, 2));
INT_OPT opt_vivify_max_eff("vivifymaxeff", "maximum number of clauses to scan in subsume", 50, INT32R(0, INT32_MAX));
INT_OPT opt_vivify_min_eff("vivifymineff", "minimum number of clauses to scan in subsume", 2e5, INT32R(0, INT32_MAX));
INT_OPT opt_vivify_rel_eff("vivifyreleff", "relative subsume efficiency per mille", 2, INT32R(0, INT32_MAX));
INT_OPT opt_lbd_tier1("lbdtier1", "lbd value of tier 1 learnts", 2, INT32R(1, INT32_MAX));
INT_OPT opt_lbd_tier2("lbdtier2", "lbd value of tier 2 learnts", 6, INT32R(3, INT32_MAX));
INT_OPT opt_lbd_fast("lbdfast", "initial lbd fast window", 33, INT32R(1, 100));
INT_OPT opt_lbd_slow("lbdslow", "initial lbd slow window", 1e5, INT32R(100, INT32_MAX));
INT_OPT opt_luby_inc("lubyinc", "luby increment value based on conflicts", 1 << 10, INT32R(1, INT32_MAX));
INT_OPT opt_luby_max("lubymax", "luby sequence maximum value", 1 << 20, INT32R(1, INT32_MAX));
INT_OPT opt_learntsub_max("subsumelearntmax", "maximum learnt clauses to subsume", 20, INT32R(0, INT32_MAX));
INT64_OPT opt_conflictout("conflictout", "set out-of-conflicts limit (must be enabled by \"boundsearch\")", INT64_MAX, INT64R(0, INT64_MAX));
INT64_OPT opt_decisionout("decisionout", "set out-of-decisions limit (must be enabled by \"boundsearch\")", INT64_MAX, INT64R(0, INT64_MAX));
DOUBLE_OPT opt_stable_rate("stablerestartrate", "stable restart increase rate", 1.0, FP64R(1, 5));
DOUBLE_OPT opt_lbd_rate("lbdrate", "slow rate in firing lbd restarts", 1.1, FP64R(1, 10));
DOUBLE_OPT opt_ternary_perc("ternaryperc", "percentage of maximum hyper clauses to add", 0.2, FP64R(0, 1));
DOUBLE_OPT opt_map_perc("mapperc", "minimum percentage of variables to map", 0.2, FP64R(0, 1));
DOUBLE_OPT opt_reduce_perc("reduceperc", "percentage of learnt clauses to reduce", 0.75, FP64R(0.1, 1));
DOUBLE_OPT opt_var_inc("varinc", "VSIDS increment value", 1.0, FP64R(1, 10));
DOUBLE_OPT opt_var_decay("vardecay", "VSIDS decay value", 0.95, FP64R(0, 1));
DOUBLE_OPT opt_garbage_perc("garbageperc", "collect garbage if its percentage exceeds this value", 0.25, FP64R(0, 1));
STRING_OPT opt_proof_out("proofout", "output file to write binary proof", "proof.out");

OPTION::OPTION() {
	RESETSTRUCT(this);
	int MAXLEN = 256;
	proof_path = pfcalloc<char>(MAXLEN);
}

OPTION::~OPTION() {
	if (proof_path != NULL) {
		std::free(proof_path);
		proof_path = NULL;
	}
}

void OPTION::init() 
{
	bumpreason_en = opt_bumpreason_en;
	chrono_en = opt_chrono_en;
	chronoreuse_en = opt_chronoreuse_en;
	chrono_min = opt_chrono_min;
	conflict_out = opt_conflictout;
	decision_out = opt_decisionout;
	debinary_en = opt_debinary_en;
	decompose_en = opt_decompose_en;
	decompose_min = opt_decompose_min;
	decompose_limit = opt_decompose_limit;
	decompose_min_eff = opt_decompose_min_eff;
	model_en = opt_model_en;
	modelprint_en = opt_modelprint_en;
	modelverify_en = opt_modelverify_en;
	mode_inc = opt_mode_inc;
	minimize_depth = opt_minimize_depth;
	mdmvsidsonly_en = opt_mdmvsidsonly_en;
	mdmassume_en = opt_mdmassume_en;
	mdmfusem_en = opt_mdmfusem_en;
	mdmfuses_en = opt_mdmfuses_en;
	mdm_mcv_en = !opt_mdmlcv_en;
	mdm_vsids_pumps = opt_mdm_vsidspumps;
	mdm_vmfq_pumps = opt_mdm_vmfqpumps;
	mdm_rounds = opt_mdm_rounds;
	mdm_freq = opt_mdm_freq;
	mdm_inc = opt_mdm_minc;
	mdm_sinc = opt_mdm_sinc;
	mdm_div = opt_mdm_div;
	map_perc = opt_map_perc;
	nap = opt_nap;
	memcpy(proof_path, opt_proof_out, opt_proof_out.length());
	parseonly_en = opt_parseonly_en;
	proof_en = opt_proof_en;
	probe_en = opt_probe_en;
	probe_sleep_en = opt_probe_sleep_en;
	probehbr_en = opt_probehbr_en;
	probe_inc = opt_probe_inc;
	probe_min = opt_probe_min;
	probe_min_eff = opt_probe_min_eff;
	probe_max_eff = opt_probe_max_eff;
	probe_rel_eff = opt_probe_rel_eff;
	prograte = opt_progress;
	polarity = opt_polarity;
	targetphase_en = opt_targetphase_en;
	ternary_en = opt_ternary_en;
	ternary_sleep_en = opt_ternary_sleep_en;
	ternary_priorbins = opt_ternary_priorbins;
	ternary_min_eff = opt_ternary_min_eff;
	ternary_max_eff = opt_ternary_max_eff;
	ternary_rel_eff = opt_ternary_rel_eff;
	ternary_perc = opt_ternary_perc;
	transitive_en = opt_transitive_en;
	transitive_min_eff = opt_transitive_min_eff;
	transitive_max_eff = opt_transitive_max_eff;
	transitive_rel_eff = opt_transitive_rel_eff;
	vsids_en = opt_vsids_en;
	vsidsonly_en = opt_vsidsonly_en;
	var_inc = opt_var_inc;
	var_decay = opt_var_decay;
	vivify_en = opt_vivify_en;
	vivify_priorbins = opt_vivify_priorbins;
	vivify_min_eff = opt_vivify_min_eff;
	vivify_max_eff = opt_vivify_max_eff;
	vivify_rel_eff = opt_vivify_rel_eff;
	report_en = opt_report_en && !quiet_en;
	reduce_en = opt_reduce_en;
	reduce_perc = opt_reduce_perc;
	reduce_inc = opt_reduce_inc;
	rephase_en = opt_rephase_en;
	rephase_inc = opt_rephase_inc;
	reusetrail_en = opt_reusetrail_en;
	restart_inc = opt_restart_inc;
	stable_en = opt_stable_en;
	stable_rate = opt_stable_rate;
	sigma_en = opt_sigpre_en;
	sigma_live_en = opt_siglive_en;
	sigma_sleep_en = opt_sigsleep_en;
	sigma_inc = opt_sigma_inc;
	sigma_min = opt_sigma_min;
	sigma_priorbins = opt_sigma_priorbins;
	subsume_en = opt_subsume_en;
	subsume_inc = opt_subsume_inc;
	subsume_priorbins = opt_subsume_priorbins;
	subsume_min_occs = opt_subsume_min_occs;
	subsume_rel_eff = opt_subsume_rel_eff;
	subsume_min_eff = opt_subsume_min_eff;
	subsume_max_eff = opt_subsume_max_eff;
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
	gc_perc = opt_garbage_perc;
	// initialize simplifier options
	if (sigma_en || sigma_live_en) {
		all_en = opt_all_en;
		aggr_cnf_sort = opt_aggr_cnf_sort;
		bce_en = opt_bce_en;
		bce_limit = opt_bce_max_occurs;
		ere_en = opt_ere_en;
		ere_limit = opt_ere_max_occurs;
		hse_en = opt_hse_en && ve_plus_en;
		hse_limit = opt_hse_max_occurs;
		lcve_min = opt_lcve_min;
		lits_min = opt_ve_phase_min;
		solve_en = opt_solve_en;
		shrink_rate = opt_cnf_free;
		mu_pos = opt_mu_pos;
		mu_neg = opt_mu_neg;
		phases = opt_phases;
		profile_simp = opt_profile_simp_en;
		ve_plus_en = opt_ve_plus_en;
		ve_en = opt_ve_en || ve_plus_en;
		ve_lbound_en = opt_ve_lbound_en;
		ve_clause_limit = opt_ve_clause_max;
		xor_max_arity = opt_xor_max_arity;
		if (all_en) ve_en = 1, ve_plus_en = 1, bce_en = 1, ere_en = 1;
		if (!phases && (ve_en || hse_en || bce_en)) phases = 1; // at least 1 phase needed
		if (phases && !(ve_en || hse_en || bce_en)) phases = 0;
		if (phases > 1 && !ve_en) phases = 1;
		ngpus = opt_gpus;
		gc_gpu = opt_gc_gpu_en;
		nstreams = opt_streams;
		lits_mul = opt_lits_mul;
		atomic_ve = opt_atomic_ve_en;
		sync_always = opt_sync_always_en;
		unified_access = opt_unified_access_en;
	}
}