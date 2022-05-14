/***********************************************************************[options.h]
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

#ifndef __SOLVER_OPTS_
#define __SOLVER_OPTS_

#include <cstring>
#include "datatypes.h"
#include "input.h"

namespace pFROST {

	struct OPTION {
		//==========================================//
		//             Solver options               //
		//==========================================//
		LIT_ST	polarity;
		//------------------------------------------//
		char*	proof_path;
		//------------------------------------------//
		int64	learntsub_max;
		//------------------------------------------//
		int64	sigma_min, sigma_inc;
		uint64	conflict_out, decision_out;
		//------------------------------------------//
		double	var_inc, var_decay;
		double	stable_rate;
		double	lbd_rate;
		double	gc_perc;
		double	map_perc;
		double	reduce_perc;
		double	ternary_perc;
		//------------------------------------------//
		int		nap;
		int		seed;
		int		prograte;
		int		mode_inc;
		int		chrono_min;
		int		reduce_inc;
		int		restart_inc;
		int		rephase_inc;
		int		decompose_min;
		int		decompose_limit;
		int		decompose_min_eff;
		int		sigma_priorbins;
		int		minimize_depth;
		int		minimize_min;
		int		minimize_lbd;
		int		luby_inc, luby_max;
		int		lbd_tier2, lbd_tier1, lbd_fast, lbd_slow;
		int		mdm_rounds, mdm_inc, mdm_vsids_pumps, mdm_vmtf_pumps;
		int		subsume_priorbins, subsume_inc, subsume_max_occs, subsume_min_eff, subsume_max_eff, subsume_rel_eff, subsume_max_csize;
		int		probe_inc, probe_min, probe_min_eff, probe_max_eff, probe_rel_eff;
		int		ternary_priorbins, ternary_min_eff, ternary_max_eff, ternary_rel_eff;
		int		transitive_min_eff, transitive_max_eff, transitive_rel_eff;
		int		vivify_priorbins, vivify_min_eff, vivify_max_eff, vivify_rel_eff;
		int		walk_priorbins, walk_min_eff, walk_max_eff, walk_rel_eff;
		//------------------------------------------//
		bool	report_en;
		bool	chrono_en;
		bool	stable_en;
		bool	reduce_en;
		bool	vivify_en;
		bool	subsume_en;
		bool	rephase_en;
		bool	bumpreason_en;
		bool	boundsearch_en;
		bool	decompose_en;
		bool	debinary_en;
		bool	transitive_en;
		bool	targetonly_en, targetphase_en;
		bool	ternary_en, ternary_sleep_en;
		bool	autarky_en, autarky_sleep_en;
		bool	proof_en, proof_nonbinary_en;
		bool	parseonly_en, parseincr_en;
		bool	vsids_en, vsidsonly_en;
		bool	probe_en, probehbr_en, probe_sleep_en;
		bool	model_en, modelprint_en, modelverify_en;
		bool	mdm_walk_en, mdm_mcv_en, mdmassume_en;
		//==========================================//
		//             Simplifier options           //
		//==========================================//
		bool	sub_en;
		bool	bce_en;
		bool	ere_en;
		bool	all_en;
		bool	solve_en;
		bool	profile_simp;
		bool	aggr_cnf_sort;
		bool	sigma_en, sigma_live_en, sigma_sleep_en;
		bool	ve_en, ve_plus_en, ve_lbound_en;
		//------------------------------------------//
		int		phases;
		int		shrink_rate;
		int		xor_max_arity;
		int		ve_clause_limit;
		int		sub_limit, bce_limit, ere_limit;
		int		ere_max_resolvent;
		//------------------------------------------//
		uint32	lcve_min;
		uint32	lits_min;
		uint32	mu_pos, mu_neg;
		//------------------------------------------//
		OPTION();
		~OPTION();
		void init();
	};

}

#endif

