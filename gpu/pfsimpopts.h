/***********************************************************************[pfsimpopts.cpp]
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
#define MAX_GPU_COUNT (31 + 1) // index starts from GPU 1 (first slave)
#define MAX_STREAMS 32

BOOL_OPT opt_pcnf_en("print-cnf", "print GPU-resident CNF on stdout", false);
BOOL_OPT opt_pot_en("print-ot", "print GPU-resident occurrence table on stdout", false);
BOOL_OPT opt_ve_en("ve", "enable bounded variable elimination (BVE)", true);
BOOL_OPT opt_sub_en("sub", "enable hybrid subsumption elimination (HSE) with high bounds in LCVE", false);
BOOL_OPT opt_ve_plus_en("ve+", "enable (BVE + HSE) untill no literals can be removed", true);
BOOL_OPT opt_bce_en("bce", "enable blocked clause elimination", false);
BOOL_OPT opt_hre_en("hre", "enable hidden redundancy elimination", false);
BOOL_OPT opt_all_en("all", "enable all simplifications", false);
INT_OPT opt_mu_pos("mu-pos", "set the positive freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_mu_neg("mu-neg", "set the negative freezing temperature in LCVE", 32, INT32R(10, INT32_MAX));
INT_OPT opt_phases("phases", "set the number of phases in stage-1 reductions", 2, INT32R(0, INT32_MAX));
INT_OPT opt_cnf_free("cnf-free-freq", "set the frequency of CNF memory shrinkage in SIGmA", 3, INT32R(0, 5));
INT_OPT opt_gpus("ngpus", "number of GPU(s) to be activated", 1, INT32R(1, MAX_GPU_COUNT));
INT_OPT opt_streams("nstreams", "number of GPU streams to be created", 4, INT32R(1, MAX_STREAMS));

#endif