/***********************************************************************[options.cuh]
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

#ifndef __GPU_OPTIONS_
#define __GPU_OPTIONS_

#include "datatypes.hpp"
#include "constants.cuh"
#include "definitions.cuh"

namespace ParaFROST {

	struct KOptions {
		bool ve_fun_en : 1;
		bool ve_lbound_en : 1;
		bool ere_extend_en : 1;
		bool ere_extend_all_en : 1;
		bool proof_en : 1;

		int	sub_clause_max;
		int	ere_clause_max;

		uint32 ve_clause_max;
		uint32 xor_max_arity;
		uint32 sub_max_occurs;
		uint32 ere_max_occurs;
		uint32 bce_max_occurs;
	};

	struct GOPTION {
		GOPTION();

		void init(const bool& opt_proof_en);

		bool unified_access;
		bool profile_gpu;
		bool sync_always;

		uint32 ve_min_threads;
		uint32 sub_min_threads;
		uint32 ere_min_threads;

		double ve_min_blocks;
		double sub_min_blocks;
		double ere_min_blocks;

		KOptions hostKOpts;
	};

	extern GOPTION gopts;

	extern __constant__ 
	KOptions kOpts[1];

	void initDevOpts();

}

#endif