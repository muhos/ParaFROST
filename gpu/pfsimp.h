/***********************************************************************[pfsimp.h]
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

#ifndef __SIGMA_SIMP_
#define __SIGMA_SIMP_

#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/counting_iterator.h>
#include "pfsimptypes.h"
#include "pfcuconst.h"

namespace pFROST {

	namespace SIGmA {
		//======================================================//
		//                GPU Wrappers Declaration              //
		//======================================================//
		void cuMemSetAsync(addr_t, const Byte&, const size_t&);
		void copyIf(uint32*, CNF*, GSTATS*);
		void calcScores(VARS*, uint32*);
		void calcScores(VARS*, uint32*, OT*);
		void calcSigCNFAsync(CNF*, const uint32&, const uint32&, const cudaStream_t&);
		void calcSigCNF(CNF*, const uint32&);
		void createOT(CNF*, OT*, const bool&);
		void createOTAsync(CNF*, OT*, const bool&);
		void reduceOTAsync(CNF*, OT*, const bool&);
		void reduceOT(CNF*, OT*, VARS*, cudaStream_t*, const bool&);
		void ve(CNF*, OT*, VARS*);
		void hse(CNF*, OT*, VARS*);
		void bce(CNF*, OT*, VARS*);
		void hre(CNF*, OT*, VARS*);
		void evalReds(CNF*, VARS*, cudaStream_t*);
		void countCls(CNF*, GSTATS*);
		void countLits(CNF*, GSTATS*);
		void countDelVars(VARS* vars);

	}
}

#endif 
