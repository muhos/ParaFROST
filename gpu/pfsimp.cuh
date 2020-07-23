/***********************************************************************[pfsimp.cuh]
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
#include <thrust/system/cuda/execution_policy.h>
#include "pfsimptypes.h"
#include "pfcuconst.h"

namespace pFROST {

	namespace SIGmA {

		// Special GPU comparators
		struct GPU_LCV_CMP {
			uint32* scores;
			_PFROST_H_D_ GPU_LCV_CMP(uint32* _scores) : scores(_scores) { assert(scores != NULL); }
			_PFROST_D_ bool operator () (uint32& x, uint32& y) {
				return (scores[x] != scores[y]) ? scores[x] < scores[y] : x < y;
			}
		};
		struct GPU_MCV_CMP {
			uint32* scores;
			_PFROST_H_D_ GPU_MCV_CMP(uint32* _scores) : scores(_scores) { assert(scores != NULL); }
			_PFROST_D_ bool operator () (uint32& x, uint32& y) {
				return (scores[x] != scores[y]) ? scores[x] > scores[y] : x > y;
			}
		};
		struct COMPACT_CMP
		{
			CNF* cnf;
			_PFROST_H_D_ COMPACT_CMP(CNF* _cnf) : cnf(_cnf) {}
			_PFROST_D_ bool operator()(const S_REF& ref) {
				return !cnf->cref(ref)->deleted();
			}
		};
		struct CNF_CMP_SZ {
			CNF* cnf;
			_PFROST_H_D_ CNF_CMP_SZ(CNF* _cnf) : cnf(_cnf) { assert(cnf != NULL); }
			_PFROST_D_ bool operator () (S_REF& a, S_REF& b) {
				SCLAUSE& x = (*cnf)[a], &y = (*cnf)[b];
				uint32 xsize = x.size();
				if (xsize != y.size()) return x.size() < y.size();
				else if (*x != *y) return *x < *y;
				else if (x[1] != y[1]) return x[1] < y[1];
				else if (xsize > 2 && x.back() != y.back()) return x.back() < y.back();
				else return x.sig() < y.sig();
			}
		};
		//======================================================//
		//                GPU Wrappers Declaration              //
		//======================================================//
		void cuMemSetAsync(addr_t, const Byte&, const size_t&);
		void copyIf(uint32*, CNF*, GSTATS*);
		void calcScores(VARS*, uint32*);
		void calcScores(VARS*, uint32*, OT*);
		void calcSigCNFAsync(CNF*, const uint32&, const uint32&, const cudaStream_t&);
		void calcSigCNF(CNF*, const uint32&);
		void createOTAsync(CNF*, OT*, const bool&);
		void reduceOTAsync(CNF*, OT*, const bool&);
		void sortOTAsync(CNF*, OT*, VARS* vars, cudaStream_t*);
		void veAsync(CNF*, OT*, VARS*, cudaStream_t*, const cuHist&, const uint32&, const uint32&, const int&, const bool&);
		void hseAsync(CNF*, OT*, VARS*, const uint32&);
		void bceAsync(CNF*, OT*, VARS*, const uint32&);
		void hreAsync(CNF*, OT*, VARS*, const uint32&);
		void evalReds(CNF*, GSTATS*, LIT_ST*);
		void countFinal(CNF*, GSTATS*, LIT_ST*);
		void countCls(CNF*, GSTATS*);
		void countAll(CNF*, GSTATS*);
		void countLits(CNF*, GSTATS*);
		void countMelted(LIT_ST*);

	}
}

#endif 
