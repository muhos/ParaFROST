/***********************************************************************[occurrence.cuh]
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

#ifndef __OCCURRENCE_
#define __OCCURRENCE_

#include "cnf.cuh"
#include "table.cuh"
#include "grid.cuh"

namespace ParaFROST {

	__global__ void assign_scores(
		uint32* __restrict__ eligible,
		uint32* __restrict__ scores,
		const uint32* __restrict__ hist,
		uint32 size)
	{
		grid_t tid = global_tx;
		while (tid < size) {
			const uint32 v = tid + 1;
			const uint32 p = V2L(v), ps = hist[p], ns = hist[NEG(p)];
			eligible[tid] = v;
			scores[v] = ps * ns;
			tid += stride_x;
		}
	}

	__global__ void assign_scores(
		uint32* __restrict__ eligible,
		uint32* __restrict__ scores,
		uint32* __restrict__ hist,
		const OT* __restrict__ ot,
		uint32 size)
	{
		grid_t tid = global_tx;
		while (tid < size) {
			const uint32 v = tid + 1;
			const uint32 p = V2L(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
			hist[p] = ps, hist[n] = ns;
			eligible[tid] = v;
			scores[v] = ps * ns;
			tid += stride_x;
		}
	}

	__global__ void reduce_ot(const CNF* __restrict__ cnfptr, OT* __restrict__ ot)
	{
		grid_t tid = global_tx;
		while (tid < ot->size()) {
			OL& ol = (*ot)[tid];
			if (ol.size()) {
				const CNF& cnf = *cnfptr;
				S_REF* j = ol;
				forall_occurs(ol, i) {
					const S_REF ref = *i;
					if (!cnf[ref].deleted())
						*j++ = ref;
				}
				ol.resize(j - ol);
			}
			tid += stride_x;
		}
	}

	__global__ void reset_ot_k(OT* ot)
	{
		grid_t tid = global_tx;
		while (tid < ot->size()) {
			(*ot)[tid].clear();
			tid += stride_x;
		}
	}

	__global__ void create_ot_k(CNF* __restrict__ cnf, OT* __restrict__ ot_ptr)
	{
		grid_t tid = global_tx;
		while (tid < cnf->size()) {
			const S_REF r = cnf->ref(tid);
			SCLAUSE& c = (*cnf)[r];
			if (c.original() || c.learnt()) {
				OT& ot = *ot_ptr;
				forall_clause(c, lit) {
					ot[*lit].insert(r);
				}
			}
			tid += stride_x;
		}
	}

}

#endif