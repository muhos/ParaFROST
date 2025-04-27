/***********************************************************************[occurrence.cu]
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

#include "solver.hpp"
#include "options.cuh"
#include "timer.cuh"
#include "grid.cuh"

namespace ParaFROST {

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

	void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
			LASTERR("Occurrence table reduction failed");
			SYNCALL;
			if (p) {
				LOGRULER('=', 30);
				LOG0("\toccurrence table");
				ot->print();
				LOGRULER('=', 30);
			}
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->rot += cutimer->gpuTime();
	}

	void resetOTAsync(CNF* cnf, OT* ot)
	{
		assert(cnf);
		assert(ot);
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reset_ot_k << <nBlocks, BLOCK1D >> > (ot);
		if (gopts.sync_always) {
			LASTERR("Occurrence table reset failed");
			SYNCALL;
			assert(ot->accViolation(inf.maxVar));
		}
	}

	void createOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		resetOTAsync(cnf, ot);
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		create_ot_k << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
			LASTERR("Occurrence table creation failed");
			SYNCALL;
			assert(ot->accViolation(inf.maxVar));
			if (p) {
				LOGRULER('=', 30);
				LOG0("\toccurrence table");
				ot->print();
				LOGRULER('=', 30);
			}
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->cot += cutimer->gpuTime();
	}


	bool Solver::reallocOT(const cudaStream_t& stream)
	{
		assert(inf.nLiterals);
		if (!flattenCNF(inf.nLiterals)) { simpstate = OTALLOC_FAIL; return false; }
		histSimp(inf.nLiterals);
		if (!cumm.resizeOTAsync(ot, inf.nLiterals, stream)) { simpstate = OTALLOC_FAIL; return false; }
		return true;
	}

}