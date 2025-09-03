/***********************************************************************[blocked.cuh]
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

#ifndef __BLOCKED_
#define __BLOCKED_

#include "elimination.cuh"

namespace ParaFROST {

_PFROST_D_ void blocked_x(
    const uint32& x,
    CNF& cnf,
    OL& poss,
    OL& negs,
    cuVecU* resolved,
    cuVecB* proof,
    uint32* sh_c) {
#pragma unroll
    forall_occurs(negs, i) { // start with negs
        SCLAUSE& ci = cnf[*i];
        if (ci.deleted() || ci.learnt()) continue;
        bool allTautology = true;
        int c_size = ci.size();
        if (c_size <= SH_MAX_BCE_IN) { // use shared memory
            ci.shareTo(sh_c);
#pragma unroll
            forall_occurs(poss, j) { // block with poss
                SCLAUSE& cj = cnf[*j];
                if (cj.deleted() || cj.learnt()) continue;
                if (!isTautology(x, cj, sh_c, c_size)) {
                    allTautology = false;
                    break;
                }
            }
        } else { // use global memory
#pragma unroll
            forall_occurs(poss, j) { // block with poss
                SCLAUSE& cj = cnf[*j];
                if (cj.deleted() || cj.learnt()) continue;
                if (!isTautology(x, ci, cj)) {
                    allTautology = false;
                    break;
                }
            }
        }
        if (allTautology) {
            const int size = ci.size() + 1;
            assert(size > 2);
            uint32* saved = resolved->jump(size);
            saveClause(saved, ci, NEG(V2L(x)));
            if (proof) {
                uint32 bytes = 0;
                countProofBytes(ci, bytes);
                addr_t pdata = proof->jump(bytes);
                saveProofClause(pdata, ci, PROOF_DELETED);
            }
            ci.markDeleted();
        }
    }
}

// kernel
__global__ void bce_k(
    CNF* __restrict__ cnf,
    OT* __restrict__ ot,
    cuVecB* __restrict__ proof,
    cuVecU* __restrict__ resolved,
    const cuVecU* __restrict__ pVars,
    const Byte* __restrict__ eliminated) {
    __shared__ uint32 sh_cls[BLBCE * SH_MAX_BCE_IN];
    for_parallel_x(tid, pVars->size()) {
        const uint32 x = (*pVars)[tid];
        assert(x);
        assert(!ELIMINATED(eliminated[x]));
        const uint32 p = V2L(x), n = NEG(p);
        if ((*ot)[p].size() <= kOpts->bce_max_occurs && (*ot)[n].size() <= kOpts->bce_max_occurs)
            blocked_x(x, *cnf, (*ot)[p], (*ot)[n], resolved, proof, sh_cls + threadIdx.x * SH_MAX_BCE_IN);
    }
}

} // namespace ParaFROST

#endif