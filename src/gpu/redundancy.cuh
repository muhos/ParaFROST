/***********************************************************************[redundancy.cuh]
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

#ifndef __SIGMA_ERE_
#define __SIGMA_ERE_

#include "device.cuh"

namespace pFROST {

	namespace SIGmA {

		_PFROST_D_ int merge_ere(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(!c1.deleted());
			assert(!c2.deleted());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1++];
				if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
			}
			while (it2 < n2) {
				lit2 = c2[it2++];
				if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
			}
			assert(len <= SH_MAX_ERE_OUT);
			return len;
		}

		_PFROST_D_ void forward_equ(CNF& cnf, OT& ot, uint32* m_c, const int& m_len, const CL_ST& type)
		{
			assert(m_len > 1);
			assert(type != DELETED);
			uint32 best = *m_c, m_sig = MAPHASH(best);
			assert(best > 1);
			int minsize = ot[best].size();
			for (int k = 1; k < m_len; k++) {
				int lsize = ot[m_c[k]].size();
				if (lsize < minsize) minsize = lsize, best = m_c[k];
				m_sig |= MAPHASH(m_c[k]);
			}
			OL& minList = ot[best];
			for (S_REF* i = threadIdx.x + minList; i < minList.end(); i += blockDim.x) {
				SCLAUSE& c = cnf[*i];
				if (m_len == c.size() && (c.learnt() || (c.status() == type)) &&
					sub(m_sig, c.sig()) && isEqual(c, m_c, m_len)) {
					c.markDeleted();
					break;
				}
			}
		}

	} // sigma namespace
} // parafrost namespace


#endif