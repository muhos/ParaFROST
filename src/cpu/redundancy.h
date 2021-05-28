/***********************************************************************[redundancy.h]
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

#ifndef __ERE_
#define __ERE_

#include "simplify.h" 
using namespace pFROST;

namespace SIGmA {

	#define MAX_ERE_OUT 190

	inline void forward_equ(Lits_t& m_c, OT& ot, const CL_ST& type)
	{
		pfrost->getStats().sigma.ere.tried++;
		int msize = m_c.size();
		assert(msize > 1);
		uint32 best = *m_c, m_sig = MAPHASH(best);
		assert(best > 1);
		int minsize = ot[best].size();
		for (int k = 1; k < msize; k++) {
			int lsize = ot[m_c[k]].size();
			if (lsize < minsize) minsize = lsize, best = m_c[k];
			m_sig |= MAPHASH(m_c[k]);
		}
		OL& minList = ot[best];
		for (int i = 0; i < minList.size(); i++) {
			S_REF c = minList[i];
			if (msize == c->size() && (c->learnt() || (c->status() == type)) &&
				sub(m_sig, c->sig()) && isEqual(*c, m_c)) {
				if (c->learnt()) pfrost->getStats().sigma.ere.learnts++;
				else pfrost->getStats().sigma.ere.orgs++;
				pfrost->removeClause(c);
				break;
			}
		}
	}

	inline bool merge_ere(const uint32& elim_var, const S_REF c1, const S_REF c2, Lits_t& out_c)
	{
		assert(elim_var);
		assert(!c1->deleted());
		assert(!c2->deleted());
		out_c.clear();
		int it1 = 0, it2 = 0;
		uint32 lit1, lit2, v1, v2;
		while (it1 < c1->size() && it2 < c2->size()) {
			lit1 = c1->lit(it1); lit2 = c2->lit(it2);
			v1 = ABS(lit1); v2 = ABS(lit2);
			if (v1 == elim_var) { it1++; }
			else if (v2 == elim_var) { it2++; }
			else if ((lit1 ^ lit2) == NEG_SIGN) return false; // tautology
			else if (v1 < v2) { it1++; out_c.push(lit1); }
			else if (v2 < v1) { it2++; out_c.push(lit2); }
			else { // repeated literal
				assert(lit1 == lit2);
				it1++; it2++;
				out_c.push(lit1);
			}
		}
		while (it1 < c1->size()) {
			lit1 = c1->lit(it1);
			if (ABS(lit1) == elim_var) it1++;
			else { it1++; out_c.push(lit1); }
		}
		while (it2 < c2->size()) {
			lit2 = c2->lit(it2);
			if (ABS(lit2) == elim_var) it2++;
			else { it2++; out_c.push(lit2); }
		}
		return true;
	}

}

#endif


