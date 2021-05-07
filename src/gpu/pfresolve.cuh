/***********************************************************************[pfresolve.cuh]
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

#ifndef __SIGMA_RESOLVE_
#define __SIGMA_RESOLVE_

#include "pfdevice.cuh"

namespace pFROST {

	namespace SIGmA {

		_PFROST_D_ bool resolve(const uint32& x, const uint32& nOrgCls, CNF& cnf, OL& poss, OL& negs, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(x);
			assert(checkMolten(cnf, poss, negs));
			// check resolvability
			nAddedCls = 0, nAddedLits = 0;
			if (countResolvents(x, nOrgCls, cnf, poss, negs, nAddedCls, nAddedLits)) return false;
#if VE_DBG
			printf("c  Resolving(%d) ==> added = %d, deleted = %d\n", x, nAddedCls, poss.size() + negs.size());
			pClauseSet(cnf, poss, negs);
#endif		
			return true; 
		}

	} // sigma namespace
} // parafrost namespace


#endif