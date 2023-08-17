/***********************************************************************[primitives.cuh]
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

#ifndef __GPU_MODEL_
#define __GPU_MODEL_

#include "primitives.cuh"

namespace ParaFROST {

	#define MODEL_DBG 0

	_PFROST_D_ void	saveWitness(uint32*& saved, const uint32& witness)
	{
		ORIGINIZELIT(orgWitness, witness);
		*saved++ = orgWitness;
		*saved++ = 1;
	}

	_PFROST_D_ void	saveClause(uint32*& saved, SCLAUSE& c, const uint32& witlit)
	{
		uint32* first = saved, * witness = NULL;
		assert(c.original());
		forall_clause(c, k) {
			const uint32 lit = *k;
			if (lit == witlit) {
				witness = saved;
			}
			ORIGINIZELIT(orgLit, lit);
			*saved++ = orgLit;
		}
		assert(witness >= first);
		if (witness != first)
			devSwap(*first, *witness);
		else
			assert(*witness == *first);
		*saved++ = c.size();
	}

} 

#endif