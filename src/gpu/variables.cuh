/***********************************************************************[variables.cuh]
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

#ifndef __GPU_VARS_
#define __GPU_VARS_

#include "vector.cuh"

namespace ParaFROST {

	struct VARS {
		cuVecU* elected,	 * units,       * resolved;
		uint32* electedData, * electedSize;
		uint32* unitsData, * unitsSize;
		uint32* scores, * varcore, * eligible, * cachedUnits;
		Byte*	eliminated,  * cachedEliminated;
		cuVecU  tmpUnits;
		uint32  numElected, nUnits, currMelted;
		bool	isEliminatedCached;
		VARS() :
			elected(NULL)
			, units(NULL)
			, resolved(NULL)
			, electedData(NULL)
			, electedSize(NULL)
			, unitsData(NULL)
			, unitsSize(NULL)
			, scores(NULL)
			, varcore(NULL)
			, eligible(NULL)
			, cachedUnits(NULL)
			, eliminated(NULL)
			, cachedEliminated(NULL)
			, numElected(0), nUnits(0), currMelted(0)
			, isEliminatedCached(false)
		{ }
	};

}

#endif
