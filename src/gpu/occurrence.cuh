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

	void resetOTAsync(CNF* cnf, OT* ot);
	void createOTAsync(CNF*, OT*, const bool&);
	void reduceOTAsync(CNF*, OT*, const bool&);

}

#endif