/***********************************************************************[histogram.cu]
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

#include "solve.h"

using namespace pFROST;
using namespace SIGmA;

void ParaFROST::calcOccurs(const uint32& numLits)
{
	assert(numLits);
	uint32* literals = cumm.resizeLits(numLits);
	PFLOGN2(2, " Copying survived literals..");
	copyIf(literals, cnf, vars->gstats);
	assert(vars->gstats->numLits == numLits);
	PFLENDING(2, 5, "(%d copied)", numLits);
	histSimp(numLits);
}

void ParaFROST::histSimp(const uint32& numLits)
{
	PFLOGN2(2, " Computing histogram on %d elements..", numLits);
	assert(numLits);
	cuLits& culits = cumm.literals();
	assert(culits.size >= numLits);
	if (profile_gpu) cutimer->start();
	t_iptr& thrust_lits = culits.thrust_lits;
	t_iptr& thrust_hist = cuhist.thrust_hist;
	thrust::sort(thrust::cuda::par(tca), thrust_lits, thrust_lits + numLits);
	thrust::counting_iterator<size_t> search_begin(0);
	thrust::upper_bound(thrust::cuda::par(tca), thrust_lits, thrust_lits + numLits, search_begin, search_begin + inf.nDualVars, thrust_hist);
	thrust::adjacent_difference(thrust::cuda::par(tca), thrust_hist, thrust_hist + inf.nDualVars, thrust_hist);
	if (profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
	PFLDONE(2, 5);
}