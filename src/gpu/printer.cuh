/***********************************************************************[printer.cuh]
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

#ifndef __CU_PRINTER_
#define __CU_PRINTER_

#include "simptypes.cuh"

namespace ParaFROST {

	_PFROST_H_D_ void pLit(const uint32& l) { printf("%6d", SIGN(l) ? -int(ABS(l)) : ABS(l)); }

	_PFROST_H_D_ void pSharedClause(const uint32* c, const int& sz)
	{
		printf("(");
		for (int k = 0; k < sz; k++) pLit(c[k]), printf(" ");
		printf(") shared\n");
	}

	_PFROST_H_D_ void pClauseSet(const CNF& cnf, const OT& ot, const uint32& v)
	{
		const OL& poss = ot[V2L(v)], &negs = ot[NEG(V2L(v))];
		for (uint32 i = 0; i < poss.size(); i++) {
			printf("c  "), cnf[poss[i]].print();
		}
		for (uint32 i = 0; i < negs.size(); i++) {
			printf("c  "), cnf[negs[i]].print();
		}
	}

	_PFROST_H_D_ void pClauseSet(const CNF& cnf, const OL& poss, const OL& negs)
	{
		for (uint32 i = 0; i < poss.size(); i++) {
			printf("c  "), cnf[poss[i]].print();
		}
		for (uint32 i = 0; i < negs.size(); i++) {
			printf("c  "), cnf[negs[i]].print();
		}
	}

	_PFROST_H_D_ void pClauseSet(const CNF& cnf, const OL& ol)
	{
		for (uint32 i = 0; i < ol.size(); i++) {
			printf("c  c(%lld)->", uint64(ol[i]));
			cnf[ol[i]].print();
		}
	}

	_PFROST_H_D_ void printGate(const CNF& cnf, const OL& poss, const OL& negs)
	{
		for (uint32 i = 0; i < poss.size(); i++) {
			if (cnf[poss[i]].molten()) {
				printf("c  ");
				cnf[poss[i]].print();
			}
		}
		for (uint32 i = 0; i < negs.size(); i++) {
			if (cnf[negs[i]].molten()) {
				printf("c  ");
				cnf[negs[i]].print();
			}
		}
	}

}

#endif