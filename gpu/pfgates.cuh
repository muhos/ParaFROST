/***********************************************************************[pfgates.cuh]
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

#ifndef __SIGMA_GATES_
#define __SIGMA_GATES_

#include "pfdevice.cuh"

namespace pFROST {

	namespace SIGmA {

#define COUNTFLIPS(X) \
	{ \
		while (__popc(++X) & 1); \
	}

		_PFROST_D_ void toblivion(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			bool which = pOrgs > nOrgs;
			const uint32 n = NEG(p);
			if (which) {
				uint32 nsLits = 0;
				countLitsBefore(cnf, negs, nsLits);
				uint32 nElems = nOrgs + nsLits + 2;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (S_REF* i = negs; i != negs.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, n);
					c.markDeleted();
				}
				saveResolved(saved, p);
			}
			else {
				uint32 psLits = 0;
				countLitsBefore(cnf, poss, psLits);
				uint32 nElems = pOrgs + psLits + 2;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (S_REF* i = poss; i != poss.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, p);
					c.markDeleted();
				}
				saveResolved(saved, n);
			}
			OL& other = which ? poss : negs;
#pragma unroll
			for (S_REF* i = other; i != other.end(); i++) cnf[*i].markDeleted();
			poss.clear(true), negs.clear(true);
		}

		_PFROST_D_ void saveResolved(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			bool which = pOrgs > nOrgs;
			const uint32 n = NEG(p);
			if (which) {
				uint32 nsLits = 0;
				countLitsBefore(cnf, negs, nsLits);
				uint32 nElems = nOrgs + nsLits + 2;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (S_REF* i = negs; i != negs.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, n);
				}
				saveResolved(saved, p);
			}
			else {
				uint32 psLits = 0;
				countLitsBefore(cnf, poss, psLits);
				uint32 nElems = pOrgs + psLits + 2;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (S_REF* i = poss; i != poss.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, p);
				}
				saveResolved(saved, n);
			}
		}

		_PFROST_D_ void countSubstituted(const uint32& x, CNF& cnf, OL& me, OL& other, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(x);
			assert(!nAddedCls);
			assert(!nAddedLits);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() != cj.molten() && !isTautology(x, ci, cj))
						nAddedCls++, nAddedLits += ci.size() + cj.size() - 2;
				}
			}
		}

		_PFROST_D_ void countSubstituted(const uint32& x, CNF& cnf, OL& me, OL& other, uint32& nAddedCls)
		{
			assert(x);
			assert(!nAddedCls);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() != cj.molten() && !isTautology(x, ci, cj))
						nAddedCls++;
				}
			}
		}

		_PFROST_D_ void substitute_single(const uint32& dx, const uint32& def, SCLAUSE& org, cuVecU* units)
		{
			assert(dx > 1);
			assert(def != dx);
			assert(org.original());
#if VE_DBG
			printf("c | Clause ");
			org.print();
#endif
			int n = 0;
#pragma unroll
			for (int i = 0; i < org.size(); i++)
				org[n++] = org[i] == dx ? def : org[i];
			org.resize(n);
			devSort(org.data(), org.size());
			calcSig(org);
#if VE_DBG
			printf("c | Substituted to ");
			org.print();
#endif
			if (org.size() == 1) units->push(*org);
		}

		_PFROST_D_ void substitute_single(const uint32& p, const uint32& def, CNF& cnf, OL& poss, OL& negs, cuVecU* units)
		{
			assert(def > 1);
			assert(!SIGN(p));
			uint32 n = NEG(p), def_f = FLIP(def);
			// substitute negatives 
#pragma unroll
			for (S_REF* i = negs; i != negs.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.learnt() || c.has(def)) c.markDeleted(); // learnt or tautology
				else substitute_single(n, def_f, c, units);
			}
			// substitute positives
#pragma unroll
			for (S_REF* i = poss; i != poss.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.learnt() || c.has(def_f)) c.markDeleted(); // learnt or tautology
				else substitute_single(p, def, c, units);
			}
		}

		_PFROST_D_ void freeze_binaries(CNF& cnf, OL& list)
		{
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.size() == 2) c.freeze();
			}
		}

		_PFROST_D_ void freeze_arities(CNF& cnf, OL& list)
		{
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.size() > 2 && c.molten()) c.freeze();
			}
		}

		_PFROST_D_ int find_fanin(const uint32& gate_out, CNF& cnf, OL& list, uint32* out_c, uint32& sig)
		{
			assert(gate_out > 1);
			sig = 0;
			uint32 imp = 0;
			int nImps = 0;
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.learnt()) continue;
				assert(!c.molten());
				if (c.size() == 2) {
					imp = FLIP(c[0] ^ c[1] ^ gate_out);
					out_c[nImps++] = imp;
					sig |= MAPHASH(imp);
					c.melt(); // mark as gate clause
				}
			}
			return nImps;
		}

		_PFROST_D_ uint32 find_sfanin(const uint32& gate_out, CNF& cnf, OL& list)
		{
			assert(gate_out > 1);
			uint32 imp = 0;
			int nImps = 0;
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.size() == 2) {
					imp = FLIP(c[0] ^ c[1] ^ gate_out);
					nImps++;
				}
				if (nImps > 1) return 0; // cannot be a single-input gate
			}
			return imp;
		}

		_PFROST_D_ S_REF fast_equality_check(CNF& cnf, OT& ot, uint32 x, uint32 y, uint32 z) {
			if (ot[y].size() > ot[z].size()) devSwap(y, z);
			if (ot[x].size() > ot[y].size()) devSwap(x, y);
			OL& list = ot[x];
			sort3(x, y, z, DEFAULT_CMP<uint32>());
			assert(x <= y && y <= z && x <= z);
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.molten()) continue;
				assert(c.isSorted());
				if (c.original() && c.size() == 3 &&
					c[0] == x && c[1] == y && c[2] == z) return *i;
			}
			return GNOREF;
		}

		_PFROST_D_ void freeze_arities(CNF& cnf, OL& me, OL& other)
		{
			freeze_arities(cnf, me);
			freeze_arities(cnf, other);
		}

		_PFROST_D_ void shareXORClause(SCLAUSE& c, uint32* shared)
		{

			uint32* lit = c, * cend = c.end();
#pragma unroll
			while (lit != cend)
				*shared++ = *lit++;
		}

		_PFROST_D_ bool checkArity(SCLAUSE& c, uint32* literals, const int& size)
		{
			assert(literals);
			const uint32* end = literals + size;
			uint32* i = c, * cend = c.end();
#pragma unroll
			while (i != cend) {
				const uint32 lit = *i++;
				uint32* j;
				for (j = literals; j != end; j++) {
					if (lit == *j) {
						break;
					}
				}
				if (j == end) return false;
			}
			return true;
		}

		_PFROST_D_ bool makeArity(CNF& cnf, OT& ot, uint32& parity, uint32* literals, const int& size)
		{
			const uint32 oldparity = parity;
			COUNTFLIPS(parity);
			for (int k = 0; k < size; k++) {
				const uint32 bit = (1UL << k);
				if ((parity & bit) ^ (oldparity & bit))
					literals[k] = FLIP(literals[k]);
			}
			// search for an arity clause
			assert(size > 2);
			uint32 best = *literals;
			assert(best > 1);
			int minsize = ot[best].size();
			for (int k = 1; k < size; k++) {
				const uint32 lit = literals[k];
				assert(lit > 1);
				int lsize = ot[lit].size();
				if (lsize < minsize) {
					minsize = lsize;
					best = literals[k];
				}
			}
			OL& minlist = ot[best];
			for (S_REF* i = minlist; i != minlist.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.size() == size && checkArity(c, literals, size)) {
					assert(c.original());
					c.melt();
					return true;
				}
			}
			return false;
		}

	}
}


#endif