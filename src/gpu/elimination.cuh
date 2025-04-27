/***********************************************************************[elimination.cuh]
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

#ifndef __ELIMINATION_
#define __ELIMINATION_

#include "cnf.cuh"
#include "table.cuh"
#include "model.cuh"
#include "options.cuh"
#include "printer.cuh"
#include "proofutils.cuh"
#include "primitives.cuh"
#include "definitions.hpp"

namespace ParaFROST {

	#define RES_UNBOUNDED 0

	#define SUBSIG(A,B) (!((A) & ~(B)))

	_PFROST_H_D_ bool checkMolten(const CNF& cnf, const OL& poss, const OL& negs)
	{
		for (uint32 i = 0; i < poss.size(); i++)
			if (cnf[poss[i]].molten()) return false;
		for (uint32 i = 0; i < negs.size(); i++)
			if (cnf[negs[i]].molten()) return false;
		return true;
	}

	_PFROST_H_D_ bool checkDeleted(const CNF& cnf, const OL& poss, const OL& negs)
	{
		for (uint32 i = 0; i < poss.size(); i++)
			if (cnf[poss[i]].deleted()) return false;
		for (uint32 i = 0; i < negs.size(); i++)
			if (cnf[negs[i]].deleted()) return false;
		return true;
	}

	_PFROST_D_ bool isEqual(SCLAUSE& c1, const uint32* c2, const int& size)
	{
		assert(!c1.deleted());
		assert(size > 1);
		assert(c1.size() == size);
		const uint32* c1_data = c1.data();
		for (int it = 0; it < size; it++) 
			if (NEQUAL(c1_data[it], c2[it])) 
				return false;
		return true;
	}

	_PFROST_D_ bool selfsub(const uint32& A, const uint32& B)
	{
		uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
		return !(A & ~B_tmp);
	}

	_PFROST_D_ void freezeBinaries(CNF& cnf, OL& list)
	{
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original() && c.size() == 2) c.freeze();
		}
	}

	_PFROST_D_ void freezeClauses(CNF& cnf, OL& poss, OL& negs)
	{
		forall_occurs(poss, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original() && c.molten())
				c.freeze();
		}
		forall_occurs(negs, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original() && c.molten())
				c.freeze();
		}
	}

	_PFROST_D_ void appendUnits(const CNF& cnf, OL& ol, uint32*& units)
	{
		assert(units);
		forall_occurs(ol, i) {
			const SCLAUSE& c = cnf[*i];
			if (c.size() == 1) {
				assert(!c.deleted());
				*units++ = c[0];
			}
		}
	}

	_PFROST_D_ void reduceOL(const CNF& cnf, OL& ol)
	{
		if (ol.empty()) return;
		S_REF* j = ol;	
		forall_occurs(ol, i) {
			const S_REF ref = *i;
			if (!cnf[ref].deleted())
				*j++ = ref;
		}
		ol.resize(j - ol);
	}

	_PFROST_D_ bool isTautology(const uint32& x, SCLAUSE& c1, SCLAUSE& c2)
	{
		assert(x);
		assert(!c1.deleted());
		assert(!c2.deleted());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		uint32* n1 = c1.end(), * n2 = c2.end();
		uint32* lit1 = c1, * lit2 = c2, v1, v2;
		while (lit1 != n1 && lit2 != n2) {
			v1 = ABS(*lit1), v2 = ABS(*lit2);
			if (v1 == x) lit1++;
			else if (v2 == x) lit2++;
			else if (IS_TAUTOLOGY(*lit1, *lit2)) return true;
			else if (v1 < v2) lit1++;
			else if (v2 < v1) lit2++;
			else { lit1++, lit2++; }
		}
		return false;
	}

	_PFROST_D_ bool isTautology(const uint32& x, const SCLAUSE& c1, const uint32* c2, const int& n2)
	{
		assert(x);
		assert(!c1.deleted());
		assert(c1.size() > 1);
		assert(n2 > 1);
		const int n1 = c1.size();
		int it1 = 0, it2 = 0;
		uint32 lit1, lit2, v1, v2;
		while (it1 < n1 && it2 < n2) {
			lit1 = c1[it1], lit2 = c2[it2];
			v1 = ABS(lit1), v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if (IS_TAUTOLOGY(lit1, lit2)) return true;
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { it1++; it2++; }
		}
		return false;
	}

	_PFROST_D_ int mergeProof(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, uint32& bytes)
	{
		assert(kOpts->proof_en);
		assert(x);
		assert(c1.original());
		assert(c2.original());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		assert(DC_PTRS->d_lbyte);
		const addr_t lbyte = DC_PTRS->d_lbyte;
		uint32 local = 2; // prefix + suffix
		const int n1 = c1.size(), n2 = c2.size();
		int it1 = 0, it2 = 0;
		int len = n1 + n2 - 2;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
			else if (v1 < v2) { 
				it1++;
				local += lbyte[lit1]; 
			}
			else if (v2 < v1) { 
				it2++;
				local += lbyte[lit2];
			}
			else { // repeated literal
				it1++, it2++;
				assert(len > 1);
				len--;
				local += lbyte[lit1];
			}
		}
		assert(len > 0);
		while (it1 < n1) {
			const uint32 lit1 = c1[it1++];
			if (NEQUAL(ABS(lit1), x)) {
				local += lbyte[lit1];
			}
		}
		while (it2 < n2) {
			const uint32 lit2 = c2[it2++];
			if (NEQUAL(ABS(lit2), x)) {
				local += lbyte[lit2];
			}
		}
		assert(len > 0);
		assert(local > 2);
		bytes += local; 
		return len;
	}

	_PFROST_D_ int merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2)
	{
		assert(x);
		assert(c1.original());
		assert(c2.original());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		const int n1 = c1.size(), n2 = c2.size();
		int it1 = 0, it2 = 0;
		int len = n1 + n2 - 2;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { // repeated literal
				it1++, it2++;
				assert(len > 1);
				len--;
			}
		}
		assert(len > 0);
		return len;
	}

	_PFROST_D_ int merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, uint32& unit)
	{
		assert(x);
		assert(c1.original());
		assert(c2.original());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		const int n1 = c1.size(), n2 = c2.size();
		int it1 = 0, it2 = 0;
		int len = n1 + n2 - 2;
		unit = 0;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { // repeated literal
				it1++, it2++;
				assert(len > 1);
				len--;
				if (len == 1) unit = lit1;
			}
		}
		assert(len > 0);
		return len;
	}

	_PFROST_D_ void merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, SCLAUSE* out_c)
	{
		assert(x);
		assert(c1.original());
		assert(c2.original());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		assert(out_c->empty());
		const int n1 = c1.size(), n2 = c2.size();
		int it1 = 0, it2 = 0;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if (v1 < v2) { it1++, out_c->push(lit1); }
			else if (v2 < v1) { it2++, out_c->push(lit2); }
			else { // repeated literal
				it1++, it2++;
				out_c->push(lit1);
			}
		}
		while (it1 < n1) {
			const uint32 lit1 = c1[it1++];
			if (NEQUAL(ABS(lit1), x)) out_c->push(lit1);
		}
		while (it2 < n2) {
			const uint32 lit2 = c2[it2++];
			if (NEQUAL(ABS(lit2), x)) out_c->push(lit2);
		}
		calcSig(*out_c);
		assert(out_c->isSorted());
		assert(out_c->hasZero() < 0);
	}

	_PFROST_D_ int merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, uint32* out_c)
	{
		assert(x);
		assert(c1.original());
		assert(c2.original());
		assert(c1.size() > 1);
		assert(c2.size() > 1);
		const int n1 = c1.size(), n2 = c2.size();
		int it1 = 0, it2 = 0;
		int len = 0;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
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
			const uint32 lit1 = c1[it1++];
			if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
		}
		while (it2 < n2) {
			const uint32 lit2 = c2[it2++];
			if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
		}
		assert(len > 0);
		return len;
	}

	_PFROST_D_ int merge(const uint32& x, const uint32* c1, const int& n1, const SCLAUSE& c2, uint32* out_c)
	{
		assert(x);
		assert(n1 > 1);
		assert(c2.original());
		assert(c2.size() > 1);
		const int n2 = c2.size();
		int it1 = 0, it2 = 0;
		int len = 0;
		while (it1 < n1 && it2 < n2) {
			const uint32 lit1 = c1[it1], lit2 = c2[it2];
			const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
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
			const uint32 lit1 = c1[it1++];
			if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
		}
		while (it2 < n2) {
			const uint32 lit2 = c2[it2++];
			if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
		}
		assert(len > 0);
		return len;
	}

	_PFROST_D_ void countOrgs(CNF& cnf, OL& list, uint32& orgs)
	{
		assert(!orgs);
		forall_occurs(list, i) {
			if (cnf[*i].original()) 
				orgs++;
		}
	}

	_PFROST_D_ void countOrgs(CNF& cnf, OL& list, uint32& nClsBefore, uint32& nLitsBefore)
	{
		assert(!nClsBefore);
		assert(!nLitsBefore);
		forall_occurs(list, i) {
			const SCLAUSE& c = cnf[*i];
			if (c.original()) {
				nClsBefore++;
				nLitsBefore += c.size();
			}
		}
	}

	_PFROST_D_ void countLitsBefore(CNF& cnf, OL& list, uint32& nLitsBefore)
	{
		forall_occurs(list, i) {
			const SCLAUSE& c = cnf[*i];
			if (c.original()) 
				nLitsBefore += c.size();
		}
	}

	_PFROST_D_ bool countSubstituted(
		const uint32& x,
		const uint32& nClsBefore, 
		CNF& cnf,
		OL& me,
		OL& other, 
		uint32& nElements,
		uint32& nAddedCls, 
		uint32& nAddedLits)
	{
		assert(x);
		assert(!nElements);
		assert(!nAddedCls);
		assert(!nAddedLits);
		const int rlimit = kOpts->ve_clause_max;
		// check if proof bytes has to be calculated
		if (kOpts->proof_en) {
			uint32 proofBytes = 0;
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				const bool ci_m = ci.molten();
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.original() && NEQUAL(ci_m, cj.molten())) {
						const int rsize = mergeProof(x, ci, cj, proofBytes);
						if (rsize == 1) 
							nElements++;
						else if (rsize) {
							if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}
			
			// GUARD for compressed proof size and #units
			if (nElements > ADDEDCLS_MAX || proofBytes > ADDEDPROOF_MAX) return true;
			nElements = ENCODEPROOFINFO(nElements, proofBytes);

			#if PROOF_DBG
			printf("c  Variable %d: counted %d units and %d proof bytes\n", x, nElements, proofBytes);
			#endif
		}
		else {
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				const bool ci_m = ci.molten();
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.original() && NEQUAL(ci_m, cj.molten())) {
						const int rsize = merge(x, ci, cj);
						if (rsize == 1) 
							nElements++;
						else if (rsize) {
							if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}
		}
		
		// GUARD for compressed variable limits
		if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX) return true;
		
		// check bound on literals
		if (kOpts->ve_lbound_en) {
			uint32 nLitsBefore = 0;
			countLitsBefore(cnf, me, nLitsBefore);
			countLitsBefore(cnf, other, nLitsBefore);
			if (nAddedLits > nLitsBefore) return true;
		}
		return false;
	}

	_PFROST_D_ void toblivion(
		const uint32& p,
		const uint32& n,
		const uint32& pOrgs,
		const uint32& nOrgs,
		CNF& cnf,
		OL& poss, 
		OL& negs,
		cuVecU* resolved)
	{
		const bool which = pOrgs > nOrgs;
		if (which) {
			uint32 nsLits = 0;
			countLitsBefore(cnf, negs, nsLits);
			uint32* saved = resolved->jump(nOrgs + nsLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), nOrgs + nsLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(negs, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, n);
				c.markDeleted();
			}
			saveWitness(saved, p);
		}
		else {
			uint32 psLits = 0;
			countLitsBefore(cnf, poss, psLits);
			uint32* saved = resolved->jump(pOrgs + psLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), pOrgs + psLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(poss, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, p);
				c.markDeleted();
			}

			saveWitness(saved, n);
		}

		OL& other = which ? poss : negs;
		forall_occurs(other, i) 
			cnf[*i].markDeleted();

		poss.clear(true), negs.clear(true);
	}

	_PFROST_D_ void toblivion(CNF& cnf, OL& poss, OL& negs)
	{
		forall_occurs(poss, i) 
			cnf[*i].markDeleted();

		forall_occurs(negs, i) 
			cnf[*i].markDeleted();

		poss.clear(true), negs.clear(true);
	}

	_PFROST_D_ void saveResolved(
		const uint32& p,
		const uint32& n,
		CNF& cnf, 
		OL& poss,
		OL& negs, 
		cuVecU* resolved)
	{
		const bool which = poss.size() > negs.size();
		if (which) {
			uint32 nsCls = 0, nsLits = 0;
			countOrgs(cnf, negs, nsCls, nsLits);
			uint32* saved = resolved->jump(nsCls + nsLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), nsCls + nsLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(negs, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, n);
			}

			saveWitness(saved, p);
		}
		else {
			uint32 psCls = 0, psLits = 0;
			countOrgs(cnf, poss, psCls, psLits);
			uint32* saved = resolved->jump(psCls + psLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), psCls + psLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(poss, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, p);
			}

			saveWitness(saved, n);
		}
	}

	_PFROST_D_ void saveResolved(
		const uint32& p,
		const uint32& n,
		const uint32& pOrgs, 
		const uint32& nOrgs,
		CNF& cnf,
		OL& poss,
		OL& negs,
		cuVecU* resolved)
	{
		const bool which = pOrgs > nOrgs;
		if (which) {
			uint32 nsLits = 0;
			countLitsBefore(cnf, negs, nsLits);
			uint32* saved = resolved->jump(nOrgs + nsLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), nOrgs + nsLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(negs, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, n);
			}

			saveWitness(saved, p);
		}
		else {
			uint32 psLits = 0;
			countLitsBefore(cnf, poss, psLits);
			uint32* saved = resolved->jump(pOrgs + psLits + 2);

			#if MODEL_DBG
			printf("c  saving witness(%d) of length %d at position %d\n",
				ABS(p), pOrgs + psLits + 2, uint32(saved - resolved->data()));
			#endif

			forall_occurs(poss, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveClause(saved, c, p);
			}

			saveWitness(saved, n);
		}
	}

} 

#endif