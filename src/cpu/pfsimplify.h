/***********************************************************************[pfsimplify.h]
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

#ifndef __SIMPLIFY_
#define __SIMPLIFY_

#include "pfsolve.h" 
#include "pfhse.h" 
using namespace pFROST;

namespace SIGmA {

	#define VE_DBG 0

	inline void printGate(OL& poss, OL& negs)
	{
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->molten()) {
				PFLOGN0(" ");
				poss[i]->print();
			}
		}
		for (int i = 0; i < negs.size(); i++) {
			if (negs[i]->molten()) {
				PFLOGN0(" ");
				negs[i]->print();
			}
		}
	}

	inline bool checkMolten(OL& poss, OL& negs)
	{
		for (int i = 0; i < poss.size(); i++)
			if (poss[i]->molten()) return false;
		for (int i = 0; i < negs.size(); i++)
			if (negs[i]->molten()) return false;
		return true;
	}

	inline void cswap(uint32& x, uint32& y)
	{
		uint32 ta = std::min(x, y);
		uint32 tb = std::max(x, y);
		x = ta, y = tb;
	}

	inline void sort3(uint32& x, uint32& y, uint32& z)
	{
		cswap(y, z);
		cswap(x, z);
		cswap(x, y);
	}

	inline bool isTautology(const uint32& x, const S_REF c1, const S_REF c2)
	{
		assert(x > 0);
		assert(!c1->deleted());
		assert(!c2->deleted());
		int it1 = 0, it2 = 0;
		while (it1 < c1->size() && it2 < c2->size()) {
			uint32 v1 = ABS(c1->lit(it1)), v2 = ABS(c2->lit(it2));
			if (v1 == x) { it1++; }
			else if (v2 == x) { it2++; }
			else if ((c1->lit(it1) ^ c2->lit(it2)) == NEG_SIGN) return true; // tautology detected ==> abort
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { it1++; it2++; }
		}
		return false;
	}

	inline void merge(const uint32& x, const S_REF c1, const S_REF c2, Lits_t& out_c)
	{
		assert(x > 0);
		assert(c1->original());
		assert(c2->original());
		out_c.clear();
		int it1 = 0, it2 = 0;
		uint32 lit1, lit2, v1, v2;
		while (it1 < c1->size() && it2 < c2->size()) {
			lit1 = c1->lit(it1);
			lit2 = c2->lit(it2);
			v1 = ABS(lit1);
			v2 = ABS(lit2);
			if (v1 == x) { it1++; }
			else if (v2 == x) { it2++; }
			else if (v1 < v2) { it1++; out_c.push(lit1); }
			else if (v2 < v1) { it2++; out_c.push(lit2); }
			else { // repeated literal
				assert(lit1 == lit2);
				it1++; it2++;
				out_c.push(lit1);
			}
		}
		while (it1 < c1->size()) {
			if (ABS(c1->lit(it1)) == x) it1++;
			else { out_c.push(c1->lit(it1)); it1++; }
		}
		while (it2 < c2->size()) {
			if (ABS(c2->lit(it2)) == x) it2++;
			else { out_c.push(c2->lit(it2)); it2++; }
		}
	}

	inline int merge(const uint32& x, const S_REF c1, const S_REF c2)
	{
		assert(x > 0);
		assert(c1->original());
		assert(c2->original());
		const int n1 = c1->size(), n2 = c2->size();
		int it1 = 0, it2 = 0;
		uint32 lit1, lit2, v1, v2;
		int len = n1 + n2 - 2;
		while (it1 < n1 && it2 < n2) {
			lit1 = c1->lit(it1);
			lit2 = c2->lit(it2);
			v1 = ABS(lit1);
			v2 = ABS(lit2);
			if (v1 == x) it1++;
			else if (v2 == x) it2++;
			else if ((lit1 ^ lit2) == NEG_SIGN) return 0;
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { // repeated literal
				assert(lit1 == lit2);
				it1++, it2++;
				assert(len > 0);
				len--;
			}
		}
		return len;
	}

	inline void countOrgs(OL& list, int& orgs)
	{
		assert(!orgs);
		for (S_REF* i = list; i != list.end(); i++)
			if ((*i)->original()) orgs++;
	}

	inline void countLitsBefore(OL& list, int& nLitsBefore)
	{
		for (S_REF* i = list; i != list.end(); i++) {
			if ((*i)->original()) nLitsBefore += (*i)->size();
		}
	}

	inline bool countSubstituted(const uint32& x, const int& clsbefore, OL& me, OL& other, int& nAddedCls)
	{
		assert(!nAddedCls);
		int nAddedLits = 0;
		const int rlimit = pfrost->opts.ve_clause_limit;
		for (int i = 0; i < me.size(); i++) {
			if (me[i]->original()) {
				const bool a = me[i]->molten();
				for (int j = 0; j < other.size(); j++) {
					if (other[j]->original()) {
						const bool b = other[j]->molten();
						int rsize;
						if (a != b && (rsize = merge(x, me[i], other[j]))) {
							if (++nAddedCls > clsbefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}
		}
		if (pfrost->opts.ve_lbound_en) {
			int nLitsBefore = 0;
			countLitsBefore(me, nLitsBefore);
			countLitsBefore(other, nLitsBefore);
			if (nAddedLits > nLitsBefore) return true;
		}
		return false;
	}

	inline bool countResolvents(const uint32& x, const int& clsbefore, OL& poss, OL& negs, int& nAddedCls)
	{
		assert(!nAddedCls);
		int nAddedLits = 0;
		const int rlimit = pfrost->opts.ve_clause_limit;
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->original()) {
				for (int j = 0; j < negs.size(); j++) {
					if (negs[j]->original()) {
						const int rsize = merge(x, poss[i], negs[j]);
						if (rsize && ((++nAddedCls > clsbefore) || (rlimit && rsize > rlimit))) return true;
						nAddedLits += rsize;
					}
				}
			}
		}
		if (pfrost->opts.ve_lbound_en) {
			int nLitsBefore = 0;
			countLitsBefore(poss, nLitsBefore);
			countLitsBefore(negs, nLitsBefore);
			if (nAddedLits > nLitsBefore) return true;
		}
		return false;
	}

	inline void freezeBinaries(OL& list)
	{
		for (S_REF* i = list; i != list.end(); i++) {
			SCLAUSE& c = **i;
			if (c.original() && c.size() == 2)
				c.freeze();
		}
	}

	inline void toblivion(OL& list)
	{
		for (int i = 0; i < list.size(); i++) pfrost->removeClause(list[i]);
		list.clear(true);
	}

	inline void toblivion(const uint32& p, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, MODEL& model)
	{
		CHECKLIT(p);
		const uint32 n = NEG(p);
		PFLOG2(4, " saving clauses of eliminated(%d) as witness", ABS(p));
		if (pOrgs > nOrgs) {
			for (int i = 0; i < negs.size(); i++) {
				SCLAUSE& c = *negs[i];
				if (c.original())
					model.saveClause(c, c.size(), n);
				c.markDeleted();
			}
			model.saveWitness(p);
			negs.clear(true);
			toblivion(poss);
		}
		else {
			for (int i = 0; i < poss.size(); i++) {
				SCLAUSE& c = *poss[i];
				if (c.original())
					model.saveClause(c, c.size(), p);
				c.markDeleted();
			}
			model.saveWitness(n);
			poss.clear(true);
			toblivion(negs);
		}
	}

	inline void substitute_x(const uint32& x, OL& poss, OL& negs, Lits_t& out_c)
	{
		assert(x);
#if VE_DBG
		PFLOG1(" Substituting(%d):", x);
		pfrost->printOL(poss), pfrost->printOL(negs);
#endif
		out_c.clear();
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->original()) {
				for (int j = 0; j < negs.size(); j++) {
					if (negs[j]->original()) {
						bool a = poss[i]->molten(), b = negs[j]->molten();
						if (a != b && !isTautology(x, poss[i], negs[j])) {
							merge(x, poss[i], negs[j], out_c);
							S_REF added = new SCLAUSE(out_c);
							pfrost->newResolvent(added);
#if VE_DBG
							PFLCLAUSE(1, (*added), " Added ");
#endif
						}
					}
				}
			}
		}
	}

	inline void resolve_x(const uint32& x, OL& poss, OL& negs, Lits_t& out_c)
	{
		assert(x);
#if VE_DBG
		PFLOG1(" Resolving(%d):", x);
		pfrost->printOL(poss), pfrost->printOL(negs);
#endif
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->original()) {
				for (int j = 0; j < negs.size(); j++) {
					if (negs[j]->original()) {
						if (!isTautology(x, poss[i], negs[j])) {
							merge(x, poss[i], negs[j], out_c);
							S_REF added = new SCLAUSE(out_c);
							pfrost->newResolvent(added);
#if VE_DBG
							PFLCLAUSE(1, (*added), " Resolvent ");
#endif
						}
					}
				}
			}
		}
	}

	inline bool mayResolve_x(const uint32& x, const int& nOrgCls, OL& poss, OL& negs, int& nAddedCls)
	{
		assert(x);
		assert(checkMolten(poss, negs));
		// check resolvability
		nAddedCls = 0;
		if (countResolvents(x, nOrgCls, poss, negs, nAddedCls)) return false;
		return true;
	}

	inline void blocked_x(const uint32& x, OL& poss, OL& negs, MODEL& model)
	{
		const uint32 p = V2L(x), n = NEG(p);
		// start with negs
		for (int i = 0; i < negs.size(); i++) {
			if (negs[i]->original()) {
				bool allTautology = true;
				for (int j = 0; j < poss.size(); j++) {
					if (poss[j]->original() && !isTautology(x, negs[i], poss[j])) { 
						allTautology = false;
						break; 
					}
				}
				if (allTautology) {
					SCLAUSE& neg = *negs[i];
					assert(neg.original());
					model.saveClause(neg, neg.size(), n);
					negs[i]->markDeleted();
				}
			}
		}
	}

}

#endif