/***********************************************************************[simplify.h]
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

#include "solve.h" 
using namespace ParaFROST;

#define RESOLUTION 1
#define SUBSTITUTION 2
#define CORESUBSTITUTION 4

struct CNF_CMP_KEY {
	inline bool operator () (const S_REF x, const S_REF y) {
		const int xsize = x->size(), ysize = y->size();
		if (xsize < ysize) return true;
		if (xsize > ysize) return false;
		const uint32 lit0 = **x, lit1 = **y;
		if (lit0 < lit1) return true;
		if (lit0 > lit1) return false;
		const uint32 xback = x->back(), yback = y->back();
		if (xback < yback) return true;
		if (xback > yback) return false;
		return x->sig() < y->sig();
	}
};

inline void printGate(const OL& poss, const OL& negs)
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

inline bool checkDeleted(OL& poss, OL& negs)
{
	for (int i = 0; i < poss.size(); i++)
		if (poss[i]->deleted()) return false;
	for (int i = 0; i < negs.size(); i++)
		if (negs[i]->deleted()) return false;
	return true;
}

inline bool checkMolten(const OL& poss, const OL& negs)
{
	for (int i = 0; i < poss.size(); i++) {
		S_REF c = poss[i];
		if (!c->deleted() && c->molten())
			return false;
	}
	for (int i = 0; i < negs.size(); i++) {
		S_REF c = negs[i];
		if (!c->deleted() && c->molten())
			return false;
	}
	return true;
}

inline bool sub(const uint32& A, const uint32& B) { return !(A & ~B); }

inline bool selfsub(const uint32& A, const uint32& B)
{
	uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
	return !(A & ~B_tmp);
}

inline bool isEqual(const SCLAUSE& c1, const Lits_t& c2)
{
	assert(!c1.deleted());
	assert(c1.size() == c2.size());
	int it = 0;
	while (it < c2.size()) {
		if (NEQUAL(c1[it], c2[it]))
			return false;
		else it++;
	}
	return true;
}

inline void cswap(uint32& x, uint32& y)
{
	const uint32 ta = MIN(x, y);
	const uint32 tb = MAX(x, y);
	x = ta, y = tb;
}

inline void sort3(uint32& x, uint32& y, uint32& z)
{
	cswap(y, z);
	cswap(x, z);
	cswap(x, y);
}

inline bool isTautology(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2)
{
	CHECKVAR(x);
	assert(c1.original());
	assert(c2.original());
	const int n1 = c1.size(), n2 = c2.size();
	int it1 = 0, it2 = 0;
	while (it1 < n1 && it2 < n2) {
		const uint32 lit1 = c1[it1], lit2 = c2[it2];
		const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
		if (v1 == x) { it1++; }
		else if (v2 == x) { it2++; }
		else if ((lit1 ^ lit2) == NEG_SIGN) return true; // tautology detected ==> abort
		else if (v1 < v2) it1++;
		else if (v2 < v1) it2++;
		else { it1++; it2++; }
	}
	return false;
}

inline bool merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2, Lits_t& out_c)
{
	CHECKVAR(x);
	assert(c1.original());
	assert(c2.original());
	out_c.clear();
	const int n1 = c1.size(), n2 = c2.size();
	int it1 = 0, it2 = 0;
	uint32 lit1, lit2, v1, v2;
	while (it1 < n1 && it2 < n2) {
		lit1 = c1[it1];
		lit2 = c2[it2];
		v1 = ABS(lit1);
		v2 = ABS(lit2);
		if (v1 == x) { it1++; }
		else if (v2 == x) { it2++; }
		else if ((lit1 ^ lit2) == NEG_SIGN) return false;
		else if (v1 < v2) { it1++; out_c.push(lit1); }
		else if (v2 < v1) { it2++; out_c.push(lit2); }
		else { // repeated literal
			assert(lit1 == lit2);
			it1++, it2++;
			out_c.push(lit1);
		}
	}
	while (it1 < n1) {
		lit1 = c1[it1++];
		if (NEQUAL(ABS(lit1), x)) out_c.push(lit1);
	}
	while (it2 < n2) {
		lit2 = c2[it2++];
		if (NEQUAL(ABS(lit2), x)) out_c.push(lit2);
	}
	return true;
}

inline int merge(const uint32& x, const SCLAUSE& c1, const SCLAUSE& c2)
{
	CHECKVAR(x);
	assert(c1.original());
	assert(c2.original());
	const int n1 = c1.size(), n2 = c2.size();
	int it1 = 0, it2 = 0;
	uint32 lit1, lit2, v1, v2;
	int len = n1 + n2 - 2;
	while (it1 < n1 && it2 < n2) {
		lit1 = c1[it1];
		lit2 = c2[it2];
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

inline void freezeBinaries(OL& list)
{
	forall_occurs(list, i) {
		SCLAUSE& c = **i;
		if (c.original() && c.size() == 2)
			c.freeze();
	}
}

inline void freezeAll(OL& list)
{
	forall_occurs(list, i) {
		SCLAUSE& c = **i;
		if (c.original()) c.freeze();
	}
}

inline void countOrgs(OL& list, int& orgs)
{
	assert(!orgs);
	forall_occurs(list, i) {
		if ((*i)->original()) orgs++;
	}
}

inline void countLitsBefore(OL& list, int& nLitsBefore)
{
	forall_occurs(list, i) {
		S_REF c = *i;
		if (c->original()) nLitsBefore += c->size();
	}
}

inline bool countSubstituted(const uint32& x, const int& clsbefore, OL& me, OL& other, int& nAddedCls)
{
	assert(!nAddedCls);
	int nAddedLits = 0;
	const int rlimit = solver->opts.ve_clause_limit;
	for (int i = 0; i < me.size(); i++) {
		const SCLAUSE& ci = *me[i];
		if (ci.original()) {
			const bool a = ci.molten();
			for (int j = 0; j < other.size(); j++) {
				const SCLAUSE& cj = *other[j];
				if (cj.original()) {
					const bool b = cj.molten();
					int rsize;
					if (a != b && (rsize = merge(x, ci, cj)) > 1) {
						if (++nAddedCls > clsbefore || (rlimit && rsize > rlimit)) return true;
						nAddedLits += rsize;
					}
				}
			}
		}
	}
	if (solver->opts.ve_lbound_en) {
		int nLitsBefore = 0;
		countLitsBefore(me, nLitsBefore);
		countLitsBefore(other, nLitsBefore);
		if (nAddedLits > nLitsBefore) return true;
	}
	return false;
}

inline bool countResolvents(const uint32& x, const int& clsbefore, OL& me, OL& other, int& nAddedCls)
{
	assert(!nAddedCls);
	int nAddedLits = 0;
	const int rlimit = solver->opts.ve_clause_limit;
	for (int i = 0; i < me.size(); i++) {
		const SCLAUSE& ci = *me[i];
		if (ci.original()) {
			for (int j = 0; j < other.size(); j++) {
				const SCLAUSE& cj = *other[j];
				if (cj.original()) {
					int rsize;
					if ((rsize = merge(x, ci, cj)) > 1) {
						if (++nAddedCls > clsbefore || (rlimit && rsize > rlimit)) return true;
						nAddedLits += rsize;
					}
				}
			}
		}
	}
	if (solver->opts.ve_lbound_en) {
		int nLitsBefore = 0;
		countLitsBefore(me, nLitsBefore);
		countLitsBefore(other, nLitsBefore);
		if (nAddedLits > nLitsBefore) return true;
	}
	return false;
}

inline void toblivion(OL& list)
{
	forall_occurs(list, i) {
		(*i)->markDeleted();
	}
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

#endif