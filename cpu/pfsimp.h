/***********************************************************************[pfsimp.h]
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

#ifndef __SIMP_
#define __SIMP_

#include "pfsort.h"
#include "pfsolve.h" 
using namespace pFROST;

namespace SIGmA {

#define VE_DBG 0
#define HSE_DBG 0
#define SUB_DBG 0
#define CE_POS_LMT 512
#define CE_NEG_LMT 512
#define	MAX_ERE_OUT 350
#define HSE_MAX_CL_SIZE 1000
#ifdef __linux__ 
#define COUNTFLIPS(X) while (__builtin_parity(++X)); 
#elif _WIN32
#define COUNTFLIPS(X) while (__popcnt(++X) & 1); 
#endif

	// OT sorting comparator
	struct CNF_CMP_KEY {
		bool operator () (S_REF x, S_REF y) {
			if (x->size() != y->size()) return x->size() < y->size();
			else if (x->lit(0) != y->lit(0)) return x->lit(0) < y->lit(0);
			else if (x->lit(1) != y->lit(1)) return x->lit(1) < y->lit(1);
			else if (x->size() > 2 && x->back() != y->back()) return x->back() < y->back();
			else return x->sig() < y->sig();
		}
	};
	struct CNF_CMP_SZ {
		bool operator () (S_REF x, S_REF y) {
			return x->size() < y->size();
		}
	};

	// Elimination sub-routines 

	inline bool checkMolten(OL& poss, OL& negs)
	{
		for (int i = 0; i < poss.size(); i++)
			if (poss[i]->molten()) return false;
		for (int i = 0; i < negs.size(); i++)
			if (negs[i]->molten()) return false;
		return true;
	}

	inline bool checkDeleted(OL& poss, OL& negs)
	{
		for (int i = 0; i < poss.size(); i++)
			if (poss[i]->deleted()) return false;
		for (int i = 0; i < negs.size(); i++)
			if (negs[i]->deleted()) return false;
		return true;
	}

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

	inline bool equal3(S_REF d, uint32 a, uint32 b, uint32 c) {
		assert(d->original());
		int found = 0;
		for (int i = 0; i < d->size(); i++) {
			uint32 lit = d->lit(i);
			if (a != lit && b != lit && c != lit) return false;
			found++;
		}
		return found == 3;
	}

	inline bool isTautology(const uint32& elim_v, const S_REF c1, const S_REF c2)
	{
		assert(elim_v > 0);
		assert(!c1->deleted());
		assert(!c2->deleted());
		int it1 = 0, it2 = 0;
		while (it1 < c1->size() && it2 < c2->size()) {
			uint32 v1 = ABS(c1->lit(it1)), v2 = ABS(c2->lit(it2));
			if (v1 == elim_v) { it1++; }
			else if (v2 == elim_v) { it2++; }
			else if ((c1->lit(it1) ^ c2->lit(it2)) == NEG_SIGN) return true; // tautology detected ==> abort
			else if (v1 < v2) it1++;
			else if (v2 < v1) it2++;
			else { it1++; it2++; }
		}
		return false;
	}

	inline void merge(const uint32& elim_v, const S_REF c1, const S_REF c2, Lits_t& out_c)
	{
		assert(elim_v > 0);
		assert(c1->original());
		assert(c2->original());
		out_c.clear();
		int it1 = 0, it2 = 0;
		register uint32 lit1, lit2, v1, v2;
		while (it1 < c1->size() && it2 < c2->size()) {
			lit1 = c1->lit(it1);
			lit2 = c2->lit(it2);
			v1 = ABS(lit1);
			v2 = ABS(lit2);
			if (v1 == elim_v) { it1++; }
			else if (v2 == elim_v) { it2++; }
			else if (v1 < v2) { it1++; out_c.push(lit1); }
			else if (v2 < v1) { it2++; out_c.push(lit2); }
			else { // repeated literal
				assert(lit1 == lit2);
				it1++; it2++;
				out_c.push(lit1);
			}
		}
		while (it1 < c1->size()) {
			if (ABS(c1->lit(it1)) == elim_v) it1++;
			else { out_c.push(c1->lit(it1)); it1++; }
		}
		while (it2 < c2->size()) {
			if (ABS(c2->lit(it2)) == elim_v) it2++;
			else { out_c.push(c2->lit(it2)); it2++; }
		}
	}

	inline bool merge_ere(const uint32& elim_var, const S_REF c1, const S_REF c2, Lits_t& out_c)
	{
		assert(elim_var);
		assert(!c1->deleted());
		assert(!c2->deleted());
		out_c.clear();
		int it1 = 0, it2 = 0;
		uint32 lit1, lit2, v1, v2;
		while (it1 < c1->size() && it2 < c2->size()) {
			lit1 = c1->lit(it1); lit2 = c2->lit(it2);
			v1 = ABS(lit1); v2 = ABS(lit2);
			if (v1 == elim_var) { it1++; }
			else if (v2 == elim_var) { it2++; }
			else if ((lit1 ^ lit2) == NEG_SIGN) return false; // tautology
			else if (v1 < v2) { it1++; out_c.push(lit1); }
			else if (v2 < v1) { it2++; out_c.push(lit2); }
			else { // repeated literal
				assert(lit1 == lit2);
				it1++; it2++;
				out_c.push(lit1);
			}
		}
		while (it1 < c1->size()) {
			lit1 = c1->lit(it1);
			if (ABS(lit1) == elim_var) it1++;
			else { it1++; out_c.push(lit1); }
		}
		while (it2 < c2->size()) {
			lit2 = c2->lit(it2);
			if (ABS(lit2) == elim_var) it2++;
			else { it2++; out_c.push(lit2); }
		}
		return true;
	}

	inline void substitute_single(const uint32& dx, SCLAUSE& org, const uint32& def)
	{
		assert(dx > 1);
		assert(def != dx);
		assert(org.original());
#if VE_DBG
		PFLCLAUSE(1, org, " Clause ");
#endif
		int n = 0;
		for (int i = 0; i < org.size(); i++) {
			if (org[i] == dx) org[n++] = def;
			else if (org[i] != def) org[n++] = org[i];
		}
		org.resize(n);
		Sort(org.data(), org.size(), LESS<uint32>());
		org.calcSig();
		assert(org.isSorted());
#if VE_DBG
		PFLCLAUSE(1, org, " Substituted to ");
#endif
		if (org.size() == 1 && pfrost->unassigned(*org)) pfrost->enqueue(*org);
	}

	inline void substitute_single(const uint32& p, const uint32& def, OL& poss, OL& negs)
	{
		assert(def > 1);
		assert(!SIGN(p));
		// substitute negatives 
		uint32 n = NEG(p);
		for (int i = 0; i < negs.size(); i++) {
			if (negs[i]->learnt() || negs[i]->has(def)) negs[i]->markDeleted(); // learnt or tautology
			else substitute_single(n, *negs[i], FLIP(def));
		}
		// substitute positives
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->learnt() || poss[i]->has(FLIP(def))) poss[i]->markDeleted();
			else substitute_single(p, *poss[i], def);
		}
	}

	inline void substitute_x(const uint32& x, OL& poss, OL& negs, Lits_t& out_c)
	{
		assert(x);
		out_c.clear();
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->learnt()) continue;
			for (int j = 0; j < negs.size(); j++) {
				if (negs[j]->learnt()) continue;
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

	inline bool subset_sig(const uint32& A, const uint32& B) { return !(A & ~B); }

	inline bool selfSubset_sig(const uint32& A, const uint32& B)
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
			if (c1[it] != c2[it]) return false;
			else it++;
		}
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

	inline S_REF fast_equality_check(OT& ot, uint32 x, uint32 y, uint32 z) {
		if (ot[y].size() > ot[z].size()) swap(y, z);
		if (ot[x].size() > ot[y].size()) swap(x, y);
		OL& list = ot[x];
		sort3(x, y, z);
		assert(x <= y && y <= z && x <= z);
		for (S_REF* i = list; i != list.end(); i++) {
			SCLAUSE& c = **i;
			if (c.molten()) continue;
			assert(c.isSorted());
			if (c.original() && c.size() == 3 &&
				c[0] == x && c[1] == y && c[2] == z) return *i;
		}
		return NULL;
	}

	inline void forward_equ(Lits_t& m_c, OT& ot, const CL_ST& type)
	{
		pfrost->getStats().n_triedreduns++;
		int msize = m_c.size();
		assert(msize > 1);
		uint32 best = *m_c, m_sig = MAPHASH(best);
		assert(best > 1);
		int minsize = ot[best].size();
		for (int k = 1; k < msize; k++) {
			int lsize = ot[m_c[k]].size();
			if (lsize < minsize) minsize = lsize, best = m_c[k];
			m_sig |= MAPHASH(m_c[k]);
		}
		OL& minList = ot[best];
		for (int i = 0; i < minList.size(); i++) {
			CL_ST st = minList[i]->status();
			if (msize == minList[i]->size() && ((st & LEARNT) || (st & type)) &&
				subset_sig(m_sig, minList[i]->sig()) && isEqual(*minList[i], m_c)) {
				minList[i]->markDeleted();  //  HR found --> eliminate
				if (st & LEARNT) pfrost->getStats().n_lrnreduns++;
				else pfrost->getStats().n_orgreduns++;
				break;
			}
		}
	}

	inline void updateOL(OL& ol)
	{
		if (ol.empty()) return;
		S_REF* i, * j, * rend = ol.end();
		for (i = ol, j = i; i != rend; i++) {
			SCLAUSE& c = **i;
			if (c.molten()) c.freeze();
			else if (!c.deleted()) *j++ = *i;
		}
		ol.shrink(int(rend - j));
	}

	inline bool subset(const S_REF sm, const S_REF lr)
	{
		assert(!sm->deleted());
		assert(!lr->deleted());
		assert(sm->size() > 1);
		assert(lr->size() > 1);
		assert(sm->size() <= lr->size());
		int it1 = 0, it2 = 0, sub = 0;
		while (it1 < sm->size() && it2 < lr->size()) {
			if (sm->lit(it1) < lr->lit(it2)) it1++;
			else if (lr->lit(it2) < sm->lit(it1)) it2++;
			else { sub++; it1++; it2++; }
		}
		if (sub == sm->size())
			return true;
		return false;
	}

	inline bool selfSubset(const uint32& x, const S_REF sm, const S_REF lr)
	{
		assert(!sm->deleted());
		assert(!lr->deleted());
		assert(sm->size() > 1);
		assert(lr->size() > 1);
		assert(sm->size() <= lr->size());
		int it1 = 0, it2 = 0, sub = 0;
		bool self = false;
		while (it1 < sm->size() && it2 < lr->size()) {
			if (sm->lit(it1) == FLIP(x)) it1++;
			else if (lr->lit(it2) == x) { self = true; it2++; }
			else if (sm->lit(it1) < lr->lit(it2)) it1++;
			else if (lr->lit(it2) < sm->lit(it1)) it2++;
			else { sub++; it1++; it2++; }
		}
		if ((sub + 1) == sm->size()) {
			if (self) return true;
			else {
				while (it2 < lr->size()) {
					if (lr->lit(it2) == x) return true;
					it2++;
				}
			}
		}
		return false;
	}

	inline void freeze_binaries(OL& list)
	{
		for (S_REF* i = list; i != list.end(); i++)
			if ((*i)->size() == 2) (*i)->freeze();
	}

	inline void freeze_arities(OL& list)
	{
		for (S_REF* i = list; i != list.end(); i++)
			if ((*i)->size() > 2 && (*i)->molten()) (*i)->freeze();
	}

	inline void freeze_arities(OL& me, OL& other)
	{
		for (S_REF* i = me; i != me.end(); i++) {
			S_REF c = *i;
			if (c->size() > 2 && c->molten())
				c->freeze();
		}
		for (S_REF* i = other; i != other.end(); i++) {
			S_REF c = *i;
			if (c->size() > 2 && c->molten())
				c->freeze();
		}
	}

	inline void copyClause(SCLAUSE& c, Lits_t& out_c)
	{
		out_c.clear();
		for (uint32* k = c; k != c.end(); k++) {
			out_c.push(*k);
		}
	}

	inline bool checkArity(SCLAUSE& c, uint32* literals, const int& size)
	{
		assert(literals);
		const uint32* end = literals + size;
		for (uint32* i = c; i != c.end(); i++) {
			const uint32 lit = *i;
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

	inline bool makeArity(OT& ot, uint32& parity, uint32* literals, const int& size)
	{
		const uint32 oldparity = parity;
		COUNTFLIPS(parity);
		for (int k = 0; k < size; k++) {
			const uint32 bit = (1UL << k);
			if ((parity & bit) != (oldparity & bit))
				literals[k] = FLIP(literals[k]);
		}
		// search for an arity clause
		assert(size > 2);
		uint32 best = *literals;
		CHECKLIT(best);
		int minsize = ot[best].size();
		for (int k = 1; k < size; k++) {
			const uint32 lit = literals[k];
			CHECKLIT(lit);
			int lsize = ot[lit].size();
			if (lsize < minsize) {
				minsize = lsize;
				best = literals[k];
			}
		}
		OL& minlist = ot[best];
		for (S_REF* i = minlist; i != minlist.end(); i++) {
			SCLAUSE& c = **i;
			if (c.original() && c.size() == size && checkArity(c, literals, size)) {
				assert(c.original());
				c.melt();
				return true;
			}
		}
		return false;
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

	inline void countSubstituted(const uint32& x, OL& me, OL& other, int& nAddedCls, int& nAddedLits)
	{
		for (int i = 0; i < me.size(); i++) {
			if (me[i]->learnt()) continue;
			for (int j = 0; j < other.size(); j++) {
				if (other[j]->learnt()) continue;
				bool a = me[i]->molten(), b = other[j]->molten();
				if (a != b && !isTautology(x, me[i], other[j]))
					nAddedCls++, nAddedLits += me[i]->size() + other[j]->size() - 2;
			}
		}
	}

	inline void countSubstituted(const uint32& x, OL& me, OL& other, int& nAddedCls)
	{
		for (int i = 0; i < me.size(); i++) {
			if (me[i]->learnt()) continue;
			for (int j = 0; j < other.size(); j++) {
				if (other[j]->learnt()) continue;
				bool a = me[i]->molten(), b = other[j]->molten();
				if (a != b && !isTautology(x, me[i], other[j]))
					nAddedCls++;
			}
		}
	}

	inline void countResolvents(const uint32& x, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, int& nAddedCls, int& nAddedLits)
	{
		int nTs = 0;
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->learnt()) continue;
			for (int j = 0; j < negs.size(); j++) {
				if (negs[j]->learnt()) continue;
				if (isTautology(x, poss[i], negs[j])) nTs++;
				else nAddedLits += poss[i]->size() + negs[j]->size() - 2;
			}
		}
		assert(pOrgs * nOrgs >= nTs);
		nAddedCls = pOrgs * nOrgs - nTs;
	}

	inline void	saveResolved(uVec1D& resolved, const uint32& lit) { resolved.push(lit), resolved.push(1); }

	inline void	saveResolved(uVec1D& resolved, SCLAUSE& c, const uint32& x)
	{
		assert(c.original());
		const uint32 last = resolved.size();
		int pos = -1;
		for (int i = 0; i < c.size(); i++) {
			const uint32 lit = c[i];
			if (lit == x)
				pos = i;
			resolved.push(c[i]);
		}
		assert(pos >= 0);
		if (pos)
			swap(resolved[pos + last], resolved[last]);
		resolved.push(c.size());
	}

	inline void saveResolved(const uint32& p, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, uVec1D& resolved)
	{
		assert(p > 1);
		bool which = pOrgs > nOrgs;
		if (which) {
			for (int i = 0; i < negs.size(); i++)
				if (negs[i]->original()) saveResolved(resolved, *negs[i], NEG(p));
			saveResolved(resolved, p);
		}
		else {
			for (int i = 0; i < poss.size(); i++)
				if (poss[i]->original()) saveResolved(resolved, *poss[i], p);
			saveResolved(resolved, NEG(p));
		}
	}

	inline void toblivion(const uint32& p, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, uVec1D& resolved)
	{
		assert(p > 1);
		bool which = pOrgs > nOrgs;
		if (which) {
			for (int i = 0; i < negs.size(); i++)
				if (negs[i]->original()) saveResolved(resolved, *negs[i], NEG(p));
			saveResolved(resolved, p);
		}
		else {
			for (int i = 0; i < poss.size(); i++)
				if (poss[i]->original()) saveResolved(resolved, *poss[i], p);
			saveResolved(resolved, NEG(p));
		}
		for (int i = 0; i < poss.size(); i++) poss[i]->markDeleted();
		for (int i = 0; i < negs.size(); i++) negs[i]->markDeleted();
		poss.clear(true), negs.clear(true);
	}

	inline uint32 find_sfanin(const uint32& gate_out, OL& list)
	{
		assert(gate_out > 1);
		uint32 imp = 0;
		int nImps = 0;
		for (S_REF* i = list; i != list.end(); i++) {
			SCLAUSE& c = **i;
			if (c.original() && c.size() == 2) {
				imp = FLIP(c[0] ^ c[1] ^ gate_out);
				nImps++;
			}
			if (nImps > 1) return 0; // cannot be a single-input gate
		}
		return imp;
	}

	inline uint32 find_BN_gate(const uint32& p, OL& poss, OL& negs)
	{
		assert(p > 1);
		assert(!SIGN(p));
		assert(checkMolten(poss, negs));
		uint32 n = NEG(p);
		uint32 first = find_sfanin(p, poss);
		if (first) {
			uint32 second = n, def = first;
			if (second < first) first = second, second = def;
			for (int i = 0; i < negs.size(); i++) {
				SCLAUSE& c = *negs[i];
				if (c.learnt()) continue;
				if (c.size() == 2 && c[0] == first && c[1] == second) {
#if VE_DBG
					PFLOG1(" Gate %d = -/+%d found", ABS(p), ABS(def));
					pfrost->printOL(poss), pfrost->printOL(negs);
#endif
					return def;
				}
			}
		}
		return 0;
	}

	inline void find_fanin(const uint32& gate_out, OL& list, Lits_t& out_c, uint32& sig)
	{
		assert(gate_out > 1);
		out_c.clear();
		sig = 0;
		uint32 imp = 0;
		for (S_REF* i = list; i != list.end(); i++) {
			SCLAUSE& c = **i;
			if (c.learnt()) continue;
			assert(!c.molten());
			if (c.size() == 2) {
				imp = FLIP(c[0] ^ c[1] ^ gate_out);
				out_c.push(imp);
				sig |= MAPHASH(imp);
				c.melt(); // mark as gate clause
			}
		}
	}

	inline bool find_AO_gate(const uint32& dx, const int& nOrgCls, OT& ot, Lits_t& out_c)
	{
		assert(dx > 1);
		assert(checkMolten(ot[dx], ot[FLIP(dx)]));
		out_c.clear();
		uint32 sig, x = ABS(dx);
		// (-) ==> look for AND , (+) ==> look for OR
#if VE_DBG
		const char* type = SIGN(dx) ? "AND" : "OR";
#endif
		OL& itarget = ot[dx]; 
		find_fanin(dx, itarget, out_c, sig);
		if (out_c.size() > 1) {
			uint32 f_dx = FLIP(dx);
			out_c.push(f_dx);
			sig |= MAPHASH(f_dx);
			Sort(out_c, LESS<uint32>());
			OL& otarget = ot[f_dx];
			for (int i = 0; i < otarget.size(); i++) {
				SCLAUSE& c = *otarget[i];
				if (c.learnt()) continue;
				if (c.size() == out_c.size() && subset_sig(c.sig(), sig) && isEqual(c, out_c)) {
					c.melt(); // mark as fanout clause
					// check resolvability
					int nAddedCls = 0;
					countSubstituted(x, itarget, otarget, nAddedCls);
					if (nAddedCls > nOrgCls) { c.freeze(); break; }
					// can be substituted
#if VE_DBG
					PFLOGN1(" Gate %d = %s(", ABS(dx), type);
					for (int k = 0; k < out_c.size(); k++) {
						if (ABS(out_c[k]) == ABS(dx)) continue;
						fprintf(stdout, " %d", ABS(out_c[k]));
						if (k < out_c.size() - 1) putc(',', stdout);
					}
					fprintf(stdout, " ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
					printGate(itarget, otarget);
#endif
					substitute_x(x, itarget, otarget, out_c);
					return true;
				}
			}
		}
		freeze_binaries(itarget);
		return false;
	}

	inline bool find_ITE_gate(const uint32& dx, const int& nOrgCls, OT& ot, Lits_t& out_c)
	{
		assert(checkMolten(ot[dx], ot[FLIP(dx)]));
		OL& itarget = ot[dx];
		for (S_REF* i = itarget; i != itarget.end(); i++) {
			SCLAUSE& ci = **i;
			if (ci.learnt() || ci.size() < 3 || ci.size() > 3) continue;
			assert(ci.original());
			uint32 xi = ci[0], yi = ci[1], zi = ci[2];
			if (yi == dx) swap(xi, yi);
			if (zi == dx) swap(xi, zi);
			assert(xi == dx);
			for (S_REF* j = i + 1; j != itarget.end(); j++) {
				SCLAUSE& cj = **j;
				if (cj.learnt() || cj.size() < 3 || cj.size() > 3) continue;
				assert(cj.original());
				uint32 xj = cj[0], yj = cj[1], zj = cj[2];
				if (yj == dx) swap(xj, yj);
				if (zj == dx) swap(xj, zj);
				assert(xj == dx);
				if (ABS(yi) == ABS(zj)) swap(yj, zj);
				if (ABS(zi) == ABS(zj)) continue;
				if (yi != FLIP(yj)) continue;
				uint32 f_dx = FLIP(dx);
				S_REF d1 = fast_equality_check(ot, f_dx, yi, FLIP(zi));
				if (!d1) continue;
				S_REF d2 = fast_equality_check(ot, f_dx, yj, FLIP(zj));
				if (!d2) continue;
				assert(d1->original());
				assert(d2->original());
				// mark gate clauses
				ci.melt(), cj.melt();
				d1->melt(), d2->melt();
				// check resolvability
				uint32 v = ABS(dx);
				int nAddedCls = 0;
				OL& otarget = ot[f_dx];
				countSubstituted(v, itarget, otarget, nAddedCls);
				if (nAddedCls > nOrgCls) {
					ci.freeze(), cj.freeze();
					d1->freeze(), d2->freeze();
					return false;
				}
				// can be substituted
#if VE_DBG
				PFLOG1(" Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d", l2i(dx), ABS(yi), ABS(zi), ABS(zj), nAddedCls, itarget.size() + otarget.size());
				printGate(itarget, otarget);
#endif
				substitute_x(v, itarget, otarget, out_c);
				return true;
			}
		}
		return false;
	}

	inline bool find_XOR_gate(const uint32& dx, const int& nOrgCls, OT& ot, Lits_t& out_c)
	{
		const uint32 fx = FLIP(dx), v = ABS(dx);
		assert(checkMolten(ot[dx], ot[fx]));
		OL& itarget = ot[dx];
		OL& otarget = ot[fx];
		for (S_REF* i = itarget; i != itarget.end(); i++) {
			SCLAUSE& ci = **i;
			if (ci.original()) {
				const int size = ci.size();
				const int arity = size - 1; // XOR arity
				if (size < 3 || arity > pfrost->opts.xor_max_arity) continue;
				// share to out_c
				copyClause(ci, out_c);
				// find arity clauses
				uint32 parity = 0;
				int itargets = (1 << arity);
				while (--itargets && makeArity(ot, parity, out_c, size));
				assert(parity < (1UL << size)); // overflow check
				assert(itargets >= 0);
				if (itargets)
					freeze_arities(itarget, otarget);
				else {
					ci.melt();
					// check resolvability
					int nAddedCls = 0;
					countSubstituted(v, itarget, otarget, nAddedCls);
					if (nAddedCls > nOrgCls) {
						freeze_arities(itarget, otarget);
						break;
					}
					// can be substituted
					if (verbose >= 4) {
						PFLOGN1(" Gate %d = XOR(", l2i(dx));
						for (int k = 0; k < out_c.size(); k++) {
							printf(" %d", ABS(out_c[k]));
							if (k < out_c.size() - 1) putchar(',');
						}
						printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
						printGate(itarget, otarget);
					}
					substitute_x(v, itarget, otarget, out_c);
					return true;
				}
			} // original block
		}
		return false;
	}

	inline bool resolve_x(const uint32& x, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, Lits_t& out_c, const bool& bound)
	{
		assert(x);
		assert(checkMolten(poss, negs));
		out_c.clear();
		uint32 p = V2L(x);
		// check resolvability
		int nAddedCls = 0, nAddedLits = 0;
		countResolvents(x, pOrgs, nOrgs, poss, negs, nAddedCls, nAddedLits);
		if (nAddedCls == 0) return true; // No resolvents to add
		if (nAddedCls > pOrgs + nOrgs) return false;
		// count literals before elimination
		if (bound) {
			int lits_before = 0;
			countLitsBefore(poss, lits_before);
			countLitsBefore(negs, lits_before);
			if (nAddedLits > lits_before) return false;
		}
		// can be eliminated
#if VE_DBG
		PFLOG1(" Resolving(%d) ==> added = %d, deleted = %d", x, nAddedCls, poss.size() + negs.size());
		pfrost->printOL(poss), pfrost->printOL(negs);
#endif
		for (int i = 0; i < poss.size(); i++) {
			if (poss[i]->learnt()) continue;
			for (int j = 0; j < negs.size(); j++) {
				if (negs[j]->learnt()) continue;
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
		return true;
	}

	inline void self_sub_x(const uint32& x, OL& poss, OL& negs)
	{
		assert(checkMolten(poss, negs));
		// positives vs negatives
		for (int i = 0; i < poss.size(); i++) {
			S_REF pos = poss[i];
			if (pos->size() > HSE_MAX_CL_SIZE) break;
			if (pos->deleted()) continue;
			// self-subsumption check
			for (int j = 0; j < negs.size(); j++) {
				S_REF neg = negs[j];
				if (neg->size() > pos->size()) break;
				if (neg->deleted()) continue;
				if (neg->size() > 1 && selfSubset_sig(neg->sig(), pos->sig()) && selfSubset(x, neg, pos)) {
#if HSE_DBG
					PFLCLAUSE(1, (*pos), " Clause ");
					PFLCLAUSE(1, (*neg), " Strengthened by ");
#endif 
					pfrost->strengthen(pos, x);
					pos->melt(); // mark for fast recongnition in ot update 
					break; // cannot strengthen "pos" anymore, 'x' already removed
				}
			}
			// subsumption check
			for (int j = 0; j < i; j++) {
				S_REF sm_c = poss[j];
				if (sm_c->deleted()) continue;
				if (pos->molten() && sm_c->size() > pos->size()) continue;
				if (sm_c->size() > 1 && subset_sig(sm_c->sig(), pos->sig()) && subset(sm_c, pos)) {
					if (sm_c->learnt() && pos->original()) sm_c->set_status(ORIGINAL);
#if HSE_DBG
					PFLCLAUSE(1, (*pos), " Clause ");
					PFLCLAUSE(1, (*sm_c), " Subsumed by ");
#endif 
					pos->markDeleted();
					break;
				}
			}
		}
		updateOL(poss);
		// negatives vs positives
		for (int i = 0; i < negs.size(); i++) {
			S_REF neg = negs[i];
			if (neg->size() > HSE_MAX_CL_SIZE) break;
			if (neg->deleted()) continue;
			// self-subsumption check
			for (int j = 0; j < poss.size(); j++) {
				S_REF pos = poss[j];
				if (pos->size() >= neg->size()) break;
				if (pos->deleted()) continue;
				if (pos->size() > 1 && selfSubset_sig(pos->sig(), neg->sig()) && selfSubset(NEG(x), pos, neg)) {
#if HSE_DBG
					PFLCLAUSE(1, (*neg), " Clause ");
					PFLCLAUSE(1, (*pos), " Strengthened by ");
#endif 
					pfrost->strengthen(neg, NEG(x));
					neg->melt();
					break;
				}
			}
			// subsumption check
			for (int j = 0; j < i; j++) {
				S_REF sm_c = negs[j];
				if (sm_c->deleted()) continue;
				if (neg->molten() && sm_c->size() > neg->size()) continue;
				if (sm_c->size() > 1 && subset_sig(sm_c->sig(), neg->sig()) && subset(sm_c, neg)) {
					if (sm_c->learnt() && neg->original()) sm_c->set_status(ORIGINAL);
#if HSE_DBG
					PFLCLAUSE(1, (*neg), " Clause ");
					PFLCLAUSE(1, (*sm_c), " Subsumed by ");
#endif 
					neg->markDeleted();
					break;
				}
			}
		}
		updateOL(negs);
		assert(checkMolten(poss, negs));
		assert(checkDeleted(poss, negs));
	}

	inline void blocked_x(const uint32& x, OL& poss, OL& negs)
	{
		// start with negs
		for (int i = 0; i < negs.size(); i++) {
			if (negs[i]->deleted() || negs[i]->learnt()) continue;
			bool allTautology = true;
			for (int j = 0; j < poss.size(); j++) {
				if (poss[j]->deleted() || poss[j]->learnt()) continue;
				if (!isTautology(x, negs[i], poss[j])) { allTautology = false; break; }
			}
			if (allTautology) negs[i]->markDeleted();
		}
	}

}

#endif