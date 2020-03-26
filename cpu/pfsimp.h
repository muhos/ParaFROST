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

#define LIT_REM_THR 10
#define MIN_PARALLEL_VARS 10
#define CE_POS_LMT 256
#define CE_NEG_LMT 256

/* Elimination sub-routines */

inline bool isTautology(const uint32& elim_v, const S_REF c1, const S_REF c2)
{
	assert(elim_v > 0);
	assert(c1->status() != DELETED);
	assert(c2->status() != DELETED);
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

inline bool isTautology_and(const S_REF org_c, const uVector1D& defs)
{
	assert(org_c->status() != DELETED);
	for (int i = 0; i < defs.size(); i++) {
		if (org_c->has(defs[i])) return true;
	}
	return false;
}

inline bool isTautology_or(const S_REF org_c, const uVector1D& defs)
{
	assert(org_c->status() != DELETED);
	for (int i = 0; i < defs.size(); i++) {
		if (org_c->has(FLIP(defs[i]))) return true;
	}
	return false;
}

inline void clause_extend_and(const uint32& neg_x, const S_REF org_c, const uVector1D& defs, uVector1D& out_c)
{
	assert(org_c->status() != DELETED);
	for (int i = 0; i < org_c->size(); i++) {
		if (org_c->lit(i) == neg_x) {
			for (int k = 0; k < defs.size(); k++) out_c.push(FLIP(defs[k]));
		}
		else out_c.push(org_c->lit(i));
	}
}

inline void clause_extend_or(const uint32& x, const S_REF org_c, const uVector1D& defs, uVector1D& out_c)
{
	assert(org_c->status() != DELETED);
	for (int i = 0; i < org_c->size(); i++) {
		if (org_c->lit(i) == x) {
			for (int k = 0; k < defs.size(); k++) out_c.push(defs[k]);
		}
		else out_c.push(org_c->lit(i));
	}
}

inline void clause_split_and(const uint32& x, const S_REF org_c, const uint32& def, S_REF added)
{
	assert(org_c->status() != DELETED);
	for (int k = 0; k < org_c->size(); k++) {
		if (org_c->lit(k) == def) continue; // repeated literal
		if (org_c->lit(k) == x) added->push(def);
		else added->push(org_c->lit(k));
	}
	Sort(added->d_ptr(), added->size());
	assert(added->isSorted());
	assert(added->hasZero() < 0);
}

inline void clause_split_or(const uint32& neg_x, const S_REF org_c, const uint32& def, S_REF added)
{
	assert(org_c->status() != DELETED);
	for (int k = 0; k < org_c->size(); k++) {
		if (org_c->lit(k) == FLIP(def)) continue; // repeated literal
		if (org_c->lit(k) == neg_x) added->push(FLIP(def));
		else added->push(org_c->lit(k));
	}
	Sort(added->d_ptr(), added->size());
	assert(added->isSorted());
	assert(added->hasZero() < 0);
}

inline void merge(const uint32& elim_v, const S_REF c1, const S_REF c2, S_REF added)
{
	assert(elim_v > 0);
	assert(c1->status() != DELETED);
	assert(c2->status() != DELETED);
	int it1 = 0, it2 = 0;
	register uint32 lit1, lit2, v1, v2;
	while (it1 < c1->size() && it2 < c2->size()) {
		lit1 = c1->lit(it1);
		lit2 = c2->lit(it2);
		v1 = ABS(lit1);
		v2 = ABS(lit2);
		if (v1 == elim_v) { it1++; } // skip pv->melted var 
		else if (v2 == elim_v) { it2++; } // skip pv->melted var 
		else if (v1 < v2) { it1++; added->push(lit1); }
		else if (v2 < v1) { it2++; added->push(lit2); }
		else { // repeated literal
			assert(lit1 == lit2);
			it1++; it2++;
			added->push(lit1);
		}
	}
	while (it1 < c1->size()) {
		if (ABS(c1->lit(it1)) == elim_v) it1++;
		else { added->push(c1->lit(it1)); it1++; }
	}
	while (it2 < c2->size()) {
		if (ABS(c2->lit(it2)) == elim_v) it2++;
		else { added->push(c2->lit(it2)); it2++; }
	}
	assert(added->isSorted());
	assert(added->hasZero() < 0);
}

inline bool merge_hre(const uint32& elim_var, const S_REF c1, const S_REF c2, uVector1D& out_c)
{
	assert(elim_var);
	assert(c1->status() != DELETED);
	assert(c2->status() != DELETED);
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

inline bool subset_sig(const uint32& A, const uint32& B)
{
	return (A & ~B) == 0;
}

inline bool selfSubset_sig(const uint32& A, const uint32& B)
{
	uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
	return (A & ~B_tmp) == 0;
}

inline bool subset(const S_REF sm, const S_REF lr)
{
	assert(sm->status() != DELETED);
	assert(lr->status() != DELETED);
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

inline bool isEqual(const SCLAUSE& c1, const uVector1D& c2)
{
	assert(c1.status() != DELETED);
	assert(c1.size() == c2.size());
	int it = 0;
	while (it < c2.size()) {
		if (c1[it] != c2[it]) return false;
		else it++;
	}
	return true;
}

inline bool selfSubset(const uint32& x, const S_REF sm, const S_REF lr)
{
	assert(sm->status() != DELETED);
	assert(lr->status() != DELETED);
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

inline bool forward_equ(const uVector1D& m_c, const uint32& m_sig, OL& x_ol)
{
	for (int i = 0; i < x_ol.size(); i++) {
		if (x_ol[i]->status() != DELETED && m_c.size() == x_ol[i]->size() &&
			subset_sig(m_sig, x_ol[i]->sig()) && isEqual(*x_ol[i], m_c)) {
			x_ol[i]->set_status(DELETED);  //  HR found --> eliminate
			return true;
		}
	}
	return false;
}

inline void resolve_x(const uint32& x, OL& p_ol, OL& n_ol)
{
	assert(x);
	for (int i = 0; i < p_ol.size(); i++) {
		for (int j = 0; j < n_ol.size(); j++) {
			if (!isTautology(x, p_ol[i], n_ol[j])) {
				S_REF added = new SCLAUSE();
				merge(x, p_ol[i], n_ol[j], added);
				g_pFrost->attachClause(added);
			}
		}
	}
	for (int i = 0; i < p_ol.size(); i++) p_ol[i]->set_status(DELETED);
	for (int j = 0; j < n_ol.size(); j++) n_ol[j]->set_status(DELETED);
}

inline void updateOL(OL& ol)
{
	if (ol.size() == 0) return;
	int idx = 0, ol_sz = ol.size();
	while (idx < ol_sz) {
		S_REF c_ref = ol[idx];
		if (c_ref->molten()) {
			c_ref->freeze();
			ol[idx] = ol[ol_sz - 1];
			ol_sz--;
		}
		else idx++;
	}
	int remOccurs = ol.size() - ol_sz;
	if (remOccurs) ol.shrink(remOccurs);
}

bool substitute_AND(const uint32& x, const uVector1D& defs, OL& p_ol, OL& n_ol, uVector1D& out_c)
{
	int numAddedClauses = 0;
	// count number of added clauses & literals for negatives
	for (int i = 0; i < n_ol.size(); i++) {
		if (!isTautology_and(n_ol[i], defs))
			numAddedClauses++;
	}
	// count number of added clauses & literals for positives
	for (int k = 0; k < defs.size(); k++) {
		for (int i = 0; i < p_ol.size(); i++) {
			if (!(p_ol[i]->has(FLIP(defs[k]))))
				numAddedClauses++;
		}
	}
	if (numAddedClauses > p_ol.size() + n_ol.size()) return false;
	// substitute negatives 
	for (int i = 0; i < n_ol.size(); i++) {
		if (!isTautology_and(n_ol[i], defs)) {
			assert(out_c.empty());
			clause_extend_and(NEG(x), n_ol[i], defs, out_c);
			Sort(out_c.d_ptr(), out_c.size());
			S_REF added = new SCLAUSE();
			added->push(out_c[0]);
			for (LIT_POS k = 1; k < out_c.size(); k++) {
				if (out_c[k - 1] == out_c[k]) continue; // skip repeated literal
				added->push(out_c[k]);
			}
			assert(added->isSorted());
			assert(added->hasZero() < 0);
			g_pFrost->attachClause(added);
			out_c.clear();
		}
	}
	// substitute positives
	for (int d = 0; d < defs.size(); d++) {
		for (int i = 0; i < p_ol.size(); i++) {
			if (!(p_ol[i]->has(FLIP(defs[d])))) {
				S_REF added = new SCLAUSE();
				clause_split_and(x, p_ol[i], defs[d], added);
				//assert(added->size() > 1);
				g_pFrost->attachClause(added);
			}
		}
	}
	// delete old occurs
	for (int i = 0; i < p_ol.size(); i++) p_ol[i]->set_status(DELETED);
	for (int j = 0; j < n_ol.size(); j++) n_ol[j]->set_status(DELETED);
	return true; // AND-substitution successful
}

bool substitute_OR(const uint32& x, const uVector1D& defs, OL& p_ol, OL& n_ol, uVector1D& out_c)
{
	// count number of added clauses & literals for positives
	int numAddedClauses = 0;
	for (int i = 0; i < p_ol.size(); i++)
		if (!isTautology_or(p_ol[i], defs))
			numAddedClauses++;
	// count number of added clauses & literals for negatives
	for (int d = 0; d < defs.size(); d++)
		for (int i = 0; i < n_ol.size(); i++)
			if (!(n_ol[i]->has(defs[d])))
				numAddedClauses++;
	if (numAddedClauses > p_ol.size() + n_ol.size()) return false;
	// substitute positives
	for (int i = 0; i < p_ol.size(); i++) {
		if (!isTautology_or(p_ol[i], defs)) {
			assert(out_c.empty());
			clause_extend_or(x, p_ol[i], defs, out_c);
			Sort(out_c.d_ptr(), out_c.size());
			S_REF added = new SCLAUSE();
			added->push(out_c[0]);
			for (int k = 1; k < out_c.size(); k++) {
				if (out_c[k - 1] == out_c[k]) continue; // skip repeated literal
				added->push(out_c[k]);
			}
			assert(added->isSorted());
			assert(added->hasZero() < 0);
			g_pFrost->attachClause(added);
			out_c.clear();
		}
	}
	// substitute negatives
	for (int d = 0; d < defs.size(); d++) {
		for (int i = 0; i < n_ol.size(); i++) {
			if (!(n_ol[i]->has(defs[d]))) {
				S_REF added = new SCLAUSE();
				clause_split_or(NEG(x), n_ol[i], defs[d], added);
				g_pFrost->attachClause(added);
			}
		}
	}
	// delete old occurs
	for (int i = 0; i < p_ol.size(); i++) p_ol[i]->set_status(DELETED);
	for (int j = 0; j < n_ol.size(); j++) n_ol[j]->set_status(DELETED);
	return true; // OR-substitution successful
}

bool gateReasoning_x(const uint32& x, OL& p_ol, OL& n_ol)
{
	uint32 imp;
	uint32 sig = 0;
	uVector1D out_c;
	// AND-gate Reasoning
	for (int i = 0; i < n_ol.size(); i++) {
		SCLAUSE& c = *n_ol[i];
		if (c.size() == 2) {
			if ((c[0] ^ x) == NEG_SIGN) { // found x with negative sign
				imp = FLIP(c[1]); // toggle implied literal sign
				out_c.push(imp);
				sig |= MAPHASH(imp);
			}
			else if ((c[1] ^ x) == NEG_SIGN) {
				imp = FLIP(c[0]); // toggle implied literal sign
				out_c.push(imp);
				sig |= MAPHASH(imp);
			}
		}
	}
	if (out_c.size() > 1) {
		out_c.push(x);
		sig |= MAPHASH(x);
		Sort(out_c.d_ptr(), out_c.size());
		for (int i = 0; i < p_ol.size(); i++) {
			SCLAUSE& c = *p_ol[i];
			if (c.size() == out_c.size() && subset_sig(c.sig(), sig) && isEqual(c, out_c)) {
				uVector1D defs(out_c.size() - 1);
				int nDefs = 0;
				for (LIT_POS k = 0; k < out_c.size(); k++)
					if (out_c[k] != x) defs[nDefs++] = FLIP(out_c[k]);
				assert(nDefs == out_c.size() - 1);
				out_c.clear();
				if (substitute_AND(x, defs, p_ol, n_ol, out_c)) {
					defs.clear(true);
					out_c.clear(true);
					return true;
				}
			}
		}
	}
	// OR-gate Reasoning
	out_c.clear();
	sig = 0;
	for (int i = 0; i < p_ol.size(); i++) {
		SCLAUSE& c = *p_ol[i];
		if (c.size() == 2) {
			if (c[0] == x) { // found x with positive sign
				imp = FLIP(c[1]); // toggle implied literal sign
				out_c.push(imp);
				sig |= MAPHASH(imp);
			}
			else if (c[1] == x) {
				imp = FLIP(c[0]); // toggle implied literal sign
				out_c.push(imp);
				sig |= MAPHASH(imp);
			}
		}
	}
	if (out_c.size() > 1) {
		out_c.push(FLIP(x));
		sig |= MAPHASH(FLIP(x));
		Sort(out_c.d_ptr(), out_c.size());
		for (int i = 0; i < n_ol.size(); i++) {
			SCLAUSE& c = *n_ol[i];
			if (c.size() == out_c.size() && subset_sig(c.sig(), sig) && isEqual(c, out_c)) {
				uVector1D defs(out_c.size() - 1);
				int nDefs = 0;
				for (LIT_POS k = 0; k < out_c.size(); k++)
					if (out_c[k] != FLIP(x)) defs[nDefs++] = out_c[k];
				assert(nDefs == out_c.size() - 1);
				out_c.clear();
				if (substitute_OR(x, defs, p_ol, n_ol, out_c)) {
					defs.clear(true);
					out_c.clear(true);
					return true;
				}
			}
		}
	}
	return false;
}

bool mayResolve_x(const uint32& x, OL& p_ol, OL& n_ol)
{
	assert(x);
	int numTautologies = 0;
	int64 numAddedLiterals = 0;
	// check resolvability
	for (int i = 0; i < p_ol.size(); i++)
		for (int j = 0; j < n_ol.size(); j++)
			if (isTautology(x, p_ol[i], n_ol[j]))
				numTautologies++;
			else
				numAddedLiterals += (int64)p_ol[i]->size() + (int64)n_ol[j]->size() - 2;
	int resolvents = p_ol.size() * n_ol.size() - numTautologies;
	if (resolvents == 0) { // No resolvents to add
		for (int i = 0; i < p_ol.size(); i++) p_ol[i]->set_status(DELETED);
		for (int j = 0; j < n_ol.size(); j++) n_ol[j]->set_status(DELETED);
		return true;
	}
	if (resolvents > p_ol.size() + n_ol.size()) return false;
	// count literals before elimination
	int64 lits_before = 0;
	for (int i = 0; i < p_ol.size(); i++) lits_before += (int64)p_ol[i]->size();
	for (int j = 0; j < n_ol.size(); j++) lits_before += (int64)n_ol[j]->size();
	if (numAddedLiterals > lits_before) return false;
	// can be eliminated
	resolve_x(x, p_ol, n_ol);
	return true; 
}

void self_sub_x(const uint32& x, OL& p_ol, OL& n_ol)
{
	// positives vs negatives
	for (int i = 0; i < p_ol.size(); i++) {
		S_REF pos_c = p_ol[i];
		if (pos_c->status() == DELETED) continue;
		// self-subsumption check
		for (int j = 0; j < n_ol.size(); j++) {
			S_REF neg_c = n_ol[j];
			if (neg_c->status() == DELETED) continue;
			if (neg_c->size() > 1 && neg_c->size() < pos_c->size() &&
				selfSubset_sig(neg_c->sig(), pos_c->sig()) && selfSubset(x, neg_c, pos_c)) {
				g_pFrost->strengthen(pos_c, x);
				pos_c->melt(); // mark for fast recongnition in ot update 
			}
			else if (neg_c->size() > 1 && neg_c->size() == pos_c->size() &&
				selfSubset_sig(neg_c->sig(), pos_c->sig()) && selfSubset(x, neg_c, pos_c)) {
				// shorten the positive occur and delete the negative occur subsumed by the former
				g_pFrost->strengthen(pos_c, x);
				pos_c->melt(); 
				neg_c->set_status(DELETED);
			}
		}
		// subsumption check
		for (int j = 0; j < p_ol.size(); j++) {
			S_REF sm_c = p_ol[j];
			if (sm_c->status() != DELETED && sm_c->size() < pos_c->size() &&
				subset_sig(sm_c->sig(), pos_c->sig()) && subset(sm_c, pos_c)) {
				pos_c->set_status(DELETED);
				break; // mission accomplished!
			}
		}
	}
	updateOL(p_ol);
	// negatives vs positives
	for (int i = 0; i < n_ol.size(); i++) {
		S_REF neg_c = n_ol[i];
		if (neg_c->status() == DELETED) continue;
		// self-subsumption check
		for (int j = 0; j < p_ol.size(); j++) {
			S_REF pos_c = p_ol[j];
			if (pos_c->status() == DELETED) continue;
			if (pos_c->size() > 1 && pos_c->size() < neg_c->size() &&
				selfSubset_sig(pos_c->sig(), neg_c->sig()) && selfSubset(NEG(x), pos_c, neg_c)) {
				g_pFrost->strengthen(neg_c, NEG(x));
				neg_c->melt(); 
			}
		}
		// subsumption check
		for (int j = 0; j < n_ol.size(); j++) {
			S_REF sm_c = n_ol[j];
			if (sm_c->status() != DELETED && sm_c->size() < neg_c->size() &&
				subset_sig(sm_c->sig(), neg_c->sig()) && subset(sm_c, neg_c)) {
				neg_c->set_status(DELETED);
				break; // mission accomplished!
			}
		}
	}
	updateOL(n_ol);
}

void blocked_x(const uint32& x, OL& p_ol, OL& n_ol)
{
	uint32 numTautologies = 0;
	if (p_ol.size() <= n_ol.size()) {  // start with positives
		for (int i = 0; i < p_ol.size(); i++) {
			numTautologies = 0;
			if (p_ol[i]->status() == DELETED) continue;
			for (int j = 0; j < n_ol.size(); j++) {
				if (n_ol[j]->status() == DELETED) continue;
				if (isTautology(x, p_ol[i], n_ol[j])) numTautologies++;
			}
			if (numTautologies == n_ol.size()) p_ol[i]->set_status(DELETED);
		}
	}
	else { // start with negatives
		for (int i = 0; i < n_ol.size(); i++) {
			numTautologies = 0;
			if (n_ol[i]->status() == DELETED) continue;
			for (int j = 0; j < p_ol.size(); j++) {
				if (p_ol[j]->status() == DELETED) continue;
				if (isTautology(x, n_ol[i], p_ol[j])) numTautologies++;
			}
			if (numTautologies == p_ol.size()) n_ol[i]->set_status(DELETED);
		}
	}
}

#endif