/***********************************************************************[equivalence.h]
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

#ifndef __EQU_
#define __EQU_

#include "simplify.h" 
using namespace pFROST;

inline uint32 substitute_single(const uint32& dx, SCLAUSE& org, const uint32& def)
{
	CHECKLIT(dx);
	assert(def != dx);
	assert(org.original());
	PFLCLAUSE(4, org, " Clause ");
	int n = 0;
	for (int i = 0; i < org.size(); i++) {
		const uint32 lit = org[i];
		if (lit == dx) org[n++] = def;
		else if (NEQUAL(lit, def)) org[n++] = lit;
	}
	org.resize(n);
	if (org.size() == 1) return *org;
	Sort(org.data(), org.size(), LESS<uint32>());
	org.calcSig();
	PFLCLAUSE(4, org, " Substituted to ");
	return 0;
}

inline bool substitute_single(const uint32& p, const uint32& def, OT& ot)
{
	CHECKLIT(def);
	assert(!SIGN(p));
	const uint32 n = NEG(p), def_f = FLIP(def);
	OL& poss = ot[p], & negs = ot[n];
	const bool proofEN = pfrost->opts.proof_en;
	// substitute negatives 
	for (int i = 0; i < negs.size(); i++) {
		SCLAUSE& neg = *negs[i];
		if (neg.learnt() || neg.molten() || neg.has(def))
			neg.markDeleted();
		else if (neg.original()) {
			uint32 unit = substitute_single(n, neg, def_f);
			if (unit) {
				const LIT_ST val = pfrost->litvalue(unit);
				if (UNASSIGNED(val))
					pfrost->enqueueUnit(unit);
				else if (!val) return true;
			}
			else if (proofEN)
				pfrost->proof.addResolvent(neg);
		}
	}
	// substitute positives
	for (int i = 0; i < poss.size(); i++) {
		SCLAUSE& pos = *poss[i];
		if (pos.learnt() || pos.molten() || pos.has(def_f))
			pos.markDeleted();
		else if (pos.original()) {
			uint32 unit = substitute_single(p, pos, def);
			if (unit) {
				const LIT_ST val = pfrost->litvalue(unit);
				if (UNASSIGNED(val))
					pfrost->enqueueUnit(unit);
				else if (!val) return true;
			}
			else if (proofEN)
				pfrost->proof.addResolvent(pos);
		}
	}
	return false; 
}

inline uint32 find_sfanin(const uint32& gate_out, OL& list)
{
	CHECKLIT(gate_out);
	uint32 imp = 0;
	int nImps = 0;
	forall_occurs(list, i) {
		SCLAUSE& c = **i;
		if (c.original() && c.size() == 2) {
			imp = FLIP(c[0] ^ c[1] ^ gate_out);
			c.melt(); // mark as gate clause
			nImps++;
		}
		if (nImps > 1) {
			// cannot be a single-input gate	
			return 0;
		}
	}
	return imp;
}

inline uint32 find_BN_gate(const uint32& p, OL& poss, OL& negs)
{
	CHECKLIT(p);
	assert(!SIGN(p));
	assert(checkMolten(poss, negs));
	const uint32 n = NEG(p);
	uint32 first = find_sfanin(p, poss);
	if (first) {
		uint32 second = n, def = first;
		if (second < first) first = second, second = def;
		for (int i = 0; i < negs.size(); i++) {
			SCLAUSE& c = *negs[i];
			if (c.original() && c.size() == 2 && c[0] == first && c[1] == second) {
				c.melt(); // mark as fanout clause
				PFLOG2(3, " Gate %d = -/+%d found", ABS(p), ABS(def));
				PFLOCCURS(pfrost, 4, ABS(p));
				return def;
			}
		}
	}
	freezeBinaries(poss);
	return 0;
}

inline void save_BN_gate(const uint32& p, const int& pOrgs, const int& nOrgs, OL& poss, OL& negs, MODEL& model)
{
	CHECKLIT(p);
	PFLOG2(4, " saving buffer/inverter clauses as witness");
	const uint32 n = NEG(p);
	if (pOrgs > nOrgs) {
		for (int i = 0; i < negs.size(); i++) {
			SCLAUSE& c = *negs[i];
			if (c.original())
				model.saveClause(c, c.size(), n);
		}
		model.saveWitness(p);
	}
	else {
		for (int i = 0; i < poss.size(); i++) {
			SCLAUSE& c = *poss[i];
			if (c.original())
				model.saveClause(c, c.size(), p);
		}
		model.saveWitness(n);
	}
}

#endif