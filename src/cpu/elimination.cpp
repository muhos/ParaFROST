/***********************************************************************[elim.cpp]
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

#include "equivalence.h"
#include "ifthenelse.h"
#include "redundancy.h"
#include "and.h"
#include "xor.h"

using namespace pFROST;
using namespace SIGmA;

bool ParaFROST::propClause(S_REF c, const uint32& me)
{
	uint32 sig = 0;
	int n = 0;
	for (int k = 0; k < c->size(); k++) {
		const uint32 lit = c->lit(k);
		if (lit != me) {
			if (isTrue(lit)) return true;
			(*c)[n++] = lit;
			sig |= MAPHASH(lit);
		}
	}
	assert(n == c->size() - 1);
	assert(c->hasZero() < 0);
	assert(c->isSorted());
	c->set_sig(sig);
	c->pop();
	return false;
}

bool ParaFROST::prop()
{
	uint32 propagated = nForced = sp->propagated;
	while (propagated < trail.size()) {
		const uint32 assign = trail[propagated++];
		CHECKLIT(assign);
		for (int i = 0; i < ot[assign].size(); i++)
			ot[assign][i]->markDeleted();
		ot[assign].clear(true);
	}
	while (sp->propagated < trail.size()) { 
		const uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		CHECKLIT(assign);
		// reduce unsatisfied
		for (int i = 0; i < ot[f_assign].size(); i++) {
			S_REF c = ot[f_assign][i];
			assert(c->size());
			if (c->deleted()) continue;
			if (propClause(c, f_assign))
				removeClause(c); // clause satisfied by assigned unit
			else {
				if (c->empty()) { learnEmpty(); return false; }
				if (c->size() == 1) {
					const uint32 unit = **c;
					CHECKLIT(unit);
					if (unassigned(unit)) enqueue(unit);
					else { learnEmpty(); return false; }
				}
			}
		}
		// remove satisfied
		for (int i = 0; i < ot[assign].size(); i++)
			removeClause(ot[assign][i]);
		ot[assign].clear(true), ot[f_assign].clear(true);
	}
	nForced = sp->propagated - nForced;
	if (nForced) {
		PFLREDALL(this, 2, "BCP Reductions");
		nForced = 0;
		if (!opts.sub_en) reduceOT();
	}
	return true;
}

void ParaFROST::bve()
{
	if (interrupted()) killSolver();
	if (opts.profile_simp) timer.pstart();
	Lits_t out_c;
	out_c.reserve(INIT_CAP);
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32& v = PVs[i];
		assert(v);
		assert(!sp->vstate[v].state);
		const uint32 p = V2L(v), n = NEG(p);
		OL& poss = ot[p], & negs = ot[n];
		int pOrgs = 0, nOrgs = 0;
		countOrgs(poss, pOrgs), countOrgs(negs, nOrgs);
		// pure-literal
		if (!pOrgs || !nOrgs) {
			toblivion(p, pOrgs, nOrgs, poss, negs, model);
			sp->vstate[v].state = MELTED_M, v = 0;
		}
		// Equiv/NOT-gate Reasoning
		else if (uint32 def = find_BN_gate(p, poss, negs)) {
			save_BN_gate(p, pOrgs, nOrgs, poss, negs, model);
			if (substitute_single(p, def, ot)) {
				learnEmpty();
				killSolver();
			}
			sp->vstate[v].state = MELTED_M, v = 0;
		}
		else {
			assert(pOrgs && nOrgs);
			out_c.clear();
			const int nOrgCls = pOrgs + nOrgs;
			int nAddedCls = 0;
			// AND/OR-gate Reasoning
			if (find_AO_gate(n, nOrgCls, ot, out_c, nAddedCls)) {
				if (nAddedCls) substitute_x(v, negs, poss, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
			else if (find_AO_gate(p, nOrgCls, ot, out_c, nAddedCls)) {
				if (nAddedCls) substitute_x(v, poss, negs, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
			// ITE-gate Reasoning
			else if (find_ITE_gate(p, nOrgCls, ot, out_c, nAddedCls)) {
				if (nAddedCls) substitute_x(v, poss, negs, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
			else if (find_ITE_gate(n, nOrgCls, ot, out_c, nAddedCls)) {
				if (nAddedCls) substitute_x(v, negs, poss, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
			// XOR-gate Reasoning
			else if (find_XOR_gate(p, nOrgCls, ot, out_c, nAddedCls)) {
				if (nAddedCls) substitute_x(v, poss, negs, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
			// n-by-m resolution
			else if (!nAddedCls && mayResolve_x(v, nOrgCls, poss, negs, nAddedCls)) {
				if (nAddedCls) resolve_x(v, poss, negs, out_c);
				toblivion(p, pOrgs, nOrgs, poss, negs, model);
				sp->vstate[v].state = MELTED_M, v = 0;
			}
		}
		if (!v) {
			assert(inf.unassigned);
			inf.unassigned--;
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.ve += timer.pcpuTime();
}

void ParaFROST::VE()
{
	if (opts.ve_en) {
		PFLOGN2(2, "  Eliminating variables..");
		bve();
		PFLDONE(2, 5);
		PFLREDALL(this, 2, "VE Reductions");
	}
}

void ParaFROST::SUB()
{
	if (opts.sub_en || opts.ve_plus_en) {
		if (interrupted()) killSolver();
		PFLOGN2(2, "  Eliminating (self)-subsumptions..");
		if (opts.profile_simp) timer.pstart();
		for (uint32 i = 0; i < PVs.size(); i++) {
			const uint32 v = PVs[i];
			assert(v);
			assert(!sp->vstate[v].state);
			uint32 p = V2L(v), n = NEG(p);
			if (ot[p].size() <= opts.sub_limit && ot[n].size() <= opts.sub_limit)
				self_sub_x(p, ot[p], ot[n]);
		}
		if (opts.profile_simp) timer.pstop(), timer.sub += timer.pcpuTime();
		PFLDONE(2, 5);
		PFLREDALL(this, 2, "SUB Reductions");
	}
}

void ParaFROST::BCE()
{
	if (opts.bce_en) {
		if (interrupted()) killSolver();
		PFLOGN2(2, " Eliminating blocked clauses..");
		if (opts.profile_simp) timer.pstart();
		for (uint32 i = 0; i < PVs.size(); i++) {
			const uint32 v = PVs[i];
			if (!v) continue;
			uint32 p = V2L(v), n = NEG(p);
			if (ot[p].size() <= opts.bce_limit && ot[n].size() <= opts.bce_limit)
				blocked_x(v, ot[p], ot[n], model);
		}
		if (opts.profile_simp) timer.pstop(), timer.bce += timer.pcpuTime();
		PFLDONE(2, 5);
		PFLREDALL(this, 2, "BCE Reductions");
	}
}

void ParaFROST::ERE()
{
	if (!opts.ere_en) return;
	if (interrupted()) killSolver();
	PFLOGN2(2, " Eliminating redundances..");
	if (opts.profile_simp) timer.pstart();
	Lits_t m_c;
	m_c.reserve(INIT_CAP);
	for (uint32 n = 0; n < PVs.size(); n++) {
		assert(PVs[n]);
		const uint32 p = V2L(PVs[n]);
		OL& poss = ot[p], & negs = ot[NEG(p)];
		if (ot[p].size() <= opts.ere_limit && ot[n].size() <= opts.ere_limit) {
			// do merging and apply forward equality check (on-the-fly) over resolvents
			for (int i = 0; i < poss.size(); i++) {
				if (poss[i]->deleted()) continue;
				for (int j = 0; j < negs.size(); j++) {
					if (negs[j]->deleted() || (poss[i]->size() + negs[j]->size() - 2) > MAX_ERE_OUT) continue;
					if (merge_ere(PVs[n], poss[i], negs[j], m_c)) {
						if (m_c.size() > 1) {
							CL_ST type;
							if (poss[i]->learnt() || negs[j]->learnt()) type = LEARNT;
							else type = ORIGINAL;
							forward_equ(m_c, ot, type);
						}
					}
				}
			}
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.ere += timer.pcpuTime();
	PFLDONE(2, 5);
	PFLREDCL(this, 2, "ERE Reductions");
}

void ParaFROST::newResolvent(S_REF s)
{
	assert(*s != NULL);
	assert(s->size());
	assert(s->hasZero() < 0);
	assert(s->original());
	assert(s->isSorted());
	// unit clauses must be added to be checked for conflicts later
	if (s->size() == 1) {
		const uint32 unit = **s;
		enqueueUnit(unit);
	}
	else if (opts.proof_en) proof.addResolvent(*s);
	s->calcSig();
	s->markAdded();
	scnf.push(s); 
}

void ParaFROST::strengthen(S_REF c, const uint32& me)
{
	if (opts.proof_en) proof.shrinkClause(*c, me);
	uint32 sig = 0;
	int n = 0;
	for (int k = 0; k < c->size(); k++) {
		const uint32 lit = c->lit(k);
		if (NEQUAL(lit, me)) {
			(*c)[n++] = lit;
			sig |= MAPHASH(lit);
		}
	}
	assert(n == c->size() - 1);
	assert(c->hasZero() < 0);
	assert(c->isSorted());
	c->set_sig(sig);
	c->pop();
	if (n == 1) {
		const uint32 unit = **c;
		enqueueUnit(unit);
	}
	else if (c->learnt()) bumpShrunken(c);
}

inline void ParaFROST::bumpShrunken(S_REF c)
{
	assert(c->learnt());
	assert(c->size() > 1);
	const int old_lbd = c->lbd();
	if (old_lbd <= opts.lbd_tier1) return; // always keep Tier1 value
	int new_lbd = std::min(c->size() - 1, old_lbd);
	if (new_lbd >= old_lbd) return;
	c->set_lbd(new_lbd);
	c->set_usage(USAGET3);
	PFLCLAUSE(4, (*c), " Bumping shrunken clause with (lbd:%d, usage:%d) ", new_lbd, c->usage());
}