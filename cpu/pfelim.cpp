/***********************************************************************[pfelim.cpp]
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

#include "pfsimp.h"

using namespace pFROST;
using namespace SIGmA;

bool ParaFROST::prop()
{
	nForced = sp->propagated;
	while (sp->propagated < trail.size()) { // propagate units
		uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		assert(assign > 1);
		// remove satisfied
		for (int i = 0; i < ot[assign].size(); i++) ot[assign][i]->markDeleted();
		// reduce unsatisfied
		for (int i = 0; i < ot[f_assign].size(); i++) {
			S_REF c = ot[f_assign][i];
			//c->print();
			assert(c->size());
			if (c->deleted() || propClause(c, f_assign)) continue; // clause satisfied
			// clause is unit or conflict
			// Note: attach & strengthen don't check for conflict before enqueue
			if (c->size() == 0) { cnfstate = UNSAT; return false; }
			if (c->size() == 1) {
				assert(**c > 1);
				if (unassigned(**c)) enqueueOrg(**c);
				else { cnfstate = UNSAT; return false; }
			}
		}
		// delete assign lists
		ot[assign].clear(true), ot[f_assign].clear(true);
	}
	nForced = sp->propagated - nForced;
	if (nForced) {
		PFLREDALL(this, 2, "BCP Reductions");
		nForced = 0;
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
		assert(sp->vstate[v] == ACTIVE);
		uint32 p = V2L(v), n = NEG(p);
		OL& poss = ot[p], & negs = ot[n];
		int pOrgs = 0, nOrgs = 0;
		countOrgs(poss, pOrgs), countOrgs(negs, nOrgs);
		// pure-literal
		if (!pOrgs || !nOrgs) {
			toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
			sp->vstate[v] = MELTED, v = 0;
		}
		else {
			assert(pOrgs && nOrgs);
			out_c.clear();
			uint32 def;
			// Equiv/NOT-gate Reasoning
			if (def = find_BN_gate(p, poss, negs)) {
				saveResolved(p, pOrgs, nOrgs, poss, negs, model.resolved);
				substitute_single(p, def, poss, negs);
				sp->vstate[v] = MELTED, v = 0;
			}
			// AND-gate Reasoning
			else if (find_AO_gate(n, pOrgs + nOrgs, ot, out_c)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
			// OR-gate Reasoning
			else if (find_AO_gate(p, pOrgs + nOrgs, ot, out_c)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
			// ITE-gate Reasoning
			else if (find_ITE_gate(p, pOrgs + nOrgs, ot, out_c)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
			else if (find_ITE_gate(n, pOrgs + nOrgs, ot, out_c)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
			// XOR-gate Reasoning
			else if (find_XOR_gate(p, pOrgs + nOrgs, ot, out_c, false)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
			// n-by-m resolution
			else if (resolve_x(v, pOrgs, nOrgs, poss, negs, out_c, false)) {
				toblivion(p, pOrgs, nOrgs, poss, negs, model.resolved);
				sp->vstate[v] = MELTED, v = 0;
			}
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

void ParaFROST::HSE()
{
	if (opts.hse_en || opts.ve_plus_en) {
		if (interrupted()) killSolver();
		PFLOGN2(2, "  Eliminating (self)-subsumptions..");
		if (opts.profile_simp) timer.pstart();
		for (uint32 i = 0; i < PVs.size(); i++) {
			uint32 v = PVs[i];
			assert(v);
			assert(sp->vstate[v] == ACTIVE);
			uint32 p = V2L(v), n = NEG(p);
			if (ot[p].size() <= opts.hse_limit && ot[n].size() <= opts.hse_limit)
				self_sub_x(p, ot[p], ot[n]);
		}
		if (opts.profile_simp) timer.pstop(), timer.hse += timer.pcpuTime();
		PFLDONE(2, 5);
		PFLREDALL(this, 2, "HSE Reductions");
	}
}

void ParaFROST::BCE()
{
	if (opts.bce_en) {
		if (interrupted()) killSolver();
		PFLOGN2(2, " Eliminating blocked clauses..");
		if (opts.profile_simp) timer.pstart();
		for (uint32 i = 0; i < PVs.size(); i++) {
			uint32 v = PVs[i];
			if (!v) continue;
			uint32 p = V2L(v), n = NEG(p);
			if (ot[p].size() <= opts.bce_limit && ot[n].size() <= opts.bce_limit)
				blocked_x(v, ot[p], ot[n]);
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
		uint32 p = V2L(PVs[n]);
		OL& poss = ot[p], & negs = ot[NEG(p)];
		if (ot[p].size() <= opts.ere_limit && ot[n].size() <= opts.ere_limit) {
			// do merging and apply forward equality check (on-the-fly) over resolvents
			for (int i = 0; i < poss.size(); i++) {
				if (poss[i]->deleted()) continue;
				for (int j = 0; j < negs.size(); j++) {
					if (negs[j]->deleted() || (poss[i]->size() + negs[j]->size() - 2) > MAX_ERE_OUT) continue;
					if (merge_ere(PVs[n], poss[i], negs[j], m_c)) {
						CL_ST type;
						if (poss[i]->learnt() || negs[j]->learnt()) type = LEARNT;
						else type = ORIGINAL;
						if (m_c.size() > 1) forward_equ(m_c, ot, type);
					}
				}
			}
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.ere += timer.pcpuTime();
	PFLDONE(2, 5);
	PFLREDCL(this, 2, "ERE Reductions");
}

bool ParaFROST::propClause(S_REF c, const uint32& f_assign)
{
	uint32 sig = 0;
	int n = 0;
	bool check = false;
	for (int k = 0; k < c->size(); k++) {
		uint32 lit = c->lit(k);
		if (lit != f_assign) {
			if (isTrue(lit)) return true;
			(*c)[n++] = lit;
			sig |= MAPHASH(lit);
		}
		else check = true;
	}
	assert(check);
	assert(n == c->size() - 1);
	assert(c->hasZero() < 0);
	assert(c->isSorted());
	c->set_sig(sig);
	c->pop();
	return false;
}

void ParaFROST::strengthen(S_REF c, const uint32& self_lit)
{
	uint32 sig = 0;
	int n = 0;
	bool check = false;
	for (int k = 0; k < c->size(); k++) {
		uint32 lit = c->lit(k);
		if (lit != self_lit) {
			(*c)[n++] = lit;
			sig |= MAPHASH(lit);
		}
		else check = true;
	}
	assert(check);
	assert(n == c->size() - 1);
	assert(c->hasZero() < 0);
	assert(c->isSorted());
	c->set_sig(sig);
	c->pop();
	if (opts.proof_en) {
		wrProof('a');
		wrProof(*c, n);
		wrProof(0);
	}
	if (n == 1 && unassigned(**c)) enqueueOrg(**c);
	else if (n > 1 && c->learnt()) bumpShrunken(c);
}

inline void ParaFROST::bumpShrunken(S_REF c)
{
	assert(c->learnt());
	assert(c->size() > 1);
	int old_lbd = c->lbd();
	if (old_lbd <= opts.lbd_tier1) return; // always keep Tier1 value
	int new_lbd = std::min(c->size() - 1, old_lbd);
	if (new_lbd >= old_lbd) return;
	c->set_lbd(new_lbd);
	c->set_usage(USAGET3);
	PFLCLAUSE(4, (*c), " Bumping shrunken clause with (lbd:%d, usage:%d) ", new_lbd, c->usage());
}