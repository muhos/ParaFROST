/***********************************************************************[elimination.cpp]
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

#include "and.h"
#include "xor.h"
#include "equivalence.h"
#include "ifthenelse.h"
#include "redundancy.h"
#include "subsume.h" 

using namespace ParaFROST;

inline void Solver::bumpShrunken(SCLAUSE& c)
{
	assert(c.learnt());
	assert(c.size() > 1);
	const int old_lbd = c.lbd();
	if (old_lbd <= opts.lbd_tier1) return; // always keep Tier1 value
	int new_lbd = MIN(c.size() - 1, old_lbd);
	if (new_lbd >= old_lbd) return;
	c.set_lbd(new_lbd);
	c.set_usage(USAGET3);
	PFLCLAUSE(4, c, " Bumping shrunken clause with (lbd:%d, usage:%d) ", new_lbd, c.usage());
}

inline bool Solver::propClause(const LIT_ST* values, const uint32& lit, SCLAUSE& c)
{
	uint32 sig = 0;
	uint32* j = c;
	forall_clause(c, i) {
		const uint32 other = *i;
		if (NEQUAL(other, lit)) {
			if (values[other] > 0) return true;
			*j++ = other;
			sig |= MAPHASH(other);
		}
		else
			assert(!values[other]);
		
	}
	assert(int(j - c) == c.size() - 1);
	assert(c.hasZero() < 0);
	c.set_sig(sig);
	c.pop();
	assert(c.isSorted());
	return false;
}

bool Solver::prop()
{
	nForced = sp->propagated;
	LIT_ST* values = sp->value;
	while (sp->propagated < trail.size()) { 
		const uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		CHECKLIT(assign);
		// reduce unsatisfied
		forall_occurs(ot[f_assign], i) {
			SCLAUSE& c = **i;
			assert(c.size());
			if (c.deleted()) continue;
			if (propClause(values, f_assign, c))
				c.markDeleted(); // clause satisfied by an assigned unit
			else {
				const int size = c.size();
				if (!size) { learnEmpty(); return false; }
				if (size == 1) {
					const uint32 unit = *c;
					CHECKLIT(unit);
					if (UNASSIGNED(values[unit])) enqueueUnit(unit);
					else { learnEmpty(); return false; }
				}
			}
		}
		ot[f_assign].clear(true);
		toblivion(ot[assign]);
	}
	nForced = sp->propagated - nForced;
	if (nForced) {
		PFLREDALL(this, 2, "BCP Reductions");
		nForced = 0;
		if (!opts.sub_en) reduceOT();
	}
	return true;
}

void Solver::VE()
{
	if (opts.ve_en) {
		PFLOG2(2, " Eliminating variables..");
		bve();
		PFLREDALL(this, 2, "VE Reductions");
	}
}

void Solver::SUB()
{
	if (opts.sub_en || opts.ve_plus_en) {
		if (interrupted()) killSolver();
		PFLOG2(2, " Eliminating (self)-subsumptions..");
		if (opts.profile_simp) timer.pstart();
		SUBSTATS& substats = stats.sigma.sub;
		for (uint32 i = 0; i < PVs.size(); i++) {
			const uint32 v = PVs[i];
			assert(v);
			assert(!sp->vstate[v].state);
			const uint32 p = V2L(v), n = NEG(p);
			OL& poss = ot[p], &negs = ot[n];
			if (poss.size() <= opts.sub_limit && negs.size() <= opts.sub_limit)
				self_sub_x(p, poss, negs, substats);
		}
		if (opts.profile_simp) timer.pstop(), timer.sub += timer.pcpuTime();
		PFLREDALL(this, 2, "SUB Reductions");
	}
}

void Solver::BCE()
{
	if (opts.bce_en) {
		if (interrupted()) killSolver();
		PFLOG2(2, " Eliminating blocked clauses..");
		if (opts.profile_simp) timer.pstart();
		for (uint32 i = 0; i < PVs.size(); i++) {
			const uint32 v = PVs[i];
			if (!v) continue;
			const uint32 p = V2L(v), n = NEG(p);
			OL& poss = ot[p], &negs = ot[n];
			if (poss.size() <= opts.bce_limit && negs.size() <= opts.bce_limit) {
				// start with negs
				for (int i = 0; i < negs.size(); i++) {
					SCLAUSE& neg = *negs[i];
					if (neg.original()) {
						bool allTautology = true;
						for (int j = 0; j < poss.size(); j++) {
							SCLAUSE& pos = *poss[j];
							if (pos.original() && !isTautology(v, neg, pos)) {
								allTautology = false;
								break;
							}
						}
						if (allTautology) {
							assert(neg.original());
							model.saveClause(neg, neg.size(), n);
							negs[i]->markDeleted();
						}
					}
				}
			}
		}
		if (opts.profile_simp) timer.pstop(), timer.bce += timer.pcpuTime();
		PFLREDALL(this, 2, "BCE Reductions");
	}
}

void Solver::ERE()
{
	if (!opts.ere_en) return;
	if (interrupted()) killSolver();
	PFLOG2(2, " Eliminating redundances..");
	if (opts.profile_simp) timer.pstart();
	ERESTATS& erestats = stats.sigma.ere;
	const int maxsize = opts.ere_max_resolvent;
	for (uint32 n = 0; n < PVs.size(); n++) {
		assert(PVs[n]);
		uint32 dx = V2L(PVs[n]), fx = NEG(dx);
		if (ot[dx].size() > ot[fx].size()) swap(dx, fx);
		OL& me = ot[dx], &other = ot[fx];
		if (me.size() <= opts.ere_limit && other.size() <= opts.ere_limit) {
			// do merging and apply forward equality check (on-the-fly) over resolvents
			for (int i = 0; i < me.size(); i++) {
				const SCLAUSE& ci = *me[i];
				if (ci.deleted()) continue;
				for (int j = 0; j < other.size(); j++) {
					const SCLAUSE& cj = *other[j];
					if (cj.deleted()) continue;
					forward_equ(PVs[n], ci, cj, ot, maxsize, erestats);
				}
			}
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.ere += timer.pcpuTime();
	PFLREDCL(this, 2, "ERE Reductions");
}

void Solver::strengthen(SCLAUSE& c, const uint32& me)
{
	uint32 sig = 0;
	int n = 0;
	for (int k = 0; k < c.size(); k++) {
		const uint32 lit = c[k];
		if (NEQUAL(lit, me)) {
			c[n++] = lit;
			sig |= MAPHASH(lit);
		}
	}
	assert(n == c.size() - 1);
	assert(c.hasZero() < 0);
	c.set_sig(sig);
	c.pop();
	if (n == 1) {
		const uint32 unit = *c;
		const LIT_ST val = sp->value[unit];
		if (UNASSIGNED(val)) {
			enqueueUnit(unit);
			toblivion(ot[unit]); // safe to remove clauses of 'unit' in 'SUB'
		}
		else if (!val) { 
			PFLOG2(2, "  SUB proved a contradiction");
			learnEmpty();
			killSolver();
		}
	}
	else {
		assert(c.isSorted());
		if (opts.proof_en) 
			proof.addResolvent(c);
		if (c.learnt()) 
			bumpShrunken(c);
	}
}