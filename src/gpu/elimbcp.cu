/***********************************************************************[elimbcp.cu]
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

#include "solver.hpp"

using namespace ParaFROST;

void Solver::createOTHost(HOT& hot)
{
	assert(hcnf != NULL);
	assert(vars->nUnits);
	assert(reallocFailed());
	for (uint32 i = 0; i < hcnf->size(); i++) {
		const S_REF r = hcnf->ref(i);
		SCLAUSE& c = (*hcnf)[r];
		if (c.learnt() || c.original()) {
			assert(c.size());
			forall_clause(c, k) {
				CHECKLIT(*k);
				hot[*k].push(r);
			}
		}
	}
}

bool Solver::propFailed()
{
	if (vars->nUnits) {
		SYNCALL; // sync 'cacheCNF'
		assert(reallocFailed());
		assert(hcnf != NULL);
		assert(inf.nClauses == hcnf->size());
		HOT hot(inf.nDualVars);
		createOTHost(hot);
		// start proping on host
		nForced = sp->propagated;
		assert(vars->cachedUnits != NULL);
		LIT_ST* values = sp->value;
		uint32* t = vars->cachedUnits + vars->nUnits;
		for (uint32* u = vars->cachedUnits; u != t; u++) {
			const uint32 unit = *u;
			CHECKLIT(unit);
			const LIT_ST val = values[unit];
			if (UNASSIGNED(val)) enqueueDevUnit(unit);
			else if (!val) return false; // early conflict detection
		}
		if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
		else LOGN2(2, " Propagating pending forced units..");
		int64 bclauses = inf.nClauses, bliterals = inf.nLiterals;
		CNF& hcnf = *this->hcnf;
		while (sp->propagated < trail.size()) { // propagate units
			uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
			CHECKLIT(assign);
			forall_occurs(hot[assign], i) {
				hcnf[*i].markDeleted();
			}
			forall_occurs(hot[f_assign], i) {
				SCLAUSE& c = hcnf[*i];
				assert(c.size());
				if (c.deleted()) continue;
				if (propClause(values, f_assign, c))
					c.markDeleted();
				else {
					assert(c.size()); // cannot be empty at this point
					if (c.size() == 1) {
						const uint32 unit = *c;
						CHECKLIT(unit);
						if (UNASSIGNED(values[unit])) enqueueUnit(unit);
						else { learnEmpty(); return false; }
					}
				}
			}
		}
		LOGDONE(2, 5);
		nForced = sp->propagated - nForced;
		LOGREDALLHOST(this, 2, "BCP Reductions");
		countAll(1);
		inf.nClauses = inf.n_cls_after;
		inf.nLiterals = inf.n_lits_after;
		stats.units.forced += nForced;
		stats.sigma.all.clauses += bclauses - int64(inf.nClauses);
		stats.sigma.all.literals += bliterals - int64(inf.nLiterals);
		vars->nUnits = nForced = 0;
	}
	return true;
}

bool Solver::prop()
{
	if (!enqueueCached(streams[3])) { learnEmpty(); return false; }
	LIT_ST* values = sp->value;
	while (sp->propagated < trail.size()) {
		uint32 assign = trail[sp->propagated++], f_assign = FLIP(assign);
		CHECKLIT(assign);
		OL& ol = (*ot)[assign];
		OL& f_ol = (*ot)[f_assign];
		CNF& cnf = *this->cnf;
		forall_occurs(f_ol, i) { // reduce unsatisfied 
			SCLAUSE& c = cnf[*i];
			assert(c.size());
			if (c.deleted()) continue;
			if (propClause(values, f_assign, c))
				c.markDeleted();
			else {
				assert(c.size()); // cannot be empty at this point
				if (c.size() == 1) {
					const uint32 unit = *c;
					CHECKLIT(unit);
					if (UNASSIGNED(values[unit])) enqueueUnit(unit);
					else { learnEmpty(); return false; }
				}
			}
		}
		forall_occurs(ol, i) { 
			SCLAUSE& c = cnf[*i];
			assert(c.size());
			c.markDeleted();
		}
		ol.clear(true), f_ol.clear(true);
	}
	cleanProped();
	return true;
}

inline bool Solver::propClause(const LIT_ST* values, const uint32& lit, SCLAUSE& c)
{
	assert(c.size() > 1);
	uint32 sig = 0;
	uint32* j = c;
	forall_clause(c, i) {
		const uint32 other = *i;
		if (NEQUAL(other, lit)) {
			if (values[other] > 0) return true;
			*j++ = other;
			sig |= MAPHASH(other);
		}
	}
	assert(int(j - c) == c.size() - 1);
	assert(c.hasZero() < 0);
	c.set_sig(sig);
	c.pop();
	assert(c.isSorted());
	return false;
}

inline bool Solver::enqueueCached(const cudaStream_t& stream) {
	if (vars->nUnits) {
		nForced = sp->propagated;
		SYNC(stream); // sync units copy
		assert(vars->cachedUnits != NULL);
		LIT_ST* values = sp->value;
		uint32* t = vars->cachedUnits + vars->nUnits;
		for (uint32* u = vars->cachedUnits; u != t; u++) {
			const uint32 unit = *u;
			CHECKLIT(unit);
			const LIT_ST val = values[unit];
			if (UNASSIGNED(val)) enqueueDevUnit(unit);
			else if (!val) return false; // early conflict detection
		}
		if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
		else LOGN2(2, " Propagating forced units..");
		SYNCALL; // sync ot creation
	}
	return true;
}

inline void	Solver::cleanProped() {
	if (vars->nUnits) {
		LOGDONE(2, 5);
		nForced = sp->propagated - nForced;
		LOGREDALL(this, 2, "BCP Reductions");
		nForced = 0, vars->tmpUnits.clear();
		assert(vars->tmpUnits.data() == cumm.unitsdPtr());
		if (!opts.sub_en) reduceOTAsync(cnf, ot, 0);
		CHECK(cudaMemcpyAsync(vars->units, &vars->tmpUnits, sizeof(cuVecU), cudaMemcpyHostToDevice));
	}
}