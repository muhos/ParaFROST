/***********************************************************************[assume.cpp]
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
#include "version.hpp"

using namespace ParaFROST;

inline void printOriginal(Lits_t& clause) {
    LOGN0("  adding original clause( ");
    forall_clause(clause, k) {
        PRINT("%d  ", SIGN(*k) ? -int(ABS(*k)) : int(ABS(*k)));
    }
    PRINT(")\n");
}

inline bool verifyMarkings(Vec<LIT_ST>& marks, Lits_t& clause) {
    forall_clause(clause, k) {
        const uint32 orglit = *k;
        assert(orglit > 1 && orglit < NOVAR);
        const uint32 orgvar = ABS(orglit);
        while (orgvar >= marks.size())
            marks.push(UNDEFINED);
        if (marks[orgvar] != UNDEFINED)
            return false;
    }
    return true;
}

uint32 Solver::iadd() {
    inf.unassigned++;
    const uint32 v = inf.orgVars = ++inf.maxVar;
    LOG2(3, "  adding new variable %d (%d unassigned)..", v, inf.unassigned);
    const uint32 lit = V2L(v);
    inf.maxDualVars = lit + 2;
    // Expanding non-SP structures.
    wt.expand(lit + 2);
    ifrozen.expand(v + 1, 0);
    bumps.expand(v + 1, 0);
    activity.expand(v + 1, 0.0);
    model.maxVar = v;
    model.lits.expand(v + 1);
    model.lits[v] = lit;
    vorg.expand(v + 1);
    vorg[v] = v;
    vmtf.init(v);
    vmtf.update(v, (bumps[v] = ++bumped));
    vsids.insert(v);
    if (sp == NULL) {
        sp = new SP();
        if (opts.proof_en)
            proof.init(sp);
    }
    if (!sp->selfallocated) {
        // Build mode: expanding SP structures the rebind.
        ivalue.expand(lit + 2, UNDEFINED);
        ilevel.expand(v + 1, UNDEFINED);
        iphase.expand(v + 1, opts.polarity);
        isource.expand(v + 1, NOREF);
        ivstate.expand(v + 1);
        ivstate[v] = VSTATE();
        sp->value = ivalue;
        sp->level = ilevel;
        sp->source = isource;
        sp->vstate = ivstate;
        sp->psaved = iphase;
    }
    return v;
}

bool Solver::itoClause(Lits_t& c, Lits_t& org) {
    if (org.empty()) {
        learnEmpty();
        LOG2(2, "  original clause is empty.");
        return false;
    }
    assert(c.empty());
    assert(verifyMarkings(imarks, org));
    bool satisfied = false;
    if (verbose >= 3) printOriginal(org);
    LOGN2(3, "  adding mapped clause  ( ");
    forall_clause(org, k) {
        const uint32 orglit = *k;
        assert(orglit > 1 && orglit < NOVAR);
        const uint32 orgvar = ABS(orglit);
        const LIT_ST sign = SIGN(orglit);
        imarks.expand(orgvar + 1, UNDEFINED);
        LIT_ST marker = imarks[orgvar];
        if (UNASSIGNED(marker)) {
            assert(sign >= 0);
            imarks[orgvar] = sign;
            uint32 mlit = imap(orgvar);
            CHECKLIT(mlit);
            uint32 mvar = ABS(mlit);
            mlit = V2DEC(mvar, sign);
            PRINT2(3, 5, "%d  ", SIGN(mlit) ? -int(mvar) : int(mvar));
            LIT_ST val = sp->value[mlit];
            if (UNASSIGNED(val))
                c.push(mlit);
            else if (val)
                satisfied = true; // satisfied unit
        } else if (NEQUAL(marker, sign))
            satisfied = true; // tautology
    }
    PRINT2(3, 5, ")\n");
    forall_clause(org, k) {
        const uint32 orglit = *k;
        assert(orglit > 1 && orglit < NOVAR);
        imarks[ABS(orglit)] = UNDEFINED;
    }
    if (satisfied) {
        if (opts.proof_en) proof.deleteClause(org);
    } else {
        int newsize = c.size();
        if (!newsize) {
            learnEmpty();
            LOG2(2, "  original clause became empty after parsing.");
            return false;
        } else if (newsize == 1) {
            const uint32 unit = *c;
            CHECKLIT(unit);
            LIT_ST val = sp->value[unit];
            if (UNASSIGNED(val))
                enqueueUnit(unit), formula.units++;
            else if (!val) {
                LOG2(2, "  unit clause(%d) is conflicting.", l2i(unit));
                return false;
            }
        } else if (newsize) {
            if (newsize == 2)
                formula.binaries++;
            else if (newsize == 3)
                formula.ternaries++;
            else
                assert(newsize > 3), formula.large++;
            if (newsize > formula.maxClauseSize)
                formula.maxClauseSize = newsize;
#ifdef LOGGING
            const C_REF newref =
#endif
                newClause(c, false);

            LOGCLAUSE(3, cm[newref], "  adding new clause");
        }
        if (opts.proof_en && newsize < org.size()) {
            proof.addClause(c);
            proof.deleteClause(org);
            org.clear();
        }
    }
    c.clear(), org.clear();
    return true;
}

bool Solver::ifailed(const uint32& v)
{
	const uint32 mlit = imap(v);
	if (!mlit) return false;
	CHECKLIT(mlit);
	const int size = iconflict.size();
	assert(size);
	for (int i = 0; i < size; i++) {
		if (ABS(mlit) == v)
			return true;
	}
	return false;
}

bool Solver::ieliminated(const uint32& v)
{
	const uint32 mlit = imap(v);
	if (!mlit) return true;
	CHECKLIT(mlit);
	return  MELTED(sp->vstate[ABS(mlit)].state) 
        ||  SUBSTITUTED(sp->vstate[ABS(mlit)].state)
        ||  FROZEN(sp->vstate[ABS(mlit)].state);
}

void Solver::ifreeze(const uint32& v) 
{
	const uint32 mlit = imap(v);
	CHECKLIT(mlit);
	const uint32 mvar = ABS(mlit);
	ifrozen[mvar] = 1;
	LOG2(3, "  freezing original variable %d (mapped to %d)..", v, mvar);
}

void Solver::iunfreeze(const uint32& v)
{
	const uint32 mlit = imap(v);
	CHECKLIT(mlit);
	const uint32 mvar = ABS(mlit);
	ifrozen[mvar] = 0;
	LOG2(3, "  melting original variable %d (mapped to %d)..", v, mvar);
}

void Solver::iassume(Lits_t& assumptions)
{
	assert(inf.maxVar);
	assert(stats.clauses.original);
	int assumed = assumptions.size();
	if (!assumed) return;
	stats.decisions.assumed += uint64(assumed);
	LOGN2(2, " Adding %d assumptions..", assumed);
	this->assumptions.reserve(assumed);
	forall_clause(assumptions, k) {
		const uint32 a = *k, v = ABS(a);
		assert(!ieliminated(v));
		uint32 mlit = imap(v);
		if (a != mlit && SIGN(a))
			mlit = FLIP(mlit);
		CHECKLIT(mlit);
		const uint32 mvar = ABS(mlit);
		ifrozen[mvar] = 1;
		this->assumptions.push(mlit);
		LOG2(4, "  assuming %d after mapping original assumption %d", l2i(mlit), l2i(a));
	}
	LOGDONE(2, 3);
}

void Solver::iunassume()
{
	assert(inf.maxVar);
	assert(stats.clauses.original);
	LOGN2(2, " Resetting %d assumptions and solver state..", assumptions.size());
	if (assumptions.size()) {
		forall_clause(assumptions, k) {
			const uint32 a = *k;
			CHECKLIT(a);
			ifrozen[ABS(a)] = 0;
		}
		assumptions.clear(true);
	}
	cnfstate = UNSOLVED;
	LOGDONE(2, 5);
	backtrack();
	iconflict.clear();
}

void Solver::idecide()
{
	assert(inf.unassigned);
	assert(isPropagated());
	assert(conflict == NOREF);
	assert(IS_UNSOLVED(cnfstate));
	int level = DL();
	uint32 dec = 0;
	while (level < assumptions.size()) {
		const uint32 a = assumptions[level];
		CHECKLIT(a);
		assert(ifrozen[ABS(a)]);
		const LIT_ST val = sp->value[a];
		if (UNASSIGNED(val)) {
			dec = a;
			break;
		}
		else if (val) incDL(), level = DL();
		else {
			assert(!val);
			ianalyze(FLIP(a));
			cnfstate = UNSAT;
			return;
		}
	}
	if (!dec) {
		const uint32 cand = vsidsEnabled() ? nextVSIDS() : nextVMFQ();
		dec = makeAssign(cand, useTarget());
	}
	enqueueDecision(dec);
	stats.decisions.single++;
}

void Solver::ianalyze(const uint32& failed)
{
	assert(vorg);
	LOG2(3, " Analyzing conflict on failed assumption (%d):", l2i(failed));
	LOGTRAIL(this, 3);
	iconflict.clear();
	iconflict.push(failed);
	if (!DL()) return;
	sp->seen[ABS(failed)] = ANALYZED_M;
	assert(trail.size());
	uint32* tail = trail.end() - 1, *pivot = trail + dlevels[1];
	while (tail >= pivot) {
		const uint32 parent = *tail--;
		CHECKLIT(parent);
		const uint32 parentv = ABS(parent);
		if (sp->seen[parentv]) {
			const C_REF r = sp->source[parentv];
			if (REASON(r)) {
				CLAUSE& c = cm[r];
				LOGCLAUSE(4, c, "  analyzing %d reason", l2i(parent));
				forall_clause(c, k) {
					const uint32 other = *k;
					if (other == parent) continue;
					const uint32 v = ABS(other);
					CHECKVAR(v);
					if (sp->level[v]) 
						sp->seen[v] = ANALYZED_M;
				}
			}
			else {
				assert(sp->level[parentv]);
				iconflict.push(FLIP(parent));
			}
			sp->seen[parentv] = 0;
		}
	}
	sp->seen[ABS(failed)] = 0;
}