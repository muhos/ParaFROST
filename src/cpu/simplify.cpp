/***********************************************************************[simp.cpp]
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

#include "simplify.h"
#include "control.h"

using namespace pFROST;

inline bool	ParaFROST::checkMem(const string& _name, const size_t& size)
{
	const size_t sysMemCons = size_t(sysMemUsed()) + size;
	if (sysMemCons > size_t(stats.sysmem)) { // to catch memout problems before exception does
		PFLOGW("not enough memory for %s (free: %lld, used: %zd), sigma will terminate",
			_name.c_str(), stats.sysmem / MBYTE, sysMemCons / MBYTE);
		return false;
	}
	return true;
}

void ParaFROST::createOT(const bool& reset)
{
	if (opts.profile_simp) timer.pstart();
	if (reset) {
		forall_variables(v) { 
			const uint32 p = V2L(v); 
			ot[p].clear(); 
			ot[NEG(p)].clear();
		}
	}
	forall_vector(S_REF, scnf, i) {
		S_REF c = *i;
		if (c->learnt() || c->original()) {
			assert(c->size());
			forall_clause((*c), k) { 
				CHECKLIT(*k);
				ot[*k].push(*i);
			}
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.cot += timer.pcpuTime();
}

void ParaFROST::reduceOL(OL& ol)
{
	if (ol.empty()) return;
	S_REF *j = ol;
	forall_occurs(ol, i) {
		S_REF c = *i;
		if (c->deleted()) continue;
		*j++ = c;
	}
	ol.resize(int(j - ol));
}

void ParaFROST::reduceOT()
{
	if (opts.profile_simp) timer.pstart();
	forall_variables(v) {
		const uint32 p = V2L(v), n = NEG(p);
		reduceOL(ot[p]);
		reduceOL(ot[n]);
	}
	if (opts.profile_simp) timer.pstop(), timer.rot += timer.pcpuTime();
}

void ParaFROST::sortOT()
{
	if (opts.profile_simp) timer.pstart();
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32 v = PVs[i];
		CHECKVAR(v);
		const uint32 p = V2L(v), n = NEG(p);
		OL& poss = ot[p], &negs = ot[n];
		if (poss.size() > 1) Sort(poss, CNF_CMP_KEY());
		if (negs.size() > 1) Sort(negs, CNF_CMP_KEY());
	}
	if (opts.profile_simp) timer.pstop(), timer.sot += timer.pcpuTime();
}

void ParaFROST::extract(BCNF& cnf)
{
	assert(hc_isize == sizeof(uint32));
	assert(hc_scsize == sizeof(SCLAUSE));
	forall_cnf(cnf, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		const CLAUSE& c = cm[ref];
		const int size = c.size();
		const size_t bytes = hc_scsize + (size - 1) * hc_isize;
		S_REF s = (S_REF) new Byte[bytes];
		s->init(c);
		assert(s->size() == size);
		s->calcSig();
		rSort(s->data(), size);
		scnf[inf.nClauses++] = s;
		inf.nLiterals += size;
	}
}

void ParaFROST::sigmify()
{
	if (!opts.phases && !(opts.all_en || opts.ere_en)) return;
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));
	assert(stats.clauses.original);
	stats.sigma.calls++;
	sigmifying();
	INCREASE_LIMIT(this, sigma, stats.sigma.calls, nlognlogn, true);
	last.sigma.reduces = stats.reduces + 1;
	if (opts.phases > 2) {
		opts.phases--;
		PFLOG2(2, "  sigmify phases decreased to %d", opts.phases);
	}
}

void ParaFROST::awaken()
{ 
	initSimp();
	PFLOGN2(2, " Allocating memory..");
	const size_t numCls = size_t(maxClauses()), numLits = size_t(maxLiterals());
	const size_t ot_cap = inf.nDualVars * sizeof(OL) + numLits * sizeof(S_REF);
	const size_t scnf_cap = numCls * sizeof(S_REF) + numLits * sizeof(uint32);
	if (!checkMem("ot", ot_cap) || !checkMem("scnf", scnf_cap)) {
		simpstate = AWAKEN_FAIL; 
		return;
	}
	ot.resize(inf.nDualVars), scnf.resize(numCls);
	PFLENDING(2, 5, "(%.1f MB used)", double(ot_cap + scnf_cap) / MBYTE);
	PFLOGN2(2, " Extracting clauses to simplifying CNF..");
	printStats(1, '-', CGREEN0);
	inf.nClauses = inf.nLiterals = 0;
	wt.clear(true);
	extract(orgs), orgs.clear(true);
	extract(learnts), learnts.clear(true);
	scnf.resize(inf.nClauses);
	cm.destroy();
	PFLENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
	PFLMEMCALL(this, 2);
	return;
}

void ParaFROST::sigmifying()
{
	/********************************/
	/*         awaken sigma         */
	/********************************/
	SLEEPING(sleep.sigma, opts.sigma_sleep_en);
	rootify();
	shrinkTop(false);
	if (orgs.empty()) {
		recycleWT();
		return;
	}
	timer.stop();
	timer.solve += timer.cpuTime();
	if (!opts.profile_simp) timer.start();
	awaken();
	if (simpstate == AWAKEN_FAIL) {
		recycle();
		return;
	}
	if (interrupted()) killSolver();
	/********************************/
	/*      V/C Eliminations        */
	/********************************/
	assert(!phase && !mu_inc);
	int64 bmelted = inf.maxMelted, bclauses = inf.nClauses, bliterals = inf.nLiterals;
	int64 litsbefore = inf.nLiterals, diff = INT64_MAX;
	while (litsbefore) {
		resizeCNF();
		createOT();
		if (!prop()) killSolver();
		if (!LCVE()) break;
		sortOT();
		if (stop(diff)) { ERE(); break; }
		SUB(), VE(), BCE();
		countAll(), filterPVs();
		inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
		diff = litsbefore - inf.nLiterals, litsbefore = inf.nLiterals;
		phase++, mu_inc++;
		mu_inc += phase == opts.phases;
	}
	/********************************/
	/*          Write Back          */
	/********************************/
	// prop. remaining units if formula is empty
	// where recreating OT is not needed as there
	// are nothing to add
	if (!prop()) killSolver(); 
	assert(sp->propagated == trail.size());
	if (interrupted()) killSolver();
	occurs.clear(true), ot.clear(true);
	countFinal();
	shrinkSimp();
	assert(inf.nClauses == scnf.size());
	bool success = (bliterals != inf.nLiterals);
	stats.sigma.all.variables += int64(inf.maxMelted) - bmelted;
	stats.sigma.all.clauses += bclauses - int64(inf.nClauses);
	stats.sigma.all.literals += bliterals - int64(inf.nLiterals);
	last.shrink.removed = stats.shrunken;
	if (inf.maxFrozen > sp->simplified) stats.units.forced += inf.maxFrozen - sp->simplified;
	if (!inf.unassigned || !inf.nClauses) { 
		PFLOG2(2, " Formula is SATISFIABLE by elimination");
		cnfstate = SAT; 
		printStats(1, 's', CGREEN); 
		return;
	}
	if (canMap()) map(true); 
	else newBeginning();
	rebuildWT(opts.sigma_priorbins);
	if (retrail()) PFLOG2(2, " Propagation after sigmify proved a contradiction");
	UPDATE_SLEEPER(this, sigma, success);
	printStats(1, 's', CGREEN);
	if (!opts.profile_simp) timer.stop(), timer.simp += timer.cpuTime();
	if (!opts.solve_en) killSolver();
	timer.start();
}

void ParaFROST::shrinkSimp() 
{
	if (opts.profile_simp) timer.start();
	S_REF* j = scnf;
	forall_vector(S_REF, scnf, i) {
		S_REF c = *i;
		if (c->deleted()) deleteClause(c);
		else *j++ = c;
	}
	scnf.resize(uint32(j - scnf));
	if (opts.profile_simp) timer.stop(), timer.gc += timer.cpuTime();
}

void ParaFROST::newBeginning() 
{
	assert(opts.sigma_en || opts.sigma_live_en);
	assert(wt.empty());
	assert(orgs.empty());
	assert(learnts.empty());
	assert(inf.nClauses == scnf.size());
	const int64 litsCap = (inf.nLiterals - (inf.nClauses << 1)) * sizeof(uint32);
	assert(litsCap >= 0);
	const C_REF bytes = inf.nClauses * sizeof(CLAUSE) + size_t(litsCap);
	cm.init(bytes);
	stats.literals.original = stats.literals.learnt = 0;
	if (opts.aggr_cnf_sort) std::stable_sort(scnf.data(), scnf.data() + scnf.size(), CNF_CMP_KEY());
	forall_vector(S_REF, scnf, s) { newClause(**s); }
	stats.clauses.original = orgs.size();
	stats.clauses.learnt = learnts.size();
	assert(maxClauses() == int64(scnf.size()));
	scnf.clear(true);
}