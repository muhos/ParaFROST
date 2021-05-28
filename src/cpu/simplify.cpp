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

using namespace pFROST;
using namespace SIGmA;

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
	for (uint32 i = 0; i < scnf.size(); i++) {
		const SCLAUSE& c = *scnf[i];
		if (c.learnt() || c.original()) {
			assert(c.size());
			for (int k = 0; k < c.size(); k++) { 
				CHECKLIT(c[k]);
				ot[c[k]].push(scnf[i]);
			}
		}
	}
	if (opts.profile_simp) timer.pstop(), timer.cot += timer.pcpuTime();
}

void ParaFROST::reduceOL(OL& ol)
{
	if (ol.empty()) return;
	int n = 0;
	for (int i = 0; i < ol.size(); i++) {
		if (ol[i]->deleted()) continue;
		ol[n++] = ol[i];
	}
	ol.resize(n);
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
		std::sort(poss.data(), poss.data() + poss.size(), CNF_CMP_KEY());
		std::sort(negs.data(), negs.data() + negs.size(), CNF_CMP_KEY());
	}
	if (opts.profile_simp) timer.pstop(), timer.sot += timer.pcpuTime();
}

void ParaFROST::newSClause(S_REF s)
{
	assert(*s != NULL);
	assert(s->size() > 1);
	s->calcSig();
	rSort(s->data(), s->size());
	assert(s->isSorted());
	scnf[inf.nClauses++] = s;
	inf.nLiterals += s->size();
}

void ParaFROST::extract(const BCNF& cnf)
{
	for (uint32 i = 0; i < cnf.size(); i++) {
		const CLAUSE& c = cm[cnf[i]];
		if (c.deleted()) continue;
		newSClause(new SCLAUSE(c));
	}
}

void ParaFROST::sigmify()
{
	if (!opts.phases && !(opts.all_en || opts.ere_en)) return;
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
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
	size_t numCls = maxClauses(), numLits = maxLiterals();
	size_t ot_cap = inf.nDualVars * sizeof(OL) + numLits * sizeof(S_REF);
	size_t scnf_cap = numCls * sizeof(S_REF) + numLits * sizeof(uint32);
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
	assert(cnfstate == UNSOLVED);
	PFLOGN2(2, " Shrinking CNF before sigmification..");
	if (sp->simplified < inf.maxFrozen) sp->simplified = inf.maxFrozen;
	int64 beforeCls = maxClauses(), beforeLits = maxLiterals();
	shrinkTop(orgs), shrinkTop(learnts);
	PFLSHRINKALL(this, 2, beforeCls, beforeLits);
	assert(stats.clauses.original == orgs.size());
	assert(stats.clauses.learnt == learnts.size());
	timer.stop(), timer.solve += timer.cpuTime();
	if (!opts.profile_simp) timer.start();
	awaken();
	if (simpstate == AWAKEN_FAIL) return;
	if (interrupted()) killSolver();
	/********************************/
	/*      V/C Eliminations        */
	/********************************/
	assert(!phase && !mu_inc);
	int64 bmelted = inf.maxMelted, bclauses = inf.nClauses, bliterals = inf.nLiterals;
	int64 litsbefore = inf.nLiterals, diff = INT64_MAX;
	while (true) {
		resizeCNF();
		createOT();
		if (!prop()) killSolver();
		if (!LCVE()) break;
		sortOT();
		if (stop(diff)) { ERE(); break; }
		HSE(), VE(), BCE();
		countAll(), filterPVs();
		inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
		diff = litsbefore - inf.nLiterals, litsbefore = inf.nLiterals;
		phase++, mu_inc++;
		mu_inc += phase == opts.phases;
	}
	/********************************/
	/*          Write Back          */
	/********************************/
	assert(sp->propagated == trail.size());
	if (interrupted()) killSolver();
	occurs.clear(true), ot.clear(true);
	countFinal();
	shrinkSimp();
	assert(inf.nClauses == scnf.size());
	if (!inf.unassigned || !inf.nClauses) {
		cnfstate = SAT;
		printStats(1, 's', CGREEN);
		return;
	}
	bool success = (bliterals != inf.nLiterals);
	stats.sigma.all.variables += int64(inf.maxMelted) - bmelted;
	stats.sigma.all.clauses += bclauses - int64(inf.nClauses);
	stats.sigma.all.literals += bliterals - int64(inf.nLiterals);
	last.shrink.removed = stats.shrunken;
	if (inf.maxFrozen > sp->simplified) stats.units.forced += inf.maxFrozen - sp->simplified;
	if (canMap()) map(true); 
	else newBeginning();
	rebuildWT(opts.sigma_priorbins);
	if (BCP()) {
		PFLOG2(1, " Propagation after sigmify proved a contradiction");
		learnEmpty();
	}
	UPDATE_SLEEPER(this, sigma, success);
	printStats(1, 's', CGREEN);
	if (!opts.profile_simp) timer.stop(), timer.simp += timer.cpuTime();
	if (!opts.solve_en) killSolver();
	timer.start();
}

void ParaFROST::shrinkSimp() {
	if (opts.profile_simp) timer.start();
	uint32 n = 0;
	for (uint32 i = 0; i < scnf.size(); i++) {
		S_REF c = scnf[i];
		if (c->deleted()) deleteClause(c);
		else scnf[n++] = c;
	}
	scnf.resize(n);
	if (opts.profile_simp) timer.stop(), timer.gc += timer.cpuTime();
}

void ParaFROST::newBeginning() {
	assert(opts.sigma_en || opts.sigma_live_en);
	assert(wt.empty());
	assert(orgs.empty());
	assert(learnts.empty());
	cm.init(scnf.size());
	stats.literals.original = stats.literals.learnt = 0;
	if (opts.aggr_cnf_sort) std::stable_sort(scnf.data(), scnf.data() + scnf.size(), CNF_CMP_KEY());
	for (S_REF* s = scnf; s != scnf.end(); s++) newClause(**s);
	stats.clauses.original = orgs.size();
	stats.clauses.learnt = learnts.size();
	assert(maxClauses() == scnf.size());
	scnf.clear(true);
}