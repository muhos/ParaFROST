/***********************************************************************[pfsimp.cpp]
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
#include "pfsimpopts.h"
using namespace pFROST;
using namespace SIGmA;

void ParaFROST::optSimp()
{
	assert(sigma_en);
	ve_en = opt_ve_en || opt_ve_plus_en;
	ve_plus_en = opt_ve_plus_en;
	sub_en = opt_sub_en;
	bce_en = opt_bce_en;
	hre_en = opt_hre_en;
	all_en = opt_all_en;
	phases = opt_phases;
	mu_pos = opt_mu_pos;
	mu_neg = opt_mu_neg;
	cnf_free_freq = opt_cnf_free;
	if (all_en) ve_en = 1, ve_plus_en = 1, sub_en = 1, bce_en = 1, hre_en = 1;
	if (!phases && ve_en) phases = 1; // at least 1 phase needed for BVE(+)
	if (phases && !ve_en) phases = 0;
}

void ParaFROST::varReorder()
{
	PFLOGN2(2, " Finding eligible variables for LCVE..");
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	assert(!scnf.empty());
	histSimp(scnf, true);
	uint32* scores = sp->tmp_stack;
	for (uint32 v = 1; v <= inf.maxVar; v++) eligible[v - 1] = v, scores[v] = rscore(v);
	Sort(eligible, LCV_CMP(scores));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible variables:");
		for (uint32 i = 0; i < eligible.size(); i++) {
			uint32 v = eligible[i];
			PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", i, v, occurs[v].ps, occurs[v].ns, scores[v]);
		}
	}
}

void ParaFROST::newClause(S_REF c, const bool& added)
{
	assert(*c != NULL);
	assert(c->size());
	assert(c->hasZero() < 0);
	if (added) {
		assert(c->status() == 0);
		assert(c->isSorted());
		if (c->size() == 1 && unassigned(**c)) enqueue(**c);
		else {
			c->set_status(ORIGINAL);
			c->calcSig();
			scnf.push(c);
		}
		for (int l = 0; l < c->size(); l++) ot[(*c)[l]].push(c); // attach to OT
		if (proof_en) {
			wrProof('a');
			wrProof(*c, c->size());
			wrProof(0);
		}
	}
	else {
		assert(c->size() > 1);
		c->calcSig();
		Sort(c->data(), c->size());
		scnf[inf.nClauses++] = c;
		inf.nLiterals += c->size();
	}
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
	if (proof_en) {
		wrProof('a');
		wrProof(*c, c->size());
		wrProof(0);
	}
	if (c->size() == 1 && unassigned(**c)) enqueue(**c);
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

bool ParaFROST::prop()
{
	while (sp->propagated < trail.size()) { // propagate units
		uint32 assign = trail[sp->propagated++], f_assign = flip(assign);
		assert(assign > 0);
		PFLBCP(this, 4, assign);
		// remove satisfied
		for (int i = 0; i < ot[assign].size(); i++) ot[assign][i]->markDeleted();
		// reduce unsatisfied
		for (int i = 0; i < ot[f_assign].size(); i++) {
			S_REF c = ot[f_assign][i];
			//c->print();
			assert(c->size());
			if (c->status() == DELETED || propClause(c, f_assign)) continue; // clause satisfied
			// clause is unit or conflict
			// Note: attach & strengthen don't check for conflict before enqueue
			if (c->size() == 0) { cnfstate = UNSAT; return false; }
			if (c->size() == 1) { 
				assert(**c > 1);
				if (unassigned(**c)) enqueue(**c);
				else { cnfstate = UNSAT; return false; } 
			}
		}
		// delete assign lists
		ot[assign].clear(true), ot[f_assign].clear(true);
	}
	return true;
}

void ParaFROST::create_ot(const bool& rst)
{
	// reset ot
	if (rst) {
		for (uint32 v = 1; v <= inf.maxVar; v++) { 
			uint32 p = v2l(v); 
			ot[p].clear(); 
			ot[neg(p)].clear();
		}
	}
	// create ot
	for (uint32 i = 0; i < scnf.size(); i++) {
		SCLAUSE& c = *scnf[i];
		if (c.status() != DELETED) {
			assert(c.status() > 0);
			assert(c.size() > 1);
			for (int k = 0; k < c.size(); k++) { 
				assert(c[k] != 0);
				ot[c[k]].push(scnf[i]);
			}
		}
	}
	assert(consistent(scnf, ot));
}

void ParaFROST::reduce_ol(OL& ol)
{
	if (ol.empty()) return;
	int new_sz = 0;
	for (int i = 0; i < ol.size(); i++)
		if (ol[i]->status() != DELETED)
			ol[new_sz++] = ol[i];
	ol.resize(new_sz);
}

void ParaFROST::reduce_ot()
{
	for (uint32 v = 1; v <= inf.maxVar; v++) {
		uint32 p = v2l(v), n = neg(p);
		reduce_ol(ot[p]);
		reduce_ol(ot[n]);
	}
	assert(consistent(scnf, ot));
}

void ParaFROST::extract(BCNF& cnf)
{
	for (uint32 i = 0; i < cnf.size(); i++) {
		CLAUSE& c = cm[cnf[i]];
		assert(!c.deleted());
		S_REF sc = new SCLAUSE(c);
		sc->set_status(c.status());
		newClause(sc, false);
	}
	cnf.clear(true);
}

bool ParaFROST::awaken()
{ 
	// deal with any remained facts at root level
	PFLOG2(2, " Propagating any remaining facts before eliminations..");
	C_REF cref = BCP();
	assert(cref == NOREF); // dare to prove?!
	PFLOG2(2, " All good.");
	assert(DL() == ROOT_LEVEL);
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	assert(sp->propagated == trail.size());
	if (sigma_live_en) {
		PFLOG2(2, " Shrinking CNF before eliminations..");
		shrinkWT();
		shrink(orgs);
		shrink(learnts);
		inf.nDelVars += (trail.size() - sp->simplified);
		PFLENDING(2, 5, " (-%d variables)", trail.size() - sp->simplified);
		sp->simplified = trail.size();
	}
	// alloc simplifier memory 
	PFLOGN2(2, " Allocating memory..");
	size_t ot_cap = inf.nDualVars * sizeof(OL);
	size_t scnf_cap = inf.nOrgCls * sizeof(S_REF);
	if (!checkMem("ot", ot_cap) || !checkMem("scnf", scnf_cap)) return false;
	ot.resize(inf.nDualVars);
	scnf.resize(inf.nOrgCls);
	PFLENDING(2, 5, "(%.1f MB used)", double(ot_cap + scnf_cap) / MBYTE);
	// append clauses to scnf
	PFLOGN2(2, " Extracting clauses..");
	BCNF bins;
	bins.reserve(inf.nOrgBins + inf.nLearntBins);
	wtBin.recycle();
	extractBins(bins);
	inf.nClauses = inf.nLiterals = 0;
	inf.nOrgBins = inf.nLearntBins = 0;
	extract(bins);
	extract(orgs);
	extract(learnts);
	scnf.resize(inf.nClauses);
	PFLENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
	// cleaning up
	wtBin.clear(true), wt.clear(true), cm.destroy(); 
	printPStats();
	PFLMEMCALL(this, 2);
	return true;
}

void ParaFROST::shrinkSimp()
{
	inf.nClauses = 0;
	inf.nLiterals = 0;
	for (uint32 i = 0; i < scnf.size(); i++) {
		S_REF c = scnf[i];
		if (c->status() == DELETED) removeClause(c);
		else {
			assert(c->size() > 1);
			scnf[inf.nClauses++] = c;
			inf.nLiterals += c->size();
		}
	}
	scnf.resize(inf.nClauses);
}

void ParaFROST::preprocess()
{
	if (!phases && !sub_en && !bce_en && !hre_en) return;
	/********************************/
	/*         awaken sigma         */
	/********************************/
	timer.start();
	assert(cnfstate == UNSOLVED);
	assert(conflict == NOREF);
	int64 lits_before, lits_diff = INT64_MAX;
	int phase = 0;
	if (!awaken()) goto writeBack;
	if (interrupted()) killSolver();
	/********************************/
	/*      1st-stage reduction     */
	/********************************/
	mu_inc = 0, lits_before = inf.nLiterals;
	create_ot(false); // created only initially
	while (lits_diff > LIT_REM_THR && phase < phases) {
		PFLOG2(2, "\t\tPhase-%d Variable Elections", phase);
		LCVE(); 
		VE();
		countCls(); 
		inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
		lits_diff = lits_before - inf.nLiterals, lits_before = inf.nLiterals;
		if (phase + 1 != phases && ((phase + 1) % cnf_free_freq == 0)) { shrinkSimp(); create_ot(); }
		else reduce_ot();
		phase++, mu_inc++;
	} 
	/********************************/
	/*      2nd-stage reduction     */
	/********************************/
	if (sub_en | hre_en | bce_en) {
		PFLOGN2(2, " Initiating clause eliminations..");
		int t_p = mu_pos, t_n = mu_neg;
		while (t_p < CE_POS_LMT && t_n < CE_NEG_LMT) mu_inc++, t_p <<= mu_inc, t_n <<= mu_inc;
		PFLDONE(2, 5);
		LCVE();
		if (sub_en) SUB();
		if (bce_en) BCE();
		if (hre_en) HRE();
	}
	/********************************/
	/*           Write Back         */
	/********************************/
writeBack:
	if (interrupted()) killSolver();
	occurs.clear(true), ot.clear(true);
	shrinkSimp(), countMelted();
	printPStats();
	if (canMap()) map();
	else assert(!mapped), newBeginning();
	timer.stop();
	timer.pre = timer.cpuTime();
}

void ParaFROST::newBeginning() {
	assert(sigma_en || sigma_live_en);
	assert(!scnf.empty());
	assert(wtBin.empty()), assert(wt.empty());
	assert(orgs.empty()), assert(learnts.empty());
	inf.nOrgBins = inf.nLearntBins = 0;
	inf.nLiterals = inf.nLearntLits = 0;
	assert(inf.maxVar > vmap.numVars());
	uint32 tableSize = mapped ? v2l(vmap.size()) : inf.nDualVars;
	wtBin.resize(tableSize), wt.resize(tableSize);
	cm.init(scnf.size());
	if (!mapped) assert(vmap.empty()), sp->lockMelted();
	createWT(), copyWatched(), copyNonWatched(); // must follow this order
	assert(orgs.size() + learnts.size() + inf.nOrgBins + inf.nLearntBins == scnf.size());
	inf.nOrgCls = inf.nClauses = orgs.size();
	inf.nOrgLits = inf.nLiterals;
	scnf.clear(true);
}

C_REF ParaFROST::newClause(SCLAUSE& in_c)
{
	C_REF r = in_c.ref();
	if (r == NOREF) {
		int sz = in_c.size();
		assert(sz > 1);
		r = cm.alloc(sz);
		in_c.set_ref(r); // new ref overwrites simplified clause sig.
		CLAUSE& new_c = cm[r];
		if (mapped) vmap.mapClause(new_c, in_c);
		else new_c.copyLitsFrom(in_c);
		assert(sz == new_c.size());
		assert(new_c[0] > 0 && new_c[1] > 0);
		assert(new_c[0] <= UINT32_MAX && new_c[1] <= UINT32_MAX);
		new_c.set_status(in_c.status());
		if (sz == 2) {
			if (in_c.status() == ORIGINAL) inf.nOrgBins++;
			else assert(in_c.status() == LEARNT), inf.nLearntBins++;
		}
		else {
			assert(sz > 2);
			if (in_c.status() == LEARNT) {
				new_c.set_LBD(new_c.size());
				learnts.push(r);
				inf.nLearntLits += sz;
			}
			orgs.push(r);
			inf.nLiterals++;
		}
	}
	return r;
}

void ParaFROST::depFreeze(const OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	for (int i = 0; i < ol.size(); i++) {
		SCLAUSE& c = *ol[i];
		for (int k = 0; k < c.size(); k++) {
			register uint32 v = l2a(c[k]);
			if (v != cand && (occurs[v].ps < p_temp || occurs[v].ns < n_temp)) sp->frozen[v] = 1;
		}
	}
}

void ParaFROST::LCVE()
{
	// reorder variables
	varReorder();
	// extended LCVE
	PFLOGN2(2, " Electing variables..");
	PVs.clear();
	for (uint32 i = 0; i < eligible.size(); i++) {
		uint32 cand = eligible[i];
		if (sp->vstate[cand] == FROZEN || sp->vstate[cand] == MELTED) continue;
		if (sp->frozen[cand]) continue;
		uint32 p = v2l(cand), n = neg(p);
		uint32 poss_sz = (uint32)ot[p].size(), negs_sz = (uint32)ot[n].size();
		assert(poss_sz == occurs[cand].ps);
		assert(negs_sz == occurs[cand].ns);
		assert(poss_sz || negs_sz);
		uint32 pos_temp = mu_pos << mu_inc, neg_temp = mu_neg << mu_inc;
		if (poss_sz >= pos_temp && negs_sz >= neg_temp) break;
		assert(sp->vstate[cand] == ACTIVE);
		PVs.push(cand);
		depFreeze(ot[p], cand, pos_temp, neg_temp);
		depFreeze(ot[n], cand, pos_temp, neg_temp);
	}
	assert(verifyLCVE());
	PFLENDING(2, 4, "(%d elected)", PVs.size());
	if (verbose >= 3) { PFLOGN0(" PLCVs "); printVars(PVs, PVs.size(), 'v'); }
	memset(sp->frozen, 0, inf.maxVar + 1ULL);
}

void ParaFROST::bve()
{
	if (interrupted()) killSolver();
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32 v = PVs[i];
		assert(v);
		assert(sp->vstate[v] == ACTIVE);
		uint32 p = v2l(v), n = neg(p);
		if (ot[p].size() == 0 || ot[n].size() == 0) {
			toblivion(model.resolved, ot[p], ot[n], v);
			sp->vstate[v] = MELTED;
			
		}
		else if (ot[p].size() == 1 || ot[n].size() == 1) { 
			resolve_x(v, ot[p], ot[n]); 
			toblivion(model.resolved, ot[p], ot[n], v);
			sp->vstate[v] = MELTED;
		}
		else if (gateReasoning_x(p, ot[p], ot[n]) || mayResolve_x(v, ot[p], ot[n])) {
			toblivion(model.resolved, ot[p], ot[n], v);
			sp->vstate[v] = MELTED;
		}
	}
	if (!prop()) killSolver();
}

void ParaFROST::VE()
{
	if (PVs.size() <= MIN_PARALLEL_VARS) return;
	int64 lits_before = inf.nLiterals, lits_removed = INT64_MAX;
	int round = 0;
	while (lits_removed > LIT_REM_THR && filterPVs()) {
		PFLOG2(2, " Elimination round %d:", round++);
		if (ve_plus_en) {
			// HSE 
			PFLOGN2(2, "  1) HSE-ing variables..");
			HSE();
			countLits(); // count remained literals
			lits_removed = lits_before - inf.n_lits_after;
			lits_before = inf.n_lits_after;
			assert(lits_removed >= 0);
			PFLENDING(2, 5, "(Literals removed : %lld)", -lits_removed);
		}
		// BVE 
		PFLOGN2(2, "  2) Eliminating variables..");
		bve();
		countLits(); // count remained literals
		lits_removed = lits_before - inf.n_lits_after;
		PFLENDING(2, 5, "(Literals removed : %c%lld)", lits_removed < 0 ? '+' : ' ', -lits_removed);
		lits_before = inf.n_lits_after;
	}
	PFLREDALL(this, 2, "BVE(+) Reductions");
}

void ParaFROST::HSE()
{
	if (interrupted()) killSolver();
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32 v = PVs[i];
		assert(v);
		assert(sp->vstate[v] == ACTIVE);
		uint32 p = v2l(v), n = neg(p);
		self_sub_x(p, ot[p], ot[n]);
	}
	if (!prop()) killSolver();
}

void ParaFROST::SUB()
{
	if (interrupted()) killSolver();
	PFLOGN2(2, " SUB-ing variables..");
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32 v = PVs[i];
		assert(v);
		assert(sp->vstate[v] == ACTIVE);
		uint32 p = v2l(v), n = neg(p);
		self_sub_x(p, ot[p], ot[n]);
	}
	if (!prop()) killSolver();
	PFLDONE(2, 5);
	PFLREDCL(this, 2, "SUB Reductions");
}

void ParaFROST::BCE()
{
	if (interrupted()) killSolver();
	PFLOGN2(2, " Eliminating blocked clauses..");
	for (uint32 i = 0; i < PVs.size(); i++) {
		uint32 v = PVs[i];
		assert(v);
		assert(sp->vstate[v] == ACTIVE);
		uint32 p = v2l(v), n = neg(p);
		blocked_x(v, ot[p], ot[n]);
	}
	PFLDONE(2, 5);
	PFLREDALL(this, 2, "BCE Reductions");
}

void ParaFROST::HRE()
{
	if (interrupted()) killSolver();
	PFLOGN2(2, " Eliminating hidden redundances..");
	Lits_t m_c;
	m_c.reserve(INIT_CAP);
	for (uint32 n = 0; n < PVs.size(); n++) {
		assert(PVs[n]);
		assert(sp->vstate[PVs[n]] == ACTIVE);
		uint32 posY = v2l(PVs[n]);
		OL& y_poss = ot[posY], &y_negs = ot[neg(posY)];
		// do merging and apply forward equality check (on-the-fly) over resolvents
		for (int i = 0; i < y_poss.size(); i++) {
			if (y_poss[i]->status() == DELETED) continue;
			for (int j = 0; j < y_negs.size(); j++) {
				if (y_negs[j]->status() == DELETED) continue;
				if (merge_hre(PVs[n], y_poss[i], y_negs[j], m_c)) {
					uint32 m_sig = 0;
					for (int k = 0; k < m_c.size(); k++) m_sig |= MAPHASH(m_c[k]);
					for (int k = 0; k < m_c.size(); k++)
						if (forward_equ(m_c, m_sig, ot[m_c[k]])) break;
				} // end-if (merge)
			} // end-for (negs_list)
		} // end-for (poss_list)
	}
	PFLDONE(2, 5);
	PFLREDCL(this, 2, "HRE Reductions");
}

bool ParaFROST::consistent(const SCNF& cnf, const OT& ot)
{
	for (uint32 i = 0; i < cnf.size(); i++) {
		S_REF c = cnf[i];
		assert(c->size());
		//c->print();
		for (int k = 0; k < c->size(); k++) {
			assert(c->lit(k));
			const OL& ol = ot[c->lit(k)];
			if (c->status() != DELETED) assert(ol.size());
			bool found = false;
			for (int j = 0; j < ol.size(); j++) {
				S_REF ref = ol[j];
				assert(ref->size());
				if (ref == c) {
					found = true;
					break;
				}
			}
			if (c->status() == DELETED && found) {
				printf(" c");
				c->print();
				printf(" found DELETED in list[%d]:", l2i(c->lit(k)));
				printOL(ol);
				return false;
			}
			if (c->status() != DELETED && !found) {
				printf(" c");
				c->print();
				printf(" NOT found in list[%d]:", l2i(c->lit(k)));
				printOL(ol);
				return false;
			}
		}
	}
	return true;
}