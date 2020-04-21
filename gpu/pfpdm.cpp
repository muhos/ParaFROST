/***********************************************************************[pfpdm.cpp]
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

#include "pfsort.h"
#include "pfsolve.h" 

void ParaFROST::varOrder()
{
	scores.resize(nOrgVars());
	occurs.resize(nOrgVars());
	hist(orgs, true);
	hist(bins);
	bins.clear(true);
	for (uint32 v = 0; v < nOrgVars(); v++) {
		scores[v].v = v;
		scores[v].sc = (!occurs[v].ps || !occurs[v].ns) ? occurs[v].ps | occurs[v].ns : occurs[v].ps * occurs[v].ns;
	}
	Sort(scores, VAR_CMP());
	if (verbose == 4) {
		printf("c | Variable occurrences:\n");
		for (uint32 v = 0; v < nOrgVars(); v++)
			printf("c | Var[%d]->(v: %d, p: %d, n: %d, s: %d)\n", v, scores[v].v + 1, occurs[scores[v].v].ps, occurs[scores[v].v].ns, scores[v].sc);
	}
}

void ParaFROST::hist(BCNF& cnf, const bool& rst)
{
	if (cnf.size() == 0) return;
	if (rst) { for (int i = 0; i < occurs.size(); i++) occurs[i].ps = 0, occurs[i].ns = 0; }
	for (int i = 0; i < cnf.size(); i++) clHist(cnf[i]);
}

void ParaFROST::hist(LCNF& cnf, const bool& rst)
{
	if (cnf.size() == 0) return;
	if (rst) { for (int i = 0; i < occurs.size(); i++) occurs[i].ps = 0, occurs[i].ns = 0; }
	for (int i = 0; i < cnf.size(); i++) clHist(cnf[i]);
}

void ParaFROST::pumpFrozen()
{
	int nPDs = sp->trail_size - sp->trail_head;
	for (int i = sp->trail_size - 1; i >= sp->trail_head; i--) {
		uint32 v_idx = V2X(sp->trail[i]);
		WATCH* w_i = wt[sp->trail[i]], * w_end = w_i + wt[sp->trail[i]].size();
		double norm_act = (double)sol->level(v_idx) / nPDs;
		while (w_i != w_end) {
			C_REF c = (C_REF)w_i->c_ref;
			uint32 v0 = V2X(**c), v1 = V2X(*(*c + 1));
			uint32 frozen_v = (v0 == v_idx) ? v1 : v0;
			if (var_heap->varActivity(frozen_v) == 0.0) var_heap->varBumpAct(frozen_v, norm_act);
			w_i++;
		}
	}
}

void ParaFROST::PDMInit()
{
	if (!pdm_rounds) return;
	R = pdm_rounds;
	if (verbose >= 3) printf("c | Electing decisions (R=%d)..\n", R);
	stats.pdm_calls++;
	varOrder(); // initial variable ordering
	if (mcv_en) { // prefer the most-constrained variables
		for (int v = (int)nOrgVars() - 1; v >= 0; v--) {
			uint32 cand = scores[v].v;
			assert(cand < nOrgVars());
			if (sp->frozen[cand] || sp->lock[cand]) continue;
			uint32 p = V2D(cand + 1), n = NEG(p);
			uint32 dec = wt[p].size() >= wt[n].size() ? p : n;
			if (depFreeze_init(wt[dec], cand)) {
				incDL();
				enqueue(dec, DL());
				sp->seen[cand] = true;
				if (var_heap->has(cand)) var_heap->remove(cand);
			}
		}
	}
	else { // prefer the least constrained
		for (int v = 0; v < (int)nOrgVars(); v++) {
			uint32 cand = scores[v].v;
			assert(cand < nOrgVars());
			if (sp->frozen[cand] || sp->lock[cand]) continue;
			uint32 p = V2D(cand + 1), n = NEG(p);
			uint32 dec = wt[p].size() >= wt[n].size() ? p : n;
			if (depFreeze_init(wt[dec], cand)) {
				incDL();
				enqueue(dec, DL());
				sp->seen[cand] = true;
				if (var_heap->has(cand)) var_heap->remove(cand);
			}
		}
	}
	scores.clear(true);
	int nPDs = sp->trail_size - sp->trail_head;
	for (int i = 0; i < nPDs; i++) {
		uint32 v_idx = V2X(sp->trail[i + sp->trail_head]);
		sp->seen[v_idx] = 0;
		if (sp->frozen[v_idx]) {
			printf("c | Var(%d) is elected and frozen.\n", v_idx + 1);
			printWatched(v_idx);
		}
		assert(!sp->frozen[v_idx]);
	}
	bool* f = sp->frozen, * f_end = f + nOrgVars();
	while (f != f_end) *f++ = 0;
	ref_vars = nOrgVars() - nPDs, R--;
	stats.n_pds += nPDs;
	if (fdp_en) pumpFrozen(); // FUD Prioritization
}

void ParaFROST::PDM()
{
	if (verbose >= 3) printf("c | Electing decisions (R=%d)..\n", R);
	stats.pdm_calls++;
	eligibility();
	assert(scores.size() > 0);
	uint32 dec, cand;
	for (int i = 0; i < scores.size(); i++) {
		cand = scores[i].v;
		assert(cand < nOrgVars());
		if (sp->frozen[cand] || sp->lock[cand]) continue;
		dec = (polarity != 0) ? (sp->pol[cand] ? NEG(V2D(cand + 1)) : V2D(cand + 1)) : (drand() < 0.5 ? NEG(V2D(cand + 1)) : V2D(cand + 1));
		if (depFreeze(wt[dec], cand)) {
			incDL();
			enqueue(dec, DL());
			sp->seen[cand] = true;
			if (var_heap->has(cand)) var_heap->remove(cand);
		}
	}
	int nPDs = sp->trail_size - sp->trail_head;
	for (int i = 0; i < nPDs; i++) {
		uint32 v_idx = V2X(sp->trail[i + sp->trail_head]);
		sp->seen[v_idx] = 0;
		if (sp->frozen[v_idx]) {
			printf("c | Var(%d) is elected and frozen.\n", v_idx + 1);
			printWatched(v_idx);
		}
		assert(!sp->frozen[v_idx]);
	}
	bool* f = sp->frozen, * f_end = f + nOrgVars();
	while (f != f_end) *f++ = 0;
	ref_vars = nOrgVars() - nPDs, R--;
	stats.n_pds += nPDs;
}

bool ParaFROST::depFreeze_init(WL& ws, const uint32& cand) {
	WATCH* w_i = ws, * w_end = w_i + ws.size();
	while (w_i != w_end) { // check binary clauses
		assert(((B_REF)w_i->c_ref)->size() != 0);
		if (((B_REF)w_i->c_ref)->size() <= 2) return false;
		w_i++;
	}
	// freeze assignments
	w_i = ws; w_end = w_i + ws.size();
	while (w_i != w_end) {
		B_REF c = (B_REF)w_i->c_ref;
		assert(c->size());
		uint32 v0 = V2X(**c), v1 = V2X(*(*c + 1));
		uint32 frozen_v = (v0 == cand) ? v1 : v0;
		if (sp->seen[frozen_v]) return false; // v elected before, then cand cannot be elected 
		sp->frozen[frozen_v] = true;
		for (LIT_POS k = 2; k < c->size(); k++) {
			uint32 v = V2X((*c)[k]);
			if (sp->seen[v]) return false;
			sp->frozen[v] = true;
		}
		w_i++;
	}
	return true;
}

bool ParaFROST::depFreeze(WL& ws, const uint32& cand) {
	WATCH* w_end = ws + ws.size();
	for (WATCH* w_i = ws; w_i != w_end; w_i++) { // check binary clauses
		B_REF c = (B_REF)w_i->c_ref;
		//c->print();
		assert(c->isWatchedVar(cand));
		assert(c->size() > 1);
		assert(w_i->imp != 0);
		if (assigns[V2X(w_i->imp)] == !ISNEG(w_i->imp)) continue; // clause satisfied
		if (c->size() == 2) return false;
		bool satisfied = false, unAssigned = false;
		LIT_POS k = 2;
		while (k < c->size() && !satisfied && !unAssigned) {
			register LIT_ST h = assigns[V2X((*c)[k])];
			if (h == !ISNEG((*c)[k])) satisfied = true; // clause satisfied
			else if (h == UNDEFINED) unAssigned = true;
			k++;
		}
		if (!satisfied && !unAssigned) return false;
	}
	// freeze assignments
	for (WATCH* w_i = ws; w_i != w_end; w_i++) {
		B_REF c = (B_REF)w_i->c_ref;
		//c->print();
		if (assigns[V2X(w_i->imp)] == !ISNEG(w_i->imp)) continue; // clause satisfied
		uint32 v0 = V2X(**c), v1 = V2X(*(*c + 1));
		uint32 frozen_v = (v0 == cand) ? v1 : v0;
		if (sp->seen[frozen_v]) return false; // v elected before, then cand cannot be elected 
		sp->frozen[frozen_v] = true;
		for (LIT_POS k = 2; k < c->size(); k++) {
			assert((*c)[k] != 0);
			uint32 v = V2X((*c)[k]);
			if (sp->seen[v]) return false;
			sp->frozen[v] = true;
		}
	}
	return true;
}

void ParaFROST::eligibility()
{
	int nVars = var_heap->size();
	uint32* heap = var_heap->h_ptr();
	double* act = var_heap->act_ptr();
	scores.resize(nVars);
	if (pdm_order == 0) { // favor non-zero activity
		int h = 0, t = nVars - 1;
		for (int i = 0; i < nVars; i++) {
			if (act[heap[i]] > 1e-6) scores[h++].v = heap[i];
			else scores[t--].v = heap[i];
		}
		assert(h == t + 1);
	}
	else if (pdm_order == -1) { // favor zero-activity
		int h = 0, t = nVars - 1;
		for (int i = 0; i < nVars; i++) {
			if (act[heap[i]] <= 1e-6) scores[h++].v = heap[i];
			else scores[t--].v = heap[i];
		}
		assert(h == t + 1);
	}
	else if (pdm_order == 1) { // sort w.r.t activity and histogram
		hist(orgs, true);
		hist(learnts);
		for (int i = 0; i < nVars; i++) {
			uint32 v = heap[i];
			scores[i].v = v;
			scores[i].sc = (!occurs[v].ps || !occurs[v].ns) ? occurs[v].ps | occurs[v].ns : occurs[v].ps * occurs[v].ns;
		}
		Sort(scores, KEY_CMP_ACTIVITY<double>(act));
		/*for (int i = 0; i < nVars; i++)
			printf(" v[%d] = (h:%d, a:%f)\n", scores[i].v, scores[i].sc, act[scores[i].v]);
		printf("\n");*/
	}
}
