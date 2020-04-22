

#include "pfsimpopts.h"
#include "pfsolve.h"
#include "pfsort.h"

void ParaFROST::free_master()
{
	
}

void ParaFROST::free_slaves()
{
	
}

void ParaFROST::optSimp()
{
	assert(pre_en);
	ngpus = opt_gpus;
	nstreams = opt_streams;
	p_cnf_en = opt_pcnf_en;
	p_ot_en = opt_pot_en;
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
	if (ngpus > devCount) ngpus = devCount;
}

void ParaFROST::_simplify()
{
	G_REF ref = BCP();
	assert(ref == NULL); // at this point BCP cannot return a conflict (dare to prove?!)
	if (sp->trail_size == sp->trail_offset) return;
	assert(sp->trail_size > 0 && DL() == ROOT_LEVEL);
	assert(sp->trail_head == sp->trail_size);
	assert(sol->assign(V2X(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	if (verbose >= 2) printf("c | Simplifying CNF..");
	// simplify watch table
	for (int i = sp->trail_offset; i < sp->trail_size; i++) {
		uint32 assign = sp->trail[i], assign_f = FLIP(assign);
		wt.collect(assign);
		wt.collect(assign_f);
		WL& ws = wt[assign_f];
		for (int j = 0; j < ws.size(); j++) {
			B_REF c = (B_REF)ws[j].c_ref;
			if (!c->garbage() && c->size() == 2) // original or learnt (undeleted) binary
				removeClause(c);
		}
		uint32 assign_idx = V2X(assign);
		if (var_heap->has(assign_idx)) var_heap->remove(assign_idx); // remove assign from the heap
		removed.push(assign); // save in case preprocessing mapped variables
	}
	// simplify input CNF
	cnf_stats.global_n_cls = simplify(orgs);
	orgs.resize(cnf_stats.global_n_cls);
	cnf_stats.global_n_del_vars += (sp->trail_size - sp->trail_offset);
	sp->trail_offset = sp->trail_size;
	// clean any garbage clauses
	wt.recycle();
	gcr.recycle();
	if (verbose >= 2) printf(" ==> done\n");
}

void ParaFROST::extractBins()
{
	assert(nClauses() == 0);
	assert(nLiterals() == 0);
	assert(cnf->size() == 0);
	assert(cnf->numLits() == 0);
	if (pre_delay) {
		if (nOrgBins()) {
			for (uint32 v = 0; v < nOrgVars(); v++) {
				uint32 p = V2D(v + 1), n = NEG(p);
				WL& poss = wt[n];
				for (int j = 0; j < poss.size(); j++) {
					assert(poss[j].c_ref != NULL);
					B_REF c = (B_REF)poss[j].c_ref;
					assert(c->size() > 1);
					if (!c->garbage() && c->size() == 2) {
						if (c->orgBin()) { assert(c->status() == UNKNOWN); cnf->attach(*c, 2); }
						removeClause(c);
					}
				}
				WL& negs = wt[p];
				for (int j = 0; j < negs.size(); j++) {
					assert(negs[j].c_ref != NULL);
					B_REF c = (B_REF)negs[j].c_ref;
					assert(c->size() > 1);
					if (!c->garbage() && c->size() == 2) {
						if (c->orgBin()) { assert(c->status() == UNKNOWN); cnf->attach(*c, 2); }
						removeClause(c);
					}
				}
			}
			gcr.recycle();
		}
	}
	else {
		for (int i = 0; i < bins.size(); i++) {
			assert(bins[i]->orgBin());
			assert(bins[i]->status() == UNKNOWN);
			cnf->attach(*bins[i], 2);
			collectClause(bins[i], false);
			assert(bins[i] == NULL);
		}
		bins.clear(true);
	}
	wt.clear(true);
}

void ParaFROST::GPU_CNF_fill() {
	// assign clause pointers and fill CNF with raw data
	cnf_stats.global_n_cls = 0;
	cnf_stats.global_n_lits = 0;
	extractBins();
	// append k-size clauses 
	for (int i = 0; i < orgs.size(); i++) {
		cnf->attach(*orgs[i], orgs[i]->size());
		collectClause(orgs[i], false);
		assert(orgs[i] == NULL);
	}
	cnf_stats.global_n_cls = cnf->size();
	cnf_stats.global_n_lits = cnf->numLits();
	// free cnfs
	orgs.clear(true);
	for (int i = 0; i < learnts.size(); i++) { collectClause(learnts[i], false); assert(learnts[i] == NULL); }
	learnts.clear();
	calc_sig(cnf, 0, nClauses());
	if (p_cnf_en) cnf->print();
}

bool ParaFROST::awaken()
{
	// simplify any remained facts at root level
	_simplify();
	// alloc fixed-size memory
	pv = new PV();
	if (!cuMem.allocStats(gstats, nOrgVars()) ||
		!cuMem.allocPV(pv, nOrgVars()) ||
		!cuMem.allocVO(d_occurs, d_scores, nOrgVars())) return false;
	// alloc memory for scnf
	if (!cuMem.resizeCNF(cnf, nClauses() + nBins(), nLiterals() + ((int64)nBins() << 1))) return false;
	// append orgs/bins to scnf
	GPU_CNF_fill();
	// alloc mmeory for variable stats
	occurs.resize(nOrgVars()), scores.resize(nOrgVars());
	histogram.resize(V2D(nOrgVars() + 1ULL)), raw_hist = thrust::raw_pointer_cast(histogram.data());
	// alloc memory for occur. table
	GPU_occurs(nLiterals(), true), ot = cuMem.resizeOT(raw_hist, nOrgVars(), nLiterals());
	if (ot == NULL) return false;
	// create streams
	createStreams();
	// reset vars states
	for (uint32 v = 0; v < nOrgVars(); v++) { sp->lock[v] = false; sp->seen[v] = false; sol->init(v); }
	sp->reset_trail();
	printPStats();
	return true;
}

void ParaFROST::GPU_preprocess()
{
	/********************************/
	/*         awaken sigma         */
	/********************************/
	timer->start();
	if (!awaken()) { timer->stop(); timer->pre = timer->cpuTime(); pre_en = false; return; } // not enough memory
	uint32 simpVars = nRemVars();
	/********************************/
	/*      1st-stage reduction     */
	/********************************/
	int64 lits_before = nLiterals(), lits_diff = INT64_MAX;
	int phase = 0;
	while (lits_diff > LIT_REM_THR && phase < phases) {
		if (interrupted()) killSolver();
		if (verbose > 1) printf("c |\t\tPhase-%d Variable Elections\n", phase);
		if (!LCVE()) goto writeBack;
		calc_added(cnf, ot, pv, gstats);
		if (!cuMem.resizeCNF(cnf, nClauses() + maxAddedCls(), nLiterals() + maxAddedLits(), phase)) goto writeBack;
		if (phase) create_ot(cnf, ot, p_ot_en);
		CNF_STATE veRet = GPU_VE(); 
		if (veRet == TERMINATE) goto writeBack;
		countCls(cnf, gstats); 
		cnf_stats.global_n_cls = cnf_stats.n_cls_after, cnf_stats.global_n_lits = cnf_stats.n_lits_after;
		lits_diff = lits_before - nLiterals(), lits_before = nLiterals();
		GPU_occurs(nLiterals()), ot = cuMem.resizeOT(raw_hist, nOrgVars(), nLiterals());
		if (ot == NULL) goto writeBack;
		phase++, pv->mu_inc++;
		pv->sol->assigns->print(true);
	}
	/********************************/
	/*      2nd-stage reduction     */
	/********************************/
	if (sub_en | hre_en | bce_en) {
		if (verbose > 1) printf("c | 2nd-Stage Clause Eliminations..\n");
		int t_p = mu_pos, t_n = mu_neg;
		while (t_p < CE_POS_LMT && t_n < CE_NEG_LMT) pv->mu_inc++, t_p <<= pv->mu_inc, t_n <<= pv->mu_inc;
		if (!LCVE()) goto writeBack;
		if (sub_en) GPU_SUB();
		if (bce_en) GPU_BCE();
		if (hre_en) GPU_HRE();
	}
	/********************************/
	/*           Write Back         */
	/********************************/
writeBack:
	if (interrupted()) killSolver();
	timer->stop();
	timer->pre = timer->cpuTime();
	printf("Pre time = %f\n", timer->cpuTime());
	destroyStreams();
	histogram.clear();
	exit(0);
	// call initial PDM again 
	PDMInit();
	pre_en = false;
	if (pre_delay) progRate += (initProgRate >> 2);
}

void ParaFROST::GPU_VO()
{
	assert(raw_hist != NULL);
	calc_vscores(d_occurs, d_scores, raw_hist);
	thrust::sort(d_scores, d_scores + nOrgVars(), VAR_CMP());
	if (verbose == 4) {
		printf("c | Variable occurrences:\n");
		for (uint32 v = 0; v < nOrgVars(); v++) {
			uint32 x = d_scores[v].v;
			printf("c | Var[%d]->(v: %d, p: %d, n: %d, s: %d)\n", v, x + 1, d_occurs[x].ps, d_occurs[x].ns, d_scores[v].sc);
		}
	}
}

void ParaFROST::GPU_occurs(const int64& numLits, const bool& init)
{
	assert(numLits > 0);
	rawLits.resize(numLits);
	if (init) copy(thrust::raw_pointer_cast(rawLits.data()), cnf, numLits);
	else {
		copyIf(thrust::raw_pointer_cast(rawLits.data()), cnf, gstats);
		assert(gstats->numLits == numLits);
	}
	GPU_hist();
	rawLits.clear();
}

void ParaFROST::GPU_hist()
{
	assert(rawLits.size() > 0);
	size_t n_bins = V2D(nOrgVars() + 1);
	assert(histogram.size() == n_bins);
	thrust::sort(rawLits.begin(), rawLits.end());
	thrust::counting_iterator<size_t> search_begin(0);
	thrust::upper_bound(rawLits.begin(), rawLits.end(), search_begin, search_begin + n_bins, histogram.begin());
	thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());
}

CNF_STATE ParaFROST::GPU_VE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating variables for round-0." << endl;
	ve(cnf, ot, pv);
	if (verbose > 1) { evalReds(cnf, gstats); printf("c |\t\tBVE Reductions\n"); logReductions(); }
	if (ve_plus_en) {
		countLits(cnf, gstats);
		int64 lits_before = cnf_stats.n_lits_after;
		int64 lits_removed = nLiterals() - cnf_stats.n_lits_after, lits_removed_pre = lits_removed;
		while (lits_removed > LIT_REM_THR) {
			// discard value variables
			uint32 n = 0;
			for (uint32 v = 0; v < pv->numPVs; v++) {
				uint32 x = pv->pVars->at(v);
				if (pv->sol->value[x] == UNDEFINED) (*pv->pVars)[n++] = x;
			}
			if (n == 0) { if (verbose > 1) printf("c | WARNING - Nothing left to eliminate\n"); break; }
			pv->numPVs = n, pv->pVars->resize(n);
			// Resizing OT (discard deleted, add resolvents)
			GPU_occurs(lits_before), ot = cuMem.resizeOT(raw_hist, nOrgVars(), lits_before);
			if (ot == NULL) return TERMINATE;
			create_ot(cnf, ot, p_ot_en);
			// HSE 
			GPU_SUB();
			if (verbose > 1) printf("c | Eliminating variables in new round..");
			ve(cnf, ot, pv);
			countLits(cnf, gstats); 
			lits_removed = (lits_before > cnf_stats.n_lits_after) ? lits_before - cnf_stats.n_lits_after : cnf_stats.n_lits_after - lits_before;
			if (verbose > 1) {
				if (lits_before > cnf_stats.n_lits_after) printf("(Literals reduction: -%lld) ==> done\n", lits_removed);
				else printf("(Literals reduction: +%lld) ==> done\n", lits_removed);
			}
			lits_before = cnf_stats.n_lits_after;
		}
		if (lits_removed_pre != lits_removed && verbose > 1) { evalReds(cnf, gstats); printf("c |\t\tBVE+ Reductions\n"); logReductions(); }
	}
	return UNSOLVED;
}

void ParaFROST::GPU_SUB()
{
	if (interrupted()) killSolver();
	if (verbose > 1) printf("c | HSE-ing parallel variables..");
	hse(cnf, ot, pv);
	if (verbose > 1) printf(" ==> done\n");
	if (verbose > 1) { evalReds(cnf, gstats); printf("c |\t\tHSE Reductions\n"); logReductions(); }
}

void ParaFROST::GPU_BCE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating blocked clauses..";
	bce(cnf, ot, pv);
	if (verbose > 1) printf(" ==> done\n");
	if (verbose > 1) { countCls(cnf, gstats); printf("c |\t\tBCE Reductions\n"); logReductions(); }
}

void ParaFROST::GPU_HRE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating hidden redundances..";
	hre(cnf, ot, pv);
	if (verbose > 1) printf(" ==> done\n");
	if (verbose > 1) { countCls(cnf, gstats); printf("c |\t\tHRE Reductions\n"); logReductions(); }
}

void ParaFROST::depFreeze(const OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	for (uint32 i = 0; i < ol.size(); i++) {
		SCLAUSE& c = (*cnf)[ol[i]];
		for (LIT_POS k = 0; k < c.size(); k++) {
			register uint32 v = V2X(c[k]);
			if (v != cand && (occurs[v].ps < p_temp || occurs[v].ns < n_temp)) sp->frozen[v] = true;
		}
	}
}

bool ParaFROST::LCVE()
{
	GPU_VO();
	cudaMemcpyAsync(occurs.data(), d_occurs, nOrgVars() * sizeof(OCCUR), cudaMemcpyDeviceToHost, streams[0]);
	cudaMemcpyAsync(scores.data(), d_scores, nOrgVars() * sizeof(SCORE), cudaMemcpyDeviceToHost, streams[1]);
	create_ot(cnf, ot, p_ot_en);
	pv->numPVs = 0, pv->pVars->clear();
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 cand = scores[v].v;
		if (sp->frozen[cand]) continue;
		uint32 p = V2D(cand + 1), n = NEG(p);
		uint32 poss_sz = (*ot)[p].size(), negs_sz = (*ot)[n].size();
		assert(poss_sz == occurs[cand].ps);
		assert(negs_sz == occurs[cand].ns);
		if (poss_sz == 0 && negs_sz == 0) continue;
		uint32 pos_temp = mu_pos << pv->mu_inc, neg_temp = mu_neg << pv->mu_inc;
		if (poss_sz >= pos_temp && negs_sz >= neg_temp) break;
		assert(pv->sol->value[cand] == UNDEFINED);
		pv->pVars->_push(cand);
		depFreeze((*ot)[p], cand, pos_temp, neg_temp);
		depFreeze((*ot)[n], cand, pos_temp, neg_temp);
	}
	pv->numPVs = pv->pVars->size();
	assert(verifyLCVE());
	if (verbose >= 3) {
		printf("c | PLCVs (n = %d): \n", pv->numPVs);
		printVars(*pv->pVars, pv->numPVs);
	}
	bool* f = sp->frozen, * f_end = f + nOrgVars();
	while (f != f_end) *f++ = 0;
	if (pv->numPVs < MIN_PVARS) {
		if (verbose > 1) printf("c | WARNING - parallel variables not enough, simplifications will be terminated\n");
		return false;
	}
	return true;
}