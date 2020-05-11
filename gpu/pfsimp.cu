

#include "pfsimpopts.h"
#include "pfsolve.h"
#include "pfsort.h"

void ParaFROST::masterFree()
{
	if (h_hist != NULL) delete[] h_hist;
	if (h_cnf != NULL) h_cnf->~CNF(), h_cnf = NULL;
	if (cnf != NULL) cnf->~CNF(), cnf = NULL;
	if (ot != NULL) ot->~OT(), ot = NULL;
	if (pv != NULL) pv->~PV(), pv = NULL;
	d_hist = NULL;
	histogram.clear(), rawLits.clear();
	histogram.shrink_to_fit(), rawLits.shrink_to_fit();
	destroyStreams();
	cuMem.~cuMM();
}

void ParaFROST::slavesFree()
{
	
}

void ParaFROST::optSimp()
{
	assert(pre_en);
	clsTile = opt_tile_size;
	litsCopyR = opt_lits_cop_r;
	ngpus = opt_gpus;
	nstreams = opt_streams;
	solve_en = opt_solve_en;
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
	cls_en = all_en || sub_en || bce_en || hre_en;
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
	assert(sol->value(V2X(sp->trail[sp->trail_size - 1])) != UNDEFINED);
	if (verbose > 1) printf("c | Simplifying CNF..");
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
	if (verbose > 1) printf(" ==> done\n");
}

void ParaFROST::extractBins()
{
	assert(cnf->size() == 0);
	assert(cnf->numLits() == 0);
	if (pre_delay) {
		if (nOrgBins()) {
			bins.clear(true), rawBins.reserve((nOrgBins() << 1));
			for (uint32 v = 0; v < nOrgVars(); v++) {
				uint32 p = V2D(v + 1), n = NEG(p);
				WL& poss = wt[n];
				for (int j = 0; j < poss.size(); j++) {
					assert(poss[j].c_ref != NULL);
					B_REF c = (B_REF)poss[j].c_ref;
					assert(c->size() > 1);
					if (!c->garbage() && c->size() == 2) {
						if (c->orgBin()) { assert(c->status() == UNKNOWN); rawBins.appendFrom(*c, 2); }
						removeClause(c);
					}
				}
				WL& negs = wt[p];
				for (int j = 0; j < negs.size(); j++) {
					assert(negs[j].c_ref != NULL);
					B_REF c = (B_REF)negs[j].c_ref;
					assert(c->size() > 1);
					if (!c->garbage() && c->size() == 2) {
						if (c->orgBin()) { assert(c->status() == UNKNOWN); rawBins.appendFrom(*c, 2); }
						removeClause(c);
					}
				}
			}
			CHECK(cudaMemcpyAsync(cuMem.litsLEA(), rawBins, rawBins.size() * sizeof(uint32), cudaMemcpyHostToDevice, streams[0]));
			gcr.recycle();
		}
	}
	else {
		rawBins.reserve((bins.size() << 1));
		for (int i = 0; i < bins.size(); i++) {
			assert(bins[i]->orgBin());
			assert(bins[i]->status() == UNKNOWN);
			rawBins.appendFrom(*bins[i], 2);
			collectClause(bins[i], false);
			assert(bins[i] == NULL);
		}
		CHECK(cudaMemcpyAsync(cuMem.litsLEA(), rawBins, rawBins.size() * sizeof(uint32), cudaMemcpyHostToDevice, streams[0]));
		bins.clear(true);
	}
	wt.clear(true);
	cnf_stats.global_n_bins = uint32(rawBins.size() >> 1);
}

void ParaFROST::GPU_CNF_fill() {
	// append binaries
	extractBins();
	assert(nBins() == (rawBins.size() >> 1));
	// append k-size clauses 
	rawOrgs.reserve(nLiterals());
	for (int i = 0; i < orgs.size(); i++) rawOrgs.appendFrom(*orgs[i], orgs[i]->size());
	// sync binaries copy
	CHECK(cudaStreamSynchronize(streams[0])); 
	uint32* dBin = cuMem.litsLEA(), *dOrg = dBin + rawBins.size();
	for (uint32 i = 0; i < nBins(); i++) (*cnf)[i].attach(dBin + (i << 1), 2);
	CHECK(cudaMemcpyAsync(dOrg, rawOrgs, rawOrgs.size() * sizeof(uint32), cudaMemcpyHostToDevice, streams[1]));
	// free learnts
	for (int i = 0; i < learnts.size(); i++) { collectClause(learnts[i], false); assert(learnts[i] == NULL); }
	learnts.clear();
	// sync orgs copy
	CHECK(cudaStreamSynchronize(streams[1]));
	cnf_stats.global_n_lits = 0;
	for (int i = 0; i < orgs.size(); i++) {
		(*cnf)[i + cnf_stats.global_n_bins].attach(dOrg + cnf_stats.global_n_lits, orgs[i]->size());
		cnf_stats.global_n_lits += orgs[i]->size();
	}
	cnf_stats.global_n_lits += rawBins.size();
	cnf_stats.global_n_cls = nBins() + orgs.size();
	cnf_stats.global_n_bins = 0;
	cnf->resize(cnf_stats.global_n_cls), cnf->resizeData(cnf_stats.global_n_lits);
	cuMem.prefetchCNF(streams[0]);
	// free the rest
	rawBins.clear(true), rawOrgs.clear(true);
	for (int i = 0; i < orgs.size(); i++) { collectClause(orgs[i], false); assert(orgs[i] == NULL); }
	orgs.clear(true);
	// calc clauses signatures
	calcSigAsync(cnf, 0, nClauses(), streams[0]);
}

bool ParaFROST::awaken()
{
	_simplify();
	createStreams();
	pv = new PV();
	// TODO make an error code for return to handle the original cnf on host
	uint32 numCls = nClauses() + nBins();
	uint64 numLits = nLiterals() + (nBins() << 1ULL);
	if (phases) {
		cnf_stats.max_added_cls = numCls, cnf_stats.max_added_lits = numLits;
		if (verbose > 1) PFLOG(" Max. added clauses/literals = %d/%lld", maxAddedCls(), maxAddedLits());
		numCls += cnf_stats.max_added_cls, numLits += cnf_stats.max_added_lits;
	}
	if (!cuMem.allocPV(pv) || !cuMem.resizeCNF(cnf, numCls, numLits)) return false;
	GPU_CNF_fill();
	histogram.resize(nDualVars());
	d_hist = thrust::raw_pointer_cast(histogram.data()), h_hist = new uint32[nDualVars()];
	cuMem.calcOTBlocks();
	GPU_occurs(nLiterals(), true);
	if (!cuMem.resizeOTAsync(ot, d_hist, numLits)) return false;
	simpVal.resize(nOrgVars(), UNDEFINED);
	for (int i = 0; i < removed.size(); i++) simpVal[V2X(removed[i])] = DEFINED;
	sp->reset_trail();
	printPStats();
	return true;
}

void ParaFROST::GPU_preprocess()
{
	if (!phases && !cls_en) return;
	/********************************/
	/*         awaken sigma         */
	/********************************/
	timer->start();
	if (!awaken()) { timer->stop(), pre_en = false; return; }
	/********************************/
	/*      1st-stage reduction     */
	/********************************/
	if (p_cnf_en) cnf->print();
	int64 lits_before = nLiterals(), lits_diff = INT64_MAX;
	int phase = 0;
	while (lits_diff > LIT_REM_THR && phase < phases) {
		if (interrupted()) killSolver();
		if (verbose > 1) printf("c |\t\tPhase-%d Variable Elections\n", phase);
		CHECK(cudaStreamSynchronize(streams[0]));
		createOT(cnf, ot, p_ot_en);
		if (prop() == UNSAT) killSolver(UNSAT); 
		if (!LCVE()) goto writeBack;
		if (phase + 1 != phases && ((phase + 1) % cnf_free_freq == 0)) {
			cnf_stats.max_added_cls = nClauses(), cnf_stats.max_added_lits = nLiterals();
			if (verbose > 1) printf("c | Max. added clauses/literals = %d/%lld\n", maxAddedCls(), maxAddedLits());
			if (!cuMem.resizeCNF(cnf, nClauses() + maxAddedCls(), nLiterals() + maxAddedLits())) goto writeBack;
			createOTAsync(cnf, ot, p_ot_en);
		}
		GPU_VE(); 
		countCls(cnf, pv->gsts);
		cnf_stats.global_n_cls = cnf_stats.n_cls_after, cnf_stats.global_n_lits = cnf_stats.n_lits_after;
		lits_diff = lits_before - nLiterals(), lits_before = nLiterals();
		if (pv->numGUnits) CHECK(cudaMemcpyAsync(sp->trail + sp->trail_head, cuMem.unitsLEA(),
			pv->numGUnits * sizeof(uint32), cudaMemcpyDeviceToHost, streams[3]));
		cuMem.resetOTCapAsync(ot, streams[0]), GPU_occurs(nLiterals());
		if (!cuMem.resizeOTAsync(ot, d_hist, nLiterals(), streams[0])) goto writeBack;
		phase++, pv->mu_inc++;
	}
	if (pv->numGUnits || cls_en) createOT(cnf, ot, p_ot_en);
	if (pv->numGUnits && prop() == UNSAT) killSolver(UNSAT); 
	/********************************/
	/*      2nd-stage reduction     */
	/********************************/
	if (cls_en) {
		if (verbose > 1) printf("c | 2nd-Stage Clause Eliminations..\n");
		int t_p = mu_pos, t_n = mu_neg;
		while (t_p <= CE_POS_LMT && t_n <= CE_NEG_LMT) pv->mu_inc++, t_p <<= pv->mu_inc, t_n <<= pv->mu_inc;
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
	h_cnf = cuMem.allocHostCNF();
	size_t copySize = h_cnf->numLits(), copyTile = copySize * litsCopyR;
	CHECK(cudaMemcpyAsync(h_cnf->data(), cuMem.litsLEA(), copyTile * sizeof(uint32), cudaMemcpyDeviceToHost, streams[0]));
	countDelVars(cnf, pv->gsts, cuMem.dvarsLEA());
	CHECK(cudaMemcpyAsync(h_cnf->data(copyTile), cuMem.litsLEA(copyTile), (copySize - copyTile) * sizeof(uint32), cudaMemcpyDeviceToHost, streams[1]));
	cnf_stats.global_n_cls = 0;
	cnf_stats.global_n_lits = 0;
	if (cnf_stats.n_del_vars - nVarsDeleted()) { // map
		mapped = true;
		mappedVars.resize(nOrgVars() + 1, 0);
		reverseVars.resize(nOrgVars() + 1, 0);
		cnf_stats.global_n_del_vars = 0;
		cnf_stats.n_org_vars -= cnf_stats.n_del_vars;
		cnf_stats.n_dual_vars = V2D(cnf_stats.n_org_vars + 1LL);
		cleanSlate();
	}
	else { assert(mapped == false); cleanSlate(); }
	printPStats();
	timer->stop();
	timer->pre = timer->cpuTime();
	if (!solve_en) killSolver();
	PDMInit(); // call initial PDM again 
	pre_en = false;
	if (pre_delay) progRate += (initProgRate >> 2);
}

void ParaFROST::cleanSlate()
{
	assert(wt.empty());
	wt.allocMem();
	uint32 copySize = h_cnf->size(), tile = clsTile;
	uint32 nTiles = (copySize + tile - 1) / tile;
	assert(nTiles > 0);
	assert(tile * (nTiles - 1) <= copySize);
	uint32 lastTile = copySize - tile * (nTiles - 1);
	size_t* tiles = new size_t[nTiles], *t = tiles, *t_end = t + nTiles - 1, off = 0;
	while (t < t_end) *t++ = tile * sizeof(SCLAUSE); 
	*t = lastTile * sizeof(SCLAUSE);
	assert(t - tiles == nTiles - 1);
	t = tiles, t_end = t + nTiles;
	CHECK(cudaStreamSynchronize(streams[0]));
	CHECK(cudaStreamSynchronize(streams[1]));
	// issue initial chunck
	CHECK(cudaMemcpy(*h_cnf, cuMem.clsLEA(), *t++, cudaMemcpyDeviceToHost));
	if (mapped) {
		while (t < t_end) {
			CHECK(cudaMemcpyAsync(*h_cnf + off + tile, cuMem.clsLEA(off + tile), *t++, cudaMemcpyDeviceToHost));
			mapTile(off, tile), off += tile;
			if (t < t_end) CHECK(cudaStreamSynchronize(0));
		}
		removed.appendFrom(sp->trail, sp->trail_size, 1);
		delete[] h_hist, tiles; delete pv;
		pv = NULL, h_hist = NULL, d_hist = NULL;
		allocSolver(true), resetSolver();
		uint32 oldVars = cnf_stats.n_org_vars + cnf_stats.n_del_vars;
		Vec<Byte> seen(oldVars + 1, 0);
		if (nTiles > 1) CHECK(cudaStreamSynchronize(0)); // sync last tile
		CHECK(cudaMemcpyAsync(seen, cuMem.seenLEA(), oldVars + 1, cudaMemcpyDeviceToHost));
		mapTile(off, lastTile);
		mappedVars.clear(true);
		cuMem.freeHostCNF(), h_cnf = NULL;
		CHECK(cudaStreamSynchronize(0));
		for (uint32 v = 1; v <= oldVars; v++) if (simpVal[v - 1] == UNDEFINED && !seen[v]) removed.push(V2D(v));
		seen.clear(true);
	}
	else {
		while (t < t_end) {
			CHECK(cudaMemcpyAsync(*h_cnf + off + tile, cuMem.clsLEA(off + tile), *t++, cudaMemcpyDeviceToHost));
			copTile(off, tile), off += tile;
			if (t < t_end) CHECK(cudaStreamSynchronize(0));
		}
		delete[] h_hist, tiles; delete pv;
		pv = NULL, h_hist = NULL, d_hist = NULL;
		resetSolver(false);
		if (nTiles > 1) CHECK(cudaStreamSynchronize(0)); // sync last tile
		copTile(off, lastTile);
		cuMem.freeHostCNF(), h_cnf = NULL;
	}
	assert(consistent(orgs, wt));
	simpVal.clear(true);
	histogram.clear(), rawLits.clear();
}

void ParaFROST::GPU_VO()
{
	assert(d_hist != NULL);
	if (pv->numGUnits) calcVarScores(pv->scores, d_hist, ot); // update d_hist & calc scores
	else calcVarScores(pv->scores, d_hist);
	CHECK(cudaMemcpyAsync(h_hist, d_hist, nDualVars() * sizeof(uint32), cudaMemcpyDeviceToHost, streams[2]));
	thrust::sort(pv->scores, pv->scores + nOrgVars(), VAR_CMP());
	CHECK(cudaStreamSynchronize(streams[2]));
	pv->numGUnits = 0;
	if (verbose == 4) {
		printf("c | Variable occurrences:\n");
		cudaMemcpy(h_hist, d_hist, nDualVars() * sizeof(uint32), cudaMemcpyDeviceToHost);
		for (uint32 v = 0; v < nOrgVars(); v++) {
			uint32 x = pv->scores[v].v, p = V2D(x), n = NEG(p);
			printf("c | Var[%d]->(v: %d, p: %d, n: %d, s: %d)\n", v, x, h_hist[p], h_hist[n], pv->scores[v].sc);
		}
	}
}

void ParaFROST::GPU_occurs(const int64& numLits, const bool& init)
{
	assert(numLits > 0);
	rawLits.resize(numLits);
	if (init) copy(thrust::raw_pointer_cast(rawLits.data()), cnf, numLits);
	else {
		copyIf(thrust::raw_pointer_cast(rawLits.data()), cnf, pv->gsts);
		assert(pv->gsts->numLits == numLits);
	}
	GPU_hist();
	rawLits.clear();
}

void ParaFROST::GPU_hist()
{
	assert(rawLits.size() > 0);
	assert(histogram.size() == V2D(nOrgVars() + 1ULL));
	thrust::sort(rawLits.begin(), rawLits.end());
	thrust::counting_iterator<size_t> search_begin(0);
	thrust::upper_bound(rawLits.begin(), rawLits.end(), search_begin, search_begin + nDualVars(), histogram.begin());
	thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());
}

CNF_STATE ParaFROST::prop()
{
	if (pv->numGUnits) {
		CHECK(cudaStreamSynchronize(streams[3]));
		if (verbose == 4) printTrail();
	}
	while (sp->trail_head < sp->trail_size) { // propagate units
		uint32 unit = sp->trail[sp->trail_head++], unitIdx = V2X(unit), f_unit = FLIP(unit);
		assert(unit);
		LIT_ST val = !ISNEG(unit);
		if (simpVal[unitIdx] != val) { // not propagated before
			assert(simpVal[unitIdx] == UNDEFINED);
			simpVal[unitIdx] = val;
			if (verbose == 4) printf("c | Propagating unit("), printLit(unit), printf("):\n"), ot->printClauseSet(*cnf, unitIdx);
			OL& ol = (*ot)[unit], &ol_f = (*ot)[f_unit];
			for (uint32 i = 0; i < ol.size(); i++) (*cnf)[ol[i]].markDeleted(); // remove satisfied
			for (uint32 i = 0; i < ol_f.size(); i++) { // reduce unsatisfied 
				SCLAUSE& c = (*cnf)[ol_f[i]];
				assert(c.size());
				if (c.status() == DELETED || propClause(c, f_unit)) continue; // clause satisfied
				if (c.size() == 0) return UNSAT; // conflict on top level
				if (c.size() == 1) {
					assert(*c > 1);
					if (simpVal[V2X(*c)] == UNDEFINED) sp->trail[sp->trail_size++] = *c; // new unit
					else return UNSAT;  // conflict on top level
				}
			}
			if (verbose == 4) ot->printClauseSet(*cnf, unitIdx);
			(*ot)[unit].clear(true), (*ot)[f_unit].clear(true);
		}
	}
	if (pv->numGUnits) {
		if (verbose > 1) { evalReds(cnf, pv->gsts); printf("c |\t\tBCP Reductions\n"); logReductions(); }
		pv->units->clear(), filterPVs(), reduceOTAsync(cnf, ot, p_ot_en);
	}
	return UNSOLVED;
}

void ParaFROST::GPU_VE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating variables for round-0." << endl;
	ve(cnf, ot, pv);
	if (verbose > 1) { evalReds(cnf, pv->gsts); printf("c |\t\tBVE Reductions\n"); logReductions(); }
	if (ve_plus_en) {
		countLits(cnf, pv->gsts);
		int64 lits_before = cnf_stats.n_lits_after, lits_removed = nLiterals() - cnf_stats.n_lits_after;
		while (lits_removed > LIT_REM_THR && filterPVs()) {
			// HSE
			if (verbose > 1) printf("c | Applying HSE..");
			hse(cnf, ot, pv);
			// count remaining literals
			countLits(cnf, pv->gsts);
			lits_removed = lits_before - cnf_stats.n_lits_after;
			if (verbose > 1) printf("(Literals reduction: -%lld) ==> done\n", lits_removed);
			if (lits_removed <= LIT_REM_THR) break;
			// VE
			if (verbose > 1) printf("c | Applying VE..");
			ve(cnf, ot, pv);
			// count remaining literals
			countLits(cnf, pv->gsts);
			lits_removed = (lits_before > cnf_stats.n_lits_after) ? lits_before - cnf_stats.n_lits_after : cnf_stats.n_lits_after - lits_before;
			if (verbose > 1) {
				if (lits_before > cnf_stats.n_lits_after) printf("(Literals reduction: -%lld) ==> done\n", lits_removed);
				else printf("(Literals reduction: +%lld) ==> done\n", lits_removed);
			}
			lits_before = cnf_stats.n_lits_after;
		}
	}
	pv->numGUnits = pv->units->size(), sp->trail_size += pv->numGUnits;
}

void ParaFROST::GPU_SUB()
{
	if (interrupted()) killSolver();
	if (verbose > 1) printf("c | HSE-ing parallel variables..");
	hse(cnf, ot, pv);
	pv->numGUnits = pv->units->size(), sp->trail_size += pv->numGUnits;
	if (pv->numGUnits) CHECK(cudaMemcpyAsync(sp->trail + sp->trail_head, cuMem.unitsLEA(),
		pv->numGUnits * sizeof(uint32), cudaMemcpyDeviceToHost, streams[3]));
	if (prop() == UNSAT) killSolver(UNSAT);
	if (verbose > 1) {
		printf(" ==> done\n");
		CHECK(cudaDeviceSynchronize()), evalReds(cnf, pv->gsts);
		printf("c |\t\tHSE Reductions\n"), logReductions();
	}
}

void ParaFROST::GPU_BCE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating blocked clauses..";
	bce(cnf, ot, pv);
	if (verbose > 1)
		printf(" ==> done\n"), countCls(cnf, pv->gsts), printf("c |\t\tBCE Reductions\n"), logReductions(1);
}

void ParaFROST::GPU_HRE()
{
	if (interrupted()) killSolver();
	if (verbose > 1) cout << "c | Eliminating hidden redundances..";
	hre(cnf, ot, pv);
	if (verbose > 1) 
		printf(" ==> done\n"), countCls(cnf, pv->gsts), printf("c |\t\tHRE Reductions\n"), logReductions(1);
}

void ParaFROST::depFreeze(const OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	for (uint32 i = 0; i < ol.size(); i++) {
		SCLAUSE& c = (*cnf)[ol[i]];
		for (LIT_POS k = 0; k < c.size(); k++) {
			register uint32 v = ABS(c[k]), p = V2D(v), n = NEG(p);
			if (v != cand && (h_hist[p] < p_temp || h_hist[n] < n_temp)) sp->frozen[v] = true;
		}
	}
}

bool ParaFROST::LCVE()
{
	GPU_VO();
	pv->numPVs = 0, pv->pVars->clear();
	for (uint32 v = 0; v < nOrgVars(); v++) {
		uint32 cand = pv->scores[v].v;
		assert(cand && cand <= nOrgVars());
		if (sp->frozen[cand]) continue;
		uint32 p = V2D(cand), n = NEG(p);
		assert((*ot)[p].size() == h_hist[p]);
		assert((*ot)[n].size() == h_hist[n]);
		if (h_hist[p] == 0 && h_hist[n] == 0) continue;
		assert(simpVal[cand - 1] == UNDEFINED);
		uint32 pos_temp = mu_pos << pv->mu_inc, neg_temp = mu_neg << pv->mu_inc;
		if (h_hist[p] >= pos_temp && h_hist[n] >= neg_temp) break;
		pv->pVars->_push(cand);
		depFreeze((*ot)[p], cand, pos_temp, neg_temp);
		depFreeze((*ot)[n], cand, pos_temp, neg_temp);
	}
	pv->numPVs = pv->pVars->size();
	assert(verifyLCVE());
	if (verbose >= 3) printf("c | PLCVs "), printVars(*pv->pVars, pv->numPVs);
	bool* f = sp->frozen, * f_end = f + nOrgVars() + 1;
	while (f != f_end) *f++ = 0;
	if (pv->numPVs < MIN_PVARS) {
		if (verbose > 1) PFLOGW("parallel variables not enough -> skip simp.");
		return false;
	}
	return true;
}