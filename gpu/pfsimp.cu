/***********************************************************************[pfsimp.cu]
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

#include "pfsimpopts.h"
#include "pfsolve.h"
#include "pfsort.h"

namespace pFROST {

	using namespace SIGmA;

	void ParaFROST::masterFree()
	{
		if (h_hist != NULL) delete[] h_hist;
		if (cnf != NULL) cnf->~CNF(), cnf = NULL;
		if (ot != NULL) ot->~OT(), ot = NULL;
		if (vars != NULL) vars->~VARS(), vars = NULL;
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
		assert(sigma_en || sigma_live_en);
		ngpus = opt_gpus;
		nstreams = opt_streams;
		solve_en = opt_solve_en;
		ve_en = opt_ve_en || opt_ve_plus_en;
		ve_plus_en = opt_ve_plus_en;
		sub_en = opt_sub_en;
		bce_en = opt_bce_en;
		hre_en = opt_hre_en;
		all_en = opt_all_en;
		phases = opt_phases;
		mu_pos = opt_mu_pos;
		mu_neg = opt_mu_neg;
		lcve_min = opt_lcve_min;
		lits_rem_min = opt_lits_rem_min;
		shrink_rate = opt_cnf_free;
		cls_en = all_en || sub_en || bce_en || hre_en;
		if (all_en) ve_en = 1, ve_plus_en = 1, sub_en = 1, bce_en = 1, hre_en = 1;
		if (!phases && ve_en) phases = 1; // at least 1 phase needed for BVE(+)
		if (phases && !ve_en) phases = 0;
		if (ngpus > devCount) ngpus = devCount;
	}

	void ParaFROST::extract(BCNF& bcnf)
	{
		for (uint32 i = 0; i < bcnf.size(); i++) {
			CLAUSE& c = cm[bcnf[i]];
			assert(!c.deleted());
			cnf->newClause(c, c.size());
			inf.nClauses++, inf.nLiterals += c.size();
		}
		bcnf.clear(true);
	}

	void ParaFROST::awaken()
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
			reduce();
			PFLOGN2(2, " Shrinking CNF before eliminations..");
			shrinkWT();
			shrink(orgs);
			shrink(learnts);
			inf.nDelVars += (trail.size() - sp->simplified);
			PFLENDING(2, 5, " (-%d variables)", trail.size() - sp->simplified);
			sp->simplified = trail.size();
		}
		initSimp();
		// alloc simplifier memory 
		uint32 numCls = maxOrgs() + maxLearnts(), numLits = maxLiterals();
		if (phases) {
			inf.maxAddedCls = maxOrgs(), inf.maxAddedLits = maxOrgLits();
			PFLOG2(2, " Maximum added clauses/literals = %d/%d", inf.maxAddedCls, inf.maxAddedLits);
			numCls += inf.maxAddedCls, numLits += inf.maxAddedLits;
		}
		if (!cuMem.allocPV(vars) || !cuMem.resizeCNF(cnf, numCls, numLits)) { errCode = CNFALLOC_FAIL; return; }
		createStreams();
		PFLOGN2(2, " Extracting clauses..");
		BCNF bins;
		bins.reserve(inf.nOrgBins + inf.nLearntBins);
		wtBin.recycle();
		extractBins(bins);
		inf.nClauses = inf.nLiterals = 0;
		inf.nOrgBins = inf.nLearntBins = 0;
		extract(bins);
		// TODO: extract clauses using WT then use orgs & learnts
		extract(orgs);
		extract(learnts);
		PFLENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
		// prefetch cnf 
		cuMem.prefetchCNF(streams[0]);
		// initialize vstates on device with the current vstates on host
		cacheVarSt(cudaMemcpyHostToDevice, streams[1]);
		// cleaning up
		wtBin.clear(true), wt.clear(true), cm.destroy();
		// compute clauses signatures
		calcSigCNFAsync(cnf, 0, inf.nClauses, streams[0]);
		// compute histogram on GPU & alloc OT
		histogram.resize(inf.nDualVars);
		d_hist = thrust::raw_pointer_cast(histogram.data()), h_hist = new uint32[inf.nDualVars];
		if (!reallocOT()) return;
		printPStats();
		sync(streams[1]);
	}

	void ParaFROST::varReorder()
	{
		PFLOGN2(2, " Finding eligible variables for LCVE..");
		assert(d_hist != NULL);
		if (vars->nUnits) calcScores(vars, d_hist, ot); // update d_hist & calc scores
		else calcScores(vars, d_hist);
		CHECK(cudaMemcpyAsync(h_hist, d_hist, inf.nDualVars * sizeof(uint32), cudaMemcpyDeviceToHost, streams[2]));
		thrust::sort(vars->eligible, vars->eligible + inf.maxVar, LCV_CMP(vars->scores));
		sync(streams[2]);
		PFLDONE(2, 5);
		vars->nUnits = 0;
		if (verbose == 4) {
			PFLOG0(" Eligible variables:");
			cudaMemcpy(h_hist, d_hist, inf.nDualVars * sizeof(uint32), cudaMemcpyDeviceToHost);
			for (uint32 v = 0; v < inf.maxVar; v++) {
				uint32 x = vars->eligible[v], p = v2l(x), n = neg(p);
				PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, h_hist[p], h_hist[n], vars->scores[x]);
			}
		}
	}

	void ParaFROST::calcOccurs(const uint32& numLits)
	{
		assert(numLits > 0);
		rawLits.resize(numLits);
		copyIf(thrust::raw_pointer_cast(rawLits.data()), cnf, vars->gstats);
		assert(vars->gstats->numLits == numLits);
		histSimp();
		rawLits.clear();
	}

	void ParaFROST::histSimp()
	{
		assert(rawLits.size() > 0);
		assert(histogram.size() == v2l(inf.maxVar + 1ULL));
		thrust::sort(rawLits.begin(), rawLits.end());
		thrust::counting_iterator<size_t> search_begin(0);
		thrust::upper_bound(rawLits.begin(), rawLits.end(), search_begin, search_begin + inf.nDualVars, histogram.begin());
		thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());
	}

	void ParaFROST::preprocess()
	{
		if (!phases && !cls_en) return;
		/********************************/
		/*         awaken sigma         */
		/********************************/
		assert(cnfstate == UNSOLVED);
		assert(conflict == NOREF);
		timer.start();
		int phase;
		int64 lits_before, lits_diff;
		awaken();
		if (errCode == CNFALLOC_FAIL) { cuMem.destroy(); return; }
		if (errCode == OTALLOC_FAIL) { cuMem.destroy(); goto writeBack; }
		if (interrupted()) killSolver();
		/********************************/
		/*      1st-stage reduction     */
		/********************************/
		phase = 0, lits_before = inf.nLiterals, lits_diff = INT64_MAX;
		while (lits_diff > lits_rem_min && phase < phases) {
			if (interrupted()) killSolver();
			sync(streams[0]);
			createOT(cnf, ot, false);
			if (!prop()) killSolver();
			PFLOG2(2, "\t\tPhase-%d Variable Elections", phase);
			if (!LCVE()) goto writeBack;
			if (!reallocCNF(phase + 1)) goto writeBack;
			VE();
			countCls();
			inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
			lits_diff = lits_before - inf.nLiterals, lits_before = inf.nLiterals;
			cacheUnits(streams[3]);
			if (!reallocOT(streams[0])) goto writeBack;
			phase++, vars->mu_inc++;
		}
		if (vars->nUnits || cls_en) createOT(cnf, ot, false);
		if (!prop()) killSolver();
		/********************************/
		/*      2nd-stage reduction     */
		/********************************/
		if (cls_en) {
			if (verbose > 1) printf("c | 2nd-Stage Clause Eliminations..\n");
			int t_p = mu_pos, t_n = mu_neg;
			while (t_p <= CE_POS_LMT && t_n <= CE_NEG_LMT) vars->mu_inc++, t_p <<= vars->mu_inc, t_n <<= vars->mu_inc;
			if (!LCVE()) goto writeBack;
			if (sub_en) SUB();
			if (bce_en) BCE();
			if (hre_en) HRE();
		}
		/********************************/
		/*           Write Back         */
		/********************************/
	writeBack:
		if (interrupted()) killSolver();
		cacheCNF(streams[0], streams[1]);
		countDelVars();
		printPStats();
		sigmified = true;
		if (canMap()) map();
		else assert(!mapped), newBeginning();
		timer.stop();
		timer.pre = timer.cpuTime();
		if (!solve_en) killSolver();
	}

	C_REF ParaFROST::newClause(SCLAUSE& s)
	{
		assert(!s.deleted());
		C_REF r = s.ref();
		if (r == NOREF) {
			int sz = s.size();
			assert(sz > 1);
			r = cm.alloc(sz);
			s.set_ref(r); // new ref overwrites simplified clause sig.
			CLAUSE& new_c = cm[r];
			if (mapped) vmap.mapClause(new_c, s);
			else new_c.copyLitsFrom(s);
			assert(sz == new_c.size());
			assert(new_c[0] > 0 && new_c[1] > 0);
			assert(new_c[0] <= UINT32_MAX && new_c[1] <= UINT32_MAX);
			new_c.set_status(s.status());
			if (sz == 2) {
				if (s.status() == ORIGINAL) inf.nOrgBins++;
				else assert(s.status() == LEARNT), inf.nLearntBins++;
			}
			else {
				assert(sz > 2);
				if (s.status() == LEARNT) {
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

	void ParaFROST::newBeginning() {
		assert(sigma_en || sigma_live_en);
		assert(sigmified);
		assert(!hcnf->empty());
		assert(wtBin.empty()), assert(wt.empty());
		assert(orgs.empty()), assert(learnts.empty());
		inf.nOrgBins = inf.nLearntBins = 0;
		inf.nLiterals = inf.nLearntLits = 0;
		assert(inf.maxVar > vmap.numVars());
		uint32 tableSize = mapped ? v2l(vmap.size()) : inf.nDualVars;
		wtBin.resize(tableSize), wt.resize(tableSize);
		cm.init(hcnf->data().size);
		sync(streams[0]), sync(streams[1]); // sync CNF caching
		if (!mapped) cacheVarSt(cudaMemcpyDeviceToHost);
		createWT(), copyWatched(), copyNonWatched(); // must follow this order
		inf.nOrgCls = inf.nClauses = orgs.size();
		inf.nOrgLits = inf.nLiterals;
		cuMem.breakMirror(); // free host CNF
		if (!mapped) assert(vmap.empty()), sync(), sp->lockMelted();
	}

	void ParaFROST::depFreeze(const OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
	{
		for (uint32 i = 0; i < ol.size(); i++) {
			SCLAUSE& c = (*cnf)[ol[i]];
			for (int k = 0; k < c.size(); k++) {
				register uint32 v = l2a(c[k]), p = v2l(v), n = neg(p);
				if (v != cand && (h_hist[p] < p_temp || h_hist[n] < n_temp)) sp->frozen[v] = 1;
			}
		}
	}

	bool ParaFROST::LCVE()
	{
		// reorder variables
		varReorder();
		// extended LCVE
		PFLOGN2(2, " Electing variables..");
		vars->numPVs = 0, vars->pVars->clear();
		for (uint32 v = 0; v < inf.maxVar; v++) {
			uint32 cand = vars->eligible[v];
			assert(cand && cand <= inf.maxVar);
			if (sp->vstate[cand] == FROZEN) continue;
			if (sp->frozen[cand]) continue;
			uint32 p = v2l(cand), n = neg(p);
			assert((*ot)[p].size() == h_hist[p]);
			assert((*ot)[n].size() == h_hist[n]);
			if (h_hist[p] == 0 && h_hist[n] == 0) continue;
			uint32 pos_temp = mu_pos << vars->mu_inc, neg_temp = mu_neg << vars->mu_inc;
			if (h_hist[p] >= pos_temp && h_hist[n] >= neg_temp) break;
			assert(sp->vstate[cand] == ACTIVE);
			vars->pVars->_push(cand);
			depFreeze((*ot)[p], cand, pos_temp, neg_temp);
			depFreeze((*ot)[n], cand, pos_temp, neg_temp);
		}
		vars->numPVs = vars->pVars->size();
		assert(verifyLCVE());
		PFLENDING(2, 5, "(%d elected)", vars->numPVs);
		if (verbose == 4) { PFLOGN0(" PLCVs "); printVars(*vars->pVars, vars->numPVs, 'v'); }
		memset(sp->frozen, 0, inf.maxVar + 1ULL);
		if (vars->numPVs < lcve_min) {
			if (verbose > 1) PFLOGW("parallel variables not enough -> skip BVE");
			return false;
		}
		return true;
	}

	bool ParaFROST::propClause(SCLAUSE& c, const uint32& f_assign)
	{
		uint32 sig = 0;
		int n = 0;
		bool check = false;
		for (int k = 0; k < c.size(); k++) {
			uint32 lit = c[k];
			if (lit != f_assign) {
				if (isTrue(lit)) return true;
				c[n++] = lit;
				sig |= MAPHASH(lit);
			}
			else check = true;
		}
		assert(check);
		assert(n == c.size() - 1);
		assert(c.hasZero() < 0);
		assert(c.isSorted());
		c.set_sig(sig);
		c.pop();
		return false;
	}

	bool ParaFROST::enqeueCached(const cudaStream_t& stream) {
		if (vars->nUnits) {
			nForced = sp->propagated;
			sync(stream);
			assert(vars->cachedUnits != NULL);
			uint32* t = vars->cachedUnits + vars->nUnits;
			for (uint32* u = vars->cachedUnits; u != t; u++) {
				LIT_ST val = value(*u);
				if (val == UNDEFINED) enqueue(*u);
				else if (!val) return false; // early conflict detection
			}
			if (trail.size() == sp->propagated) vars->nUnits = nForced = 0; // duplicate units
			else PFLTRAIL(this, 3);
		}
		return true;
	}

	bool ParaFROST::prop()
	{
		enqeueCached(streams[3]);
		while (sp->propagated < trail.size()) { // propagate units
			uint32 assign = trail[sp->propagated++], f_assign = flip(assign);
			assert(assign);
			PFLBCP(this, 4, assign);
			OL& ol = (*ot)[assign], & f_ol = (*ot)[f_assign];
			for (uint32 i = 0; i < ol.size(); i++) (*cnf)[ol[i]].markDeleted(); // remove satisfied
			for (uint32 i = 0; i < f_ol.size(); i++) { // reduce unsatisfied 
				SCLAUSE& c = (*cnf)[f_ol[i]];
				assert(c.size());
				if (c.deleted() || propClause(c, f_assign)) continue; // clause satisfied
				assert(c.size()); // cannot be empty at this point
				if (c.size() == 1) {
					assert(*c > 1);
					if (unassigned(*c)) enqueue(*c); 
					else { cnfstate = UNSAT; return false; }  // conflict on top level
				}
			}
			(*ot)[assign].clear(true), (*ot)[f_assign].clear(true);
		}
		cleanProped();
		return true;
	}

	void ParaFROST::VE()
	{
		if (interrupted()) killSolver();
		int64 lits_before = inf.nLiterals, lits_removed = INT64_MAX;
		int round = 0;
		while (lits_removed > lits_rem_min) {
			PFLOG2(2, " Elimination round %d:", round++);
			if (ve_plus_en) {
				PFLOGN2(2, "  1) HSE-ing variables..");
				hse(cnf, ot, vars), countLits();
				lits_removed = lits_before - inf.n_lits_after;
				lits_before = inf.n_lits_after;
				assert(lits_removed >= 0);
				PFLENDING(2, 5, "(Literals removed : %lld)", -lits_removed);
			}
			PFLOGN2(2, "  2) Eliminating variables..");
			ve(cnf, ot, vars), countLits();
			lits_removed = lits_before - inf.n_lits_after;
			PFLENDING(2, 5, "(Literals removed : %c%lld)", lits_removed < 0 ? '+' : '-', abs(lits_removed));
			lits_before = inf.n_lits_after;
			if (filterPVs()) break;
		}
		PFLREDALL(this, 2, "BVE(+) Reductions");
		cacheNumUnits(streams[3]);
	}

	void ParaFROST::SUB()
	{
		if (interrupted()) killSolver();
		PFLOGN2(2, " SUB-ing variables..");
		hse(cnf, ot, vars);
		cacheNumUnits(streams[3]);
		cacheUnits(streams[3]);
		if (!prop()) killSolver();
		PFLDONE(2, 5);
		PFLREDCL(this, 2, "SUB Reductions");
	}

	void ParaFROST::BCE()
	{
		if (interrupted()) killSolver();
		PFLOGN2(2, " Eliminating blocked clauses..");
		bce(cnf, ot, vars);
		PFLDONE(2, 5);
		PFLREDALL(this, 2, "BCE Reductions");
	}

	void ParaFROST::HRE()
	{
		if (interrupted()) killSolver();
		PFLOGN2(2, " Eliminating hidden redundances..");
		hre(cnf, ot, vars);
		PFLDONE(2, 5);
		PFLREDCL(this, 2, "HRE Reductions");
	}

}