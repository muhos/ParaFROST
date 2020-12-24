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

#include "pfsolve.h"
#include <cub/device/device_select.cuh>
using namespace cub;

namespace pFROST {
	using namespace SIGmA;
	cuTIMER *cutimer = NULL;
	bool unified_access = false;
	bool profile_gpu = false;
	bool sync_always = false;
	bool atomic_ve = false;
	bool gc_par = false;

	void ParaFROST::masterFree()
	{
		syncAll();
		cumm.freeFixed();
		cumm.freePinned();
		tca.destroy();
		cleanDynamic();
		destroyStreams();
		cumm.breakMirror(), hcnf = NULL;
	}

	void ParaFROST::slavesFree()
	{

	}

	void ParaFROST::optSimp()
	{
		culimit.limits[0] = opts.hse_limit;
		culimit.limits[1] = opts.bce_limit;
		culimit.limits[2] = opts.ere_limit;
		culimit.limits[3] = opts.xor_max_arity;
		if (opts.ngpus > devCount) opts.ngpus = devCount;
		if (profile_gpu = opts.profile_simp) cutimer = new cuTIMER;
		initLimits(culimit);
	}

	void ParaFROST::calcOccurs(const uint32& numLits)
	{
		assert(numLits);
		PFLOGN2(2, " Copying survived literals..");
		copyIf(cuhist.d_lits, cnf, vars->gstats);
		assert(vars->gstats->numLits == numLits);
		PFLENDING(2, 5, "(%d copied)", numLits);
		histSimp(numLits);
	}

	void ParaFROST::histSimp(const uint32& numLits)
	{
		assert(numLits);
		PFLOGN2(2, " Computing histogram on %d elements..", numLits);
		if (profile_gpu) cutimer->start();
		thrust::sort(thrust::cuda::par(tca), cuhist.thrust_lits, cuhist.thrust_lits + numLits);
		thrust::counting_iterator<size_t> search_begin(0);
		thrust::upper_bound(thrust::cuda::par(tca), cuhist.thrust_lits, cuhist.thrust_lits + numLits, search_begin, search_begin + inf.nDualVars, cuhist.thrust_hist);
		thrust::adjacent_difference(thrust::cuda::par(tca), cuhist.thrust_hist, cuhist.thrust_hist + inf.nDualVars, cuhist.thrust_hist);
		if (profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		PFLDONE(2, 5);
	}

	void ParaFROST::extract(CNF* dest, BCNF& src)
	{
		for (uint32 i = 0; i < src.size(); i++) {
			CLAUSE& c = cm[src[i]];
			if (c.deleted()) continue;
			dest->newClause(c);
			inf.nClauses++, inf.nLiterals += c.size();
		}
	}

	void ParaFROST::awaken(const bool& strict)
	{
		assert(conflict == NOREF);
		assert(cnfstate == UNSOLVED);
		assert(sp->propagated == trail.size());
		initSimp();
		if (stats.sigmifications && strict) reduceTop(strict);
		if (orgs.empty()) { sigState = AWAKEN_FAIL; return; }
		// alloc simplifier memory 
		size_t numCls = maxClauses(), numLits = maxLiterals();
		if (opts.phases) {
			size_t maxAddedCls = opts.ve_en ? orgs.size() : 0;
			size_t maxAddedLits = opts.ve_en ? inf.nLiterals * opts.lits_mul : 0;
			PFLOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
			numCls += maxAddedCls, numLits += maxAddedLits;
		}
		assert(inf.nDualVars);
		if (!cumm.allocHist(cuhist, numLits) ||
			!cumm.allocAux(numCls) ||
			!cumm.allocVars(vars, numLits) ||
			!cumm.resizeCNF(cnf, numCls, numLits))
		{
			sigState = CNFALLOC_FAIL;
			return;
		}
		if (unified_access) {
			PFLOGN2(2, " Extracting clauses directly to device..");
			if (profile_gpu) cutimer->start();
			wt.clear(true);
			printStats(1, '-', CGREEN0), inf.nClauses = inf.nLiterals = 0;
			extract(cnf, orgs), orgs.clear(true);
			extract(cnf, learnts), learnts.clear(true);
			cm.destroy();
			assert(inf.nClauses == cnf->size());
			if (profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
		}
		else {
			PFLOGN2(2, " Extracting clauses heterogeneously to device..");
			if (profile_gpu) cutimer->start(streams[0]);
			wt.clear(true);
			cumm.createMirror(hcnf, maxClauses(), maxLiterals());
			printStats(1, '-', CGREEN0), inf.nClauses = inf.nLiterals = 0;
			extract(hcnf, orgs), reflectCNF(streams[0], streams[1]), orgs.clear(true);
			extract(hcnf, learnts), reflectCNF(streams[0], streams[1]), learnts.clear(true);
			cm.destroy(); 
			sync(streams[0]), sync(streams[1]);
			cumm.resizeCNFAsync(cnf, hcnf->data().size, hcnf->size());
			assert(hcnf->data().size == dataoff);
			assert(inf.nClauses == hcnf->size() && hcnf->size() == csoff);
			if (profile_gpu) cutimer->stop(streams[0]), cutimer->io += cutimer->gpuTime();
		}
		PFLENDING(2, 5, "(%d clauses extracted)", inf.nClauses);
		sync(); // sync device CNF resize
		prepareCNFAsync(cnf, streams[0]);
	}

	void ParaFROST::sigmify()
	{
		/********************************/
		/*		Getting ready...        */
		/********************************/
		if (alldisabled()) return;
		backtrack();
		if (BCP()) { cnfstate = UNSAT; return; }
		shrink(orgs), shrink(learnts);
		/********************************/
		/*         awaken sigma         */
		/********************************/
		timer.stop();
		timer.solve += timer.cpuTime();
		timer.start();
		awaken();
		if (sleeping()) return;
		/********************************/
		/*      reduction phases        */
		/********************************/
		assert(!phase && !mu_inc);
		int64 litsbefore = inf.nLiterals, diff = INT64_MAX;
		while (!interrupted()) {
			if (!reallocOT(streams[0])) break;
			sync(streams[0]);
			if (!reallocCNF()) break;
			createOTAsync(cnf, ot, 0);
			if (!prop()) killSolver();
			if (!LCVE()) break;
			sortOTAsync(cnf, ot, vars, streams);
			if (stop(diff)) { ERE(); break; }
			HSE(), VE(), BCE();
			countAll();
			cacheNumUnits(streams[3]);
			inf.nClauses = inf.n_cls_after, inf.nLiterals = inf.n_lits_after;
			diff = litsbefore - inf.nLiterals, litsbefore = inf.nLiterals;
			cacheUnits(streams[3]);
			phase++, mu_inc++;
			mu_inc += phase == opts.phases;
		}
		/********************************/
		/*          Write Back          */
		/********************************/
		assert(sp->propagated == trail.size());
		cumm.updateMaxCap(), tca.updateMaxCap(); // for reporting GPU memory only
		if (interrupted()) killSolver();
		if (reallocFailed()) syncAll();
		cacheCNF(streams[0], streams[1]);
		stats.sigmifications++;
		lrn.elim_lastmarked = lrn.elim_marked;
		if (ereCls > inf.nClauses) stats.n_reduns += ereCls - inf.nClauses;
		if (inf.maxFrozen > sp->simplified) stats.n_forced += inf.maxFrozen - sp->simplified;
		if (satisfied() || !inf.nClauses)
			cnfstate = SAT, printStats(1, 'p', CGREEN);
		else {
			if (maxInactive() > opts.map_min) map(true);
			else assert(!mapped), newBeginning();
			sigmaDelay();
		}
		timer.stop(), timer.simp += timer.cpuTime();
		if (!opts.solve_en) killSolver();
		timer.start();
	}

	void ParaFROST::varReorder()
	{
		PFLOGN2(2, " Finding eligible variables for LCVE..");
		assert(cuhist.d_hist != NULL);
		// NOTE: OT creation will be synced in calcScores call
		if (vars->nUnits) calcScores(vars, cuhist.d_hist, ot); // update d_hist & calc scores
		else calcScores(vars, cuhist.d_hist);
		cuhist.cacheHist(streams[2]);
		if (profile_gpu) cutimer->start(streams[3]);
		thrust::sort(thrust::cuda::par(tca).on(streams[3]), vars->eligible, vars->eligible + inf.maxVar, GPU_LCV_CMP(vars->scores));
		PFLDONE(2, 5);
		vars->nUnits = 0;
		sync(streams[2]);
		if (profile_gpu) cutimer->stop(streams[3]), cutimer->vo += cutimer->gpuTime();
		if (verbose == 4) {
			PFLOG0(" Eligible variables:");
			for (uint32 v = 0; v < inf.maxVar; v++) {
				uint32 x = vars->eligible[v], p = V2L(x), n = NEG(p);
				PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", v, x, cuhist[p], cuhist[n], vars->scores[x]);
			}
		}
	}

	bool ParaFROST::LCVE()
	{
		// reorder variables
		varReorder();
		// extended LCVE
		PFLOGN2(2, " Electing variables (p-mu: %d, n-mu: %d)..", opts.mu_pos << mu_inc, opts.mu_neg << mu_inc);
		vars->numPVs = 0, vars->pVars->clear();
		OT& ot = *this->ot; // cache 'ot' reference on host
		for (uint32 i = 0; i < inf.maxVar; i++) {
			uint32 cand = vars->eligible[i];
			assert(cand && cand <= inf.maxVar);
			if (sp->vstate[cand] == FROZEN || sp->vstate[cand] == MELTED) continue;
			if (sp->frozen[cand]) continue;
			uint32 p = V2L(cand), n = NEG(p);
			assert(ot[p].size() >= cuhist[p]);
			assert(ot[n].size() >= cuhist[n]);
			if (!cuhist[p] && !cuhist[n]) continue;
			uint32 pos_temp = opts.mu_pos << mu_inc, neg_temp = opts.mu_neg << mu_inc;
			if (cuhist[p] >= pos_temp && cuhist[n] >= neg_temp) break;
			assert(sp->vstate[cand] == ACTIVE);
			vars->pVars->_push(cand);
			depFreeze(ot[p], cand, pos_temp, neg_temp);
			depFreeze(ot[n], cand, pos_temp, neg_temp);
		}
		vars->numPVs = vars->pVars->size();
		assert(verifyLCVE());
		memset(sp->frozen, 0, inf.maxVar + 1ULL);
		if (vars->numPVs) {
			uint32 mcv = vars->pVars->back(), pmcv = V2L(mcv);
			PFLENDING(2, 5, "(%d elected, mcv: %d, pH: %d, nH: %d)", vars->numPVs, mcv, cuhist[pmcv], cuhist[NEG(pmcv)]);
			if (verbose == 4) { PFLOGN0(" PLCVs "); printVars(*vars->pVars, vars->numPVs, 'v'); }
		}
		if (vars->numPVs < opts.lcve_min) {
			PFLDONE(2, 5);
			if (verbose > 1) PFLOGW("parallel variables not enough -> skip SIGmA");
			return false;
		}
		return true;
	}

	inline void ParaFROST::depFreeze(OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
	{
		CNF& cnf = *this->cnf; // cache 'cnf' reference on host
		for (S_REF* i = ol; i != ol.end(); i++) {
			SCLAUSE& c = cnf[*i];
			if (c.deleted()) continue;
			for (uint32* k = c; k != c.end(); k++) {
				register uint32 lit = *k, v = ABS(lit);
				register uint32 p = POS(lit), n = NEG(lit);
				if (v != cand && (cuhist[p] < p_temp || cuhist[n] < n_temp)) sp->frozen[v] = 1;
			}
		}
	}

	void ParaFROST::newBeginning() {
		assert(opts.sigma_en || opts.sigma_live_en);
		assert(wt.empty());
		assert(orgs.empty());
		assert(learnts.empty());
		assert(inf.maxVar > vmap.numVars());
		assert(inf.nClauses <= hcnf->size());
		if (unified_access) assert(hcnf == NULL), hcnf = cnf;
		if (!mapped) assert(vmap.empty()), sp->lockMelted(inf.maxVar);
		assert(!hcnf->empty());
		cm.init(hcnf->data().size);
		if (!unified_access) {
			sync(streams[0]), sync(streams[1]); // sync CNF caching
			if (profile_gpu) cutimer->stop(streams[1]), cutimer->io += cutimer->gpuTime();
		}
		cacheResolved(streams[2]);
		writeBack();
		syncAll();
		if (unified_access) {
			hcnf = NULL;
			if (profile_gpu) cutimer->stop(), cutimer->io += cutimer->gpuTime();
		}
		else cumm.breakMirror(), hcnf = NULL;
		printStats(1, 'p', CGREEN);
	}

	inline void ParaFROST::writeBack() {
		if (!propHost()) killSolver(); // prop any pending units via host memory if 'realloc' failed
		inf.nLiterals = inf.nLearntLits = 0;
		for (uint32 i = 0; i < hcnf->size(); i++) newClause(hcnf->clause(i));
		assert(size_t(orgs.size() + learnts.size()) == inf.nClauses);
		inf.nOrgCls = orgs.size(), inf.nOrgLits = inf.nLiterals;
		wt.resize(mapped ? V2L(vmap.size()) : inf.nDualVars);
		rebuildWT(opts.priorbins_en);
	}

	inline void	ParaFROST::initSimp() {
		cleanDynamic();
		cumm.freeFixed();
		tca.destroy();
		sigState = AWAKEN_SUCC;
		dataoff = csoff = 0;
		phase = mu_inc = 0;
		ereCls = 0;
		nForced = 0;
		compacted = false;
		assert(hcnf == NULL);
		if (stats.sigmifications) {
			// memory must be cleaned before this inquiry
			size_t free = 0, tot = 0;
			CHECK(cudaMemGetInfo(&free, &tot));
			cumm.reset(free);
		}
	}

	inline void	ParaFROST::sigmaDelay() {
		if (opts.sigma_live_en) {
			// update sigma trigger (inspired by Cadical) 
			// but we decrease phases too
			double current_inc = scale(opts.sigma_inc * (phase + 1.0));
			lrn.sigma_conf_max = nConflicts + current_inc;
			PFLOG2(2, " inprocessing limit increased to %lld conflicts by a weight of %.2f", lrn.sigma_conf_max, current_inc);
			if (opts.phases > 2) {
				opts.phases--;
				PFLOG2(2, " inprocessing phases decreased to %d", opts.phases);
			}
		}
	}

	inline void ParaFROST::cacheCNF(const cudaStream_t& s1, const cudaStream_t& s2)
	{
		// NOTE: if there are units propagated at the last phase,
		// deleted clauses will be left (not compacted), 
		// thus cnf->size() or hcnf->size() must always be used after step 1
		// 0) count clauses on GPU and variables on CPU
		countFinal();
		if (!inf.nClauses) return; // all eliminated
		if (sigState == CNFALLOC_FAIL || sigState == OTALLOC_FAIL) opts.aggr_cnf_sort = false;
		if (unified_access) {
			if (profile_gpu) cutimer->start();
			// 1) compact cs w.r.t clause status on gpu
			if (!compacted) {
				uint32* ts = cuhist.d_lits;
				size_t tb = 0;
				DeviceSelect::If(NULL, tb, cnf->csData(), cnf->csData(), ts, cnf->size(), COMPACT_CMP(cnf)), assert(tb <= cumm.literalsCap());
				DeviceSelect::If(ts + 1, tb, cnf->csData(), cnf->csData(), ts, cnf->size(), COMPACT_CMP(cnf), s1);
				sync();
				cnf->resize(inf.nClauses); 
			}
			// 2) sort cs w.r.t clause size on gpu (user-enabled)
			if (opts.aggr_cnf_sort) thrust::stable_sort(thrust::cuda::par(tca).on(s1), cnf->csData(), cnf->csData() + cnf->size(), CNF_CMP_KEY(cnf));
		}
		else {
			cumm.mirrorCNF(hcnf);
			if (profile_gpu) cutimer->start(s2);
			// 1) compact cs w.r.t  clause status on gpu
			if (!compacted) {
				uint32* ts = cuhist.d_lits;
				size_t tb = 0;
				DeviceSelect::If(NULL, tb, cumm.csMem(), cumm.csMem(), ts, hcnf->size(), COMPACT_CMP(cnf)), assert(tb <= cumm.literalsCap());
				DeviceSelect::If(ts + 1, tb, cumm.csMem(), cumm.csMem(), ts, hcnf->size(), COMPACT_CMP(cnf), s1);
				hcnf->resize(inf.nClauses);
			}
			// 2) copy actual cnf data async.
			CHECK(cudaMemcpyAsync(hcnf->data().mem, cumm.cnfMem(), hcnf->data().size * hc_bucket, cudaMemcpyDeviceToHost, s2));
			// 3) sort cs w.r.t clause size on gpu (user-enabled)
			if (opts.aggr_cnf_sort) thrust::stable_sort(thrust::cuda::par(tca).on(s1), cumm.csMem(), cumm.csMem() + hcnf->size(), CNF_CMP_KEY(cnf));
			// 4) copy compact cs async.
			CHECK(cudaMemcpyAsync(hcnf->csData(), cumm.csMem(), hcnf->size() * sizeof(S_REF), cudaMemcpyDeviceToHost, s1));
		}
	}

	inline void	ParaFROST::cleanDynamic() {
		if (vars != NULL) delete vars;
		cumm.freeVars(), cumm.freeCNF(), cumm.freeOT();
		vars = NULL, cnf = NULL, ot = NULL;
	}

}