/***********************************************************************[solver.hpp]
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

#ifndef __SOLVE_
#define __SOLVE_

#include "sort.hpp"
#include "heap.hpp"
#include "vmap.hpp"
#include "walk.hpp"
#include "model.hpp"
#include "proof.hpp"
#include "limit.hpp"
#include "queue.hpp"
#include "dimacs.hpp"
#include "random.hpp"
#include "memory.hpp"
#include "restart.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "solvertypes.hpp"
#include "thrustalloc.cuh"
#include "variables.cuh"
#include "histogram.cuh"
#include "memory.cuh"
#include "cache.cuh"
#include "proof.cuh"
#include "table.cuh"
#include "cnf.cuh"
#include "key.cuh"
#include "timer.cuh"

namespace ParaFROST {
	/*****************************************************/
	/*  Name:     Solver                              */
	/*  Usage:    global handler for solver/simplifier   */
	/*  Scope:    host only                              */
	/*  Memory:   system memory                          */
	/*****************************************************/
	class Solver {
	protected:
		FORMULA			formula;
		TIMER			timer;
		cuTIMER			cutimer;
		CMM				cm;
		WT				wt;
		SP				*sp;
		LIMIT			limit;
		LAST			last;
		STATS			stats;
		SLEEP			sleep;
		BCNF			orgs;
		BCNF 			learnts;
		BCNF 			reduced;
		VMAP			vmap;
		Lits_t			learntC;
		CLAUSE			subbin;
		QUEUE			vmtf;
		HEAP<VSIDS_CMP>	vsids;
		HEAP<SCORS_CMP>	vschedule;
		Vec<uint64>		bumps;
		Vec<double>		activity;
		Vec<CSIZE>		scheduled;
		Vec<OCCUR>		occurs;
		Vec<DWATCH>		dwatches;
		Vec<WOL>		wot;
		Vec<BOL>		bot;	
		Vec1D			lbdlevels;
		uVec1D			vorg;
		uVec1D 			vhist;
		uVec1D			eligible;
		uVec1D			probes;
		uVec1D			dlevels;
		uVec1D			trail;
		uVec1D			analyzed;
		uVec1D 			minimized;
		LBDREST			lbdrest;
		LUBYREST		lubyrest;
		RANDOM			random;
		WALK			tracker;
		uint64			bumped;
		C_REF			conflict;
		C_REF 			ignore;
		size_t			solLineLen;
		string			solLine;
		CNFState		cnfstate;
		bool			intr;
		bool 			stable;
		bool			probed;
		bool 			incremental;

	public:
		OPTION			opts;
		MODEL			model;
		PROOF			proof;

		//============== inline methods ===============
		inline int		calcLBD				(CLAUSE&);
		inline void		bumpClause			(CLAUSE&);
		inline void		countoccurs			(CLAUSE&);
		inline bool		subsumeCheck		(CLAUSE*, uint32&);
		inline CL_ST	subsumeClause		(CLAUSE&, const C_REF&);
		inline void		removeSubsumed		(CLAUSE&, const C_REF&, CLAUSE*);
		inline void	    strengthenOTF       (CLAUSE&, const C_REF&, const uint32&);
		inline void		strengthen			(CLAUSE&, const uint32&);
		inline LIT_ST	sortClause			(CLAUSE&, const int&, const int&, const bool&);
		inline void		moveClause			(C_REF&, CMM&);
		inline void		moveWatches			(WL&, CMM&);
		inline uint32	minReachable		(WL&, DFS*, const DFS&);
		inline bool		depFreeze			(WL&, const uint32&);
		inline bool		valid				(WL&);
		inline void		recycleWL			(const uint32&);
		inline bool		findBinary			(uint32, uint32);
		inline bool		findTernary			(uint32, uint32, uint32);
		inline void		analyzeLit			(const uint32&, int&, int&);
		inline uint32	analyzeReason		(const C_REF&, const uint32&);
		inline bool		analyzeReason		(const C_REF&, const uint32&, int&);
		inline bool		isBinary			(const C_REF&, uint32&, uint32&);
		inline uint32	propAutarkClause	(const bool&, const C_REF&, CLAUSE&, const LIT_ST*, LIT_ST*);
		inline bool		proplarge			(const uint32&, const bool&);
		inline bool		propbinary			(const uint32&);
		inline void		cancelAssign		(const uint32&);
		inline void		cancelAutark		(const bool&, const uint32&, LIT_ST*);
		inline void		pumpFrozenHeap		(const uint32&);
		inline void		pumpFrozenQue		(const uint32&);
		inline void		bumpVariable		(const uint32&);
		inline void		bumpReason			(const uint32&);
		inline void		bumpReasons			(const uint32&);
		inline void		bumpReasons			();
		inline void		bumpVariables		();
		inline void		savePhases			();
		inline void		varOrgPhase			();
		inline void		varInvPhase			();
		inline void		varFlipPhase		();
		inline void		varBestPhase		();
		inline void		varWalkPhase		();
		inline void		varRandPhase		();
		inline bool		verifyMDM			();
		inline bool		verifySeen			();
		inline void		clearMDM			();
		inline int		calcLBD				();
		inline void		resetoccurs			();
		inline void		interrupt			() { intr = true; }
		inline void		nointerrupt			() { intr = false; }
		inline void		incDL				() { dlevels.push(trail.size()); }
		inline bool		interrupted			() const { return intr; }
		inline CNFState	status				() const { return cnfstate; }
		inline int		DL					() const { assert(dlevels.size()); return (int)dlevels.size() - 1; }
		inline int64	maxClauses			() const { return stats.clauses.original + stats.clauses.learnt; }
		inline int64	maxLiterals			() const { return stats.literals.original + stats.literals.learnt; }
		inline bool		useTarget			() const { return (stable && opts.targetphase_en) || opts.targetonly_en; }
		inline bool		vsidsOnly			() const { return (stable && opts.vsidsonly_en); }
		inline bool		vsidsEnabled		() const { return (stable && opts.vsids_en); }
		inline bool		varsEnough			() const { assert(trail.size() < inf.maxVar); return (inf.maxVar - trail.size()) > last.mdm.unassigned; }
		inline bool		canPreSimplify		() const { return opts.sigma_en && stats.clauses.original; }
		inline bool		canRephase			() const { return opts.rephase_en && stats.conflicts > limit.rephase; }
		inline bool		canReduce			() const { return opts.reduce_en && stats.clauses.learnt && stats.conflicts >= limit.reduce; }
		inline bool		canCollect			() const { return cm.garbage() > (cm.size() * opts.gc_perc); }
		inline bool		canProbe			() const {
			if (!opts.probe_en) return false;
			if (last.probe.reduces > stats.reduces) return false;
			return limit.probe <= stats.conflicts;
		}
		inline bool		canSubsume			() const {
			if (!opts.subsume_en) return false;
			if (stats.conflicts < limit.subsume) return false;
			return true;
		}
		inline bool		canMap				() const {
			if (DL()) return false;
			const uint32 inactive = maxInactive();
			assert(inactive <= inf.maxVar);
			return inactive > (opts.map_perc * inf.maxVar);
		}
		inline bool		canSigmify			() const {
			if (!opts.sigma_live_en) return false;
			if (!stats.clauses.original) return false;
			if (last.sigma.reduces > stats.reduces) return false;
			if (limit.sigma > stats.conflicts) return false;
			if (sp->simplified >= opts.sigma_min) return true;
			return ((stats.shrunken - last.shrink.removed) > (opts.sigma_min << 4));
		}
		inline bool		runningout			() const { 
			return interrupted()
				|| (termCallback && termCallback(termCallbackState))
				|| (opts.boundsearch_en && (stats.conflicts >= opts.conflict_out || stats.decisions.single >= opts.decision_out));
		}
		inline bool		isPropagated		() const { return sp->propagated == trail.size(); }
		inline void		rootify				() {
			assert(IS_UNSOLVED(cnfstate));
			assert(isPropagated());
			backtrack();
			if (BCP()) learnEmpty();
		}
		inline bool		retrail				() {
			assert(!DL());
			sp->propagated = 0;
			if (BCP()) {
				learnEmpty();
				return true;
			}
			return false;
		}
		inline void		markLearnt			() {
			forall_clause(learntC, k) {
				const uint32 lit = *k;
				sp->marks[ABS(lit)] = SIGN(lit);
			}
		}
		inline void		unmarkLearnt		() {
			forall_clause(learntC, k) {
				assert(l2marker(*k) >= 0);
				sp->marks[ABS(*k)] = UNDEFINED;
			}
		}
		inline void		initQueue			() {
			if (verbose == 4) LOG2(2, " Initializing VMFQ Queue with %d variables..", inf.maxVar);
			else LOGN2(2, " Initializing VMFQ Queue with %d variables..", inf.maxVar);
			forall_variables(v) { vmtf.init(v), vmtf.update(v, (bumps[v] = ++bumped)); }
			LOGDONE(2, 4);
		}
		inline void		initHeap			() {
			LOGN2(2, " Initializing VSIDS Heap with %d variables..", inf.maxVar);
			forall_variables(v) { vsids.insert(v); }
			LOGDONE(2, 5);
		}
		inline void		initVars			() {
			LOGN2(2, " Initializing original variables array with %d variables..", inf.maxVar);
			vorg.resize(inf.maxVar + 1);
			vorg[0] = 0;
			forall_variables(v) { vorg[v] = v; }
			LOGDONE(2, 5);
		}
		inline void		updateQueue			()
		{
			const uint32 last = vmtf.last();
			vmtf.update(last, bumps[last]);
		}
		inline void		updateHeap			()
		{
			forall_variables(v) {
				if (!sp->vstate[v].state && !vsids.has(v))
					vsids.insert(v);
			}
		}
		inline void		learnEmpty			() {
			if (opts.proof_en) proof.addEmpty();
			cnfstate = UNSAT;
		}
		inline void		clearLevels			() {
			forall_vector(int, lbdlevels, dl) {
				sp->vstate[*dl].dlcount = 0;
			}
			lbdlevels.clear();
		}
		inline void		clearMinimized		() {
			forall_vector(uint32, minimized, v) {
				sp->seen[*v] = 0;
			}
			minimized.clear();
		}
		inline void		clearAnalyzed		() {
			forall_vector(uint32, analyzed, a) {
				sp->seen[*a] = 0;
			}
			analyzed.clear();
		}
		inline void		clearVivify			() {
			forall_vector(uint32, analyzed, a) {
				CHECKLIT(*a);
				sp->seen[ABS(*a)] = 0;
			}
			analyzed.clear();
			learntC.clear();
		}
		inline void		clearFrozen			() {
			assert(sp->stacktail - sp->tmpstack <= inf.maxVar);
			uint32* start = sp->tmpstack, *tail = sp->stacktail;
			while (start != tail)
				sp->frozen[*start++] = 0;
			assert(verifyFrozen());
		}
		inline bool		verifyFrozen		() {
			for (uint32 v = 0; v <= inf.maxVar; v++) {
				if (sp->frozen[v]) {
					LOG0("");
					LOGERRORN("frozen(%d) is not melted", v);
					printWatched(v);
					return false;
				}
			}
			return true;
		}
		inline void		attachDelayed		() {
			forall_dwatches(dwatches, d) {
				wt[FLIP(d->lit)].push(WATCH(d->ref, d->size, d->imp));
			}
			dwatches.clear();
		}
		inline uint64	cacheLines			(const size_t& len, const size_t& unit) {
			return (len * unit) >> 6; // 64-byte cache line is assumed
		}
		inline void		markLit				(const uint32& lit) { CHECKLIT(lit); sp->marks[ABS(lit)] = SIGN(lit); }
		inline void		unmarkLit			(const uint32& lit) { CHECKLIT(lit); sp->marks[ABS(lit)] = UNDEFINED; }
		inline void		markProbe			(const uint32& lit) { CHECKLIT(lit); sp->vstate[ABS(lit)].probe = 1; }
		inline void		unmarkProbe			(const uint32& lit) { CHECKLIT(lit); sp->vstate[ABS(lit)].probe = 0; }
		inline void		markSubsume			(const uint32& lit) { CHECKLIT(lit); sp->vstate[ABS(lit)].subsume = 1; }
		inline void		unmarkSubsume		(const uint32& lit) { CHECKLIT(lit); sp->vstate[ABS(lit)].subsume = 0; }
		inline bool		markedSubsume		(const uint32& lit) const { CHECKLIT(lit); return sp->vstate[ABS(lit)].subsume; }
		inline bool		markedProbe			(const uint32& lit) const { CHECKLIT(lit); return sp->vstate[ABS(lit)].probe; }
		inline bool		active				(const uint32& lit) const { CHECKLIT(lit); return !sp->vstate[ABS(lit)].state; }
		inline bool		inactive			(const uint32& lit) const { CHECKLIT(lit); return sp->vstate[ABS(lit)].state; }
		inline bool		subsumed			(const uint32& lit) const { CHECKLIT(lit); return sp->marks[ABS(lit)] == SIGN(lit); }
		inline bool		notsubsumed			(const uint32& lit) const { CHECKLIT(lit); return sp->marks[ABS(lit)] ^ SIGN(lit); }
		inline bool		selfsubsumed		(const uint32& lit) const { CHECKLIT(lit); return sp->marks[ABS(lit)] == !SIGN(lit); }
		inline LIT_ST	l2marker			(const uint32& lit) const { CHECKLIT(lit); return sp->marks[ABS(lit)]; }
		inline int		l2dl				(const uint32& lit) const { CHECKLIT(lit); return sp->level[ABS(lit)]; }
		inline C_REF	l2r					(const uint32& lit) const { CHECKLIT(lit); return sp->source[ABS(lit)]; }
		inline LIT_ST	unassigned			(const uint32& lit) const { CHECKLIT(lit); return UNASSIGNED(sp->value[lit]); }
		inline LIT_ST	isFalse				(const uint32& lit) const { CHECKLIT(lit); return !sp->value[lit]; }
		inline LIT_ST	isTrue				(const uint32& lit) const { 
			CHECKLIT(lit);
			LIT_ST val = sp->value[lit];
			return val && !UNASSIGNED(val);
		}
		inline void		attachClause		(const C_REF& ref, const CLAUSE& c) {
			assert(ref < NOREF);
			const int size = c.size();
			assert(size > 1);
			for (int i = 0; i < size; i++) {
				CHECKLIT(c[i]);
				wot[c[i]].push(ref);
			}
		}
		inline void		attachBinary		(const C_REF& ref, const CLAUSE& c) {
			assert(ref < NOREF);
			assert(c.binary());
			wot[c[0]].push(ref);
			wot[c[1]].push(ref);
		}
		inline void		attachWatch			(const C_REF& ref, const CLAUSE& c) {
			assert(ref < NOREF);
			CHECKLIT(c[0]);
			CHECKLIT(c[1]);
			const int size = c.size();
			assert(size > 1);
			wt[FLIP(c[0])].push(WATCH(ref, size, c[1]));
			wt[FLIP(c[1])].push(WATCH(ref, size, c[0]));
		}
		inline void		attachWatch			(const uint32& lit, const uint32& imp, const C_REF& ref, const int& size) {
			CHECKLIT(lit);
			CHECKLIT(imp);
			assert(lit != imp);
			assert(ref < NOREF);
			assert(size > 1);
			wt[FLIP(lit)].push(WATCH(ref, size, imp));
		}
		inline void		delayWatch			(const uint32& lit, const uint32& imp, const C_REF& ref, const int& size) {
			CHECKLIT(lit);
			CHECKLIT(imp);
			assert(lit != imp);
			assert(ref < NOREF);
			assert(size > 1);
			dwatches.push(DWATCH(lit, imp, ref, size));
		}
		inline void		detachWatch			(const uint32& lit, const C_REF& ref) {
			CHECKLIT(lit);
			assert(ref < NOREF);
			WL& ws = wt[lit];
			if (ws.empty()) return;
			WATCH *j = ws;
			forall_watches(ws, i) {
				const WATCH w = *i;
				if (w.ref != ref) *j++ = w;
			}
			assert(j + 1 == ws.end());
			ws.resize(int(j - ws));
		}
		inline void		enqueueDecision		(const uint32& lit) {
			CHECKLIT(lit);
			assert(inf.unassigned);
			incDL();
			enqueue(lit, DL());
		}
		inline void		enqueue				(const uint32& lit, const int& level = 0, const C_REF src = NOREF) {
			CHECKLIT(lit);
			assert(unassigned(lit));
			const uint32 v = ABS(lit);
			if (!probed) sp->psaved[v] = SIGN(lit);
			sp->source[v] = src;
			sp->level[v] = level;
			sp->value[lit] = 1;
			sp->value[FLIP(lit)] = 0;
			trail.push(lit);
			assert(inf.unassigned);
			inf.unassigned--;
			if (!level) learnUnit(lit, v);
#ifdef LOGGING
			LOGNEWLIT(this, 4, src, lit);
#endif
			if (wt.size()) {
				WL& ws = wt[lit];
				if (ws.size()) {
#ifdef _WIN32
					PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, &ws[0]);
#else
					__builtin_prefetch(&ws[0], 0, 1);
#endif
				}
			}
		}
		inline void		enqueueUnit			(const uint32& lit) {
			CHECKLIT(lit);
			assert(active(lit));
			assert(!DL());
			const uint32 v = ABS(lit);
			learnUnit(lit, v);
			sp->level[v] = 0;
			sp->value[lit] = 1;
			sp->value[FLIP(lit)] = 0;
			trail.push(lit);
			assert(inf.unassigned);
			inf.unassigned--;
#ifdef LOGGING
			LOGNEWLIT(this, 3, NOREF, lit);
#endif
		}
		inline void		enqueueImp			(const uint32& lit, const C_REF& src) {
			CHECKLIT(lit);
			enqueue(lit, forcedLevel(lit, src), src);
		}
		inline int		forcedLevel			(const uint32& lit, const C_REF& src) {
			CHECKLIT(lit);
			assert(REASON(src));
			int fl = 0;
			CLAUSE& c = cm[src];
			forall_clause(c, k) {
				const uint32 other = *k;
				if (NEQUAL(other, lit)) {
					assert(!unassigned(other));
					const int otherlevel = l2dl(other);
					if (otherlevel > fl) fl = otherlevel;
				}
			}
			return fl;
		}
		inline void		learnUnit			(const uint32& lit, const uint32& v) {
			if (opts.proof_en) proof.addUnit(lit);
			assert(ABS(lit) == v);
			markFrozen(v);
		}
		inline void		varBumpHeap			(const uint32& v, const double& norm_act) {
			CHECKVAR(v);
			assert(norm_act);
			assert(activity[v] == 0.0);
			activity[v] = norm_act;
			vsids.bump(v);
		}
		inline void		varBumpHeap			(const uint32& v) {
			CHECKVAR(v);
			assert(last.vsids.inc);
			const double min = 1e-150, max = 1e150;
			assert(activity[v] <= max);
			double newAct = activity[v] + last.vsids.inc;
			if (newAct > max) {
				for (uint32 i = 1; i <= inf.maxVar; i++) activity[i] *= min;
				last.vsids.scale(min);
				newAct = activity[v] + last.vsids.inc;
			}
			assert(newAct <= max);
			activity[v] = newAct;
			vsids.bump(v);
		}
		inline void		varBumpQueue		(const uint32& v) {
			CHECKVAR(v);
			if (!vmtf.next(v)) return;
			vmtf.toFront(v);
			assert(bumped && bumped < UINT64_MAX);
			bumps[v] = ++bumped;
			if (UNASSIGNED(sp->value[V2L(v)])) vmtf.update(v, bumps[v]);
		}
		inline void		varBumpQueueNU		(const uint32& v) {
			CHECKVAR(v);
			if (!vmtf.next(v)) return;
			assert(bumped && bumped < UINT64_MAX);
			vmtf.toFront(v);
			bumps[v] = ++bumped;
#ifdef LOGGING
			LOG2(4, " Variable %d moved to queue front & bumped to %lld", v, bumped);
#endif
		}
		inline uint32	prescore			(const uint32& v) {
			CHECKVAR(v);
			return occurs[v].ps * occurs[v].ns;
		}
		inline uint32	livescore			(const uint32& v) {
			CHECKVAR(v);
			const uint32 p = V2L(v), n = NEG(p);
			const int ps = wot[p].size(), ns = wot[n].size();
			return uint32(ps) * ns;
		}
		inline void		markFrozen			(const uint32& v) {
			CHECKVAR(v);
			sp->vstate[v].state = FROZEN_M;
			inf.maxFrozen++;
		}
		inline void		markEliminated		(const uint32& v) {
			CHECKVAR(v);
			assert(!sp->vstate[v].state);
			sp->vstate[v].state = MELTED_M;
			assert(inf.unassigned);
			inf.unassigned--;
		}
		inline void		markSubstituted		(const uint32& v) {
			CHECKVAR(v);
			assert(!sp->vstate[v].state);
			sp->vstate[v].state = SUBSTITUTED_M;
			assert(inf.unassigned);
			inf.unassigned--;
			inf.maxSubstituted++;
		}
		inline LIT_ST*	getValues			() { return sp->value; }
		inline int*		getLevels			() { return sp->level; }
		//==============================================
		bool			keeping				(CLAUSE&);
		CL_ST			rooted				(CLAUSE&);
		CL_ST			rootedTop			(CLAUSE&);
		void			markLits			(CLAUSE&);
		void			unmarkLits			(CLAUSE&);
		void			sortClause			(CLAUSE&);
		void			histClause			(CLAUSE&);
		void			markSubsume			(CLAUSE&);
		void			bumpShrunken		(CLAUSE&);
		int				removeRooted		(CLAUSE&);
		void			removeClause		(CLAUSE&, const C_REF&);
		uint32			hyper2Resolve		(CLAUSE&, const uint32&);
		bool			hyper3Resolve		(CLAUSE&, CLAUSE&, const uint32&);
		bool			analyzeVivify		(CLAUSE&, bool&);
		bool			learnVivify			(CLAUSE&, const C_REF&, const int&, const bool&);
		void			shrinkClause		(CLAUSE&, const int&);
		void			shrinkClause		(const C_REF&);
		bool			vivifyClause		(const C_REF&);
		void			markSubsume			(SCLAUSE&);
		void			newClause			(SCLAUSE&);
		C_REF			newClause			(const Lits_t&, const bool&);
		void			newClause			(const C_REF&, CLAUSE&, const bool&);
		bool			toClause			(Lits_t&, Lits_t&, int&);
		bool			toClause			(Lits_t&, Lits_t&, char*&);
		void			backtrack			(const int& jmplevel = 0);
		void			recycle				(CMM&);
		void			filter				(BCNF&, CMM&);
		void			filter				(BCNF&);
		void			shrink				(BCNF&);
		void			shrinkTop			(BCNF&);
		void			histBins			(BCNF&);
		void			sortviv				(BCNF&);
		void			schedule2sub		(BCNF&);
		void			schedule2viv		(BCNF&, const bool& tier2, const bool& learnt);
		void			histCNF				(BCNF&, const bool& reset = false);
		void			histCNFFlat			(BCNF&, const bool& reset = false);
		void			attachBins			(BCNF&, const bool& hasElim = false);
		void			attachNonBins		(BCNF&, const bool& hasElim = false);
		void			attachClauses		(BCNF&, const bool& hasElim = false);
		bool			substitute			(BCNF&, uint32*);
		void			attachTernary		(BCNF&, LIT_ST*);
		void			scheduleTernary		(LIT_ST*);
		uint32			autarkReasoning		(LIT_ST*);
		uint32			useAutarky			(LIT_ST*);
		uint32			propAutarky			(const LIT_ST*, LIT_ST*);
		void			vivifying			(const CL_ST&);
		void			ternarying			(const uint64&, const uint64&);
		void			transiting			(const uint32&, const uint64&, uint64&, uint32&);
		void			ternaryResolve		(const uint32&, const uint64&);
		void			subsumeLearnt		(const C_REF&);
		void			analyzeFailed		(const uint32&);
		uint32			makeAssign			(const uint32&, const bool& tphase = false);
		bool			minimize			(const uint32&, const int& depth = 0);
		void			rebuildWT			(const CL_ST& priorbins = 0);
		void			binarizeWT			(const bool& keeplearnts);
		void			detachClauses		(const bool& keepbinaries);
		bool			canELS				(const bool&);
		void			ELS					(const bool&);
		void			shrinkTop			(const bool&);
		void			newHyper3			(const bool&);
		void			newHyper2			();
		bool			shrink				();
		bool			decompose			();
		void			debinary			();
		bool			canTernary			();
		void			ternary				();
		void			transitive			();
		bool			canVivify			();
		void			vivify				();
		void			sortWT				();
		void			pumpFrozen			();
		void			allocSolver			();
		void			initLimits			();
		void			initSolver			();
		void			killSolver			();
		void			varOrder			();
		void			markReasons         ();
		void			unmarkReasons       ();
		void			recycleWT			();
		void			recycle				();
		void			reduce				();
		void			reduceLearnts		();
		void			rephase				();
		void			autarky				();
		void			filterAutarky		();
		void			subsume				();
		bool			subsumeAll			();
		void			filterOrg			();
		void			minimize			();
		void			minimizebin			();
		int				reuse				();
		bool			canRestart			();
		void			vibrate				();
		void			updateModeLimit		();
		void			updateUnstableLimit	();
		void			stableMode			();
		void			unstableMode		();
		void			restart				();
		void			probe				();
		void			FLE					();
		void			scheduleProbes		();
		uint32			nextProbe			();
		int				where				();
		C_REF			backjump			();
		void			analyze				();
		bool			finduip				();
		bool			chronoAnalyze		();
		bool			chronoHasRoot		();
		bool			BCPVivify			();
		bool			BCPProbe			();
		bool			BCP					();
		uint32			nextVSIDS			();
		uint32			nextVMFQ			();
		void			MDMInit				();
		void			MDM					();
		bool			canMMD				();
		void			eligibleVSIDS		();
		void			eligibleVMFQ		();
		void			decide				();
		void			report				();
		void			wrapup				();
		void			solve				();
		bool			parse				();
		bool			parse				(const string&);
		void			map					(BCNF&);
		void			map					(WL&);
		void			map					(WT&);
		void			map					(const bool& sigmified = false);
		void 			initialize			(const bool& banner);
		explicit 		Solver				(const std::string& path);  				// incremental = false, takes formula
						Solver				();                       					// incremental = false
						~Solver				();
		//==========================================//
		//                Simplifier                //
		//==========================================//
	protected:
		VARS			*vars;
		OT				*ot;
		CNF				*cnf, *hcnf;
		CACHER			cacher;
		TCA				tca;
		cuMM			cumm;
		cuHist			cuhist;
		cuPROOF			cuproof;
		S_REF			dataoff;
		uint32			csoff, multiplier, ereCls;
		cudaStream_t	*streams;
		COMPACT_CMP		compact_cmp;
		OLIST_CMP		olist_cmp;
		bool			mapped, compacted, flattened;
		int				phase, nForced, simpstate, devCount;
	
	public:

		// Getters
		inline VARS*	getVars				()		{ return vars; }
		inline OT*		getDeviceTable		()		{ return ot; }
		inline CNF*		getDeviceCNF		()		{ return cnf; }
		inline CNF*		getHostCNF			()		{ return hcnf; }
		inline CACHER&	getCacher			()		{ return cacher; }
		inline TCA&		getThrustAllocator	()		{ return tca; }
		inline cuMM&	getDeviceAllocator	()		{ return cumm; }
		inline cuHist&	getDeviceHistogram	()		{ return cuhist; }
		inline cuPROOF&	getDeviceProver		()		{ return cuproof; }
		inline int		getDeviceCount		()		{ return devCount; }
		inline int		getSimplifierState	()		{ return simpstate; }
		inline bool		isCompacted			()		{ return compacted; }
		// Helpers
		inline void		enqueueDevUnit		(const uint32& lit) {
			CHECKLIT(lit);
			assert(active(lit));
			assert(!DL());
			const uint32 v = ABS(lit);
			markFrozen(v);
			sp->level[v] = 0;
			sp->value[lit] = 1;
			sp->value[FLIP(lit)] = 0;
			trail.push(lit);
			assert(inf.unassigned);
			inf.unassigned--;
#ifdef LOGGING
			LOGNEWLIT(this, 3, NOREF, lit);
#endif
		}
		inline bool		alldisabled         () {
			return !opts.phases && !(opts.all_en | opts.ere_en);
		}
		inline bool		reallocFailed       () {
			return (simpstate == OTALLOC_FAIL || simpstate == CNFALLOC_FAIL);
		}
		inline void		createStreams       () {
			if (streams == NULL) {
				LOGN2(2, " Allocating GPU streams..");
				streams = new cudaStream_t[opts.nstreams];
				for (int i = 0; i < opts.nstreams; i++) cudaStreamCreate(streams + i);
				LOGDONE(2, 5);
			}
		}
		inline void		destroyStreams      () {
			if (streams != NULL) {
				for (int i = 0; i < opts.nstreams; i++) cudaStreamDestroy(streams[i]);
				delete[] streams;
			}
		}
		inline void		updateNumPVs		();
		inline bool		verifyLCVE			();
		void			logReductions		();
		uint32			updateNumElims		();
		void			countCls			(const bool& host = 0);
		void			countLits			(const bool& host = 0);
		void			countAll			(const bool& host = 0);
		inline bool		stop                (const int64 cr, const int64 lr) {
			return (phase == opts.phases)
				|| (simpstate == CNFALLOC_FAIL)
				|| (!cr && !lr) 
				|| (phase > 2 && lr <= opts.phase_lits_min);
		}
		//===========================================//
		void			optSimp             ();
		void			freeSimp			();
		void			varReorder          ();
		void			newBeginning        ();
		void			awaken              ();
		bool			LCVE                ();
		void			sortOT				();
		void			postVE              ();
		void			VE                  ();
		void			SUB                 ();
		void			ERE                 ();
		void			BCE                 ();
		bool			prop                ();
		bool			propFailed          ();
		void			writeBackCNF		();
		bool			reallocCNF			();
		bool			reallocCNF			(const bool& realloc);
		void			simplifying         (const bool& keep_gpu_mem = false);
		void			simplify            (const bool& keep_gpu_mem = false);
		void			extractCNF			(CNF*, BCNF&);
		uint32*			flattenCNF          (const uint32&);
		void			reflectCNF			(const cudaStream_t&, const cudaStream_t&);
		void			cacheCNF			(const cudaStream_t&, const cudaStream_t&);
		void			histSimp            (const uint32&);
		void 			resetOTAsync		(const cudaStream_t& stream = 0);
		void 			reduceOTAsync		(const bool& print = false, const cudaStream_t& stream = 0);
		void			createOTAsync		(const bool& print = false, const cudaStream_t& stream = 0);
		void			createOTHost		(HOT&);
		void			cacheUnits          (const cudaStream_t&);
		void			cacheNumUnits       (const cudaStream_t&);
		void			cacheResolved       (const cudaStream_t&);
		void			cacheEliminated     (const cudaStream_t&);
		void			markEliminated		(const cudaStream_t&);
		bool			reallocOT			(const cudaStream_t& stream = 0);
		inline bool		depFreeze           (OL& ol, LIT_ST* frozen, uint32*& tail, const uint32& cand, const uint32& pmax, const uint32& nmax);
		inline bool		propClause          (const LIT_ST*, const uint32&, SCLAUSE&);
		inline bool		enqueueCached       (const cudaStream_t&);
		inline void		mapFrozen			();
		inline void		cleanProped         ();
		inline void		cleanManaged        ();
		inline void		initSimplifier      ();
			   void		initSharedMem		();
		//==========================================//
		//             Local search                 //
		//==========================================//
		inline bool		popUnsat			(const uint32&, const uint32&, Vec<CINFO>&);
		inline void		saveTrail			(const LIT_ST*, const bool&);
		inline void		saveAll				(const LIT_ST*);
		inline uint32	breakValue			(const uint32&);
		inline void		makeClauses			(const uint32&);
		inline void		breakClauses		(const uint32&);
		inline void		walkassign			();
		uint32			promoteLit			();
		uint32			ipromoteLit			();
		bool			walkschedule		();
		void			walkinit			();
		void			walkstep			();
		void			walkstop			();
		void			walking				();
		void			walk				();
		//==========================================//
		//          Incremental Solving             //
		//==========================================//
	protected:
		Vec<LIT_ST>		ifrozen, ivalue, iphase, imarks;
		Vec<VSTATE>		ivstate;
		Vec<C_REF>		isource;
		Vec1D			ilevel;
		Lits_t			assumptions, iconflict;

	public:
		void*			termCallbackState;
        void*			learnCallbackState;
        int*			learnCallbackBuffer;
        int				learnCallbackLimit;
		explicit 		Solver				(const bool& inc);   						// sets incremental
		explicit 		Solver				(const bool& inc, const std::string& path); // sets incremental, takes formula
		void			iunassume           ();
		void			iallocspace         ();
		uint32			iadd                ();
		void			idecide             ();
		void			ianalyze            (const uint32&);
		bool			iclause             (Lits_t&, Lits_t&);
		void			iassume             (Lits_t&);
		void			isolve              (Lits_t&);
		bool		    ifailed             (const uint32& v);
		void		    ifreeze             (const uint32& v);
		void		    iunfreeze           (const uint32& v);
		bool		    ieliminated         (const uint32& v);
		void 			resetextended   	() { if (model.extended) model.extended = false; }
		inline uint32   imap                (const uint32& v) const {
			assert(v && v < NOVAR);
			assert(model.lits.size() > v);
			return model.lits[v];
		}
		inline bool		iassumed            (const uint32& v) const {
			CHECKVAR(v);
			return incremental && ifrozen[v];
		}
		int				(*termCallback)		(void* state);
        void			(*learnCallback)	(void* state, int* clause);
        void			setTermCallback		(void* state, int (*termCallback)(void*)) {
            this->termCallbackState = state;
            this->termCallback = termCallback;
        }
        void			setLearnCallback	(void* state, int maxlength, void (*learn)(void* state, int* clause)) {
            this->learnCallbackState = state;
            this->learnCallbackLimit = maxlength;
            this->learnCallback = learn;
            pfralloc<int>(this->learnCallbackBuffer, (maxlength + 1) * sizeof(int));
        }


		//==========================================//
		//			       Printers                 //
		//==========================================//
		void			printStats			(const bool& _p = true, const Byte& _t = ' ', const char* _c = CNORMAL);
		void			printVars			(const uint32* arr, const uint32& size, const LIT_ST& type = 'x');
		void			printClause			(const Lits_t&);
		void			printTrail			(const uint32& off = 0);
		void			printCNF			(const BCNF&, const int& off = 0);
		void			printOL				(const OL& list, const bool& host = 0);
		void			printOL				(const uint32& lit, const bool& host = 0);
		void			printWL				(const uint32&, const bool& bin = 0);
		void			printWL				(const WL&, const bool& bin = 0);
		void			printWatched		(const uint32& v);
		void			printBinaries		(const uint32& v);
		void			printSortedStack	(const int&);
		void			printTable			();
		void			printWT				();
		void			printHeap			();
		void			printSource			();
		void			printLearnt			();
		
	};
}

#endif 