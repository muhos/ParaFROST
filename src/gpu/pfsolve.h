/***********************************************************************[pfsolve.h]
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

#ifndef __SOLVE_
#define __SOLVE_

#include "pfsimplify.cuh"
#include "pfmodel.h"
#include "pfproof.h"
#include "pfsort.h"
#include "pfvector.h"
#include "pfinput.h"
#include "pfheap.h"
#include "pfqueue.h"
#include "pfrestart.h"
#include "pfvmap.h"
#include "pfsolvertypes.h"
#include "pflimit.h"
#include "pfstats.h"
#include "pfoptions.h"
#include "pfrandom.h"
#include "pfcontrol.h"

namespace pFROST {
	/*****************************************************/
	/*  Name:     ParaFROST                              */
	/*  Usage:    global handler for solver/simplifier   */
	/*  Scope:    host only                              */
	/*  Memory:   system memory                          */
	/*****************************************************/
	using namespace SIGmA;

	class ParaFROST {
	protected:
		FORMULA			formula;
		TIMER			timer;
		CMM				cm;
		WT				wt;
		SP				*sp;
		LIMIT			limit;
		LAST			last;
		STATS			stats;
		SLEEP			sleep;
		BCNF			orgs, learnts, reduced;
		VMAP			vmap;
		C_REF			conflict, ignore;
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
		uVec1D			vorg, vhist;
		uVec1D			eligible;
		uVec1D			probes;
		uVec1D			dlevels;
		uVec1D			trail, trailindex;
		uVec1D			analyzed, minimized;
		LBDREST			lbdrest;
		LUBYREST		lubyrest;
		RANDOM			random;
		CNF_ST			cnfstate;
		uint64			bumped;
		size_t			solLineLen;
		string			solLine;
		bool			intr, stable, probed, incremental;
	public:
		OPTION			opts;
		MODEL			model;
		PROOF			proof;
		//============== inline methods ===============
		inline int		calcLBD				(CLAUSE&);
		inline void		bumpClause          (CLAUSE&);
		inline void		countoccurs         (CLAUSE&);
		inline bool		subsumeCheck		(CLAUSE*, uint32&);
		inline CL_ST	subsumeClause		(CLAUSE&, const C_REF&);
		inline void		removeSubsumed		(CLAUSE&, const C_REF&, CLAUSE*);
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
		inline void		analyzeLit			(const uint32&, int&);
		inline uint32	analyzeReason		(const C_REF&, const uint32&);
		inline void		analyzeReason		(const C_REF&, const uint32&, int&);
		inline bool		isBinary			(const C_REF&, uint32&, uint32&);
		inline bool		proplarge			(const uint32&, const bool&);
		inline bool		propbinary			(const uint32&);
		inline void		cancelAssign		(const uint32&);
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
		inline void		varRandPhase		();
		inline bool		verifyMDM			();
		inline bool		verifySeen			();
		inline void		clearMDM			();
		inline void		resetoccurs         ();
		//==============================================
		inline			~ParaFROST			() { }
		inline void		interrupt			() { intr = true; }
		inline void		nointerrupt			() { intr = false; }
		inline void		incDL				() { dlevels.push(trail.size()); }
		inline STATS&	getStats			() { return stats; }
		inline bool		interrupted			() const { return intr; }
		inline CNF_ST	status				() const { return cnfstate; }
		inline int		DL					() const { assert(dlevels.size()); return (int)dlevels.size() - 1; }
		inline int64	maxClauses			() const { return stats.clauses.original + stats.clauses.learnt; }
		inline int64	maxLiterals			() const { return stats.literals.original + stats.literals.learnt; }
		inline bool		useTarget			() const { return stable && opts.targetphase_en; }
		inline bool		vsidsOnly			() const { return stable && opts.vsidsonly_en; }
		inline bool		vsidsEnabled		() const { return stable && opts.vsids_en; }
		inline bool		varsEnough			() const { assert(trail.size() < inf.maxVar); return (inf.maxVar - trail.size()) > last.mdm.unassigned; }
		inline bool		canPreSigmify		() const { return opts.sigma_en; }
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
			uint32 inactive = maxInactive();
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
				|| (opts.boundsearch_en && (stats.conflicts >= opts.conflict_out || stats.decisions.single >= opts.decision_out));
		}
		inline void		rootify				() {
			assert(cnfstate);
			assert(sp->propagated == trail.size());
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
			if (verbose == 4) PFLOG2(2, " Initializing VMFQ Queue with %d variables..", inf.maxVar);
			else PFLOGN2(2, " Initializing VMFQ Queue with %d variables..", inf.maxVar);
			forall_variables(v) { vmtf.init(v), vmtf.update(v, (bumps[v] = ++bumped)); }
			PFLDONE(2, 4);
		}
		inline void		initHeap			() {
			PFLOGN2(2, " Initializing VSIDS Heap with %d variables..", inf.maxVar);
			forall_variables(v) { vsids.insert(v); }
			PFLDONE(2, 5);
		}
		inline void		initVars            () {
			PFLOGN2(2, " Initializing original variables array with %d variables..", inf.maxVar);
			vorg.resize(inf.maxVar + 1);
			vorg[0] = 0;
			forall_variables(v) { vorg[v] = v; }
			PFLDONE(2, 5);
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
				const uint32 v = *a;
				sp->seen[v] = 0;
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
		inline void		clearFrozen         () {
			assert(sp->stacktail - sp->tmp_stack <= inf.maxVar);
			uint32* start = sp->tmp_stack, * tail = sp->stacktail;
			while (start != tail)
				sp->frozen[*start++] = 0;
			assert(verifyFrozen());
		}
		inline bool		verifyFrozen        () {
			for (uint32 v = 0; v <= inf.maxVar; v++) {
				if (sp->frozen[v]) {
					PFLOG0("");
					PFLOGEN("frozen(%d) is not melted", v);
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
		inline uint64	scale				(const uint64& increase)
		{
			const uint64 factor = logn(stats.clauses.original);
			assert(factor >= 1);
			const uint64 dfactor = factor * factor;
			assert(dfactor >= 1);
			uint64 scaled = dfactor * increase;
			assert(increase <= scaled);
			return scaled;
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
			PFLNEWLIT(this, 4, src, lit);
			if (wt.size()) {
				WL& ws = wt[lit];
				if (ws.size()) {
#if _WIN32
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
			PFLNEWLIT(this, 3, NOREF, lit);
		}
		inline void		enqueueImp			(const uint32& lit, const C_REF& src) {
			CHECKLIT(lit);
			const int level = forcedLevel(lit, src);
			enqueue(lit, level, src);
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
			PFLOG2(4, " Variable %d moved to queue front & bumped to %lld", v, bumped);
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
		inline void		markSubstituted     (const uint32& v) {
			CHECKVAR(v);
			assert(!sp->vstate[v].state);
			sp->vstate[v].state = SUBSTITUTED_M;
			assert(inf.unassigned);
			inf.unassigned--;
			inf.maxSubstituted++;
		}
		//==============================================
		bool	keeping				(CLAUSE&);
		CL_ST	rooted				(CLAUSE&);
		void	markLits			(CLAUSE&);
		void	unmarkLits			(CLAUSE&);
		void	sortClause			(CLAUSE&);
		void	histClause          (CLAUSE&);
		void	markSubsume			(CLAUSE&);
		void	bumpShrunken		(CLAUSE&);
		int		removeRooted		(CLAUSE&);
		void	removeClause		(CLAUSE&, const C_REF&);
		void	hyper2Resolve		(CLAUSE&, const uint32&);
		bool	hyper3Resolve		(CLAUSE&, CLAUSE&, const uint32&);
		bool	analyzeVivify		(CLAUSE&, bool&);
		bool	learnVivify			(CLAUSE&, const C_REF&, const int&, const bool&);
		void	shrinkClause		(CLAUSE&, const int&);
		void	shrinkClause		(const C_REF&);
		bool	vivifyClause		(const C_REF&);
		void	markSubsume			(SCLAUSE&);
		void	newClause			(SCLAUSE&);
		void	newHyper2			(const bool&);
		void	newHyper3			(const bool&);
		C_REF	newClause			(const Lits_t&, const bool&);
		void	newClause			(const C_REF&, CLAUSE&, const bool&);
		bool	toClause			(Lits_t&, Lits_t&, char*&);
		void	backtrack			(const int& jmplevel = 0);
		void	recycle				(CMM&);
		void	filter				(BCNF&, CMM&);
		void	filter				(BCNF&);
		void	shrink				(BCNF&);
		void	shrinkTop			(BCNF&);
		void	histBins			(BCNF&);
		void	attachBins			(BCNF&);
		void	attachNonBins		(BCNF&);
		void	attachClauses		(BCNF&);
		void	sortviv				(BCNF&);
		void	schedule2sub		(BCNF&);
		void	schedule2viv		(BCNF&, const bool& tier2, const bool& learnt);
		void	histCNF             (BCNF&, const bool& reset = false);
		void	histCNFFlat         (BCNF&, const bool& reset = false);
		bool	substitute			(BCNF&, uint32*);
		void	attachTernary		(BCNF&, LIT_ST*);
		void	scheduleTernary		(LIT_ST*);
		void	vivifying			(const CL_ST&);
		void	ternarying			(const uint64&, const uint64&);
		void	transiting			(const uint32&, const uint64&, uint64&, uint32&);
		void	ternaryResolve		(const uint32&, const uint64&);
		void	subsumeLearnt		(const C_REF&);
		void	analyzeFailed		(const uint32&);
		uint32	makeAssign			(const uint32&, const bool& tphase = false);
		bool	minimize			(const uint32&, const int& depth = 0);
		void	rebuildWT			(const CL_ST& priorbins = 0);
		bool	canELS				(const bool&);
		void	ELS					(const bool&);
		bool	decompose			();
		void	debinary			();
		bool	canTernary			();
		void	ternary				();
		void	transitive			();
		bool	canVivify			();
		void	vivify				();
		void	sortWT				();
		void	pumpFrozen			();
		void	allocSolver			();
		void	initLimits			();
		void	initSolver			();
		void	killSolver			();
		void	varOrder			();
		bool	shrink				();
		void	markReasons         ();
		void	unmarkReasons       ();
		void	recycleWT			();
		void	recycle				();
		void	reduce				();
		void	reduceLearnts		();
		void	rephase				();
		void	subsume				();
		bool	subsumeAll			();
		void	filterOrg			();
		void	minimize			();
		int		reuse				();
		bool	canRestart			();
		void	vibrate				();
		void	updateModeLimit		();
		void	updateUnstableLimit	();
		void	stableMode			();
		void	unstableMode		();
		void	restart				();
		void	probe				();
		void	FLE					();
		void	scheduleProbes		();
		uint32	nextProbe			();
		int		where				();
		C_REF	backjump			();
		void	analyze				();
		void	finduip				();
		bool	chronoAnalyze		();
		bool	chronoHasRoot		();
		bool	BCPVivify			();
		bool	BCPProbe			();
		bool	BCP					();
		uint32	nextVSIDS			();
		uint32	nextVMFQ			();
		void	MDMFuseMaster		();
		void	MDMFuseSlave		();
		void	MDMInit				();
		void	MDM					();
		bool	canMMD				();
		void	eligibleVSIDS		();
		void	eligibleVMFQ		();
		void	decide				();
		void	report				();
		void	wrapup				();
		bool	parser				();
		void	solve				();
		void	map					(BCNF&);
		void	map					(WL&);
		void	map					(WT&);
		void	map					(const bool& sigmified = false);
				ParaFROST			(const string&);
		//==========================================//
		//                Simplifier                //
		//==========================================//
	protected:
		VARS			*vars;
		TCA				tca;
		cuMM			cumm;
		cuHist			cuhist;
		OT				*ot;
		CNF				*cnf, *hcnf;
		S_REF			dataoff;
		uint32			csoff, mu_inc, ereCls;
		cudaStream_t	*streams;
		std::ofstream	outputFile;
		bool			mapped, compacted;
		int				phase, nForced, sigState, devCount;
	public:
		cuLimit			culimit;
	public:
		//============= inline methods ==============//
		inline bool		alldisabled			() {
			return !opts.phases && !(opts.all_en | opts.ere_en);
		}
		inline bool		reallocFailed		() {
			return (sigState == OTALLOC_FAIL || sigState == CNFALLOC_FAIL);
		}
		inline void		createStreams		() {
			if (streams == NULL) {
				PFLOGN2(2, " Allocating GPU streams..");
				streams = new cudaStream_t[opts.nstreams];
				for (int i = 0; i < opts.nstreams; i++) cudaStreamCreate(streams + i);
				PFLDONE(2, 5);
			}
		}
		inline void		destroyStreams		() {
			if (streams != NULL) {
				for (int i = 0; i < opts.nstreams; i++) cudaStreamDestroy(streams[i]);
				delete[] streams;
			}
		}
		inline bool		verifyLCVE			() {
			for (uint32 v = 0; v < vars->numPVs; v++)
				if (sp->frozen[vars->pVars->at(v)])
					return false;
			return true;
		}
		inline void		logReductions		() {
			int64 varsRemoved = int64(inf.n_del_vars_after) + nForced;
			int64 clsRemoved = int64(inf.nClauses) - inf.n_cls_after;
			int64 litsRemoved = int64(inf.nLiterals) - inf.n_lits_after;
			const char* header = "  %s%-10s  %-10s %-10s %-10s%s";
			PFLOG1(header, CREPORT, " ", "Variables", "Clauses", "Literals", CNORMAL);
			const char* rem = "  %s%-10s: %s%-9lld  %c%-8lld  %c%-8lld%s";
			const char* sur = "  %s%-10s: %s%-9d  %-9d  %-9d%s";
			PFLOG1(rem, CREPORT, "Removed", CREPORTVAL,
				-varsRemoved,
				clsRemoved < 0 ? '+' : '-', abs(clsRemoved),
				litsRemoved < 0 ? '+' : '-', abs(litsRemoved), CNORMAL);
			PFLOG1(sur, CREPORT, "Survived", CREPORTVAL,
				maxActive(),
				inf.n_cls_after,
				inf.n_lits_after, CNORMAL);
		}
		inline void		countMelted			() {
			SIGmA::countMelted(sp->vstate);
		}
		inline void		countFinal			(const bool& host = 0) {
			if (host) countMelted(), countAll(1);
			else SIGmA::countFinal(cnf, vars->gstats, sp->vstate);
			inf.nClauses = inf.n_cls_after;
			inf.nLiterals = inf.n_lits_after;
		}
		inline void		countCls			(const bool& host = 0) {
			if (host) {
				assert(!hcnf->empty());
				inf.n_cls_after = 0;
				for (uint32 i = 0; i < hcnf->size(); i++) {
					SCLAUSE& c = hcnf->clause(i);
					if (c.original() || c.learnt())
						inf.n_cls_after++;
				}
			}
			else SIGmA::countCls(cnf, vars->gstats);
		}
		inline void		countLits			(const bool& host = 0) {
			if (host) {
				assert(!hcnf->empty());
				inf.n_lits_after = 0;
				for (uint32 i = 0; i < hcnf->size(); i++) {
					SCLAUSE& c = hcnf->clause(i);
					if (c.original() || c.learnt())
						inf.n_lits_after += c.size();
				}
			}
			else SIGmA::countLits(cnf, vars->gstats);
		}
		inline void		countAll			(const bool& host = 0) {
			if (host) {
				assert(!hcnf->empty());
				inf.n_cls_after = 0, inf.n_lits_after = 0;
				for (uint32 i = 0; i < hcnf->size(); i++) {
					SCLAUSE& c = hcnf->clause(i);
					if (c.original() || c.learnt())
						inf.n_cls_after++, inf.n_lits_after += c.size();
				}
			}
			else SIGmA::countAll(cnf, vars->gstats);
		}
		inline void		evalReds			(const bool& host = 0) {
			if (host) countMelted(), countAll(1);
			else SIGmA::evalReds(cnf, vars->gstats, sp->vstate);
		}
		inline void		cacheUnits			(const cudaStream_t& stream) {
			sync(stream);
			if ((vars->nUnits = vars->tmpObj.size())) 
				CHECK(cudaMemcpyAsync(vars->cachedUnits, cumm.unitsdPtr(), vars->nUnits * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
			if (sync_always) sync(stream);
		}
		inline void		cacheNumUnits		(const cudaStream_t& stream) {
			CHECK(cudaMemcpyAsync(&vars->tmpObj, vars->units, sizeof(cuVecU), cudaMemcpyDeviceToHost, stream));
			if (sync_always) sync(stream);
		}
		inline void		cacheResolved		(const cudaStream_t& stream) {
			uint32* devStart = *vars->resolved, devSize = vars->resolved->size();
			if (devSize == 0) return;
			uint32 off = model.resolved.size();
			model.resolved.resize(off + devSize);
			uint32* start = model.resolved + off;
			CHECK(cudaMemcpyAsync(start, devStart, devSize * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
			if (sync_always || unified_access) sync(stream);
		}
		inline bool		stop				(const int64 cr, const int64 lr) {
			return (phase == opts.phases) || (!cr && !lr) || (phase > 2 && lr <= opts.lits_min);
		}
		//===========================================//
		void			optSimp				();
		void			varReorder			();
		void			newBeginning		();
		void			sigmifying			();
		void			sigmify				();
		void			awaken				();
		bool			LCVE				();
		void			postVE				();
		void			VE					();
		void			HSE					();
		void			ERE					();
		void			BCE					();
		bool			prop				();
		bool			propFailed			();
		void			masterFree			();
		void			slavesFree			();
		void			createOTHost		(HOT&);
		void			calcOccurs			(const uint32&);
		void			histSimp			(const uint32&);
		void			extract				(CNF*, BCNF&);
		void			cacheCNF			(const cudaStream_t&, const cudaStream_t&);
		void			reflectCNF			(const cudaStream_t&, const cudaStream_t&);
		inline void		depFreeze			(OL&, const uint32&, const uint32&, const uint32&);
		inline bool		propClause			(const LIT_ST*, const uint32&, SCLAUSE&);
		inline bool		enqueueCached		(const cudaStream_t&);
		inline bool		reallocOT			(const cudaStream_t& stream = 0);
		inline bool		reallocCNF			();
		inline void		cleanProped			();
		inline void		cleanManaged		();
		inline void		initSimplifier		();
		inline void		writeback			();
		//==========================================//
		//          Incremental Solving             //
		//==========================================//
	protected:
		Vec<LIT_ST>		ifrozen, ivalue, imarks;
		Vec<VSTATE>		ivstate;
		Vec1D			ilevel;
		Lits_t			assumptions, iconflict;
	public:
						ParaFROST			();
		void			iunassume			();
		void			iallocSpace			();
		uint32			incVariable			();
		void			idecide				();
		void			ianalyze			(const uint32&);
		bool			itoClause			(Lits_t&, Lits_t&);
		void			iassume				(Lits_t&);
		void			isolve				(Lits_t&);
		bool		    ifailed             (const uint32& v);
		void		    ifreeze             (const uint32& v);
		void		    iunfreeze           (const uint32& v);
		bool		    ieliminated         (const uint32& v);
		inline uint32   imap                (const uint32& v) const {
			assert(v && v < NOVAR);
			assert(model.lits.size() > v);
			return model.lits[v];
		}
		inline bool		iassumed            (const uint32& v) const {
			CHECKVAR(v);
			return incremental && ifrozen[v];
		}
		//==========================================//
		//			       Printers                 //
		//==========================================//
		void printStats			(const bool& _p = true, const Byte& _t = ' ', const char* _c = CNORMAL);
		void printVars			(const uint32* arr, const uint32& size, const LIT_ST& type = 'x');
		void printClause		(const Lits_t&);
		void printTrail			(const uint32& off = 0);
		void printCNF			(const BCNF&, const int& off = 0);
		void printOL			(const OL& list, const bool& host = 0);
		void printOL			(const uint32& lit, const bool& host = 0);
		void printWL			(const uint32&, const bool& bin = 0);
		void printWatched		(const uint32&);
		void printBinaries		(const uint32&);
		void printSortedStack	(const int&);
		void printTable			();
		void printWT			();
		void printHeap			();
		void printSource		();
		void printLearnt		();
		
	};
	extern ParaFROST* pfrost;
}

#endif 