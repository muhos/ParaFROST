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

#include "pfsimp.cuh"
#include "pfmodel.h"
#include "pfsort.h"
#include "pfvec.h"
#include "pfargs.h"
#include "pfheap.h"
#include "pfqueue.h"
#include "pfrestart.h"
#include "pfvmap.h"
#include "pfsolvertypes.h"
#include "pfoptions.h"
#include "pfcontrol.h"

namespace pFROST {

	using namespace SIGmA;
	/*****************************************************/
	/*  Name:     ParaFROST                              */
	/*  Usage:    global handler for solver/simplifier   */
	/*  Scope:    host only                              */
	/*  Memory:   system memory                          */
	/*****************************************************/
	class ParaFROST {
	protected:
		string			path;
		TIMER			timer;
		CMM				cm;
		SP				*sp;
		LEARN			lrn;
		STATS			stats;
		WT				wt;
		BCNF			orgs, learnts, reduced;
		VMAP			vmap;
		MODEL			model;
		C_REF			conflict;
		Lits_t			learntC;
		CLAUSE			subbin;
		Vec<CSIZE>		scheduled;
		Vec<OCCUR>		occurs;
		Vec<int64>		bumps;
		Vec<double>		activity;
		Vec<WOL>		wot;
		Vec<BOL>		bot;
		uVec1D			trail, dlevels, eligible, analyzed, minimized;
		uVec1D          subhist, vorg;
		HEAP<HEAP_CMP>	vsids;
		QUEUE			vmfq;
		LBDREST			lbdrest;
		LUBYREST		lubyrest;
		int64			nConflicts, subleftovers;
		uint32			starts;
		CNF_ST			cnfstate;
		size_t			solLineLen;
		string			solLine;
		std::ofstream	proofFile;
		bool			intr;
	public:
		OPTION			opts;
		//============== inline methods ===============
		inline void		strengthen			(CLAUSE&, const uint32&);
		inline int		removeRooted		(CLAUSE&);
		inline void		removeSubsumed		(CLAUSE&, const C_REF&, CLAUSE*, const C_REF&);
		inline bool		depFreeze			(WL&, const uint32&);
		inline bool		valid				(WL&);
		inline void		recycleWL			(WL&, CMM&);
		inline void		recycleWL			(const uint32&);
		inline void		reduceWeight		(double&);
		inline void		savePhases			(const int&);
		inline CL_ST	subsumeClause		(const C_REF&, CLAUSE&, BCNF&);
		inline bool		subsumeCheck		(CLAUSE*, uint32&);
		inline void		bumpClause			(CLAUSE&);
		inline void		moveClause			(C_REF&, CMM&);
		inline void		analyzeLit			(const uint32&, int&);
		inline void		analyzeReason		(const C_REF&, const uint32&, int&);
		inline void		cancelAssign		(const uint32&);
		inline int		calcLBD				(const CLAUSE&);
		inline void		pumpFrozenHeap		(const uint32&);
		inline void		pumpFrozenQue		(const uint32&);
		inline void		bumpVariable		(const uint32&);
		inline bool		bumpReason			(const uint32&);
		inline void		bumpReasons			(const uint32&, const int&);
		inline void		bumpReasons			();
		inline void		bumpVariables		();
		inline int		whereToJump			();
		inline int		calcLBD				();
		inline void		clearAnalyzed		();
		inline void		clearMinimized		();
		inline void		varOrgPhase			();
		inline void		varInvPhase			();
		inline void		varFlipPhase		();
		inline void		varBestPhase		();
		inline void		varRandPhase		();
		inline bool		verifyMDM			();
		inline bool		verifySeen			();
		//==============================================
		inline			~ParaFROST			() { }
		inline void		interrupt			() { intr = true; }
		inline void		incDL				() { dlevels.push(trail.size()); }
		inline void		decayVarAct			() { lrn.var_inc *= (1.0 / lrn.var_decay); }
		inline STATS&	getStats			() { return stats; }
		inline int64	maxClauses			() const { return int64(orgs.size()) + learnts.size(); }
		inline uint32	DL					() const { return dlevels.size() - 1; }
		inline double	C2VRatio			() const { return ratio(orgs.size(), maxActive()); }
		inline bool		interrupted			() const { return intr; }
		inline bool		useTarget			() const { return lrn.stable && opts.target_phase_en; }
		inline bool		vsidsOnly			() const { return lrn.stable && opts.vsidsonly_en; }
		inline bool		vsidsEnabled		() const { return lrn.stable && opts.vsids_en; }
		inline bool		varsEnough			() const { assert(trail.size() < inf.maxVar); return (inf.maxVar - trail.size()) > lrn.nRefVars; }
		inline bool		canPreSigmify		() const { return opts.sigma_en; }
		inline bool		canMMD				() const { return lrn.rounds && varsEnough(); }
		inline bool		canRephase			() const { return opts.rephase_en && nConflicts > lrn.rephase_conf_max; }
		inline bool		canReduce			() const { return opts.reduce_en && learnts.size() && nConflicts >= lrn.reduce_conf_max; }
		inline bool		canCollect			() const { return cm.garbage() > (cm.size() * opts.gc_perc); }
		inline bool		canSubsume			() const {
			if (!opts.subsume_en) return false;
			if (nConflicts != lrn.lastreduce || nConflicts < lrn.subsume_conf_max) return false;
			return true;
		}
		inline bool		canMap				() const {
			if (DL()) return false;
			if (nConflicts < lrn.map_conf_max) return false;
			uint32 inactive = inf.maxFrozen + inf.maxMelted;
			assert(inactive <= inf.maxVar);
			if (inactive < opts.map_min) return false;
			return inactive >= (opts.map_perc * inf.maxVar);
		}
		inline bool		canSigmify			() const {
			if (!opts.sigma_live_en) return false;
			if (nConflicts <= lrn.sigma_conf_max) return false;
			if (sp->simplified >= opts.sigma_min) return true;
			return ((lrn.elim_marked - lrn.elim_lastmarked) > (int64(opts.sigma_min) << 4));
		}
		inline bool		satisfied			() const { 
			uint32 assigned = trail.size();
			if (sp->propagated < assigned) return false;
			return (assigned == inf.maxVar - inf.maxMelted);
		}
		inline void		markLearnt			() {
			assert(learntC.size() > 1);
			for (uint32* k = learntC + 1; k != learntC.end(); k++) {
				uint32 lit = *k;
				sp->marks[ABS(lit)] = SIGN(lit);
			}
		}
		inline void		unmarkLearnt		() {
			assert(learntC.size() > 1);
			for (uint32* k = learntC + 1; k != learntC.end(); k++)
				sp->marks[ABS(*k)] = UNDEFINED;
		}
		inline void		scaleVarAct			() {
			for (uint32 v = 1; v <= inf.maxVar; v++) activity[v] *= 1e-150;
			lrn.var_inc *= 1e-150;
		}
		inline void		initQueue			() {
			if (verbose == 4) PFLOG2(2, "  Initializing VMFQ Queue with %d variables..", inf.maxVar);
			else PFLOGN2(2, "  Initializing VMFQ Queue with %d variables..", inf.maxVar);
			for (uint32 v = 1; v <= inf.maxVar; v++)
				vmfq.init(v), vmfq.update(v, (bumps[v] = ++lrn.bumped));
			PFLDONE(2, 4);
		}
		inline void		initHeap			() {
			PFLOGN2(2, "  Initializing VSIDS Heap with %d variables..", inf.maxVar);
			for (uint32 v = 1; v <= inf.maxVar; v++)
				vsids.insert(v);
			PFLDONE(2, 5);
		}
		inline void		initVars			() {
			size_t maxSize = (size_t)inf.maxVar + 1;
			memset(sp->value, UNDEFINED, inf.nDualVars);
			memset(sp->marks, UNDEFINED, maxSize);
			memset(sp->ptarget, UNDEFINED, maxSize);
			memset(sp->pbest, UNDEFINED, maxSize);
			memset(sp->psaved, opts.polarity, maxSize);
			vorg.resize((uint32)maxSize);
			vorg[0] = 0;
			for (uint32 v = 1; v <= inf.maxVar; v++) {
				sp->level[v] = UNDEFINED, sp->source[v] = NOREF;
				vorg[v] = v;
			}
		}
		inline bool		canRestart			() {
			if (DL() < 2) return false;
			if (vibrate()) return lubyrest;
			if (nConflicts <= lrn.restarts_conf_max) return false;
			if (opts.mdmfusem_en && !lrn.rounds) MDMFuseMaster();
			return lbdrest.restart();
		}
		inline double	scale				(const double& val) {
			double c2v = C2VRatio();
			double factor = (c2v <= 2) ? 1.0 : log(c2v) / log(2);
			double newval = val * factor;
			if (newval < 1) newval = 1;
			return newval;
		}
		inline void		markLit				(const uint32& lit) { assert(lit > 1); sp->marks[ABS(lit)] = SIGN(lit); }
		inline void		unmarkLit			(const uint32& lit) { assert(lit > 1); sp->marks[ABS(lit)] = UNDEFINED; }
		inline bool		subsumed			(const uint32& lit) const { assert(lit > 1); return sp->marks[ABS(lit)] == SIGN(lit); }
		inline bool		selfsubsumed		(const uint32& lit) const { assert(lit > 1); return sp->marks[ABS(lit)] == !SIGN(lit); }
		inline LIT_ST	value				(const uint32& lit) const { assert(lit > 1); return sp->value[lit]; }
		inline LIT_ST	unassigned			(const uint32& lit) const { assert(lit > 1); return UNASSIGNED(sp->value[lit]); }
		inline LIT_ST	isFalse				(const uint32& lit) const { assert(lit > 1); return !sp->value[lit]; }
		inline LIT_ST	l2marker			(const uint32& lit) const { assert(lit > 1); return sp->marks[ABS(lit)]; }
		inline int		l2dl				(const uint32& lit) const { assert(lit > 1); return sp->level[ABS(lit)]; }
		inline C_REF	l2r					(const uint32& lit) const { assert(lit > 1); return sp->source[ABS(lit)]; }
		inline LIT_ST	isTrue				(const uint32& lit) const { 
			assert(lit > 1);
			LIT_ST val = sp->value[lit];
			return val && !UNASSIGNED(val);
		}
		inline void		varBumpHeap			(const uint32& v, const double& norm_act) {
			assert(v && v <= inf.maxVar);
			assert(norm_act);
			assert(activity[v] == 0.0);
			activity[v] = norm_act;
			vsids.bump(v);
		}
		inline void		varBumpHeap			(const uint32& v) {
			assert(v && v <= inf.maxVar);
			assert(lrn.var_inc);
			assert(activity[v] <= 1e150);
			double newAct = activity[v] + lrn.var_inc;
			if (newAct > 1e150) {
				scaleVarAct();
				assert(lrn.var_inc <= 1e150);
				newAct = activity[v] + lrn.var_inc;
			}
			assert(newAct <= 1e150);
			activity[v] = newAct;
			vsids.bump(v);
		}
		inline void		varBumpQueue		(const uint32& v) {
			assert(v && v <= inf.maxVar);
			if (!vmfq.next(v)) return;
			vmfq.toFront(v);
			assert(lrn.bumped && lrn.bumped != INT64_MAX);
			bumps[v] = ++lrn.bumped;
			PFLOG2(4, " Variable %d moved to queue front & bumped to %lld", v, lrn.bumped);
			if (!sp->locked[v]) vmfq.update(v, bumps[v]);
		}
		inline void		varBumpQueueNU		(const uint32& v) {
			assert(v && v <= inf.maxVar);
			if (!vmfq.next(v)) return;
			assert(lrn.bumped && lrn.bumped != INT64_MAX);
			vmfq.toFront(v);
			bumps[v] = ++lrn.bumped;
			PFLOG2(4, " Variable %d moved to queue front & bumped to %lld", v, lrn.bumped);
		}
		inline uint32	rscore				(const uint32& v) {
			assert(v && v <= inf.maxVar);
			return (!occurs[v].ps || !occurs[v].ns) ? occurs[v].ps | occurs[v].ns : occurs[v].ps * occurs[v].ns;
		}
		inline void		attachWatch			(const C_REF& ref, const CLAUSE& c) {
			assert(c[0] > 1);
			assert(c[1] > 1);
			assert(ref < NOREF);
			int sz = c.size();
			assert(sz > 1);
			wt[FLIP(c[0])].push(WATCH(ref, sz, c[1]));
			wt[FLIP(c[1])].push(WATCH(ref, sz, c[0]));
		}
		inline void		attachWatch			(const uint32& lit, const uint32& imp, const C_REF& ref, const int& size) {
			assert(lit != imp);
			assert(ref < NOREF);
			assert(size > 1);
			wt[FLIP(lit)].push(WATCH(ref, size, imp));
		}
		inline void		detachWatch			(const uint32& lit, const C_REF& ref) {
			assert(lit > 1);
			assert(ref < NOREF);
			WL& ws = wt[lit];
			if (ws.empty()) return;
			WATCH* i = ws, *end = ws.end();
			for (WATCH* j = i; j != end; j++) {
				const WATCH& w = *i++ = *j;
				if (w.ref == ref) i--;
			}
			assert(i + 1 == end);
			ws.resize(int(i - ws));
		}
		inline void		enqueue				(const uint32& lit, const int& pLevel = 0, const C_REF src = NOREF) {
			assert(lit > 1);
			uint32 v = ABS(lit);
			if (!pLevel) sp->vstate[v] = FROZEN, inf.maxFrozen++;
			sp->psaved[v] = SIGN(lit);
			sp->source[v] = src;
			sp->locked[v] = 1;
			sp->level[v] = pLevel;
			sp->value[lit] = 1;
			sp->value[FLIP(lit)] = 0;
			trail.push(lit);
			if (!wt.empty()) {
				WL& ws = wt[lit];
				if (!ws.empty()) {
#if _WIN32
					PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, &ws[0]);
#else
					__builtin_prefetch(&ws[0], 0, 1);
#endif
				}
			}
			PFLNEWLIT(this, 3, src, lit);
		}
		inline void		enqueueImp			(const uint32& lit, const C_REF src) {
			int assignLevel = opts.chrono_en ? forcedLevel(lit, src) : DL();
			enqueue(lit, assignLevel, src);
		}
		inline void		enqueueOrg			(const uint32& lit) {
			assert(lit > 1);
			assert(!DL());
			uint32 v = ABS(lit);
			sp->vstate[v] = FROZEN, inf.maxFrozen++;
			sp->locked[v] = 1;
			sp->level[v] = 0;
			sp->value[lit] = 1;
			sp->value[FLIP(lit)] = 0;
			trail.push(lit);
			PFLNEWLIT(this, 3, NOREF, lit);
		}
		inline int		forcedLevel			(const uint32& lit, const C_REF& src) {
			assert(opts.chrono_en);
			assert(REASON(src));
			int fl = 0;
			CLAUSE& c = cm[src];
			for (uint32* k = c; k != c.end(); k++) {
				uint32 other = *k;
				if (other == lit) continue;
				assert(!unassigned(other));
				int otherLevel = l2dl(other);
				if (otherLevel > fl) fl = otherLevel;
			}
			return fl;
		}
		inline void		hist				(CLAUSE& c) {
			for (uint32* k = c; k != c.end(); k++) {
				uint32 lit = *k;
				assert(lit > 1);
				if (SIGN(lit)) occurs[ABS(lit)].ns++;
				else occurs[ABS(lit)].ps++;
			}
		}
		inline void		hist				(const BCNF& cnf, const bool& rst = false) {
			if (cnf.empty()) return;
			if (rst) for (uint32 i = 0; i < occurs.size(); i++) occurs[i] = { 0 , 0 };
			for (uint32 i = 0; i < cnf.size(); i++) hist(cm[cnf[i]]);
			assert(occurs[0].ps == 0 && occurs[0].ns == 0);
		}
		inline void		savePhases			(LIT_ST* to)
		{
			for (uint32 v = 1; v <= inf.maxVar; v++)
				to[v] = sp->psaved[v];
		}
		inline void		wrProof				(const Byte& byte) { proofFile << byte; }
		inline void		wrProof				(uint32* lits, const int& len) {
			assert(len > 0);
			uint32* lit = lits, * end = lits + len;
			while (lit != end) {
				register uint32 b = 0;
				/*if (mapped) b = vmap.revLit(*lit++);
				else*/ b = *lit++;
				while (b > 127) { wrProof(Byte(128 | (b & 127))); b >>= 7; }
				wrProof(Byte(b));
			}
		}
		//==============================================
		bool	keeping				(CLAUSE&);
		CL_ST	rooted				(CLAUSE&);
		void	markLits			(CLAUSE&);
		void	unmarkLits			(CLAUSE&);
		void	markSubsume			(CLAUSE&);
		void	bumpShrunken		(CLAUSE&);
		void	shrinkClause		(CLAUSE&, const int&);
		void	shrinkClause		(const C_REF&);
		C_REF	newClause			(const Lits_t&, const CL_ST& type = ORIGINAL);
		void	newClause			(SCLAUSE&);
		void	markSubsume			(SCLAUSE&);
		bool	toClause			(Lits_t&, Lits_t&, char*&);
		void	removeClause		(const C_REF&);
		void	backtrack			(const int& bt_level = 0);
		C_REF	backjump			(const int&);
		void	recycle				(CMM&);
		void	filter				(BCNF&, CMM&);
		void	filter				(BCNF& s, BCNF& d, const CL_ST& t);
		void	filter				(BCNF&);
		void	shrink				(BCNF&);
		void	schedule			(BCNF&);
		void	attachBins			(BCNF&);
		void	attachNonBins		(BCNF&);
		void	attachClauses		(BCNF&);
		void	subsumeLearnt		(const C_REF&);
		uint32	makeAssign			(const uint32&, const bool& tphase = false);
		bool	minimize			(const uint32&, const int& depth = 0);
		void	reduceLearnts		(const bool& sizeonly = false);
		void	reduceTop			(const bool&);
		void	rebuildWT			(const bool&);
		void	pumpFrozen			();
		void	allocSolver			();
		void	resetSolver			();
		void	initSolver			();
		void	killSolver			();
		void	varOrder			();
		bool	shrink				();
		void	protectReasons		();
		void	unprotectReasons	();
		void	recycleWT			();
		void	recycle				();
		void	reduce				();
		void	rephase				();
		void	subsume				();
		bool	subsumeAll			();
		void	minimizeBin			();
		void	minimize			();
		void	solve				();
		int		reuse				();
		bool	vibrate				();
		void	restart				();
		void	analyze				();
		bool	chronoAnalyze		();
		bool	BCPChronoRoot		();
		bool	BCP					();
		uint32	nextVSIDS			();
		uint32	nextVMFQ			();
		void	MDMFuseMaster		();
		void	MDMFuseSlave		();
		void	MDMInit				();
		void	MDM					();
		void	eligibleVSIDS		();
		void	eligibleVMFQ		();
		void	decide				();
		void	report				();
		void	wrapup				();
		bool	parser				();
		void	map					(BCNF&);
		void	map					(WL&);
		void	map					(WT&);
		void	map					(const bool& = false);
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
			return !opts.phases && !(opts.all_en || opts.ere_en);
		}
		inline bool		sleeping			() {
			return (cnfstate == UNSAT || sigState == AWAKEN_FAIL || sigState == CNFALLOC_FAIL);
		}
		inline bool		reallocFailed		() {
			return (sigState == OTALLOC_FAIL || sigState == CNFALLOC_FAIL);
		}
		inline bool		reallocCNF			() {
			int times = phase + 1;
			if (times > 1 && times != opts.phases && (times % opts.shrink_rate) == 0) {
				size_t maxAddedCls = opts.ve_en ? inf.nClauses : 0;
				size_t maxAddedLits = opts.ve_en ? inf.nLiterals * opts.lits_mul : 0;
				PFLOG2(2, " Maximum added clauses/literals = %zd/%zd", maxAddedCls, maxAddedLits);
				if (!cumm.resizeCNF(cnf, inf.nClauses + maxAddedCls, inf.nLiterals + maxAddedLits)) {
					sigState = CNFALLOC_FAIL, compacted = false;
					return false;
				}
				compacted = true;
			}
			else cumm.cacheCNFPtr(cnf), compacted = false;
			return true;
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
			for (uint32 v = 0; v < vars->numPVs; v++) if (sp->frozen[vars->pVars->at(v)]) return false;
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
			if (host) countMelted(), countCls(1);
			else SIGmA::countFinal(cnf, vars->gstats, sp->vstate);
			inf.nClauses = inf.n_cls_after;
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
		inline bool		reallocOT			(const cudaStream_t& stream = 0) {
			assert(inf.nLiterals);
			calcOccurs(inf.nLiterals);
			if (!cumm.resizeOTAsync(ot, inf.nLiterals, stream)) { sigState = OTALLOC_FAIL; return false; }
			return true;
		}
		inline void		reflectCNF			(const cudaStream_t& s1, const cudaStream_t& s2) {
			S_REF len1 = hcnf->data().size - dataoff;
			if (!len1) return;
			uint32 len2 = hcnf->size() - csoff;
			CHECK(cudaMemcpyAsync(cumm.cnfMem() + dataoff, hcnf->data().mem + dataoff, len1 * hc_bucket, cudaMemcpyHostToDevice, s1));
			CHECK(cudaMemcpyAsync(cumm.csMem() + csoff, hcnf->csData() + csoff, len2 * sizeof(S_REF), cudaMemcpyHostToDevice, s2));
			dataoff = hcnf->data().size, csoff = hcnf->size();
		}
		inline void		cacheUnits			(const cudaStream_t& stream) {
			sync(stream);
			if (vars->nUnits = vars->tmpObj.size()) CHECK(cudaMemcpyAsync(vars->cachedUnits, cumm.unitsdPtr(),
				vars->nUnits * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
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
		inline bool		stop				(const int64 lr) {
			return (phase == opts.phases) || (lr <= opts.lits_min && phase > 2);
		}
		//===========================================//
		void			optSimp				();
		void			varReorder			();
		void			newBeginning		();
		void			sigmify				();
		bool			LCVE				();
		void			postVE				();
		void			VE					();
		void			ERE					();
		void			BCE					();
		void			HSE					();
		void			masterFree			();
		void			slavesFree			();
		bool			prop				();
		bool			propHost			();
		void			awaken				(const bool& strict = 0);
		void			histSimp			(const uint32&);
		void			calcOccurs			(const uint32&);
		void			extract				(CNF*, BCNF&);
		void			createOTHost		(HOT&);
		inline bool		propClause			(SCLAUSE&, const uint32&);
		inline void		depFreeze			(OL&, const uint32&, const uint32&, const uint32&);
		inline void		cacheCNF			(const cudaStream_t&, const cudaStream_t&);
		inline bool		enqeueCached		(const cudaStream_t&);
		inline void		cleanProped			();
		inline void		cleanDynamic		();
		inline void		initSimp			();
		inline void		sigmaDelay			();
		inline void		writeBack			();
		//==========================================//
		//			       Printers                 //
		//==========================================//
		void printStats(const bool& _p = true, const Byte& _t = ' ', const char* _c = CNORMAL);
		void printVars(const uint32* arr, const uint32& size, const LIT_ST& type = 'x');
		void printClause(const Lits_t&);
		void printTrail(const uint32& off = 0);
		void printCNF(const BCNF&, const int& off = 0);
		void printOL(const OL& list, const bool& host = 0);
		void printOL(const uint32& lit, const bool& host = 0);
		void printWL(const uint32&, const bool& bin = 0);
		void printWatched(const uint32&);
		void printBinaries(const uint32&);
		void printTable();
		void printWT();
		void printOT();
		void printHeap();
		void printSource();
		void printLearnt();
	};
	extern ParaFROST* pfrost;
}

#endif 