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
#include "pfsort.h"
#include "pfvec.h"
#include "pfargs.h"
#include "pfheap.h"
#include "pfqueue.h"
#include "pfrestart.h"
#include "pfvmap.h"
#include "pfsolvertypes.h"
#include "pfmodel.h"

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
		WT				wt, wtBin;
		BCNF			orgs, learnts, reduced;
		VMAP			vmap;
		MODEL			model;
		C_REF			conflict;
		Lits_t			learntC;
		Vec<OCCUR>		occurs;
		Vec<int64>		varBumps;
		Vec<double>		varAct;
		uVec1D			trail, trail_lens, eligible, analyzed, toSimp;
		HEAP<HEAP_CMP>	vHeap;
		QUEUE			vQueue;
		LBDREST			lbdrest;
		LUBYREST		lubyrest;
		int64			nConflicts;
		uint32			starts, minLevel;
		CNF_ST			cnfstate;
		size_t			solLineLen;
		string			solLine;
		std::ofstream	proofFile;
		bool			intr, units_f, search_guess;
	public:
		//============== inline methods ===============
						~ParaFROST			() { if (sigma_en || sigma_live_en) masterFree(), slavesFree(), CHECK(cudaDeviceReset()); }
		inline void		interrupt			() { intr = true; }
		inline void		incDL				() { trail_lens.push(trail.size()); }
		inline void		decayVarAct			() { lrn.var_inc *= (1.0 / lrn.var_decay); }
		inline uint32	maxOrgs				() const { return orgs.size() + inf.nOrgBins; }
		inline uint32	maxLearnts			() const { return learnts.size() + inf.nLearntBins; }
		inline uint32	DL					() const { return trail_lens.size(); }
		inline double	drand				() const { return ((double)rand() / (double)RAND_MAX); }
		inline double	C2VRatio			() const { return ratio(maxOrgs(), maxActive()); }
		inline bool		interrupted			() const { return intr; }
		inline bool		satisfied			() const { assert(sp->propagated == trail.size()); return (trail.size() == inf.maxVar - inf.maxMelted); }
		inline bool		useTarget			() const { return lrn.stable && target_phase_en; }
		inline bool		vsidsOnly			() const { return lrn.stable && vsidsonly_en; }
		inline bool		vsids				() const { return lrn.stable && vsids_en; }
		inline bool		varsEnough			() const { assert(trail.size() < inf.maxVar); return (inf.maxVar - trail.size()) > lrn.nRefVars; }
		inline bool		canPreSigmify		() const { return sigma_en; }
		inline bool		canMMD				() const { return lrn.rounds && varsEnough(); }
		inline bool		canRephase			() const { return rephase_en && nConflicts > lrn.rephase_conf_max; }
		inline bool		canReduce			() const { return reduce_en && learnts.size() && nConflicts >= lrn.reduce_conf_max; }
		inline bool		canMap				() const { return inf.maxFrozen > map_min || inf.maxMelted > map_min; }
		inline bool		canShrink			() const { 
			return DL() == ROOT_LEVEL && int(trail.size() - sp->simplified) >= shrink_min;
		}
		inline bool		canSigmify			() const { 
			return sigma_live_en && DL() == ROOT_LEVEL && nConflicts > lrn.sigma_conf_max && lrn.lastsimplified >= sigma_min;
		}
		inline bool		abortSub			(const uint32& level) const { return (MAPHASH(level) & minLevel) == 0; }
		inline LIT_ST	value				(const uint32& lit) const { assert(lit > 1); return sp->value[lit]; }
		inline LIT_ST	unassigned			(const uint32& lit) const { assert(lit > 1); return sp->value[lit] & NOVAL_MASK; }
		inline LIT_ST	isTrue				(const uint32& lit) const { assert(lit > 1); return !unassigned(lit) && sp->value[lit]; }
		inline LIT_ST	isFalse				(const uint32& lit) const { assert(lit > 1); return !sp->value[lit]; }
		inline int		l2dl				(const uint32& lit) const { assert(lit > 1); return sp->level[l2a(lit)]; }
		inline C_REF	l2r					(const uint32& lit) const { assert(lit > 1); return sp->source[l2a(lit)]; }
		inline uint32	l2hl				(const uint32& lit) const { return MAPHASH(l2dl(lit)); }
		inline int		calcLBD				(Lits_t& c) {
			stats.marker++;
			register int lbd = 0;
			for (int i = 0; i < c.size(); i++) {
				int litLevel = l2dl(c[i]);
				if (sp->board[litLevel] != stats.marker) { sp->board[litLevel] = stats.marker; lbd++; }
			}
			return lbd;
		}
		inline int		calcLBD				(CLAUSE& c) {
			stats.marker++;
			register int lbd = 0;
			for (int i = 0; i < c.size(); i++) {
				int litLevel = l2dl(c[i]);
				if (sp->board[litLevel] != stats.marker) { sp->board[litLevel] = stats.marker; lbd++; }
			}
			return lbd;
		}
		inline uint32	rscore				(const uint32& v) {
			assert(v && v <= inf.maxVar);
			return (!occurs[v].ps || !occurs[v].ns) ? occurs[v].ps | occurs[v].ns : occurs[v].ps * occurs[v].ns;
		}
		inline uint32	what				(const uint32& v, const bool& tphase = false) {
			assert(v < NOVAR);
			LIT_ST pol = UNDEFINED;
			if (pol == UNDEFINED && tphase) pol = sp->ptarget[v];
			if (pol == UNDEFINED) pol = sp->psaved[v];
			if (pol == UNDEFINED) pol = polarity;
			assert(pol >= 0);
			return v2dec(v, pol);
		}
		inline void		bumpVariable		(const uint32& v) {
			assert(v && v <= inf.maxVar);
			if (vsids()) varBumpHeap(v);
			else varBumpQueue(v);
		}
		inline void		varBumpQueue		(const uint32& v) {
			assert(v && v <= inf.maxVar);
			if (!vQueue.next(v)) return;
			assert(lrn.bumped && lrn.bumped != INT64_MAX);
			vQueue.toFront(v);
			varBumps[v] = ++lrn.bumped;
			PFLOG2(4, " Variable %d moved to queue front & bumped to %lld", v, lrn.bumped);
			if (!sp->locked[v]) vQueue.update(v, varBumps[v]);
		}
		inline void		varBumpQueueNU		(const uint32& v) {
			assert(v && v <= inf.maxVar);
			if (!vQueue.next(v)) return;
			assert(lrn.bumped && lrn.bumped != INT64_MAX);
			vQueue.toFront(v);
			varBumps[v] = ++lrn.bumped;
			PFLOG2(4, " Variable %d moved to queue front & bumped to %lld", v, lrn.bumped);
		}
		inline void		varBumpHeap			(const uint32& v, double norm_act) {
			assert(v && v <= inf.maxVar);
			assert(norm_act);
			assert(varAct[v] == 0.0);
			varAct[v] += norm_act;
			vHeap.bump(v);
		}
		inline void		varBumpHeap			(const uint32& v) {
			assert(v && v <= inf.maxVar);
			assert(lrn.var_inc);
			if ((varAct[v] += lrn.var_inc) > 1e150) scaleVarAct();
			assert(varAct[v] <= 1e150);
			assert(lrn.var_inc <= 1e150);
			vHeap.bump(v);
		}
		inline void		enqueue				(const uint32& lit, const int& pLevel = ROOT_LEVEL, const uint32 src = NOREF) {
			assert(lit > 1);
			uint32 v = l2a(lit);
			if (!search_guess) sp->psaved[v] = sign(lit);
			if (!pLevel) sp->vstate[v] = FROZEN, inf.maxFrozen++;
			if (src != NOREF) cm[src].markReason();
			sp->source[v] = src;
			sp->locked[v] = 1;
			sp->level[v] = pLevel;
			sp->value[lit] = 1;
			sp->value[flip(lit)] = 0;
			trail.push(lit);
			PFLNEWLIT(this, 4, src, lit);
		}
		inline void		cancelAssign		(const uint32& lit) {
			assert(lit > 1);
			uint32 v = l2a(lit);
			C_REF& r = sp->source[v];
			sp->value[lit] = UNDEFINED;
			sp->value[flip(lit)] = UNDEFINED;
			sp->locked[v] = 0;
			if (!vHeap.has(v)) vHeap.insert(v);
			if (vQueue.bumped() < varBumps[v]) vQueue.update(v, varBumps[v]);
			if (r != NOREF) cm[r].initReason(), r = NOREF;
			PFLOG2(4, " Literal %d@%d cancelled", l2i(lit), sp->level[v]);
			sp->level[v] = UNDEFINED;
		}
		inline void		analyzeLit			(const uint32& lit, int& track) {
			assert(lit > 1);
			uint32 v = l2a(lit);
			if (sp->seen[v]) return;
			int confLevel = DL(), litLevel = sp->level[v];
			if (litLevel > ROOT_LEVEL) {
				assert(litLevel <= confLevel);
				if (litLevel == confLevel) track++;
				else if (litLevel < confLevel) learntC.push(lit);
				analyzed.push(v);
				sp->seen[v] = 1;
			}
		}
		inline void		analyzeReason		(const C_REF& r, const uint32& parent, int& track) {
			CLAUSE& c = cm[r];
			PFLCLAUSE(4, c, " Analyzing reason");
			if (parent && c.size() == 2 && isFalse(*c)) {
				assert(isTrue(c[1]));
				c.swapWatched();
			}
			if (c.learnt() && c.size() > 2) bumpClause(c);
			for (int j = (parent == 0) ? 0 : 1; j < c.size(); j++)
				analyzeLit(c[j], track);
		}
		inline void		bumpClause			(CLAUSE& c) {
			assert(c.learnt());
			assert(c.size() > 2);
			int old_lbd = c.lbd();
			if (old_lbd <= lbd_tier1) return; // always keep Tier1 value
			int new_lbd = calcLBD(c);
			if (new_lbd < old_lbd) { // try to update old LBD
				if (old_lbd > lbd_tier2 && new_lbd <= lbd_tier2) c.initTier2();
				else c.initTier3(); // set as Tier 3 
				c.set_lbd(new_lbd);
				PFLCLAUSE(4, c, " Bumping clause with LBD %d ", new_lbd);
			}
			else if (c.usage() && old_lbd <= lbd_tier2) c.initTier2(); // set as Tier2
			else c.initTier3(); // set as Tier 3 
		}
		inline void		bumpShrinked		(CLAUSE& c) {
			assert(c.learnt());
			assert(c.size() > 2);
			int old_lbd = c.lbd();
			if (old_lbd <= lbd_tier1) return; // always keep Tier1 value
			int new_lbd = std::min(c.size() - 1, old_lbd);
			if (new_lbd >= old_lbd) return;
			if (old_lbd > lbd_tier2 && new_lbd <= lbd_tier2) c.initTier2();
			c.set_lbd(new_lbd);
			PFLCLAUSE(4, c, " Bumping clause with LBD %d ", new_lbd);
		}
		inline void		depleteSource		(CLAUSE& c) {
			assert(c.deleted());
			if (c.reason()) {
				if (c.size() == 2 && isFalse(*c)) {
					assert(isTrue(c[1]));
					assert(l2dl(c[1]) == ROOT_LEVEL);
					sp->source[l2a(c[1])] = NOREF;
				}
				else {
					assert(isTrue(*c));
					assert(l2dl(*c) == ROOT_LEVEL);
					sp->source[l2a(*c)] = NOREF;
				}
			}
		}
		inline void		hist				(const C_REF& r) {
			CLAUSE& c = cm[r];
			for (int k = 0; k < c.size(); k++) {
				assert(c[k] > 0);
				if (sign(c[k])) occurs[l2a(c[k])].ns++;
				else occurs[l2a(c[k])].ps++;
			}
		}
		inline void		hist				(BCNF& cnf, const bool& rst = false) {
			if (cnf.size() == 0) return;
			if (rst) for (uint32 i = 0; i < occurs.size(); i++) occurs[i] = { 0 , 0 };
			for (uint32 i = 0; i < cnf.size(); i++) hist(cnf[i]);
			assert(occurs[0].ps == 0 && occurs[0].ns == 0);
		}
		inline void		removeSat			(BCNF& cnf) {
			if (cnf.size() == 0) return;
			uint32 n = 0;
			for (uint32 i = 0; i < cnf.size(); i++) {
				C_REF r = cnf[i];
				assert(!cm[r].moved());
				if (satisfied(cm[r])) removeClause(r);
				else cnf[n++] = r;
			}
			cnf.resize(n);
		}
		inline void		shrink				(BCNF& cnf) {
			if (cnf.size() == 0) return;
			uint32 n = 0;
			for (uint32 i = 0; i < cnf.size(); i++) {
				C_REF r = cnf[i];
				assert(!cm[r].moved());
				if (satisfied(cm[r])) removeClause(r);
				else {
					shrinkClause(r);
					if (cm[r].size() > 2) cnf[n++] = r;
				}
			}
			cnf.resize(n);
		}
		inline void		extractBins			(BCNF& dest) {
			for (uint32 lit = 2; lit < wtBin.size(); lit++) {
				WL& ws = wtBin[lit];
				if (ws.empty()) continue;
				for (WATCH* w = ws; w != ws.end(); w++) {
					C_REF r = w->ref;
					if (cm[r].moved()) continue;
					dest.push(r);
					cm[r].markMoved();
				}
			}
		}
		inline void		mapBins				(BCNF& cnf) {
			assert(!vmap.empty());
			for (uint32 i = 0; i < cnf.size(); i++) {
				CLAUSE& c = cm[cnf[i]];
				assert(c.moved());
				c.initMoved();
				vmap.mapBinClause(c);
			}
		}
		inline void		map					(BCNF& cnf) {
			assert(!vmap.empty());
			for (uint32 i = 0; i < cnf.size(); i++) {
				CLAUSE& c = cm[cnf[i]];
				vmap.mapClause(c);
			}
		}
		inline void		map					(WL& ws) {
			if (ws.empty()) return;
			for (WATCH* w = ws; w != ws.end(); w++)
				w->imp = vmap.mapLit(w->imp);
		}
		inline void		reduceWeight		(double& val) {
			double orgSize = maxOrgs();
			if (orgSize > 1e5) {
				val *= log(orgSize / 1e4) / log(10);
				if (val < 1.0) val = 1.0;
			}
		}
		inline void		sigmaWeight			(double& val) {
			double c2v = C2VRatio();
			val *= (c2v <= 2) ? 1.0 : log(c2v) / log(2);
			if (val < 1.0) val = 1.0;
		}
		inline void		bumpVariables		() {
			bool vsidsEn = vsids();
			if (!vsidsEn) Sort(analyzed, ANALYZE_CMP(varBumps));
			for (uint32 i = 0; i < analyzed.size(); i++)
				bumpVariable(analyzed[i]);
			if (vsidsEn) decayVarAct();
			analyzed.clear();
		}
		inline void		scaleVarAct			() {
			for (uint32 v = 1; v <= inf.maxVar; v++) varAct[v] *= 1e-150;
			lrn.var_inc *= 1e-150;
		}
		inline void		varOrgPhase			() {
			memset(sp->psaved, polarity, inf.maxVar + 1ULL);
			lrn.lastrephased = ORGPHASE;
		}
		inline void		varFlipPhase		() {
			LIT_ST* end = sp->psaved + inf.maxVar + 1;
			for (LIT_ST* s = sp->psaved; s != end; s++)
				*s ^= NEG_SIGN;
			lrn.lastrephased = FLIPPHASE;
		}
		inline void		varBestPhase		() {
			for (uint32 v = 1; v <= inf.maxVar; v++)
				if (sp->pbest[v] != UNDEFINED)
					sp->psaved[v] = sp->pbest[v];
			lrn.lastrephased = BESTPHASE;
		}
		inline void		MDMFuseMaster		() {
			if (mdm_rounds && nConflicts >= lrn.mdm_conf_max) {
				lrn.rounds = mdm_rounds;
				lrn.mdm_conf_max = nConflicts + mdm_minc;
				mdm_minc <<= 1;
				PFLOG2(2, " MDM limit increased by (conflicts: %lld, Increment: %d) to %lld\n", nConflicts, mdm_minc, lrn.mdm_conf_max);
			}
		}
		inline void		MDMFuseSlave		() {
			if (mdm_rounds && mdm_div) {
				int q = int(nConflicts / mdm_div) + 1;
				if (starts % q == 0) {
					lrn.nRefVars = 0, lrn.rounds = mdm_rounds;
					PFLOG2(2, " Starts: %d, conflicts: %lld, dividor: %d, q: %d, rounds: %d", starts, nConflicts, mdm_div, q, lrn.rounds);
				}
				if (nConflicts % mdm_freq == 0) mdm_div += mdm_sinc;
			}
		}
		inline bool		canRestart			() {
			if (vibrate()) return lubyrest.restart();
			if (nConflicts <= lrn.restarts_conf_max) return false;
			if (mdmfusem_en && !lrn.rounds) MDMFuseMaster();
			return lbdrest.restart();
		}
		inline bool		verifyMDM			() {
			for (uint32 i = sp->propagated; i < trail.size(); i++) {
				uint32 v = l2a(trail[i]);
				if (sp->frozen[v]) {
					PFLOG0("");
					PFLOGEN("decision(%d) is elected and frozen", v);
					printWatched(v);
					return false;
				}
			}
			return true;
		}
		inline bool		verifySeen			() {
			for (uint32 v = 0; v <= inf.maxVar; v++) {
				if (sp->seen[v]) {
					PFLOG0("");
					PFLOGEN("seen(%d) is not unseen", v);
					printWatched(v);
					return false;
				}
			}
			return true;
		}
		inline bool		guessing			(const uint32& guess) {
			assert(guess && guess <= inf.nDualVars + 1);
			if (unassigned(guess)) {
				incDL();
				enqueue(guess, DL());
				if (BCP() != NOREF) {
					cancelAssigns();
					return false;
				}
			}
			return true;
		}
		inline bool		satisfied			(const CLAUSE& c) {
			if (c.deleted()) return true;
			for (int k = 0; k < c.size(); k++) if (isTrue(c[k]))  return true;
			return false;
		}
		inline void		recycleWL			(WL& ws) {
			for (WATCH* w = ws; w != ws.end(); w++)
				w->ref = newClause((*hcnf)[w->ref]);
		}
		inline void		recycleWL			(WL& ws, CMM& new_cm) {
			for (WATCH* w = ws; w != ws.end(); w++)
				cm.move(w->ref, new_cm);
		}
		inline void		detachClause		(C_REF& r, const bool& gc = true) {
			CLAUSE& c = cm[r];
			WT& ws = c.size() == 2 ? wtBin : wt;
			if (gc) ws.collect(flip(c[0])), ws.collect(flip(c[1]));
			else ws.remWatch(flip(c[0]), r), ws.remWatch(flip(c[1]), r);
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
		inline void		printTrail			(const uint32& off = 0) {
			assert(off < trail.size());
			uint32* x = trail + off;
			PFLOGN1(" Trail (size = %d)->[", trail.size());
			for (uint32 i = 0; i < trail.size(); i++) {
				printf("%6d@%-5d", l2i(x[i]), l2dl(x[i]));
				if (i && i < trail.size() - 1 && i % 8 == 0) { putc('\n', stdout); PFLOGN0("\t\t\t"); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		inline void		printCNF			(const BCNF& cnf, const int& off = 0) {
			PFLOG1("\tHost CNF(size = %d)", cnf.size());
			for (uint32 c = off; c < cnf.size(); c++) {
				if (cm[cnf[c]].size() > 0) {
					PFLCLAUSE(1, cm[cnf[c]], " C(%d)->", c);
				}
			}
		}
		inline void		printWatched		(const uint32& v) {
			uint32 p = v2l(v);
			wtBin.print(p), wtBin.print(neg(p));
			wt.print(p), wt.print(neg(p));
		}
		inline void		printStats			(const bool& _p = true, const Byte& type = ' ') {
			if (verbose == 1 && _p) {
				int l2c = learnts.size() == 0 ? 0 : int(inf.nLearntLits / learnts.size());
				solLine[0] = type;
				PFLOGN1(solLine.c_str(), maxActive(), orgs.size(), inf.nOrgBins, inf.nLiterals,
					nConflicts, starts - 1,
					learnts.size(), inf.nGlues + inf.nLearntBins, l2c);
				REPCH(' ', RULELEN - solLineLen), putc('|', stdout), putc('\n', stdout);
			}
		}
		inline void		printTable			() {
			const char* header = " Progress ";
			size_t len = strlen(header);
			if (RULELEN < len) PFLOGE("ruler length is smaller than the table title");
			size_t gap = (RULELEN - strlen(header)) / 2;
			PFLOGN0(""); REPCH('-', gap), fprintf(stdout, "%s", header);
			REPCH('-', gap); putc('|', stdout), putc('\n', stdout);
			string h = "";
			const char* leq = u8"\u2264";
#ifdef _WIN32
			SetConsoleOutputCP(65001);
#endif
			h = "                    ORG                    Conflicts   Restarts            Learnt";
			if (RULELEN < h.size()) PFLOGE("ruler length is smaller than the table header");
			PFLOGN0(h.c_str()); REPCH(' ', RULELEN - h.size()); putc('|', stdout), putc('\n', stdout);
			h = "       V        C          B          L                            C       BG(%s2)     L/C";
			if (RULELEN < h.size()) PFLOGE("ruler length is smaller than the table header");
			PFLOGN1(h.c_str(), leq); REPCH(' ', RULELEN - h.size() + 1); putc('|', stdout), putc('\n', stdout);
			solLine = "  %9d %9d %9d %10d %10d %8d %9d %9d %6d", solLineLen = 89;
			if (RULELEN < solLineLen) PFLOGE("ruler length is smaller than the progress line");
			PFLOGR('-', RULELEN);
		}
		inline void		printHeap			() {
			PFLOG1(" Heap (size = %d):", vHeap.size());
			for (uint32 i = 0; i < vHeap.size(); i++)
				PFLOG1(" h(%d)->(v: %d, a: %g)", i, vHeap[i], varAct[vHeap[i]]);
		}
		inline void		printLearnt			() {
			PFLOGN0(" Learnt(");
			for (int i = 0; i < learntC.size(); i++)
				printf("%d@%-6d", l2i(learntC[i]), l2dl(learntC[i]));
			putc(')', stdout), putc('\n', stdout);
		}
		//==============================================
		void	backJump			(const int&);
		void	cancelAssigns		(const int& = ROOT_LEVEL);
		void	resetSolver			(const bool& rst = true);
		void	histBins			(const bool& rst = false);
		void	savePhases			(const int&);
		void	savePhases			(LIT_ST*);
		void	savePhases			(LIT_ST*, const uint32&, const uint32&);
		void	newClause			(Lits_t&, const CL_ST& type = ORIGINAL);
		void	removeClause		(C_REF&, const bool& gc = true);
		void	shrinkClause		(C_REF&);
		void	recycle				(CMM&);
		bool	selfsub				(const uint32&, uint32*&);
		bool	depFreeze			(WL&, const uint32&);
		bool	valid				(WL&, WL&, const uint32&);
		void	pumpFrozenHeap		(const uint32&);
		void	pumpFrozenQue		(const uint32&);
		void	pumpFrozen			();
		void	initQueue			();
		void	initHeap			();
		void	initSolver			();
		void	killSolver			();
		void	recycle				();
		void	varOrder			();
		void	shrinkWT			();
		void	shrinkTop			();
		void	shrink				();
		void	reduce				();
		void	reduceLearnts		();
		void	rephase				();
		void	whereToJump			();
		void	selfsubBin			();
		void	selfsub				();
		void	analyze				();
		void	solve				();
		void	restart				();
		bool	vibrate				();
		C_REF	BCP					();
		uint32	nextVSIDS			();
		uint32	nextVMFQ			();
		void	MDMInit				();
		void	MDM					();
		void	eligibleVSIDS		();
		void	eligibleVMFQ		();
		void	decide				();
		void	guess				();
		void	allNegs				();
		void	allPoss				();
		void	forwardNegs			();
		void	forwardPoss			();
		void	backwardNegs		();
		void	backwardPoss		();
		void	printReport			();
		void	allocVars			();
		void	wrapup				();
		int64	sysMemUsed			();
		CNF_ST	parser				();
		void	mapBins				();
		void	map					(WT&);
		void	map					(const bool& sigmified = 0);
				ParaFROST			(const string&);
		//==========================================//
		//             Solver options               //
		//==========================================//
		string	proof_path;
		int64	stabrestart_r, stabrestart_inc;
		double	var_inc, var_decay;
		double	lbdrate, gc_perc, reduce_perc;
		uint32	map_min;
		int		seed, timeout, prograte;
		int		cbt_dist, cbt_conf_max;
		int		mdm_rounds, mdm_freq, mdm_div;
		int		mdm_minc, mdm_sinc;
		int		mdm_vsids_pumps, mdm_vmfq_pumps;
		int		shrink_min, sigma_min, sigma_inc;
		int		polarity, rephase_inc;
		int		luby_inc, luby_max, restart_base, restart_inc;
		int		lbd_reduce_min, lbd_csize_min;
		int		lbd_tier2, lbd_tier1, reduce_inc;
		bool	parse_only_en, report_en, model_en, proof_en;
		bool	vsids_en, mdmvsidsonly_en, vsidsonly_en, stable_en;
		bool	target_phase_en, rephase_en, guess_en, cbt_en;
		bool	reduce_en, mdmfusem_en, mdmfuses_en, mcv_en;
		//==========================================//
		//                Simplifier                //
		//==========================================//
	protected:
		VARS			*vars;
		cuMM			cumem;
		cuHist			cuhist;
		OT				*ot;
		CNF				*cnf, *hcnf;
		uint32			off1, off2;
		cudaStream_t	*streams;
		std::ofstream	outputFile;
		bool			mapped;
		int				nForced, sigState, devCount;
	public:
		//============= inline methods ==============//
		inline void		createStreams		() {
			if (streams == NULL) {
				PFLOGN2(2, " Allocating GPU streams..");
				streams = new cudaStream_t[nstreams];
				for (int i = 0; i < nstreams; i++) cudaStreamCreate(streams + i);
				PFLDONE(2, 5);
			}
		}
		inline void		destroyStreams		() {
			if (streams != NULL) {
				for (int i = 0; i < nstreams; i++) cudaStreamDestroy(streams[i]);
				delete[] streams;
			}
		}
		inline bool		verifyLCVE			() {
			for (uint32 v = 0; v < vars->numPVs; v++) if (sp->frozen[vars->pVars->at(v)]) return false;
			return true;
		}
		inline void		createWT			() {
			// construct watches based on simplified hcnf (scnf)
			assert(!hcnf->empty());
			for (uint32 i = 0; i < hcnf->size(); i++) {
				S_REF href = hcnf->ref(i);
				SCLAUSE& s = (*hcnf)[href];
				assert(!s.deleted());
				assert(s.size() > 1);
				s.set_ref(NOREF);
				WT& ws = (s.size() == 2) ? wtBin : wt;
				if (mapped) ws.newWatch(href, vmap.mapLit(s[0]), vmap.mapLit(s[1]));
				else ws.newWatch(href, s[0], s[1]);
			}

		}
		inline void		copyWatched			() {
			// copy clauses based on watches then update watches
			assert(!hcnf->empty());
			uint32 maxVar = mapped ? vmap.numVars() : inf.maxVar;
			assert(maxVar);
			for (uint32 v = 1; v <= maxVar; v++) {
				uint32 p = v2l(v), n = neg(p);
				recycleWL(wtBin[p]), recycleWL(wt[p]);
				recycleWL(wtBin[n]), recycleWL(wt[n]);
			}
		}
		inline void		copyNonWatched		() {
			assert(!hcnf->empty());
			for (uint32 i = 0; i < hcnf->size(); i++)
				newClause(hcnf->clause(i));
		}
		inline void		countMelted			() {
			SIGmA::countMelted(sp->vstate);
		}
		inline void		countFinal			() {
			SIGmA::countFinal(cnf, vars->gstats, sp->vstate);
			inf.nClauses = inf.n_cls_after;
		}
		inline void		countCls			() {
			SIGmA::countCls(cnf, vars->gstats);
		}
		inline void		countLits			() {
			SIGmA::countLits(cnf, vars->gstats);
		}
		inline void		countAll			() {
			SIGmA::countAll(cnf, vars->gstats);
		}
		inline void		evalReds			() { 
			SIGmA::evalReds(cnf, vars->gstats, sp->vstate);
		}
		inline void		logReductions		() {
			int64 varsRemoved	= int64(inf.n_del_vars_after) + nForced;
			int64 clsRemoved	= int64(inf.nClauses)	- inf.n_cls_after;
			int64 litsRemoved	= int64(inf.nLiterals)	- inf.n_lits_after;
			const char* header	= "  %-10s  %-10s %-10s %-10s";
			PFLOG1(header, " ", "Variables", "Clauses", "Literals");
			const char* rem = "  %-10s: %-9lld  %c%-8lld  %c%-8lld";
			const char* sur = "  %-10s: %-9d  %-9d  %-9d";
			PFLOG1(rem, "Removed",
				-varsRemoved,
				clsRemoved < 0 ? '+' : '-', abs(clsRemoved),
				litsRemoved < 0 ? '+' : '-', abs(litsRemoved));
			PFLOG1(sur, "Survived",
				maxActive(),
				inf.n_cls_after,
				inf.n_lits_after);
		}
		inline void		filterPVs			() {
			int idx = 0, n = vars->pVars->size();
			while (idx < n) {
				uint32& x = (*vars->pVars)[idx];
				if (ELIMINATED(x)) 
					toblivion(x), x = (*vars->pVars)[--n];
				else idx++;
			}
			vars->pVars->resize(n), vars->numPVs = n;
		}
		inline void		printOT				() {
			PFLOG0(" Positive occurs:");
			for (uint32 v = 1; v <= inf.maxVar; v++) printOL(v2l(v));
			PFLOG0(" Negative occurs:");
			for (uint32 v = 1; v <= inf.maxVar; v++) printOL(neg(v2l(v)));
		}
		inline void		printOL				(const OL& list) {
			for (uint32 i = 0; i < list.size(); i++) {
				PFLOGN0(" ");
				(*cnf)[list[i]].print();
			}
		}
		inline void		printOL				(const uint32& lit) {
			if ((*ot)[lit].size() == 0) return;
			PFLOG1(" List(%d):", l2i(lit));
			printOL((*ot)[lit]);
		}
		inline void		toblivion			(const uint32 x) {
			sp->vstate[RECOVERVAR(x)] = MELTED;
		}
		inline bool		reallocCNF			(const int& times) {
			if (times > 1 && times != phases && (times % shrink_rate) == 0) {
				inf.maxAddedCls = inf.nClauses, inf.maxAddedLits = inf.nLiterals;
				PFLOG2(2, " Maximum added clauses/literals = %d/%d", inf.maxAddedCls, inf.maxAddedLits);
				if (!cumem.resizeCNF(cnf, inf.nClauses + inf.maxAddedCls, inf.nLiterals + inf.maxAddedLits)) return false;
			}
			else cumem.cacheCNFPtr(cnf);
			return true;
		}
		inline bool		reallocOT			(const cudaStream_t& stream = 0) {
			assert(inf.nLiterals);
			calcOccurs(inf.nLiterals);
			if (!cumem.resizeOTAsync(ot, inf.nLiterals + inf.maxAddedLits, stream)) { sigState = OTALLOC_FAIL; return false; }
			return true;
		}
		inline void		reflectCNF			(const cudaStream_t& s1, const cudaStream_t& s2) {
			uint32 len1 = hcnf->data().size - off1;
			if (!len1) return;
			uint32 len2 = hcnf->size() - off2;
			CHECK(cudaMemcpyAsync(cumem.cnfMem() + off1, hcnf->data().mem + off1, len1 * sizeof(S_REF), cudaMemcpyHostToDevice, s1));
			CHECK(cudaMemcpyAsync(cumem.csMem() + off2, hcnf->csData() + off2, len2 * sizeof(S_REF), cudaMemcpyHostToDevice, s2));
			off1 = hcnf->data().size, off2 = hcnf->size();
		}
		inline void		cacheUnits			(const cudaStream_t& stream) {
			sync(stream);
			if (vars->nUnits = vars->tmpObj.size()) CHECK(cudaMemcpyAsync(vars->cachedUnits, cumem.unitsdPtr(),
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
			uint32 *start = model.resolved + off;
			CHECK(cudaMemcpyAsync(start, devStart, devSize * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
			if (sync_always || unified_access) sync(stream);
		}
		//===========================================//
		void			optSimp				();
		void			varReorder			();
		void			awaken				();
		void			newBeginning		();
		void			sigmify				();
		bool			LCVE				();
		bool			prop				();
		void			VE					();
		void			SUB					();
		void			HRE					();
		void			BCE					();
		void			masterFree			();
		void			slavesFree			();
		void			histSimp			(const uint32&);
		void			calcOccurs			(const uint32&);
		void			extract				(CNF*, WT&);
		void			extract				(CNF*, BCNF&);
		bool			propClause			(SCLAUSE&, const uint32&);
		C_REF			newClause			(SCLAUSE&);
		inline void		depFreeze			(const OL&, const uint32&, const uint32&, const uint32&);
		inline void		cacheCNF			(const cudaStream_t&, const cudaStream_t&);
		inline bool		enqeueCached		(const cudaStream_t&);
		inline void		cleanProped			();
		inline void		cleanFixed			();
		inline void		cleanDynamic		();
		inline void		initSimp			();
		//================ options ================//
		bool	sigma_en, sigma_live_en, solve_en, sort_cnf_en;
		bool	ve_en, ve_plus_en, sub_en, bce_en, hre_en, cls_en, all_en;
		int		phases, shrink_rate, xor_limit, ngpus, nstreams;
		uint32	lcve_min, ve_phase_min, mu_pos, mu_neg;
		uint32	hse_limit, bce_limit, hre_limit;
	};
	extern ParaFROST* pfrost;
}

#endif 