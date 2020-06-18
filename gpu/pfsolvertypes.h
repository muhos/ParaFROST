/***********************************************************************[pfsolvertypes.h]
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

#ifndef __SOLVER_TYPES_
#define __SOLVER_TYPES_

#include "pfclause.h"
#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {
	/*****************************************************/
	/*  Usage:   store clause index per literal          */
	/*  Dependency:  none                                */
	/*****************************************************/
	struct WATCH {
		C_REF	ref;
		uint32	imp;
		inline ~WATCH		() { ref = NOREF; }
		inline WATCH		() { ref = NOREF; imp = 0; }
		inline WATCH		(const C_REF& cref, const uint32& lit) { ref = cref; imp = lit; }
	};
	/*****************************************************/
	/*  Global vector types for CNF/Watch lists          */
	/*****************************************************/
	typedef Vec<WATCH, int> WL;
	typedef Vec<C_REF> BCNF;
	/*****************************************************/
	/*  Usage:   store watch lists per literal           */
	/*  Dependency:  Watch                               */
	/*****************************************************/
	class WT {
		Vec<WL> wt;
		Vec<bool> collected;
		Vec<uint32> garbage;
		const CMM& cm;
	public:
		inline				WT			(const CMM& _cm) : cm(_cm) {}
		inline void			resize		(const uint32& size) { assert(size); wt.resize(size), collected.resize(size, 0); }
		inline WL&			getClean	(const uint32& lit) { assert(lit > 1); if (collected[lit]) recycle(lit); return wt[lit]; }
		inline WL&			operator[]	(const uint32& lit) { assert(lit > 1); return wt[lit]; }
		inline const WL&	operator[]	(const uint32& lit) const { assert(lit > 1); return wt[lit]; }
		inline void			destroy		(const uint32& lit) { wt[lit].clear(true); }
		inline void			collect		(const uint32& lit) { assert(lit > 1); if (!collected[lit]) { collected[lit] = 1; garbage.push(lit); } }
		inline void			recycle		(const uint32& lit) {
			WL& ws = wt[lit];
			if (ws.empty()) { ws.clear(true); collected[lit] = 0; return; }
			WATCH* w_i, * w_j, * end = ws.end();
			for (w_i = ws, w_j = ws; w_i != end; w_i++) {
				assert(w_i->ref != NOREF);
				if (!cm[w_i->ref].deleted()) *w_j++ = *w_i;
			}
			ws.shrink(int(w_i - w_j));
			if (ws.empty()) ws.clear(true);
			collected[lit] = 0;
		}
		inline void			newWatch	(const C_REF& r, const uint32& l0, const uint32& l1) {
			wt[flip(l0)].push(WATCH(r, l1));
			wt[flip(l1)].push(WATCH(r, l0));
		}
		inline void			remWatch	(const uint32& lit, const C_REF cref) {
			WL& ws = wt[lit];
			if (ws.size() == 0) return;
			if (ws.size() > 1) {
				int c_idx = 0;
				while (ws[c_idx].ref != cref) c_idx++;
				assert(c_idx < ws.size());
				while (c_idx < ws.size() - 1) { ws[c_idx] = ws[c_idx + 1]; c_idx++; }
			}
			ws.pop();
		}
		inline void			print		(const uint32& lit) {
			WL& ws = wt[lit];
			if (ws.size()) PFLOG1(" list(%d):", l2i(lit));
			for (int i = 0; i < ws.size(); i++) {
				PFLCLAUSE(1, cm[ws[i].ref], " W(r: %-4d, b: %-4d)->", ws[i].ref, l2i(ws[i].imp));
			}
		}
		inline void			clear		(bool _free = true) {
			wt.clear(_free);
			collected.clear(_free);
			garbage.clear(_free);
		}
		inline void			print		() {
			PFLOG0(" Watches:");
			for (uint32 lit = 2; lit < wt.size(); lit++) 
				print(lit);
		}
		inline void			recycle		() {
			for (uint32 i = 0; i < garbage.size(); i++) if (collected[garbage[i]]) recycle(garbage[i]);
			garbage.clear();
		}
		inline void			shrinkCap	() { wt.shrinkCap(), collected.shrinkCap(); }
		inline int64		size		() const { return wt.size(); }
		inline bool			empty		() const { return wt.size() == 0; }

	};
	/*****************************************************/
	/*  Usage:    Information of search space            */
	/*  Dependency: none                                 */
	/*****************************************************/
	class SP {
		addr_t		_mem;
		size_t		_sz, _cap;
		template <class T>
		size_t		calcBytes		(const uint32& sz, const uint32& nVecs) { 
			assert(sz); return size_t(sz) * nVecs * sizeof(T);
		}
	public:
		int			*level, *board, bt_level, cbt_level;
		LIT_ST		*value, *locked, *seen, *frozen, *vstate, *psaved, *ptarget, *pbest;
		uint32		*tmp_stack, propagated, simplified, learnt_lbd;
		C_REF		*source;
		bool		max1Found;
					SP				(const uint32& size) :
						_mem(NULL)
						, _sz(size)
						, _cap(0)
						, max1Found(false)
						, propagated(0)
						, simplified(0)
						, learnt_lbd(0)
						, bt_level(ROOT_LEVEL)
						, cbt_level(UNDEFINED) {
						size_t vec1Bytes = calcBytes<LIT_ST>(size, 9);
						size_t vec4Bytes = calcBytes<uint32>(size, 4);
						_cap = vec1Bytes + vec4Bytes;
						assert(_cap);
						pfalloc(_mem, _cap);
						assert(_mem != NULL);
						memset(_mem, 0, _cap);
						// 1-byte arrays
						value = (LIT_ST*)_mem;
						locked = value + _sz + _sz;
						frozen = locked + _sz;
						seen = frozen + _sz;
						psaved = seen + _sz;
						ptarget = psaved + _sz;
						pbest = ptarget + _sz;
						vstate = pbest + _sz;
						// 4-byte arrays
						level = (int*)(_mem + vec1Bytes);
						board = level + _sz;
						tmp_stack = (uint32*)(board + _sz);
						source = tmp_stack + _sz;
					}
		void		printStates		(const uint32& size) {
						PFLOGN1(" States->[");
						for (uint32 v = 1; v <= size; v++) {
							fprintf(stdout, "%5d:%d ", v, vstate[v]);
							if (v > 1 && v < size - 1 && v % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
						}
						putc(']', stdout), putc('\n', stdout);
					}
		void		printValues		(const uint32& size) {
			PFLOGN1(" Values->[");
			for (uint32 lit = 2; lit <= v2l(size); lit++) {
				fprintf(stdout, "%5d:%d ", l2i(lit), value[lit]);
				if (lit > 2 && lit < size - 1 && lit % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		printLevels		(const uint32& size) {
			PFLOGN1(" Levels->[");
			for (uint32 v = 1; v <= size; v++) {
				fprintf(stdout, "%5d@%d ", v, level[v]);
				if (v > 1 && v < size - 1 && v % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		printSource		(const uint32& size) {
			PFLOGN1(" Sources->[");
			int noref = -1;
			for (uint32 v = 1; v <= size; v++) {		
				fprintf(stdout, "%5d:%-5d ", v, source[v] == NOREF ? noref : source[v]);
				if (v > 1 && v < size - 1 && v % 10 == 0) { putc('\n', stdout); PFLOGN0("\t "); }
			}
			putc(']', stdout), putc('\n', stdout);
		}
		void		lockMelted		() {
			for (uint32 v = 1; v <= _sz; v++)
				if (vstate[v] == MELTED)
					locked[v] = 1;
		}
		void		reset			() { 
			bt_level = ROOT_LEVEL, cbt_level = UNDEFINED, max1Found = false; 
		}
		void		destroy			() { if (_mem != NULL) std::free(_mem); }
					~SP				() { destroy(); }
	};
	struct LEARN {
		int64 simpProps, reductions, bumped;
		int64 mdm_conf_max, reduce_conf_max;
		int64 restarts_conf_max, stable_conf_max;
		int64 rephased[2], rephase_conf_max, rephase_last_max;
		double var_inc, var_decay;
		float cl_inc, cl_decay;
		uint32 numMDs, nRefVars;
		uint32 best, target;
		int	rounds;
		LIT_ST lastrephased;
		bool stable;
	};
	struct MODEL {
		Vec<LIT_ST> value;
		uVec1D lits, resolved;
		uint32 maxVar;
					MODEL			() : maxVar(0) {}
		void		init			() {
			assert(inf.maxVar);
			PFLOG2(2, " Initially mapping original variables to literals..");
			maxVar = inf.maxVar;
			lits.resize(maxVar + 1), lits[0] = 0;
			for (uint32 v = 1; v <= inf.maxVar; v++) lits[v] = v2l(v);
		}
		void		print			() {
			PFLMH('v');
			for (uint32 v = 1; v <= maxVar; v++)
				PFLMLIT(v, value[v]);
			putc('\n', stdout);
			if (!quiet_en) PFLOG0("");
		}
		void		extend			(LIT_ST*);
		__forceinline
		bool		satisfied		(const uint32& lit) { return value[l2a(lit)] == !sign(lit); }
	};
	struct STATS {
		int64 sysMemAvail;
		int64 n_rephs;
		int64 n_fuds, n_pds;
		int64 n_units, n_props;
		int64 tot_lits, max_lits;
		int marker, mdm_calls;
		int mappings, shrinkages;
		int reduces, recyclings;
		int stab_restarts, ncbt, cbt;
		bool guess_succ;
		const char* guess_who;
		void reset() {
			marker = 0;
			n_rephs = 0;
			reduces = 0;
			mappings = 0;
			mdm_calls = 0;
			shrinkages = 0;
			recyclings = 0;
			ncbt = cbt = 0;
			stab_restarts = 0;
			n_pds = n_fuds = 0;
			n_props = n_units = 0;
			max_lits = tot_lits = 0;
			guess_succ = false;
			guess_who = "none succeeded";
		}
	};
	struct HEAP_CMP {
		const Vec<double>& act;
		HEAP_CMP(const Vec<double>& _act) : act(_act) {}
		bool operator () (uint32& a, uint32& b) { return act[a] < act[b]; }
	};
	struct LEARNT_CMP {
		const CMM& cm;
		LEARNT_CMP(const CMM& _cm) : cm(_cm) {}
		bool operator () (C_REF& a, C_REF& b) {
			if (cm[a].LBD() != cm[b].LBD()) return cm[a].LBD() > cm[b].LBD();
			else return cm[a].activity() != cm[b].activity() ? cm[a].activity() < cm[b].activity() : cm[a].size() > cm[b].size();
		}
	};
	struct ANALYZE_CMP {
		Vec<int64>& bumped;
		ANALYZE_CMP(Vec<int64>& _bumped) : bumped(_bumped) {}
		bool operator () (uint32& a, uint32& b) { return bumped[a] < bumped[b]; }
	};
}

#endif