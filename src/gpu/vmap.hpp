/***********************************************************************[vmap.h]
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

#ifndef __VMAP_
#define __VMAP_

#include "solvertypes.hpp"

namespace ParaFROST {

	class VMAP {
		SP* sp;
		uVec1D _mapped;
		uint32 newVars, firstDL0, mappedFirstDL0;
		LIT_ST valFirstDL0;
		inline void			map					(const uint32& old) { 
			CHECKVAR(old); 
			_mapped[old] = ++newVars;
		}
	public:
							~VMAP				() { destroy(); }
							VMAP				() : sp(NULL), newVars(0), firstDL0(0), mappedFirstDL0(0), valFirstDL0(UNDEFINED) {}
		inline uint32*		operator*			() { return _mapped; }
		inline bool			empty				() const { return !newVars; }
		inline uint32		size				() const { return newVars + 1; }
		inline uint32		numVars				() const { return newVars; }
		inline uint32		firstL0				() const { return firstDL0; }
		inline uint32		mapped				(const uint32& old) const { return _mapped[old]; }
		inline uint32		mapLit				(const uint32& lit) {
			assert(!_mapped.empty());
			CHECKLIT(lit);
			uint32 oldVar = ABS(lit), newVar = mapped(oldVar);
			assert(newVar <= newVars);
			if (newVar && NEQUAL(oldVar, firstDL0)) {
				assert(UNASSIGNED(sp->value[lit]));
				assert(!sp->vstate[oldVar].state);
				return V2DEC(newVar, SIGN(lit));
			}
			LIT_ST val = sp->value[lit];
			if (UNASSIGNED(val)) return 0;
			assert(val >= 0);
			uint32 newLitDL0 = V2L(mappedFirstDL0);
			if (valFirstDL0 != val) newLitDL0 = FLIP(newLitDL0);
			return newLitDL0;
		}
		inline void			initiate			(SP* _sp) {
			assert(inf.maxVar);
			assert(_sp != NULL);
			sp = _sp;
			_mapped.resize(inf.maxVar + 1, 0);
			forall_variables(old) {
				if (!sp->vstate[old].state) map(old);
				else if (!firstDL0 && !UNASSIGNED(sp->value[V2L(old)])) {
					firstDL0 = old, valFirstDL0 = sp->value[V2L(firstDL0)];
					map(firstDL0), mappedFirstDL0 = newVars;
				}
			}
			assert(newVars <= inf.maxVar);
			PFLOG2(2, " Mapped %d to %d, first frozen/autartic literal \"%d\"", inf.maxVar, newVars,
				firstDL0 ? (valFirstDL0 ? firstDL0 : -int(firstDL0)) : 0);
		}
		inline void			mapTransitive		(uint32& lit) {
			if (lit <= 2) return;
			CHECKLIT(lit);
			uint32 v = ABS(lit);
			uint32 mlit = 0;
			if (sp->vstate[v].state) {
				while (v <= inf.maxVar && sp->vstate[v].state) v++;
				if (v <= inf.maxVar) mlit = mapLit(V2L(v));
			}
			else mlit = mapLit(lit);
			if (!mlit) mlit = 2;
			lit = mlit;
		}
		inline void			mapSP				(SP* to) {
			// map all arrays
			forall_variables(v) {
				const uint32 mVar = mapped(v);
				if (mVar) {
					const uint32 p = V2L(v), n = NEG(p);
					const uint32 mpos = V2L(mVar), mneg = NEG(mpos);
					// map 'value'
					to->value[mpos] = sp->value[p];
					to->value[mneg] = sp->value[n];
					// map others
					to->source[mVar] = sp->source[v];
					to->level[mVar] = sp->level[v];
					to->board[mVar] = sp->board[v];
					to->pbest[mVar] = sp->pbest[v];
					to->psaved[mVar] = sp->psaved[v];
					to->ptarget[mVar] = sp->ptarget[v];
					to->marks[mVar] = sp->marks[v];
					to->vstate[mVar] = sp->vstate[v];
				}
			}
			// update counters
			to->propagated = sp->propagated;
			to->simplified = firstDL0 ? 1 : 0;
			sp = NULL; // nullify local reference
		}
		inline void			mapOrgs             (Lits_t& lits) {
			forall_vector(uint32, lits, i) {
				const uint32 lit = *i;
				if (lit) {
					CHECKLIT(lit);
					*i = mapLit(lit);
					PFLOG2(4, " Literal %d mapped to %d", l2i(lit), *i ? l2i(*i) : 0);
				}
			}
		}
		inline void			mapOrgs             (uVec1D& lits) {
			forall_vector(uint32, lits, i) {
				const uint32 lit = *i;
				if (lit) {
					CHECKLIT(lit);
					*i = mapLit(lit);
					PFLOG2(4, " Literal %d mapped to %d", l2i(lit), *i ? l2i(*i) : 0);
				}
			}
		}
		template <class T>
		inline void			mapShrinkVars		(Vec<T>& vars) {
			assert(inf.maxVar >= newVars);
			forall_variables(v) {
				uint32 mVar = mapped(v);
				if (mVar) vars[mVar] = vars[v];
			}
			vars.resize(size());
			vars.shrinkCap();
		}
		template <class T>
		inline void			mapShrinkDualVars	(Vec<T>& vars) {
			assert(inf.maxVar >= newVars);
			forall_variables(v) {
				uint32 mVar = mapped(v);
				if (mVar) {
					uint32 p = V2L(v), n = NEG(p);
					uint32 mpos = V2L(mVar), mneg = NEG(mpos);
					vars[mpos] = vars[p];
					vars[mneg] = vars[n];
				}
			}
			vars.resize(V2L(size()));
			vars.shrinkCap();
		}
		template <class T>
		inline void			mapShrinkLits		(Vec<T>& lits) {
			uint32 *d = lits;
			forall_vector(uint32, lits, s) {
				const uint32 srclit = *s;
				CHECKLIT(srclit);
				uint32 mVar = mapped(ABS(srclit));
				if (mVar) *d++ = V2DEC(mVar, SIGN(srclit));
			}
			lits.resize(uint32(d - lits));
			lits.shrinkCap();
		}
		template <class SRC, class DEST>
		inline void			mapClause			(DEST& dest, SRC& src) {
			assert(src.size() > 1);
			assert(!src.deleted());
			PFLCLAUSE(4, src, " Clause    ");
			for (int i = 0; i < src.size(); i++) {
				CHECKLIT(src[i]);
				assert(UNASSIGNED(sp->value[src[i]]));
				dest[i] = mapLit(src[i]);
			}
			PFLCLAUSE(4, dest, " mapped to ");
		}
		template <class C>
		inline void			mapClause			(C& c) {
			assert(c.size() > 1);
			assert(!c.deleted());
			assert(!c.moved());
			PFLCLAUSE(4, c, " Clause    ");
			forall_clause(c, i) {
				CHECKLIT(*i);
				assert(UNASSIGNED(sp->value[*i]));
				*i = mapLit(*i);
			}
			PFLCLAUSE(4, c, " mapped to ");
		}
		inline void			destroy				() {
			sp = NULL, _mapped.clear(true);
			newVars = 0, firstDL0 = mappedFirstDL0 = 0;
			valFirstDL0 = UNDEFINED;
		}
	};

}

#endif