/***********************************************************************[pfvmap.h]
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

#include "pfsolvertypes.h"

namespace pFROST {

	class VMAP {
		SP* sp;
		uVec1D _mapped;
		uint32 newVars, firstDL0, mappedFirstDL0;
		LIT_ST valFirstDL0;
	public:
							~VMAP			() { destroy(); }
							VMAP			() : sp(NULL), newVars(0), firstDL0(0), mappedFirstDL0(0), valFirstDL0(UNDEFINED) {}
		inline uint32*		operator*		() { return _mapped; }
		inline bool			empty			() const { return !newVars; }
		inline uint32		size			() const { return newVars + 1; }
		inline uint32		numVars			() const { return newVars; }
		inline uint32		firstL0			() const { return firstDL0; }
		inline uint32		mapped			(const uint32& old) const {  return _mapped[old]; }
		inline void			map				(const uint32& old) { assert(old && old <= inf.maxVar); _mapped[old] = ++newVars; }
		inline uint32		mapLit			(const uint32& lit) {
			assert(!_mapped.empty());
			assert(lit > 1);
			uint32 oldVar = ABS(lit), newVar = mapped(oldVar);
			assert(newVar <= newVars);
			if (newVar && oldVar != firstDL0) {
				assert(UNASSIGNED(sp->value[lit]));
				assert(sp->vstate[oldVar] == ACTIVE);
				return V2DEC(newVar, SIGN(lit));
			}
			LIT_ST val = sp->value[lit];
			if (UNASSIGNED(val)) return 0;
			assert(val >= 0);
			assert(sp->vstate[oldVar] == FROZEN);
			uint32 newLitDL0 = V2L(mappedFirstDL0);
			if (valFirstDL0 != val) newLitDL0 = FLIP(newLitDL0);
			return newLitDL0;
		}
		inline void			initiate		(SP* _sp) {
			assert(inf.maxVar);
			assert(_sp != NULL);
			sp = _sp;
			uint32 oldVars = inf.maxVar;
			_mapped.resize(oldVars + 1, 0);
			for (uint32 old = 1; old <= oldVars; old++) {
				if (sp->vstate[old] == ACTIVE) map(old);
				else if (sp->vstate[old] == FROZEN && !firstDL0) {
					firstDL0 = old, valFirstDL0 = sp->value[V2L(firstDL0)];
					map(firstDL0), mappedFirstDL0 = newVars;
				}
			}
			assert(newVars <= inf.maxVar);
			PFLOG2(2, " Mapped %d to %d, first frozen literal \"%d\"", oldVars, newVars,
				firstDL0 ? (valFirstDL0 ? firstDL0 : -int(firstDL0)) : 0);
		}
		inline void			mapSP			(SP* to) {
			// map values
			for (uint32 v = 1; v <= inf.maxVar; v++) {
				uint32 mVar = mapped(v);
				if (mVar) {
					uint32 p = V2L(v), n = NEG(p);
					uint32 mpos = V2L(mVar), mneg = NEG(mpos);
					to->value[mpos] = sp->value[p];
					to->value[mneg] = sp->value[n];
				}
			}
			// map lock, variable state, subsume, and phases
			mapVars(to->locked	,	sp->locked);
			mapVars(to->vstate	,	sp->vstate);
			mapVars(to->subsume ,	sp->subsume);
			mapVars(to->pbest	,	sp->pbest);
			mapVars(to->psaved	,	sp->psaved);
			mapVars(to->ptarget	,	sp->ptarget);
			// map source and level
			mapVars(to->source	,	sp->source);
			mapVars(to->level	,	sp->level);
			// update counters
			to->propagated = sp->propagated;
			to->simplified = firstDL0 ? 1 : 0;
			sp = NULL; // nullify local reference
		}
		inline void			mapOrgs			(uVec1D& lits) {
			uint32* end = lits.end();
			for (uint32* i = lits; i != end; i++) {
				uint32 lit = *i;
				if (lit) {
					assert(lit > 1);
					*i = mapLit(lit);
					PFLOG2(4, " Literal %d mapped to %d", l2i(lit), *i ? l2i(*i) : 0);
				}
			}
		}
		template <class T>
		inline void			mapVars			(Vec<T>& vars) {
			assert(inf.maxVar - newVars >= 1);
			for (uint32 v = 1; v <= inf.maxVar; v++) {
				uint32 mVar = mapped(v);
				if (mVar) vars[mVar] = vars[v];
			}
		}
		template <class T>
		inline void			mapShrinkVars	(Vec<T>& vars) {
			assert(inf.maxVar >= newVars);
			for (uint32 v = 1; v <= inf.maxVar; v++) {
				uint32 mVar = mapped(v);
				if (mVar) vars[mVar] = vars[v];
			}
			vars.resize(size());
			vars.shrinkCap();
		}
		template <class T>
		inline void			mapShrinkLits	(Vec<T>& lits) {
			uint32 *s = lits, *end = lits.end(), *d = s;
			while (s != end) {
				assert(*s > 1);
				uint32 mVar = mapped(ABS(*s));
				if (mVar) *d++ = V2DEC(mVar, SIGN(*s));
				s++;
			}
			lits.resize(uint32(d - lits));
			lits.shrinkCap();
		}
		template <class SRC, class DEST>
		inline void			mapVars			(DEST& dest, SRC& src) {
			assert(inf.maxVar >= newVars);
			for (uint32 v = 1; v <= inf.maxVar; v++) {
				uint32 mVar = mapped(v);
				if (mVar) dest[mVar] = src[v];
			}
		}
		template <class SRC, class DEST>
		inline void			mapClause		(DEST& dest, SRC& src) {
			assert(src.size() > 1);
			assert(!src.deleted());
			PFLCLAUSE(4, src, " Clause    ");
			for (int i = 0; i < src.size(); i++) {
				assert(src[i] > 1);
				assert(UNASSIGNED(sp->value[src[i]]));
				dest[i] = mapLit(src[i]);
			}
			PFLCLAUSE(4, dest, " mapped to ");
		}
		template <class C>
		inline void			mapClause		(C& c) {
			assert(c.size() > 1);
			assert(!c.deleted());
			assert(!c.moved());
			PFLCLAUSE(4, c, " Clause    ");
			uint32* end = c.end();
			for (uint32 *i = c; i != end; i++) {
				assert(*i > 1);
				assert(UNASSIGNED(sp->value[*i]));
				*i = mapLit(*i);
			}
			PFLCLAUSE(4, c, " mapped to ");
		}
		inline void			destroy			() {
			sp = NULL, _mapped.clear(true);
			newVars = 0, firstDL0 = mappedFirstDL0 = 0;
			valFirstDL0 = UNDEFINED;
		}
	};

}

#endif