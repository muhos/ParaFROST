/***********************************************************************[pfbve.cuh]
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

#ifndef __SIGMA_BVE_
#define __SIGMA_BVE_

#include "pfdevice.cuh"
#include "pfgates.cuh"

namespace pFROST {

	namespace SIGmA {

#define VE_DBG 0 // set to serialize BVE

		_PFROST_D_ void toblivion(const uint32& dx, const uint32& nSaved, CNF& cnf, OL& toSave, OL& other, cuVecU* resolved)
		{
			uint32* saved = resolved->jump(nSaved);
#pragma unroll
			for (S_REF* i = toSave; i != toSave.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveResolved(saved, c);
				c.markDeleted();
			}
			saveResolved(saved, dx);
#pragma unroll
			for (S_REF* i = other; i != other.end(); i++) cnf[*i].markDeleted();
			toSave.clear(true), other.clear(true);
		}

		_PFROST_D_ void saveResolved(const uint32& dx, const uint32& nSaved, CNF& cnf, OL& toSave, cuVecU* resolved)
		{
			uint32* saved = resolved->jump(nSaved);
#pragma unroll
			for (S_REF* i = toSave; i != toSave.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) saveResolved(saved, c);
			}
			saveResolved(saved, dx);
		}

		_PFROST_D_ void calcResolvents(const uint32& x, CNF& cnf, OL& me, OL& other, const uint32& pOrgs, const uint32& nOrgs, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(x);
			assert(pOrgs);
			assert(nOrgs);
			assert(!nAddedLits);
			uint32 nTs = 0;
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (isTautology(x, ci, cj)) nTs++;
					else nAddedLits += ci.size() + cj.size() - 2;
				}
			}
			assert(pOrgs * nOrgs >= nTs);
			nAddedCls = pOrgs * nOrgs - nTs;
		}

		_PFROST_D_ void resolve_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			int nbuckets = dc_nbuckets;
			uint32 checksum = 0;
			S_REF ref;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && checksum < nAddedCls; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					int max_sz = ci.size() + cj.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (!isTautology(x, ci, cj)) {
							SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
							merge(x, ci, cj, added), max_sz = added->size();
							added->markAdded();
							if (max_sz == 1) units->push(**added);
							cs[checksum++] = ref, ref += (max_sz - 1) + nbuckets;
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							added->markAdded();
							if (merged_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (merged_sz - 1) + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
				}
			}
			assert(checksum == nAddedCls);
		}

		_PFROST_D_ void substitute_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			int nbuckets = dc_nbuckets;
			uint32 checksum = 0;
			S_REF ref;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && checksum < nAddedCls; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() == cj.molten()) continue;
					int max_sz = ci.size() + cj.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (!isTautology(x, ci, cj)) {
							SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
							merge(x, ci, cj, added), max_sz = added->size();
							added->markAdded();
							if (max_sz == 1) units->push(**added);
							cs[checksum++] = ref, ref += (max_sz - 1) + nbuckets;
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							added->markAdded();
							if (merged_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (merged_sz - 1) + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
				}
			}
			assert(checksum == nAddedCls);
		}

		_PFROST_D_ bool resolve(const uint32& x, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, cuVecU* resolved, uint32* out_c, const bool& bound)
		{
			assert(x);
			assert(checkMolten(cnf, poss, negs));
			uint32 p = V2L(x);
			uint32 nAddedCls = 0, nAddedLits = 0, psLits = 0, nsLits = 0;
			// count clauses/literals before and after
			calcResolvents(x, cnf, poss, negs, pOrgs, nOrgs, nAddedCls, nAddedLits);
			if (nAddedCls) {
				if (nAddedCls > pOrgs + nOrgs) return false;
				countLitsBefore(cnf, poss, psLits), countLitsBefore(cnf, negs, nsLits);
				if (bound && nAddedLits > psLits + nsLits) return false;
			}
#if VE_DBG
			printf("c | Resolving(%d) ==> added = %d, deleted = %d\n", x, nAddedCls, poss.size() + negs.size());
			pClauseSet(cnf, poss, negs);
#endif
			// resolve
			if (nAddedCls) resolve_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, out_c);
			// save resolved
			uint32 nSaved, lit;
			bool which = pOrgs > nOrgs;
			if (which) nSaved = nOrgs + nsLits + 2, lit = p;
			else nSaved = pOrgs + psLits + 2, lit = NEG(p);
			toblivion(lit, nSaved, cnf, which ? negs : poss, which ? poss : negs, resolved);
			return true; // resolution successful
		}

		_PFROST_D_ bool find_equ_gate(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, cuVecU* resolved)
		{
			assert(p > 1);
			assert(!SIGN(p));
			assert(checkMolten(cnf, poss, negs));
			assert(pOrgs && nOrgs);
			uint32 first;
			if (first = find_sfanin(p, cnf, poss)) {
				uint32 second = NEG(p), def = first;
				if (second < def) first = second, second = def;
				assert(first < second);
				for (S_REF* i = negs; i != negs.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original() && c.size() == 2 && c[0] == first && c[1] == second) {
#if VE_DBG
						printf("c | Gate %d = -/+%d found\n", ABS(p), ABS(def));
						pClauseSet(cnf, poss, negs);
#endif
						saveResolved(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
						substitute_single(p, def, cnf, poss, negs, units);
						return true;
					}
				}
			}
			return false;
		}

		_PFROST_D_ bool find_ao_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, cuVecU* units, cuVecU* resolved, uint32* out_c)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			uint32 sig, x = ABS(dx);
			// (-) ==> look for AND , (+) ==> look for OR
#if VE_DBG
			const char* type = SIGN(dx) ? "AND" : "OR";
#endif
			OL& itarget = ot[dx];
			int nImps = find_fanin(dx, cnf, itarget, out_c, sig);
			if (nImps > 1) {
				uint32 f_dx = FLIP(dx);
				out_c[nImps++] = f_dx;
				sig |= MAPHASH(f_dx);
				devSort(out_c, nImps);
				OL& otarget = ot[f_dx];
				for (S_REF* i = otarget; i != otarget.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.learnt()) continue;
					if (c.size() == nImps && sub(c.sig(), sig) && isEqual(c, out_c, nImps)) {
						c.melt(); // mark as fanout clause
						// check resolvability
						uint32 nAddedCls = 0, nAddedLits = 0;
						countSubstituted(x, cnf, itarget, otarget, nAddedCls, nAddedLits);
						if (nAddedCls > nOrgCls) { c.freeze(); break; }
						// can be substituted
#if VE_DBG
						printf("c | Gate %d = %s(", ABS(dx), type);
						for (int k = 0; k < nImps; k++) {
							if (ABS(out_c[k]) == ABS(dx)) continue;
							printf(" %d", ABS(out_c[k]));
							if (k < nImps - 1) printf(",");
						}
						printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
						pClauseSet(cnf, itarget, otarget);
#endif
						if (nAddedCls) substitute_x(x, nAddedCls, nAddedLits, cnf, itarget, otarget, units, out_c);
						return true;
					}
				}
			}
			freeze_binaries(cnf, itarget);
			return false;
		}

		_PFROST_D_ bool find_ite_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, cuVecU* units, cuVecU* resolved, uint32* out_c)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			OL& itarget = ot[dx];
			for (S_REF* i = itarget; i != itarget.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt() || ci.size() < 3 || ci.size() > 3) continue;
				uint32 xi = ci[0], yi = ci[1], zi = ci[2];
				if (yi == dx) devSwap(xi, yi);
				if (zi == dx) devSwap(xi, zi);
				assert(xi == dx);
				for (S_REF* j = i + 1; j != itarget.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt() || cj.size() < 3 || cj.size() > 3) continue;
					uint32 xj = cj[0], yj = cj[1], zj = cj[2];
					if (yj == dx) devSwap(xj, yj);
					if (zj == dx) devSwap(xj, zj);
					assert(xj == dx);
					if (ABS(yi) == ABS(zj)) devSwap(yj, zj);
					if (ABS(zi) == ABS(zj)) continue;
					if (yi != FLIP(yj)) continue;
					uint32 f_dx = FLIP(dx);
					S_REF r1 = fast_equality_check(cnf, ot, f_dx, yi, FLIP(zi));
					if (r1 == GNOREF) continue;
					S_REF r2 = fast_equality_check(cnf, ot, f_dx, yj, FLIP(zj));
					if (r2 == GNOREF) continue;
					assert(cnf[r1].original());
					assert(cnf[r2].original());
					// mark gate clauses
					ci.melt(), cj.melt();
					cnf[r1].melt(), cnf[r2].melt();
					// check resolvability
					uint32 v = ABS(dx);
					uint32 nAddedCls = 0, nAddedLits = 0;
					OL& otarget = ot[f_dx];
					countSubstituted(v, cnf, itarget, otarget, nAddedCls, nAddedLits);
					if (nAddedCls > nOrgCls) {
						ci.freeze(), cj.freeze();
						cnf[r1].freeze(), cnf[r2].freeze();
						return false;
					}
					// can be substituted
#if VE_DBG
					printf("c | Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d", ABS(dx), ABS(yi), ABS(zi), ABS(zj), nAddedCls, itarget.size() + otarget.size());
					pClauseSet(cnf, itarget, otarget);
#endif
					if (nAddedCls) substitute_x(v, nAddedCls, nAddedLits, cnf, itarget, otarget, units, out_c);
					return true;
				}
			}
			return false;
		}

		_PFROST_D_ bool find_xor_gate(const uint32& dx, const uint32& nOrgCls, const int& xor_max_arity, CNF& cnf, OT& ot, cuVecU* units, cuVecU* resolved, uint32* out_c)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			OL& itarget = ot[dx];
			for (S_REF* i = itarget; i != itarget.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				int size = ci.size();
				int arity = size - 1; // XOR arity
				if (ci.learnt() || size < 3 || arity > xor_max_arity || arity > SH_MAX_BVE_OUT) continue;
				// share to out_c except dx
				shareXORClause(dx, ci, out_c);
				// find arity clauses
				int itargets = 0;
				for (int j = 0; j < arity; j++) {
					S_REF r = find_fanin(dx, j, cnf, itarget, out_c, arity);
					if (r == GNOREF) break;
					cnf[r].melt(), itargets++;
				}
				if (itargets < arity) {
					freeze_arities(cnf, itarget);
					continue;
				}
				// find all +/-  
				uint32 f_dx = FLIP(dx);
				S_REF r1 = find_all(f_dx, cnf, ot, out_c, arity);
				if (r1 == GNOREF) break;
				flip_all(out_c, arity);
				S_REF r2 = find_all(f_dx, cnf, ot, out_c, arity);
				if (r2 == GNOREF) break;
				assert(cnf[r1].original());
				assert(cnf[r2].original());
				cnf[r1].melt(), cnf[r2].melt();
				// check resolvability
				uint32 v = ABS(dx);
				OL& otarget = ot[f_dx];
				uint32 nAddedCls = 0, nAddedLits = 0;
				countSubstituted(v, cnf, itarget, otarget, nAddedCls, nAddedLits);
				if (nAddedCls > nOrgCls) {
					cnf[r1].freeze(), cnf[r2].freeze();
					break;
				}
				// can be substituted
#if VE_DBG
				printf("c | Gate %d = XOR(", ABS(dx));
				for (int k = 0; k < arity; k++) {
					printf(" %d", ABS(out_c[k]));
					if (k < arity - 1) printf(",");
				}
				printf(" ) found ==> added = %d, deleted = %d\f_dx", nAddedCls, itarget.size() + otarget.size());
				pClauseSet(cnf, itarget, otarget);
#endif
				// substitute
				if (nAddedCls) substitute_x(v, nAddedCls, nAddedLits, cnf, itarget, otarget, units, out_c);
				return true;
			}
			freeze_arities(cnf, itarget);
			return false;
		}

		// 3-phase approach
		_PFROST_D_ void toblivion(CNF& cnf, OL& poss, OL& negs)
		{
#pragma unroll
			for (S_REF* i = poss; i != poss.end(); i++) cnf[*i].markDeleted();
#pragma unroll
			for (S_REF* i = negs; i != negs.end(); i++) cnf[*i].markDeleted();
			poss.clear(true), negs.clear(true);
		}

		_PFROST_D_ void resolve_x(const uint32& x, const uint32& nAddedCls, uint32& addedPos, S_REF& addedRef, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			assert(nAddedCls);
#if VE_DBG
			printf("c | Resolving %d, added = %d, addedPos = %d, addedRef = %lld:\n", x, nAddedCls, addedPos, addedRef);
#endif
			int nbuckets = dc_nbuckets;
			uint32 checksum = addedPos + nAddedCls;
			S_REF* cs = cnf.csData();
#pragma unroll
			for (S_REF* i = me; i != me.end() && addedPos < checksum; i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && addedPos < checksum; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					int max_sz = ci.size() + cj.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT2) {
						// must use "isTautology" check first to avoid data racing on global memory
						if (!isTautology(x, ci, cj)) {
							SCLAUSE* added = new (cnf.cref(addedRef)) SCLAUSE();
							merge(x, ci, cj, added), max_sz = added->size();
							added->markAdded();
							if (max_sz == 1) units->push(**added);
							cs[addedPos++] = addedRef, addedRef += (max_sz - 1) + nbuckets;
#if VE_DBG
							printf("c | C(%d, r: %lld)->", addedPos - 1, addedRef - added->blockSize()), added->print();
#endif
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT2);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(addedRef);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							added->markAdded();
							if (merged_sz == 1) units->push(*out_c);
							cs[addedPos++] = addedRef, addedRef += (merged_sz - 1) + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
#if VE_DBG
							printf("c | C(%d, r: %lld)->", addedPos - 1, addedRef - added->blockSize()), added->print();
#endif
						}
					}
				}
			}
			assert(checksum == addedPos);
			// delete resolved
			toblivion(cnf, me, other);
		}

		_PFROST_D_ void substitute_x(const uint32& x, const uint32& nAddedCls, uint32& addedPos, S_REF& addedRef, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			assert(nAddedCls);
#if VE_DBG
			printf("c | Substituting %d, added = %d, addedPos = %d, addedRef = %d:\n", x, nAddedCls, addedPos, addedRef);
#endif
			int nbuckets = dc_nbuckets;
			uint32 checksum = addedPos + nAddedCls;
			S_REF* cs = cnf.csData();
#pragma unroll
			for (S_REF* i = me; i != me.end() && addedPos < checksum; i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && addedPos < checksum; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() == cj.molten()) continue;
					int max_sz = ci.size() + cj.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT2) {
						// must use "isTautology" check first to avoid data racing on global memory
						if (!isTautology(x, ci, cj)) { 
							SCLAUSE* added = new (cnf.cref(addedRef)) SCLAUSE();
							merge(x, ci, cj, added), max_sz = added->size();
							added->markAdded();
							if (max_sz == 1) units->push(**added);
							cs[addedPos++] = addedRef, addedRef += (max_sz - 1) + nbuckets;
#if VE_DBG
							printf("c | C(%d, r: %d)->", addedPos - 1, addedRef - added->blockSize()), added->print();
#endif
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT2);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(addedRef);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							added->markAdded();
							if (merged_sz == 1) units->push(*out_c);
							cs[addedPos++] = addedRef, addedRef += (merged_sz - 1) + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
#if VE_DBG
							printf("c | C(%d, r: %d)->", addedPos - 1, addedRef - added->blockSize()), added->print();
#endif
						}
					}
				}
			}
			assert(checksum == addedPos);
			// delete resolved
			toblivion(cnf, me, other);
		}

		_PFROST_D_ bool resolve(const uint32& x, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved, uint32& nAddedCls, uint32& nAddedLits, const bool& bound)
		{
			assert(x);
			assert(checkMolten(cnf, poss, negs));
			uint32 p = V2L(x);
			uint32 psLits = 0, nsLits = 0;
			// count clauses/literals before and after
			nAddedCls = 0, nAddedLits = 0;
			calcResolvents(x, cnf, poss, negs, pOrgs, nOrgs, nAddedCls, nAddedLits);
			if (nAddedCls) {
				if (nAddedCls > pOrgs + nOrgs) return false;
				countLitsBefore(cnf, poss, psLits);
				countLitsBefore(cnf, negs, nsLits);
				if (bound && nAddedLits > psLits + nsLits) return false;
			}
#if VE_DBG
			printf("c | Resolving(%d) ==> added = %d, deleted = %d\n", x, nAddedCls, poss.size() + negs.size());
			pClauseSet(cnf, poss, negs);
#endif
			// save deleted
			uint32 nSaved, lit;
			bool which = pOrgs > nOrgs;
			if (which) nSaved = nOrgs + nsLits + 2, lit = p;
			else nSaved = pOrgs + psLits + 2, lit = NEG(p);
			if (nAddedCls) saveResolved(lit, nSaved, cnf, which ? negs : poss, resolved);
			else toblivion(lit, nSaved, cnf, which ? negs : poss, which ? poss : negs, resolved);
#if VE_DBG
			printf("c | Resolved"), resolved->print();
#endif
			return true; // resolution successful
		}

		_PFROST_D_ uint32 find_equ_gate(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			assert(p > 1);
			assert(!SIGN(p));
			assert(checkMolten(cnf, poss, negs));
			assert(pOrgs && nOrgs);
			uint32 first;
			if (first = find_sfanin(p, cnf, poss)) {
				uint32 second = NEG(p), def = first;
				if (second < def) first = second, second = def;
				assert(first < second);
				for (S_REF* i = negs; i != negs.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original() && c.size() == 2 && c[0] == first && c[1] == second) {
						assert(def == first || def == second);
#if VE_DBG
						printf("c | Gate %d = -/+%d found\n", ABS(p), ABS(def));
						pClauseSet(cnf, poss, negs);
#endif
						saveResolved(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
#if VE_DBG
						printf("c | Resolved"), resolved->print();
#endif
						return def;
					}
				}
			}
			return 0;
		}

		_PFROST_D_ bool find_ao_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, cuVecU* resolved, uint32* out_c, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			uint32 sig, x = ABS(dx);
			// (-) ==> look for AND , (+) ==> look for OR
#if VE_DBG
			const char* type = SIGN(dx) ? "AND" : "OR";
#endif
			OL& itarget = ot[dx];
			int nImps = find_fanin(dx, cnf, itarget, out_c, sig);
			if (nImps > 1) {
				uint32 f_dx = FLIP(dx);
				out_c[nImps++] = f_dx;
				sig |= MAPHASH(f_dx);
				devSort(out_c, nImps);
				OL& otarget = ot[f_dx];
				for (S_REF* i = otarget; i != otarget.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.learnt()) continue;
					if (c.size() == nImps && sub(c.sig(), sig) && isEqual(c, out_c, nImps)) {
						c.melt(); // mark as fanout clause
						// check resolvability
						nAddedCls = 0, nAddedLits = 0;
						countSubstituted(x, cnf, itarget, otarget, nAddedCls, nAddedLits);
						if (nAddedCls > nOrgCls) { c.freeze(); break; }
						// can be substituted
#if VE_DBG
						printf("c | Gate %d = %s(", ABS(dx), type);
						for (int k = 0; k < nImps; k++) {
							if (ABS(out_c[k]) == ABS(dx)) continue;
							printf(" %d", ABS(out_c[k]));
							if (k < nImps - 1) printf(",");
						}
						printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
						printGate(cnf, itarget, otarget);
#endif
						return true;
					}
				}
			}
			freeze_binaries(cnf, itarget);
			return false;
		}

		_PFROST_D_ bool find_ite_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, cuVecU* resolved, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			OL& itarget = ot[dx];
			for (S_REF* i = itarget; i != itarget.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt() || ci.size() < 3 || ci.size() > 3) continue;
				assert(ci.original());
				uint32 xi = ci[0], yi = ci[1], zi = ci[2];
				if (yi == dx) devSwap(xi, yi);
				if (zi == dx) devSwap(xi, zi);
				assert(xi == dx);
				for (S_REF* j = i + 1; j != itarget.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt() || cj.size() < 3 || cj.size() > 3) continue;
					assert(cj.original());
					uint32 xj = cj[0], yj = cj[1], zj = cj[2];
					if (yj == dx) devSwap(xj, yj);
					if (zj == dx) devSwap(xj, zj);
					assert(xj == dx);
					if (ABS(yi) == ABS(zj)) devSwap(yj, zj);
					if (ABS(zi) == ABS(zj)) continue;
					if (yi != FLIP(yj)) continue;
					uint32 f_dx = FLIP(dx);
					S_REF r1 = fast_equality_check(cnf, ot, f_dx, yi, FLIP(zi));
					if (r1 == GNOREF) continue;
					S_REF r2 = fast_equality_check(cnf, ot, f_dx, yj, FLIP(zj));
					if (r2 == GNOREF) continue;
					assert(cnf[r1].original());
					assert(cnf[r2].original());
					// mark gate clauses
					ci.melt(), cj.melt();
					cnf[r1].melt(), cnf[r2].melt();
					// check resolvability
					uint32 v = ABS(dx);
					OL& otarget = ot[f_dx];
					nAddedCls = 0, nAddedLits = 0;
					countSubstituted(v, cnf, itarget, otarget, nAddedCls, nAddedLits);
					if (nAddedCls > nOrgCls) {
						ci.freeze(), cj.freeze();
						cnf[r1].freeze(), cnf[r2].freeze();
						return false;
					}
					// can be substituted
#if VE_DBG
					printf("c | Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d", ABS(dx), ABS(yi), ABS(zi), ABS(zj), nAddedCls, itarget.size() + otarget.size());
					printGate(cnf, itarget, otarget);
#endif
					return true;
				}
			}
			return false;
		}

		_PFROST_D_ bool find_xor_gate(const uint32& dx, const uint32& nOrgCls, const int& xor_max_arity, CNF& cnf, OT& ot, cuVecU* resolved, uint32* out_c, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			OL& itarget = ot[dx];
			for (S_REF* i = itarget; i != itarget.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				int size = ci.size();
				int arity = size - 1; // XOR arity
				if (ci.learnt() || size < 3 || arity > xor_max_arity || arity > SH_MAX_BVE_OUT1) continue;
				assert(ci.original());
				// share to out_c except dx
				shareXORClause(dx, ci, out_c);
				// find arity clauses
				int itargets = 0;
				for (int j = 0; j < arity; j++) {
					S_REF r = find_fanin(dx, j, cnf, itarget, out_c, arity);
					if (r == GNOREF) break;
					assert(cnf[r].original());
					cnf[r].melt(), itargets++;
				}
				if (itargets < arity) {
					freeze_arities(cnf, itarget);
					continue;
				}
				// find all +/-  
				uint32 f_dx = FLIP(dx);
				S_REF r1 = find_all(f_dx, cnf, ot, out_c, arity);
				if (r1 == GNOREF) break;
				flip_all(out_c, arity);
				S_REF r2 = find_all(f_dx, cnf, ot, out_c, arity);
				if (r2 == GNOREF) break;
				assert(cnf[r1].original());
				assert(cnf[r2].original());
				cnf[r1].melt(), cnf[r2].melt();
				// check resolvability
				uint32 v = ABS(dx);
				OL& otarget = ot[f_dx];
				nAddedCls = 0, nAddedLits = 0;
				countSubstituted(v, cnf, itarget, otarget, nAddedCls, nAddedLits);
				if (nAddedCls > nOrgCls) {
					cnf[r1].freeze(), cnf[r2].freeze();
					break;
				}
				// can be substituted
#if VE_DBG
				printf("c | Gate %d = XOR(", ABS(dx));
				for (int k = 0; k < arity; k++) {
					printf(" %d", ABS(out_c[k]));
					if (k < arity - 1) printf(",");
				}
				printf(" ) found ==> added = %d, deleted = %d\f_dx", nAddedCls, itarget.size() + otarget.size());
				printGate(cnf, itarget, otarget);
#endif
				return true;
			}
			freeze_arities(cnf, itarget);
			return false;
		}

	} // sigma namespace
} // parafrost namespace


#endif