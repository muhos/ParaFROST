/***********************************************************************[pfbve.h]
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

namespace pFROST {

	namespace SIGmA {

#define VE_DBG 0

		_PFROST_D_ void	saveResolved(uint32*& saved, const uint32& lit)
		{
			saved[0] = 0, saved[1] = lit, saved[2] = 0;
			saved += 3;
		}

		_PFROST_D_ void	saveResolved(uint32*& saved, SCLAUSE& c, const uint32& x)
		{
			assert(c.original());
			saveResolved(saved, x);
			uint32* lit = c, * cend = c.end();
			while (lit != cend)
				*saved++ = *lit++;
		}

		_PFROST_D_ bool isTautology_and(SCLAUSE& org, uint32* defs, const int& nDefs)
		{
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
				if (org.has(*d)) return true;
			}
			return false;
		}

		_PFROST_D_ bool isTautology_or(SCLAUSE& org, uint32* defs, const int& nDefs)
		{
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
				if (org.has(FLIP(*d))) return true;
			}
			return false;
		}

		_PFROST_D_ bool clause_extend_and(const uint32& neg_x, SCLAUSE& org, uint32* defs, const int& nDefs, SCLAUSE& out_c)
		{
			assert(neg_x);
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
			out_c.clear();
			if (isTautology_and(org, defs, nDefs)) return false; // tautology
			// extend 
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == neg_x) {
#pragma unroll
					for (uint32* d = defs; d != defs + nDefs; d++) out_c.push(FLIP(*d));
				}
				else out_c.push(lit);
			}
			// attach 
			devSort(out_c, out_c.size());
			out_c.filter();
			out_c.set_status(ORIGINAL);
			assert(out_c.hasZero() < 0);
			assert(out_c.isSorted());
			return true;
		}

		_PFROST_D_ int clause_extend_and(const uint32& neg_x, SCLAUSE& org, uint32* defs, const int& nDefs, uint32* out_c)
		{
			assert(neg_x);
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
			if (isTautology_and(org, defs, nDefs)) return 0; // tautology
			int len = 0;
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == neg_x) {
#pragma unroll
					for (uint32* d = defs; d != defs + nDefs; d++) out_c[len++] = FLIP(*d);
				}
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_extend_or(const uint32& x, SCLAUSE& org, uint32* defs, const int& nDefs, SCLAUSE& out_c)
		{
			assert(x);
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
			out_c.clear();
			if (isTautology_or(org, defs, nDefs)) return false; // tautology
			// extend
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == x) {
#pragma unroll
					for (uint32* d = defs; d != defs + nDefs; d++) out_c.push(*d);
				}
				else out_c.push(lit);
			}
			// attach 
			devSort(out_c, out_c.size());
			out_c.filter();
			out_c.set_status(ORIGINAL);
			assert(out_c.hasZero() < 0);
			assert(out_c.isSorted());
			return true;
		}

		_PFROST_D_ int clause_extend_or(const uint32& x, SCLAUSE& org, uint32* defs, const int& nDefs, uint32* out_c)
		{
			assert(x);
			assert(org.original());
			assert(org.size() > 1);
			assert(nDefs > 1);
			if (isTautology_or(org, defs, nDefs)) return 0; // tautology
			int len = 0;
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == x) {
#pragma unroll
					for (uint32* d = defs; d != defs + nDefs; d++) out_c[len++] = *d;
				}
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_split_and(const uint32& x, const uint32& def, SCLAUSE& org, SCLAUSE& out_c)
		{
			assert(x);
			assert(def);
			assert(org.original());
			assert(org.size() > 1);
			out_c.clear();
			if (org.has(FLIP(def))) return false;  // tautology
			// split
			uint32 sig = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == def) continue; // repeated literal
				if (lit == x) { out_c.push(def); sig |= MAPHASH(def); }
				else { out_c.push(lit); sig |= MAPHASH(lit); }
			}
			// attach
			devSort(out_c, out_c.size());
			out_c.set_sig(sig);
			out_c.set_status(ORIGINAL);
			assert(out_c.isSorted());
			assert(out_c.hasZero() < 0);
			return true;
		}

		_PFROST_D_ int clause_split_and(const uint32& x, const uint32& def, SCLAUSE& org, uint32* out_c)
		{
			assert(x);
			assert(def);
			assert(org.original());
			assert(org.size() > 1);
			if (org.has(FLIP(def))) return 0;
			int len = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == def) continue; // repeated literal
				if (lit == x) out_c[len++] = def;
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_split_or(const uint32& neg_x, const uint32& def, SCLAUSE& org, SCLAUSE& out_c)
		{
			assert(neg_x);
			assert(def);
			assert(org.original());
			assert(org.size() > 1);
			out_c.clear();
			if (org.has(def)) return false;  // tautology
			// split
			uint32 sig = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == FLIP(def)) continue; // repeated literal
				if (lit == neg_x) { out_c.push(FLIP(def)); sig |= MAPHASH(FLIP(def)); }
				else { out_c.push(lit); sig |= MAPHASH(lit); }
			}
			// attach
			devSort(out_c, out_c.size());
			out_c.set_sig(sig);
			out_c.set_status(ORIGINAL);
			assert(out_c.isSorted());
			assert(out_c.hasZero() < 0);
			return true;
		}

		_PFROST_D_ int clause_split_or(const uint32& neg_x, const uint32& def, SCLAUSE& org, uint32* out_c)
		{
			assert(neg_x);
			assert(def);
			assert(org.original());
			assert(org.size() > 1);
			if (org.has(def)) return 0;  // tautology
			uint32 len = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == FLIP(def)) continue; // repeated literal
				if (lit == neg_x) out_c[len++] = FLIP(def);
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ void countLitsBefore(CNF& cnf, OL& list, uint32& nLitsBefore)
		{
#pragma unroll
			for (uint32 *i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) nLitsBefore += c.size();
			}
		}

		_PFROST_D_ void calcResolvents(const uint32& x, CNF& cnf, OL& poss, OL& negs, const uint32& pOrgs, const uint32& nOrgs, uint32& addedCls, uint32& addedLits)
		{
			assert(x);
			assert(pOrgs);
			assert(nOrgs);
			assert(!addedLits);
			uint32 nTs = 0;
#pragma unroll
			for (uint32* i = poss; i != poss.end(); i++) {
				SCLAUSE& pos = cnf[*i];
				if (pos.learnt()) continue;
#pragma unroll
				for (uint32* j = negs; j != negs.end(); j++) {
					SCLAUSE& neg = cnf[*j];
					if (neg.learnt()) continue;
					if (isTautology(x, pos, neg)) nTs++;
					else addedLits += pos.size() + neg.size() - 2;
				}
			}
			assert(pOrgs * nOrgs >= nTs);
			addedCls = pOrgs * nOrgs - nTs;
		}

		_PFROST_D_ void toblivion(CNF& cnf, OL& poss, OL& negs, cuVecU* resolved, const uint32& pOrgs, const uint32& nOrgs, const uint32& psLits, const uint32& nsLits, const uint32& p)
		{
			bool which = pOrgs > nOrgs;
			if (which) {
				uint32 nElems = nOrgs * 3 + nsLits;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (uint32* i = negs; i != negs.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, NEG(p));
					c.markDeleted();
				}
			}
			else {
				uint32 nElems = pOrgs * 3 + psLits;
				uint32* saved = resolved->jump(nElems);
#pragma unroll
				for (uint32* i = poss; i != poss.end(); i++) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveResolved(saved, c, p);
					c.markDeleted();
				}
			}
			OL& other = which ? poss : negs;
#pragma unroll
			for (uint32* i = other; i != other.end(); i++) cnf[*i].markDeleted();
			poss.clear(true), negs.clear(true);
#if VE_DBG
			printf("c | Resolved"), resolved->print();
#endif
		}

		_PFROST_D_ bool substitute_AND(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* defs, const int& nDefs, uint32* out_c)
		{
			uint32 nAddedCls = 0, nAddedLits = 0;
			// count number of added cls & literals for negatives
#pragma unroll
			for (uint32* i = negs; i != negs.end(); i++) {
				SCLAUSE& neg = cnf[*i];
				if (neg.learnt()) continue;
				if (!isTautology_and(neg, defs, nDefs)) {
					nAddedCls++;
					nAddedLits += neg.size() + nDefs - 1;
				}
			}
			// count number of added cls & literals for positives
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
#pragma unroll
				for (uint32* i = poss; i != poss.end(); i++) {
					SCLAUSE& pos = cnf[*i];
					if (pos.learnt()) continue;
					if (!pos.has(FLIP(*d))) {
						nAddedCls++;
						nAddedLits += pos.size();
					}
				}
			}
			if (nAddedCls == 0) return true; // No substitutions to add
			if (nAddedCls > pOrgs + nOrgs) return false;
#if VE_DBG
			printf("c | AND(%d): added = %d, deleted = %d\n", ABS(p), nAddedCls, poss.size() + negs.size()), pClauseSet(cnf, poss, negs);
#endif
			S_REF ref, checksum = 0, n = NEG(p);
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
			// substitute negatives 
#pragma unroll
			for (uint32* i = negs; i != negs.end() && checksum < nAddedCls; i++) {
				SCLAUSE& neg = cnf[*i];
				if (neg.learnt()) continue;
				int max_sz = neg.size() + nDefs - 1;
				if (max_sz > SH_MAX_BVE_OUT) { // use global memory
					if (clause_extend_and(n, neg, defs, nDefs, cnf[ref])) {
						SCLAUSE& added = cnf[ref];
						if (added.size() == 1) units->push(*added);
						cs[checksum++] = ref, ref += added.blockSize();
					}
				}
				else { // use shared memory
					int ext_sz = 0;
					if ((ext_sz = clause_extend_and(n, neg, defs, nDefs, out_c))) {
						assert(ext_sz <= SH_MAX_BVE_OUT);
						SCLAUSE* added = cnf.cref(ref);
						added = new (added) SCLAUSE();
						added->filterCopy(out_c, ext_sz);
						if (ext_sz == 1) units->push(*out_c);
						cs[checksum++] = ref, ref += (ext_sz - 1) + sizeof(S_REF);
					}
				}
#if VE_DBG
				printf("c | Added %d substitutions\n", checksum);
				printf("c | Units"), units->print();
#endif
			}
			// substitute positives
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
#pragma unroll
				for (uint32* i = poss; i != poss.end() && checksum < nAddedCls; i++) {
					SCLAUSE& pos = cnf[*i];
					if (pos.learnt()) continue;
					int max_sz = pos.size();
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (clause_split_and(p, *d, pos, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int split_sz = 0;
						if (split_sz = clause_split_and(p, *d, pos, out_c)) {
							assert(split_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, split_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, split_sz);
							added->set_sig(sig);
							if (split_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (split_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
#if VE_DBG
					printf("c | Added %d substitutions\n", checksum);
					printf("c | Units"), units->print();
#endif
				}
			}
			assert(checksum == nAddedCls);
			return true; // AND-substitution successful
		}

		_PFROST_D_ bool substitute_OR(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* defs, const int& nDefs, uint32* out_c)
		{
			uint32 nAddedCls = 0, nAddedLits = 0;
			// count number of added cls & literals for positives
#pragma unroll
			for (uint32* i = poss; i != poss.end(); i++) {
				SCLAUSE& pos = cnf[*i];
				if (pos.learnt()) continue;
				if (!isTautology_or(pos, defs, nDefs)) {
					nAddedCls++;
					nAddedLits += pos.size() + nDefs - 1;
				}
			}
			// count number of added cls & literals for negatives
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
#pragma unroll
				for (uint32* i = negs; i != negs.end(); i++) {
					SCLAUSE& neg = cnf[*i];
					if (neg.learnt()) continue;
					if (!neg.has(*d)) {
						nAddedCls++;
						nAddedLits += neg.size();
					}
				}
			}
			if (nAddedCls == 0) return true; // No substitutions to add
			if (nAddedCls > pOrgs + nOrgs) return false;
#if VE_DBG
			printf("c | OR(%d): added = %d, deleted = %d\n", ABS(p), nAddedCls, poss.size() + negs.size()), pClauseSet(cnf, poss, negs);
#endif
			S_REF ref, checksum = 0;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
			// substitute positives
#pragma unroll
			for (uint32* i = poss; i != poss.end() && checksum < nAddedCls; i++) {
				SCLAUSE& pos = cnf[*i];
				if (pos.learnt()) continue;
				int max_sz = pos.size() + nDefs - 1;
				if (max_sz > SH_MAX_BVE_OUT) { // use global memory
					if (clause_extend_or(p, pos, defs, nDefs, cnf[ref])) {
						SCLAUSE& added = cnf[ref];
						if (added.size() == 1) units->push(*added);
						cs[checksum++] = ref, ref += added.blockSize();
					}
				}
				else { // use shared memory
					int ext_sz = 0;
					if ((ext_sz = clause_extend_or(p, pos, defs, nDefs, out_c))) {
						assert(ext_sz <= SH_MAX_BVE_OUT);
						SCLAUSE* added = cnf.cref(ref);
						added = new (added) SCLAUSE();
						added->filterCopy(out_c, ext_sz);
						if (ext_sz == 1) units->push(*out_c);
						cs[checksum++] = ref, ref += (ext_sz - 1) + sizeof(S_REF);
					}
				}
#if VE_DBG
				printf("c | Added %d substitutions\n", checksum);
				printf("c | Units"), units->print();
#endif
			}
			// substitute negatives
			uint32 n = NEG(p);
#pragma unroll
			for (uint32* d = defs; d != defs + nDefs; d++) {
#pragma unroll
				for (uint32* i = negs; i != negs.end() && checksum < nAddedCls; i++) {
					SCLAUSE& neg = cnf[*i];
					if (neg.learnt()) continue;
					int max_sz = neg.size();
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (clause_split_or(n, *d, neg, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int split_sz = 0;
						if (split_sz = clause_split_or(n, *d, neg, out_c)) {
							assert(split_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, split_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, split_sz);
							added->set_sig(sig);
							if (split_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (split_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
#if VE_DBG
					printf("c | Added %d substitutions\n", checksum);
					printf("c | Units"), units->print();
#endif
				}
			}
			assert(checksum == nAddedCls);
			return true; // OR-substitution successful
		}

		_PFROST_D_ int find_fanin(const uint32& gate_out, CNF& cnf, OL& list, uint32* out_c, uint32& sig)
		{
			assert(gate_out > 1);
			sig = 0;
			uint32 imp = 0;
			int nImps = 0;
#pragma unroll
			for (uint32* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.learnt()) continue;
				if (c.size() == 2) {
					if (c[0] == gate_out) { // found gate output 
						imp = FLIP(c[1]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
					else if (c[1] == gate_out) {
						imp = FLIP(c[0]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
				}
			}
			return nImps;
		}

		_PFROST_D_ bool gateReasoning_x(const uint32& p, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, cuVecU* resolved, uint32* defs, uint32* out_c)
		{
			uint32 sig, n = NEG(p);
			// AND-gate Reasoning
			if (nOrgs > (FAN_LMT - 1)) return false; // shared memory GUARD
			int nImps = find_fanin(n, cnf, negs, out_c, sig);
			if (nImps > 1) {
				assert(nImps <= FAN_LMT - 1);
				out_c[nImps++] = p;
				sig |= MAPHASH(p);
				devSort(out_c, nImps);
#pragma unroll
				for (uint32* i = poss; i != poss.end(); i++) {
					SCLAUSE& pos = cnf[*i];
					if (pos.learnt()) continue;
					int pos_size = pos.size();
					if (pos_size == nImps && sub(pos.sig(), sig) && isEqual(pos, out_c, nImps)) {
						int nDefs = 0;
#pragma unroll
						for (int j = 0; j < nImps; j++) {
							uint32 lit = out_c[j];
							if (lit != p) defs[nDefs++] = FLIP(lit);
						}
						assert(nDefs < FAN_LMT);
						assert(nDefs == pos_size - 1);
						if (substitute_AND(p, pOrgs, nOrgs, cnf, poss, negs, units, defs, nDefs, out_c)) {
							uint32 psLits = 0, nsLits = 0;
							countLitsBefore(cnf, poss, psLits), countLitsBefore(cnf, negs, nsLits);
							toblivion(cnf, poss, negs, resolved, pOrgs, nOrgs, psLits, nsLits, p);
							return true;
						}
					}
				}
			}
			// OR-gate Reasoning
			if (pOrgs > (FAN_LMT - 1)) return false; // shared memory GUARD
			nImps = find_fanin(p, cnf, poss, out_c, sig);
			if (nImps > 1) {
				assert(nImps <= FAN_LMT - 1);
				out_c[nImps++] = n;
				sig |= MAPHASH(n);
				devSort(out_c, nImps);
#pragma unroll
				for (uint32* i = negs; i != negs.end(); i++) {
					SCLAUSE& neg = cnf[*i];
					if (neg.learnt()) continue;
					int neg_size = neg.size();
					if (neg_size == nImps && sub(neg.sig(), sig) && isEqual(neg, out_c, nImps)) {
						int nDefs = 0;
#pragma unroll
						for (int j = 0; j < nImps; j++) {
							uint32 lit = out_c[j];
							if (lit != n) defs[nDefs++] = lit;
						}
						assert(nDefs < FAN_LMT);
						assert(nDefs == neg_size - 1);
						if (substitute_OR(p, pOrgs, nOrgs, cnf, poss, negs, units, defs, nDefs, out_c)) {
							uint32 psLits = 0, nsLits = 0;
							countLitsBefore(cnf, poss, psLits), countLitsBefore(cnf, negs, nsLits);
							toblivion(cnf, poss, negs, resolved, pOrgs, nOrgs, psLits, nsLits, p);
							return true;
						}
					}
				}
			}
			return false;
		}

		_PFROST_D_ bool resolve_x(const uint32& x, const uint32& pOrgs, const uint32& nOrgs, CNF& cnf, OL& poss, OL& negs, cuVecU* units, cuVecU* resolved, uint32* out_c)
		{
			uint32 p = V2D(x);
			uint32 nAddedCls = 0, nAddedLits = 0, psLits = 0, nsLits = 0;
			// count clauses/literals before and after
			calcResolvents(x, cnf, poss, negs, pOrgs, nOrgs, nAddedCls, nAddedLits);
			if (nAddedCls == 0) { // No resolvents to add
				countLitsBefore(cnf, poss, psLits), countLitsBefore(cnf, negs, nsLits);
				toblivion(cnf, poss, negs, resolved, pOrgs, nOrgs, psLits, nsLits, p);
				return true;
			}
			if (nAddedCls > pOrgs + nOrgs) return false;
			countLitsBefore(cnf, poss, psLits), countLitsBefore(cnf, negs, nsLits);
			if (nAddedLits > psLits + nsLits) return false;
#if VE_DBG
			printf("c | Resolving %d -> adding %d, deleting %d\n", x, nAddedCls, poss.size() + negs.size()), pClauseSet(cnf, poss, negs);
#endif
			// resolve x
			S_REF ref, checksum = 0;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (uint32* i = poss; i != poss.end(); i++) {
				SCLAUSE& pos = cnf[*i];
				if (pos.learnt()) continue;
#pragma unroll
				for (uint32* j = negs; j != negs.end() && checksum < nAddedCls; j++) {
					SCLAUSE& neg = cnf[*j];
					if (neg.learnt()) continue;
					int max_sz = pos.size() + neg.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (merge(x, pos, neg, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if ((merged_sz = merge(x, pos, neg, out_c))) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							if (merged_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (merged_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
				}
			}
#if VE_DBG
			printf("c | Added %d resolvents\n", checksum);
			printf("c | Units"), units->print();
#endif
			assert(checksum == nAddedCls);
			toblivion(cnf, poss, negs, resolved, pOrgs, nOrgs, psLits, nsLits, p);
			return true; // resolution successful
		}
	}
}


#endif