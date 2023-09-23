/***********************************************************************[function.cuh]
Copyright(c) 2021, Muhammad Osama - Anton Wijs,
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

#ifndef __FUN_
#define __FUN_

#include "elimination.cuh"

namespace ParaFROST {

	// inspired by the function table reassoning in Lingeling

	#define FUN_DBG 0

	__constant__ uint64 MAGICCONSTS[6] = {
	  0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL, 0xf0f0f0f0f0f0f0f0ULL,
	  0xff00ff00ff00ff00ULL, 0xffff0000ffff0000ULL, 0xffffffff00000000ULL,
	};

	constexpr int	 MAXFUNVAR = 12;
	constexpr uint32 FUNTABLEN = 64;
	constexpr uint64 ALLONES = ~0ULL;

	typedef uint64 Fun[FUNTABLEN];

	_PFROST_D_ void falsefun(Fun f)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			f[i] = 0ULL;
	}

	_PFROST_D_ void truefun(Fun f)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			f[i] = ALLONES;
	}

	_PFROST_D_ bool isfalsefun(const Fun f)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			if (f[i]) return false;
		return true;
	}

	_PFROST_D_ bool istruefun(const Fun f)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			if (f[i] != ALLONES) return false;
		return true;
	}

	_PFROST_D_ void orfun(Fun a, const Fun b)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			a[i] |= b[i];
	}

	_PFROST_D_ void andfun(Fun a, const Fun b)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			a[i] &= b[i];
	}

	_PFROST_D_ void copyfun(Fun a, const Fun b)
	{
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			a[i] = b[i];
	}

	_PFROST_D_ void clause2fun(const int& v, const bool& sign, Fun f)
	{
		assert(v >= 0 && v < MAXFUNVAR);
		if (v < 6) {
			uint64 val = MAGICCONSTS[v];
			if (sign) val = ~val;
			#pragma unroll
			for (uint32 i = 0; i < FUNTABLEN; ++i)
				f[i] |= val;
		}
		else {
			uint64 val = sign ? ALLONES : 0ULL;
			int j = 0;
			int sv = 1 << (v - 6);
			#pragma unroll
			for (uint32 i = 0; i < FUNTABLEN; ++i) {
				f[i] |= val;
				if (++j >= sv) {
					val = ~val;
					j = 0;
				}
			}
		}
	}

	_PFROST_D_ bool buildfuntab
	(
		const uint32& lit,
		const uint32* varcore,
		CNF& cnf,
		OT& ot, 
		Fun cls, // shared memory
		Fun f    // local memory
	)
	{
		assert(lit > 1);
		truefun(f);
		OL& list = ot[lit];
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.learnt()) continue;
			assert(!c.deleted());
			falsefun(cls);
			forall_clause(c, k) {
				const uint32 other = *k;
				if (other == lit) continue;
				assert(other != FLIP(lit));
				// mvar[0 : MAXFUNVAR - 1] here has to be unique 
				// and not repeated or in other words mapped 
				// monotonically from 0 to max. frozen variable
				const uint32 mvar = varcore[ABS(other)]; 
				if (mvar >= MAXFUNVAR) return false; 
				clause2fun(mvar, SIGN(other), cls);
			}
			assert(!isfalsefun(cls));
			assert(!istruefun(cls));
			andfun(f, cls);
		}
		return true;
	}

	_PFROST_D_ void buildfuntab(
		const uint32& lit,
		const uint32* varcore,
		const int& tail, 
		const OL& ol, 
		CNF& cnf, 
		Fun cls, // shared memory
		Fun fun, // local memory
		bool& core)
	{
		assert(lit > 1);
		for (int j = 0; j < tail; ++j) {
			SCLAUSE& c = cnf[ol[j]];
			if (c.learnt()) continue;
			assert(!c.deleted());
			falsefun(cls);
			forall_clause(c, k) {
				const uint32 other = *k;
				if (other == lit) continue;
				assert(other != FLIP(lit));
				const uint32 mvar = varcore[ABS(other)];
				assert(mvar < MAXFUNVAR);
				clause2fun(mvar, SIGN(other), cls);
			}
			assert(!isfalsefun(cls));
			assert(!istruefun(cls));
			andfun(fun, cls);
		}
		if (isfalsefun(fun)) {
			cnf[ol[tail]].melt(); // non-gate clause (not resolved with its kind)
			core = true;
		}
	}

	_PFROST_D_ uint64 collapsefun(const Fun b, const Fun c)
	{
		uint64 allzero = 0;
		#pragma unroll
		for (uint32 i = 0; i < FUNTABLEN; ++i)
			allzero |= (b[i] & c[i]);
		return allzero;
	}

	_PFROST_D_ bool countCoreSubstituted(
		const uint32& x,
		const uint32& nClsBefore,
		CNF& cnf,
		OL& me,
		OL& other,
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(x);
		assert(!nElements);
		assert(!nAddedCls);
		assert(!nAddedLits);
		const int rlimit = kOpts->ve_clause_max;
		// check if proof bytes has to be calculated
		if (kOpts->proof_en) {
			uint32 proofBytes = 0;
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				const bool ci_m = ci.molten();
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.original() && (!ci_m || !cj.molten())) {
						const int rsize = mergeProof(x, ci, cj, proofBytes);
						if (rsize == 1)
							nElements++;
						else if (rsize) {
							if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}

			// GUARD for compressed proof size and #units
			if (nElements > ADDEDCLS_MAX || proofBytes > ADDEDPROOF_MAX) return true;

			nElements = ENCODEPROOFINFO(nElements, proofBytes);

			#if FUN_DBG
			printf("c  Variable %d: counted %d units and %d proof bytes\n", x, nElements, proofBytes);
			#endif
		}
		else {
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				const bool ci_m = ci.molten();
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.original() && (!ci_m || !cj.molten())) {
						const int rsize = merge(x, ci, cj);
						if (rsize == 1)
							nElements++;
						else if (rsize) {
							if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}
		}

		// GUARD for compressed variable limits
		if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX) return true;

		// check bound on literals
		if (kOpts->ve_lbound_en) {
			uint32 nLitsBefore = 0;
			countLitsBefore(cnf, me, nLitsBefore);
			countLitsBefore(cnf, other, nLitsBefore);
			if (nAddedLits > nLitsBefore) return true;
		}

		return false;
	}

	_PFROST_D_ bool find_fun_gate(
		const uint32& p,
		const uint32& n,
		const uint32& nOrgCls,
		const uint32* varcore,
		CNF& cnf, 
		OT& ot,
		uint32* shmem, 
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(p > 1);
		assert(varcore);
		assert(checkMolten(cnf, ot[p], ot[n]));

		uint64* cls = (uint64*)shmem;

		Fun pos, neg;
		if (buildfuntab(p, varcore, cnf, ot, cls, pos) && buildfuntab(n, varcore, cnf, ot, cls, neg)) {
			if (!collapsefun(pos, neg)) {
				// core optimization (minimality not guaranteed)
				Fun& fun = pos;
				OL& poss = ot[p];
				bool core = false;
				for (int i = poss.size() - 1; i >= 0; i--) {
					copyfun(fun, neg);
					if (cnf[poss[i]].original())
						buildfuntab(p, varcore, i, poss, cnf, cls, fun, core);
				}
				OL& negs = ot[n];
				for (int i = negs.size() - 1; i >= 0; i--) {
					truefun(fun);
					if (cnf[negs[i]].original())
						buildfuntab(n, varcore, i, negs, cnf, cls, fun, core);
				}
				// check resolvability
				nElements = 0, nAddedCls = 0, nAddedLits = 0;
				if (countCoreSubstituted(ABS(p), nOrgCls, cnf, poss, negs, nElements, nAddedCls, nAddedLits)) {
					if (core) freezeClauses(cnf, poss, negs);
					return false;
				}
				// can be substituted
				#if FUN_DBG
				printf("c  Fun Gate %d ", ABS(p));
				printf("found ==> added = %d, deleted = %d\n", nAddedCls, poss.size() + negs.size());
				pClauseSet(cnf, ot, ABS(p));
				#endif	
				return true;
			}
		}
		return false;
	}

} 


#endif