/***********************************************************************[proofutils.cuh]
Copyright(c) 2021, Muhammad Osama - Anton Wijs,
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

#ifndef __GPU_PROOFUTILS_
#define __GPU_PROOFUTILS_

#include "primitives.cuh"
#include "printer.cuh"

namespace pFROST {

	#define PROOF_DBG 0

	//=======================================================
	// constant proof lookup table
	#define MAXLEADINGZEROS 30
	__constant__ Byte BLUT[MAXLEADINGZEROS + 1] =
	{
		1, 1, 1, 1, 1, 1,
		2, 2, 2, 2, 2, 2, 2,
		3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5
	};
	#define COUNTBYTES(LIT) BLUT[MAXLEADINGZEROS - __clz(LIT)]
	//=======================================================

	_PFROST_D_ void countBytes(const uint32& lit, Byte& perLit, uint32& total)
	{
		// here we rely on data racing to avoid
		// counting the bytes for a duplicate
		// assuming perLit is initially 0
		if (!perLit) {
			Byte local;
			ORIGINIZELIT(orgLit, lit);
			perLit = local = COUNTBYTES(orgLit);
			total += local;
		}
		else // the more threads taking this path the better
			total += perLit;
		assert(total);
		assert(perLit >= 1 && perLit <= 5);
	}

	_PFROST_D_ void countProofBytes(const uint32 unit, uint32& bytes)
	{
		bytes += 2 + GETBYTECOUNT(unit); // 2: 1 prefix + 1 suffix
		#if PROOF_DBG
		printf("c  counted %d bytes for unit ", bytes), pLit(orgLit), printf("\n");
		#endif
	}

	_PFROST_D_ void countProofBytes(const uint32* lits, const int& size, uint32& bytes)
	{
		assert(size);
		bytes += 2; // prefix + suffix
		#pragma unroll
		for (int i = 0; i < size; i++) {
			bytes += GETBYTECOUNT(lits[i]);
		}
		#if PROOF_DBG
		printf("c  counted %d bytes for clause", bytes), pSharedClause(lits, size);
		#endif
	}

	_PFROST_D_ void countProofBytes(SCLAUSE& c, uint32& bytes)
	{
		assert(c.size());
		bytes += 2; // prefix + suffix
		#pragma unroll
		forall_clause(c, k) {
			bytes += GETBYTECOUNT(*k);
		}
		#if PROOF_DBG
		printf("c  counted %d bytes for clause", bytes), c.print();
		#endif
	}

	_PFROST_D_ void countProofOrg(CNF& cnf, OL& ol, uint32& bytes)
	{
		assert(dc_ptrs->d_lbyte);
		#pragma unroll
		forall_occurs(ol, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original()) 
				countProofBytes(c, bytes);
		}
	}

	_PFROST_D_ void countProofSub(CNF& cnf, OL& ol, uint32& bytes)
	{
		assert(dc_ptrs->d_lbyte);
		#pragma unroll
		forall_occurs(ol, i) {
			SCLAUSE& c = cnf[*i];
			const bool molten = c.molten();
			const bool deleted = c.deleted();
			uint32 local = 0;
			if (molten || deleted) countProofBytes(c, local);
			if (molten) {
				assert(local > 2);
				bytes += local;
			}
			if (deleted) {
				assert(local > 2);
				bytes += local;
			}
		}
	}

	_PFROST_D_ void	saveProofLiteral(addr_t& saved, const uint32& lit)
	{
		ORIGINIZELIT(orgLit, lit);

		#if PROOF_DBG
		printf("c   proof saving literal "), pLit(orgLit), printf(": ");
		#endif

		while (ISLARGE(orgLit)) {
			*saved++ = L2B(orgLit);
			orgLit >>= 7;

			#if PROOF_DBG
			printf("%02X ", *(saved - 1));
			#endif
		}
		assert(orgLit <= BYTEMAX);
		*saved++ = Byte(orgLit);

		#if PROOF_DBG
		printf("%02X\n", orgLit);
		#endif
	}

	_PFROST_D_ void	saveProofUnit(addr_t& saved, const uint32 unit)
	{
		assert(saved);
		*saved++ = PROOF_ADDED;
		saveProofLiteral(saved, unit);
		*saved++ = 0;
	}

	_PFROST_D_ void	saveProofClause(addr_t& saved, const uint32* lits, const int& size, const Byte& state)
	{
		assert(size);
		assert(lits);
		assert(saved);
		assert(state == PROOF_ADDED || state == PROOF_DELETED);

		#if PROOF_DBG
		printf("c  proof saving clause("), pSharedClause(lits, size);
		printf("c   proof saving state '%c'\n", state);
		#endif

		*saved++ = state;
		#pragma unroll
		for (int i = 0; i < size; i++) {
			saveProofLiteral(saved, lits[i]);
		}
		*saved++ = 0;
	}

	_PFROST_D_ void	saveProofClause(addr_t& saved, SCLAUSE& c, const Byte& state)
	{
		assert(c.size());
		assert(saved);
		assert(state == PROOF_ADDED || state == PROOF_DELETED);

		#if PROOF_DBG
		printf("c  proof saving clause("), c.print();
		printf("c   proof saving state '%c'\n", state);
		#endif

		*saved++ = state;
		#pragma unroll
		forall_clause(c, k) {
			saveProofLiteral(saved, *k);
		}
		*saved++ = 0;
	}

	_PFROST_D_ void	addProof(addr_t& saved, CNF& cnf, OL& ol)
	{
		#pragma unroll
		forall_occurs(ol, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original())
				saveProofClause(saved, c, PROOF_ADDED);
		}
	}

	_PFROST_D_ void saveProof(CNF& cnf, OL& ol, addr_t& saved)
	{
		if (ol.empty()) 
			return;
		S_REF* j = ol;
		#pragma unroll
		forall_occurs(ol, i) {
			SCLAUSE& c = cnf[*i];
			const bool deleted = c.deleted();
			if (c.molten()) {
				saveProofClause(saved, c, PROOF_ADDED);
				c.freeze();
			}
			else if (!deleted) *j++ = *i;
			if (deleted) {
				assert(c.size() > 1);
				saveProofClause(saved, c, PROOF_DELETED);
			}
		}
		ol.resize(j - ol);
	}

} 

#endif