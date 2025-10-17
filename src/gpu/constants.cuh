/***********************************************************************[constants.cuh]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
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

#ifndef __GL_CU_MACROS_
#define __GL_CU_MACROS_

#include "datatypes.hpp"

namespace ParaFROST {

	// Global
	constexpr int MASTER_GPU 	= 0;
	constexpr int AWAKEN_SUCC 	= 0;
	constexpr int AWAKEN_FAIL 	= 1;
	constexpr int CNFALLOC_FAIL = 2;
	constexpr int OTALLOC_FAIL 	= 3;

	// BVE bit-masking
	constexpr Byte	 MELTING_MASK	= 0x01;
	constexpr Byte   ADDING_MASK	= 0x02;
	constexpr Byte   FORCED_MASK	= 0x04;
	constexpr Byte   MAX_ELIM_MASK	= FORCED_MASK;
	constexpr uint32 RES_MASK		= 0x00000001U;
	constexpr uint32 AOIX_MASK		= 0x00000002U;
	constexpr uint32 CORE_MASK		= 0x00000003U;
	constexpr uint32 TYPE_MASK		= 0x00000003U;
	constexpr uint32 ADDEDCLS_MASK	= 0x0000FFFCU;
	constexpr uint32 ADDEDLITS_MASK = 0xFFFF0000U;
	constexpr uint32 ADDEDCLS_MAX	= 0x00003FFFU;
	constexpr uint32 ADDEDLITS_MAX	= 0x0000FFFFU;
	constexpr uint32 ADDEDPROOF_MAX = 0x0003FFFFU;

	#define ELIMINATE(x)					(x |= MELTING_MASK)
	#define MARKADDING(x)					(x |= ADDING_MASK)
	#define MARKFORCED(x)					(x |= FORCED_MASK)
	#define IS_TAUTOLOGY(x,y)				(((x) ^ (y)) == NEG_SIGN)	
	#define IS_RES(x)						((x) == RES_MASK)
	#define IS_AOIX(x)						((x) == AOIX_MASK)
	#define IS_CORE(x)						((x) == CORE_MASK)
	#define IS_ADDING(x)					((x) & ADDING_MASK)
	#define IS_FORCED(x)					((x) & FORCED_MASK)
	#define ELIMINATED(x)					((x) & MELTING_MASK)
	#define RECOVERTYPE(x)					((x) & TYPE_MASK)
	#define RECOVERADDEDCLS(x)				(((x) & ADDEDCLS_MASK) >> 2)
	#define RECOVERADDEDLITS(x)				(((x) & ADDEDLITS_MASK) >> 16)
	#define RECOVERADDEDUNITS(x)			((x) & ADDEDCLS_MAX)
	#define RECOVERADDEDPROOF(x)			((x) >> 14)
	#define ENCODEVARINFO(T, CLS, LITS)		((T) | ((CLS) << 2) | ((LITS) << 16))
	#define ENCODEPROOFINFO(UNITS, BYTES)	((UNITS) | ((BYTES) << 14))

	// Kernel configuration & GUARDs
	#define LBD_TIER1		2
	#define BLOCK1D			256
	#define BLVE1			64
	#define BLVE2			128
	#define BLSUB			64
	#define BLBCE			128
	#define BLEREY			8
	#define BLEREX			64
#if defined(EXTSHMEM)
	#define SH_MAX_BVE_OUT1	250 
	#define SH_MAX_BVE_OUT2	120
	#define SH_MAX_SUB_IN	250
	#define SH_MAX_ERE_OUT	250
#else
	#define SH_MAX_BVE_OUT1	190
	#define SH_MAX_BVE_OUT2	95
	#define SH_MAX_SUB_IN	190
	#define SH_MAX_ERE_OUT	190
#endif
	#define SH_MAX_BVE_OUT	190
	#define SH_MAX_BCE_IN	95 
	#define SUB_MAX_CL_SIZE 1000
	
}
#endif