/***********************************************************************[constants.cuh]
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

#ifndef __GL_CU_MACROS_
#define __GL_CU_MACROS_

namespace ParaFROST {

	extern size_t hc_bucket;
	extern size_t hc_nbuckets;

	// Global
	#define MASTER_GPU		0
	#define AWAKEN_SUCC		0
	#define AWAKEN_FAIL		1
	#define CNFALLOC_FAIL	2
	#define OTALLOC_FAIL	3
	#define NLIMITS			7

	// BVE bit-masking
	#define MAXFUNVAR			12
	#define FUNTABLEN			64
	#define MELTING_MASK		(Byte) 0x01
	#define ADDING_MASK	    	(Byte) 0x02
	#define RES_MASK			0x00000001U
	#define AOIX_MASK			0x00000002U
	#define CORE_MASK			0x00000003U
	#define TYPE_MASK			0x00000003U
	#define ADDEDCLS_MASK		0x0000FFFCU
	#define ADDEDLITS_MASK		0xFFFF0000U
	#define ADDEDCLS_MAX		0x00003FFFU
	#define ADDEDLITS_MAX		0x0000FFFFU
	#define ADDEDPROOF_MAX		0x0003FFFFU
	#define ELIMINATE(x)					(x |= MELTING_MASK)
	#define MARKADDING(x)					(x |= ADDING_MASK)
	#define IS_TAUTOLOGY(x,y)				(((x) ^ (y)) == NEG_SIGN)	
	#define IS_RES(x)						((x) == RES_MASK)
	#define IS_AOIX(x)						((x) == AOIX_MASK)
	#define IS_CORE(x)						((x) == CORE_MASK)
	#define IS_ADDING(x)					((x) & ADDING_MASK)
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
	#define CE_POS_LMT		512
	#define CE_NEG_LMT		512
	#define BLSORT			256 
	#define BLVE			64
	#define BLVE1			64
	#define BLVE2			128
	#define BLSUB			64
	#define BLBCE			128
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