/***********************************************************************[constants.hpp]
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

#ifndef __GL_MACROS_
#define __GL_MACROS_

#include "datatypes.hpp"

extern bool quiet_en;
extern int verbose;

namespace ParaFROST {

	#define MBYTE			0x00100000
	#define KBYTE			0x00000400
	#define GBYTE			0x40000000
	#define NOREF			UINT64_MAX
	#define GNOREF			UINT64_MAX
	#define NOVAR			UINT32_MAX
	#define INIT_CAP		32
	#define UNDEFINED		-1
	#define ORGPHASE		1
	#define INVPHASE		2
	#define FLIPPHASE		3
	#define BESTPHASE		4
	#define RANDPHASE		5
	#define WALKPHASE		6
	#define UNSAT			0
	#define SAT				1
	#define UNSOLVED_M		2
	#define PROOF_ADDED	    97  // 'a'
	#define PROOF_DELETED	100 // 'd'
	#define UNSOLVED(x)		((x) & UNSOLVED_M)
	#define RESETSTATE(x)	(x = UNSOLVED_M)

	#define NEG_SIGN		0x00000001
	#define HASH_MASK		0x0000001F
	#define MAX_DLC			0x00000003
	#define MAX_LBD			0x04000000UL
	#define MAX_LBD_M		0x03FFFFFFUL
	#define BYTEMAX			0x00000080UL
	#define BYTEMASK		0x0000007FUL
	#define IBYTEMAX		0xFFFFFF80UL
	#define NOVAL_MASK		(LIT_ST)-2
	#define VAL_MASK		(LIT_ST) 1
	#define MELTED_M		(LIT_ST)0x01
	#define FROZEN_M		(LIT_ST)0x02
	#define SUBSTITUTED_M	(LIT_ST)0x04
	#define ANALYZED_M		(LIT_ST)0x01
	#define REMOVABLE_M		(LIT_ST)0x02
	#define POISONED_M		(LIT_ST)0x04
	#define KEEP_M			(LIT_ST)0x08
	#define USAGET3			(CL_ST)0x01
	#define USAGET2			(CL_ST)0x02
	#define USAGET1			(CL_ST)0x03
	#define ORIGINAL		(CL_ST)0x00
	#define LEARNT			(CL_ST)0x01
	#define DELETED			(CL_ST)0x02
	#define POS(x)			((x) & 0xFFFFFFFE)
	#define ABS(x)			((x) >> 1)
	#define V2L(x)			((x) << 1)
	#define SIGN(x)			(LIT_ST)((x) & NEG_SIGN)
	#define NEG(x)			((x) | NEG_SIGN)
	#define V2DEC(x,s)		(V2L(x) | (s))
	#define ISLARGE(x)		((x) & IBYTEMAX)
	#define L2B(x)			(((x) & BYTEMASK) | BYTEMAX)
	#define FLIP(x)			((x) ^ NEG_SIGN)
	#define HASH(x)			((x) & HASH_MASK)
	#define MAPHASH(x)		(1UL << HASH(x))
	#define MELTED(x)		((x) & MELTED_M)
	#define FROZEN(x)		((x) & FROZEN_M)
	#define SUBSTITUTED(x)	((x) & SUBSTITUTED_M)
	#define ANALYZED(x)		((x) & ANALYZED_M)
	#define	REMOVABLE(x)	((x) & REMOVABLE_M)	
	#define	POISONED(x)		((x) & POISONED_M)
	#define	KEPT(x)			((x) & KEEP_M)
	#define UNASSIGNED(x)	((x) & NOVAL_MASK)
	#define REASON(x)		((x) ^ NOREF)
	#define DECISION(x)		(!REASON(x))
	#define NEQUAL(x,y)		((x) ^ (y))
	#define MIN(x,y)		((x) < (y) ? (x) : (y))
	#define MAX(x,y)		((x) > (y) ? (x) : (y))
	#define RESETSTRUCT(MEMPTR) \
		memset(MEMPTR, 0, sizeof(*MEMPTR));
	
}

#endif