/***********************************************************************[pfcuconst.h]
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

namespace pFROST {

	// Global device-specific switches
	extern bool unified_access;
	extern bool profile_gpu;
	extern bool sync_always;
	extern bool atomic_ve;
	extern bool gc_par;

	namespace SIGmA {

		// Global
		#define MASTER_GPU		0
		#define AWAKEN_SUCC		0
		#define AWAKEN_FAIL		1
		#define CNFALLOC_FAIL	2
		#define OTALLOC_FAIL	3

		// BVE marking
		#define MELTING_MASK	0x80000000U
		#define UNMELT_MASK		0x7FFFFFFFU
		#define RES_MASK		0x40000000U
		#define AOIX_MASK		0x80000000U
		#define TYPE_MASK		0xC0000000U
		#define ADDED_MASK		0x3FFFFFFFU
		#define IS_RES(x)		(x & RES_MASK)
		#define IS_AOIX(x)		(x & AOIX_MASK)
		#define ELIMINATED(x)	(x & MELTING_MASK)
		#define RECOVERVAR(x)	(x & UNMELT_MASK)
		#define RECOVERADDED(x)	(x & ADDED_MASK)
		#define RECOVERTYPE(x)	(x & TYPE_MASK)


		// Kernel configuration & GUARDs
		#define BLOCK1D			256
		#define CE_POS_LMT		512
		#define CE_NEG_LMT		512
		#define BLVE			64
		#define BLVE1			64
		#define BLVE2			128
		#define SH_MAX_BVE_OUT	190
		#define SH_MAX_BVE_OUT1	190
		#define SH_MAX_BVE_OUT2	95
		#define BLBCE			128
		#define SH_MAX_BCE_IN	95 
		#define BLHSE			64
		#define HSE_MAX_CL_SIZE 1000
		#define SH_MAX_HSE_IN	190
		#define SH_MAX_HRE_OUT	350

	}
}
#endif