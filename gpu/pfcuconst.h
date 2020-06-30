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

	namespace SIGmA {

		// Global
		#define MASTER_GPU		0
		#define AWAKEN_SUCC		0
		#define AWAKEN_FAIL		1
		#define CNFALLOC_FAIL	2
		#define OTALLOC_FAIL	3
		#define LCVE_FAIL		4
		#define BLOCK1D			256
		#define FULLWARP		0xFFFFFFFFU

		// BVE marking
		#define MELTING_MASK	0x80000000U
		#define UNMELT_MASK		0x7FFFFFFFU
		#define ELIMINATED(x)	(x & MELTING_MASK)
		#define RECOVERVAR(x)	(x & UNMELT_MASK)

		// Kernel configuration & GUARDs
		#define CE_POS_LMT		512
		#define CE_NEG_LMT		512
		#define FAN_LMT			64
		#define SH_MAX_BVE_OUT	125
		#define SH_MAX_BCE_IN	95 
		#define SH_MAX_HSE_IN	180
		#define SH_MAX_HRE_IN	100
		#define SH_MAX_HRE_OUT	250
		#define BLUB			256
		#define BLVE			64
		#define BLHSE			64
		#define BLBCE			128

	}
}
#endif