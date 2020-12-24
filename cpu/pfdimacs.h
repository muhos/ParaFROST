/***********************************************************************[pfdimacs.h]
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

#ifndef __DIMACS_
#define __DIMACS_

#include "pfsort.h"
#include "pfdefs.h"
#include <fcntl.h>
#include <sys/stat.h>

namespace pFROST {

	inline bool isDigit(const char& ch) { return (ch ^ '0') <= 9; }

	inline void eatWS(char*& str) { while ((*str >= 9 && *str <= 13) || *str == 32) str++; }

	inline void eatLine(char*& str) { while (*str) if (*str++ == '\n') return; }

	inline uint32 toInteger(char*& str, uint32& sign)
	{
		eatWS(str);
		sign = 0;
		if (*str == '-') sign = 1, str++;
		else if (*str == '+') str++;
		if (!isDigit(*str)) PFLOGE("expected a digit but ASCII(%d) is found", *str);
		uint32 n = 0;
		while (isDigit(*str)) n = n * 10 + (*str++ - '0');
		return n;
	}

}

#endif 