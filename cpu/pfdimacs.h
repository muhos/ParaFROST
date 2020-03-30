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

inline bool isDigit(const char& ch) { return (ch ^ '0') <= 9; }

inline void skipWS(char*& str) { while ((*str >= 9 && *str <= 13) || *str == 32) str++; }

static void skipLine(char*& str) { while (*str) if (*str++ == '\n') return; }

inline int toInteger(char*& str)
{
	skipWS(str);
	int n = 0;
	bool sign = false;
	if (*str == '-') sign = true, str++;
	else if (*str == '+') str++;
	if (!isDigit(*str)) { printf("Error - expected a digit but ASCII(%d) is found\n", *str), exit(EXIT_FAILURE); }
	while (isDigit(*str)) n = n * 10 + (*str++ - '0');
	return sign ? -n : n;
}

inline void toClause(uVector1D& c, char*& str)
{
	c.clear();
	int v = 0;
	while ((v = toInteger(str)) != 0) {
		uint32 abs_v = abs(v);
		if (abs_v > nOrgVars()) { printf("Error - too many variables\n"), exit(EXIT_FAILURE); }
		c.push((v > 0) ? V2D(abs_v) : NEG(V2D(abs_v)));
	}
}

inline bool checkClause(uVector1D& c)
{
	if (c.size() == 1) return true;
	Sort(c.d_ptr(), c.size());
	// check if duplicates exist before doing rewriting
	int c_sz = 1;
	for (int l = 1; l < c.size(); l++) {
		if ((c[l - 1] ^ c[l]) == NEG_SIGN) return false; // c is a tautology
		if (c[l - 1] != c[l]) c[c_sz++] = c[l];
	}
	c.resize(c_sz);
	return true;
}

#endif 