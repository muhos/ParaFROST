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
#include "pfdefinitions.h"
#include <fcntl.h>
#include <sys/stat.h>

namespace pFROST {

	struct FORMULA {
		string path;
		uint64 size;
		double c2v;
		uint32 units, binaries, ternaries, large;
		int maxClauseSize;
		FORMULA() : 
			path()
			, c2v(0)
			, size(0)
			, units(0)
			, large(0)
			, binaries(0)
			, ternaries(0)
			, maxClauseSize(0) {}
		FORMULA(const string& path) :
			path(path)
			, c2v(0)
			, size(0)
			, units(0)
			, large(0)
			, binaries(0)
			, ternaries(0)
			, maxClauseSize(0) {}
	};

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

	inline bool canAccess(const char* path, struct stat& st)
	{
		if (stat(path, &st)) return false;
#ifdef _WIN32
#define R_OK 4
		if (_access(path, R_OK)) return false;
#else
		if (access(path, R_OK)) return false;
#endif
		return true;
	}

}

#endif 