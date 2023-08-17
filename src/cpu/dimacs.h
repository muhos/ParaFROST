/***********************************************************************[dimacs.h]
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

#include "sort.h"
#include "definitions.h"
#include <fcntl.h>
#include <sys/stat.h>

namespace ParaFROST {

	template <class T>
	inline bool isDigit(const T& ch) { return (ch ^ 48) <= 9; }

	template <class T>
	inline bool isSpace(const T& ch) { return (ch >= 9 && ch <= 13) || ch == 32; }

	struct FORMULA {
		string path;
		double c2v;
		uint64 size;
		uint32 units, large, binaries, ternaries;
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

		inline int get() { size++; return std::cin.get(); }

		inline void eatWS(int& ch) { while (isSpace((ch = get()))); }

		inline void eatComment(int& ch) {
			while (true) {
				ch = get();
				if (isSpace(ch)) continue;
				if (ch != 'c') break;
				while ((ch = get()) != '\n')
					if (ch == EOF)
						PFLOGE("unexpected EOF in comment");
			}
		}

		inline uint32 toInteger(int& ch, uint32& sign)
		{
			sign = 0;
			if (ch == '-') {
				sign = 1;
				ch = get();
			}
			else if (ch == '+')
				ch = get();
			if (!isDigit(ch))
				PFLOGE("expected a digit but ASCII(%d) is found", ch);
			uint32 n = 0;
			while (isDigit(ch)) {
				n = n * 10 + (ch - 48);
				ch = get();
			}
			return n;
		}
	};

	inline bool isDigit(const char& ch) { return (ch ^ '0') <= 9; }

	inline void eatWS(char*& str) { while (isSpace(*str)) str++; }

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