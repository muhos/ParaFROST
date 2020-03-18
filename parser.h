/***********************************************************************
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
************************************************************************/
#ifndef __PARSER_
#define __PARSER_

#include "Sort.h"
#include "defs.h"

inline int count_spaces(char* str)
{
	int spaces = 0;
	while (*str) if (*str++ != ' ' && *str == ' ') spaces++;
	return spaces;
}

inline uint32 toLiteral(char* str)
{
	uint32 v = 0, sign = 0;
	if (*str == '-') { sign = 1; str++; }
	while (*str) v = v * 10 + (*str++ - '0');
	if (sign && !v) {
		cout << "c | Error --> expected a digit after (-) sign." << endl;
		exit(EXIT_FAILURE);
	}
	if (cnf_stats.max_org_vars < v) cnf_stats.max_org_vars = v;
	return V2D(v) | sign;
}

inline uint32 toInteger(char* str)
{
	uint32 number = 0;
	while (*str) number = number * 10 + (*str++ - '0');
	return number;
}

inline void read_header(uVector1D& header, char* tmp, char* line)
{
	char digits = 0, numbers = 0;
	while (*line)
	{
		if (isdigit(*line)) *(tmp + digits++) = *line;
		if (isspace(*line) && digits > 0)
		{
			assert(digits <= CHAR_MAX && digits <= LIT_LEN);
			*(tmp + digits) = '\0';
			header[numbers++] = toInteger(tmp);
			digits = 0;
		}
		line++;
	}
	if (digits > 0) { *(tmp + digits) = '\0'; header[numbers++] = toInteger(tmp); }
	assert(numbers == 2);
}

inline void toLits(uint32* clause, char* tmp, char* line, int& len)
{
	char digits = 0;
	int clause_len = 0;
	while (*line)
	{
		if (isdigit(*line) || *line == '-') *(tmp + digits++) = *line;
		if (isspace(*line) && digits > 0 && *tmp != '0')
		{
			assert(digits <= CHAR_MAX && digits <= LIT_LEN);
			*(tmp + digits) = '\0';
			clause[clause_len++] = toLiteral(tmp);
			digits = 0;
		}
		line++;
	}
	assert(clause_len > 0 && clause_len <= len);
	len = clause_len;
}

inline bool checkClause(uint32* clause, int& len)
{
	if (len == 1) return true;
	else if (len == 2) {
		if ((*(clause + 1) ^ *clause) == NEG_SIGN) // clause is a tautology
			return false;
		if (*(clause + 1) == *clause) len = 1;
		else if (*(clause + 1) <  *clause) Swap(clause + 1, clause);
		assert(clause[0] <= clause[1]);
		return true;
	}
	Sort(clause, len);
	// check if duplicates exist before doing rewriting
	int clause_len = 1;
	for (int l = 1; l < len; l++) {
		if ((clause[l - 1] ^ clause[l]) == NEG_SIGN) // clause is a tautology
			return false;
		if (clause[l - 1] != clause[l])
			clause[clause_len++] = clause[l];
	}
	len = clause_len;
	return true;
}

#endif // !__PARSER_