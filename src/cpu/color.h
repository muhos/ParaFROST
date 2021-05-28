/***********************************************************************[color.h]
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

#ifndef __COLOR_
#define __COLOR_

#include <cstdio>

#define CNORMAL     "\x1B[0m"
#define CRED        "\x1B[31m"
#define CGREEN      "\x1B[32m"
#define CYELLOW     "\x1B[33m"
#define CBLUE       "\x1B[34m"
#define CMAGENTA    "\x1B[35m"
#define CCYAN       "\x1B[36m"
#define CWHITE      "\x1B[37m"
#define CLRED       "\x1B[91m"
#define CLGREEN     "\x1B[92m"
#define CLYELLOW    "\x1B[93m"
#define CLBLUE      "\x1B[94m"
#define CLMAGENTA   "\x1B[95m"
#define CLCYAN      "\x1B[96m"
#define CGREEN0		"\x1B[38;5;78m"
#define CGREEN1		"\x1B[38;5;77m"
#define CGREEN2		"\x1B[38;5;76m"
#define CGREEN3		"\x1B[38;5;40m"
#define CORANGE0	"\x1B[38;5;203m"
#define CORANGE1	"\x1B[38;5;202m"
#define CVIOLET0	"\x1B[38;5;168m"
#define CVIOLET1	"\x1B[38;5;170m"
#define CVIOLET2	"\x1B[38;5;165m"
#define CVIOLET3	"\x1B[38;5;164m"
#define CVIOLET4	"\x1B[38;5;163m"
#define CVIOLET5	"\x1B[38;5;125m"

#define CSOLVER		"\x1B[38;5;117m"
#define CAUTHOR		"\x1B[38;5;116m"
#define CRIGHTS		"\x1B[38;5;1m"
#define CERROR		"\x1B[38;5;1m"
#define CWARNING	"\x1B[38;5;226m"
#define CHEADER		"\x1B[38;5;109m"
#define CREPORT		"\x1B[38;5;187m"
#define CREPORTVAL	"\x1B[38;5;106m"
#define CLOGGING	"\x1B[38;5;225m"
#define CCONFLICT	"\x1B[38;5;198m"
#define CHELP		"\x1B[38;5;108m"
#define CARGDEFAULT	"\x1B[38;5;187m"
#define CARGVALUE	"\x1B[38;5;34m"
#define CARGON		"\x1B[38;5;34m"
#define CARGOFF		"\x1B[38;5;124m"


#define SETCOLOR(COLOR, STD) fprintf(STD, "%s", COLOR); 

#endif