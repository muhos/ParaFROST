/***********************************************************************[banner.cpp]
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

#include "banner.hpp"

void LOGFANCYBANNER(const char* VER) 
{
#if defined(_WIN32)
    SetConsoleOutputCP(65001);
#endif
    const char* solvername1 = u8"██████╗░░█████╗░██████╗░░█████╗░███████╗██████╗░░█████╗░░██████╗████████╗";
    const char* solvername2 = u8"██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝";
    const char* solvername3 = u8"██████╔╝███████║██████╔╝███████║█████╗░░██████╔╝██║░░██║╚█████╗░░░░██║░░░";
    const char* solvername4 = u8"██╔═══╝░██╔══██║██╔══██╗██╔══██║██╔══╝░░██╔══██╗██║░░██║░╚═══██╗░░░██║░░░";
    const char* solvername5 = u8"██║░░░░░██║░░██║██║░░██║██║░░██║██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░██║░░░";
    const char* solvername6 = u8"╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░╚═╝░░░";
    const char* name_version = u8"Copyright\u00A9 Muhammad Osama Mahmoud                                  ";
    size_t len = 75 + 1;
    if (RULELEN < len) LOGERROR("ruler length is smaller than the title (%zd)", len);
    size_t gap = (RULELEN - len - 3) / 2;
    PRINT(PREFIX);
    PUTCH('\n');
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername1, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername2, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername3, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername4, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername5, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%s\n", CSOLVER, solvername6, CNORMAL);
    PRINT(PREFIX);
    REPCH(' ', gap);
    PRINT("%s%s%6s%s\n", CSOLVER, name_version, VER, CNORMAL);
    PRINT(PREFIX);
    PUTCH('\n');
}