/***********************************************************************[version.h]
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

namespace pFROST {

	const char* version();
	const char* compiler();
	const char* compilemode();
	const char* osystem();
	const char* date();

}
#define VERSION "2.1.2"
#define COMPILER "g++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0"
#define OSYSTEM "linux nbwin1583 5.4.72-microsoft-standard-wsl2 x86_64"
#define DATE "Fri May 28 14:52:08 CEST 2021"
