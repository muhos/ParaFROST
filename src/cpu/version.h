/***********************************************************************[version.hpp]
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

namespace ParaFROST {

	const char* version();
	const char* signature();
	const char* compiler();
	const char* compilemode();
	const char* osystem();
	const char* date();

}
#define VERSION "3.2.5"
#define OSYSTEM "linux omega 6.8.0-41-generic x86_64"
#define DATE "Fri Aug 30 04:00:26 PM CEST 2024"
