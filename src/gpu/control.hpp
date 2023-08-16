/***********************************************************************[control.hpp]
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

#ifndef __CONTROL_
#define __CONTROL_

#include "datatypes.hpp"
#include "definitions.hpp"
#include "dimacs.hpp"

namespace ParaFROST {

void forceFPU();
int64 sysMemUsed();
int64 getAvailSysMem();
void getBuildInfo();
void getCPUInfo(uint64& _free);
int getGPUInfo(size_t& _free, size_t& _penalty);
void signal_handler(void h_intr(int), void h_timeout(int) = NULL);
void set_timeout(int);
void set_memoryout(int);
void handler_terminate(int);
void handler_mercy_interrupt(int);
void handler_mercy_timeout(int);
void segmentation_fault(int);
void illegal_code(int);
void arithmetic_error(int);

#define FAULT_DETECTOR                       \
    {                                        \
        signal(SIGSEGV, segmentation_fault); \
        signal(SIGILL, illegal_code);        \
        signal(SIGFPE, arithmetic_error);    \
    }
} // namespace ParaFROST

#endif