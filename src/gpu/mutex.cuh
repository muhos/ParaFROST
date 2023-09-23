/***********************************************************************[mutex.cuh]
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

#ifndef __MUTEX_
#define __MUTEX_

#include "definitions.cuh"

namespace ParaFROST {

	__managed__ int mutex = 0;

	_PFROST_D_ void lock(int* mutex) {
		while (atomicCAS(mutex, 0, 1) != 0);
	}

	_PFROST_D_ void unlock(int* mutex) {
		atomicExch(mutex, 0);
	}

	__global__ void reset_mutex_k() { mutex = 0; }

	#define RESETMUTEX reset_mutex_k<<<1,1>>>();

}

#endif