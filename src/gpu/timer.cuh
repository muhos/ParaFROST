/***********************************************************************[timer.cuh]
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

#ifndef __GPU_TIMER_
#define __GPU_TIMER_

#include "definitions.cuh"
#include "constants.hpp"

namespace ParaFROST {

	class cuTIMER {
	private:
		cudaEvent_t _start, _stop;
		float _gpuTime;
	public:
		cuTIMER() :
			_start(NULL), _stop(NULL), _gpuTime(0) 
		{
			if (!_start) cudaEventCreate(&_start);
			if (!_stop) cudaEventCreate(&_stop);
		}
		~cuTIMER() {
			if (_start) cudaEventDestroy(_start);
			if (_stop) cudaEventDestroy(_stop);
			_start = NULL, _stop = NULL;
		}
		inline void  start  (const cudaStream_t& _s = 0) { 
			assert(_start && _stop);
			cudaEventRecord(_start, _s); 
		}
		inline void  stop   (const cudaStream_t& _s = 0) { 
			assert(_start && _stop);
			cudaEventRecord(_stop, _s); 
		}
		inline float gpuTime() {
			assert(_start && _stop);
			_gpuTime = 0;
			cudaEventSynchronize(_stop);
			cudaEventElapsedTime(&_gpuTime, _start, _stop);
			return _gpuTime;
		}
	};

}

#endif
