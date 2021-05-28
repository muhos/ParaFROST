/***********************************************************************[solve.cpp]
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

#include "solve.h" 
#include "definitions.cuh"
#include "version.h"

#if defined(__linux__)
#include <fpu_control.h>
#endif

namespace pFROST {

	int64 sysMemUsed()
	{
		int64 memUsed = 0;
#if defined(__linux__) || defined(__CYGWIN__)
		FILE* file = fopen("/proc/self/status", "r");
		char line[128];
		uint32 sign = 0;
		while (fgets(line, 128, file) != NULL) {
			char* str = line;
			if (eq(str, "VmRSS:")) {
				eatWS(str);
				memUsed = toInteger(str, sign);
				break;
			}
		}
		fclose(file);
		return memUsed * KBYTE;
#elif defined(_WIN32)
		PROCESS_MEMORY_COUNTERS_EX memInfo;
		GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memInfo, sizeof(PROCESS_MEMORY_COUNTERS_EX));
		memUsed = memInfo.WorkingSetSize;
#endif
		return memUsed;
	}

	int64 getAvailSysMem()
	{
#if defined(__linux__) || defined(__CYGWIN__)
		long pages = sysconf(_SC_PHYS_PAGES);
		long page_size = sysconf(_SC_PAGE_SIZE);
		return pages * page_size;
#elif defined(_WIN32)
		MEMORYSTATUSEX memInfo;
		memInfo.dwLength = sizeof(MEMORYSTATUSEX);
		GlobalMemoryStatusEx(&memInfo);
		return memInfo.ullAvailPhys;
#endif
	}

	void set_timeout(int time_limit)
	{
#if defined(__linux__) || defined(__CYGWIN__)
		if (time_limit) {
			rlimit limit;
			getrlimit(RLIMIT_CPU, &limit);
			if (limit.rlim_max == RLIM_INFINITY || (rlim_t)time_limit < limit.rlim_max) {
				limit.rlim_cur = time_limit;
				if (setrlimit(RLIMIT_CPU, &limit) == -1) PFLOGW("timeout cannot be set");
			}
		}
#elif defined(_WIN32)
		PFLOGW("timeout not supported on Windows");
#endif
	}

	void set_memoryout(int memory_limit)
	{
#if defined(__linux__)
		if (memory_limit) {
			rlim64_t limitbytes = (rlim64_t)memory_limit * GBYTE;
			rlimit64 limit;
			getrlimit64(RLIMIT_AS, &limit);
			if (limit.rlim_max == RLIM_INFINITY || limitbytes < limit.rlim_max) {
				limit.rlim_cur = limitbytes;
				if (setrlimit64(RLIMIT_AS, &limit) == -1) PFLOGW("memoryout cannot be set");
			}
		}
#elif defined(__CYGWIN__)
		if (memory_limit) {
			rlim_t limitbytes = (rlim_t)memory_limit * GBYTE;
			rlimit limit;
			getrlimit(RLIMIT_AS, &limit);
			if (limit.rlim_max == RLIM_INFINITY || limitbytes < limit.rlim_max) {
				limit.rlim_cur = limitbytes;
				if (setrlimit(RLIMIT_AS, &limit) == -1) PFLOGW("memoryout cannot be set");
			}
		}
#elif defined(_WIN32)
		PFLOGW("memoryout not supported on Windows");
#endif
			}

	void forceFPU()
	{
#if defined(__linux__) && defined(_FPU_EXTENDED) && defined(_FPU_DOUBLE)
		fpu_control_t cw = (_FPU_DEFAULT & ~_FPU_EXTENDED) | _FPU_DOUBLE;
		_FPU_SETCW(cw);
		PFLOG2(1, " Forced FPU to use %sdouble precision%s", CREPORTVAL, CNORMAL);
#endif 
	}

	void handler_terminate(int)
	{
		fflush(stdout);
		PFLOG1("%s%s%s", CYELLOW, "INTERRUPTED", CNORMAL);
		PFLOGS("UNKNOWN");
		_exit(EXIT_FAILURE);
	}

	void handler_mercy_interrupt(int)
	{
		fflush(stdout);
		PFLOG1("%s%s%s", CYELLOW, "INTERRUPTED", CNORMAL);
		pfrost->interrupt();
	}

	void handler_mercy_timeout(int)
	{
		fflush(stdout);
		PFLOG1("%s%s%s", CYELLOW, "TIME OUT", CNORMAL);
		PFLOGS("UNKNOWN");
		pfrost->interrupt();
	}

	void signal_handler(void h_intr(int), void h_timeout(int))
	{
		signal(SIGINT, h_intr);
		signal(SIGTERM, h_intr);
#ifdef SIGXCPU
		if (h_timeout != NULL) signal(SIGXCPU, h_timeout);
#endif
	}

	void segmentation_fault(int)
	{
		fflush(stdout);
		PFLOGEN("segmentation fault detected.");
		_exit(EXIT_FAILURE);
	}

	void illegal_code(int)
	{
		fflush(stdout);
		PFLOGEN("illegal code detected.");
		_exit(EXIT_FAILURE);
	}

	void arithmetic_error(int)
	{
		fflush(stdout);
		PFLOGEN("arithmetic flaw detected.");
		_exit(EXIT_FAILURE);
	}

	void getCPUInfo(uint64& _free)
	{
#ifndef __CYGWIN__
		char cpuid[0x40] = { 0 };
#if defined(_WIN32)
		int CPUInfo[4] = { -1 };
		__cpuid(CPUInfo, 0x80000000);
#elif defined(__linux__)
		int CPUInfo[4] = { 0, 0, 0, 0 };
		__cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#endif
		uint32 nExIds = CPUInfo[0];
		for (uint32 i = 0x80000000; i <= nExIds; ++i) {
#if defined(_WIN32)
			__cpuid(CPUInfo, i);
#elif defined(__linux__)
			__cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#endif
			if (i == 0x80000002)
				memcpy(cpuid, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000003)
				memcpy(cpuid + 16, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000004)
				memcpy(cpuid + 32, CPUInfo, sizeof(CPUInfo));
		}
		PFLOG2(1, " Available CPU: %s%s%s", CREPORTVAL, cpuid, CNORMAL);
#endif
		_free = getAvailSysMem();
		PFLOG2(1, " Available system memory: %lld GB", _free / GBYTE);
	}

	void getBuildInfo()
	{
		PFLOG1(" Built on %s%s%s at %s%s%s", CREPORTVAL, osystem(), CNORMAL, CREPORTVAL, date(), CNORMAL);
		PFLOG1("       using %s%s %s%s", CREPORTVAL, compiler(), compilemode(), CNORMAL);
	}

	inline 
	int SM2Cores(int major, int minor)
	{
		typedef struct { int SM; int Cores; } SM;

		SM nCores[] = { {0x30, 192}, {0x32, 192}, {0x35, 192}, {0x37, 192}, {0x50, 128}, {0x52, 128}, {0x53, 128},
			{0x60,  64}, {0x61, 128}, {0x62, 128}, {0x70,  64}, {0x72,  64}, {0x75,  64}, {-1, -1} };

		int index = 0;
		while (nCores[index].SM != -1) {
			if (nCores[index].SM == ((major << 4) + minor)) {
				return nCores[index].Cores;
			}
			index++;
		}
		PFLOGEN("cannot map to cores for SM %d.%d. Default to use %d Cores/SM", major, minor, nCores[index - 1].Cores);
		return nCores[index - 1].Cores;
	}

	int getGPUInfo(size_t& _free, size_t& _penalty)
	{
		int devCount = 0;
		CHECK(cudaGetDeviceCount(&devCount));
		if (!devCount) return 0;
		CHECK(cudaGetDeviceProperties(&devProp, MASTER_GPU));
		assert(devProp.totalGlobalMem);
		if (devProp.warpSize != 32) PFLOGE("GPU warp size not supported");
#if defined(__linux__) || defined(__CYGWIN__) 
		_penalty = 220 * MBYTE;
#elif defined(_WIN32)
		_penalty = 300 * MBYTE;
#endif
		_free = devProp.totalGlobalMem;
		maxGPUThreads = devProp.multiProcessorCount * devProp.maxThreadsPerMultiProcessor;
		maxGPUSharedMem = devProp.sharedMemPerBlock;
		if (!quiet_en) {
			PFLOG1(" Available GPU: %d x %s%s @ %.2fGHz%s (compute cap: %d.%d)", devCount, CREPORTVAL, devProp.name, ratio((double)devProp.clockRate, 1e6), CNORMAL, devProp.major, devProp.minor);
			PFLOG1(" Available GPU Multiprocessors: %d MPs (%d cores/MP)", devProp.multiProcessorCount, SM2Cores(devProp.major, devProp.minor));
			PFLOG1(" Available Global memory: %zd GB", _free / GBYTE);
			PFLRULER('-', RULELEN);
		}
		return devCount;
	}

}

void pFROST::ParaFROST::killSolver()
{
	syncAll();
	wrapup();
	PFLOG0("");
	PFLOGN2(1, " Cleaning up..");
	this->~ParaFROST();
	pfrost = NULL;
	PFLDONE(1, 5);
	if (!quiet_en) PFLRULER('-', RULELEN);
	exit(EXIT_SUCCESS);
}