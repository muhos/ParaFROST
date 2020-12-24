/***********************************************************************[pfsolve.cpp]
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

#include "pfsolve.h" 
#include "pfcudefs.h"

namespace pFROST {

	int64 sysMemUsed()
	{
		int64 memUsed = 0;
#ifdef __linux__ 
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
#elif _WIN32
		PROCESS_MEMORY_COUNTERS_EX memInfo;
		GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memInfo, sizeof(PROCESS_MEMORY_COUNTERS_EX));
		memUsed = memInfo.WorkingSetSize;
#endif
		return memUsed;
	}

	int64 getAvailSysMem()
	{
#ifdef __linux__ 
		long pages = sysconf(_SC_PHYS_PAGES);
		long page_size = sysconf(_SC_PAGE_SIZE);
		return pages * page_size;
#elif _WIN32
		MEMORYSTATUSEX memInfo;
		memInfo.dwLength = sizeof(MEMORYSTATUSEX);
		GlobalMemoryStatusEx(&memInfo);
		return memInfo.ullAvailPhys;
#endif
	}

	void set_timeout(int time_limit)
	{
#ifdef __linux__ 
		if (time_limit != 0) {
			rlimit rl;
			getrlimit(RLIMIT_CPU, &rl);
			if (rl.rlim_max == RLIM_INFINITY || (rlim_t)time_limit < rl.rlim_max) {
				rl.rlim_cur = time_limit;
				if (setrlimit(RLIMIT_CPU, &rl) == -1) PFLOGW("timeout cannot be set");
			}
		}
#elif _WIN32
		PFLOGW("timeout not supported on Windows");
#endif
	}

	void handler_terminate(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%s%45s%s", CYELLOW, "Interrupted", CNORMAL);
			PFLOG0("");
		}
		PFLOGS("UNKNOWN");
		if (!quiet_en) {
			PFLOG0("");
			PFLRULER('-', RULELEN);
		}
		_exit(EXIT_FAILURE);
	}

	void handler_mercy_interrupt(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%s%45s%s", CYELLOW, "Interrupted", CNORMAL);
			PFLOG0("");
		}
		pfrost->interrupt();
	}

	void handler_mercy_timeout(int)
	{
		fflush(stdout);
		if (!quiet_en) {
			PFLOG0("");
			PFLOG1("%s%45s%s", CYELLOW, "Timeout", CNORMAL);
			PFLOG0("");
		}
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

	void getCPUInfo(uint64& _free)
	{
		char cpuid[0x40];
#ifdef _WIN32
		int CPUInfo[4] = { -1 };
		__cpuid(CPUInfo, 0x80000000);
#else
		int CPUInfo[4] = { 0, 0, 0, 0 };
		__cpuid(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#endif
		uint32 nExIds = CPUInfo[0];
		for (uint32 i = 0x80000000; i <= nExIds; ++i) {
#ifdef _WIN32
			__cpuid(CPUInfo, i);
#else
			__cpuid(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#endif
			if (i == 0x80000002)
				memcpy(cpuid, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000003)
				memcpy(cpuid + 16, CPUInfo, sizeof(CPUInfo));
			else if (i == 0x80000004)
				memcpy(cpuid + 32, CPUInfo, sizeof(CPUInfo));
		}
		_free = getAvailSysMem();
		PFLOG2(1, " Available CPU: %s%s%s", CREPORTVAL, cpuid, CNORMAL);
		PFLOG2(1, " Available system memory: %lld GB", _free / GBYTE);
	}

	inline int	SM2Cores(int major, int minor) {
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
		PFLOGEN("MapSMtoCores for SM %d.%d is undefined. Default to use %d Cores/SM", major, minor, nCores[index - 1].Cores);
		return nCores[index - 1].Cores;
	}

	int getGPUInfo(size_t& _free, size_t& _penalty)
	{
		int devCount = 0;
		SIGmA::CHECK(cudaGetDeviceCount(&devCount));
		if (!devCount) return 0;
		SIGmA::CHECK(cudaGetDeviceProperties(&devProp, MASTER_GPU));
		assert(devProp.totalGlobalMem);
		if (devProp.warpSize != 32) PFLOGE("GPU warp size not supported");
#ifdef __linux__ 
		_penalty = 220 * MBYTE;
#elif _WIN32
		_penalty = 300 * MBYTE;
#endif
		_free = devProp.totalGlobalMem;
		maxGPUTereads = devProp.multiProcessorCount * devProp.maxThreadsPerMultiProcessor;
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
	wrapup();
	PFLOG0("");
	PFLOGN2(1, " Cleaning up..");
	this->~ParaFROST();
	pfrost = NULL;
	PFLDONE(1, 5);
	if (!quiet_en) PFLRULER('-', RULELEN);
	exit(EXIT_SUCCESS);
}