/***********************************************************************[memory.cuh]
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

#ifndef __CU_MEMORY_
#define __CU_MEMORY_

#include "cnf.cuh"
#include "table.cuh"
#include "constants.cuh"
#include "variables.cuh"
#include "histogram.cuh"
#include "vstate.hpp"
#include "timer.cuh"
#include <cuarena/allocator.cuh>

namespace ParaFROST {
	/*****************************************************/
	/*  Usage:    GPU global memory manager              */
	/*  Dependency: all simplifier data types            */
	/*****************************************************/

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

	struct cuPool {
		addr_t mem;
		size_t cap;
		cuPool() : mem(NULL), cap(0) {}
	};

	struct cuLits {
		uint32* mem;
		size_t  size, cap;
		cuLits() : mem(NULL), size(0), cap(0) {}
	};

	class cuMM {

		// pools & arrays
		cuLits  litsPool;
		cuPool	hcnfPool,	cnfPool,	otPool;
		cuPool	varsPool,	histPool,	auxPool;
		cuPool	pinnedPool;
		CNF		* pinned_cnf;
		S_REF	* d_refs_mem,	* d_scatter,	* d_segs,	* d_occurs;
		uint32	* d_hist,		* d_cnf_mem;
		Byte	* d_stencil;
		VSTATE	* d_vstate;
		size_t	nscatters;
		// trackers
		cuTIMER cutimer;
		float 	_compacttime;
		int64	cap, dcap, penalty;
		// cuarena backend
		cuArena::DeviceArena arena;

	public:

		inline size_t 	getFreeMemory	() {
			size_t free = 0, tot = 0;
			CHECK(cudaMemGetInfo(&free, &tot));
			return free;
		}
		inline bool		hasDeviceMem	(const size_t& min_cap, const char* name,
									     const cuArena::Region& type = cuArena::Region::Stable) {
			const int64 used = dcap + min_cap;
			const size_t avail = (type == cuArena::Region::Stable) ?
				arena.gpu_stable_available() : arena.gpu_available();
			LOG2(2, " Allocating GPU-only memory for %s (used/free = %.2f/%.2f MB)", name,
				double(used) / MBYTE, double(avail) / MBYTE);
			if (min_cap > avail) { LOGWARNING("not enough memory for %s (current = %lld MB) - skip GPU simplifier", name, used / MBYTE); return false; }
			dcap = used;
			return true;
		}
		// Roll back a reservation made by hasDeviceMem when the subsequent
		// arena.allocate fails, so dcap stays consistent (no leak).
		inline void		undoDeviceMem	(const size_t& min_cap) {
			dcap  -= min_cap;
			assert(dcap >= 0);
		}
		// Frees a single device pool back to the arena.
		template <class POOL>
		inline void		freeDevice		(POOL& pool) {
			if (pool.mem) {
				assert(pool.cap);
				arena.deallocate(pool.mem);
				pool.mem = NULL;
				dcap -= pool.cap;
				assert(dcap >= 0);
				pool.cap = 0;
			}
		}
		// Device dynamic block (e.g. proof stream): arena Dynamic with compact-and-retry.
		// Freed via freeDevice.
		template <class POOL>
		inline bool		allocDynamic	(POOL& pool, const size_t& min_cap, const char* name) {
			if (!hasDeviceMem(min_cap, name, cuArena::Region::Dynamic)) return false;
			try { pool.mem = (addr_t)arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
			catch (const cuArena::gpu_memory_error&) {
				arena.compact_gpu_dynamic(nullptr);
				try { pool.mem = (addr_t)arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
				catch (const cuArena::gpu_memory_error& e2) {
					LOGWARNING("cuArena: %s alloc failed after compact (%s)", name, e2.what());
					undoDeviceMem(min_cap); return false;
				}
			}
			pool.cap = min_cap;
			return true;
		}
		// Pinned host block from the arena CPU pool.
		template <class POOL>
		inline bool		allocPinned		(POOL& pool, const size_t& min_cap) {
			try { pool.mem = (addr_t)arena.allocate_pinned<Byte>(min_cap); }
			catch (const cuArena::cpu_memory_error& e) {
				LOGWARNING("cuArena: pinned block alloc failed (%s)", e.what()); return false;
			}
			pool.cap = min_cap;
			return true;
		}
		template <class POOL>
		inline void		freePinned		(POOL& pool) {
			if (!pool.mem) return;
			arena.deallocate_pinned(pool.mem);
			pool.mem = NULL;
			pool.cap = 0;
		}
		inline cuArena::DeviceArena* deviceArena () {
			return &arena;
		}
						cuMM			();
		void			breakMirror		();
		void			freePinned		();
		void			freeDevice		();
		bool			initDeviceArena	(const size_t& numCls, const size_t& numLits, const size_t& resolvedCap, const bool& proofEnabled);
		inline bool		empty			() const { return cap == 0; }
		inline size_t	scatterCap		() const { return nscatters * sizeof(S_REF); }
		inline float 	compactTime		() const { return _compacttime; }
		inline size_t	pagedUsed		() const { return hcnfPool.cap; }
		inline S_REF*	refsMem			() { return d_refs_mem; }
		inline S_REF*	scatter			() { return d_scatter; }
		inline S_REF*	occurs			() { return d_occurs; }
		inline uint32*	cnfMem			() { return d_cnf_mem; }
		inline CNF*		pinnedCNF		() { return pinned_cnf; }
		inline VSTATE*	deviceVstate	() { return d_vstate; }
		inline cuLits&	literals		() { return litsPool; }
		inline void		init			(const int64 _penalty) { penalty = _penalty; }
		inline void		cacheCNFPtr		(const CNF* d_cnf, const cudaStream_t& _s = 0) {
			CHECK(cudaMemcpyAsync(pinned_cnf, d_cnf, sizeof(CNF), cudaMemcpyDeviceToHost, _s));
		}
		void			cuMemSetAsync	(addr_t, const Byte&, const size_t&);
		void	        resizeCNFAsync	(CNF*, const S_REF&, const uint32&);
		bool			resizeCNF		(CNF*&, const size_t&, const size_t&);
		uint32*			resizeLits		(const size_t&);
		bool			allocAux		(const size_t&);
		bool			allocPinned		(VARS*, cuHist&);
		bool			allocHist		(cuHist&, const bool& proofEnabled);
		bool			allocVars		(VARS*&, const size_t& resolvedCap);
		bool			resizeOTAsync	(OT*&, const size_t&, const cudaStream_t& _s = 0);
		void			createMirror	(CNF*&, const size_t&, const size_t&);
		void			mirrorCNF		(CNF*&);
		void			compactCNF		(CNF*, CNF*);

	};

	

#if defined(__GNUC__) && (__GNUC__ >= 8)
#pragma GCC diagnostic pop
#endif

}

#endif
