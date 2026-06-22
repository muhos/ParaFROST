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
		int64	_tot, _free, cap, dcap, maxcap, penalty;
		// cuarena backend
		cuArena::DeviceArena arena;
		bool                 arena_ready;

	public:

		inline bool		hasDeviceMem	(const size_t& min_cap, const char* name,
									     const cuArena::Region& type = cuArena::Region::Stable) {
			const int64 used = dcap + min_cap;
			if (arena_ready) {
				const size_t avail = (type == cuArena::Region::Stable) ?
					arena.gpu_stable_available() : arena.gpu_available();
				LOG2(2, " Allocating GPU-only memory for %s (used/free = %.2f/%.2f MB)", name,
					double(used) / MBYTE, double(avail) / MBYTE);
				if (min_cap > avail) { LOGWARNING("not enough memory for %s (current = %lld MB) - skip GPU simplifier", name, used / MBYTE); return false; }
				dcap = used;
				return true;
			}
			LOG2(2, " Allocating GPU-only memory for %s (used/free = %.2f/%lld MB)", name,
				double(used) / MBYTE, _free / MBYTE);
			if (int64(min_cap) >= _free) { LOGWARNING("not enough memory for %s (current = %lld MB) - skip GPU simplifier", name, used / MBYTE); return false; }
			_free -= min_cap;
			dcap = used;
			assert(_free >= 0);
			assert(_tot >= _free);
			assert(_tot > dcap);
			return true;
		}
		// Roll back a reservation made by hasDeviceMem when the subsequent
		// arena.allocate fails, so dcap/_free stay consistent (no leak).
		inline void		undoDeviceMem	(const size_t& min_cap) {
			_free += min_cap;
			dcap  -= min_cap;
			assert(dcap >= 0);
		}
		template <class POOL>
		inline void		DFREE			(POOL& pool, const bool& usearena = true) {
			if (pool.mem) {
				assert(pool.cap);
				if (arena_ready && usearena) arena.deallocate(pool.mem);
				else                           CHECK(cudaFree(pool.mem));
				pool.mem = NULL;
				dcap -= pool.cap;
				assert(dcap >= 0);
				_free += pool.cap;
				pool.cap = 0;
			}
		}
		// Device dynamic block (e.g. proof stream): arena Dynamic with compact-and-retry,
		// or raw cudaMalloc when the arena is unavailable. Freed via DFREE.
		template <class POOL>
		inline bool		allocDynamic	(POOL& pool, const size_t& min_cap, const char* name) {
			if (arena_ready) {
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
			if (!hasDeviceMem(min_cap, name)) return false;
			CHECK(cudaMalloc((void**)&pool.mem, min_cap));
			pool.cap = min_cap;
			return true;
		}
		// Frees a block obtained from allocDynamic, matching its build-config source.
		template <class POOL>
		inline void		freeDynamic		(POOL& pool) {
			DFREE(pool, true);
		}
		// Pinned host block: arena CPU pool when available, else cudaMallocHost.
		template <class POOL>
		inline bool		pinnedAlloc		(POOL& pool, const size_t& min_cap) {
			if (arena_ready && arena.cpu_capacity()) {
				try { pool.mem = (addr_t)arena.allocate_pinned<Byte>(min_cap); }
				catch (const cuArena::cpu_memory_error& e) {
					LOGWARNING("cuArena: pinned block alloc failed (%s)", e.what()); return false;
				}
				pool.cap = min_cap;
				return true;
			}
			if (cudaMallocHost((void**)&pool.mem, min_cap) != cudaSuccess) return false;
			pool.cap = min_cap;
			return true;
		}
		template <class POOL>
		inline void		pinnedFree		(POOL& pool) {
			if (!pool.mem) return;
			if (arena_ready && arena.cpu_capacity()) arena.deallocate_pinned(pool.mem);
			else CHECK(cudaFreeHost(pool.mem));
			pool.mem = NULL;
			pool.cap = 0;
		}
						cuMM			();
		bool			initDeviceArena	(const size_t& numCls, const size_t& numLits, const size_t& resolvedCap, const bool& proofEnabled);
		void			breakMirror		();
		void			freePinned		();
		void			freeFixed		();
		void			freeVars		();
		void			freeCNF			();
		void			freeOT			();
		inline void		updateMaxCap	() { if (maxcap < (cap + dcap)) maxcap = cap + dcap; }
		inline bool		empty			() const { return cap == 0; }
		inline int64	ucapacity		() const { return cap; }
		inline int64	dcapacity		() const { return dcap; }
		inline int64	maxCapacity		() const { return maxcap; }
		inline size_t	scatterCap		() const { return nscatters * sizeof(S_REF); }
		inline float 	compactTime		() const { return _compacttime; }
		inline S_REF*	refsMem			() { return d_refs_mem; }
		inline S_REF*	scatter			() { return d_scatter; }
		inline S_REF*	occurs			() { return d_occurs; }
		inline uint32*	cnfMem			() { return d_cnf_mem; }
		inline CNF*		pinnedCNF		() { return pinned_cnf; }
		inline VSTATE*	deviceVstate	() { return d_vstate; }
		inline cuArena::DeviceArena* deviceArena () { return arena_ready ? &arena : nullptr; }
		inline cuLits&	literals		() { return litsPool; }
		void			reset			(const int64 _free) {
			this->_free = _free - penalty;
		}
		void			init			(const int64 _tot, const int64 _penalty) {
			penalty = _penalty;
			this->_tot = _tot;
			_free = _tot - penalty;
		}
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
