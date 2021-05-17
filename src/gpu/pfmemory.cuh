/***********************************************************************[pfmemory.cuh]
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

#ifndef __CU_MEMORY_
#define __CU_MEMORY_

#include "pfsimptypes.h"
#include <map>

namespace pFROST {

	namespace SIGmA {
		/*****************************************************/
		/*  Usage:    GPU global memory manager              */
		/*  Dependency: all simplifier data types            */
		/*****************************************************/
		class cuMM {
			// pools & arrays
			cuLits  litsPool;
			cuPool	hcnfPool,  cnfPool,  otPool,  varsPool;
			cuPool	hhistPool, histPool, auxPool, tmpPool;
			S_REF   * d_refs_mem, * d_scatter,   * d_segs;
			uint32  * d_hist,     * d_cnf_mem;
			uint32  * d_units,    * pinned_units;
			Byte    * d_stencil;
			CNF     * pinned_cnf;
			// trackers
			int64	_tot, _free, cap, dcap, maxcap, penalty;
			size_t	nscatters;
			inline bool		hasUnifiedMem(const size_t& min_cap, const char* name) {
				const int64 used = cap + min_cap;
				PFLOG2(2, " Allocating GPU unified memory for %s (used/free = %.3f/%lld MB)", name, double(used) / MBYTE, _free / MBYTE);
				if (used >= _free) { PFLOGW("not enough memory for %s (current = %lld MB) -> skip GPU simplifier", name, used / MBYTE); return false; }
				_free -= used;
				cap = used;
				assert(_free >= 0);
				assert(_tot >= _free);
				assert(_tot > cap);
				return true;
			}
			inline bool		hasDeviceMem(const size_t& min_cap, const char* name) {
				const int64 used = dcap + min_cap;
				PFLOG2(2, " Allocating GPU-only memory for %s (used/free = %.3f/%lld MB)", name, double(used) / MBYTE, _free / MBYTE);
				if (used >= _free) { PFLOGW("not enough memory for %s (current = %lld MB) -> skip GPU simplifier", name, used / MBYTE); return false; }
				_free -= used;
				dcap = used;
				assert(_free >= 0);
				assert(_tot >= _free);
				assert(_tot > dcap);
				return true;
			}
			template <class POOL>
			inline void		FREE(POOL& pool) {
				if (pool.mem) {
					assert(pool.cap);
					CHECK(cudaFree(pool.mem));
					pool.mem = NULL;
					cap -= pool.cap;
					assert(cap >= 0);
					_free += pool.cap;
					pool.cap = 0;
				}
			}
			template <class POOL>
			inline void		DFREE(POOL& pool) {
				if (pool.mem) {
					assert(pool.cap);
					CHECK(cudaFree(pool.mem));
					pool.mem = NULL;
					dcap -= pool.cap;
					assert(dcap >= 0);
					_free += pool.cap;
					pool.cap = 0;
				}
			}

		public:
			cuMM() {
				RESETSTRUCT(this);
			}
			void			breakMirror      ();
			void			freePinned       ();
			void			freeFixed        ();
			void			freeVars         ();
			void			freeCNF          ();
			void			freeOT           ();
			inline void		updateMaxCap     () { if (maxcap < (cap + dcap)) maxcap = cap + dcap; }
			inline bool		empty            () const { return cap == 0; }
			inline int64	ucapacity        () const { return cap; }
			inline int64	dcapacity        () const { return dcap; }
			inline int64	maxCapacity      () const { return maxcap; }
			inline size_t	scatterCap       () const { return nscatters * sizeof(S_REF); }
			inline S_REF*   refsMem          () { return d_refs_mem; }
			inline S_REF*   scatter          () { return d_scatter; }
			inline uint32*  cnfMem           () { return d_cnf_mem; }
			inline uint32*  unitsdPtr        () { return d_units; }
			inline CNF*     pinnedCNF        () { return pinned_cnf; }
			inline cuLits&  literals         () { return litsPool; }
			void			reset            (const int64 _free) {
				this->_free = _free - penalty;
			}
			void			init             (const int64 _tot, const int64 _penalty) {
				penalty = _penalty;
				this->_tot = _tot;
				_free = _tot - penalty;
			}
			inline void		cacheCNFPtr      (const CNF* d_cnf, const cudaStream_t& _s = 0) {
				CHECK(cudaMemcpyAsync(pinned_cnf, d_cnf, sizeof(CNF), cudaMemcpyDeviceToHost, _s));
			}
			inline void		prefetchCNF2D    (const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching CNF to global memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(cnfPool.mem, cnfPool.cap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
			}
			inline void		prefetchCNF2H    (const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching CNF to system memory..");
					CHECK(cudaMemPrefetchAsync(cnfPool.mem, cnfPool.cap, cudaCpuDeviceId, _s));
					PFLDONE(2, 5);
				}
			}
			inline void		prefetchOT2H     (const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching OT to system memory..");
					CHECK(cudaMemPrefetchAsync(otPool.mem, otPool.cap, cudaCpuDeviceId, _s));
					PFLDONE(2, 5);
				}
			}
			void	        resizeCNFAsync   (CNF*, const S_REF&, const uint32&);
			bool			resizeCNF        (CNF*&, const size_t&, const size_t&);
			uint32*         resizeLits       (const size_t&);
			addr_t			allocTemp        (const size_t&);
			bool			allocAux         (const size_t&);
			bool			allocHist        (cuHist&);
			bool			allocVars        (VARS*&, const size_t&);
			bool			resizeOTAsync    (OT*&, const size_t&, const cudaStream_t& _s = 0);
			void			createMirror     (CNF*&, const size_t&, const size_t&);
			void			mirrorCNF        (CNF*&);
			void			scatterCNF       (CNF*, S_REF*, Byte*);
			void			compactCNF       (CNF*, CNF*);
			
		};
		/*****************************************************/
		/*  Usage:    Thrust cached memory allocator         */
		/*  Dependency: None                                 */
		/*****************************************************/
		class INVALID_PTR {
			string msg;
		public:
			INVALID_PTR(void* p) : msg() {
				msg = string(CERROR) + "pointer \"" + std::to_string((size_t)p) + "\" not allocated by Thrust cached allocator" + CNORMAL;
			}
			virtual arg_t what() const { return msg.c_str(); }
		};
		class TCA {
			typedef std::map<char*, int64>      allocBlock_t;
			typedef std::multimap<int64, char*> freeBlock_t;
			freeBlock_t     freeBlocks;
			allocBlock_t    allocBlocks;
			size_t			used, maxcap;

		public:
			typedef char value_type;
			TCA() { used = maxcap = 0; }
			void	destroy      ();
			char*   allocate     (int64);
			void	deallocate   (char*, size_t);
			size_t	maxCapacity  () const { return maxcap; }
			void	updateMaxCap () { if (maxcap < used) maxcap = used; }
		};

	}
}

#endif
