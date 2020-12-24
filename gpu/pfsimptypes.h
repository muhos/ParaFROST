/***********************************************************************[pfsimptypes.h]
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

#ifndef __SIMP_TYPES_
#define __SIMP_TYPES_

#include <thrust/system/cuda/vector.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/pair.h>
#include <thrust/device_ptr.h>
#include <sstream>
#include <map>
#include "pfsclause.h"
#include "pfclause.h"
#include "pfdefs.h"
#include "pfcuvec.h"
#include "pfcuconst.h"
#include "pfsort.h"


namespace pFROST {

	class cuTIMER {
	private:
		cudaEvent_t _start, _stop;
		float _gpuTime;
	public:
		float vo, ve, hse, bce, ere, sub, cot, rot, sot, sig, gc, io;
		cuTIMER() {
			memset(this, 0, sizeof(*this));
			cudaEventCreate(&_start);
			cudaEventCreate(&_stop);
		}
		~cuTIMER() {
			cudaEventDestroy(_start);
			cudaEventDestroy(_stop);
			_start = NULL, _stop = NULL;
		}
		inline void start(const cudaStream_t& _s = 0) { cudaEventRecord(_start, _s); }
		inline void stop(const cudaStream_t& _s = 0) { cudaEventRecord(_stop, _s); }
		inline float gpuTime() {
			_gpuTime = 0;
			cudaEventSynchronize(_stop);
			cudaEventElapsedTime(&_gpuTime, _start, _stop);
			return _gpuTime;
		}
	};
	extern cuTIMER*			cutimer;
	extern uint32			maxGPUTereads;
	extern cudaDeviceProp	devProp;

	namespace SIGmA {
		/*****************************************************/
		/*  Usage:    Global simplifier types                */
		/*  Dependency: none                                 */
		/*****************************************************/
		typedef thrust::device_ptr<uint32> t_iptr;
		typedef cuVec<uint32> cuVecU;
		typedef cuVec<S_REF> cuREF;
		typedef cuREF OL;
		typedef Vec<S_REF> HOL;
		typedef Vec<HOL> HOT;
		/*****************************************************/
		/*  Usage:    raw data memory/pointer types          */
		/*  Dependency: none                                 */
		/*****************************************************/
		struct cuPool { addr_t mem; size_t cap; };
		struct cuCNF { uint32* mem; S_REF size, cap; };
		/*****************************************************/
		/*  Usage:  Containers for LCVE, hist, stats, limits */
		/*  Dependency: none                                 */
		/*****************************************************/
		struct __align__(16) GSTATS {
			uint32 numDelVars, numClauses, numLits;
			GSTATS() : numLits(0), numDelVars(0), numClauses(0) {}
		};
		struct VARS {
			GSTATS *gstats;
			cuVecU *pVars, *units, *resolved, tmpObj;
			uint32 *scores, *eligible, *cachedUnits;
			uint32 numPVs, nUnits;
			VARS() { memset(this, 0, sizeof(*this)); }
		};
		struct cuHist {
			S_REF	*d_segs;
			uint32	*h_hist, *d_hist, *d_lits;
			t_iptr	thrust_hist, thrust_lits;
			cuHist		() : h_hist(NULL), d_hist(NULL), d_segs(NULL), d_lits(NULL) {}
			~cuHist		() { h_hist = NULL, d_hist = NULL, d_segs = NULL, d_lits = NULL; }
			inline uint32	operator[]	(const uint32& i) const { assert(h_hist && i < inf.nDualVars); return h_hist[i]; }
			inline uint32&	operator[]	(const uint32& i)		{ assert(h_hist && i < inf.nDualVars); return h_hist[i]; }
			inline void		cacheHist	(const cudaStream_t& _s = 0) {
				CHECK(cudaMemcpyAsync(h_hist, d_hist, inf.nDualVars * sizeof(uint32), cudaMemcpyDeviceToHost, _s));
			}
		};
		struct cuLimit {
			// [0] = hse_limit, [1] = bce_limit,
			// [2] = ere_limit, [3] = xor_max_arity
			uint32 limits[NLIMITS];
		};
		/*****************************************************/
		/*  Usage:    Simplifier CNF                         */
		/*  Dependency: none                                 */
		/*****************************************************/
		class CNF {
			cuCNF _data;
			cuREF cs;
			Byte _bucket;
		public:
			_PFROST_H_D_					CNF			() : _bucket((Byte)sizeof(uint32)), _data({ NULL, 0, 0 }) { }
			_PFROST_H_D_					CNF			(const S_REF& data_cap, const uint32& cs_cap) : _bucket((Byte)sizeof(uint32)) {
				assert(_bucket == sizeof(uint32));
				assert(data_cap);
				assert(cs_cap);
				S_REF* csmem = (S_REF*)(this + 1);
				cs.alloc(csmem, cs_cap);
				_data.mem = (uint32*)(csmem + cs_cap);
				_data.cap = data_cap, _data.size = 0;
			}
			_PFROST_H_D_		void		resize		(const S_REF& data_size, const uint32& cs_size) {
				assert(data_size);
				assert(cs_size);
				assert(data_size <= _data.cap);
				assert(cs_size <= cs.capacity());
				_data.size = data_size, cs.resize(cs_size);
			}
			_PFROST_H_D_		void		resize		(const uint32& cs_size) { 
				assert(cs_size);
				assert(cs_size <= cs.capacity());
				cs.resize(cs_size); 
			}
			_PFROST_H_D_		void		fixPointer	() {
				assert(_data.size);
				S_REF* csmem = (S_REF*)(this + 1);
				cs.alloc(csmem);
				_data.mem = (uint32*)(csmem + cs.size());
			}
			_PFROST_H_D_ 		cuCNF&		data		() { return _data; }
			_PFROST_H_D_ 		cuREF&		csVec		() { return cs; }
			_PFROST_H_D_		size_t		bucket		() const { return _bucket; }
			_PFROST_H_D_		uint32		size		() const { return cs.size(); }
			_PFROST_H_D_		uint32		empty		() const { return !cs.size(); }
			_PFROST_H_D_ 		S_REF*		csData		(const uint32& i = 0)	{ return cs + i; }
			_PFROST_H_D_		S_REF		ref			(const uint32& i)		{ assert(i < cs.size()); return cs[i]; }
			_PFROST_H_D_ const	S_REF&		ref			(const uint32& i) const { assert(i < cs.size()); return cs[i]; }
			_PFROST_H_D_		SCLAUSE&	clause		(const uint32& i)		{ assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
			_PFROST_H_D_ const	SCLAUSE&	clause		(const uint32& i) const { assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
			_PFROST_H_D_		SCLAUSE*	cref		(const S_REF& r)	   { return (SCLAUSE*)(_data.mem + r); }
			_PFROST_H_D_ const	SCLAUSE*	cref		(const S_REF& r) const { return (SCLAUSE*)(_data.mem + r); }
			_PFROST_H_D_		SCLAUSE&	operator[]	(const S_REF& r)	   { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
			_PFROST_H_D_ const	SCLAUSE&	operator[]	(const S_REF& r) const { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
			_PFROST_H_D_		size_t		calcSize	(const int& nLits) { return (sizeof(SCLAUSE) + (size_t(nLits) - 1) * sizeof(uint32)); }
			_PFROST_H_D_		void		print		(const uint32& off = 0, const bool& p_ref = true) {
				for (uint32 i = off; i < size(); i++) {
					SCLAUSE& c = clause(i);
					if (c.size()) {
						if (p_ref) printf("c | C(%d, r: %lld)->", i, uint64(cs[i]));
						else printf("c | C(%d)->", i);
						c.print();
					}
				}
			}
			_PFROST_H_D_		void		printAdded	(const uint32& off = 0, const bool& p_ref = true) {
				for (uint32 i = off; i < size(); i++) {
					SCLAUSE& c = clause(i);
					if (c.size() && c.added()) {
						if (p_ref) printf("c | C(%d, r: %lld)->", i, uint64(cs[i]));
						else printf("c | C(%d)->", i);
						c.print();
					}
				}
			}
			_PFROST_H_D_		void		printDeleted(const uint32& off = 0, const bool& p_ref = true) {
				for (uint32 i = off; i < size(); i++) {
					SCLAUSE& c = clause(i);
					if (c.size() && c.deleted()) {
						if (p_ref) printf("c | C(%d, r: %lld)->", i, uint64(cs[i]));
						else printf("c | C(%d)->", i);
						c.print();
					}
				}
			}
			_PFROST_H_			void		newClause	(CLAUSE& src) {
				assert(src.size() > 1);
				size_t cBytes = calcSize(src.size());
				assert(_data.size < _data.cap);
				new (cref(_data.size)) SCLAUSE(src);
				assert(cref(_data.size)->capacity() == cBytes);
				assert(src.size() == cref(_data.size)->size());
				cs._push(_data.size);
				_data.size += S_REF(cBytes / _bucket);
			}
			_PFROST_H_			void		newClause	(SCLAUSE& src) {
				assert(_data.size < _data.cap);
				assert(src.size());
				new (cref(_data.size)) SCLAUSE(src);
				assert(src.size() == cref(_data.size)->size());
				cs._push(_data.size);
				_data.size += src.blockSize();

			}
			_PFROST_H_			void		copyFrom	(CNF* src) {
				assert(_data.mem != NULL);
				assert(_data.cap);
				assert(_data.size == 0);
				assert(cs.empty());
				PFLOGN2(2, " Compacting simplified CNF (%d to %d) on CPU..", src->size(), inf.nClauses);
				for (uint32 i = 0; i < src->size(); i++) {
					SCLAUSE& s = src->clause(i);
					if (s.learnt() || s.original())
						newClause(s);
				}
				PFLDONE(2, 5);
			}
			_PFROST_D_			S_REF*		jump		(S_REF&, const uint32&, const uint32&);
		};
		/*****************************************************/
		/*  Usage:    Simplifier occurrence table            */
		/*  Dependency: none                                 */
		/*****************************************************/
		class OT {
			OL* _lists;
			S_REF* _occurs;
			uint32 maxLists;
		public:
									~OT				() { _lists = NULL, _occurs = NULL; }
			_PFROST_H_D_			OT				() : _lists(NULL), _occurs(NULL), maxLists(0) {}
			_PFROST_H_D_			OT				(const uint32& nlists) : maxLists(nlists) {
				assert(nlists);
				_lists = (OL*)(this + 1);
				_occurs = (S_REF*)(_lists + maxLists);
			}
			_PFROST_H_D_	S_REF*	data			(const uint32& i = 0)	{ return _occurs + i; }
			_PFROST_H_D_	OL&		operator []		(const uint32& i)		{ assert(i < maxLists); return _lists[i]; }
			_PFROST_H_D_	OL		operator []		(const uint32& i) const { assert(i < maxLists); return _lists[i]; }
			_PFROST_H_D_			operator OL*	() { return _lists; }
			_PFROST_H_D_	uint32	size			() const { return maxLists; }
			_PFROST_H_D_	void	print			() {
				for (uint32 v = 2; v < size(); v++) {
					int sign_v = SIGN(v) ? -int(ABS(v)) : ABS(v);
					printf("c | list[ %d ][cap = %d]", sign_v, _lists[v].capacity()), _lists[v].print();
				}
			}
			_PFROST_H_		bool	accViolation	() {
				for (uint32 v = 1; v <= inf.maxVar; v++) {
					uint32 p = V2L(v), n = NEG(p);
					if (_lists[p].size() > _lists[p].capacity()) {
						PFLOGEN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
							v, _lists[v].capacity(), _lists[v].size());
						return false;
					}
					if (_lists[n].size() > _lists[n].capacity()) {
						PFLOGEN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
							v, _lists[v].capacity(), _lists[v].size());
						return false;
					}
				}
				return true;
			}
		};
		/*****************************************************/
		/*  Usage:    GPU shared memory manager              */
		/*  Dependency: none                                 */
		/*****************************************************/
		template<class T>
		class SharedMemory
		{
		public:
			_PFROST_D_ operator T* () {
				extern __shared__ int _smem[];
				return (T*)_smem;
			}
			_PFROST_D_ operator const T* () const {
				extern __shared__ int _smem[];
				return (T*)_smem;
			}
		};
		/*****************************************************/
		/*  Usage:    GPU global memory manager              */
		/*  Dependency: all above except shared memory       */
		/*****************************************************/
		class cuMM {
			// pools & arrays
			cuPool	hcnfPool, cnfPool, otPool, varsPool;
			cuPool	hhistPool, histPool, auxPool;
			S_REF	*d_cs_mem, *d_scatter, *d_segs;
			uint32	*d_hist, *d_lits, *d_cnf_mem;
			uint32	*d_units, *pinned_units;
			Byte	*d_stencil;
			CNF		*pinned_cnf;
			// trackers
			int64	_tot, _free, cap, dcap, maxcap, penalty;
			size_t	nscatters, litsbytes;
			uint32	otBlocks;
			inline bool		hasUnifiedMem	(const char* name) {
				PFLOG2(2, " Allocating GPU unified memory for %s (used/free = %.2f/%lld MB)", name, double(cap) / MBYTE, _free / MBYTE);
				if (cap >= _free) { if (verbose > 1) PFLOGW("not enough memory for %s (current = %lld MB) -> skip SIGmA", name, cap / MBYTE); return false; }
				_free -= cap;
				assert(_free >= 0);
				assert(_tot >= _free);
				assert(_tot > cap);
				return true;
			}
			inline bool		hasDeviceMem	(const char* name) {
				PFLOG2(2, " Allocating GPU-only memory for %s (used/free = %.2f/%lld MB)", name, double(dcap) / MBYTE, _free / MBYTE);
				if (dcap >= _free) { if (verbose > 1) PFLOGW("not enough memory for %s (current = %lld MB) -> skip SIGmA", name, dcap / MBYTE); return false; }
				_free -= dcap;
				assert(_free >= 0);
				assert(_tot >= _free);
				assert(_tot > dcap);
				return true;
			}

		public:
							cuMM			() { memset(this, 0, sizeof(*this)); }
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
			inline size_t	literalsCap		() const { return litsbytes; }
			inline S_REF*	csMem			() { return d_cs_mem; }
			inline S_REF*	scatter			() { return d_scatter; }
			inline uint32*	cnfMem			() { return d_cnf_mem; }
			inline uint32*	unitsdPtr		() { return d_units; }
			inline CNF*		pinnedCNF		() { return pinned_cnf; }
			inline void		cacheCNFPtr		(CNF* d_cnf, const cudaStream_t& _s = 0) {
				CHECK(cudaMemcpyAsync(pinned_cnf, d_cnf, sizeof(CNF), cudaMemcpyDeviceToHost, _s));
			}
			inline void		prefetchCNF2D	(const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching CNF to global memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(cnfPool.mem, cnfPool.cap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
			}
			inline void		prefetchCNF2H	(const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching CNF to system memory..");
					CHECK(cudaMemPrefetchAsync(cnfPool.mem, cnfPool.cap, cudaCpuDeviceId, _s));
					PFLDONE(2, 5);
				}
			}
			inline void		prefetchOT2H	(const cudaStream_t& _s = 0) {
				if (devProp.major > 5) {
					PFLOGN2(2, " Prefetching OT to system memory..");
					CHECK(cudaMemPrefetchAsync(otPool.mem, otPool.cap, cudaCpuDeviceId, _s));
					PFLDONE(2, 5);
				}
			}
			bool			allocAux		(const size_t&);
			bool			allocVars		(VARS*&, const size_t&);
			bool			allocHist		(cuHist&, const size_t&);
			bool			resizeOTAsync	(OT*&, const size_t&, const cudaStream_t& _s = 0);
			bool			resizeCNF		(CNF*&, const size_t&, const size_t&);
			void			createMirror	(CNF*&, const size_t&, const size_t&);
			void			mirrorCNF		(CNF*&);
			void			compactCNF		(CNF*, CNF*);
			void			resizeCNFAsync	(CNF*, const S_REF&, const uint32&);
			void			init			(const int64 _tot, const int64 _penalty) { 
				penalty = _penalty;
				this->_tot = _tot;
				_free = _tot - penalty;
			}
			void			reset			(const int64 _free) { this->_free = _free - penalty; }
		};
		/*****************************************************/
		/*  Usage:    Thrust cached memory allocator         */
		/*  Dependency: None                                 */
		/*****************************************************/
		class INVALID_PTR {
			string msg;
		public:
							INVALID_PTR	(void* p) : msg() {
				std::stringstream s;
				s << "pointer " << p << " not allocated by mine.";
				msg = s.str();
			}
			virtual arg_t	what		() const { return msg.c_str(); }
		};
		class TCA {
			typedef std::map<char*, int64>      allocBlock_t;
			typedef std::multimap<int64, char*> freeBlock_t;
			freeBlock_t     freeBlocks;
			allocBlock_t    allocBlocks;
			size_t			used, maxcap;

		public:
			typedef char value_type;
					TCA			() { used = maxcap = 0; }
			void	destroy		();
			char*	allocate	(int64);
			void	deallocate	(char*, size_t);
			size_t	maxCapacity	() { return maxcap; }
			void	updateMaxCap() { if (maxcap < used) maxcap = used; }
		};

	}
}

#endif
