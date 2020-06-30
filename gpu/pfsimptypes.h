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

#include <thrust/device_vector.h>
#include "pfsclause.h"
#include "pfclause.h"
#include "pfdefs.h"
#include "pfcuvec.h"
#include "pfcuconst.h"
#include "pfsort.h"

namespace pFROST {
	extern cudaDeviceProp	devProp;
	extern uint32			maxGPUThreads;
	namespace SIGmA {
		/*****************************************************/
		/*  Usage:    Global simplifier types                */
		/*  Dependency: none                                 */
		/*****************************************************/
		typedef thrust::device_vector<uint32> uTHRVector;
		typedef cuVec<uint32> cuVecU;
		typedef cuVecU OL;
		/*****************************************************/
		/*  Usage:    raw data memory types for CNF / OT     */
		/*  Dependency: none                                 */
		/*****************************************************/
		struct cuPool { addr_t mem; size_t cap; };
		struct cuCNF { S_REF* mem, size, cap; };
		/*****************************************************/
		/*  Usage:    Containers for LCVE, stats             */
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
			uint32 numPVs, nUnits, mu_inc;
			VARS() : gstats(NULL)
					, pVars(NULL)
					, units(NULL)
					, scores(NULL)
					, resolved(NULL)
					, eligible(NULL)
					, cachedUnits(NULL)
					, numPVs(0), nUnits(0), mu_inc(0) {}
			~VARS() { 
				if (cachedUnits != NULL) delete[] cachedUnits;
				cachedUnits = NULL;
			}
		};
		/*****************************************************/
		/*  Usage:    Simplifier CNF                         */
		/*  Dependency: none                                 */
		/*****************************************************/
		typedef cuVec<S_REF> cuREF;
		class CNF {
			cuCNF _data;
			cuREF cs;
			Byte _bucket;
		public:
			_PFROST_H_D_					CNF			() : _bucket((Byte)sizeof(S_REF)), _data({ NULL, 0, 0 }) { assert(_bucket == sizeof(uint32)); }
			_PFROST_H_D_					CNF			(const S_REF& data_cap, const uint32& cs_cap) : _bucket((Byte)sizeof(S_REF)) {
				assert(_bucket == sizeof(uint32));
				assert(data_cap);
				assert(cs_cap);
				_data.cap = data_cap, _data.size = 0;
				_data.mem = (S_REF*)(this + 1);
				cs.alloc(_data.mem + data_cap, cs_cap);
			}
			_PFROST_H_D_		void		resize		(const S_REF& data_size, const uint32& cs_size) {
				assert(data_size);
				assert(cs_size);
				_data.size = data_size, cs.resize(cs_size);
			}
			_PFROST_H_D_		void		resize		(const uint32& cs_size) { 
				assert(cs_size);
				cs.resize(cs_size); 
			}
			_PFROST_H_D_		void		fixPointer	() { assert(_data.cap); _data.mem = (S_REF*)(this + 1), cs.alloc(_data.mem + _data.cap); }
			_PFROST_H_D_ 		cuCNF&		data		() { return _data; }
			_PFROST_H_D_ 		S_REF*		csData		() { return cs; }
			_PFROST_H_D_		size_t		bucket		() const { return _bucket; }
			_PFROST_H_D_		uint32		size		() const { return cs.size(); }
			_PFROST_H_D_		uint32		empty		() const { return cs.size() == 0; }
			_PFROST_H_D_		S_REF		ref			(const uint32& i)		{ assert(i < cs.size()); return cs[i]; }
			_PFROST_H_D_ const	S_REF&		ref			(const uint32& i) const { assert(i < cs.size()); return cs[i]; }
			_PFROST_H_D_		SCLAUSE&	clause		(const uint32& i)		{ assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
			_PFROST_H_D_ const	SCLAUSE&	clause		(const uint32& i) const { assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
			_PFROST_H_D_		SCLAUSE*	cref		(const S_REF& r)	   { return (SCLAUSE*)(_data.mem + r); }
			_PFROST_H_D_ const	SCLAUSE*	cref		(const S_REF& r) const { return (SCLAUSE*)(_data.mem + r); }
			_PFROST_H_D_		SCLAUSE&	operator[]	(const S_REF& r)	   { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
			_PFROST_H_D_ const	SCLAUSE&	operator[]	(const S_REF& r) const { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
			_PFROST_H_D_		size_t		calcSize	(const int& nLits) { return (sizeof(SCLAUSE) + (size_t(nLits) - 1) * sizeof(uint32)); }
			_PFROST_H_D_		void		print		(const bool& p_ref = true) {
				for (S_REF i = 0; i < size(); i++) {
					SCLAUSE& c = clause(i);
					if (c.size()) {
						if (p_ref) printf("c | C(%d, r: %d)->", i, cs[i]);
						else printf("c | C(%d)->", i);
						c.print();
					}
				}
			}
			_PFROST_H_			void		newClause	(CLAUSE& src, const bool& _sort) {
				assert(src.size() > 1);
				size_t cBytes = calcSize(src.size());
				assert(_data.size < _data.cap);
				new (cref(_data.size)) SCLAUSE(src, src.size(), src.status());
				if (_sort) Sort(cref(_data.size)->data(), src.size());
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
				uint32 size = src->size();
				for (uint32 i = 0; i < size; i++) {
					SCLAUSE& s = src->clause(i);
					if (s.learnt() || s.original())
						newClause(s);
				}
			}
			_PFROST_D_			void		insert		(SCLAUSE& src);
			_PFROST_D_			S_REF*		jump		(S_REF&, const uint32&, const uint32&);
		};
		/*****************************************************/
		/*  Usage:    Simplifier occurrence table            */
		/*  Dependency: none                                 */
		/*****************************************************/
		class OT {
			OL* lists;
			uint32* occurs;
			uint32 maxLists, maxEntries;
		public:
									~OT			() { lists = NULL, occurs = NULL; }
			_PFROST_H_D_			OT			() : lists(NULL), occurs(NULL), maxLists(0), maxEntries(0) {}
			_PFROST_H_D_			OT			(const uint32& nlists) : maxLists(nlists), maxEntries(0) {
				assert(nlists);
				lists = (OL*)(this + 1);
				occurs = (uint32*)(lists + maxLists);
			}
			_PFROST_D_		uint32* data		(const uint32&);
			_PFROST_H_D_	OL&		operator [] (const uint32& i)		{ assert(i < maxLists); return lists[i]; }
			_PFROST_H_D_	OL		operator [] (const uint32& i) const { assert(i < maxLists); return lists[i]; }
			_PFROST_H_D_			operator OL*() { return lists; }
			_PFROST_H_D_	void	resetCap	() { maxEntries = 0; }
			_PFROST_H_D_	uint32	capacity	() const { return maxEntries; }
			_PFROST_H_D_	uint32	size		() const { return maxLists; }
			_PFROST_H_D_	void	print		() {
				for (uint32 v = 2; v < size(); v++) {
					int64 sign_v = ISNEG(v) ? -int64(ABS(v)) : ABS(v);
					printf("c | list[%lld][cap = %d]", sign_v, lists[v].capacity()), lists[v].print();
				}
			}
			_PFROST_H_		bool	accViolation() {
				for (uint32 v = 0; v < inf.maxVar; v++) {
					uint32 p = V2D(v + 1), n = NEG(p);
					if (lists[p].size() > lists[p].capacity()) {
						PFLOGEN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
							v + 1, lists[v].capacity(), lists[v].size());
						return false;
					}
					if (lists[n].size() > lists[n].capacity()) {
						PFLOGEN("list(%d) size exceeded allocated capacity (cap: %d, sz: %d):",
							v + 1, lists[v].capacity(), lists[v].size());
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
		/*  Dependency: none                                 */
		/*****************************************************/
		class cuMM {
			cuPool hcnfPool, cnfPool, otPool; // dynamic pools
			cuPool varsPool; // fixed pool
			size_t _tot, _free, cap, maxcap;
			S_REF  *d_cs_mem, *d_cnf_mem;
			uint32 *d_units, otBlocks;
			inline bool	hasFreeMem		(const char* name) {
				PFLOG2(2, " Allocating GPU memory for %s (used/free = %.2f/%zd MB)", name, double(cap) / MBYTE, _free / MBYTE);
				if (cap > _free) { PFLOGW("not enough memory (current = %zd MB) -> skip SIGmA", cap / MBYTE); return false; }
				_free -= cap;
				assert(_tot >= _free);
				assert(_tot > cap);
				return true;
			}

		public:
						cuMM			() : 
										cnfPool({ NULL, 0 })
										, hcnfPool({ NULL, 0 })
										, varsPool({ NULL, 0 })
										, otPool({ NULL, 0 })
										, d_units(NULL)
										, d_cs_mem(NULL)
										, d_cnf_mem(NULL)
										, cap(0)
										, _tot(0)
										, _free(0)
										, maxcap(0)
										, otBlocks(0) {}
				void	destroy			();
				void	breakMirror		();
		inline 	void	updateMaxCap	() { if (maxcap < cap) maxcap = cap; }
		inline	bool	empty			() const { return cap == 0; }
		inline	size_t	capacity		() const { return cap; }
		inline	size_t	maxCapacity		() const { return maxcap; }
		inline	S_REF*	cnfClsdPtr		() { return d_cs_mem; }
		inline	S_REF*	cnfDatadPtr		() { return d_cnf_mem; }
		inline	uint32*	unitsdPtr		() { return d_units; }
		inline	void	prefetchCNF		(const cudaStream_t& _s = (cudaStream_t)0) {
					   if (devProp.major > 5) {
						   PFLOGN2(2, " Prefetching CNF to global to memory..");
						   CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
						   CHECK(cudaMemPrefetchAsync(cnfPool.mem, cnfPool.cap, MASTER_GPU, _s));
						   PFLDONE(2, 5);
					   }
				   }
				bool	allocFixed		(VARS*&, const uint32&);
				void	mirrorCNF		(CNF*&);
				void	resizeHostCNF	(CNF*&, const uint32&, const uint32&);
				bool	resizeCNF		(CNF*&, const uint32&, const uint32&);
				void	resizeCNFAsync	(CNF*, CNF*);
				bool	resizeOTAsync	(OT*&, uint32*, const uint32&, const cudaStream_t& _s = (cudaStream_t)0);
				void	resetOTCapAsync	(OT*, const cudaStream_t& _s = (cudaStream_t)0);
				void	init			(const size_t _free) { this->_free = _tot = _free; }
		};

	}
}

#endif
