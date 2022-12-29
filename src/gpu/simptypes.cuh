/***********************************************************************[simptypes.cuh]
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

#include <thrust/device_ptr.h>
#include "vector.cuh"
#include "definitions.cuh"
#include "constants.cuh"
#include "definitions.hpp"
#include "sclause.hpp"
#include "clause.hpp"
#include "sort.hpp"

namespace ParaFROST {

	class cuTIMER {
	private:
		cudaEvent_t _start, _stop;
		float _gpuTime;
	public:
		float vo, ve, sub, bce, ere, cot, rot, sot, sig, gc, io;
		cuTIMER() {
			RESETSTRUCT(this);
			cudaEventCreate(&_start);
			cudaEventCreate(&_stop);
		}
		~cuTIMER() {
			cudaEventDestroy(_start);
			cudaEventDestroy(_stop);
			_start = NULL, _stop = NULL;
		}
		inline void  start  (const cudaStream_t& _s = 0) { cudaEventRecord(_start, _s); }
		inline void  stop   (const cudaStream_t& _s = 0) { cudaEventRecord(_stop, _s); }
		inline float gpuTime() {
			_gpuTime = 0;
			cudaEventSynchronize(_stop);
			cudaEventElapsedTime(&_gpuTime, _start, _stop);
			return _gpuTime;
		}
	};
	extern cuTIMER* cutimer;
	/*****************************************************/
	/*  Usage:    Global simplifier types                */
	/*  Dependency: none                                 */
	/*****************************************************/
	typedef thrust::device_ptr<uint32> t_iptr;
	typedef cuVec<uint32> cuVecU;
	typedef cuVec<Byte> cuVecB;
	typedef cuVec<S_REF> cuREF;
	typedef cuREF OL;
	typedef Vec<S_REF> HOL;
	typedef Vec<HOL> HOT;
	/*****************************************************/
	/*  Usage:    raw data memory/pointer types          */
	/*  Dependency: none                                 */
	/*****************************************************/
	struct cuPool {
		addr_t mem;
		size_t cap;
	};
	struct cuCNF {
		uint32* mem;
		S_REF size, cap;
	};
	struct cuLits {
		uint32* mem;
		t_iptr  thrust_lits;
		size_t  size, cap;
	};
	/*****************************************************/
	/*  Usage:  Containers for LCVE, hist, stats, limits */
	/*  Dependency: none                                 */
	/*****************************************************/
	struct VARS {
		cuVecU* pVars, * units, * resolved;
		uint32* pVarsData, * pVarsSize, * scores, * varcore, * eligible, * cachedUnits;
		Byte*   eliminated, *cachedEliminated;
		cuVecU  tmpUnits;
		uint32  numPVs, nUnits, nMelted;
		bool	isEliminatedCached;
		VARS() : 
			pVars(NULL)
			, units(NULL)
			, resolved(NULL)
			, pVarsData(NULL)
			, pVarsSize(NULL)
			, scores(NULL)
			, varcore(NULL)
			, eligible(NULL)
			, cachedUnits(NULL)
			, eliminated(NULL)
			, cachedEliminated(NULL)
			, numPVs(0), nUnits(0), nMelted(0)
			, isEliminatedCached(false)
		{}
	};
	struct cuHist {
		S_REF   * d_segs;
		uint32  * d_hist, * h_hist;
		uint32	* d_vorg;
		Byte	* d_lbyte;		// only used if proof is enabled
		t_iptr	thrust_hist;
		cuHist		() : 
			d_segs(NULL)
			, d_hist(NULL)
			, h_hist(NULL)
			, d_vorg(NULL) 
			, d_lbyte(NULL)
		{}
		~cuHist		() 
		{
			h_hist = NULL, d_hist = NULL;
			d_segs = NULL, d_vorg = NULL;
			d_lbyte = NULL; 
		}
		inline uint32	operator[]	(const uint32& i) const { assert(h_hist && i < inf.nDualVars); return h_hist[i]; }
		inline uint32&	operator[]	(const uint32& i) { assert(h_hist && i < inf.nDualVars); return h_hist[i]; }
		inline void		cacheHist	(const cudaStream_t& _s = 0) {
			CHECK(cudaMemcpyAsync(h_hist, d_hist, inf.nDualVars * sizeof(uint32), cudaMemcpyDeviceToHost, _s));
		}
		inline void		fetchVars	(const uint32* vorg, const cudaStream_t& _s = 0) {
			CHECK(cudaMemcpyAsync(d_vorg, vorg, (inf.maxVar + 1) * sizeof(uint32), cudaMemcpyHostToDevice, _s));
		}
	};
	struct cuOptions {
		// [0] = sub_limit, [1] = bce_limit,
		// [2] = ere_limit, [3] = xor_max_arity
		// [4] = ve_clause_limit, [5] = opts.ve_lbound_en;
		uint32 options[NLIMITS];
		cuOptions() {
			options[0] = 3000, options[1] = 3000;
			options[2] = 3000, options[3] = 10;
			options[4] = 100;
			options[5] = options[6] = 0;
		}
	};
	/*****************************************************/
	/*  Usage:    Simplifier CNF                         */
	/*  Dependency: none                                 */
	/*****************************************************/
	class CNF {
		cuCNF _data;
		cuREF _refs;
		Byte _bucket;
	public:
		_PFROST_H_D_					CNF			() :
			_data({ NULL, 0, 0 })
			, _bucket((Byte)sizeof(uint32))
			 { }
		_PFROST_H_D_					CNF			(const S_REF& data_cap, const uint32& cs_cap) :
			_bucket((Byte)sizeof(uint32)) 
		{
			assert(_bucket == sizeof(uint32));
			assert(data_cap);
			assert(cs_cap);
			S_REF* csmem = (S_REF*)(this + 1);
			_refs.alloc(csmem, cs_cap);
			_data.mem = (uint32*)(csmem + cs_cap);
			_data.cap = data_cap, _data.size = 0;
		}
		_PFROST_H_D_		void		resize		(const S_REF& data_size, const uint32& cs_size) {
			assert(data_size);
			assert(cs_size);
			assert(data_size <= _data.cap);
			assert(cs_size <= _refs.capacity());
			_data.size = data_size, _refs.resize(cs_size);
		}
		_PFROST_H_D_		void		resize		(const uint32& cs_size) {
			assert(cs_size);
			assert(cs_size <= _refs.capacity());
			_refs.resize(cs_size);
		}
		_PFROST_H_D_		void		fixPointer	() {
			assert(_data.size);
			S_REF* csmem = (S_REF*)(this + 1);
			_refs.alloc(csmem);
			_data.mem = (uint32*)(csmem + _refs.size());
		}
		_PFROST_H_D_ 		cuCNF&		data		() { return _data; }
		_PFROST_H_D_ 		cuREF&		refs		() { return _refs; }
		_PFROST_H_D_		size_t		bucket		() const { return _bucket; }
		_PFROST_H_D_		uint32		size		() const { return _refs.size(); }
		_PFROST_H_D_		uint32		empty		() const { return !_refs.size(); }
		_PFROST_H_D_ 		S_REF*		refsData    (const uint32& i = 0) { return _refs + i; }
		_PFROST_H_D_		S_REF		ref			(const uint32& i) { assert(i < _refs.size()); return _refs[i]; }
		_PFROST_H_D_ const	S_REF&		ref			(const uint32& i) const { assert(i < _refs.size()); return _refs[i]; }
		_PFROST_H_D_		SCLAUSE&	clause		(const uint32& i) { assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
		_PFROST_H_D_ const	SCLAUSE&	clause		(const uint32& i) const { assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
		_PFROST_H_D_		SCLAUSE*	cref		(const S_REF& r) { return (SCLAUSE*)(_data.mem + r); }
		_PFROST_H_D_ const	SCLAUSE*	cref		(const S_REF& r) const { return (SCLAUSE*)(_data.mem + r); }
		_PFROST_H_D_		SCLAUSE&	operator[]	(const S_REF& r) { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
		_PFROST_H_D_ const	SCLAUSE&	operator[]	(const S_REF& r) const { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
		_PFROST_H_D_		size_t		calcSize	(const int& nLits) { return sizeof(SCLAUSE) + nLits * sizeof(uint32); }
		_PFROST_H_D_		void		print		(const uint32& off = 0, const bool& p_ref = true) {
			for (uint32 i = off; i < size(); i++) {
				SCLAUSE& c = clause(i);
				if (c.size()) {
					if (p_ref) printf("c  C(%d, r: %lld)->", i, uint64(_refs[i]));
					else printf("c  C(%d)->", i);
					c.print();
				}
			}
		}
		_PFROST_H_D_		void		printAdded	(const uint32& off = 0, const bool& p_ref = true) {
			for (uint32 i = off; i < size(); i++) {
				SCLAUSE& c = clause(i);
				if (c.size() && c.added()) {
					if (p_ref) printf("c  C(%d, r: %lld)->", i, uint64(_refs[i]));
					else printf("c  C(%d)->", i);
					c.print();
				}
			}
		}
		_PFROST_H_D_		void		printDeleted(const uint32& off = 0, const bool& p_ref = true) {
			for (uint32 i = off; i < size(); i++) {
				SCLAUSE& c = clause(i);
				if (c.size() && c.deleted()) {
					if (p_ref) printf("c  C(%d, r: %lld)->", i, uint64(_refs[i]));
					else printf("c  C(%d)->", i);
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
			_refs._push(_data.size);
			_data.size += S_REF(cBytes / _bucket);
		}
		_PFROST_H_			void		newClause	(SCLAUSE& src) {
			assert(_data.size < _data.cap);
			assert(src.size());
			new (cref(_data.size)) SCLAUSE(src);
			assert(src.size() == cref(_data.size)->size());
			_refs._push(_data.size);
			_data.size += src.blockSize();

		}
		_PFROST_H_			void		copyFrom	(CNF* src) {
			assert(_data.mem != NULL);
			assert(_data.cap);
			assert(_data.size == 0);
			assert(_refs.empty());
			PFLOGN2(2, " Compacting simplified CNF (%d to %d) on CPU..", src->size(), inf.nClauses);
			for (uint32 i = 0; i < src->size(); i++) {
				SCLAUSE& s = src->clause(i);
				if (s.learnt() || s.original())
					newClause(s);
			}
			PFLDONE(2, 5);
		}
		_PFROST_IN_D_			S_REF*		jump	(S_REF&, const uint32&, const uint32&);
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
		~OT() { _lists = NULL, _occurs = NULL; }
		_PFROST_H_D_			OT				() : _lists(NULL), _occurs(NULL), maxLists(0) {}
		_PFROST_H_D_			OT				(const uint32& nlists) : maxLists(nlists) {
			assert(nlists);
			_lists = (OL*)(this + 1);
			_occurs = (S_REF*)(_lists + maxLists);
		}
		_PFROST_H_D_	S_REF*	data			(const uint32& i = 0) { return _occurs + i; }
		_PFROST_H_D_	OL&		operator []		(const uint32& i) { assert(i < maxLists); return _lists[i]; }
		_PFROST_H_D_	OL		operator []		(const uint32& i) const { assert(i < maxLists); return _lists[i]; }
		_PFROST_H_D_			operator OL*	() { return _lists; }
		_PFROST_H_D_	uint32	size			() const { return maxLists; }
		_PFROST_H_D_	void	print			() {
			for (uint32 v = 2; v < size(); v++) {
				int sign_v = SIGN(v) ? -int(ABS(v)) : ABS(v);
				printf("c  list[ %d ][cap = %d]", sign_v, _lists[v].capacity()), _lists[v].print();
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
	/*  Usage:    GPU shared memory template             */
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
	/*  Usage:    GPU comparators                        */
	/*  Dependency: CNF, dtypes                          */
	/*****************************************************/
	template<class T>
	struct GPU_DEFAULT_CMP {
		_PFROST_D_ bool operator () (const T& x, const T& y) {
			return x < y;
		}
	};
	struct GPU_LCV_CMP {
		const uint32* scores;
		_PFROST_H_D_ GPU_LCV_CMP(const uint32* _scores) : scores(_scores) { assert(scores != NULL); }
		_PFROST_D_ bool operator () (const uint32& a, const uint32& b) const {
			const uint32 x = scores[a], y = scores[b];
			if (x < y) return true;
			if (x > y) return false;
			return a < b;
		}
	};
	struct GPU_MCV_CMP {
		const uint32* scores;
		_PFROST_H_D_ GPU_MCV_CMP(const uint32* _scores) : scores(_scores) { assert(scores != NULL); }
		_PFROST_D_ bool operator () (const uint32& a, const uint32& b) const {
			const uint32 x = scores[a], y = scores[b];
			if (x > y) return true;
			if (x < y) return false;
			return a > b;
		}
	};
	struct COMPACT_VARS {
		_PFROST_D_ bool operator()(const uint32& var) {
			return var;
		}
	};
	struct COMPACT_CMP {
		const CNF* cnf;
		_PFROST_H_D_ COMPACT_CMP() : cnf(NULL) {}
		_PFROST_H_D_ COMPACT_CMP(const CNF* _cnf) : cnf(_cnf) { assert(cnf != NULL); }
		_PFROST_H_D_ void init(const CNF* _cnf) { cnf = _cnf; }
		_PFROST_D_ bool operator()(const S_REF& ref) {
			return !cnf->cref(ref)->deleted();
		}
	};
	struct OLIST_CMP {
		const CNF* cnf;
		_PFROST_H_D_ OLIST_CMP() : cnf(NULL) {}
		_PFROST_H_D_ OLIST_CMP(const CNF* _cnf) : cnf(_cnf) { assert(cnf != NULL); }
		_PFROST_H_D_ void init(const CNF* _cnf) { cnf = _cnf; }
		_PFROST_D_ bool operator () (const S_REF& a, const S_REF& b) {
			const SCLAUSE& x = (*cnf)[a];
			const SCLAUSE& y = (*cnf)[b];
			uint32 xdata = uint32(x.size()), ydata = uint32(y.size());
			if (NEQUAL(xdata, ydata)) return xdata < ydata;
			xdata = x[0], ydata = y[0];
			if (NEQUAL(xdata, ydata)) return xdata < ydata;
			xdata = x.back(), ydata = y.back();
			if (NEQUAL(xdata, ydata)) return xdata < ydata;
			xdata = x.sig(), ydata = y.sig();
			if (NEQUAL(xdata, ydata)) return xdata < ydata;
			return a < b;
		}
	};

}

#endif
