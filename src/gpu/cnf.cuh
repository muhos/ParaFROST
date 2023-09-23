/***********************************************************************[cnf.cuh]
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

#ifndef __GPU_CNF_
#define __GPU_CNF_

#include "vector.cuh"
#include "sclause.cuh"
#include "definitions.cuh"

namespace ParaFROST {

	typedef cuVec<S_REF> cuREF;

	struct cuCNF {
		uint32* mem;
		S_REF size, cap;
	};

	class CNF {

		cuCNF _data;
		cuREF _refs;

	public:
		_PFROST_H_D_					CNF			() : _data({ NULL, 0, 0 }) { }
		_PFROST_H_D_					CNF			(const S_REF& data_cap, const uint32& cs_cap)
		{
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
		_PFROST_H_D_		uint32		size		() const { return _refs.size(); }
		_PFROST_H_D_		uint32		empty		() const { return _refs.empty(); }
		_PFROST_H_D_ 		S_REF*		refsData    (const uint32& i = 0) { return _refs + i; }
		_PFROST_H_D_		S_REF		ref			(const uint32& i)		{ assert(i < _refs.size()); return _refs[i]; }
		_PFROST_H_D_ const	S_REF&		ref			(const uint32& i) const { assert(i < _refs.size()); return _refs[i]; }
		_PFROST_H_D_		SCLAUSE&	clause		(const uint32& i)		{ assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
		_PFROST_H_D_ const	SCLAUSE&	clause		(const uint32& i) const { assert(ref(i) < _data.size); return (SCLAUSE&)_data.mem[ref(i)]; }
		_PFROST_H_D_		SCLAUSE*	cref		(const S_REF& r)	   { return (SCLAUSE*)(_data.mem + r); }
		_PFROST_H_D_ const	SCLAUSE*	cref		(const S_REF& r) const { return (SCLAUSE*)(_data.mem + r); }
		_PFROST_H_D_		SCLAUSE&	operator[]	(const S_REF& r)	   { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
		_PFROST_H_D_ const	SCLAUSE&	operator[]	(const S_REF& r) const { assert(r < _data.size); return (SCLAUSE&)_data.mem[r]; }
		_PFROST_H_			void		newClause	(CLAUSE& src) {
			assert(src.size() > 1);
			assert(_data.size < _data.cap);
			new (cref(_data.size)) SCLAUSE(src);
			assert(src.size() == cref(_data.size)->size());
			_refs._push(_data.size);
			_data.size += SCBUCKETS(src.size());
		}
		_PFROST_H_			void		newClause	(SCLAUSE& src) {
			assert(_data.size < _data.cap);
			assert(src.size());
			new (cref(_data.size)) SCLAUSE(src);
			assert(src.size() == cref(_data.size)->size());
			_refs._push(_data.size);
			_data.size += SCBUCKETS(src.size());

		}
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
		_PFROST_IN_D_		S_REF*		jump		(S_REF&, const uint32&, const uint32&);
	};

}

#endif
