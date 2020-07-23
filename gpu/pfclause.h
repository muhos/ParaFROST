/***********************************************************************[pfclause.h]
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

#ifndef __CLAUSE_
#define __CLAUSE_

#include "pfconst.h"
#include "pfalloc.h"
#include "pfvec.h"
#include <cassert>
#include <string>

namespace pFROST {
	/*****************************************************/
	/*  Usage:   simple structure for CNF clause storage */
	/*  Dependency:  none                                */
	/*****************************************************/
	class CLAUSE {
		CL_ST	_st, _r, _f, _m;
		int		_sz, _pos, _lbd;
		union { uint32 _lits[2]; C_REF _ref; };
	public:
		size_t				capacity		() const { return (_sz - 2) * sizeof(uint32) + sizeof(*this); }
		inline				CLAUSE			() {
			_sz = 0, _pos = 2, _lbd = 0;
			_st = 0, _r = NOREASON, _f = CFREEZE, _m = STILL;
		}
		inline				CLAUSE			(const int& size) {
			assert(size > 1);
			_sz = size, _pos = 2, _r = NOREASON, _f = CFREEZE, _m = STILL;
			_st = 0, _lbd = 0;
		}
		inline				CLAUSE			(const Lits_t& lits) {
			assert(lits.size() > 1);
			_sz = lits.size(), copyLitsFrom(lits);
			_pos = 2, _lbd = 0;
			_st = 0, _r = NOREASON, _f = CFREEZE, _m = STILL;
		}
		inline				CLAUSE			(const CLAUSE& src) {
			assert(src.size() > 1);
			_sz = src.size(), copyLitsFrom(src);
			_pos = src.pos(), _st = src.status(), _r = src.reason(), _m = STILL;
			_lbd = src.lbd(), _f = src.usage();
		}
		template <class SRC>
		inline	void		copyLitsFrom	(const SRC& src) {
			assert(_sz > 1);
			for (int k = 0; k < _sz; k++) _lits[k] = src[k];
		}
		inline	uint32&		operator[]		(const int& i) { assert(i < _sz); return _lits[i]; }
		inline	uint32		operator[]		(const int& i) const { assert(i < _sz); return _lits[i]; }
		inline	operator	uint32*			() { assert(_sz != 0); return _lits; }
		inline	uint32*		data			() { return _lits; }
		inline	uint32*		mid				() { assert(_pos < _sz); return _lits + _pos; }
		inline	uint32*		end				() { return _lits + _sz; }
		inline	void		pop				() { _sz--; }
		inline	int			size			() const { return _sz; }
		inline	int			pos				() const { return _pos; }
		inline	C_REF		ref				() const { return _ref; }
		inline	int			lbd				() const { return _lbd; }
		inline	CL_ST		status			() const { return _st; }
		inline	bool		deleted			() const { return _st & DELETED; }
		inline	bool		original		() const { return _st & ORIGINAL; }
		inline	bool		learnt			() const { return _st & LEARNT; }
		inline	bool		reason			() const { return _r; }
		inline	bool		moved			() const { return _m; }
		inline	CL_ST		usage			() const { return _f; }
		inline	void		initTier2		() { _f = USAGET2; }
		inline	void		initTier3		() { _f = USAGET3; }
		inline	void		warm			() { _f--; }
		inline  void		initMoved		() { _m = STILL; }
		inline	void		initReason		() { _r = NOREASON; }
		inline	void		markReason		() { _r = REASON; }
		inline	void		markDeleted		() { _st = DELETED; }
		inline	void		markMoved		() { _m = MOVED; }
		inline	void		shrink			(const int& n) { _sz -= n; if (_pos >= _sz) _pos = 2; }
		inline	void		resize			(const int& newSz) { _sz = newSz; }
		inline	void		set_ref			(const C_REF& r) { _m = MOVED, _ref = r; }
		inline	void		set_pos			(const int& newPos) { assert(newPos >= 2); _pos = newPos; }
		inline	void		set_lbd			(const uint32& lbd) { _lbd = lbd; }
		inline	void		set_status		(const CL_ST& status) { _st = status; }
		inline	void		swap			(const int& p1, const int& p2) { // swap two literals
			std::swap(_lits[p1], _lits[p2]);
		}
		inline	void		swapWatched		() { // swap watched literals
			std::swap(*_lits, _lits[1]);
		}
		inline	void		print			() const {
			fprintf(stdout, "(");
			for (int l = 0; l < _sz; l++) {
				uint32 lit = _lits[l];
				fprintf(stdout, "%4d ", SIGN(lit) ? -int(ABS(lit)) : int(ABS(lit)));
			}
			char st = 'U';
			if (deleted()) st = 'X';
			else if (original()) st = 'O';
			else if (learnt()) st = 'L';
			fprintf(stdout, ") %c:%d, used=%d, lbd=%d\n", st, reason(), usage(), lbd());
		}
	};
	/*****************************************************/
	/*  Usage: memory manager for CNF clauses            */
	/*  Dependency:  CLAUSE                              */
	/*****************************************************/
	class CMM : public SMM<C_REF>
	{
	public:
		CMM() { assert(sizeof(CLAUSE) % bucket() == 0); }
		explicit				CMM			(const C_REF& init_cap) : SMM<C_REF>(init_cap) {}
		inline void				init		(const C_REF& init_cap) { SMM<C_REF>::init(init_cap * C_REF(sizeof(CLAUSE) / bucket())); }
		inline		 CLAUSE&	operator[]	(const C_REF& r) { return (CLAUSE&)SMM<C_REF>::operator[](r); }
		inline const CLAUSE&	operator[]	(const C_REF& r) const { return (CLAUSE&)SMM<C_REF>::operator[](r); }
		inline		 CLAUSE*	ref			(const C_REF& r) { return (CLAUSE*)address(r); }
		inline const CLAUSE*	ref			(const C_REF& r) const { return (CLAUSE*)address(r); }
		inline void				collect		(const C_REF& r) { SMM<C_REF>::collect(C_REF(calcSize(ref(r)->size()) / bucket())); }
		inline void				collect		(const int& size) { SMM<C_REF>::collect(size); }
		inline size_t			calcSize	(const int& size) { return (sizeof(CLAUSE) + (size_t(size) - 2) * sizeof(C_REF)); }
		inline void				migrate		(CMM& newblock) { SMM<C_REF>::migrate(newblock); }
		template <class SRC>
		inline C_REF			alloc		(const SRC& src) {
			assert(src.size() > 1);
			assert(sizeof(C_REF) == bucket());
			size_t cBytes = calcSize(src.size());
			C_REF r = SMM<C_REF>::alloc(C_REF(cBytes / bucket()));
			new (ref(r)) CLAUSE(src);
			assert(ref(r)->capacity() == cBytes);
			assert(src.size() == ref(r)->size());
			return r;
		}
		inline C_REF			alloc		(const int& size) {
			assert(size > 1);
			assert(sizeof(C_REF) == bucket());
			size_t cBytes = calcSize(size);
			C_REF r = SMM<C_REF>::alloc(C_REF(cBytes / bucket()));
			new (ref(r)) CLAUSE(size);
			assert(ref(r)->capacity() == cBytes);
			assert(size == ref(r)->size());
			return r;
		}
		inline void				move		(C_REF& r, CMM& newBlock) {
			CLAUSE& c = operator[](r);
			if (c.moved()) { r = c.ref(); return; }
			r = newBlock.alloc(c);
			c.set_ref(r);
		}
		inline void				destroy		() { dealloc(); }
	};

}

#endif
