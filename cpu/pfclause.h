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
		CL_ST	_st, _k, _b, _r, _u, _m, _s;
		int		_sz, _pos, _lbd;
		union { uint32 _lits[2]; C_REF _ref; };
	public:
		size_t				capacity	() const { return (size_t(_sz) - 2) * sizeof(uint32) + sizeof(*this); }
		inline				CLAUSE		() {
			_sz = 0, _b = 0, _pos = 2, _lbd = 0;
			_st = 0, _r = 0, _u = 0, _m = 0, _s = 0, _k = 1;
		}
		inline				CLAUSE		(const int& size) {
			assert(size > 1);
			_sz = size, _b = size == 2, _pos = 2, _lbd = 0;
			_st = 0, _r = 0, _u = 0, _m = 0, _s = 0, _k = 1;
		}
		inline				CLAUSE		(const Lits_t& lits) {
			assert(lits.size() > 1);
			_sz = lits.size(), copyLitsFrom(lits);
			_b = _sz == 2, _pos = 2, _lbd = 0;
			_st = 0, _r = 0, _u = 0, _m = 0, _s = 0, _k = 1;
		}
		inline				CLAUSE		(const CLAUSE& src) {
			assert(src.size() > 1);
			_sz = src.size(), copyLitsFrom(src);
			_b = _sz == 2, _pos = src.pos(), _st = src.status(), _r = src.reason(), _m = 0;
			_lbd = src.lbd(), _u = src.usage(), _s = src.subsume(), _k = src.keep();
		}
		template <class SRC>
		inline	void		copyLitsFrom(const SRC& src) {
			assert(_sz > 1);
			for (int k = 0; k < _sz; k++) {
				assert(src[k] > 1);
				_lits[k] = src[k];
			}
		}
		inline	uint32&		operator[]	(const int& i) { assert(i < _sz); return _lits[i]; }
		inline	uint32		operator[]	(const int& i) const { assert(i < _sz); return _lits[i]; }
		inline	operator	uint32*		() { assert(_sz != 0); return _lits; }
		inline	uint32*		data		() { return _lits; }
		inline	uint32*		mid			() { assert(_pos < _sz); return _lits + _pos; }
		inline	uint32*		end			() { return _lits + _sz; }
		inline	void		pop			() { _sz--, _b = _sz == 2; }
		inline	int			size		() const { return _sz; }
		inline	int			pos			() const { return _pos; }
		inline	C_REF		ref			() const { return _ref; }
		inline	int			lbd			() const { return _lbd; }
		inline	CL_ST		status		() const { return _st; }
		inline	bool		deleted		() const { return _st & DELETED; }
		inline	bool		original	() const { return _st & ORIGINAL; }
		inline	bool		learnt		() const { return _st & LEARNT; }
		inline	bool		binary		() const { return _b; }
		inline	bool		reason		() const { return _r; }
		inline	bool		moved		() const { return _m; }
		inline	bool		subsume		() const { return _s; }
		inline	bool		keep		() const { return _k; }
		inline	CL_ST		usage		() const { return _u; }
		inline	void		initTier1	() { _u = USAGET2; }
		inline	void		initTier2	() { _u = USAGET2; }
		inline	void		initTier3	() { _u = USAGET3; }
		inline	void		warm		() { _u--; }
		inline  void		initMoved	() { _m = 0; }
		inline	void		initReason	() { _r = 0; }
		inline	void		initSubsume	() { _s = 0; }
		inline	void		markReason	() { _r = 1; }
		inline	void		markMoved	() { _m = 1; }
		inline	void		markSubsume	() { _s = 1; }
		inline	void		markDeleted	() { _st = DELETED; }
		inline	void		shrink		(const int& n) { 
			_sz -= n; 
			if (_pos >= _sz) _pos = 2;
			_b = _sz == 2;
		}
		inline	void		resize		(const int& newSz) {
			_sz = newSz;
			if (_pos >= _sz) _pos = 2;
			_b = _sz == 2;
		}
		inline	void		set_ref		(const C_REF& r) { _m = 1, _ref = r; }
		inline	void		set_pos		(const int& newPos) { assert(newPos >= 2); _pos = newPos; }
		inline	void		set_lbd		(const int& lbd) { _lbd = lbd; }
		inline	void		set_status	(const CL_ST& status) { _st = status; }
		inline	void		set_usage	(const CL_ST& usage) { _u = usage; }
		inline	void		set_keep	(const bool& keep) { _k = keep; }
		inline	void		print		() const {
			fprintf(stdout, "(");
			for (int l = 0; l < _sz; l++) {
				int lit = int(ABS(_lits[l]));
				lit = (SIGN(_lits[l])) ? -lit : lit;
				fprintf(stdout, "%4d ", lit);
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
	typedef SMM<Byte, C_REF> CTYPE;
	class CMM : public CTYPE
	{
	public:
								CMM			() {}
		explicit				CMM			(const C_REF& init_cap) : CTYPE(init_cap) {}
		inline void				init		(const C_REF& init_cap) { CTYPE::init(init_cap * sizeof(CLAUSE)); }
		inline		 CLAUSE&	operator[]	(const C_REF& r) { return (CLAUSE&)CTYPE::operator[](r); }
		inline const CLAUSE&	operator[]	(const C_REF& r) const { return (CLAUSE&)CTYPE::operator[](r); }
		inline		 CLAUSE*	clause		(const C_REF& r) { return (CLAUSE*)address(r); }
		inline const CLAUSE*	clause		(const C_REF& r) const { return (CLAUSE*)address(r); }
		inline void				collect		(const C_REF& r) { CTYPE::collect(calcSize(clause(r)->size())); }
		inline void				collect		(const int& size) { CTYPE::collect(size); }
		inline size_t			calcSize	(const int& size) { return (sizeof(CLAUSE) + (size_t(size) - 2) * sizeof(uint32)); }
		inline void				migrate		(CMM& newblock) { CTYPE::migrate(newblock); }
		template <class SRC>
		inline C_REF			alloc		(const SRC& src) {
			assert(src.size() > 1);
			size_t cBytes = calcSize(src.size());
			C_REF r = CTYPE::alloc(cBytes);
			new (clause(r)) CLAUSE(src);
			assert(clause(r)->capacity() == cBytes);
			assert(src.size() == clause(r)->size());
			return r;
		}
		inline C_REF			alloc		(const int& size) {
			assert(size > 1);
			size_t cBytes = calcSize(size);
			C_REF r = CTYPE::alloc(cBytes);
			new (clause(r)) CLAUSE(size);
			assert(clause(r)->capacity() == cBytes);
			assert(size == clause(r)->size());
			return r;
		}
		inline void				destroy		() { dealloc(); }
	};

}

#endif
