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
#include "pfvector.h"
#include <cassert>
#include <string>

namespace pFROST {

	/*****************************************************/
	/*  Usage:   structure for CNF clause storage        */
	/*  Dependency:  none                                */
	/*****************************************************/
	class CLAUSE {
		bool	_l : 1;
		bool	_d : 1;
		bool	_k : 1;
		bool	_b : 1;
		bool	_r : 1;
		bool	_m : 1;
		bool	_s : 1;
		bool	_h : 1;
		bool	_v;
		CL_ST	_used; 
		int		_sz, _pos, _lbd;
		union { uint32 _lits[2]; C_REF _ref; };
	public:
		size_t capacity() const { return (size_t(_sz) - 2) * sizeof(uint32) + sizeof(*this); }
		inline CLAUSE() 
			: _sz(0)
			, _pos(2)
			, _lbd(0)
			, _used(0)
			, _l(false)
			, _d(false)
			, _r(false)
			, _m(false)
			, _s(false)
			, _h(false)
			, _k(true)
			, _v(false)
			, _b(false)
		{ }
		inline CLAUSE(const int& size) 
			: _sz(size)
			, _pos(2)
			, _lbd(0)
			, _used(0)
			, _l(false)
			, _d(false)
			, _r(false)
			, _m(false)
			, _s(false)
			, _h(false)
			, _v(false)
			, _k(true)
		{ assert(_sz > 1); _b = _sz == 2; }
		inline CLAUSE(const Lits_t& lits) 
			: _sz(lits.size())
			, _pos(2)
			, _lbd(0)
			, _used(0)
			, _l(false)
			, _d(false)
			, _r(false)
			, _m(false)
			, _s(false)
			, _h(false)
			, _v(false)
			, _k(true)
		{ assert(_sz > 1); _b = _sz == 2; copyLitsFrom(lits); }
		inline CLAUSE(const CLAUSE& src) 
			: _sz(src.size())
			, _pos(src.pos())
			, _lbd(src.lbd())
			, _used(src.usage())
			, _l(src.learnt())
			, _d(src.deleted())
			, _r(src.reason())
			, _s(src.subsume())
			, _v(src.vivify())
			, _h(src.hyper())
			, _k(src.keep())
			, _m(false)
		{ assert(_sz > 1); _b = _sz == 2; copyLitsFrom(src); }
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
		inline	operator	uint32*		() { assert(_sz); return _lits; }
		inline	uint32*		data		() { return _lits; }
		inline	uint32*		mid			() { assert(_pos < _sz); return _lits + _pos; }
		inline	uint32*		end			() { return _lits + _sz; }
		inline	void		pop			() { assert(_sz); _sz--, _b = _sz == 2; }
		inline	C_REF		ref			() const { assert(_m); return _ref; }
		inline	int			pos			() const { assert(_pos > 1); return _pos; }
		inline	CL_ST		usage		() const { assert(_used <= USAGET2); return _used; }
		inline	int			size		() const { return _sz; }
		inline	int			lbd			() const { return _lbd; }
		inline	bool		original	() const { return !_l; }
		inline	bool		learnt		() const { return _l; }
		inline	bool		deleted		() const { return _d; }
		inline	bool		reason		() const { return _r; }
		inline	bool		moved		() const { return _m; }
		inline	bool		subsume		() const { return _s; }
		inline	bool		vivify		() const { return _v; }
		inline	bool		keep		() const { return _k; }
		inline	bool		hyper		() const { return _h; }
		inline	bool		binary		() const { assert(_sz > 2 || _b); return _b; }
		inline	void		initTier1	() { _used = USAGET2; }
		inline	void		initTier2	() { _used = USAGET2; }
		inline	void		initTier3	() { _used = USAGET3; }
		inline	void		warm		() { assert(_used); _used--; }
		inline  void		initMoved	() { _m = false; }
		inline	void		initReason	() { _r = false; }
		inline	void		initSubsume	() { _s = false; }
		inline	void		initVivify	() { _v = false; }
		inline	void		markOriginal() { _l = false; }
		inline	void		markLearnt	() { _l = true; }
		inline	void		markReason	() { _r = true; }
		inline	void		markMoved	() { _m = true; }
		inline	void		markHyper	() { _h = true; }
		inline	void		markSubsume	() { _s = true; }
		inline	void		markVivify	() { _v = true; }
		inline	void		markDeleted	() { _d = true; }
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
		inline	void		set_usage	(const CL_ST& usage) { _used = usage; }
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
			fprintf(stdout, ") %c:%d[u:%d k:%d h:%d lbd:%d]\n", 
				st, reason(), usage(), keep(), hyper(), lbd());
		}
	};
	/*****************************************************/
	/*  Usage: memory manager for CNF clauses            */
	/*  Dependency:  CLAUSE, SMM                         */
	/*****************************************************/
	typedef SMM<Byte, C_REF> CTYPE;
	class CMM : public CTYPE
	{
		Vec<bool, C_REF> stencil;
	public:
								CMM				() { assert(CTYPE::bucket() == 1); }
		explicit				CMM				(const C_REF& init_cap) : CTYPE(init_cap), stencil(init_cap, 0) { assert(CTYPE::bucket() == 1); }
		inline void				init			(const C_REF& init_cap) { CTYPE::init(init_cap * sizeof(CLAUSE)), stencil.resize(init_cap, 0); }
		inline		 CLAUSE&	operator[]		(const C_REF& r) { return (CLAUSE&)CTYPE::operator[](r); }
		inline const CLAUSE&	operator[]		(const C_REF& r) const { return (CLAUSE&)CTYPE::operator[](r); }
		inline		 CLAUSE*	clause			(const C_REF& r) { return (CLAUSE*)address(r); }
		inline const CLAUSE*	clause			(const C_REF& r) const { return (CLAUSE*)address(r); }
		inline bool				deleted			(const C_REF& r) const { assert(check(r)); return stencil[r]; }
		inline void				collectClause	(const C_REF& r, const int& size) { CTYPE::collect(bytes(size)); assert(check(r)); stencil[r] = true; }
		inline void				collectLiterals	(const int& size) { CTYPE::collect(size * sizeof(uint32)); }
		inline void				migrateTo		(CMM& dest) { 
			CTYPE::migrateTo(dest);
			stencil.migrateTo(dest.stencil);
		}
		template <class SRC>
		inline C_REF			alloc			(const SRC& src) {
			assert(src.size() > 1);
			size_t cBytes = bytes(src.size());
			C_REF r = CTYPE::alloc(cBytes);
			new (clause(r)) CLAUSE(src);
			assert(clause(r)->capacity() == cBytes);
			assert(src.size() == clause(r)->size());
			stencil.expand(r + 1, 0);
			return r;
		}
		inline C_REF			alloc			(const int& size) {
			assert(size > 1);
			size_t cBytes = bytes(size);
			C_REF r = CTYPE::alloc(cBytes);
			new (clause(r)) CLAUSE(size);
			assert(clause(r)->capacity() == cBytes);
			assert(size == clause(r)->size());
			stencil.expand(r + 1, 0);
			return r;
		}
		inline size_t			bytes			(const int& size) {
			assert(size > 1); 
			return (sizeof(CLAUSE) + (size_t(size) - 2) * sizeof(uint32));
		}
		inline void				destroy			() { dealloc(), stencil.clear(true); }
	};

}

#endif
