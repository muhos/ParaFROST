/***********************************************************************[clause.h]
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

#ifndef __SCLAUSE_
#define __SCLAUSE_

#include "clause.h"

namespace pFROST {
	/*****************************************************/
	/*  Usage:    abstract clause for simp. on host      */
	/*  Dependency:  none                                */
	/*****************************************************/
	class SCLAUSE {
		unsigned _st : 2, _f : 1, _a : 1, _u : 2;
		unsigned _lbd : 26;
		uint32 _sig;
		int _sz;
		uint32* _lits;
	public:
		SCLAUSE		() :
			_st(ORIGINAL)
			, _f(0)
			, _a(0)
			, _u(0)
			, _lbd(0)
			, _sig(0)
			, _sz(0)
			, _lits(NULL)
		{}
		~SCLAUSE	() { clear(true); }
		SCLAUSE(const CLAUSE& src) { init(src); }
		SCLAUSE	(const Lits_t& src) { init(src); }
		inline void		init		(const Lits_t& src) {
			_st = ORIGINAL;
			_sz = src.size();
			_lbd = 0;
			_sig = 0;
			_f = 0;
			_a = 0;
			_u = 0;
			_lits = new uint32[_sz];
			copyLitsFrom(src);
		}
		inline void		init		(const CLAUSE& src) {
			_st = src.learnt();
			_sz = src.size();
			_sig = 0;
			_f = 0;
			_a = 0;
			assert(original() == !src.learnt());
			assert(!src.deleted());
			if (learnt()) {
				_lbd = src.lbd() & MAX_LBD_M;
				_u = src.usage();
			}
			else { _lbd = 0, _u = 0; }
			_lits = new uint32[_sz];
			copyLitsFrom(src);
		}
		template <class SRC>
		inline void		copyLitsFrom(const SRC& src) {
			assert(_sz);
			for (int k = 0; k < _sz; k++) {
				assert(src[k] > 1);
				_lits[k] = src[k];
			}
		}
		inline void		set_lbd		(const unsigned& lbd) { assert(_lbd < MAX_LBD); _lbd = lbd; }
		inline void		set_sig		(const uint32& sig) { _sig = sig; }
		inline void		set_usage	(const CL_ST& usage) { _u = usage; }
		inline void		set_status	(const CL_ST& status) { _st = status; }
		inline void		shrink		(const int& n) { _sz -= n; }
		inline void		resize		(const int& n) { _sz = n; }
		inline uint32	lit			(const int& i) { assert(i < _sz); return _lits[i]; }
		inline uint32&	operator [] (const int& i) { assert(i < _sz); return _lits[i]; }
		inline uint32	operator [] (const int& i) const { assert(i < _sz); return _lits[i]; }
		inline operator uint32*		() { assert(_sz != 0); return _lits; }
		inline uint32*	data		() { return _lits; }
		inline uint32*	end			() { return _lits + _sz; }
		inline uint32	back		() { return _lits[_sz - 1]; }
		inline void		pop			() { _sz--; }
		inline void		freeze		() { _f = 0; }
		inline void		melt		() { _f = 1; }
		inline void		markAdded	() { _a = 1; }
		inline void		markDeleted	() { _st = DELETED; }
		inline CL_ST	usage		() const { return _u; }
		inline bool		molten		() const { return _f; }
		inline bool		added		() const { return _a; }
		inline bool		empty		() const { return !_sz; }
		inline bool		original	() const { return !_st; }
		inline bool		deleted		() const { return _st & DELETED; }
		inline bool		learnt		() const { return _st & LEARNT; }
		inline CL_ST	status		() const { return _st; }
		inline int		size		() const { return _sz; }
		inline unsigned	lbd			() const { return _lbd; }
		inline uint32	sig			() { return _sig; }
		inline int		hasZero		() {
			for (int i = 0; i < _sz; i++)
				if (!_lits[i]) 
					return i;
			return -1;
		}
		inline bool		isSorted	() {
			for (int i = 0; i < _sz; i++) {
				if (i > 0 && _lits[i] < _lits[i - 1]) return false;
			}
			return true;
		}
		inline void		calcSig		(const uint32& init_sig = 0) {
			_sig = init_sig;
			for (int i = 0; i < _sz; i++)
				_sig |= MAPHASH(_lits[i]);
		}
		inline bool		has			(const uint32& lit) {
			if (_sz == 2) {
				if (_lits[0] == lit || _lits[1] == lit) return true;
				else return false;
			}
			else {
				assert(this->isSorted());
				int low = 0, high = _sz - 1, mid;
				uint32 first = _lits[low], last = _lits[high];
				while (first <= lit && last >= lit) {
					mid = (low + high) >> 1;
					uint32 m = _lits[mid];
					if (m < lit) first = _lits[low = mid + 1];
					else if (m > lit) last = _lits[high = mid - 1];
					else return true; // found
				}
				if (_lits[low] == lit) return true; // found
				else return false; // Not found
			}
		}
		inline void		clear		(bool _free = false) {
			if (_free && _lits != NULL) { delete[] _lits; _lits = NULL; }
			_sz = 0;
		}
		inline void		print		() const {
			printf("(");
			for (int l = 0; l < _sz; l++) {
				int lit = int(ABS(_lits[l]));
				lit = (SIGN(_lits[l])) ? -lit : lit;
				printf("%4d ", lit);
			}
			char st = 'U';
			if (deleted()) st = 'X';
			else if (added()) st = 'A';
			else if (original()) st = 'O';
			else if (learnt()) st = 'L';
			printf(") %c:%d, used=%d, lbd=%d, s=0x%X\n", st, molten(), usage(), _lbd, _sig);
		}
	};
	typedef SCLAUSE* S_REF;

}

#endif
