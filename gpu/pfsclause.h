/***********************************************************************[pfsclause.h]
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

#include "pfcudefs.h"
#include "pfconst.h"

namespace pFROST {

	namespace SIGmA {
		/*****************************************************/
		/*  Usage:      abstract clause for GPU simplifier   */
		/*  Dependency:  none                                */
		/*****************************************************/
		typedef uint32 S_REF;
		class SCLAUSE {
			CL_ST _st, _f;
			int _sz;
			uint32 _sig;
			union { uint32 _lits[1]; };
		public:
			_PFROST_H_D_ S_REF		blockSize			() const { assert(_sz); return (_sz - 1) + sizeof(S_REF); }
			_PFROST_H_D_ size_t		capacity			() const { assert(_sz); return size_t(_sz - 1) * sizeof(uint32) + sizeof(*this); }
			_PFROST_H_D_			SCLAUSE				() { _sz = 0, _st = ORIGINAL, _f = CFREEZE; }
			_PFROST_H_D_			SCLAUSE				(const int& size) { _sz = size, _st = ORIGINAL, _f = CFREEZE; }
			_PFROST_H_D_			SCLAUSE				(uint32* lits, const int& size) {
				_sig = 0, _st = ORIGINAL, _f = CFREEZE;
				_sz = size, copyLitsFrom(lits);
			}
			_PFROST_H_D_			SCLAUSE				(uint32* lits, const int& size, const CL_ST& state) {
				_sig = 0, _f = CFREEZE;
				_st = state, _sz = size, copyLitsFrom(lits);
			}
			_PFROST_H_D_			SCLAUSE				(SCLAUSE& src) {
				assert(!src.deleted());
				_f = src.molten();
				_st = src.status();
				_sig = src.sig();
				_sz = src.size();
				copyLitsFrom(src);
			}
			_PFROST_H_D_ void		copyLitsFrom		(uint32* src) {
				assert(_sz);
				for (int k = 0; k < _sz; k++) _lits[k] = src[k];
			}
			_PFROST_H_D_ void		resize				(const int& size) { _sz = size; }
			_PFROST_H_D_ void		push				(const uint32& lit) { _lits[_sz++] = lit; }
			_PFROST_H_D_ void		set_sig				(const uint32& sig) { _sig = sig; }
			_PFROST_H_D_ void		set_ref				(const C_REF& _ref) { _sig = _ref; }
			_PFROST_H_D_ void		set_status			(const CL_ST& status) { _st = status; }
			_PFROST_H_D_ uint32&	operator[]			(const int& i)		 { assert(i < _sz); return _lits[i]; }
			_PFROST_H_D_ uint32		operator[]			(const int& i) const { assert(i < _sz); return _lits[i]; }
			_PFROST_H_D_ uint32*	data				(const int& i = 0) { assert(i < _sz); return _lits + i; }
			_PFROST_H_D_ uint32*	end					() { return _lits + _sz; }
			_PFROST_H_D_ uint32		back				() { assert(_sz); return _lits[_sz - 1]; }
			_PFROST_H_D_ operator	uint32*				() { assert(_sz); return _lits; }
			_PFROST_H_D_ void		pop					() { _sz--; }
			_PFROST_H_D_ void		clear				() { _sz = 0; }
			_PFROST_H_D_ void		markDeleted			() { _st = DELETED; }
			_PFROST_H_D_ void		melt				() { _f = CMELT; }
			_PFROST_H_D_ void		freeze				() { _f = CFREEZE; }
			_PFROST_H_D_ int		size				() const { return _sz; }
			_PFROST_H_D_ uint32		sig					() const { return _sig; }
			_PFROST_H_D_ uint32		ref					() const { return _sig; }
			_PFROST_H_D_ bool		deleted				() const { return _st & DELETED; }
			_PFROST_H_D_ bool		learnt				() const { return _st & LEARNT; }
			_PFROST_H_D_ bool		original			() const { return _st & ORIGINAL; }
			_PFROST_H_D_ CL_ST		status				() const { return _st; }
			_PFROST_H_D_ bool		molten				() const { return _f; }
			_PFROST_H_D_ void		filter				() {
				_sig = 0;
				int newSz = 1;
				for (int k = 1; k < _sz; k++) {
					uint32 next = _lits[k];
					if (_lits[k - 1] != next) {
						_lits[newSz++] = next;
						_sig |= MAPHASH(next);
					}
				}
				_sz = newSz;
			}
			_PFROST_H_D_ void		filterCopy			(uint32* src, const int& size) {
				_sig = 0, _sz = 1;
				*_lits = *src;
				for (int k = 1; k < size; k++) {
					uint32 next = src[k];
					if (src[k - 1] != next) {
						_lits[_sz++] = next;
						_sig |= MAPHASH(next);
					}
				}
				assert(isSorted());
				assert(hasZero() < 0);
			}
			_PFROST_H_D_ void		shareTo				(uint32* dest) {
				assert(_sz > 1);
				uint32* s = _lits, * e = s + _sz;
				while (s != e) *dest++ = *s++;
			}
			_PFROST_H_D_ bool		has					(const uint32& lit) const { // binary search
				assert(_sz);
				if (_sz == 2) {
					if (_lits[0] == lit || _lits[1] == lit) return true;
					else return false;
				}
				else {
					assert(isSorted());
					int low = 0, high = _sz - 1, mid;
					uint32 first = _lits[low], last = _lits[high];
					while (first <= lit && last >= lit) {
						mid = ((low + high) >> 1);
						uint32 m = _lits[mid];
						if (m < lit) first = _lits[low = mid + 1];
						else if (m > lit) last = _lits[high = mid - 1];
						else return true; // found
					}
					if (_lits[low] == lit) return true; // found
					else return false; // Not found
				}
			}
			_PFROST_H_D_ bool		isSorted			() const {
				for (int i = 0; i < _sz; i++)
					if (i > 0 && _lits[i] < _lits[i - 1])
						return false;
				return true;
			}
			_PFROST_H_D_ int		hasZero				() const {
				for (int i = 0; i < _sz; i++)
					if (_lits[i] == 0) 
						return i;
				return -1;
			}
			_PFROST_H_D_ void		print				() {
				printf("(");
				for (int l = 0; l < _sz; l++) {
					int lit = int(ABS(_lits[l]));
					lit = (ISNEG(_lits[l])) ? -lit : lit;
					printf("%4d ", lit);
				}
				char st = 'U';
				if (deleted()) st = 'X';
				else if (original()) st = 'O';
				else if (learnt()) st = 'L';
				printf(") %c, s=0x%X\n", st, _sig);
			}
		};
	}
}

#endif