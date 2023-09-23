/***********************************************************************[key.cuh]
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

#ifndef __GPU_KEY_
#define __GPU_KEY_

#include "cnf.cuh"
#include "definitions.cuh"

namespace ParaFROST {

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
