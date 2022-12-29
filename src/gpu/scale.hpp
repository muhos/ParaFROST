/***********************************************************************[scale.hpp]
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

#ifndef __SCALE_
#define __SCALE_

#include <cmath>
#include <cassert>
#include "datatypes.hpp"

namespace ParaFROST {

	inline uint64 linear	    (const uint64& n) { return n; }
	inline uint64 quadratic	    (const uint64& n) { return n * n; }
	inline uint64 logn		    (const uint64& n) { return (uint64)log10(n + 10); }
	inline uint64 nlogn		    (const uint64& n) { return n * logn(n); }
	inline uint64 nbylogn	    (const uint64& n) { return n / logn(n); }
	inline uint64 lognlogn	    (const uint64& n) {
		double val = log10(n + 10);
		return uint64(val * val);
	}
	inline uint64 lognlognlogn  (const uint64& n) {
		double val = log10(n + 10);
		return uint64(val * val * val);
	}
	inline uint64 nlognlogn	    (const uint64& n) { return n * lognlogn(n); }
	inline uint64 nlognlognlogn (const uint64& n) { return n * lognlognlogn(n); }
	inline uint64 relscale	    (const uint64& n, const uint64& increase)
	{
		const uint64 factor = logn(n);
		assert(factor >= 1);
		const uint64 dfactor = factor * factor;
		assert(dfactor >= 1);
		uint64 scaled = dfactor * increase;
		assert(increase <= scaled);
		return scaled;
	}

}

#endif 