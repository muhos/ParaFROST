/***********************************************************************[restart.h]
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

#ifndef __RESTART_
#define __RESTART_

#include "definitions.hpp"

namespace ParaFROST {

    // reference: https://doi.org/10.29007/89dw
    struct EMA {
        double val, a, b;
        int64 c, p;

        inline EMA() : val(0), a(0), b(0), c(0), p(0) {}
        inline EMA(const int& window) : val(0), b(1.0), c(0), p(0) {
            assert(window >= 1);
            a = 1.0 / (double)window;
            assert(a >= 0), assert(b >= a), assert(b <= 1);
        }
        inline operator double() const { return val; }
        inline void update(const double& x) {
            val += b * (x - val);
            if (b <= a || c--) return;
            c = p = ((p + 1) << 1) - 1;
            b *= 0.5;
            if (b < a) b = a;
        }
    };
    struct LBDFS { EMA fast, slow; };

    class LBDREST {
        LBDFS average[2];
        double rate;
        int lbdfast, lbdslow;
    public:
        LBDREST() : rate(0), lbdfast(0), lbdslow(0) {}
        inline void init(const double& r, const int& f, const int& s) {
            rate = r, lbdfast = f, lbdslow = s;
        }
        inline void reset() {
            assert(rate && lbdfast && lbdslow);
            average[0].fast = EMA(lbdfast);
            average[0].slow = EMA(lbdslow);
            average[1].fast = EMA(lbdfast);
            average[1].slow = EMA(lbdslow);
        }
        inline void update(const bool& stable, const double& x) {
            average[stable].fast.update(x);
            average[stable].slow.update(x);
        }
        inline bool restart(const bool& stable) const { 
            return (rate * average[stable].slow) <= average[stable].fast;
        }
    };

    // reference: MiniSat, SAT4j, CadiCal
    class LUBYREST {
        uint64 u, v, limit;
        uint64 period, countdown;
        bool restart, limited;

    public:
        LUBYREST() : period(0), restart(false) {}
        void enable(uint64 period, uint64 limit) {
            if (limit && period > limit) period = limit;
            this->period = period;
            this->limit = limit;
            countdown = period;
            limited = (limit > 0);
            restart = false;
            u = v = 1;
            
        }
        void disable() { period = 0, restart = false; }
        void update();
        operator bool() {
            if (!restart) return false;
            restart = false;
            return true;
        }
    };

}

#endif