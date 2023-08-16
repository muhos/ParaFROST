/***********************************************************************[key.hpp]
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

#ifndef __SORT_KEY_
#define __SORT_KEY_

#include "definitions.hpp"
#include "vector.hpp"

namespace ParaFROST {

class Solver;

//============================//
//  Default Comparators       //
//============================//
template <class T>
struct DEFAULT_RANK {
    T operator()(const T& val) { return val; }
};

struct PTR_RANK {
    size_t operator()(void* ptr) { return (size_t)ptr; }
};

template <class T>
struct LESS {
    bool operator()(const T& x, const T& y) const {
        return x < y;
    }
};

template <class T>
struct GREATER {
    bool operator()(const T& x, const T& y) const {
        return x > y;
    }
};
//============================//
//  Custom Comparators        //
//============================//
struct LCV_CMP {
    const uint32* scores;
    LCV_CMP(const uint32* _scores) : scores(_scores) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const uint32 x = scores[a], y = scores[b];
        if (x < y) return true;
        if (x > y) return false;
        return a < b;
    }
};
struct MCV_CMP {
    const uint32* scores;
    MCV_CMP(const uint32* _scores) : scores(_scores) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const uint32 x = scores[a], y = scores[b];
        if (x > y) return true;
        if (x < y) return false;
        return a > b;
    }
};
struct SCORS_CMP {
    Solver* solver;
    SCORS_CMP(Solver* _solver) : solver(_solver) {}
    inline bool operator()(const uint32& a, const uint32& b) const;
};
struct VSIDS_CMP {
    const Vec<double>& act;
    VSIDS_CMP(const Vec<double>& _act) : act(_act) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const double xact = act[a], yact = act[b];
        if (xact < yact) return true;
        if (xact > yact) return false;
        return a > b;
    }
};
struct QUEUE_CMP {
    const Vec<uint64>& bumped;
    QUEUE_CMP(const Vec<uint64>& _bumped) : bumped(_bumped) {}
    bool operator()(const uint32& a, const uint32& b) const {
        return bumped[a] < bumped[b];
    }
};
struct HIST_LCV_CMP {
    const uVec1D& hist;
    HIST_LCV_CMP(const uVec1D& _hist) : hist(_hist) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const uint32 xh = hist[a], yh = hist[b];
        if (xh < yh) return true;
        if (xh > yh) return false;
        return a < b;
    }
};
struct HIST_MCV_CMP {
    const uVec1D& hist;
    HIST_MCV_CMP(const uVec1D& _hist) : hist(_hist) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const uint32 xh = hist[a], yh = hist[b];
        if (xh > yh) return true;
        if (xh < yh) return false;
        return a < b;
    }
};
struct KEY_CMP_ACTIVITY {
    const Vec<double>& acts;
    const uint32* scores;
    KEY_CMP_ACTIVITY(Vec<double>& _acts, const uint32* _scores) : acts(_acts), scores(_scores) {}
    bool operator()(const uint32& a, const uint32& b) const {
        const double dx = acts[a], dy = acts[b];
        if (dx > dy) return true;
        if (dx < dy) return false;
        const uint32 x = scores[a], y = scores[b];
        if (x > y) return true;
        if (x < y) return false;
        return a > b;
    }
};
struct KEY_CMP_BUMP {
    const Vec<uint64>& bumped;
    KEY_CMP_BUMP(const Vec<uint64>& _bumped) : bumped(_bumped) {}
    bool operator()(const uint32& x, const uint32& y) const {
        return bumped[x] > bumped[y];
    }
};
//============================//
//  Custom Rankers	          //
//============================//
struct QUEUE_RANK {
    const Vec<uint64>& bumped;
    QUEUE_RANK(const Vec<uint64>& _bumped) : bumped(_bumped) {}
    uint64 operator()(const uint32& a) const {
        return bumped[a];
    }
};
struct KEY_RANK_BUMP {
    const Vec<uint64>& bumped;
    KEY_RANK_BUMP(const Vec<uint64>& _bumped) : bumped(_bumped) {}
    uint64 operator()(const uint32& x) const {
        uint64 b = bumped[x];
        return ~b;
    }
};
struct LCV_RANK {
    const uint32* scores;
    LCV_RANK(const uint32* _scores) : scores(_scores) {}
    uint32 operator()(const uint32& a) const {
        return scores[a];
    }
};
struct MCV_RANK {
    const uint32* scores;
    MCV_RANK(const uint32* _scores) : scores(_scores) {}
    uint32 operator()(const uint32& a) const {
        return ~scores[a];
    }
};

} // namespace ParaFROST

#endif