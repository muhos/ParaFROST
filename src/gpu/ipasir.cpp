/***********************************************************************[ipasir.cpp]
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

#include "solve.hpp"
#include "version.hpp"
#include "ipasir.hpp"

using namespace ParaFROST;

class ipasir_t: public Solver {
    Lits_t filtered, clause;
    bool nomodel;

	inline uint32 abs(int lit) {
        return uint32(lit < 0 ? -lit : lit);
	}
	inline uint32 import(int lit) {
        uint32 v = abs(lit);
        while (v > inf.maxVar) iadd();
        return V2DEC(v, (lit < 0));
	}

    public:

    ipasir_t() : clause(INIT_CAP), nomodel(false) {}
    ~ipasir_t() {}

    void add(int lit) {
        nomodel = true;
        if (lit)
            clause.push(import(lit));
        else {
            itoClause(filtered, clause);
            clause.clear();
        }
    }
    void assume(int32_t lit) {
        nomodel = true;
        assumptions.push(import(lit));
    }
    int solve() {
        isolve(assumptions);
        assumptions.clear();
        nomodel = (cnfstate != SAT);
        return UNSOLVED(cnfstate) ? 0 : (cnfstate == SAT ? 10 : 20);
    }
    int val(int lit) {
        if (nomodel) return 0;
        return model.satisfied(import(lit)) ? lit : -lit;
    }
    int failed(int lit) {
        return ifailed(import(lit));
    }
};

extern "C" {

#include "ipasir.hpp"

#define IPASIR(S) ((ipasir_t*)S)

const char* ipasir_signature() { return signature(); }
void* ipasir_init() { return new ipasir_t(); }
void ipasir_release(void* s) { delete IPASIR(s); }
int ipasir_solve(void* s) { return IPASIR(s)->solve(); }
void ipasir_add(void* s, int l) { IPASIR(s)->add(l); }
void ipasir_assume(void* s, int l) { IPASIR(s)->assume(l); }
int ipasir_val(void* s, int l) { return IPASIR(s)->val(l); }
int ipasir_failed(void* s, int l) { return IPASIR(s)->failed(l); }
void ipasir_set_terminate(void* s, void* state, int (*callback)(void* state)) { IPASIR(s)->setTermCallback(state, callback); }
void ipasir_set_learn(void* s, void* state, int max_length, void (*learn)(void* state, int* clause)) { IPASIR(s)->setLearnCallback(state, max_length, learn); }

};