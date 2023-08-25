/***********************************************************************[ipasir.hpp]
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

#ifndef __IPASIR_
#define __IPASIR_

#ifdef __cplusplus
extern "C" {
#endif


// Return solver name + version.
const char* ipasir_signature();

// Create a new solver pointer and return it.
// The solver pointer is used for the all
// IPASIR Functions.
void* ipasir_init();

// Release the solver memory (destructor). The input
// pointer (solver) cannot be used after this call.
void ipasir_release(void* solver);

// Add a literal to the currently added clause or 
// finalize a clause with 0.
void ipasir_add(void* solver, int lit_or_zero);

// Add an assumption for the next ipasir_solve call.
// Call it k times to add k assumptions.
void ipasir_assume(void* solver, int lit);

// Solve the formula under the given assumptions.  
// If the formula is satisfiable the function returns 10
// and the state of the solver is changed to (SAT).
// If the formula is unsatisfiable the function returns 20 
// and the state of the solver is changed to (UNSAT).
// If the search is interrupted (see ipasir_set_terminate) 
// 0 is returned and the state of the solver becomes (INPUT).
int ipasir_solve(void* solver);

// Get the truth value of the given literal.
// Returns lit if true, -lit if false, 0 if
// lit and -lit are both satisfying.
int ipasir_val(void* solver, int lit);


// Check if the given assumption literal (lit) is part
// of the conflicting clause that proved the formula
// to be UNSATISFIABLE. Return 1 if so, 0 otherwise.
int ipasir_failed(void* solver, int lit);

// Set a callback function used to indicate a termination 
// signal to the solver.
void ipasir_set_terminate(void* solver, void* data, int (*terminate)(void* data));

// Set a callback function used to extract learned 
// clauses up to a given length from the solver.
void ipasir_set_learn(void* solver, void* data, int max_length, void (*learn)(void* data, int* clause));


#ifdef __cplusplus
} 
#endif


#endif
