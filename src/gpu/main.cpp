/***********************************************************************[main.cpp]
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

#include "control.hpp"
#include "solve.hpp"
#include "version.hpp"

using namespace ParaFROST;

bool quiet_en = false;
int  verbose  = -1;

int main(int argc, char **argv)
{             
	BOOL_OPT opt_quiet_en("quiet", "enable quiet mode, same as verbose=0", false);
	INT_OPT opt_verbose("verbose", "set the verbosity", 1, INT32R(0, 4));
	INT_OPT opt_timeout("timeout", "set out-of-time limit in seconds", 0, INT32R(0, INT32_MAX));
	INT_OPT opt_memoryout("memoryout", "set out-of-memory limit in gigabytes", 0, INT32R(0, 256));
	OPTION_VEC& options = ARG::opts();
	try {
		bool parsed = parseArguments(argc, argv);
		quiet_en = opt_quiet_en, verbose = opt_verbose;
		if (quiet_en) verbose = 0;
		else if (!verbose) quiet_en = true;
		if (!quiet_en && verbose) {
			PFNAME("ParaFROST (Parallel Formal ReasOning about Satisfiability)", version());
			PFAUTHORS("Muhammad Osama Mahmoud");
			PFLRULER('-', RULELEN);
			if (parsed) {
				PFLOGN0(" Embedded options: ");
				for (int i = 0, j = 0; i < options.size(); i++) {
					if (options[i]->isParsed()) {
						options[i]->printArgument();
						if (++j % 4 == 0) { PUTCH('\n'); PFLOGN0("\t\t      "); }
					}
				}
				PUTCH('\n'); PFLRULER('-', RULELEN);
			}
		}
		signal_handler(handler_terminate);
		string formula = argc > 1 ? argv[1] : "";
		Solver* parafrost = new Solver(formula);
		solver = parafrost;
		if (opt_timeout > 0) set_timeout(opt_timeout);
		if (opt_memoryout > 0) set_memoryout(opt_memoryout);
		signal_handler(handler_mercy_interrupt, handler_mercy_timeout);
		parafrost->solve();
		if (!quiet_en) PFLOG0("");
		PFLOGN2(1, " Cleaning up..");
		solver = NULL;
		delete parafrost;
		PFLDONE(1, 5);
		if (!quiet_en) PFLRULER('-', RULELEN);
		return EXIT_SUCCESS;
	}
	catch (std::bad_alloc&) {
		PRINT("%s%s%s", CYELLOW, "ARE YOU SERIOUS NOW?\n", CNORMAL);
		PFLOGS("UNKNOWN");
		return EXIT_FAILURE;
	}
	catch (MEMOUTEXCEPTION&) {
		syncAll();
		PRINT("%s%s%s", CYELLOW, "SOLVER MEMORY OUT\n", CNORMAL);
		PFLOGS("UNKNOWN");
		return EXIT_FAILURE;
	}
	catch (CACHEMEMOUT&) {
		syncAll();
		PRINT("%s%s%s", CYELLOW, "CACHER MEMORY OUT\n", CNORMAL);
		PFLOGS("UNKNOWN");
		return EXIT_FAILURE;
	}
}