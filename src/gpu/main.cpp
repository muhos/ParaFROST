/***********************************************************************[main.cpp]
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

/*

	██████╗░░█████╗░██████╗░░█████╗░███████╗██████╗░░█████╗░░██████╗████████╗
	██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝
	██████╔╝███████║██████╔╝███████║█████╗░░██████╔╝██║░░██║╚█████╗░░░░██║░░░
	██╔═══╝░██╔══██║██╔══██╗██╔══██║██╔══╝░░██╔══██╗██║░░██║░╚═══██╗░░░██║░░░
	██║░░░░░██║░░██║██║░░██║██║░░██║██║░░░░░██║░░██║╚█████╔╝██████╔╝░░░██║░░░
	╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░░░░╚═╝░░░

*/

#include "control.hpp"
#include "banner.hpp"
#include "solve.hpp"

using namespace ParaFROST;

bool quiet_en = false;
int  verbose  = -1;

int main(int argc, char **argv)
{   
	BOOL_OPT opt_quiet_en("quiet", "enable quiet mode, same as verbose=0", false);
	INT_OPT opt_verbose("verbose", "set the verbosity", 1, INT32R(0, 4));
	INT_OPT opt_timeout("timeout", "set the timeout in seconds", 0, INT32R(0, INT32_MAX));
	INT_OPT opt_memoryout("memoryout", "set the memout in gigabytes", 0, INT32R(0, 256));
	OPTION_VEC& options = ARG::opts();
	if (argc == 1) LOGERROR("no input file specified");
	try {
		parseArguments(argc, argv);
		quiet_en = opt_quiet_en, verbose = opt_verbose;
		if (quiet_en) verbose = 0;
		else if (!verbose) quiet_en = true;
		if (!quiet_en && verbose) {
			LOGHEADER(1, 5, "Banner");
			LOGFANCYBANNER(version());
			if (argc > 2) {
                LOGHEADER(1, 5, "Options");
                int i = 0, j = 0;
                const int MAX_PER_LINE = 4;
				for (i = 0, j = 0; i < options.size(); i++) {
					if (options[i]->isParsed()) {
                        if (j++ % MAX_PER_LINE == 0) PUTCH('c');                        
                        options[i]->printArgument();
                        if (j % MAX_PER_LINE == 0) PUTCH('\n');
					}
				}
                if (j < MAX_PER_LINE) PUTCH('\n');
			}
		}
		signal_handler(handler_terminate);
		string formula = argv[1];
		Solver* parafrost = new Solver(formula);
		solver = parafrost;
		if (opt_timeout > 0) set_timeout(opt_timeout);
		if (opt_memoryout > 0) set_memoryout(opt_memoryout);
		signal_handler(handler_mercy_interrupt, handler_mercy_timeout);
		parafrost->solve();
		if (!quiet_en) LOG0("");
		LOGHEADER(1, 5, "Exit");
		solver = NULL;
		delete parafrost;
        if (verbose > 1) LOGRULER('-', RULELEN);
		return EXIT_SUCCESS;
	}
	catch (std::bad_alloc&) {
		CHECK(cudaDeviceReset());
		PRINT("%s%s%s", CYELLOW, "ARE YOU SERIOUS NOW?\n", CNORMAL);
		LOGSAT("UNKNOWN");
		return EXIT_FAILURE;
	}
	catch (MEMOUTEXCEPTION&) {
		CHECK(cudaDeviceReset());
		PRINT("%s%s%s", CYELLOW, "SOLVER MEMORY OUT\n", CNORMAL);
		LOGSAT("UNKNOWN");
		return EXIT_FAILURE;
	}
}