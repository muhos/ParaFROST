/***********************************************************************
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
************************************************************************/

#include "solver.h"

int main(int argc, char **argv)
{             
	bool quiet = false;
	if (argc == 1) { cout << "No input file specified." << endl; exit(EXIT_FAILURE); }
	try {
		parseArguments(argc, argv);
		if (!isQuiet()) {
			cout << "c |--------------------------------------------------------------------------------------|" << endl;
			cout << "c |                              ParaFROST SAT Solver                                    |" << endl;
			cout << "c | Technische Universiteit Eindhoven (TU/e), all rights reserved.                       |" << endl;
			cout << "c |--------------------------------------------------------------------------------------|" << endl;
			cout << "c | Embedded options: ";
			for (int i = 0, j = 0; i < options.size(); i++) {
				if (options[i]->isParsed()) {
					options[i]->printArgument();
					if (++j % 4 == 0) cout << "\nc |                   ";
				}
			}
			cout << "\nc |--------------------------------------------------------------------------------------|" << endl;
		}
		string formula = argv[1];
		if (formula.find(".cnf") == -1 && formula.find(".dimacs") == -1) {
			cout << "Input file not recognizable." << endl; 
			exit(EXIT_FAILURE); 
		}
		sig_handler(handler_terminate);
		ParaFROST* pFrost = new ParaFROST(formula);
		g_pFrost = pFrost;
		if (g_pFrost->timeout > 0) set_timeout(pFrost->timeout);
		sig_handler(handler_mercy_intr, handler_mercy_timeout);
		pFrost->solve();
		delete pFrost;
		return EXIT_SUCCESS;
	}
	catch (MEMOUTEXCEPTION&) {
		printf("c |\n");
		printf("c |%45s\n", "Memoryout");
		printf("c |\n");
		printf("s UNKNOWN\n");
		exit(EXIT_FAILURE);
	}
}