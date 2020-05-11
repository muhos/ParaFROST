/***********************************************************************[pfmain.cpp]
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

#include "pfsolve.h"
#include "pfargs.h"

int main(int argc, char** argv)
{
	Vec<ARG*>& options = ARG::opts();
	if (argc == 1) { PFLOGE("No input file specified"); exit(EXIT_FAILURE); }
	try {
		parseArguments(argc, argv);
		if (!isQuiet()) {
			PFNAME("ParaFROST");
			PFLOGN(" Embedded options: ");
			for (int i = 0, j = 0; i < options.size(); i++) {
				if (options[i]->isParsed()) {
					options[i]->printArgument();
					if (++j % 4 == 0) { putc('\n', stdout); PFLOGN("\t\t\t"); }
				}
			}
			putc('\n', stdout); PFLOGR('-', RULELEN);
		}
		string formula = argv[1];
		if (formula.find(".cnf") == -1 && formula.find(".dimacs") == -1) {
			cout << "Input file not recognizable." << endl;
			exit(EXIT_FAILURE);
		}
		sig_handler(handler_terminate);
		ParaFROST* pFrost = new ParaFROST(formula);
		gpfrost = pFrost;
		if (gpfrost->timeout > 0) set_timeout(pFrost->timeout);
		sig_handler(handler_mercy_intr, handler_mercy_timeout);
		pFrost->solve();
		gpfrost = NULL;
		delete pFrost;
		exit(EXIT_SUCCESS);
	}
	catch (MEMOUTEXCEPTION&) {
		PFLOG("");
		PFLOG("%45s", "Memoryout");
		PFLOG("");
		PFLOGS("UNKNOWN");
		exit(EXIT_SUCCESS);
	}
}