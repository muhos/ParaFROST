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
	cout << "c |--------------------------------------------------------------------------------------|" << endl;
	cout << "c |                              ParaFROST SAT Solver                                    |" << endl;
	cout << "c | Technische Universiteit Eindhoven (TU/e), all rights reserved.                       |" << endl;
	cout << "c |--------------------------------------------------------------------------------------|" << endl;
	try {
		if (argc == 1) { cout << "No input file specified." << endl; exit(EXIT_FAILURE); }
		parseArguments(argc, argv);
		cout << "c | Embedded options: ";
		for (int i = 0, j = 0; i < options.size(); i++) {
			if (options[i]->isParsed()) {
				options[i]->printArgument();
				if (++j % 4 == 0) cout << "\nc |                   ";
			}
		}
		cout << "\nc |--------------------------------------------------------------------------------------|" << endl;
		string path = argv[1];
		ParaFROST* pFrost = new ParaFROST(path);
		g_pFrost = pFrost;
		if (pFrost->timeout > 0) set_timeout(pFrost->timeout);
		sig_handler(handler_exit);
		pFrost->solve();
		delete pFrost;
		cout << "c |--------------------------------------------------------------------------------------|" << endl;
		return EXIT_SUCCESS;
	}
	catch (MEMOUTEXCEPTION&) {
		cout << "c |\ns UNKNOWN" << endl;
		exit(EXIT_SUCCESS);
	}
}