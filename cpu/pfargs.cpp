/***********************************************************************[pfargs.cpp]
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

#include "pfsort.h"
#include "pfargs.h"

namespace pFROST {

    Vec<ARG*, int> options; // container for all options available

    void printUsage(int argc, char** argv, bool verbose)
    {
        PFNAME("ParaFROST");
        PFLOG0(" Usage: parafrost [<option> ...][<infile>.<cnf>][<option> ...]");
        Sort(options.data(), options.size(), ARG::ARG_CMP());
        arg_t prev_type = NULL;
        PFLOG0("");
        PFLOG0(" Options (simplification + solve):");
        for (int i = 0; i < options.size(); i++) {
            if (options[i]->type != prev_type) fprintf(stdout, "c |\n");
            options[i]->help(verbose);
            prev_type = options[i]->type;
        }
        PFLOG0("");
        PFLOG0("  -h or --help  print available options.");
        PFLOG0("  --help-more   print available options with verbose message.");
        PFLOG0("");
        PFLOGR('-', RULELEN);
        exit(EXIT_SUCCESS);
    }

    void parseArguments(int& argc, char** argv)
    {
        int i, j;
        for (i = j = 1; i < argc; i++) {
            if (string(argv[i]) == "-h") printUsage(argc, argv);
            char* arg = argv[i];
            if (eq(arg, "--") && eq(arg, "help")) {
                if (*arg == '\0')
                    printUsage(argc, argv);
                else if (eq(arg, "-more"))
                    printUsage(argc, argv, true);
            }
            else {
                int k = 0;
                bool parsed = false;
                while (k < options.size() && !(parsed = options[k++]->parse(argv[i])));
                if (!parsed) {
                    if (eq(argv[i], "--"))
                        PFLOGE("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                    else
                        argv[j++] = argv[i];
                }
            }
        }
        argc -= (i - j);
    }

    void ARG::insert(ARG* opt) { options.push(this); }

}