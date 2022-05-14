/***********************************************************************[args.cpp]
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

#include "sort.h"
#include "input.h"
#include "version.h"

namespace pFROST {

    void printUsage(int argc, char** argv, bool verbose)
    {
        PFNAME("ParaFROST (Parallel Formal Reasoning On Satisfiability)", version());
        PFAUTHORS("Muhammad Osama Mahmoud");
        PFRIGHTS("Technische Universiteit Eindhoven (TU/e)");
        PFLRULER('-', RULELEN);
        PFLOG0("");
        PFLOG1(" %sUsage: parafrost [<option> ...][<infile>.<cnf>][<option> ...]%s", CLGREEN, CNORMAL);
        PFLOG0("");
        OPTION_VEC& options = ARG::opts();
        Sort(options.data(), options.size(), ARG::ARG_CMP());
        arg_t prev_type = NULL;
        PFLOG0("");
        PFLOG1(" %sOptions (simplification + solve):%s", CLBLUE, CNORMAL);
        for (int i = 0; i < options.size(); i++) {
            if (options[i]->type != prev_type) fprintf(stdout, "c \n");
            options[i]->help(verbose);
            prev_type = options[i]->type;
        }
        PFLOG0("");
        PFLOG1("  %s-h or --help  print available options.%s", CLBLUE, CNORMAL);
        PFLOG1("  %s--helpmore   print available options with verbose message.%s", CLBLUE, CNORMAL);
        PFLOG0("");
        PFLRULER('-', RULELEN);
        exit(EXIT_SUCCESS);
    }

    void parseArguments(int& argc, char** argv)
    {
        OPTION_VEC& options = ARG::opts();
        int i, j;
        for (i = j = 1; i < argc; i++) {
            if (strlen(argv[i]) == 2 && eqn(argv[i], "-h"))
                printUsage(argc, argv);
            char* arg = argv[i];
            if (eq(arg, "--") && eq(arg, "help")) {
                if (*arg == '\0')
                    printUsage(argc, argv);
                else if (eq(arg, "more"))
                    printUsage(argc, argv, true);
            }
            else {
                int k = 0;
                bool parsed = false;
                while (k < options.size() && !(parsed = options[k++]->parse(argv[i])));
                if (!parsed) {
                    if (eqn(argv[i], "--") || eqn(argv[i], "-"))
                        PFLOGE("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                    else
                        argv[j++] = argv[i];
                }
            }
        }
        argc -= (i - j);
    }

}