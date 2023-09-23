/***********************************************************************[args.cpp]
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

#include "sort.hpp"
#include "input.hpp"
#include "banner.hpp"
#include "control.hpp"

namespace ParaFROST {

    size_t  RADIXBUFFER[RADIXWIDTH]; // defined here since we don't have 'pfsort.cpp'

    void printUsage(int argc, char** argv, bool verbose)
    {
        LOGHEADER(0, 5, "Banner");
        LOGFANCYBANNER(version());
        LOGHEADER(0, 5, "Build");
        uint64 sysmem = 0;
        getCPUInfo(sysmem);
        getBuildInfo();
        size_t _gfree = 0, _gpenalty = 0;
        getGPUInfo(_gfree, _gpenalty);
        LOGHEADER(0, 5, "Usage");
        LOG0("");
        LOG1(" %sparafrost%s [<formula>.<cnf>][<option> ...]", CSOLVER, CNORMAL);
        LOG0("");
        OPTION_VEC& options = ARG::opts();
        Sort(options.data(), options.size(), ARG::ARG_CMP());
        arg_t prev_type = NULL;
        LOG0("");
        LOG1(" Options (simplification + solve):");
        for (int i = 0; i < options.size(); i++) {
            if (options[i]->type != prev_type) fprintf(stdout, "c \n");
            options[i]->help(verbose);
            prev_type = options[i]->type;
        }
        LOG0("");
        LOG1("  %s-h or --help  print available options.%s", CHELP, CNORMAL);
        LOG1("  %s--helpmore    print available options with verbose message.%s", CHELP, CNORMAL);
        LOG0("");
        LOGRULER('-', RULELEN);
        exit(EXIT_SUCCESS);
    }

    bool parseArguments(int& argc, char** argv)
    {
        if (argc <= 1) return false;
        const char* arg = argv[1];
        int dashes = (arg[0] == '-') + (arg[1] == '-');
        if ((dashes & 1) && arg[1] == 'h')
            printUsage(argc, argv);
        else if ((dashes & 2) && hasstr(arg, "help")) {
            if (hasstr(arg, "more"))
                printUsage(argc, argv, true);
            else
                printUsage(argc, argv);
        }
        struct stat st;
        bool ispath = canAccess(arg, st);
        bool ret = false;
        OPTION_VEC& options = ARG::opts();
        for (int i = 1 + ispath; i < argc; i++) {
            const size_t arglen = strlen(argv[i]);
            if (arglen == 1)
                LOGERROR("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
            else if (arglen > 1) {
                const char* arg = argv[i];
                int dashes = (arg[0] == '-') + (arg[1] == '-');
                if (!dashes)
                    LOGERROR("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                else if ((dashes & 1) && arg[1] == 'h')
                    printUsage(argc, argv);
                if ((dashes & 2) && hasstr(arg, "help")) {
                    if (hasstr(arg, "more"))
                        printUsage(argc, argv, true);
                    else
                        printUsage(argc, argv);
                }
                else {
                    int k = 0;
                    bool parsed = false;
                    while (k < options.size() && !(parsed = options[k++]->parse(argv[i])));
                    if (!parsed)  LOGERROR("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                    if (!ret && parsed) ret = true;
                }
            }
        }
        return ret;
    }

}