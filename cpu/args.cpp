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

#include "Sort.h"
#include "args.h"

Vec<ARG*> options; // container for all options available

void ARG::insert(ARG* opt) {
    options.push(this);
}

bool isQuiet(void) 
{
    for (int i = 0; i < options.size(); i++) {
        if (options[i]->type == "<bool>" && options[i]->arg == "q" && options[i]->isParsed())
            return true;
    }
    return false;
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
                    fprintf(stderr, "ERROR - Unknown flag \"%s\". Use '-h or --help' for help.\n", argv[i]), exit(EXIT_FAILURE);
                else
                    argv[j++] = argv[i];
            }
        }
    }
    argc -= (i - j);
}

void printUsage(int argc, char** argv, bool verbose)
{
    cout << "c | Usage: parafrost [<option> ...][<infile>.<cnf>][<option> ...]" << endl;
    Sort(options, ARG::ARG_CMP());
    const char* prev_type = NULL;
    fprintf(stderr, "c |\nc | OPTIONS:\n");
    for (int i = 0; i < options.size(); i++) {
        const char* type = options[i]->type;
        if (type != prev_type) fprintf(stderr, "c |\n");
        options[i]->help(verbose);
        prev_type = options[i]->type;
    }
    fprintf(stderr, "c |\nc |  -h or --help  Print help message.\n");
    fprintf(stderr, "c |  --help-more   Print verbose help message.\n");
    fprintf(stderr, "c |\nc |--------------------------------------------------------------------------------------|\n");
    exit(EXIT_SUCCESS);
}