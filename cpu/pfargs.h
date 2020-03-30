/***********************************************************************[pfargs.h]
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

#ifndef __ARGS_
#define __ARGS_

#include "pfdefs.h"

// global
bool isQuiet(void);
void parseArguments(int& argc, char** argv);
void printUsage(int  argc, char** argv, bool verbose = false);

class ARG
{
protected:
	arg_t arg, text, type;
	bool parsed;
	struct ARG_CMP {
		bool operator()(const ARG* x, const ARG* y) {
			int type_t = strcmp(x->type, y->type);
			return type_t < 0 || type_t == 0 && strcmp(x->arg, y->arg) < 0;
		}
	};
	ARG() { arg = ""; text = ""; type = ""; parsed = false; }
	ARG(arg_t a, arg_t x, arg_t t) {
		arg = a; text = x; type = t; parsed = false;
		insert(this);
	}
public:
	// global
	friend bool isQuiet(void);
	friend void parseArguments(int& argc, char** argv);
	friend void printUsage(int  argc, char** argv, bool verbose);
	// derived methods
	virtual ~ARG() {}
	virtual bool parse(arg_t input) = 0;
	virtual void help(bool verbose = false) = 0;
	virtual void printArgument() = 0;
	// local
	void insert(ARG*);
	bool isParsed() { return parsed; }
};
extern Vec<ARG*> options;

struct INT32R {
	int h, t;
	INT32R() { h = UNDEFINED; t = UNDEFINED; }
	INT32R(int h, int t) { this->h = h; this->t = t; }
};

struct INT64R {
	int64 h, t;
	INT64R() { h = UNDEFINED; t = UNDEFINED; }
	INT64R(int64 h, int64 t) { this->h = h; this->t = t; }
};

struct FP64R {
	double h, t;
	FP64R() { h = 0.0; t = 0.0; }
	FP64R(double h, double t) { this->h = h; this->t = t; }
};

class INT_OPT : public ARG
{
protected:
	INT32R r;
	int val;

public:
	INT_OPT(arg_t a, arg_t x, int val = 0, INT32R r = INT32R(INT32_MIN, INT32_MAX))
		: ARG(a, x, "<int>"), r(r), val(val) {}

	operator int (void) const { return val; }
	operator int& (void) { return val; }
	INT_OPT& operator= (int x) { val = x; return *this; }

	virtual bool parse(arg_t input) {
		arg_t strVal = input;
		if (!eq(strVal, "--") || !eq(strVal, arg) || !eq(strVal, "="))
			return false;
		char* end;
		int tmpVal = strtol(strVal, &end, 10);
		if (end == NULL)
			return false;
		else if (tmpVal > r.t)
			printf("ERROR - maximum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		else if (tmpVal < r.h)
			printf("ERROR - minimum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		val = tmpVal;
		parsed = true;
		return true;
	}

	virtual void help(bool verbose = false) {
		printf("c |  --%-15s = %-8s [", arg, type);
		if (r.h == INT32_MIN) printf("%-5s", "-I32");
		else printf("%-5d", r.h);
		printf(" .. ");
		if (r.t == INT32_MAX) printf("%-5s", "+I32");
		else printf("%5d", r.t);
		printf("] (default: %6d)\n", val);
		if (verbose) {
			printf("c |   %s\n", text);
			printf("c |\n");
		}
	}

	virtual void printArgument() { printf(" %s<%d> ", arg, val); }
};

class INT64_OPT : public ARG
{
protected:
	INT64R r;
	int64  val;

public:
	INT64_OPT(arg_t a, arg_t x, int64 val = 0LL, INT64R r = INT64R(INT64_MIN, INT64_MAX))
		: ARG(a, x, "<int64>"), r(r), val(val) {}

	operator int64 (void) const { return val; }
	operator int64& (void) { return val; }
	INT64_OPT& operator= (int64 x) { val = x; return *this; }

	virtual bool parse(arg_t input) {
		arg_t strVal = input;

		if (!eq(strVal, "--") || !eq(strVal, arg) || !eq(strVal, "="))
			return false;

		char* end;
		int64 tmpVal = strtoll(strVal, &end, 10);

		if (end == NULL)
			return false;
		else if (tmpVal > r.t)
			printf("ERROR - maximum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		else if (tmpVal < r.h)
			printf("ERROR - minimum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		val = tmpVal;
		parsed = true;
		return true;
	}

	virtual void help(bool verbose = false) {
		printf("c |  --%-15s = %-8s [", arg, type);
		if (r.h == INT64_MIN) printf("%-5s", "-I64");
		else printf("%5lld", r.h);
		printf(" .. ");
		if (r.t == INT64_MAX) printf("%5s", "+I64");
		else printf("%5lld", r.t);
		printf("] (default: %6lld)\n", val);
		if (verbose) {
			printf("c |   %s\n", text);
			printf("c |\n");
		}
	}

	virtual void printArgument() { printf(" %s<%lld> ", arg, val); }
};

class DOUBLE_OPT : public ARG
{
	FP64R r;
	double val;

public:
	DOUBLE_OPT(arg_t a, arg_t x, double val = 0.0, FP64R r = FP64R(-INFINITY, INFINITY)) 
	: ARG(a, x, "<double>"), r(r), val(val) {}

	operator double(void) const { return val; }
	operator double& (void) { return val; }
	DOUBLE_OPT& operator=(double x) { val = x; return *this; }

	virtual bool parse(arg_t input) {
		arg_t strVal = input;
		if (!eq(strVal, "--") || !eq(strVal, arg) || !eq(strVal, "="))
			return false;
		char* end;
		double tmpVal = strtod(strVal, &end);
		if (end == NULL)
			return false;
		else if (tmpVal > r.t)
			printf("ERROR - maximum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		else if (tmpVal < r.h)
			printf("ERROR - minimum value exceeded for option \"%s\".\n", arg), exit(EXIT_FAILURE);
		val = tmpVal;
		parsed = true;
		return true;
	}

	virtual void help(bool verbose = false) {
		printf("c |  --%-15s = %-8s [", arg, type);
		if (r.h == -INFINITY) printf("%-6s", "-inf");
		else printf("%6.2f", r.h);
		printf(" .. ");
		if (r.t == INFINITY) printf("%6s", "inf");
		else printf("%6.2f", r.t);
		printf("] (default: %6.2f)\n", val);
		if (verbose) {
			printf("c |   %s\n", text);
			printf("c |\n");
		}
	}

	virtual void printArgument() { printf(" %s<%.2f> ", arg, val); }
};

class STRING_OPT : public ARG
{
	arg_t val;

public:
	STRING_OPT(arg_t a, arg_t x, arg_t val = NULL)
		: ARG(a, x, "<string>"), val(val) {}

	operator arg_t (void) const { return val; }
	operator arg_t& (void) { return val; }
	STRING_OPT& operator=(arg_t x) { val = x; return *this; }

	virtual bool parse(const char *input) {
		arg_t strVal = input;
		if (!eq(strVal, "--") || !eq(strVal, arg) || !eq(strVal, "="))
			return false;
		val = strVal;
		parsed = true;
		return true;
	}

	virtual void help(bool verbose = false) {
		printf("c |  --%-15s = %8s  (default: %s)\n", arg, type, val);
		if (verbose) {
			printf("c |   %s\n", text);
			printf("c |\n");
		}
	}

	virtual void printArgument() { printf(" %s<%s> ", arg, val); }
};

class BOOL_OPT : public ARG
{
	bool val;

public:
	BOOL_OPT(arg_t a, arg_t x, bool val = false)
		: ARG(a, x, "<bool>"), val(val) {}

	operator bool(void) const { return val; }
	operator bool& (void) { return val; }
	BOOL_OPT& operator=(bool b) { val = b; return *this; }
	bool operator!() { return !val;}

	virtual bool parse(const char *input) {
		arg_t strVal = input;
		if (eq(strVal, "-")) {
			bool bVal = !eq(strVal, "no-");
			if (strcmp(strVal, arg) == 0) {
				val = bVal;
				parsed = true;
				return true;
			}
		}
		return false;
	}

	virtual void help(bool verbose = false) {
		printf("c |  -%-10s -no-%-10s", arg, arg);
		printf("                 ");
		printf("(default: %s)\n", val ? "on" : "off");
		if (verbose) {
			printf("c |   %s\n", text);
			printf("c |\n");
		}
	}

	virtual void printArgument() { printf(" %s:%s ", arg, val ? "on" : "off"); }
};

#endif
