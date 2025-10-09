/***********************************************************************[solverinc.cpp]
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

#include "control.hpp"
#include "banner.hpp"
#include "solver.hpp" 

using namespace ParaFROST;

Solver::Solver(const bool& inc)
  : Solver()
{
    incremental = inc;              // override default=false
	initialize(true);
}

Solver::Solver(const bool& inc, const std::string& path)
  : Solver(inc)
{
	formula = path;
	LOGHEADER(1, 5, "Parser")
	if (!parse() || BCP()) { assert(cnfstate == UNSAT), killSolver(); }
	if (opts.parseonly_en) killSolver();
}

void Solver::iallocSpace()
{
	imarks.clear(true);
	if (sp->size() == size_t(inf.maxVar) + 1) return; // avoid allocation if 'maxVar' didn't change
	assert(inf.maxVar);
	LOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	assert(inf.orgVars == inf.maxVar);
	assert(vorg.size() == inf.maxVar + 1);
	assert(V2L(inf.maxVar + 1) == inf.maxDualVars);
	assert(model.lits.size() == inf.maxVar + 1);
	assert(ilevel.size() == ivstate.size());
	assert(ilevel.size() == inf.maxVar + 1);
	assert(imarks.empty());
	inf.orgCls = orgs.size();
	vorg[0] = 0;
	model.lits[0] = 0;
	model.init(vorg);
	SP* newSP = new SP(inf.maxVar + 1);
	newSP->copyFrom(sp);
	delete sp;
	sp = newSP;
	if (opts.proof_en)
		proof.init(sp, vorg);
	ilevel.clear(true);
	ivalue.clear(true);
	iphase.clear(true);
	isource.clear(true);
	ivstate.clear(true);
	LOGDONE(2, 5);
	LOGMEMCALL(this, 2);
}

void Solver::isolve(Lits_t& assumptions)
{
	FAULT_DETECTOR;
	LOGHEADER(1, 5, "Search");
	if (!stats.clauses.original) {
		assert(orgs.empty());
		LOGWARNING("Formula is already SATISFIABLE by elimination");
		return;
	}
	timer.start();
	iallocSpace();
	iunassume();
	assert(IS_UNSOLVED(cnfstate));
	if (BCP()) {
		LOG2(2, " Incremental formula has a contradiction on top level");
		learnEmpty();
	}
	else {
		initLimits();
		iassume(assumptions);
		if (verbose == 1) printTable();
		if (canPreSimplify()) simplify();
		if (IS_UNSOLVED(cnfstate)) {
			MDMInit();
			while (IS_UNSOLVED(cnfstate) && !interrupted()) {
				LOGDL(this, 3);
				if (BCP()) analyze();
				else if (!inf.unassigned) cnfstate = SAT;
				else if (canReduce()) reduce();
				else if (canRestart()) restart();
				else if (canRephase()) rephase();
				else if (canSigmify()) simplify();
				else if (canProbe()) probe();
				else if (canMMD()) MDM();
				else idecide();
			}
		}
	}
	timer.stop(), stats.time.solve += timer.cpuTime();
	wrapup();
}