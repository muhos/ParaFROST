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

Solver::Solver() :
	sp(NULL)
	, vsids(VSIDS_CMP(activity))
	, vschedule(SCORS_CMP(this))
	, bumped(0)
	, conflict(NOREF)
	, ignore(NOREF)
	, cnfstate(UNSOLVED_M)
	, intr(false)
	, stable(false)
	, probed(false)
	, incremental(true)
	, vars(NULL)
	, ot(NULL)
	, cnf(NULL)
	, hcnf(NULL)
	, cuproof(cumm, proof)
	, streams(NULL)
	, mapped(false)
	, compacted(false)
	, flattened(false)
	, phase(0)
	, nForced(0)
	, simpstate(AWAKEN_SUCC)
	, devCount(0)
	, termCallbackState(NULL)
	, learnCallbackState(NULL)
	, learnCallbackBuffer(NULL)
	, termCallback(NULL)
	, learnCallback(NULL)
{
	LOGHEADER(1, 5, "Banner");
	LOGFANCYBANNER(version());
	LOGHEADER(1, 5, "Build")
	getCPUInfo(stats.sysmem);
	getBuildInfo();
	initSolver();
	size_t _gfree = 0, _gpenalty = 0;
	if (!(devCount = getGPUInfo(_gfree, _gpenalty))) {
		LOGERRORN("no GPU(s) available that support CUDA");
		killSolver();
	}
	if (_gfree <= 200 * MBYTE) {
		LOGWARNING("not enough GPU memory (free = %zd MB): skip simplifier", _gfree / MBYTE);
		opts.sigma_en = opts.sigma_live_en = false;
	}
	else cumm.init(_gfree, _gpenalty);
	if (opts.sigma_en || opts.sigma_live_en) { optSimp(), createStreams(); }
}

void Solver::iallocSpace()
{
	imarks.clear(true);
	if (sp->size() == size_t(inf.maxVar) + 1) return; // avoid allocation if 'maxVar' didn't change
	assert(inf.maxVar);
	LOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	assert(inf.orgVars == inf.maxVar);
	assert(vorg.size() == inf.maxVar + 1);
	assert(V2L(inf.maxVar + 1) == inf.nDualVars);
	assert(model.lits.size() == inf.maxVar + 1);
	assert(ilevel.size() == ivstate.size());
	assert(ilevel.size() == inf.maxVar + 1);
	assert(imarks.empty());
	inf.nOrgCls = orgs.size();
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
	assert(UNSOLVED(cnfstate));
	if (BCP()) {
		LOG2(2, " Incremental formula has a contradiction on top level");
		learnEmpty();
	}
	else {
		initLimits();
		iassume(assumptions);
		if (verbose == 1) printTable();
		if (canPreSimplify()) simplify();
		if (UNSOLVED(cnfstate)) {
			MDMInit();
			while (UNSOLVED(cnfstate) && !interrupted()) {
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
	timer.stop(), timer.solve += timer.cpuTime();
	wrapup();
}