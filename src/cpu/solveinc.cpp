/***********************************************************************[solveinc.cpp]
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

#include "control.h"
#include "solve.h" 
#include "version.h"

using namespace pFROST;

ParaFROST::ParaFROST() :
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
	, mapped(false)
	, simpstate(AWAKEN_SUCC)
{
	PFNAME("ParaFROST (Parallel Formal Reasoning On Satisfiability)", version());
	getCPUInfo(stats.sysmem);
	getBuildInfo();
	initSolver();
}

void ParaFROST::iallocSpace()
{
	imarks.clear(true);
	if (sp->size() == size_t(inf.maxVar) + 1) return; // avoid allocation if 'maxVar' didn't change
	assert(inf.maxVar);
	assert(inf.orgVars == inf.maxVar);
	assert(vorg.size() == inf.maxVar + 1);
	assert(V2L(inf.maxVar + 1) == inf.nDualVars);
	assert(model.lits.size() == inf.maxVar + 1);
	assert(ilevel.size() == ivstate.size());
	assert(ilevel.size() == inf.maxVar + 1);
	assert(imarks.empty());
	vorg[0] = 0;
	model.lits[0] = 0;
	model.init(vorg);
	PFLOGN2(2, " Allocating fixed memory for %d variables..", inf.maxVar);
	SP* newSP = new SP(inf.maxVar + 1);
	newSP->initSaved(opts.polarity);
	newSP->copyFrom(sp);
	delete sp;
	sp = newSP;
	if (opts.proof_en)
		proof.init(sp, vorg);
	ilevel.clear(true);
	ivalue.clear(true);
	ivstate.clear(true);
	PFLDONE(2, 5);
	PFLMEMCALL(this, 2);
}

void ParaFROST::isolve(Lits_t& assumptions)
{
	timer.start();
	iallocSpace();
	iunassume();
	assert(UNSOLVED(cnfstate));
	if (BCP()) {
		PFLOG2(2, " Incremental formula has a contradiction on top level");
		learnEmpty();
	}
	else {
		initLimits();
		iassume(assumptions);
		if (verbose == 1) printTable();
		if (canPreSigmify()) sigmify();
		if (UNSOLVED(cnfstate)) {
			PFLOG2(2, "-- Incremental CDCL search started..");
			MDMInit();
			while (UNSOLVED(cnfstate) && !interrupted()) {
				PFLDL(this, 3);
				if (BCP()) analyze();
				else if (!inf.unassigned) cnfstate = SAT;
				else if (canReduce()) reduce();
				else if (canRestart()) restart();
				else if (canRephase()) rephase();
				else if (canSigmify()) sigmify();
				else if (canProbe()) probe();
				else if (canMMD()) MDM();
				else idecide();
			}
			PFLOG2(2, "-- Incremental CDCL search completed successfully");
		}
	}
	timer.stop(), timer.solve += timer.cpuTime();
	wrapup();
}