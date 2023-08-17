/***********************************************************************[statistics.cpp]
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

#include "solve.h" 
#include "control.h"
using namespace ParaFROST;

void Solver::report()
{
	if (opts.report_en) {
		PFLOG0("");
		if (opts.sigma_en || opts.sigma_live_en) {
			PFLOG1("\t\t\t%sSimplifier Report%s", CREPORT, CNORMAL);
			if (!opts.profile_simp)
				PFLOG1(" %sSimplifier time        : %s%-16.3f  sec%s", CREPORT, CREPORTVAL, timer.simp, CNORMAL);
			else {
				PFLOG1(" %s - Var ordering        : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.vo, CNORMAL);
				PFLOG1(" %s - CNF compact         : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.gc, CNORMAL);
				PFLOG1(" %s - OT  creation        : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.cot, CNORMAL);
				PFLOG1(" %s - OT  sorting         : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.sot, CNORMAL);
				PFLOG1(" %s - OT  reduction       : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.rot, CNORMAL);
				PFLOG1(" %s - BVE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.ve, CNORMAL);
				PFLOG1(" %s - HSE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.sub, CNORMAL);
				PFLOG1(" %s - BCE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.bce, CNORMAL);
				PFLOG1(" %s - ERE                 : %s%-16.2f  ms%s", CREPORT, CREPORTVAL, timer.ere, CNORMAL);
			}
			PFLOG1(" %sSigmifications         : %s%-10d%s", CREPORT, CREPORTVAL, stats.sigma.calls, CNORMAL);
#ifdef STATISTICS
			PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.variables + stats.units.forced, CNORMAL);
			PFLOG1(" %s  Resolutions          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.resolutions, CNORMAL);
			PFLOG1(" %s  Forced units         : %s%-10d%s", CREPORT, CREPORTVAL, stats.units.forced, CNORMAL);
			PFLOG1(" %s  Pure literals        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.pures, CNORMAL);
			PFLOG1(" %s  If-Then-Else         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.ites, CNORMAL);
			PFLOG1(" %s  Inverter             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.inverters, CNORMAL);
			PFLOG1(" %s  AND-OR               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.andors, CNORMAL);
			PFLOG1(" %s  Alien                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.aliens, CNORMAL);
			PFLOG1(" %s  XOR                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.bve.xors, CNORMAL);
			PFLOG1(" %s Removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.clauses, CNORMAL);
			PFLOG1(" %s  Subsumed             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.sub.subsumed, CNORMAL);
			PFLOG1(" %s  Strengthened         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.sub.strengthened, CNORMAL);
			PFLOG1(" %s Tried redundancies    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.tried, CNORMAL);
			PFLOG1(" %s  Original removed     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.orgs, CNORMAL);
			PFLOG1(" %s  Learnt removed       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.ere.learnts, CNORMAL);
			PFLOG1(" %s Removed literals      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.literals, CNORMAL);
#else
			PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.variables + stats.units.forced, CNORMAL);
			PFLOG1(" %s  Forced units         : %s%-10d%s", CREPORT, CREPORTVAL, stats.units.forced, CNORMAL);
			PFLOG1(" %s Removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.clauses, CNORMAL);
			PFLOG1(" %s Removed literals      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.sigma.all.literals, CNORMAL);
#endif
		}
		PFLOG1("\t\t\t%sSolver Report%s", CREPORT, CNORMAL);
		PFLOG1(" %sSolver time            : %s%-16.3f  sec%s", CREPORT, CREPORTVAL, timer.solve, CNORMAL);
		PFLOG1(" %sSystem memory          : %s%-16.3f  MB%s", CREPORT, CREPORTVAL, ratio(double(sysMemUsed()), double(MBYTE)), CNORMAL);
		PFLOG1(" %sFormula                : %s%-s%s", CREPORT, CREPORTVAL, formula.path.c_str(), CNORMAL);
		PFLOG1(" %s Size                  : %s%-16.3f  MB%s", CREPORT, CREPORTVAL, ratio(double(formula.size), double(MBYTE)), CNORMAL);
		PFLOG1(" %s Units                 : %s%-10d%s", CREPORT, CREPORTVAL, formula.units, CNORMAL);
		PFLOG1(" %s Binaries              : %s%-10d%s", CREPORT, CREPORTVAL, formula.binaries, CNORMAL);
		PFLOG1(" %s Ternaries             : %s%-10d%s", CREPORT, CREPORTVAL, formula.ternaries, CNORMAL);
		PFLOG1(" %s Larger                : %s%-10d%s", CREPORT, CREPORTVAL, formula.large, CNORMAL);
		PFLOG1(" %s Max clause size       : %s%-10d%s", CREPORT, CREPORTVAL, formula.maxClauseSize, CNORMAL);
		PFLOG1(" %s C2V ratio             : %s%-10.3f%s", CREPORT, CREPORTVAL, formula.c2v, CNORMAL);
		if (opts.proof_en)
			PFLOG1(" %s Proof lines           : %s%-16zd%s", CREPORT, CREPORTVAL, proof.numClauses(), CNORMAL);
		PFLOG1(" %sAutarky calls          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.autarky.calls, CNORMAL);
		PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.autarky.eliminated, CNORMAL);
		PFLOG1(" %sBacktracks             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.chrono + stats.backtrack.nonchrono, CNORMAL);
		PFLOG1(" %s Chronological         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.chrono, CNORMAL);
		PFLOG1(" %s Non-Chronological     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.backtrack.nonchrono, CNORMAL);
		PFLOG1(" %s Trail reuses          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.reuses, CNORMAL);
		PFLOG1(" %sConflicts              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.conflicts, CNORMAL);
		PFLOG1(" %s OTF strengthened      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.strengthenedfly, CNORMAL);
		PFLOG1(" %s OTF subsumed          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.subsumedfly, CNORMAL);
		PFLOG1(" %s Learnt OTF subsumes   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subtried, CNORMAL);
		PFLOG1(" %s Learnt OTF subsumed   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.learntfly, CNORMAL);
		PFLOG1(" %s Learnt units          : %s%-10d%s", CREPORT, CREPORTVAL, stats.units.learnt, CNORMAL);
		PFLOG1(" %s Learnt literals       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.literals.learnt, CNORMAL);
#ifdef STATISTICS
		PFLOG1(" %s Minimized literals    : %s%2.2f %%%s", CREPORT, CREPORTVAL, percent((double)stats.minimize.before - stats.minimize.after, (double)stats.minimize.before), CNORMAL);
#endif
		PFLOG1(" %sDeduplications         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.calls, CNORMAL);
		PFLOG1(" %s Hyper unaries         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.hyperunary, CNORMAL);
		PFLOG1(" %s Duplicated binaries   : %s%-16lld%s", CREPORT, CREPORTVAL, stats.debinary.binaries, CNORMAL);
		PFLOG1(" %sDecompositions         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.calls, CNORMAL);
		PFLOG1(" %s SCCs                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.scc, CNORMAL);
		PFLOG1(" %s Hyper unaries         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.hyperunary, CNORMAL);
		PFLOG1(" %s Removed variables     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.variables, CNORMAL);
		PFLOG1(" %s Removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decompose.clauses, CNORMAL);
		PFLOG1(" %sHyper binary resolves  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.resolutions, CNORMAL);
		PFLOG1(" %s Added binaries        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.resolvents, CNORMAL);
		PFLOG1(" %sHyper ternary resolves : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.resolutions, CNORMAL);
		PFLOG1(" %s Added binaries        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.binaries, CNORMAL);
		PFLOG1(" %s Added ternaries       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.ternaries, CNORMAL);
		PFLOG1(" %s Subsumed ternaries    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.binaries * 2, CNORMAL);
		PFLOG1(" %sReduces                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.reduces, CNORMAL);
#ifdef STATISTICS
		PFLOG1(" %s Removed binaries      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.binary.reduced, CNORMAL);
		PFLOG1(" %s Removed ternaries     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.ternary.reduced, CNORMAL);
#endif
		PFLOG1(" %sRestarts               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.restart.all, CNORMAL);
		PFLOG1(" %s Stable restarts       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.restart.stable, CNORMAL);
		PFLOG1(" %s Stable modes          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.stablemodes, CNORMAL);
		PFLOG1(" %s Unstable modes        : %s%-16lld%s", CREPORT, CREPORTVAL, stats.unstablemodes, CNORMAL);
		PFLOG1(" %sRephases               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.all, CNORMAL);
#ifdef STATISTICS
		PFLOG1(" %s Original              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.org, CNORMAL);
		PFLOG1(" %s Random                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.random, CNORMAL);
		PFLOG1(" %s Invert                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.inv, CNORMAL);
		PFLOG1(" %s Flip                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.flip, CNORMAL);
		PFLOG1(" %s Best                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.rephase.best, CNORMAL);
#endif
		PFLOG1(" %s Walk                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.calls, CNORMAL);
		PFLOG1(" %sRecyclings             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.soft + stats.recycle.hard, CNORMAL);
		PFLOG1(" %s Soft                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.soft, CNORMAL);
		PFLOG1(" %s Hard                  : %s%-16lld%s", CREPORT, CREPORTVAL, stats.recycle.hard, CNORMAL);
		PFLOG1(" %sProbes calls           : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.calls, CNORMAL);
		PFLOG1(" %s Rounds                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.rounds, CNORMAL);
		PFLOG1(" %s Probed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.probed, CNORMAL);
		PFLOG1(" %s Failed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.failed, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probeticks, CNORMAL);
		PFLOG1(" %sTransitive calls       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.probe.calls, CNORMAL);
		PFLOG1(" %s Probed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.probed, CNORMAL);
		PFLOG1(" %s Failed                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.failed, CNORMAL);
		PFLOG1(" %s Transitives           : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitive.removed, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.transitiveticks, CNORMAL);
#ifdef STATISTICS
		PFLOG1(" %sShrinks                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.calls, CNORMAL);
		PFLOG1(" %s removed clauses       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.clauses, CNORMAL);
		PFLOG1(" %s removed literals      : %s%-16lld%s", CREPORT, CREPORTVAL, stats.shrink.literals, CNORMAL);
#endif
		PFLOG1(" %sSubsume calls          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.calls, CNORMAL);
		PFLOG1(" %s Checks                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.checks, CNORMAL);
		PFLOG1(" %s Subsumed              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.subsumed, CNORMAL);
		PFLOG1(" %s Strengthened          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.subsume.strengthened, CNORMAL);
		PFLOG1(" %sSearch decisions       : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decisions.single, CNORMAL);
		PFLOG1(" %s Propagations          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.searchprops, CNORMAL);
		PFLOG1(" %s Ticks                 : %s%-16lld%s", CREPORT, CREPORTVAL, stats.searchticks, CNORMAL);
		PFLOG1(" %sMappings               : %s%-10d%s", CREPORT, CREPORTVAL, stats.mappings, CNORMAL);
		PFLOG1(" %sMDM calls              : %s%-10d%s", CREPORT, CREPORTVAL, stats.mdm.calls, CNORMAL);
		PFLOG1(" %s Walks                 : %s%-10d%s", CREPORT, CREPORTVAL, stats.mdm.walks, CNORMAL);
		PFLOG1(" %s VMTF uses             : %s%-10d%s", CREPORT, CREPORTVAL, stats.mdm.vmtf, CNORMAL);
		PFLOG1(" %s VSIDS uses            : %s%-10d%s", CREPORT, CREPORTVAL, stats.mdm.vsids, CNORMAL);
		PFLOG1(" %s All decisions         : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decisions.multiple, CNORMAL);
		PFLOG1(" %s Assumed decisions     : %s%-16lld%s", CREPORT, CREPORTVAL, stats.decisions.massumed, CNORMAL);
		PFLOG1(" %s Last-made decisions   : %s%-10d%s", CREPORT, CREPORTVAL, last.mdm.decisions, CNORMAL);
		PFLOG1(" %sVivification checks    : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.checks, CNORMAL);
		PFLOG1(" %s Vivified              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.vivified, CNORMAL);
#ifdef STATISTICS
		PFLOG1(" %s Reused                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.reused, CNORMAL);
		PFLOG1(" %s Assumed               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.assumed, CNORMAL);
		PFLOG1(" %s Implied               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.implied, CNORMAL);
		PFLOG1(" %s Subsumed              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.subsumed, CNORMAL);
		PFLOG1(" %s Strengthened          : %s%-16lld%s", CREPORT, CREPORTVAL, stats.vivify.strengthened, CNORMAL);
#endif
		PFLOG1(" %sWalk calls             : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.calls, CNORMAL);
		PFLOG1(" %s Checks                : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.checks, CNORMAL);
		PFLOG1(" %s Minimum               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.minimum, CNORMAL);
		PFLOG1(" %s Flipped               : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.flipped, CNORMAL);
		PFLOG1(" %s Improved              : %s%-16lld%s", CREPORT, CREPORTVAL, stats.walk.improved, CNORMAL);
	}
}