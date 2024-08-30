/***********************************************************************[elimination.cu]
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

#include "solver.hpp"
#include <cub/device/device_select.cuh>

using namespace cub;
using namespace ParaFROST;

void Solver::VE()
{
	if (opts.ve_en) {
		if (interrupted()) killSolver();
		LOG2(2, " Eliminating variables..");
		inf.n_del_vars_after = 0;
		veAsync(cnf, ot, vars, streams, cuproof.gpuStream(), cumm, cuhist, stats.sigma.calls > 1);
		postVE();
		LOGREDALL(this, 2, "BVE Reductions");
	}
}

void Solver::postVE()
{
	size_t bytes = 0;
	uint32* tmpmem = NULL;
	DeviceSelect::If(NULL, bytes, vars->eligible, vars->electedData, tmpmem, vars->numElected, COMPACT_VARS());
	tmpmem = (uint32*) ((bytes > cumm.scatterCap()) ? cacher.allocate(bytes) : cumm.scatter());
	if (!tmpmem) throw MEMOUTEXCEPTION();
	DeviceSelect::If(tmpmem + 1, bytes,  vars->eligible, vars->electedData, vars->electedSize, vars->numElected, COMPACT_VARS());
	veResizeCNFAsync(cnf, cuhist);
	if (bytes > cumm.scatterCap()) {
		assert(tmpmem != (uint32*)cumm.scatter());
		cacher.deallocate(tmpmem);
	}
}

void Solver::SUB()
{
	if (opts.sub_en || opts.ve_plus_en) {
		if (interrupted()) killSolver();
		LOG2(2, " Eliminating (self)-subsumptions..");
		subAsync(cnf, ot, vars, cuproof.gpuStream());
		LOGREDCL(this, 2, "SUB Reductions");
	}
}

void Solver::BCE()
{
	if (opts.bce_en) {
		if (interrupted()) killSolver();
		if (!vars->numElected) return;
		LOG2(2, " Eliminating blocked clauses..");
		bceAsync(cnf, ot, vars, cuproof.gpuStream());
		LOGREDCL(this, 2, "BCE Reductions");
	}
}

void Solver::ERE()
{
	if (opts.ere_en) {
		if (interrupted()) killSolver();
		if (!vars->numElected) return;
		LOG2(2, " Eliminating redundances..");
		ereCls = inf.nClauses;
		cacheEliminated(streams[5]);
		ereAsync(cnf, ot, vars, cuproof.gpuStream());
		LOGREDCL(this, 2, "ERE Reductions");
		cuproof.cacheProof(0);
		cuproof.writeProof(0);
	}
}