/***********************************************************************[transfer.cu]
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
#include "options.cuh"
#include "timer.cuh"

using namespace ParaFROST;

void Solver::newBeginning() 
{
	assert(opts.sigma_en || opts.sigma_live_en);
	assert(wt.empty());
	assert(orgs.empty());
	assert(learnts.empty());
	assert(inf.numClauses <= hcnf->size());
	if (gopts.unified_access) assert(hcnf == NULL), hcnf = cnf;
	assert(!hcnf->empty());
	cm.init(hcnf->data().size);
	if (!gopts.unified_access) {
		SYNC(streams[0]); SYNC(streams[1]); // sync CNF caching
		if (gopts.profile_gpu) cutimer.stop(streams[1]), stats.sigma.time.io += cutimer.gpuTime();
	}
	cacheResolved(streams[2]);
	writeBackCNF();
	SYNCALL;
	if (gopts.unified_access) {
		hcnf = NULL;
		if (gopts.profile_gpu) cutimer.stop(), stats.sigma.time.io += cutimer.gpuTime();
	}
	else cumm.breakMirror(), hcnf = NULL;
}

void Solver::markEliminated(const cudaStream_t& _s)
{
	assert(vars->isEliminatedCached);

#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
	uint32 unassigned = inf.unassigned;
#endif
	
	SYNC(_s);

	forall_variables(v) {
		if (vars->cachedEliminated[v] && !IS_FORCED(vars->cachedEliminated[v]))
			markEliminated(v);
	}

#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
	assert(unassigned >= inf.unassigned);
#endif
}

void Solver::cacheUnits(const cudaStream_t& stream) 
{
	SYNC(stream);
	if ((vars->nUnits = vars->tmpUnits.size()))
		CHECK(cudaMemcpyAsync(vars->cachedUnits, vars->unitsData, vars->nUnits * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}

void Solver::cacheEliminated(const cudaStream_t& stream)
{
	if (vars->isEliminatedCached) return;
	CHECK(cudaMemcpyAsync(vars->cachedEliminated, vars->eliminated, inf.maxVar + 1, cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
	vars->isEliminatedCached = true;
}

void Solver::cacheNumUnits(const cudaStream_t& stream)
{
	CHECK(cudaMemcpyAsync(&vars->tmpUnits, vars->units, sizeof(cuVecU), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}

void Solver::cacheResolved(const cudaStream_t& stream)
{
	cuVecU tmpObj;
	CHECK(cudaMemcpy(&tmpObj, vars->resolved, sizeof(cuVecU), cudaMemcpyDeviceToHost));
	const uint32* devStart = tmpObj.data();
	const uint32 devSize = tmpObj.size();
	if (!devSize) return;
	assert(devStart);
	const uint32 off = model.resolved.size();
	model.resolved.resize(off + devSize);
	uint32* start = model.resolved + off;
	CHECK(cudaMemcpyAsync(start, devStart, devSize * sizeof(uint32), cudaMemcpyDeviceToHost, stream));
	if (gopts.sync_always) SYNC(stream);
}