#include "pfsimp.h" 
#include "pfsolve.h"
#include "pfsort.h"
//=======================================//
//		 GEARS defined members           //
//=======================================//
void GEARS::mGPU_other_alloc(int32_t id)
{
	slaves[id].streams = new cudaStream_t[MAX_STREAMS];
	for (int s = 0; s < MAX_STREAMS; s++) CHECK(cudaStreamCreate(&slaves[id].streams[s]));
	CHECK(cudaMallocManaged((void **)&slaves[id].GPUStats, sizeof(GPU_STATS)));
	CHECK(cudaMallocManaged((void **)&slaves[id].pv->PVs, nOrgVars() * sizeof(uint32)));
	CHECK(cudaMallocManaged((void **)&slaves[id].pv->elim, nOrgVars() * sizeof(bool)));
	CHECK(cudaMallocManaged((void **)&slaves[id].nLits_after, sizeof(uint32)));
	slaves[id].gMemCons += sizeof(GPU_STATS) + nOrgVars() * 2.0;
}

void GEARS::mGPU_init_simp(int32_t id)
{
	cudaMemset(pv->isElectedCl, 0, nOrgClauses());
	slaves[id].pv->nParVars = 0;
	slaves[id].pv->nVarsPerGPU = 0;
	cudaMemset(slaves[id].pv->PVs, 0, nOrgVars() * sizeof(uint32));
	cudaMemset(slaves[id].pv->elim, 0, nOrgVars());
	// initialize CPU cnf_stats
	slaves[id].cnf_stats.n_org_vars = nOrgVars();
	slaves[id].cnf_stats.n_cls_after = 0;
	slaves[id].cnf_stats.n_lits_after = 0;
	slaves[id].cnf_stats.n_del_vars = 0;
}

void GEARS::mGPU_distribution_init(bool ve)
{
	/* distribute CNF & OT among available GPUs */
	if ((int32_t)pv->nParVars < devCount) {
		cout << "c | Cannot proceed with multiGPUs option --> reason: pv->nParVars < devCount" << endl;
		cout << "c | Action --> exit" << endl;
		exit(EXIT_FAILURE);
	}
	/* allocate memory for slaves GPUs */
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id)); // set the current context to GPU id
		mGPU_other_alloc(id);
		mGPU_init_simp(id);
		mGPU_OT_alloc(id);
	}
	/* grant access for slave GPUs */
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		cudaDeviceEnablePeerAccess(MASTER_GPU, 0);
		int accessible = 0;
		cudaDeviceCanAccessPeer(&accessible, id, MASTER_GPU);
		assert(accessible);
	}
	/* distribute elected variables (thanks to Mark van den Brand for the org. idea) */
	/* update: alternating between GPUs while distribting MCVs have better performance  */
	uVector2D mGPU_elected_vars(devCount);
	register int32_t MCV_idx = pv->nParVars - 1;
	register int32_t LCV_idx = 0;
	register int32_t v_idx = MCV_idx;
	register bool flip = false;
	while (v_idx >= devCount) {
		if (flip)
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[(devCount - 1) - id].push(pv->PVs[v_idx - id]);
		else
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[id].push(pv->PVs[v_idx - id]);
		flip = !flip;
		v_idx -= devCount;
	}
	while (v_idx >= LCV_idx) mGPU_elected_vars[MASTER_GPU].push(pv->PVs[v_idx--]);
	register uint32 checksum = 0;
	for (int32_t id = MASTER_GPU; id < devCount; id++) checksum += (uint32)mGPU_elected_vars[id].size();
	assert(checksum == pv->nParVars);
	// rewrite the master list based on redistribution
	pv->nVarsPerGPU = (uint32)mGPU_elected_vars[MASTER_GPU].size();
	for (uint32 i = 0; i < pv->nVarsPerGPU; i++) pv->PVs[i] = mGPU_elected_vars[MASTER_GPU][i];
	register uint32 offset = 0;
	for (int32_t id = 1; id < devCount; id++) {
		offset += (uint32)mGPU_elected_vars[id - 1].size();
		for (uint32 i = 0; i < mGPU_elected_vars[id].size(); i++) pv->PVs[offset + i] = mGPU_elected_vars[id][i];
	}
	// transfer the elected vars across the slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].pv->nParVars = slaves[id].pv->nVarsPerGPU = (uint32)mGPU_elected_vars[id].size();
		for (uint32 i = 0; i < slaves[id].pv->nParVars; i++) slaves[id].pv->PVs[i] = mGPU_elected_vars[id][i];
	}
	if (ve) {
		/* Calculating resolvents upper bound on master GPU */
		calc_up(slaves, cnf, ot, pv, mGPUStats, devCount, devProps, streams);
		// Verify results
		assert(up_validity(cnf, slaves, ot, pv, devCount, devProps));
	}
	else { // HSE on all pv. vars
		/* Calculating org. clauses upper bound on master GPU */
		calc_up(slaves, cnf, ot, pv, mGPUStats, devCount, devProps, streams, 1);
		// Verify results
		assert(up_validity(cnf, slaves, ot, pv, devCount, devProps, 1));
	}
	assert(cnf_stats.max_org_cls > 0 && cnf_stats.max_org_lits > 0);
	/* Allocate memory for CNF on slave GPUs */
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		// re-initialize GPU cnf_stats
		assert(slaves[id].cnf_stats.max_org_cls > 0 && slaves[id].cnf_stats.max_org_lits > 0);
		slaves[id].cnf_stats.n_org_cls = slaves[id].cnf_stats.max_org_cls;
		slaves[id].cnf_stats.n_org_lits = slaves[id].cnf_stats.max_org_lits;
		slaves[id].GPUStats->org_nClauses = slaves[id].cnf_stats.max_org_cls;
		slaves[id].GPUStats->org_nLits = slaves[id].cnf_stats.max_org_lits;
		// allocation based on #clauses & #literals calculated by upper-bound kernel
		if (ve) slaves[id].cnf_stats.max_added_cl_width = MAX_ADDED_CL_LEN; // estimated added clause width
		else {
			slaves[id].cnf_stats.max_added_cl_width = 0;
			slaves[id].cnf_stats.max_added_cls = 0;
		}
		if (!mGPU_CNF_alloc_slave(id)) {
			cout << "c | Cannot proceed with BVE on multi-GPUs." << endl;
			cudaDeviceReset();
			exit(EXIT_FAILURE);
		}
	}
	/* broadcast CNF to all slaves */
	if (verbose >= 2) cout << "c | Broadcasting CNF to all slaves...";
	broadcast_cnf(slaves, cnf, ot, devCount, devProps, timer);
	if (verbose >= 2) cout << " ==> done." << endl;
	/* Reallocate memory for CNF on master GPU */
	cudaSetDevice(MASTER_GPU);
	if (ve) cnf_stats.max_added_cl_width = MAX_ADDED_CL_LEN; // estimated added clause width
	else {
		cnf_stats.max_added_cl_width = 0;
		cnf_stats.max_added_cls = 0;
	}
	if (!mGPU_CNF_realloc_master()) {
		cout << "c | Cannot proceed with BVE on multi-GPUs." << endl;
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
	// create occurrence tables for all GPUs
	create_ot(slaves, cnf, ot, devCount, devProps, timer);
	// free system memory
	mGPU_elected_vars.clear(true);
}

bool GEARS::mGPU_distribution(int8_t phase)
{
	/* distribute CNF & OT among available GPUs */
	if ((int32_t)pv->nParVars < devCount) {
		cout << "c | Cannot proceed with multiGPUs option --> reason: pv->nParVars < devCount" << endl;
		cout << "c | Action --> exit" << endl;
		exit(EXIT_FAILURE);
	}
	/* reallocate OT memory on slaves */
	// broadcast occurrences from master GPU to all slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].v_score.clear();
		slaves[id].v_score.resize(nOrgVars());
		for (uint32 v = 0; v < nOrgVars(); v++) {
			slaves[id].v_score[v].first = v_score[v].first;
			slaves[id].v_score[v].second = v_score[v].second;
		}
		mGPU_OT_realloc_slave(id);
	}
	uVector2D mGPU_elected_vars(devCount);
	register int32_t MCV_idx = pv->nParVars - 1;
	register int32_t LCV_idx = 0;
	register int32_t v_idx = MCV_idx;
	register bool flip = false;
	while (v_idx >= devCount) {
		if (flip)
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[(devCount - 1) - id].push(pv->PVs[v_idx - id]);
		else
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[id].push(pv->PVs[v_idx - id]);
		flip = !flip;
		v_idx -= devCount;
	}
	while (v_idx >= LCV_idx) mGPU_elected_vars[MASTER_GPU].push(pv->PVs[v_idx--]);
	register uint32 checksum = 0;
	for (int32_t id = MASTER_GPU; id < devCount; id++) checksum += (uint32)mGPU_elected_vars[id].size();
	assert(checksum == pv->nParVars);
	// rewrite the master list based on redistribution
	pv->nVarsPerGPU = (uint32)mGPU_elected_vars[MASTER_GPU].size();
	for (uint32 i = 0; i < pv->nVarsPerGPU; i++) pv->PVs[i] = mGPU_elected_vars[MASTER_GPU][i];
	register uint32 offset = 0;
	for (int32_t id = 1; id < devCount; id++) {
		offset += (uint32)mGPU_elected_vars[id - 1].size();
		for (uint32 i = 0; i < mGPU_elected_vars[id].size(); i++) pv->PVs[offset + i] = mGPU_elected_vars[id][i];
	}
	// transfer the elected vars across the slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].pv->nParVars = slaves[id].pv->nVarsPerGPU = (uint32)mGPU_elected_vars[id].size();
		for (uint32 i = 0; i < slaves[id].pv->nParVars; i++) slaves[id].pv->PVs[i] = mGPU_elected_vars[id][i];
	}
	/* Calculating resolvents upper bound on master GPU */
	calc_up(slaves, cnf, ot, pv, mGPUStats, devCount, devProps, streams);
	// Verify results
	assert(up_validity(cnf, slaves, ot, pv, devCount, devProps));
	assert(cnf_stats.max_org_cls > 0 && cnf_stats.max_org_lits > 0);
	/* Allocate memory for CNF on slave GPUs */
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		// re-initialize GPU cnf_stats
		assert(slaves[id].cnf_stats.max_org_cls > 0 && slaves[id].cnf_stats.max_org_lits > 0);
		slaves[id].cnf_stats.n_org_cls = slaves[id].cnf_stats.max_org_cls;
		slaves[id].cnf_stats.n_org_lits = slaves[id].cnf_stats.max_org_lits;
		slaves[id].GPUStats->org_nClauses = slaves[id].cnf_stats.max_org_cls;
		slaves[id].GPUStats->org_nLits = slaves[id].cnf_stats.max_org_lits;
		slaves[id].cnf_stats.max_added_cl_width = MAX_ADDED_CL_LEN; // estimated added clause width
		if (!mGPU_CNF_realloc_slave(id)) return false;  
	}
	/* broadcast CNF to all slaves */
	if (verbose >= 2) cout << "c | Broadcasting CNF to all slaves...";
	broadcast_cnf(slaves, cnf, ot, devCount, devProps, timer, 1);
	if (verbose >= 2) cout << " ==> done." << endl;
	/* Reallocate memory for CNF on master GPU */
	cudaSetDevice(MASTER_GPU);
	cnf_stats.max_added_cl_width = MAX_ADDED_CL_LEN; // estimated added clause width
	if (!mGPU_CNF_realloc_master(false)) return false;
	// create occurrence tables for all GPUs
	create_ot(slaves, cnf, ot, devCount, devProps, timer);
	// free system memory
	mGPU_elected_vars.clear(true);
	return true;
}

bool GEARS::mGPU_redistribution()
{
	/* distribute CNF & OT among available GPUs */
	if ((int32_t)pv->nParVars < devCount) {
		cout << "c | Cannot proceed with multiGPUs option --> reason: pv->nParVars < devCount" << endl;
		cout << "c | Action --> exit" << endl;
		exit(EXIT_FAILURE);
	}
	/* reallocate OT memory on slaves */
	// broadcast occurrences from master GPU to all slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].v_score.clear();
		slaves[id].v_score.resize(nOrgVars());
		for (uint32 v = 0; v < nOrgVars(); v++) {
			slaves[id].v_score[v].first = v_score[v].first;
			slaves[id].v_score[v].second = v_score[v].second;
		}
		mGPU_OT_realloc_slave(id);
	}
	/* redistribute elected variables */
	uVector2D mGPU_elected_vars(devCount);
	register int32_t MCV_idx = pv->nParVars - 1;
	register int32_t LCV_idx = 0;
	register int32_t v_idx = MCV_idx;
	register bool flip = false;
	while (v_idx >= devCount) {
		if (flip)
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[(devCount - 1) - id].push(pv->PVs[v_idx - id]);
		else 
			for (int32_t id = MASTER_GPU; id < devCount; id++) mGPU_elected_vars[id].push(pv->PVs[v_idx - id]);
		flip = !flip;
		v_idx -= devCount;
	}
	while (v_idx >= LCV_idx) mGPU_elected_vars[MASTER_GPU].push(pv->PVs[v_idx--]);
	register uint32 checksum = 0;
	for (int32_t id = MASTER_GPU; id < devCount; id++) checksum += (uint32)mGPU_elected_vars[id].size();
	assert(checksum == pv->nParVars);
	// rewrite the master list based on redistribution
	pv->nVarsPerGPU = (uint32) mGPU_elected_vars[MASTER_GPU].size();
	for (uint32 i = 0; i < pv->nVarsPerGPU; i++) pv->PVs[i] = mGPU_elected_vars[MASTER_GPU][i];
	register uint32 offset = 0;
	for (int32_t id = 1; id < devCount; id++) {
		offset += (uint32)mGPU_elected_vars[id - 1].size();
		for (uint32 i = 0; i < mGPU_elected_vars[id].size(); i++) pv->PVs[offset + i] = mGPU_elected_vars[id][i];
	}
	// transfer the elected vars across the slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].pv->nParVars = slaves[id].pv->nVarsPerGPU = (uint32)mGPU_elected_vars[id].size();
		for (uint32 i = 0; i < slaves[id].pv->nParVars; i++) slaves[id].pv->PVs[i] = mGPU_elected_vars[id][i];
	}
	/* Calculating org. clauses upper bound on master GPU */
	calc_up(slaves, cnf, ot, pv, mGPUStats, devCount, devProps, streams, 1);
	// Verify results
	assert(up_validity(cnf, slaves, ot, pv, devCount, devProps, 1));
	assert(cnf_stats.max_org_cls > 0 && cnf_stats.max_org_lits > 0);
	/* Allocate memory for CNF on slave GPUs */
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		// re-initialize GPU cnf_stats
		assert(slaves[id].cnf_stats.max_org_cls > 0 && slaves[id].cnf_stats.max_org_lits > 0);
		slaves[id].cnf_stats.n_org_cls = slaves[id].cnf_stats.max_org_cls;
		slaves[id].cnf_stats.n_org_lits = slaves[id].cnf_stats.max_org_lits;
		slaves[id].GPUStats->org_nClauses = slaves[id].cnf_stats.max_org_cls;
		slaves[id].GPUStats->org_nLits = slaves[id].cnf_stats.max_org_lits;
		slaves[id].cnf_stats.max_added_cls = 0;
		slaves[id].cnf_stats.max_added_cl_width = 0;
		// allocation based on #clauses & #literals calculated by upper-bound kernel
		if (!mGPU_CNF_realloc_slave(id)) return false;
	}
	/* broadcast CNF to all slaves */
	if (verbose >= 2) cout << "c | Broadcasting CNF to all slaves...";
	broadcast_cnf(slaves, cnf, ot, devCount, devProps, timer, 1);
	if (verbose >= 2) cout << " ==> done." << endl;
	/* Reallocate memory for CNF on master GPU */
	cudaSetDevice(MASTER_GPU);
	cnf_stats.max_added_cl_width = 0; 
	cnf_stats.max_added_cls = 0;
	if (!mGPU_CNF_realloc_master(false)) return false;
	// create occurrence tables for all GPUs
	create_ot(slaves, cnf, ot, devCount, devProps, timer);
	// free system memory
	mGPU_elected_vars.clear(true);
	return true;
}
//TODO
void GEARS::mGPU_CNF_merge()
{
	if (verbose >= 2) cout << "c | Merging all CNFs to the master GPU..";
	// count the remaining of clauses & literals in all CNFs
	cnt_all(slaves, cnf, devCount, GPUStats, devProps);
	uint32 numClsAfter = acc_cls_after();
	int64 numLitsAfter = acc_lits_after();
	/* Check memory availability */
	double cnf_cap = sizeof(CNF) + (double)(numClsAfter * sizeof(CLAUSE) + numLitsAfter * sizeof(uint32)); // new master CNF
	if (verbose >= 2) cout << "(Cap.: " << cnf_cap / MBYTE << " MB)";
	// free old memory 
	double cnf_cap_old = sizeof(CNF) + (double)(cnf->size() * sizeof(CLAUSE) + cnf->nLits() * sizeof(uint32)); // old cnf
	gMemCons -= cnf_cap_old;
	if (gMemCons < 0) gMemCons = 0;
	gMemCons += cnf_cap;
	if (gMemCons > gMemMax) {
		if (verbose >= 2) cout << "\nc | Not enough memory space for CNF (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory space for CNF (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
	// allocate temp. memory for master CNF 
	cudaSetDevice(MASTER_GPU);
	size_t mem_sz = sizeof(CNF) + size_t(cnf_stats.n_cls_after) * sizeof(CLAUSE) + (size_t)cnf_stats.n_lits_after * sizeof(uint32);
	if (gMemCons + mem_sz > gMemMax) {
		if (verbose > 1) cout << "\nc | Not enough GPU memory for tmp CNF (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory for tmp CNF (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
	addr_t gMemCNF_tmp;
	CHECK(cudaMallocManaged((void**)&gMemCNF_tmp, mem_sz));
	addr_t bottom = gMemCNF_tmp + mem_sz;
	tmp_cnf = (CNF*)gMemCNF_tmp;
	gMemCNF_tmp += sizeof(CNF);
	tmp_cnf->set_nClauses(cnf_stats.n_cls_after);
	tmp_cnf->set_nLits(cnf_stats.n_lits_after);
	tmp_cnf->alloc_cls(&gMemCNF_tmp);
	tmp_cnf->alloc_lits(&gMemCNF_tmp);
	assert(gMemCNF_tmp == bottom);
	// copy master CNF to tmp. CNF
	tmp_cnf->reset_added();
	copy_cnf(cnf, tmp_cnf);
	assert(tmp_cnf->num_added_cls() == cnf_stats.n_cls_after);
	assert(tmp_cnf->num_added_lits() == cnf_stats.n_lits_after);
	// reallocate master CNF 
	mem_sz = sizeof(CNF) + (size_t)numClsAfter * sizeof(CLAUSE) + (size_t)numClsAfter * sizeof(uint32);
	if (verbose >= 2) cout << "c | New CNF mem. cap: " << mem_sz / MBYTE << " MB." << endl;
	assert(gMemCNF != NULL);
	CHECK(cudaFree(gMemCNF));
	gMemCNF = NULL;
	CHECK(cudaMallocManaged((void**)&gMemCNF, mem_sz));
	bottom = gMemCNF + mem_sz;
	cnf = (CNF*)gMemCNF;
	gMemCNF += sizeof(CNF);
	cnf->set_nClauses(numClsAfter);
	cnf->set_nLits(numClsAfter);
	cnf->alloc_cls(&gMemCNF);
	cnf->alloc_lits(&gMemCNF);
	assert(gMemCNF == bottom);
	// merge slave CNFs to the master
	merge_cnf(slaves, tmp_cnf, cnf, devCount, devProps, timer);
	cudaSetDevice(MASTER_GPU); // return context to master GPU
	assert(cnf->num_added_cls() == numClsAfter);
	assert(cnf->num_added_lits() == numLitsAfter);
	cnf->reset_added();
	// reinitialize original CNF info
	GPUStats->org_nClauses = cnf_stats.max_org_cls = numClsAfter;
	GPUStats->org_nLits = cnf_stats.max_org_lits = numLitsAfter;
	// free tmp. CNF
	CHECK(cudaFree(gMemCNF_tmp));
	if (verbose >= 2) cout << " ==> done." << endl;
	// reorder d_vars
	mGPU_VO_master(numLitsAfter);
	// realloc occur. table
	GPU_OT_alloc(0);
	if (verbose > 1 && !ve_plus_en) cout << "c |-----------------------------------------------------------------|" << endl;
}

bool GEARS::mGPU_CNF_realloc_master(bool initial)
{
	if (verbose >= 2) cout << "c | Reallocating managed memory for CNF on master GPU..";
	assert(pv->nVarsPerGPU > 0);
	register uint32 nVarsPerGPU = pv->nVarsPerGPU;
	/* Set added CNF information */
	cnf_stats.max_added_lits = cnf_stats.max_added_cl_width * cnf_stats.max_added_cls;
	/* Check memory availability */
	double cnf_cap = sizeof(CNF);
	cnf_cap += cnf_stats.max_org_cls * sizeof(CLAUSE) + cnf_stats.max_org_lits * sizeof(uint32); // original CNF
	cnf_cap += cnf_stats.max_added_cls * sizeof(CLAUSE) + cnf_stats.max_added_lits * sizeof(uint32); // added CNF
	if (verbose >= 2) cout << "(Cap.: " << cnf_cap / MBYTE << " MB)";
	// free old memory 
	double cnf_cap_old = sizeof(CNF) + (double)(cnf->size() * sizeof(CLAUSE) + cnf->nLits() * sizeof(uint32)); // old cnf
	gMemCons -= cnf_cap_old;
	if (gMemCons < 0) gMemCons = 0;
	gMemCons += cnf_cap;
	// leave some space for variable re-ordering if (ve+) is enabled 
	// variable re-ordering is done on what's left from ve, so the max. space won't exceed the original CNF
	double cap_org_lits = 0;
	if (ve_plus_en) cap_org_lits = (double)cnf_stats.max_org_lits * sizeof(uint32);
	bool mem_opt = false;
	if (gMemCons > (gMemMax - cap_org_lits)) {
		if (verbose >= 2) cout << "\nc | Not enough memory space for CNF on master GPU (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory space for CNF on master GPU (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		cout << "c | Trying to optimize memory consumption...";
		while (gMemCons > gMemMax && cnf_stats.max_added_cl_width >= MIN_ADDED_CL_LEN) {
			gMemCons -= cnf_stats.max_added_lits * sizeof(uint32);
			if (gMemCons < 0) gMemCons = 0;
			cnf_stats.max_added_cl_width >>= 1; // reduce clause width by half
			cnf_stats.max_added_lits = cnf_stats.max_added_cl_width * cnf_stats.max_added_cls;
			gMemCons += cnf_stats.max_added_lits * sizeof(uint32);
		}
		if (cnf_stats.max_added_cl_width < MIN_ADDED_CL_LEN) {
			cout << "(clause width " << cnf_stats.max_added_cl_width << " not sufficient) ==> failed." << endl;
			gMemCons -= cnf_cap;
			if (gMemCons < 0) gMemCons = 0;
			return false; 
		}
		mem_opt = true;
	}
	if (initial) {
		// verify clause index in CPU-side memory the same as GPU-side occurrence table
		// store the CNF info of non-elected variables on CPU side before recreation 
		unElected_CNF = new CNF_UE; 
		unElected_CNF->size = nOrgClauses() - cnf_stats.max_org_cls;
		unElected_CNF->nLits = nOrgLits() - cnf_stats.max_org_lits;
		for (int32_t id = 1; id < devCount; id++) {
			unElected_CNF->size -= slaves[id].cnf_stats.n_org_cls;
			unElected_CNF->nLits -= slaves[id].cnf_stats.n_org_lits;
		}
		unElected_CNF->clauses = new CLAUSE_UE[unElected_CNF->size];
		register uint32 unElectedCls = 0;
		register int64 unElectedLits = 0;
		for (uint32 c = 0; c < nOrgClauses(); c++) {
			if (!pv->isElectedCl[c]) {
				assert(unElectedCls <= unElected_CNF->size);
				unElected_CNF->clauses[unElectedCls].occur = c;
				unElected_CNF->clauses[unElectedCls].size = orgs[c]->size();
				unElectedLits += unElected_CNF->clauses[unElectedCls].size;
				unElectedCls++;
			}
		}
		assert(unElectedCls == unElected_CNF->size && unElectedLits == unElected_CNF->nLits);
		// reallocate org CNF 
		register uint32 tot_nClauses = cnf_stats.max_org_cls + unElected_CNF->size;
		register int64 tot_nLits = cnf_stats.max_org_lits + unElected_CNF->nLits;
		cnf->set_nClauses(tot_nClauses + cnf_stats.max_added_cls);
		cnf->set_nLits(tot_nLits + cnf_stats.max_added_lits);
		/*cnf->realloc_info();
		cnf->realloc_data();*/
		//cnf->reset(); 
		// copy elected part to master GPU
		register uint32 nCls = 0;
		register int64 nLits = 0;
		register uint32 x, poss_size, negs_size;
		
		for (uint32 v = 0; v < nVarsPerGPU; v++) {
			x = pv->PVs[v];
			uint32 p = V2D(x + 1), n = NEG(p);
			poss_size = ot->list(p)->size();
			negs_size = ot->list(n)->size();
			for (uint32 c = 0; c < poss_size; c++) {
				register uint32 occur = ot->list(p)->occur(c);
				B_REF h_cl = orgs[occur];
				CLAUSE* d_cl = cnf->clause(nCls++);
				d_cl->set_cl_ptr(cnf->d_ptr(nLits));
				d_cl->copyLitsFrom(h_cl->cl_ptr(), h_cl->size());
				d_cl->set_status(ORIGINAL);
				nLits += h_cl->size();
			}
			for (uint32 c = 0; c < negs_size; c++) {
				register uint32 occur = ot->list(n)->occur(c);
				B_REF h_cl = orgs[occur];
				CLAUSE* d_cl = cnf->clause(nCls++);
				d_cl->set_cl_ptr(cnf->d_ptr(nLits));
				d_cl->copyLitsFrom(h_cl->cl_ptr(), h_cl->size());
				d_cl->set_status(ORIGINAL);
				nLits += h_cl->size();
			}
		}
		assert(nCls == cnf_stats.max_org_cls && nLits == cnf_stats.max_org_lits);
		// restore the non-elected part
		for (uint32 c = 0; c < unElected_CNF->size; c++) {
			register uint32 occur = unElected_CNF->clauses[c].occur;
			assert((uint32)orgs[occur]->size() == unElected_CNF->clauses[c].size);
			B_REF h_cl = orgs[occur];
			CLAUSE* d_cl = cnf->clause(nCls++);
			d_cl->set_cl_ptr(cnf->d_ptr(nLits));
			d_cl->copyLitsFrom(h_cl->cl_ptr(), h_cl->size());
			d_cl->set_status(ORIGINAL);
			nLits += h_cl->size();
		}
		delete[] unElected_CNF;
		assert(nCls == tot_nClauses && nLits == tot_nLits);
		
		
		// calculate clauses signatures
		calc_clause_sig(cnf, 0, tot_nClauses, devProps);
		// re-initialize original cnf_stats
		GPUStats->org_nClauses = cnf_stats.max_org_cls = tot_nClauses;
		GPUStats->org_nLits = cnf_stats.max_org_lits = tot_nLits;
	}
	else {
		register uint32 unElected_nCls = cnf->size() - cnf_stats.max_org_cls;
		register int64 unElected_nLits = cnf->nLits() - cnf_stats.max_org_lits;
		for (int32_t id = 1; id < devCount; id++) {
			unElected_nCls -= slaves[id].cnf_stats.n_org_cls;
			unElected_nLits -= slaves[id].cnf_stats.n_org_lits;
		}
		register uint32 tot_nClauses = cnf_stats.max_org_cls + unElected_nCls;
		register int64 tot_nLits = cnf_stats.max_org_lits + unElected_nLits;
		// allocate tmp. CNF
		CHECK(cudaMallocManaged((void **)&tmp_cnf, sizeof(CNF)));
		tmp_cnf->set_nClauses(tot_nClauses);
		tmp_cnf->set_nLits(tot_nLits);
		/*tmp_cnf->alloc_cls();
		tmp_cnf->alloc_lits();*/
		// copy master CNF to tmp. CNF (elected + unelected)
		tmp_cnf->reset_added();
		copy_elect_cnf(cnf, tmp_cnf, ot, pv, devProps);
		assert(tmp_cnf->num_added_cls() == cnf_stats.max_org_cls);
		assert(tmp_cnf->num_added_lits() == cnf_stats.max_org_lits);
		copy_unelect_cnf(cnf, tmp_cnf, pv, devProps);
		assert(tmp_cnf->num_added_cls() == tot_nClauses);
		assert(tmp_cnf->num_added_lits() == tot_nLits);
		// reallocate master CNF
		cnf->set_nClauses(tot_nClauses + cnf_stats.max_added_cls);
		cnf->set_nLits(tot_nLits + cnf_stats.max_added_lits);
		/*cnf->realloc_info();
		cnf->realloc_data();*/
		//cnf->reset(); 
		// restore the org. CNF
		copy_cnf(tmp_cnf, cnf);
		assert(cnf->num_added_cls() == tot_nClauses);
		assert(cnf->num_added_lits() == tot_nLits);
		cnf->reset_added();
		// free tmp. CNF
		tmp_cnf->~CNF();
		CHECK(cudaFree(tmp_cnf));
		// reinitialize original CNF info
		GPUStats->org_nClauses = cnf_stats.max_org_cls = tot_nClauses;
		GPUStats->org_nLits = cnf_stats.max_org_lits = tot_nLits;
	}
	if (mem_opt) cout << "(Consumed: " << gMemCons / MBYTE << ") ==> successful." << endl;
	else if (verbose >= 2) cout << " ==> done." << endl;
	return true;
}

bool GEARS::mGPU_CNF_alloc_slave(int32_t id)
{
	if (verbose >= 2) cout << "c | Allocating managed memory for CNF on GPU " << id << "..";
	/* Set added CNF information */
	assert(slaves[id].cnf_stats.n_org_cls > 0 && slaves[id].cnf_stats.n_org_cls == slaves[id].cnf_stats.max_org_cls);
	assert(slaves[id].cnf_stats.n_org_lits > 0 && slaves[id].cnf_stats.n_org_lits == slaves[id].cnf_stats.max_org_lits);
	slaves[id].cnf_stats.max_added_lits = slaves[id].cnf_stats.max_added_cl_width * slaves[id].cnf_stats.max_added_cls;
	/* Check memory availability */
	double cnf_cap = sizeof(CNF);
	cnf_cap += slaves[id].cnf_stats.n_org_cls * sizeof(CLAUSE) + slaves[id].cnf_stats.n_org_lits * sizeof(uint32); // original CNF
	cnf_cap += slaves[id].cnf_stats.max_added_cls * sizeof(CLAUSE) + slaves[id].cnf_stats.max_added_lits * sizeof(uint32); // added CNF
	if (verbose >= 2) cout << "(Cap.: " << cnf_cap / MBYTE << " MB)";
	slaves[id].gMemCons += cnf_cap;
	double cap_org_lits = 0;
	if (ve_plus_en) cap_org_lits = (double)slaves[id].cnf_stats.n_org_lits * sizeof(uint32);
	bool mem_opt = false;
	if (slaves[id].gMemCons > (gMemMax - cap_org_lits)) {
		if (verbose >= 2) cout << "\nc | Not enough memory space for CNF on GPU " << id << " (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory space for CNF (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
		cout << "c | Trying to optimize memory consumption...";
		while (slaves[id].gMemCons > gMemMax && slaves[id].cnf_stats.max_added_cl_width >= MIN_ADDED_CL_LEN) {
			slaves[id].gMemCons -= slaves[id].cnf_stats.max_added_lits * sizeof(uint32);
			if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
			slaves[id].cnf_stats.max_added_cl_width >>= 1; // reduce clause width by half
			slaves[id].cnf_stats.max_added_lits = slaves[id].cnf_stats.max_added_cl_width * slaves[id].cnf_stats.max_added_cls;
			slaves[id].gMemCons += slaves[id].cnf_stats.max_added_lits * sizeof(uint32);
		}
		if (slaves[id].cnf_stats.max_added_cl_width < MIN_ADDED_CL_LEN) {
			cout << "(clause width " << slaves[id].cnf_stats.max_added_cl_width << " not sufficient) ==> failed." << endl;
			slaves[id].gMemCons -= cnf_cap;
			if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
			return false;
		}
		mem_opt = true;
	}
	CHECK(cudaMallocManaged((void **)&slaves[id].cnf, sizeof(CNF)));
	slaves[id].cnf->set_nClauses(slaves[id].cnf_stats.n_org_cls + slaves[id].cnf_stats.max_added_cls);
	//slaves[id].cnf->alloc_cls();
	slaves[id].cnf->set_nLits(slaves[id].cnf_stats.n_org_lits + slaves[id].cnf_stats.max_added_lits);
	//slaves[id].cnf->alloc_lits();
	slaves[id].cnf->reset_added();
	if (mem_opt) cout << "(Consumed: " << slaves[id].gMemCons / MBYTE << ") ==> successful." << endl;
	else if (verbose >= 2) cout << " ==> done." << endl;
	return true;
}

bool GEARS::mGPU_CNF_realloc_slave(int32_t id)
{
	if (verbose >= 2) cout << "c | Reallocating managed memory for CNF on GPU " << id << "..";
	/* Set added CNF information */
	assert(slaves[id].cnf_stats.n_org_cls > 0 && slaves[id].cnf_stats.n_org_cls == slaves[id].cnf_stats.max_org_cls);
	assert(slaves[id].cnf_stats.n_org_lits > 0 && slaves[id].cnf_stats.n_org_lits == slaves[id].cnf_stats.max_org_lits);
	slaves[id].cnf_stats.max_added_lits = slaves[id].cnf_stats.max_added_cl_width * slaves[id].cnf_stats.max_added_cls;
	/* Check memory availability */
	double cnf_cap = sizeof(CNF);
	cnf_cap += slaves[id].cnf_stats.n_org_cls * sizeof(CLAUSE) + slaves[id].cnf_stats.n_org_lits * sizeof(uint32); // original CNF
	cnf_cap += slaves[id].cnf_stats.max_added_cls * sizeof(CLAUSE) + slaves[id].cnf_stats.max_added_lits * sizeof(uint32); // added CNF
	if (verbose >= 2) cout << "(Cap.: " << cnf_cap / MBYTE << " MB)";
	// free old memory 
	double cnf_cap_old = sizeof(CNF) + (double)(slaves[id].cnf->size() * sizeof(CLAUSE) + slaves[id].cnf->nLits() * sizeof(uint32)); // old CNF
	slaves[id].gMemCons -= cnf_cap_old;
	if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
	slaves[id].gMemCons += cnf_cap;
	double cap_org_lits = 0;
	if (ve_plus_en) cap_org_lits = (double)slaves[id].cnf_stats.n_org_lits * sizeof(uint32);
	bool mem_opt = false;
	if (slaves[id].gMemCons > (gMemMax - cap_org_lits)) {
		if (verbose >= 2) cout << "\nc | Not enough memory space for CNF on GPU " << id << " (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory space for CNF on GPU " << id << " (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
		cout << "c | Trying to optimize memory consumption...";
		while (slaves[id].gMemCons > gMemMax && slaves[id].cnf_stats.max_added_cl_width >= MIN_ADDED_CL_LEN) {
			slaves[id].gMemCons -= slaves[id].cnf_stats.max_added_lits * sizeof(uint32);
			if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
			slaves[id].cnf_stats.max_added_cl_width >>= 1; // reduce clause width by half
			slaves[id].cnf_stats.max_added_lits = slaves[id].cnf_stats.max_added_cl_width * slaves[id].cnf_stats.max_added_cls;
			slaves[id].gMemCons += slaves[id].cnf_stats.max_added_lits * sizeof(uint32);
		}
		if (slaves[id].cnf_stats.max_added_cl_width < MIN_ADDED_CL_LEN) {
			cout << "(clause width " << slaves[id].cnf_stats.max_added_cl_width << " not sufficient) ==> failed." << endl;
			slaves[id].gMemCons -= cnf_cap;
			if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
			return false; 
		}
		mem_opt = true;
	}
	slaves[id].cnf->set_nClauses(slaves[id].cnf_stats.n_org_cls + slaves[id].cnf_stats.max_added_cls);
	slaves[id].cnf->set_nLits(slaves[id].cnf_stats.n_org_lits + slaves[id].cnf_stats.max_added_lits);
	/*slaves[id].cnf->realloc_info();
	slaves[id].cnf->realloc_data();*/
	//slaves[id].cnf->reset(); // write zero! (achieves the first touch & memory coherency)
	if (mem_opt) cout << "(Consumed: " << slaves[id].gMemCons / MBYTE << ") ==> successful." << endl;
	else if (verbose >= 2) cout << " ==> done." << endl;
	return true;
}

void GEARS::mGPU_VO_master(int64 tot_nLits)
{
	if (verbose > 1) cout << "c | Reordering variables after merging..";
	double mem_cap = (double)tot_nLits * sizeof(uint32);
	if (verbose > 1) cout << "(Cap.: " << mem_cap / MBYTE << " MB)";
	gMemCons += mem_cap;
	if (gMemCons > gMemMax) {
		if (verbose >= 2) cout << "\nc | Not enough memory space for reordering (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		else cout << "c | Not enough memory space for reordering (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}
	// calculate the histogram on CNF literals
	GPU_hist(tot_nLits);
	// set variable scores and calculate their weights (occur(x) * occur(~x))
	
	set_var_scores(d_vars, thrust::raw_pointer_cast(histogram.data()), devProps);
	// order variables w.r.t their scores
	thrust::sort(d_vars->occur(0), d_vars->occur(nOrgVars()), GPU_VAR_LCV());
	
	
	if (verbose >= 3) d_vars->print();
	// free histogram memory (no longer needed)
	histogram.clear();
	histogram.shrink_to_fit();
	gMemCons -= mem_cap;
	if (gMemCons < 0) gMemCons = 0;
	if (verbose > 1) cout << " ==> done." << endl;
}

void GEARS::mGPU_VE()
{
	int64 global_lits_before = 0;
	int64 global_lits_after = 0;
	int64 global_lits_removed = 0;
	if (pv->nParVars == 0 && verbose > 1)
		cout << "c | No variables elected --> (try another method)." << endl;
	if (pv->nParVars > 0) {
		if (verbose > 1) cout << "c | Eliminating variables for first round...";
		ve(slaves, cnf, ot, pv, devCount, GPUStats, devProps, timer);
		if (verbose > 1) cout << " ==> done." << endl;
		// evaluate reductions
		eval_reds(slaves, cnf, pv->elim, devCount, GPUStats, devProps, streams);
		cnf_stats.global_n_del_vars = acc_vars_eliminated();
		cnf_stats.global_n_cls = acc_cls_after();
		cnf_stats.global_n_lits = acc_lits_after();
		if (verbose == 1) {
			cout << "c |\nc |                        BVE Reductions          " << endl;
			cout << "c |  Parallel variables = " << pv->nParVars << endl;
			cout << "c |  Variables after = " << nOrgVars() - cnf_stats.global_n_del_vars << " (-" << cnf_stats.global_n_del_vars << ")" << endl;
			cout << "c |  Clauses after = " << cnf_stats.global_n_cls << " (-" << nOrgClauses() - cnf_stats.global_n_cls << ")" << endl;
			if (cnf_stats.global_n_lits > nOrgLits()) cout << "c |  Literals after = " << cnf_stats.global_n_lits << " (+" << cnf_stats.global_n_lits - nOrgLits() << ")" << endl;
			else cout << "c |  Literals after = " << cnf_stats.global_n_lits << " (-" << nOrgLits() - cnf_stats.global_n_lits << ")" << endl;
			if (cnf_stats.global_n_lits == 0) cout << "c | SATISFIABLE" << endl;
			cout << "c |" << endl;
		}
		if (verbose > 1) {
			cout << "c |\nc |                        BVE Reductions          " << endl;
			cout << "c |  Master GPU:" << endl;
			cout << "c |  -----------" << endl;
			cout << "c |  Variables: before = " << pv->nVarsPerGPU << ", after = " << pv->nVarsPerGPU - cnf_stats.n_del_vars << " (-" << cnf_stats.n_del_vars << ")" << endl;
			cout << "c |  Clauses  : before = " << cnf_stats.max_org_cls << ", after = " << cnf_stats.n_cls_after << " (-" << cnf_stats.max_org_cls - cnf_stats.n_cls_after << ")" << endl;
			if (cnf_stats.n_lits_after < GPUStats->org_nLits)
				cout << "c |  Literals : before = " << cnf_stats.max_org_lits << ", after = " << cnf_stats.n_lits_after << " (-" << cnf_stats.max_org_lits - cnf_stats.n_lits_after << ")" << endl;
			else cout << "c |  Literals : before = " << cnf_stats.max_org_lits << ", after = " << cnf_stats.n_lits_after << " (+" << cnf_stats.n_lits_after - cnf_stats.max_org_lits << ")" << endl;
			if (cnf_stats.n_lits_after == 0) cout << "c | SATISFIABLE" << endl;
			cout << "c |" << endl;
			for (int32_t id = 1; id < devCount; id++) {
				cout << "c |  Slave GPU " << id << ":" << endl;
				cout << "c |  ------------" << endl;
				cout << "c |  Variables: before = " << slaves[id].pv->nParVars << ", after = " << slaves[id].pv->nParVars - slaves[id].cnf_stats.n_del_vars << " (-" << slaves[id].cnf_stats.n_del_vars << ")" << endl;
				cout << "c |  Clauses  : before = " << slaves[id].cnf_stats.max_org_cls << ", after = " << slaves[id].cnf_stats.n_cls_after << " (-" << slaves[id].cnf_stats.max_org_cls - slaves[id].cnf_stats.n_cls_after << ")" << endl;
				if (slaves[id].cnf_stats.n_lits_after < slaves[id].cnf_stats.max_org_lits)
					cout << "c |  Literals : before = " << slaves[id].cnf_stats.max_org_lits << ", after = " << slaves[id].cnf_stats.n_lits_after << " (-" << slaves[id].cnf_stats.max_org_lits - slaves[id].cnf_stats.n_lits_after << ")" << endl;
				else cout << "c |  Literals : before = " << slaves[id].cnf_stats.max_org_lits << ", after = " << slaves[id].cnf_stats.n_lits_after << " (+" << slaves[id].cnf_stats.n_lits_after - slaves[id].cnf_stats.max_org_lits << ")" << endl;
				if (slaves[id].cnf_stats.n_cls_after == 0) cout << "c | SATISFIABLE" << endl;
				cout << "c |" << endl;
			}
		}
		// extend BVE untill no more reductions can obtained
		if (ve_plus_en) {
			// per GPU reductions
			cnf_stats.n_lits_before = cnf_stats.n_lits_after;
			for (int32_t id = 1; id < devCount; id++) slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
			// accumulate all reductions
			global_lits_before = cnf_stats.global_n_lits;
			global_lits_removed = nOrgLits() - global_lits_before;
			register uint32 nElectedVars_s, nElectedVars_master, nElectedVars_slave;
			while (global_lits_removed > LIT_REM_THR) {
				// master
				nElectedVars_s = 0; nElectedVars_master = 0;
				CHECK(cudaSetDevice(MASTER_GPU));
				for (uint32 v = 0; v < pv->nVarsPerGPU; v++) {
					register uint32 x = pv->PVs[v];
					if (!pv->elim[x]) pv->PVs[nElectedVars_master++] = x;
				}
				nElectedVars_s += nElectedVars_master;
				pv->nVarsPerGPU = nElectedVars_master;
				// slaves
				for (int32_t id = 1; id < devCount; id++) {
					CHECK(cudaSetDevice(id));
					nElectedVars_slave = 0;
					for (uint32 v = 0; v < slaves[id].pv->nParVars; v++) {
						register uint32 x = slaves[id].pv->PVs[v];
						if (!slaves[id].pv->elim[x]) slaves[id].pv->PVs[nElectedVars_slave++] = x;
					}
					slaves[id].pv->nParVars = slaves[id].pv->nVarsPerGPU = nElectedVars_slave;
					nElectedVars_s += nElectedVars_slave;
				}
				pv->nParVars = nElectedVars_s;
				if (pv->nParVars == 0) {
					if (verbose > 1) cout << "c | Nothing left to eliminate --> terminate procedure." << endl;
					return;
				}
				// HSE 
				if (verbose > 1) cout << "c | HSEing non-elim variables...";
				hse(slaves, cnf, ot, pv, devCount, devProps, timer);
				// count remaining literals
				cnt_lits(slaves, cnf, devCount, GPUStats, devProps);
				global_lits_after = acc_lits_after();
				global_lits_removed = global_lits_before - global_lits_after;
				global_lits_before = global_lits_after;
				if (verbose == 1) cout << "(Literals removed: -" << global_lits_removed << ") ==> done." << endl;
				if (verbose > 1) {
					cout << "\nc |  Literals removed on master GPU: -" << cnf_stats.n_lits_before - cnf_stats.n_lits_after << " ==> done." << endl;
					cnf_stats.n_lits_before = cnf_stats.n_lits_after;
					for (int32_t id = 1; id < devCount; id++) {
						cout << "c |  Literals removed on slave GPU " << id << ": -" << slaves[id].cnf_stats.n_lits_before - slaves[id].cnf_stats.n_lits_after << " ==> done." << endl;
						slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
					}
				}
				if (global_lits_removed <= LIT_REM_THR) break;
				// reallocate occurrence table
				mGPU_var_freq_nround();
				mGPU_OT_realloc();
				// execute BVE again
				if (verbose > 1) cout << "c | Eliminating variables for another round...";
				// create occurrence table to consider added clauses from previous round
				create_ot(slaves, cnf, ot, devCount, devProps, timer);
				ve(slaves, cnf, ot, pv, devCount, GPUStats, devProps, timer);
				// count remaining literals
				cnt_lits(slaves, cnf, devCount, GPUStats, devProps);
				global_lits_after = acc_lits_after();
				global_lits_removed = (global_lits_before > global_lits_after) ? global_lits_before - global_lits_after : global_lits_after - global_lits_before;
				if (verbose == 1) {
					if (global_lits_before > global_lits_after) cout << "(Literals reduction: -" << global_lits_removed << ") ==> done." << endl;
					else cout << "(Literals reduction: +" << global_lits_removed << ") ==> done." << endl;
				}
				global_lits_before = global_lits_after;
				if (verbose > 1) {
					cout << "\nc |  Literals reduction on master GPU: ";
					if (cnf_stats.n_lits_before > cnf_stats.n_lits_after) cout << "-" << cnf_stats.n_lits_before - cnf_stats.n_lits_after << " ==> done." << endl;
					else cout << "+" << cnf_stats.n_lits_after - cnf_stats.n_lits_before << " ==> done." << endl;
					cnf_stats.n_lits_before = cnf_stats.n_lits_after;
					for (int32_t id = 1; id < devCount; id++) {
						cout << "c |  Literals reduction on slave GPU " << id << ": ";
						if (slaves[id].cnf_stats.n_lits_before > slaves[id].cnf_stats.n_lits_after) cout << "-" << slaves[id].cnf_stats.n_lits_before - slaves[id].cnf_stats.n_lits_after << " ==> done." << endl;
						else cout << "-" << slaves[id].cnf_stats.n_lits_after - slaves[id].cnf_stats.n_lits_before << " ==> done." << endl;
						slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
					}
				}
			}
		}
	}
}

void GEARS::mGPU_VE_ext()
{
	int64 global_lits_before = 0;
	int64 global_lits_after = 0;
	int64 global_lits_removed = 0;
	if (pv->nParVars == 0) {
		cout << "c | No variables elected --> (try another method)." << endl;
	}
	else {
		if (verbose > 1) cout << "c | Eliminating variables for first round...";
		ve(slaves, cnf, ot, pv, devCount, GPUStats, devProps, timer);
		if (verbose > 1) cout << " ==> done." << endl;
		// evaluate reductions
		eval_reds(slaves, cnf, pv->elim, devCount, GPUStats, devProps, streams);
		uint32 numVarsElim = acc_vars_eliminated();
		uint32 numClsAfter = acc_cls_after();
		int64 numLitsAfter = acc_lits_after();
		if (verbose > 1) {
			cout << "c |\nc |                        BVE Reductions (all GPUs)         " << endl;
			cout << "c |  Parallel variables = " << pv->nParVars << endl;
			cout << "c |  Variables after = " << nOrgVars() - numVarsElim << " (-" << abs((int)(cnf_stats.global_n_del_vars - numVarsElim)) << ")" << endl;
			cout << "c |  Clauses after = " << numClsAfter << " (-" << cnf_stats.global_n_cls - numClsAfter << ")" << endl;
			if (numLitsAfter > cnf_stats.global_n_lits)
				cout << "c |  Literals after = " << numLitsAfter << " (+" << numLitsAfter - cnf_stats.global_n_lits << ")" << endl;
			else cout << "c |  Literals after = " << numLitsAfter << " (-" << cnf_stats.global_n_lits - numLitsAfter << ")" << endl;
			if (numLitsAfter == 0) cout << "c | SATISFIABLE" << endl;
			cout << "c |" << endl;
		}
		// update global counters
		cnf_stats.global_n_del_vars = numVarsElim;
		cnf_stats.global_n_cls = numClsAfter;
		cnf_stats.global_n_lits = numLitsAfter;
		// extend BVE untill no more reductions can obtained
		if (ve_plus_en) {
			// per GPU reductions
			cnf_stats.n_lits_before = cnf_stats.n_lits_after;
			for (int32_t id = 1; id < devCount; id++) slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
			// accumulate all reductions
			global_lits_before = cnf_stats.global_n_lits;
			global_lits_removed = nOrgLits() - global_lits_before;
			register uint32 nElectedVars_s, nElectedVars_master, nElectedVars_slave;
			while (global_lits_removed > LIT_REM_THR) {
				// master
				nElectedVars_s = 0; nElectedVars_master = 0;
				CHECK(cudaSetDevice(MASTER_GPU));
				for (uint32 v = 0; v < pv->nVarsPerGPU; v++) {
					register uint32 x = pv->PVs[v];
					if (!pv->elim[x]) pv->PVs[nElectedVars_master++] = x;
				}
				nElectedVars_s += nElectedVars_master;
				pv->nVarsPerGPU = nElectedVars_master;
				// slaves
				for (int32_t id = 1; id < devCount; id++) {
					CHECK(cudaSetDevice(id));
					nElectedVars_slave = 0;
					for (uint32 v = 0; v < slaves[id].pv->nParVars; v++) {
						register uint32 x = slaves[id].pv->PVs[v];
						if (!slaves[id].pv->elim[x]) slaves[id].pv->PVs[nElectedVars_slave++] = x;
					}
					slaves[id].pv->nParVars = slaves[id].pv->nVarsPerGPU = nElectedVars_slave;
					nElectedVars_s += nElectedVars_slave;
				}
				pv->nParVars = nElectedVars_s;
				if (pv->nParVars == 0) {
					if (verbose > 1) cout << "c | Nothing left to eliminate --> terminate procedure." << endl;
					return;
				}
				// HSE 
				if (verbose > 1) cout << "c | HSEing non-elim variables...";
				hse(slaves, cnf, ot, pv, devCount, devProps, timer);
				// count remaining literals
				cnt_lits(slaves, cnf, devCount, GPUStats, devProps);
				global_lits_after = acc_lits_after();
				global_lits_removed = global_lits_before - global_lits_after;
				global_lits_before = global_lits_after;
				if (verbose == 1) cout << "(Literals removed: -" << global_lits_removed << ") ==> done." << endl;
				if (verbose > 1) {
					cout << "\nc |  Literals removed on master GPU: -" << cnf_stats.n_lits_before - cnf_stats.n_lits_after << " ==> done." << endl;
					cnf_stats.n_lits_before = cnf_stats.n_lits_after;
					for (int32_t id = 1; id < devCount; id++) {
						cout << "c |  Literals removed on slave GPU " << id << ": -" << slaves[id].cnf_stats.n_lits_before - slaves[id].cnf_stats.n_lits_after << " ==> done." << endl;
						slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
					}
				}
				if (global_lits_removed <= LIT_REM_THR) break;
				// reallocate occurrence table
				mGPU_var_freq_nround();
				mGPU_OT_realloc();
				// execute BVE again
				if (verbose > 1) cout << "c | Eliminating variables for another round...";
				// create occurrence table to consider added clauses from previous round
				create_ot(slaves, cnf, ot, devCount, devProps, timer);
				ve(slaves, cnf, ot, pv, devCount, GPUStats, devProps, timer);
				// count remaining literals
				cnt_lits(slaves, cnf, devCount, GPUStats, devProps);
				global_lits_after = acc_lits_after();
				global_lits_removed = (global_lits_before > global_lits_after) ? global_lits_before - global_lits_after : global_lits_after - global_lits_before;
				if (verbose == 1) {
					if (global_lits_before > global_lits_after) cout << "(Literals reduction: -" << global_lits_removed << ") ==> done." << endl;
					else cout << "(Literals reduction: +" << global_lits_removed << ") ==> done." << endl;
				}
				global_lits_before = global_lits_after;
				if (verbose > 1) {
					cout << "\nc |  Literals reduction on master GPU: ";
					if (cnf_stats.n_lits_before > cnf_stats.n_lits_after) cout << "-" << cnf_stats.n_lits_before - cnf_stats.n_lits_after << " ==> done." << endl;
					else cout << "+" << cnf_stats.n_lits_after - cnf_stats.n_lits_before << " ==> done." << endl;
					cnf_stats.n_lits_before = cnf_stats.n_lits_after;
					for (int32_t id = 1; id < devCount; id++) {
						cout << "c |  Literals reduction on slave GPU " << id << ": ";
						if (slaves[id].cnf_stats.n_lits_before > slaves[id].cnf_stats.n_lits_after) cout << "-" << slaves[id].cnf_stats.n_lits_before - slaves[id].cnf_stats.n_lits_after << " ==> done." << endl;
						else cout << "-" << slaves[id].cnf_stats.n_lits_after - slaves[id].cnf_stats.n_lits_before << " ==> done." << endl;
						slaves[id].cnf_stats.n_lits_before = slaves[id].cnf_stats.n_lits_after;
					}
				}
			}
		}
	}
}

void GEARS::mGPU_var_freq_nround()
{
	// master
	CHECK(cudaSetDevice(MASTER_GPU));
	assert(cnf_stats.n_lits_after > 0);
	double mem_cap = (double)cnf_stats.n_lits_after * sizeof(uint32);
	gMemCons += mem_cap;
	if (gMemCons > gMemMax) {
		cout << "c | Not enough memory space for score-reassigning on master GPU (Allowed: " << gMemMax / MBYTE << ", Consumed: " << gMemCons / MBYTE << " MB)" << endl;
		// fail-safe switch
		cout << "c | Fail-safe switch is triggered." << endl;
		mGPU_CNF_merge();
		// write last simplified formula
		write_simp_cnf();
		print_reports();
		cudaDeviceReset();
		exit(EXIT_SUCCESS);
	}
	simplified.resize(cnf_stats.n_lits_after);
	*nLits_after = 0;
	// slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		double mem_cap_slave = (double)slaves[id].cnf_stats.n_lits_after * sizeof(uint32);
		slaves[id].gMemCons += mem_cap_slave;
		if (slaves[id].gMemCons > gMemMax) {
			cout << "c | Not enough memory space for score-reassigning (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
			// fail-safe switch
			cout << "c | Fail-safe switch is triggered." << endl;
			mGPU_CNF_merge();
			// write last simplified formula
			write_simp_cnf();
			print_reports();
			cudaDeviceReset();
			exit(EXIT_SUCCESS);
		}
		slaves[id].simplified.resize(slaves[id].cnf_stats.n_lits_after);
		*slaves[id].nLits_after = 0;
	}
	if (verbose > 1) cout << "c | Assigning new variables scores on all GPUs...";
	// filter CNFs from removed literals & copy to simplified (all GPUs)
	copy_if(slaves, cnf, thrust::raw_pointer_cast(simplified.data()), nLits_after, devCount, devProps);
	// calculate histogram for the simplified formulas (all GPUs)
	mGPU_hist_simp();
	// assign new scores to master GPU variables
	CHECK(cudaSetDevice(MASTER_GPU));
	v_score.clear();
	v_score.resize(nOrgVars());
	register uint32 v_idx, v_pos, v_neg;
	thrust::host_vector<uint32> h_histogram_m = histogram;
	for (uint32 v = 0; v < nOrgVars(); v++) {
		v_idx = v + 1;
		v_pos = V2D(v_idx);
		v_neg = NEG(v_pos);
		v_score[v].first = h_histogram_m[v_pos];
		v_score[v].second = h_histogram_m[v_neg];
		if (pv->elim[v]) assert(v_score[v].first == 0 && v_score[v].second == 0);
	}
	h_histogram_m.clear();
	h_histogram_m.shrink_to_fit();
	histogram.clear();
	histogram.shrink_to_fit();
	simplified.clear();
	simplified.shrink_to_fit();
	gMemCons -= mem_cap;
	if (gMemCons < 0) gMemCons = 0;
	// assign new scores to slave GPUs 
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		slaves[id].v_score.clear();
		slaves[id].v_score.resize(nOrgVars());
		thrust::host_vector<uint32> h_histogram_s = slaves[id].histogram;
		for (uint32 v = 0; v < nOrgVars(); v++) {
			v_idx = v + 1;
			v_pos = V2D(v_idx);
			v_neg = NEG(v_pos);
			slaves[id].v_score[v].first = h_histogram_s[v_pos];
			slaves[id].v_score[v].second = h_histogram_s[v_neg];
			if (slaves[id].pv->elim[v]) assert(slaves[id].v_score[v].first == 0 && slaves[id].v_score[v].second == 0);
		}
		slaves[id].histogram.clear();
		slaves[id].histogram.shrink_to_fit();
		slaves[id].simplified.clear();
		slaves[id].simplified.shrink_to_fit();
		slaves[id].gMemCons -= mem_cap;
		if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
	}
	if (verbose > 1) cout << " ==> done." << endl;
}

//TODO
void GEARS::mGPU_OT_realloc()
{
	//double vt_cap = ((double) sizeof(OT) + nOrgVars() * sizeof(OL)) * 2;
	//register uint32 pos_entries = 0, pscore;
	//register uint32 neg_entries = 0, nscore;
	//register uint32 pos_entries_idx = 0, neg_entries_idx = 0;
	//// master
	//if (verbose > 1) cout << "c | Reallocating occurrence table on master GPU..";
	//CHECK(cudaSetDevice(MASTER_GPU));
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	pscore = v_score[v].first;
	//	nscore = v_score[v].second;
	//	p_ot->list(v)->set_size(pscore);
	//	n_ot->list(v)->set_size(nscore);
	//	pos_entries += pscore;
	//	neg_entries += nscore;
	//}
	//// modify memory consumed by OT
	//gMemCons -= (double)p_ot->cap() * sizeof(uint32);
	//gMemCons -= (double)n_ot->cap() * sizeof(uint32);
	//gMemCons -= vt_cap;
	//if (gMemCons < 0) gMemCons = 0;
	//vt_cap += (pos_entries + neg_entries) * sizeof(uint32);
	//gMemCons += vt_cap;
	//if (verbose > 1) cout << "(Cap.: " << vt_cap / MBYTE << " MB)";
	//p_ot->set_nEntries(pos_entries);
	////p_ot->realloc_list_data();
	//n_ot->set_nEntries(neg_entries);
	////n_ot->realloc_list_data();
	//// assign list/tail pointers
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	uint32 *plist_ptr = p_ot->d_ptr(pos_entries_idx);
	//	uint32 *nlist_ptr = n_ot->d_ptr(neg_entries_idx);
	//	p_ot->list(v)->set_ptr(plist_ptr);
	//	n_ot->list(v)->set_ptr(nlist_ptr);
	//	pos_entries_idx += v_score[v].first;
	//	neg_entries_idx += v_score[v].second;
	//}
	//if (verbose > 1) cout << " ==> done." << endl;
	//// slaves
	//for (int32_t id = 1; id < devCount; id++) {
	//	CHECK(cudaSetDevice(id));
	//	if (verbose > 1) cout << "c | Reallocating occurrence table on slave GPU " << id << "..";
	//	vt_cap = ((double) sizeof(OT) + nOrgVars() * sizeof(OL)) * 2;
	//	pos_entries = 0; pscore = 0;
	//	neg_entries = 0; nscore = 0;
	//	for (uint32 v = 0; v < nOrgVars(); v++) {
	//		pscore = slaves[id].v_score[v].first;
	//		nscore = slaves[id].v_score[v].second;
	//		slaves[id].p_ot->list(v)->set_size(pscore);
	//		slaves[id].n_ot->list(v)->set_size(nscore);
	//		pos_entries += pscore;
	//		neg_entries += nscore;
	//	}
	//	// modify memory consumed by OT
	//	slaves[id].gMemCons -= (double)slaves[id].p_ot->cap() * sizeof(uint32);
	//	slaves[id].gMemCons -= (double)slaves[id].n_ot->cap() * sizeof(uint32);
	//	slaves[id].gMemCons -= vt_cap;
	//	if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
	//	vt_cap += (pos_entries + neg_entries) * sizeof(uint32);
	//	slaves[id].gMemCons += vt_cap;
	//	if (verbose > 1) cout << "(Cap.: " << vt_cap / MBYTE << " MB)";
	//	slaves[id].p_ot->set_nEntries(pos_entries);
	//	//slaves[id].p_ot->realloc_list_data();
	//	slaves[id].n_ot->set_nEntries(neg_entries);
	//	//slaves[id].n_ot->realloc_list_data();
	//	// assign list/tail pointers
	//	pos_entries_idx = 0; neg_entries_idx = 0;
	//	for (uint32 v = 0; v < nOrgVars(); v++) {
	//		uint32 *plist_ptr = slaves[id].p_ot->d_ptr(pos_entries_idx);
	//		uint32 *nlist_ptr = slaves[id].n_ot->d_ptr(neg_entries_idx);
	//		slaves[id].p_ot->list(v)->set_ptr(plist_ptr);
	//		slaves[id].n_ot->list(v)->set_ptr(nlist_ptr);
	//		pos_entries_idx += slaves[id].v_score[v].first;
	//		neg_entries_idx += slaves[id].v_score[v].second;
	//	}
	//	if (verbose > 1) cout << " ==> done." << endl;
	//}
}
//TODO
void GEARS::mGPU_OT_realloc_slave(int32_t id)
{
	//if (verbose > 1) cout << "c | Reallocating occurrence table on slave GPU " << id << "..";
	//double vt_cap = ((double) sizeof(OT) + nOrgVars() * sizeof(OL)) * 2;
	//register uint32 pos_entries = 0, pscore = 0;
	//register uint32 neg_entries = 0, nscore = 0;
	//register uint32 pos_entries_idx = 0, neg_entries_idx = 0;
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	pscore = slaves[id].v_score[v].first;
	//	nscore = slaves[id].v_score[v].second;
	//	slaves[id].p_ot->list(v)->set_size(pscore);
	//	slaves[id].n_ot->list(v)->set_size(nscore);
	//	pos_entries += pscore;
	//	neg_entries += nscore;
	//}
	//// modify memory consumed by OT
	//slaves[id].gMemCons -= (double)slaves[id].p_ot->cap() * sizeof(uint32);
	//slaves[id].gMemCons -= (double)slaves[id].n_ot->cap() * sizeof(uint32);
	//slaves[id].gMemCons -= vt_cap;
	//if (slaves[id].gMemCons < 0) slaves[id].gMemCons = 0;
	//vt_cap += (pos_entries + neg_entries) * sizeof(uint32);
	//slaves[id].gMemCons += vt_cap;
	//if (verbose > 1) cout << "(Cap.: " << vt_cap / MBYTE << " MB)";
	//slaves[id].p_ot->set_nEntries(pos_entries);
	////slaves[id].p_ot->realloc_list_data();
	//slaves[id].n_ot->set_nEntries(neg_entries);
	////slaves[id].n_ot->realloc_list_data();
	//// assign list/tail pointers
	//pos_entries_idx = 0; neg_entries_idx = 0;
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	uint32 *plist_ptr = slaves[id].p_ot->d_ptr(pos_entries_idx);
	//	uint32 *nlist_ptr = slaves[id].n_ot->d_ptr(neg_entries_idx);
	//	slaves[id].p_ot->list(v)->set_ptr(plist_ptr);
	//	slaves[id].n_ot->list(v)->set_ptr(nlist_ptr);
	//	pos_entries_idx += slaves[id].v_score[v].first;
	//	neg_entries_idx += slaves[id].v_score[v].second;
	//}
	//if (verbose > 1) cout << " ==> done." << endl;
}
//TODO
void GEARS::mGPU_OT_alloc(int32_t id)
{
	//if (verbose >= 2) cout << "c | Allocating managed memory for the occurrence table on GPU " << id << "..";
	//assert(nOrgVars() > 0);
	///* Check memory availability */
	//double vt_cap = ((double)sizeof(OT) + nOrgVars() * sizeof(OL)) * 2;
	///* allocate managed device memory for the occurrence list */
	//CHECK(cudaMallocManaged((void **)&slaves[id].p_ot, sizeof(OT)));
	//CHECK(cudaMallocManaged((void **)&slaves[id].n_ot, sizeof(OT)));
	//slaves[id].p_ot->set_size(nOrgVars());
	//slaves[id].n_ot->set_size(nOrgVars());
	////slaves[id].p_ot->alloc_lists();
	////slaves[id].n_ot->alloc_lists();
	//register uint32 pos_entries = 0, pscore;
	//register uint32 neg_entries = 0, nscore;
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	pscore = v_score[v].first;
	//	nscore = v_score[v].second;
	//	slaves[id].p_ot->list(v)->set_size(pscore);
	//	slaves[id].n_ot->list(v)->set_size(nscore);
	//	pos_entries += pscore;
	//	neg_entries += nscore;
	//}
	//vt_cap += (pos_entries + neg_entries) * sizeof(uint32);
	//if (verbose >= 2) cout << "(Cap.: " << vt_cap / MBYTE << " MB)";
	//slaves[id].gMemCons += vt_cap;
	//if (slaves[id].gMemCons > gMemMax) {
	//	if (verbose >= 2) cout << "\nc | Not enough memory in GPU " << id << "for OT (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
	//	else cout << "c | Not enough memory in GPU " << id << "for OT (Allowed: " << gMemMax / MBYTE << ", Consumed: " << slaves[id].gMemCons / MBYTE << " MB)" << endl;
	//	CHECK(cudaFree(slaves[id].p_ot->list(0)));
	//	CHECK(cudaFree(slaves[id].n_ot->list(0)));
	//	CHECK(cudaFree(slaves[id].p_ot));
	//	CHECK(cudaFree(slaves[id].n_ot));
	//	exit(EXIT_FAILURE);
	//}
	//slaves[id].p_ot->set_nEntries(pos_entries);
	////slaves[id].p_ot->alloc_entries();
	//slaves[id].n_ot->set_nEntries(neg_entries);
	////slaves[id].n_ot->alloc_entries();
	//// assign list/tail pointers
	//register uint32 pos_entries_idx = 0, neg_entries_idx = 0;
	//for (uint32 v = 0; v < nOrgVars(); v++) {
	//	uint32 *plist_ptr = slaves[id].p_ot->d_ptr(pos_entries_idx);
	//	uint32 *nlist_ptr = slaves[id].n_ot->d_ptr(neg_entries_idx);
	//	slaves[id].p_ot->list(v)->set_ptr(plist_ptr);
	//	slaves[id].n_ot->list(v)->set_ptr(nlist_ptr);
	//	pos_entries_idx += v_score[v].first;
	//	neg_entries_idx += v_score[v].second;
	//}
	//if (verbose >= 2) cout << " ==> done." << endl;
}

void GEARS::mGPU_hist_simp()
{
	// set max. number of histogram bins
	size_t num_bins = V2D(nOrgVars() + 1);
	// master
	CHECK(cudaSetDevice(MASTER_GPU));
	thrust::sort(simplified.begin(), simplified.end());
	histogram.resize(num_bins);
	thrust::counting_iterator<size_t> search_begin_m(0);
	thrust::upper_bound(simplified.begin(), simplified.end(), search_begin_m, search_begin_m + num_bins, histogram.begin());
	thrust::adjacent_difference(histogram.begin(), histogram.end(), histogram.begin());
	// slaves
	for (int32_t id = 1; id < devCount; id++) {
		CHECK(cudaSetDevice(id));
		thrust::sort(slaves[id].simplified.begin(), slaves[id].simplified.end());
		slaves[id].histogram.resize(num_bins);
		thrust::counting_iterator<size_t> search_begin_s(0);
		thrust::upper_bound(slaves[id].simplified.begin(), slaves[id].simplified.end(), search_begin_s, search_begin_s + num_bins, slaves[id].histogram.begin());
		thrust::adjacent_difference(slaves[id].histogram.begin(), slaves[id].histogram.end(), slaves[id].histogram.begin());
	}
}

int64 GEARS::acc_lits_after(void)
{
	// master
	register int64 numLitsAfter = cnf_stats.n_lits_after;
	// slaves
	for (int32_t id = 1; id < devCount; id++) numLitsAfter += slaves[id].cnf_stats.n_lits_after;
	return numLitsAfter;
}

uint32 GEARS::acc_cls_after(void)
{
	// master
	register uint32 numClsAfter = cnf_stats.n_cls_after;
	// slaves
	for (int32_t id = 1; id < devCount; id++) numClsAfter += slaves[id].cnf_stats.n_cls_after;
	return numClsAfter;
}

uint32 GEARS::acc_vars_eliminated(void)
{
	// master
	uint32 numVarsElim = cnf_stats.n_del_vars;
	// slaves
	for (int32_t id = 1; id < devCount; id++) numVarsElim += slaves[id].cnf_stats.n_del_vars;
	return numVarsElim;
}