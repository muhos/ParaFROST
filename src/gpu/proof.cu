/***********************************************************************[proof.cu]
Copyright(c) 2021, Muhammad Osama - Anton Wijs,
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

#include "count.cuh"
#include "proof.cuh"
#include "timer.cuh"
#include "shared.cuh"
#include "reduce.cuh"
#include "options.cuh"
#include "primitives.cuh"

using namespace ParaFROST;

//=======================================================
// constant proof lookup table
#define MAXLEADINGZEROS 30
__constant__
Byte BLUT[MAXLEADINGZEROS + 1] =
{
	1, 1, 1, 1, 1, 1,
	2, 2, 2, 2, 2, 2, 2,
	3, 3, 3, 3, 3, 3, 3,
	4, 4, 4, 4, 4, 4, 4,
	5, 5, 5, 5
};
#define COUNTBYTES(LIT) BLUT[MAXLEADINGZEROS - __clz(LIT)]
//=======================================================

__global__ void printHead(cuVecB* proof)
{
	printf("this = %p, mem = %p, size = %d, cap = %d", 
		proof, proof->data(), proof->size(), proof->capacity());
}

_PFROST_D_ void countBytes(const uint32& lit, Byte& perLit, uint32& total)
{
	// here we rely on data racing to avoid
	// counting the bytes for a duplicate
	// assuming perLit is initially 0
	if (!perLit) {
		Byte local;
		ORIGINIZELIT(orgLit, lit);
		perLit = local = COUNTBYTES(orgLit);
		total += local;
	}
	else // the more threads taking this path the better
		total += perLit;
	assert(total);
	assert(perLit >= 1 && perLit <= 5);
}

__global__ 
void cnt_proof(const uint32* __restrict__ literals, const uint32 numLits)
{
	uint32* sh_bytes = SharedMemory<uint32>();
	grid_t tid = global_tx_off;
	uint32 nbytes = 0;
	while (tid < numLits) {
		addr_t lbyte = DC_PTRS->d_lbyte;
		uint32 lit = literals[tid];
		countBytes(lit, lbyte[lit], nbytes);
		grid_t off = tid + blockDim.x;
		if (off < numLits) {
			lit = literals[off];
			countBytes(lit, lbyte[lit], nbytes);
		}
		tid += stride_x_off;
	}
	loadShared(sh_bytes, nbytes, numLits);
	sharedReduce(sh_bytes, nbytes);
	warpReduce(sh_bytes, nbytes);
	if (!threadIdx.x) devLBlocks[blockIdx.x] = nbytes;
}

__global__ 
void cnt_proof_verify(const uint32* __restrict__ literals, const uint32 numLits)
{
	grid_t tid = 0;
	while (tid < numLits) {
		addr_t lbyte = DC_PTRS->d_lbyte;
		const uint32 lit = literals[tid];
		if (lit & 0xF0000000) lbyte[lit] = 5;
		else if (lit & 0x0FE00000) lbyte[lit] = 4;
		else if (lit & 0x001FC000) lbyte[lit] = 3;
		else if (lit & 0x00003F80) lbyte[lit] = 2;
		else lbyte[lit] = 1;
		printf(" literal(%d) has %d bytes of its original\n", SIGN(lit) ? -int(ABS(lit)) : ABS(lit), lbyte[lit]);
		gcounter += lbyte[lit];
		tid++;
	}
	printf(" total = %d\n", gcounter);
}

uint32 cuPROOF::count(const uint32* literals, const uint32& numLits)
{
	if (!proof.checkFile()) LOGERROR("host proof system is not activated");
	if (!literals) return 0;
	if (!numLits) return 0;
	enabled = true;
	LOGN2(2, " Counting proof bytes..");
	OPTIMIZEBLOCKS2(numLits, BLOCK1D);
	OPTIMIZESHARED(blockSize, sizeof(uint32));
	SYNCALL; // sync any pending kernels or transfers
	if (gopts.profile_gpu) cutimer->start();
	cnt_proof << <nBlocks, blockSize, smemSize >> > (literals, numLits);
	LASTERR("Proof counting failed");
	CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
	if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	const uint32 maxcap = seqreduceBlocks(hostLBlocks, nBlocks);
	assert(maxcap && maxcap < (numLits * sizeof(uint32)));
	LOGENDING(2, 5, "(%d bytes)", maxcap);
	return maxcap;
}

void cuPROOF::writeClause(addr_t& byte)
{
	assert(enabled);
	const size_t bytes = 16;
	const char delimiter = '0';
	char line[bytes];
	char* tail = line + bytes;
	*--tail = 0;
	while (*byte) {
		assert(byte != hostStream->end());
		uint32 ulit = 0, shift = 0;
		Byte b;
		do {
			b = *byte++;
			ulit |= (b & BYTEMASK) << shift;
			shift += 7;
		} while (b & BYTEMAX);
		const LIT_ST sign = SIGN(ulit);
		if (sign) proof.write('-');
		char* nstr = tail;
		assert(!*nstr);
		int digit = ABS(ulit);
		while (digit) {
			*--nstr = (digit % 10) + delimiter;
			digit /= 10;
		}
		while (nstr != tail) proof.write(*nstr++);
		proof.write(' ');
	}
	proof.write(delimiter);
	proof.write('\n');
}

void cuPROOF::writeProof(const cudaStream_t& _s)
{
	if (!enabled) return;
	if (hostStream->empty()) return;
	SYNC(_s);
	assert(**hostStream == PROOF_ADDED || **hostStream == PROOF_DELETED);
	assert(!hostStream->back());
	size_t prevlines = deviceAdded;
	LOGN2(2, " Writing GPU proof data..");
	Byte* byte = *hostStream;
	Byte* end = hostStream->end();
	if (proof.isNonBinary()) {
		while (byte != end) {
			assert(*byte == PROOF_ADDED || *byte == PROOF_DELETED);
			if (*byte++ == PROOF_DELETED) {
				proof.write('d');
				proof.write(' ');
			}
			writeClause(byte);
			deviceAdded++;
			byte++; // skip 0
		}
	}
	else {
		while (byte != end) {
			const Byte b = *byte++;
			proof.write(b);
			if (!b) deviceAdded++;
		}
	}
	size_t lines = deviceAdded - prevlines;
	LOGENDING(2, 5, "(%zd clauses, %d bytes)", lines, hostStream->size());
	bytesWritten += hostStream->size();
	hostStream->clear();
}

bool cuPROOF::alloc(const uint32& maxcap)
{
	if (!maxcap) return false;
	assert(enabled);
	assert(hostPool.cap == devicePool.cap);
	const size_t header_size = sizeof(cuVecB);
	const size_t proof_cap   = maxcap;
	const size_t min_cap     = header_size + proof_cap;
	assert(min_cap);
	if (devicePool.cap < min_cap) {
		// device memory
		cumm.DFREE(devicePool);
		assert(devicePool.mem == NULL);
		if (!cumm.hasDeviceMem(min_cap, "Proof")) return false;
		CHECK(cudaMalloc((void**)&devicePool.mem, min_cap));
		addr_t ea = devicePool.mem;
		deviceStream = (cuVecB*)ea, ea += header_size;
		header.alloc(ea, uint32(proof_cap)), ea += proof_cap;
		CHECK(cudaMemcpyAsync(deviceStream, &header, header_size, cudaMemcpyHostToDevice));
		assert(ea == devicePool.mem + min_cap);
		devicePool.cap = min_cap;
		// pinned host memory
		if (hostPool.mem) {
			CHECK(cudaFreeHost(hostPool.mem));
			hostPool.mem = NULL;
		}
		CHECK(cudaMallocHost((void**)&hostPool.mem, min_cap));
		ea = hostPool.mem;
		hostStream = (cuVecB*)ea, ea += header_size;
		hostStream->alloc(ea, uint32(proof_cap)), ea += proof_cap;
		assert(ea == hostPool.mem + min_cap);
		hostPool.cap = min_cap;
		SYNC(0); 
		header.clear(true);
	}
	return true;
}

void cuPROOF::cacheProof(const cudaStream_t& _s) 
{
	if (!enabled) return;
	assert(hostStream);
	assert(deviceStream);
	CHECK(cudaMemcpy(&header, deviceStream, sizeof(cuVecB), cudaMemcpyDeviceToHost)); 
	const uint32 devSize = header.size();
	if (!devSize) return;
	header.clear();
	CHECK(cudaMemcpyAsync(deviceStream, &header, sizeof(cuVecB), cudaMemcpyHostToDevice, _s)); // reset gpu proof size
	assert(header.data());
	assert(hostStream->capacity() == header.capacity());
	hostStream->resize(devSize);
	if (gopts.profile_gpu) cutimer->start(_s);
	CHECK(cudaMemcpyAsync(hostStream->data(), header.data(), devSize, cudaMemcpyDeviceToHost, _s));
	if (gopts.sync_always) SYNC(_s);
	if (gopts.profile_gpu) cutimer->stop(_s), cutimer->ve += cutimer->gpuTime();
}

void cuPROOF::destroy()
{
	LOGN2(2, " Freeing up proof host-device memory..");
	cumm.DFREE(devicePool);
	if (hostPool.mem) {
		CHECK(cudaFreeHost(hostPool.mem));
		hostPool.mem = NULL, hostPool.cap = 0;
	}
	LOGDONE(2, 5);
}