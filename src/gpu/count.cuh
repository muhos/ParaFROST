/***********************************************************************[count.cuh]
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

#ifndef __GPU_COUNT_
#define __GPU_COUNT_

#include "definitions.cuh"
#include "definitions.hpp"

namespace ParaFROST {

	#define MAXREDUCEBLOCKS 1024

	extern 
	uint32 hostCBlocks[MAXREDUCEBLOCKS];
	extern
	uint32 hostLBlocks[MAXREDUCEBLOCKS];

	extern __device__
	uint32 gcounter;
	extern __device__
	uint32 devCBlocks[MAXREDUCEBLOCKS];
	extern __device__
	uint32 devLBlocks[MAXREDUCEBLOCKS];

	__global__ void reset_counter();

	__global__ void print_counter();

	__global__ void check_counter(const uint32 checksum);

	template <typename D, typename I>
	inline D seqreduceBlocks(const D* blocks, const I& n)
	{
		D finalcount = 0;
		for (I i = 0; i < n; ++i)
			finalcount += blocks[i];
		return finalcount;
	}

	template <typename D, typename I>
	inline void seqreduceBlocks(const D* CBlocks, const D* LBlocks, const I& n)
	{
		inf.numClausesSurvived = 0;
		inf.numLiteralsSurvived = 0;
		for (I i = 0; i < n; ++i) {
			inf.numClausesSurvived += CBlocks[i];
			inf.numLiteralsSurvived += LBlocks[i];
		}
	}

}

#endif