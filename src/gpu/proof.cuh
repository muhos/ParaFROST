/***********************************************************************[proof.cuh]
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

#ifndef __GPU_PROOF_
#define __GPU_PROOF_

#include "proof.hpp"
#include "memory.cuh"
#include "vector.cuh"

namespace ParaFROST {

	class cuMM;
	class PROOF;

	class cuPROOF {

		cuMM&	cumm;
		PROOF&	proof;

		cuTIMER	cutimer;
		cuPool  hostPool, devicePool;
		cuVecB  *hostStream, *deviceStream;
		cuVecB  header;
		size_t  deviceAdded, bytesWritten;
		float 	time;
		bool	enabled;

	public:

		cuPROOF(cuMM& _cumm, PROOF& _proof) :
			cumm(_cumm), proof(_proof)
			, cutimer()
			, hostPool()
			, devicePool()
			, hostStream(NULL)
			, deviceStream(NULL)
			, deviceAdded(0)
			, bytesWritten(0)
			, time(0.0f)
			, enabled(false)
		{}

		inline float 	gpuTime		() const { return time; }
		inline size_t	gpuClauses	() const { return deviceAdded; }
		inline size_t	gpuBytes	() const { return bytesWritten; }
		inline cuVecB*	gpuStream	() { return deviceStream; }	

		void	destroy			();
		void	writeClause		(addr_t&);
		void	cacheProof		(const cudaStream_t&);
		void	writeProof		(const cudaStream_t&);
		bool	alloc			(const uint32&);
		uint32	count			(const uint32*, const uint32&);
	};

} 

#endif