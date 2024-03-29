/***********************************************************************[proof.hpp]
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

#ifndef __PROOF_
#define __PROOF_

#include "vector.hpp"
#include "clause.hpp"
#include "logging.hpp"
#include "space.hpp"

namespace ParaFROST {

	class PROOF {

		FILE*	proofFile;
		SP*		sp;
		uint32*	vars;
		Lits_t	clause, tmpclause;
		size_t	added, bytesWritten;
		bool	nonbinary_en;

		#ifdef _WIN32
			#define writebyte putc
		#else
			#define writebyte putc_unlocked
		#endif

	protected:

		inline void	write			(const uint32*, const int&);
		inline void	binary			(const uint32*, const int&);
		inline void	nonbinary		(const uint32*, const int&);
		inline void	addline			(const uint32*, const int&);
		inline void	delline			(const uint32*, const int&);
		inline void	addClause		();
		inline void	deleteClause	();
	
	public:

		inline void		write		(const Byte& byte) { 
			bytesWritten++;
			writebyte(byte, proofFile); 
		}
		inline bool		isNonBinary	() const { return nonbinary_en; }
		inline size_t	clauses		() const { return added; }
		inline size_t	bytes		() const { return bytesWritten; }

				PROOF			();
				~PROOF			();
		bool	checkFile		();
		void	close			();
		int		exists			();
		int		exists			(uint32*, const int&);
		void	init			(SP*);
		void	init			(SP*, uint32*);
		void	handFile		(arg_t path, const bool&);
		void	checkInput		(arg_t input);
		void	addEmpty		();
		void	addUnit			(uint32);
		void	addClause		(Lits_t&);
		void	addClause		(CLAUSE&);
		void	deleteClause	(Lits_t&);
		void	deleteClause	(CLAUSE&);
		void	deleteClause	(uint32*, const int&);
		void	shrinkClause	(CLAUSE&);
		void	shrinkClause	(CLAUSE&, const uint32&);
		void	printClause		(const char*, const uint32*, const int&, const bool& map = false);

	};

}

#endif