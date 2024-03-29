/***********************************************************************[proof.h]
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

#ifndef __PROOF_
#define __PROOF_

#include "vector.h"
#include "definitions.h"
#include "clause.h"
#include "simptypes.h"
#include "logging.h"
#include "space.h"

namespace ParaFROST {

	class PROOF {

		FILE* proofFile;
		SP*		sp;
		uint32*	vars;
		Lits_t	clause, tmpclause;
		size_t	added;
		bool	nonbinary_en;

		inline void		write		(const Byte&);
		inline void		write		(const uint32*, const int&);
		inline void		binary		(const uint32*, const int&);
		inline void		nonbinary	(const uint32*, const int&);
		inline void		addline		(const uint32*, const int&);
		inline void		delline		(const uint32*, const int&);
		inline int		exists		(uint32*, const int&);
		inline int		exists		();
		inline void		addClause	();
		inline void		deleteClause();
		inline bool		checkFile	();
	
	public:

		PROOF	();
		~PROOF	();

		size_t numClauses		() const { return added; }
		void close				();
		void init				(SP*);
		void init				(SP*, uint32*);
		void handFile			(arg_t path, const bool&);
		void checkInput			(arg_t input);
		void addEmpty			();
		void addUnit			(uint32);
		void addClause			(Lits_t&);
		void addClause			(CLAUSE&);
		void addClause			(SCLAUSE&);
		void addResolvent		(SCLAUSE&);
		void deleteList			(OL&);
		void deleteClause		(Lits_t&);
		void deleteClause		(CLAUSE&);
		void deleteClause		(SCLAUSE&);
		void deleteClause		(SCLAUSE&, const uint32& def, const uint32& other);
		void shrinkClause		(CLAUSE&);
		void shrinkClause		(CLAUSE&, const uint32&);
		void shrinkClause		(SCLAUSE&, const uint32&);
		void resolve			(const uint32&, OL&, OL&);
		void resolve			(const uint32&, SCLAUSE&, SCLAUSE&);
		bool merge				(const uint32&, SCLAUSE&, SCLAUSE&);
		void printClause		(const char*, const uint32*, const int&, const bool& map = false);

	};

}

#endif