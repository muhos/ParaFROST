/***********************************************************************[proof.cpp]
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

#include "proof.h"
#include "sort.h"

#if defined(_WIN32) || defined(__CYGWIN__)
#define writebyte putc
#elif defined(__linux__)
#define writebyte putc_unlocked
#endif

#define BYTEMAX 128
#define IBYTEMAX -128
#define BYTEMASK 127
#define L2B(x) (((x)&BYTEMASK) | BYTEMAX)
#define ISLARGE(x) ((x)&IBYTEMAX)

using namespace ParaFROST;

PROOF::PROOF() : proofFile(NULL), sp(NULL), vars(NULL), added(0), nonbinary_en(false) {}

PROOF::~PROOF() {
  clause.clear(true);
  tmpclause.clear(true);
  close();
  vars = NULL, sp = NULL, proofFile = NULL;
}

void PROOF::close() {
  if (proofFile != NULL) {
    fclose(proofFile);
    proofFile = NULL;
  }
}

void PROOF::handFile(arg_t path, const bool& _nonbinary_en) {
  PFLOGN2(1, " Handing over \"%s%s%s\" to the proof system..", CREPORTVAL, path, CNORMAL);
  proofFile = fopen(path, "w");
  if (proofFile == NULL) PFLOGE("cannot open proof file %s", path);
  nonbinary_en = _nonbinary_en;
  PFLENDING(1, 5, "(binary %s)", nonbinary_en ? "disabled" : "enabled");
}

void PROOF::init(SP* _sp) {
  assert(_sp);
  sp = _sp;
}

void PROOF::init(SP* _sp, uint32* vorg) {
  assert(_sp);
  assert(vorg);
  sp = _sp;
  vars = vorg;
}

inline bool PROOF::checkFile() {
  if (proofFile == NULL) {
    PFLOGEN("proof file is not opened or cannot be accessed");
    return false;
  }
  return true;
}

inline void PROOF::write(const Byte& byte) { writebyte(byte, proofFile); }

inline void PROOF::write(const uint32* lits, const int& len) {
  assert(checkFile());
  assert(len > 0);
  if (nonbinary_en)
    nonbinary(lits, len);
  else
    binary(lits, len);
  added++;
}

inline void PROOF::addline(const uint32* lits, const int& len) {
  assert(checkFile());
  if (!nonbinary_en) write('a');
  write(lits, len);
}

inline void PROOF::delline(const uint32* lits, const int& len) {
  assert(checkFile());
  write('d');
  if (nonbinary_en) write(' ');
  write(lits, len);
}

inline void PROOF::binary(const uint32* lits, const int& len) {
  for (int i = 0; i < len; i++) {
    const uint32 lit = lits[i];
    CHECKLIT(lit);
    uint32 r = V2DEC(vars[ABS(lit)], SIGN(lit));
    assert(r > 1 && r < NOVAR);
    Byte b;
    while (ISLARGE(r)) {
      b = L2B(r);
      write(b);
      r >>= 7;
    }
    b = r;
    write(b);
  }
  write(0);
}

inline void PROOF::nonbinary(const uint32* lits, const int& len) {
  const size_t bytes = 16;
  const char delimiter = '0';
  char line[bytes];
  char* tail = line + bytes;
  *--tail = 0;
  for (int i = 0; i < len; i++) {
    const uint32 lit = lits[i];
    CHECKLIT(lit);
    const uint32 mvar = vars[ABS(lit)];
    const LIT_ST sign = SIGN(lit);
    if (sign) write('-');
    char* nstr = tail;
    assert(!*nstr);
    uint32 digit = mvar;
    while (digit) {
      *--nstr = (digit % 10) + delimiter;
      digit /= 10;
    }
    while (nstr != tail) write(*nstr++);
    write(' ');
  }
  write(delimiter);
  write('\n');
}

void PROOF::addEmpty() {
  assert(checkFile());
  if (nonbinary_en) {
    write('0');
  } else {
    write('a');
    write(0);
  }
  added++;
}

inline void PROOF::addClause() {
  addline(clause, clause.size());
  clause.clear();
}

inline void PROOF::deleteClause() {
  delline(clause, clause.size());
  clause.clear();
}

void PROOF::addUnit(uint32 unit) { addline(&unit, 1); }

void PROOF::addClause(Lits_t& c) { addline(c, c.size()); }

void PROOF::addClause(CLAUSE& c) { addline(c, c.size()); }

void PROOF::addClause(SCLAUSE& c) { addline(c, c.size()); }

void PROOF::deleteClause(Lits_t& c) { delline(c, c.size()); }

void PROOF::deleteClause(CLAUSE& c) { delline(c, c.size()); }

void PROOF::deleteClause(SCLAUSE& c) { delline(c, c.size()); }

void PROOF::deleteClause(SCLAUSE& c, const uint32& def, const uint32& other) {
  assert(clause.empty());
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i];
    CHECKLIT(lit);
    if (NEQUAL(lit, def))
      clause.push(lit);
    else
      clause.push(other);
  }
  assert(clause.size() == c.size());
  deleteClause();
}

void PROOF::deleteList(OL& list) {
  for (int i = 0; i < list.size(); i++)
    deleteClause(*list[i]);
}

void PROOF::addResolvent(SCLAUSE& c) {
  assert(clause.empty());
  assert(sp != NULL);
  LIT_ST* values = sp->value;
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i];
    CHECKLIT(lit);
    if (values[lit] > 0) {
      clause.clear();
      return;
    }
    clause.push(lit);
  }
  addClause();
}

void PROOF::shrinkClause(CLAUSE& c, const uint32& me) {
  assert(clause.empty());
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i];
    CHECKLIT(lit);
    if (NEQUAL(lit, me)) clause.push(lit);
  }
  assert(clause.size() == c.size() - 1);
  addClause();
  deleteClause(c);
}

void PROOF::shrinkClause(SCLAUSE& c, const uint32& me) {
  assert(clause.empty());
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i];
    CHECKLIT(lit);
    if (NEQUAL(lit, me)) clause.push(lit);
  }
  assert(clause.size() == c.size() - 1);
  addClause();
  deleteClause(c);
}

void PROOF::shrinkClause(CLAUSE& c) {
  assert(clause.empty());
  assert(sp != NULL);
  int* levels = sp->level;
  for (int i = 0; i < c.size(); i++) {
    const uint32 lit = c[i], v = ABS(lit);
    CHECKLIT(lit);
    if (!levels[v]) {
      assert(!sp->value[lit]);
      continue;
    }
    clause.push(lit);
  }
  assert(clause.size() > 1);
  addClause();
  deleteClause(c);
}

void PROOF::checkInput(arg_t input) {
  PFLOGN2(1, " Handing test clause \"%s%s%s\" to the proof system..", CREPORTVAL, input, CNORMAL);
  while (*input) {
    while (isspace(*input)) input++;
    uint32 v = 0, s = 0;
    if (*input == '-')
      s = 1, input++;
    else if (*input == '+')
      input++;
    while (isdigit(*input)) v = v * 10 + (*input++ - '0');
    if (v) {
      uint32 lit = V2DEC(v, s);
      CHECKLIT(lit);
      tmpclause.push(lit);
    }
  }
  Sort(tmpclause, LESS<uint32>());
  PFLDONE(1, 5);
}

void PROOF::printClause(const char* msg, const uint32* lits, const int& len, const bool& map) {
  PFLOGN1(" Proof: %s clause(", msg);
  for (int i = 0; i < len; i++) {
    int lit = int(ABS(lits[i]));
    PRINT("%4d ", SIGN(lits[i]) ? -lit : lit);
  }
  PRINT(")\n");
  if (map) {
    PFLOGN1(" Proof: %s mapped clause(", msg);
    for (int i = 0; i < len; i++) {
      const uint32 lit = lits[i];
      const uint32 mvar = vars[ABS(lit)];
      PRINT("%4d ", SIGN(lit) ? -int(mvar) : int(mvar));
    }
    PRINT(")\n");
  }
}

inline int PROOF::exists() {
  if (tmpclause.empty()) return -1;
  if (tmpclause.size() != clause.size()) return 1;
  Sort(clause, LESS<uint32>());
  for (int i = 0; i < clause.size(); i++)
    if (clause[i] != tmpclause[i])
      return 1;
  printClause("found", clause.data(), clause.size());
  return 0;
}

inline int PROOF::exists(uint32* lits, const int& len) {
  if (tmpclause.empty()) return -1;
  if (tmpclause.size() != len) return 1;
  Lits_t tmp(len);
  tmp.copyFrom(lits, len);
  Sort(tmp.data(), len, LESS<uint32>());
  for (int i = 0; i < len; i++)
    if (tmp[i] != tmpclause[i])
      return 1;
  printClause("found", tmp, len);
  return 0;
}

bool PROOF::merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2) {
  assert(clause.empty());
  const int n1 = c1.size(), n2 = c2.size();
  int it1 = 0, it2 = 0;
  while (it1 < n1 && it2 < n2) {
    const uint32 lit1 = c1[it1], lit2 = c2[it2];
    const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
    if (v1 == x)
      it1++;
    else if (v2 == x)
      it2++;
    else if ((lit1 ^ lit2) == NEG_SIGN)
      return false;
    else if (v1 < v2) {
      it1++, clause.push(lit1);
    } else if (v2 < v1) {
      it2++, clause.push(lit2);
    } else { // repeated literal
      assert(lit1 == lit2);
      it1++, it2++;
      clause.push(lit1);
    }
  }
  while (it1 < n1) {
    const uint32 lit1 = c1[it1++];
    if (NEQUAL(ABS(lit1), x)) clause.push(lit1);
  }
  while (it2 < n2) {
    const uint32 lit2 = c2[it2++];
    if (NEQUAL(ABS(lit2), x)) clause.push(lit2);
  }
  return true;
}

void PROOF::resolve(const uint32& x, OL& poss, OL& negs) {
  assert(clause.empty());
  assert(x);
  for (int i = 0; i < poss.size(); i++) {
    if (poss[i]->original()) {
      for (int j = 0; j < negs.size(); j++) {
        if (negs[j]->original())
          resolve(x, *poss[i], *negs[j]);
      }
    }
  }
}

void PROOF::resolve(const uint32& x, SCLAUSE& c1, SCLAUSE& c2) {
  assert(clause.empty());
  assert(x);
  assert(c1.original());
  assert(c2.original());
  if (merge(x, c1, c2))
    addClause();
  else
    clause.clear();
}