#ifndef MAGRES_h_
#define MAGRES_h_

// size of buffer to allocate for error messages
#define MAGRES_ERRM_MAXLEN 1024

#include <string>
#include <string.h>

#ifdef USE_LIBCMATRIX
#include "geometry.h"
#include "List.h"
#else
#include <list>
#endif

// The atomic properties

namespace MagRes {

#ifdef USE_LIBCMATRIX
  typedef libcmatrix::vector3 vector_t;
  typedef libcmatrix::rmatrix3 tensor_t;
#define LIST_T libcmatrix::List
  typedef libcmatrix::Failed exception_t;
  typedef libcmatrix::InvalidParameter notmagres_exception_t;
#else
  typedef double vector_t [3];
  typedef double tensor_t [3][3];
#define LIST_T std::list
  typedef std::runtime_error exception_t;
  typedef std::runtime_error notmagres_exception_t; //!< should use different type so can distinguish failure modes
#endif

  static const char line_separator='\n';
  static const char token_separator[] = " \t\r";
  static const int major_version_limit=1;

struct atomid_t {
  atomid_t(const char* speciesv, size_t indexv)
    : species(speciesv), index(indexv) {}

  std::string species;
  size_t index;
};

struct MagresLine;

struct MagresAtom {
  MagresAtom(const MagresLine&);
  size_t index;
  std::string species; 
  std::string label;
  vector_t position;

  std::ostream& printname(std::ostream& ostr) const { return ostr << species << ' ' << index; };
  friend std::ostream& operator<< (std::ostream& ostr, const MagresAtom& atom);

  bool operator== (const atomid_t& id) const 
  { return (id.index==index) && (strcmp(species.c_str(),id.species.c_str())==0); }

};

struct MagresLattice {
  tensor_t lattice;
  MagresLattice(const MagresLine&);
};

struct MagresSymmetry {
  std::string symmetry_string;
  MagresSymmetry(const MagresLine&);
};

struct MagresFile;

// The magnetic resonance properties
struct MagresMs {
  const MagresAtom *atom;
  tensor_t sigma;
  MagresMs(const MagresLine&, MagresFile&);
};

struct MagresIsc {
  const MagresAtom* atom1;
  const MagresAtom* atom2;
  tensor_t K;
  MagresIsc(const MagresLine&, MagresFile&);
};

struct MagresEfg {
  const MagresAtom* atom;
  tensor_t V;
  MagresEfg(const MagresLine&, MagresFile&);
};

enum block_t { ATOM, MAGRES };

// The whole file
struct MagresFile {
  MagresFile() : latticep(NULL) {}
  ~MagresFile() { delete latticep; }

  MagresLattice* latticep;
  
  size_t num_atoms() const { return atoms.size(); }
  LIST_T<MagresAtom> atoms;

  size_t num_symmetries() const { return symmetries.size(); }
  LIST_T<MagresSymmetry> symmetries;

  size_t num_ms() const { return ms.size(); }
  LIST_T<MagresMs> ms;

  size_t num_isc() const { return isc.size(); }
  LIST_T<MagresIsc> isc;

  size_t num_efg() const { return efg.size(); }
  LIST_T<MagresEfg> efg;

  const MagresAtom* find_atom(const atomid_t&) const;
  const MagresAtom* parse_find_atom(const MagresLine& line,size_t base) const;
  void parse_from_string(char* file);
  void parse_from_file(const char* file);
  void parse_atom(const MagresLine&);
  void parse_magres(const MagresLine&);
  void parse_lines(char *block, block_t blocktype);
  char* windforward(char*);

  friend std::ostream& operator<< (std::ostream& ostr, const MagresFile& magres_file);
};

// Generic line tokenizer
struct MagresLine {
  MagresLine(char*);

  size_t num_cols() const { return cols.size(); }
  const char* column(size_t i) const { return cols(i).c_str(); }
  void verify_columns(size_t cols, const char* label) const;
  void read_tensor(tensor_t& dest, size_t base) const;

  LIST_T<std::string> cols;
};

} //namespace MagRes

#endif
