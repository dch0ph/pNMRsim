#ifndef cmatrix_hash_set_h_
#define cmatrix_hash_set_h_

//Quick and dirty string hash set
//Not a fully encapsulated version of set

//need this for definition of stringhash
#include "cmatrix_hash_map.h"
#include <set>

namespace libcmatrix {

  template<class Hasher =stringhash> class hashed_set : public std::set<size_t> {
    Hasher hasher_;
    typedef std::set<size_t> set_t;
  public:
    size_t count(const typename Hasher::argument_type& name) const { return set_t::count(hasher_(name)); }
    std::pair<set_t::iterator, bool> insert(const typename Hasher::argument_type& name) { return set_t::insert(hasher_(name)); }
  };
  
}

#endif
