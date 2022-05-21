#ifndef NMRsim_matrixspec_h_
#define NMRsim_matrixspec_h_

#include "NMRsim.h"

class ExchangeMatrix : public Setable
{
public:
  ExchangeMatrix(const char*, const BaseList<double>&, bool, bool =false);
  //  ExchangeMatrix(const char*, const BaseList<double>&, const BaseList<double>&, bool);

  void printvariablename(std::ostream&, subsid_t) const;
  void printraw(std::ostream&) const; //!< print in input format
  void set(double v, subsid_t subsid) { set(BaseList<double>(1,&v),subsid); }
  void set(const BaseList<double>&, subsid_t);
  size_t sites() const { return matrix_.rows(); } //!< return number of exchanging sites (0 if unset)
  const Matrix<complex>& operator()() const { return matrix_; } //!< return exchange matrix
  bool isexchange() const { return !isgeneral_; }
  bool isreal() const { return isreal_; }
  //  bool isreal() const { return rawvalsi_.empty(); } //!< return true if no imaginary component (or if matrix is undefined in this context)
  bool isconstant() const { return isconst_; } //!< return true if contents are not constant
  bool operator!() const { return matrix_.empty(); } //!< return false if matrix is undefined
  void update();  //!< ensure exchange matrix is up-to-date (against rawvals)

private:
  std::string name_; //!< matrix name
  size_t sites_; //!< number of exchanging sites
  LIST<double> rawvals_; //!< unique raw elements
  bool isreal_; //!< true if last specification was purely real
  //  LIST<double> rawvalsi_; //!< unique raw elements (imag)
  bool isconst_;
  bool isgeneral_;
  Matrix<complex> matrix_; //!< full matrix

  void checkset(size_t, size_t, double);
  static ThreadWarning<> unphysical_warning;
};

std::ostream& operator<< (std::ostream&, const ExchangeMatrix&);

//! stores reference to a "matrix" of some description

//! object linking product operator specification with destination operator
struct operator_def {
  operator_def(const setableoperator_spec* specpv, BlockedOperator& opv)
    : specp(specpv), op(opv) {}
  const setableoperator_spec* specp;
  BlockedOperator& op;
};

//typedef ListList<int> filter_spec;
//typedef UnionHolder< 2, filter_spec, spinorder_spec > generalfilter_spec;

//! object linking filter specification with resulting matrix
struct filter_def {
/*   filter_def(const filter_spec& specv, const BlockedFilter* filterpv =NMRSIM_NULL) */
/*     : spec(specv), filterp(filterpv) {} */
  filter_def(const spinorder_spec& specv, const BlockedFilter* filterpv =NMRSIM_NULL)
    : spec(specv), filterp(filterpv) {}

  const spinorder_spec spec;
  const BlockedFilter* filterp;
  mutable Matrix<bool> filterfull; //!< cache for expanded
  mutable block_pattern blkspec; //!< block pattern of *expanded* filter

  bool operator!() const { return (filterp==NMRSIM_NULL); }
  bool isspinorder() const { return spec.isspinorder(); }
  // { return spec.istype(Type2Type<spinorder_spec>()); }
  bool ensure(bool needed =true); //!< return true if able to construct, and check if current
  const Matrix<bool>& ensurefull() const;
  const BlockedFilter& operator()() const;
};

/** As only storing pointers, class uses a crude pointer to void 
    and casts.  Could be more elegant, but if it works ... */
struct matrixspec {
  enum matrix_t { NONE, //!< not set
		  BOOL, //!< filter matrix
		  OPERATOR, //!< complex matrix
		  //		  DOUBLE, //!< real matrix
		  SPECIAL, //!< "special" e.g. Hamiltonian
		  EXCHANGE //!< (real) exchange matrix
  };

  const void* ptr; //!< pointer to object
  matrix_t type; //!< matrix type

  matrixspec() : ptr(NMRSIM_NULL), type(NONE) {} //!< default constructor creates unset object

  //! only useful for SPECIAL type
  matrixspec(const void* ptr_, matrix_t type_)
    : ptr(ptr_), type(type_) {
    if (type_!=SPECIAL)
      std::cerr << "Warning: explicit matrixspec construction - very dodgy!\n";
  } 

  explicit matrixspec(operator_def& ptr_) //! don't take const to emphasise that it must remain active
    : ptr(static_cast<const void*>(&ptr_)), type(OPERATOR) {} //!< create reference to operator matrix
  //  explicit matrixspec(const BlockedMatrix<complex>& ptr_)
  //  : ptr(static_cast<const void*>(&ptr_)), type(COMPLEX) {} //!< create reference to complex matrix
  //explicit matrixspec(const BlockedMatrix<double>& ptr_)
  //  : ptr(static_cast<const void*>(&ptr_)), type(DOUBLE) {} //!< create reference to real matrix
  explicit matrixspec(const ExchangeMatrix& ptr_)
    : ptr(static_cast<const void*>(&ptr_)), type(EXCHANGE) {} //!< create reference to exchange matrix
  explicit matrixspec(const filter_def& ptr_)
    : ptr(static_cast<const void*>(&ptr_)), type(BOOL) {} //!< create reference to filter (boolean) matrix

  bool isconstant() const;
  bool operator!() const; //!< return \c true if matrix is empty
  usage_t usage() const; //!< return memory usage

  void set_operator(const operator_def&); //!< assign to operator matrix
  //void set_complex(const BlockedMatrix<complex>&); //!< assign to complex matrix
  //  void set_double(const BlockedMatrix<double>&); //!< assign to real matrix
  void set_filter(const filter_def&); //!< assign to filter matrix
  void set_exchange(const ExchangeMatrix&); //!< assign to exchange matrix
  matrixspec& operator= (const matrixspec&); //!< assignment

  const operator_def& asoperator() const; //!< return (if possible) operator matrix
  //  const BlockedMatrix<complex>& ascomplex() const; //!< return (if possible) complex matrix
  //const BlockedMatrix<double>& asdouble() const; //!< return (if possible) real matrix
  const filter_def& asfilter() const; //!< return (if possible) filter matrix
  const ExchangeMatrix& asexchange() const; //!< return (if possible) exchange/transformation matrix

  void log(logfile_controller&, const char* name, int flags =0) const; //!< output to "log file"
  void update(); //!< ensure up-to-date

  //  void assign(BlockedOperator&) const;
  static ThreadWarning<> writefailed_warning; //!< failed to write (undefined?)

  template<class T> static usage_t usage_(const T& a) {
    return !a ? usage_t() : usage_t(a.row());
  }
};

//! true if two ::matrixspec refer to same matrix
inline bool operator== (const matrixspec& a, const matrixspec& b)
{ return (a.ptr==b.ptr); }

//! \c true if two ::matrixspec refer to different matrices
inline bool operator!= (const matrixspec& a, const matrixspec& b)
{ return (a.ptr!=b.ptr); }

std::ostream& operator<< (std::ostream&, const matrixspec&);

typedef MAPTYPE(matrixspec) matrixmap_type; 
extern matrixmap_type matrixmap; //!< map of matrices
			    
			    const filter_def& getfiltermatrix(const char*);
#endif
