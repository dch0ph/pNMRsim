#ifndef NMRsim_common_h_
#define NMRsim_common_h_

// 10/2/16  Removed in favour of configure --disable-CXX11
//#define NMRSIM_USE_CXX11 1

/*! \file NMRsim_common.h
 \brief  Header file required by component parsing pNMRsim expressions etc.
*/

//! default error estimate is 10%
#define NMRSIM_DEFAULT_ERROR 0.1

#include <string>
#include <cstring>
#include <iostream>
#include <map>

#include "config.h"

#ifdef NMRSIM_USE_CXX11
#define NMRSIM_USE_HASH 1
#define NMRSIM_NULL nullptr
#endif

#ifndef NMRSIM_NULL
#define NMRSIM_NULL NULL
#endif


#if NMRSIM_USE_HASH
//#include "cmatrix_hash_map.h"
//! type used for directive maps
#include <unordered_map>

struct lcm_string_hash : public std::unary_function<const char *, size_t>
{
  size_t operator()(const char * str) const
  {
    size_t hash = 5381;
    int c;
    
    while ((c = *str++))
      hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    
    return hash;
  }
};

struct lcm_string_compare {
  bool operator()(const char* str1, const char* str2) const 
  {
    for (; *str1 && *str1 == *str2; ++str1, ++str2) {}
    return (*str1 == *str2);
  }
};

#define FASTMAPTYPE(DATA) std::unordered_map<const char*, DATA, lcm_string_hash, lcm_string_compare >
#define MAPTYPE(DATA) std::unordered_map<std::string, DATA>

#else
#define FASTMAPTYPE MAPTYPE
#define MAPTYPE(DATA) std::map<std::string, DATA>
#endif

#include "List.h"
#include "ListList.h"
#include "Warnings.h"

#ifdef NMRSIM_USE_CXX11
#define NMRSIM_ISNAN std::isnan
#else
#define NMRSIM_ISNAN isnan
#endif

std::ostream& parser_printcontext(std::ostream& =std::cerr, bool printtoken =false); //!< print current parsing position (prior to warning/error)

template<typename T =libcmatrix::Failed> class ContextWarning : public libcmatrix::Warning<T> {
public:
  ContextWarning(const char* namev, libcmatrix::BaseWarning* parentpv, libcmatrix::BaseWarning::warning_t typev =libcmatrix::BaseWarning::Inherit, std::ostream& ostrv =std::cerr) 
: libcmatrix::Warning<T>(namev,parentpv,typev,ostrv) {}

  std::ostream& print_message(const char* extra) const {
    //  if (includecontext_)
    parser_printcontext(libcmatrix::Warning<T>::BaseWarning::ostr_);
    return libcmatrix::BaseWarning::print_message(extra);
  }
//private:
//  bool includecontext_;
};

std::ostream& parser_printthread(std::ostream& ostr =std::cout); //!< print thread ID (if multi-threaded)

template<typename T =libcmatrix::Failed> class ThreadWarning : public libcmatrix::Warning<T> {
public:
  ThreadWarning(const char* namev, libcmatrix::BaseWarning* parentpv, libcmatrix::BaseWarning::warning_t typev =libcmatrix::BaseWarning::Inherit, std::ostream& ostrv =std::cerr) 
: libcmatrix::Warning<T>(namev,parentpv,typev,ostrv) {}

  std::ostream& print_message(const char* extra) const {
    parser_printthread(libcmatrix::Warning<T>::BaseWarning::ostr_);
    return libcmatrix::BaseWarning::print_message(extra);
  }
};

extern libcmatrix::BaseWarning NMRsim_repeat_warning;
extern libcmatrix::BaseWarning NMRsim_once_warning;

//! type used for list/vector
/** \note define necessary as mpi.h may define its own List in the global namespace (bad, bad) */
#define LIST libcmatrix::List

// big include but required for productoperator_spec
#include "BaseMetaPropagation.h"
#include "smartptr.h"

using namespace libcmatrix;

//! maximum deviation from integral value
#define NMRSIM_ROUNDTOL 0.1

//! recurse level that triggers warning (0 if none)
#define NMRSIM_RECURSE_WARN 200
//! recurse level that forces abort (0 if none)
#define NMRSIM_RECURSE_ABORT 10000

//! parse input flags
enum { F_DENYZERO=1, //!< disallow zero
       F_ISPHASE=2, //!< flag phase parameter
       //       F_ISSHIFT=4, //!< flag shift parameter (allow p)
       F_ALLOWMISSING=8, //!< allow empty response (string)
       F_DENYEXPR=16, //!< disallow expression
       F_DENYVAR=32, //!< disallow variable
       F_ALLOWSLOT=64, //!< allow #<n> slot numbers
       F_ALLOWLIST=256, //!< allow list quantity
       F_DENYARRAY=512, //!< disallow {} array
       F_ALWAYSNEW=1024, //!< always create ::VariableBase
       F_ISINTEGER=2048, //!< flag integer quantity
       F_DENYSUM=4096, //!< disallow [] sum
       F_IGNOREDOLLAR=8192, //!< ignore unresolved variable when parsing strings (if checks enabled)
       F_REPLACEDOLLAR=16384, //!< resolve $ when parsing string
       F_IGNORESYNTAX=32768 //!< ignore ? syntax request when parsing string
};

//! attribute flags
enum { 
  A_VAR=32,
  A_ARRAY=512,
  A_SUM=4096
};

const int F_SIMPLE = F_DENYVAR | F_DENYARRAY | F_DENYSUM; // flags for 'simple' simulation
//extern int global_flags; //global settings
extern int F_defaultdataset; //!< default flags for variables structuring data set (sw, np)

//! exit error codes
/** \todo Error codes ought to be replaced by something more sophisticated
 */
enum { ERR_INVALID_INPUT=3, //!< parse failed
       ERR_FAILED=4 //!< logic error
};

//! verbose flags
enum { 
  VER_GEN=1, //!< general verbose
  VER_POWDER=8, //!< powder averaging info (also implies VER_PARALLEL)
  VER_OPTIM=16, //!< optimisation/fitting info
  VER_PARSE=32, //!< debugging of parsing
  VER_PROFILE=64, //!< enabled profiling
  VER_PARALLEL=128 //!< parallelisation info
};
extern int verbose; //!< accumulated verbose flags
extern int verbose_level; //!< verbosity level (0 quiet)
extern bool silent; //!< \c true suppresses all output not specifically enabled by verbose
extern bool nochecks; //!< \c true if checks are disabled
//inline bool allowwarnings() { return (!nochecks && !silent); } //!< \c true if warnings/checks have not been disabled

//limited number of variables/functions that the parser requires
extern int sum_index; //!< index into sum array
extern int row_index; //!< row index into data set
extern int var_index; //!< index to variable {} set
extern int evaluation_index; //!< number of simulation runs
extern LIST<size_t> suminds; //!< sum virtual indices 
extern size_t nchannels; //!< number of RF channels
extern LIST<size_t> nucids; //!< nucleus for each RF channel (from \c channels)

//interfacing functions between parser and NMR system
size_t parser_verbose_level(); //!< return verbose level for parser
void fudge_shift(double&, size_t =0); //!< convert shift from ppm to frequency
//void make_common_par_variables(); //!< variables for \c par shared by all readers
void error_abort(const char*, int =ERR_INVALID_INPUT); //!< exit with error
void error_abort(int =ERR_INVALID_INPUT); //!< abort

double handle_variable(int flags,size_t subsid =1); //!< create variable using input flags \a flags
double handle_operator_variable(int,size_t,size_t); //!< create variable co-efficient in ::setableoperator_spec
const basespin_system* get_spin_system(); //!< get spin system
bool proton_freq_isconstant(); //!< returns \c true if proton Larmor frequency is fixed during calculation
extern double curgrat; //!< current ratio of gamma/gamma1H (0 if not valid)
extern double defgrat; //!< default current ratio (0 unless homonuclear)
extern double proton_freq; //!< proton frequency (Hz, 0 if unset)
FILE* pathopen(const char* fname, const char* mode); //!< open file using global path
extern LIST<const char*> argnamestack;

//! flags for substitute_string
enum { SUB_NUMERIC=1, //!< replace substitution variables ($1 etc.) otherwise user variables will be replaced 
       SUB_ESCAPE=2, //!< replace \\$ by $
       SUB_ABORTHASH=4, //!< don't substitute $1 etc. if first character is hash i.e. line is comment
       SUB_NONCONSTWARN=8, //!< warn if formatting non-constant variables
       SUB_FULLPREC=16 //!< use full precision when evaluating numerics
};
extern size_t substitute_maxchars; //!< maximum number of output characters for a single list item
char* substitute_string(char* out, int, const char* in, int flags =0); //!< substitute $ variables
char* substitute_string(char* out, int, const char* in, int flags, int& accuses); //!< substitute $ variables

class VariableBase;

//! (abstract) base class for potentially variable parameters
class RealVariable {
public:
  RealVariable(const char* namev, bool issetv, bool isconstv, int usesv =0)
    : name_(namev), varp_(NMRSIM_NULL), issetable_(issetv), isconst_(isconstv), isused_(false), uses_(usesv) {
    if (!namev)
      throw InvalidParameter("RealVariable");
  }
  RealVariable(const char* namev, VariableBase&);
  virtual ~RealVariable() {}
  virtual const BaseList<double> get_list() const =0; //!< return value
  const char* name() const { return name_.c_str(); }
  bool isconstant() const { return isconst_; } //!< returns \c true if parameter fixed
  bool isdefined() const { return isconst_ || (uses_!=0); } //!< returns \c true if variable has been defined, either to constant value or some kind of variable
  bool issetable() const { return issetable_; } //!< returns \c true if parameter can be altered (by user)
  //! set const status
  void isconstant(bool isconstv) { 
    if (!isconstv && !issetable_)
      throw InternalError("RealVariable: can't be fixed and non-const");
    isconst_=isconstv;
  }
  bool hasbase() const { return (varp_!=NMRSIM_NULL); }
  int uses() const { return uses_; }
  bool isused() const { return isused_; } //!< returns \c true if the variable has been referred to
  void flagused() { //!< flag that variable is linked
    //, but only if variable is not const
    //    if (!isconst_)
      isused_=true;
  } 
  //! return value
  const VariableBase& value() const;
  VariableBase& value();
  bool validate() const; //!< return true if variable is well-defined in this context
  
  void reset(VariableBase&);//, bool isconstv) //!< reset to new value
  void set_error(double, size_t =0U); //!< set error value
  double get_error() const; //!< return current error
  const BaseList<double> get_errors() const; //!< return all errors

private:
  std::string name_; //!< variable name
  VariableBase* varp_; //!< pointer to associated VariableBase (NMRSIM_NULL if none)
  bool issetable_; //!< \c true if parameter can be altered
  bool isconst_; //!< \c true if parameter is fixed
  bool isused_; //!< \c true if parameter has been used
  int uses_; //!< flags indicating what type of variable
};

//! actual type used for Setable component
/*! \note This is an \c int rather than an enum for \c Setables such as ::setableoperator_spec
    where the component is an index running from 0 (-1 is then an error/trap value) */
typedef int subsid_t;

class PhasedSequence;
PhasedSequence* handle_phasedsequence(int flags, subsid_t); //!< create reference to ::PhasedSequence 
std::ostream& operator<< (std::ostream&, const PhasedSequence&);


//! round \c double to nearest integer
/** exit with ::ERR_FAILED if not within ::NMRSIM_ROUNDTOL */
int round_int(double);
inline int rawround_int(double f) { return (f<0.0) ? -int(0.5-f) : int(0.5+f); }

class multioperator_spec : public LIST<productoperator_spec> {
public:
  multioperator_spec() { reset_(); }
  multioperator_spec(const productoperator_spec& spec)
	: LIST<productoperator_spec>(1, spec) { reset_(); }

  void clear() { LIST<productoperator_spec>::clear(); reset_(); }
  void swap(multioperator_spec&);
  void setnucleus(); //!< check for common nucleus
  size_t nucleus() const { return nuc_; }
  size_t arraytag_; //!< virtual dimension tag (0 if none) - shouldn't really be public (ease of parsing)
private:
  size_t nuc_; //!< common nucleus id (NULL_NUCLEUS if none)
  void reset_() { arraytag_ = 0; nuc_ = NULL_NUCLEUS; }
};

class Expression;
class ExpressionNamedBase;

//! base class for expression components
class ExpressionBase {
public:
  ExpressionBase() {}

  //! single child component
  explicit ExpressionBase(size_t child)
    : children(1,child) {}

  //! two child component
  ExpressionBase(size_t left, size_t right)
    : children(2)
  { children.front()=left; children(1U)=right; }

  //! general number of children
  ExpressionBase(const BaseList<size_t>& childrenv)
    : children(childrenv) {}
  
  virtual ~ExpressionBase() {}
  virtual ExpressionBase* clone() const =0;
  virtual void print(std::ostream&, const Expression&) const =0; //!< stream expression from here
  virtual void get(LIST<double>&, const Expression&) const =0; //!< get value from here
  double get(const Expression& expr) const; //!< get (scalar) value from here
  virtual bool isinteger() const { return false; } //!< \c true if result is integer (rather than real)
  virtual bool isconstant(int&, const Expression&) const; //!< \c true if sub-expression is constant
  bool isconstant(const Expression& expr) const { int junk; return isconstant(junk,expr); }
  int uses(const Expression&) const; //!< return set of flags indicating what type of variables used
  virtual int uses() const { return 0; }

  friend class VariableBuilder;
  friend class Expression;

protected:
  LIST<size_t> children; //!< offsets to children

  void setchildren(size_t start, size_t end); //!< explicitly set child pointers
  // const_t isconst;
  // void setisconst(const_t constv) { isconst=constv; } //!< \internal
};

class ExpressionNamedBase : public ExpressionBase
{
public:
  ExpressionNamedBase(const char* namev, const BaseList<size_t>& offs)
    : ExpressionBase(offs), name_(namev) {}

  explicit ExpressionNamedBase(const char* namev)
    : name_(namev) {}
  
  const char* name() const { return name_; }

protected:
  const char* name_; //!< function name (pointer to spec)
};

//! class for storing expressions
class Expression : public LIST< smartptr<const ExpressionBase> >
{
public:
  typedef LIST< smartptr<const ExpressionBase> > support_type;
  Expression() {}
  Expression(const ExpressionBase&, const BaseList< LIST<double> >& args); //!< simple construction from single operator and arguments

  //! returns \c true if expression is constant
  bool isconstant() const { return back()->isconstant(*this); }
  int uses() const { return back()->uses(*this); }

  //! add node to expression true, returning offset from root
  size_t push(const ExpressionBase* op) {
    push_back(smartptr<const ExpressionBase>(op));
    return size()-1;
  }
  size_t push_original(ExpressionBase* op) {
    push_back(smartptr<const ExpressionBase>());
    back().reset(op);
    return size()-1;
  }

  void get(LIST<double>&) const; //!< return expression value (list)
  void print(std::ostream& =std::cout) const; //!< stream expression
};

inline std::ostream& operator<< (std::ostream& ostr, const Expression& a)
{ a.print(ostr); return ostr; }

//! signature for function
/** name + number of arguments (-1 multivariate) */

struct function_spec : public std::pair<const char*,int> {
  function_spec(const char* name, int nargs)
    : std::pair<const char*,int>(name,nargs) {}

  bool operator== (const function_spec& b) const
  { return ((second == b.second) && (strcmp(first,b.first)==0)); }
  
  bool operator!= (const function_spec& b) const
  { return ((second != b.second) || (strcmp(first,b.first)!=0)); }

};

//! user defined function definition
struct function_def_t {
  const ExpressionNamedBase* objp; //!< pointer to "parent" object from which copies are cloned
  bool used; //!< \c true if function is `in use'
  explicit function_def_t(const ExpressionNamedBase* objpv =NMRSIM_NULL)
    : objp(objpv), used(false) {}
  const ExpressionNamedBase& operator*() const; //!< resolve pointer to source function object
};

#if NMRSIM_USE_HASH
//! hash for function signature
class FunctionHash : public std::unary_function<function_spec,size_t> {
 public:
  //! \note Not a perfect hash, but can't confuse functions with same name or same number of args */
  size_t operator()(const function_spec& a) const;
 private:
  lcm_string_hash stringhasher_;
};
typedef std::unordered_map<function_spec,function_def_t, FunctionHash> Function_Factory_t;
#else

struct FunctionCompare {
  bool operator()(const function_spec&, const function_spec&) const;
};

typedef std::map<function_spec,function_def_t,FunctionCompare> Function_Factory_t;
#endif

extern Function_Factory_t& get_Function_Factory(); //!< returns registry for function names

//! class used for flexible storage
class VariableBase {
public:
  VariableBase() : isconst_(false) {}
  explicit VariableBase(double val, double err =0.0); //!< construct from scalar
  explicit VariableBase(const BaseList<double>& val) : value_(val), isconst_(true) {} //!< construct from list
  explicit VariableBase(const Expression& exprv) : expr_(exprv) { update(); } //!< construct from expression

    bool empty() const { return value_.empty() && expr_.empty(); } //!< contents undefined
  bool isarray() const { return !(array.empty()); } //!< returns \c true if array
  bool isexpr() const { return (!expr_.empty()); } //!< returns \c true if expression
  bool isexprarray() const { return !(arrayexprs_.empty()); } //!< returns \c true if contains array with expressions
  bool issimple() const { return array.empty() && expr_.empty(); } //!< returns \c true if simple scalar
  bool isconst() const { return isconst_; } //!< returns \c true if fixed
  void clear(); //!< clear contents

  size_t array_item_length() const { return array.empty() ? 0 : arraylength_; } 
  //! set to scalar
//   VariableBase& operator= (double val) { 
//     initialise(val);
//     return *this;
//   }
  VariableBase& operator= (const BaseList<double>&); //!< set to list
  VariableBase& operator= (const Matrix<double>&); //!< create from matrix

  //  void initialise(double initvalv); //!< initialise to double
  void initialise(double initvalv, double errv =0.0); //!< initialise to double
  void initialise(Expression&); //!< intialise to expression (DESTROYED)
  void initialise(bool isconstv); //!< initialise to empty
  void initialise(size_t rows, size_t cols, bool isconstv =false); //!< initialise to (zero) matrix

  void swapin(VariableBase&); //!< \internal initialise from temporary

  const BaseList<double>& value() const { return value_; }
  LIST<double>& value() { return value_; } //!< return value

  bool updateindex();  //!< update value based on overall index, return \a false if no update possible or varindex didn't correspond to new state (getindex)
  void update(bool allowwarnings =true); //!< update value
  void set(double); //!< set scalar value 
  void set(const BaseList<double>& v) { value_=v; } //!< set list value

  void print(std::ostream&, bool full =true) const; //!< stream value or complete quantity if \a full is \c true

  const Expression& expression() const { return expr_; } //!< return expression
  const ListList<double>& get_array() const { return array; } //!< return array
  const BaseList<double> get_row(size_t i) const; //!< return row i (only valid for array type)
  BaseList<double> get_row(size_t i); //!< return row i (only valid for array type)
  void set_current_row(const BaseList<double>&); //!< set current row of rectangular array
  size_t summationindex() const; //!< get summation block
  size_t sumtag() const { return sumtag_; } //!< sum virtual dimension
  std::pair<size_t,bool> getoffset() const; //!< return current offset into array and flag indicating whether this is a "new" value for this variable
  const BaseList<size_t>& arraytags() const { return arraytags_; } //!< array virtual dimensions
  double get_error(size_t i) const; //!< return error on index \a i
  double get_error() const; //!< return error on current index
  const BaseList<double>& get_errors() const { return errors_; } //!< return errors list
  void set_error(double, size_t =0U); //!< set error value
  int uses() const; //!< return attributes used
  int arrayuses() const; //!< return overall attributes of array *contents*

  friend class VariableBuilder;
  friend void verify_arraysizes(const VariableBase&);

  bool validate_context(const char*) const; //!< warn if variable evaluation is dodgy in this context

private:
  LIST<double> value_; //!< value
  ListList<double> array; //!< array values
  LIST<double> errors_; //!< error values
  LIST<size_t> arraytags_; //!< array virtual dimensions
  LIST<Expression*> arrayexprs_; //!< expressions in array values
  LIST<size_t> exproffsets_; //!< offsets 
  size_t sumtag_; //!< sum virtual dimension
  Expression expr_; //!< expression
  bool isconst_; //!< true if value is constant  
  size_t arraylength_; //!< length of items in {} array
};

inline std::ostream& operator<< (std::ostream& ostr, const VariableBase& a)
{ a.print(ostr); return ostr << '\n'; }

RealVariable& parse_variable_name(const char* name, bool flagused =true); //!< find variable called \a name, automatically flag used if \a flagused set

void clearevalstacks(); //!< \internal
std::ostream& dumparray(std::ostream&, const BaseList<double>&, char ='['); //!< stream array contents
char matchingbracket(char); //!< return complementing bracket

class CompSequenceBase;
typedef MAPTYPE(CompSequenceBase*) seqmap_type;
extern seqmap_type seqmap; //!< map of sequence fragments

template<typename T> void dumpmap(const T& map, const char* name, const char* mapname, bool dumpoptions =true)
{
  parser_printcontext() << mapname << ' ' << name << " not found\n";
  if (!(map.empty()) && dumpoptions) {
    std::cerr << "Defined entries are:";
    const typename T::const_iterator end(map.end());
    typename T::const_iterator start(map.begin());
    while (start!=end)
      std::cerr << ' ' << (*start++).first;
    std::cerr << '\n';
  }
}

//! find quantity \a name in map

template<typename T> typename T::mapped_type findmap(const T& map, const char* name, const char* mapname, bool dumpoptions =true)
{
  const typename T::const_iterator curp(map.find(name));
  if (curp==map.end()) {
    dumpmap(map,name,mapname,dumpoptions);
    error_abort();
  }
  return curp->second;
}

//! "Dirty" status of sequence fragment or overall simulation
/** higher dirty value implies lower ones */
enum dirty_t { DIRTY_CLEAN =0, //!< sequence has not been changed
	       DIRTY_GAMMA, //!< only gamma angle has been changed
	       DIRTY_HAMILTONIANRF, //!< system hamiltonian or RF parameters have changed
	       DIRTY_ALL //!< Sequence will need rebuilding (timing has changed)
}; 

enum context_t {
  CONTEXT_PRE=1, //!< pre simulation  
  CONTEXT_MAINLOOP, //!< within main evaluation loop / per-spectrum processing
  CONTEXT_POSTCALC, //!< processing after summation
  CONTEXT_FINALISE, //!< processing of final results
  CONTEXT_TERMINATE, //!< post all calculation
  CONTEXT_UNCHECKED //!< ignore context validity
};

extern context_t evaluation_state; //!< global flag of evaluation state
extern bool isinteractive; //!< true if in an interactive mode
extern BaseList<complex> current_data_row; //!< reference to current line of data set (empty if outside of processing);

inline bool within_evaluation() { return (evaluation_state==CONTEXT_MAINLOOP); } //!< return \c true if within main simulation loop

bool getindex(size_t& transindex, size_t tag, size_t arraylen); //!< return index into array of length \a arraylen given dimension tag \a tag

void ensure_library_version(const char*); //!< ensure library is at least given version

double getgrat(size_t i); //!< return gamma ratio for nucleus \a i
size_t get_thread_num(); //!< get thread number (always 0 if not threading)

extern size_t global_workers; //!< number of parallel workers (0 if not using parallel)

extern bool debug; //!< debug flag

int gcd(int a, int b); //!< return greatest common divisor of \a a and \a b
int lcm(int a, int b); //!< return lowest common multiple of \a a and \a b

bool try_find_sw(); //!< try to find alternative source of sw info 
bool try_find_np(); //!< try to find alternative source of np info

std::ostream& print_syntax_string(std::ostream&, const char*, size_t);
void verify_syntax(const char* arg, const char* syn, size_t =0); //!< 0 indicates no separate arguments - print whole
// fragments for constructing syntax strings
#define NMRSIM_RANGESTR "<range of spectrum as indices numbered from 1>"
#define NMRSIM_ROWCOLSTR "[<row selection>]#<column selection>"

#endif
