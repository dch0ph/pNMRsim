#ifndef NMRsim_h_
#define NMRsim_h_

/*! \file NMRsim.h
  \brief  Header file required by most pNMRsim components

 (rather an untidy mess...) */

#include "NMRsim_common.h"
#include "MetaPropagation.h"
#include "UnionHolder.h"
#include "smartptr.h"
#include "simple_counter.h"
#include "timer.h"  
#include "optim.h"

#include "config.h"

#include "lcm_binaryutils.h"

typedef unsigned int count_t; //!< counting type - avoid size_t as 64 bit integers not defined for Matlab file format

//! if set, then offset signs are reversed to give more intuitive results
#define NMRSIM_INTUITIVE_OFFSET
//! nuclear spin properties table to use by default
#ifndef NMRSIM_NUCLEUS_PROPERTIES
#define NMRSIM_NUCLEUS_PROPERTIES "HarrisIUPAC"
#endif

#if defined HAVE_LIBMPI && !defined DISABLE_MPI
#define HAVE_PARSYS
#define USEMPI
#else
#ifdef HAVE_FORK_CONTROLLER
#define HAVE_PARSYS
#else
#undef HAVE_PARSYS
#endif
#endif
#ifdef HAVE_PARSYS
#include "lcm_basethreads.h"
#endif

//! Maximum number of dimensions in nD data sets (can't generally exceed LCM_MULTIMATRIX_DIMS)
#define MAX_DIMENSIONS 4
//! Maximum number of RF channels
#define MAXCHANNELS 4
//! default propagator cache limit (M)
#define NMRSIM_CACHE_LIMIT 200
//! cache_limit over 100 Gb is assumed to be mistake!
#define NMRSIM_HUGE_CACHE_LIMIT 1e5

//! Use Chebyshev propagation when matrix is less than this density
#define DEFAULT_CHEBYSHEV_CROSSOVER 0.5

//! don't warn about ignored imaginary components below this level
#define NMRSIM_IMAG_TOLERANCE 1e-8

//! default number of steps in lineshape histograms
#define NMRSIM_DEFAULT_LINESHAPE_STEPS 5

//! default tolerance for non-zero coherence when reporting statistics (override with NMRSIM_COHERENCE_TOLERANCE environment variable)
#define NMRSIM_DEFAULT_COHERENCE_TOLERANCE 1e-8

//! default lineshape cutoff (0 none)
#define NMRSIM_DEFAULT_LINESHAPE_CUTOFF 0.0

//! default number of lineshapes to cache
#define NMRSIM_DEFAULT_LINESHAPES_CACHE NMRSIM_DEFAULT_LINESHAPE_STEPS

//! character to use instead of illegal Matlab
#define NMRSIM_MATLAB_SQUASH_CHAR 'X'

/* NOTHING TO ALTER BELOW HERE */

#define NMRSIM_INVALID -10000

typedef accumulating_timer<>::guard timer_guard_t;

#ifndef NOPERIODIC
#include "CrystalSystem.h"
#endif

using namespace libcmatrix;

#define STRNCMP(A,B) strncmp(A,B,sizeof(B)-1)
#ifdef NDEBUG
#define NMRSIM_EXPECT(A) A
#else
#define NMRSIM_EXPECT(A) assert(A)
#endif

//! Used to specify component of non-trivial Setable's
enum { S_NONE=0, //!< Flag value: \e not setable
       S_ALPHA, S_BETA, S_GAMMA, //!< Euler angle components
       S_PHASE, //!< rotor / RF phase
       S_SPEED, //!< rotor speed
       S_ANGLE, //!< pulse tip angle
       S_OFFSET, //!< RF offset
       S_RF, //!< RF nutation rate
       S_DUR, //!< pulse duration
       S_ISO, //!< isotropic component
       S_ANISO, //!< anisotropy of interaction
       S_ASYM, //!< asymmetry of interaction specified as xx-yy or eta
       S_GFRAC, S_GFRAC1, //!< LG fraction
       S_ARG1, S_ARG2 //!< generic arguments 
 };

//! mix-in class for any quantity that altered during a simulation
class Setable {
public:
  Setable() : uses_(0) {}
  virtual ~Setable() {}

  virtual void printvariablename(std::ostream&, subsid_t) const =0; //!< output variable identifier
  virtual void set(double fval,subsid_t =S_ARG1) =0; //!< set to value method (qualified by parameter type)
  virtual void set(const BaseList<double>&, subsid_t =S_ARG1); //!< set to list of values (not allowed by default)

  int setable_uses() const { return uses_; }
  void setable_uses(int usesv) { uses_= usesv; }

 private:
  int uses_; //!< overall attributes
};

//! mix of ::multioperator_spec and ::Setable to allow co-efficients to be varied
/*! \note A separate class is defined for the mix since ::multioperator_spec needs
  to be visible to the parser which does not know (or want to know) about ::Setable */
class setableoperator_spec : public multioperator_spec, public Setable {
public:
  setableoperator_spec(const operator_spec& spec_, const basespin_system& sys_) 
    :  
    multioperator_spec(productoperator_spec(spec_,sys_)),
      isconst_(true) {} //!< construct ::setableoperator_spec

  setableoperator_spec(const productoperator_spec& spec_)
    : multioperator_spec(spec_), isconst_(true) {}

  setableoperator_spec(const operator_spec& spec_) 
    : 
    multioperator_spec(productoperator_spec(spec_)),
      isconst_(true) {} //!< construct ::setableoperator_spec

  explicit setableoperator_spec(bool isconstv =true) 
    : isconst_(isconstv) {} //!< create empty (which can be added to)

  const productoperator_spec& current() const { return (*this)(currentindex()); } ; //!< return current (array-indexed) product operator specification
  productoperator_spec& current() { return (*this)(currentindex()); } 

  size_t arraytag() const { return arraytag_; }
  void set(double, subsid_t =S_NONE); //!< set coefficient method (subsid_t codes for 
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  bool isconstant() const { return (size()==1) && isconst_; } //!< \c false if operator specification is non-trivial (a list and/or contains variable co-efficients)
  bool operator!() const { return empty(); }
private:
  size_t currentindex() const;  //!< return current index
  bool isconst_; //!< \c false if operator specification contains variable co-efficients
};

//! \c true if two ::setableoperator_spec are "simple" and have a non-null intersection
inline bool arematching(const setableoperator_spec& a, const setableoperator_spec& b)
{
  if (!a.isconstant() || !b.isconstant())
    return false;
  return arematching(a.front(),b.front());
}

//! current processing state for single row/dimension
struct processing_state {
  processing_state(double,bool istd, double sfrqv =0.0, double refv =0.0); //!< initialise from supplied SW
  processing_state() : sw(0.0) {}
  double sw; //!< spectral width (0 if undefined)
  bool istimedomain; //!< \c true if data is in time domain
  double sfrq; //!< NMR frequency
  double ref; //!< centre-of-spectrum 

  bool operator== (const processing_state&) const; //!< very simple minded comparison
  bool operator!= (const processing_state& a) const { return !(*this==a); }
};

std::ostream& operator<< (std::ostream&, const processing_state&);

//! object used to store data set
class DataStore : private UnionHolder< 2, ListList<complex>, cmatrix >  {
public:
  typedef UnionHolder< 2, ListList<complex>, cmatrix > store_type; //!< data storage
  LIST<processing_state> procstates_;

  using store_type::clear;

  void reset(bool istd); //!< create empty data set and reset processing
  void create(size_t r); //!< create empty irregular data set with \a r rows
  void create(size_t r, size_t np, const complex& v); //!< create regular nD data set with \a r rows (total), \a np data points and initial value \a v (uses global sw values)
  void createrow(size_t np, const complex& v, const processing_state&); //!< add new row in irregular data set with \a np points and value \a v
  void push_back(const BaseList<complex>&, const processing_state&); //!< add new row to irregular data set
  bool transpose(); //!< transpose shape

  static Warning<> mismatchedtype_warning;
  //! apply transformation
  template<class F> void apply_ip(const F& func, const DataStore& b)
  {
    if ((rows()!=1) || (b.rows()!=1)) {
      if (isnD()==b.isnD())
	store_type::apply_ip(func,b);
      mismatchedtype_warning.raise();
    }
    BaseList<complex> arow(row()); //!< assume they match!
    func(arow,b.row());
  }

//   double sw(size_t r =0) const; //!< spectral width for given row
//   double& sw(size_t r =0);
//   //!< return list of spectral widths
//   const BaseList<double>& sws() const {
//     assert(sws_.size()==rows()); //!< quick sanity check
//     return sws_;
//   }
  LIST<processing_state>& states() { return procstates_; } //!< return raw list of processing states
  const BaseList<processing_state>& states() const { return procstates_; } 
  processing_state& state(size_t i); //!< return processing state for row or dimension \a i
  const processing_state& state(size_t) const;
  bool empty() const; //!< \c true if data set is empty
  const bool isnD() const { return (type()==2); } //!< \c true if contains nD data
  size_t rows() const; //!< total number of rows
  size_t size() const { return row().size(); } //!< total number of data points
  void print(std::ostream& =std::cout) const; //!< print contents to stream
  void print_structure(std::ostream& =std::cout) const; //!< display structure

  const BaseList<complex> row() const; //!< return data set as row vector
  BaseList<complex> row();
  const BaseList<complex> row(size_t r) const; //!< return row \c r of data set
  BaseList<complex> row(size_t);

  cmatrix trymatrix() const; //!< try to return data as a matrix (exception if not nD)
  cmatrix& matrix(); //!< return data set as simple matrix (fails if not nD data)
  const cmatrix& matrix() const;
  ListList<complex>& listlist(); //!< return data as ListList (fails if not irregular)
  const ListList<complex>& listlist() const;
  //  bool checkprocstates(double =1e-5) const;

  DataStore& operator+= (const DataStore&); //!< (in place) addition of data sets [proc states are ignored]
  DataStore& operator-= (const DataStore&); //!< (in place) subtraction
  DataStore& operator= (const complex& v); //!< set to constant value \a v
  //  DataStore& operator= (const cmatrix&); //!< set to simple matrix
  void assign(const cmatrix&, const BaseList<processing_state>&); //!< replace data with new 2D matrix and state set

  void copy_structure(LIST<size_t>&) const; //!< store row sizes
  void duplicate_structure(const DataStore&); //!< copy structure
  void duplicate_structure_processing(const DataStore&); //!< copy structure and processing

  void print_processing(std::ostream& =std::cout) const; //!< display processing state
};

bool arematching(const DataStore&, const DataStore&); //!< \c true if data sets are `compatible'

//! ensure \a a has the same structure as \a b
inline void duplicate_structure(DataStore& a, const DataStore& b)
{
  a.duplicate_structure(b);
}

inline std::ostream& operator<< (std::ostream& ostr, const DataStore& a) {
  a.print(ostr);
  return ostr;
}

//! forward declarations
struct Ucache_obj;
class MasterObj;

//! type for storing memory usage (items,bytes)
struct usage_t {
  usage_t()
    : items(0), bytes(0) {}
  usage_t(size_t itemsv, long bytesv) 
    : items(itemsv), bytes(bytesv) {}

  template<class T> usage_t(size_t itemsv, Type2Type<T>) //!< usage of \a itemsv objects of type T
    : items(itemsv), bytes(long(itemsv)*sizeof(LCM_VAL(T))) {}

  template<class T> explicit usage_t(const T& a) //!< usage of container object \a a
    : items(a.size()), bytes(long(items)*sizeof(LCM_VAL(T))) {}

  usage_t operator+ (const usage_t& a) const {
    return usage_t(items+a.items,bytes+a.bytes);
  }
  usage_t& operator+= (const usage_t& a) {
    items+=a.items; bytes+=a.bytes;
    return *this;
  }
  operator bool() const { return (bytes!=0); } //!< object is \c true if usage is non-zero

  size_t items; //!< number of ``items'' (not well defined for mixed types)
  long bytes; //!< bytes used
};

std::ostream& operator<< (std::ostream&, const usage_t&);

//! forward declarations
class logfile_controller;

//! Fake object to allow timing of operations not associated with action / processing objects
struct InternalCommand : public accumulating_timer<> {
  InternalCommand(const char* namev) : name(namev) {}

  double time() const { return operator()(); } //!< return time used
  void print(std::ostream& ostr) const { ostr << name; }
  
  const usage_t& usage() const { return usage_; } //!< return memory usage
  void usage(const usage_t& usagev) { usage_=usagev; } //!< set memory usage

  std::string name; //!< description of activity
  usage_t usage_; //!< memory usage
};

extern LIST<InternalCommand*> internal_timers; //!< list of internal timing objects
InternalCommand* create_timer(const char* name); //!< create new timing object called \a name
extern InternalCommand* setup_cp;

std::ostream& operator<< (std::ostream& ostr, const accumulating_timer<>&);

//! state of optional feature
class option {
 public:
  enum optional_t { OFF, //!< feature disabled
	 AUTO, //!< no explicit instructions - use feature if appropriate
	 ON, //!< force feature on if possible
	 TEST //!< test effect of flag
  };

  enum used_t { USED, NOTUSED, NOTKNOWN };

  option(const char* labelv, const char* descv ="", optional_t def =AUTO, used_t =NOTKNOWN);

    bool operator() () const { return isnotstate(OFF); }
    bool operator! () const { return isstate(OFF); }
    bool isenabled() const { return isstate(ON); }
    bool isdisabled() const { return isstate(OFF); }
    //   bool checkisdisabled();
    bool isauto() const { return isstate(AUTO); }
    bool istest() const { return isstate(TEST); }
    bool isused() const { return (used_==USED); }
    bool isstate(optional_t) const;
    bool isnotstate(optional_t) const;
    void setdefault(optional_t);
    const char* name() const { return label_.c_str(); }
    const char* description() const { return desc_.c_str(); }
    void check(bool used, const char* desc =NMRSIM_NULL);
    void setusage(bool);
    used_t getusage() const { return used_; }
    optional_t get() const { return value_; }
    void set(optional_t v) { value_=v; }

 private:
    std::string label_;
    std::string desc_;
    optional_t value_;
    used_t used_;
};

//bool isopton(optional_t);
//bool isoptoff(optional_t);
//bool isoptauto(optional_t);

extern LIST<const char*> testedoptions; //!< list of options to be tested (safer to use name rather than pointer into map)

//! map of optional features 
/** \note little point in hashing name since optional flags are usually only checked once */
typedef std::map<std::string,option*> optional_map_t;
void add_option(option&); //!< register option
optional_map_t& get_optional_map(); //!< get optional map
void setdefaultoptions(option::optional_t =option::AUTO); //!< set default option status 
extern option optcache; //!< cache propagators feature
extern option optpartition; //!< partitioning feature (for Chebyshev propagation)
extern option optupdate; //!< minimal updating
extern option optsync; //!< try to synchronise
extern option optgeneralisedQ; //!< generalised quadrupole treatment
extern option optparallel; //!< use parallel facilities
extern option optcombinepropagators;
extern option optsmartprop;
extern option optgamma;
extern option optforcepointbypoint;
extern option optphasemod;
extern option optclassicQ;

//void check_optional(const option& enable, const char* desc, bool isused, const char* name =NMRSIM_NULL); //!< check actual used status (\a isused) against user request (\a enable)

extern bool abortonwarning; //!< if \c true abort on precision / accuracy warning
extern bool noexecute; //!< \c true if only parsing
extern bool active2D; //!< \c true for nD simulations (any indirect dimension set)
extern double spin_rate; //!< MAS spin rate (0 if static)
extern double chebyshev_crossover; //!< density below which to use Chebyshev propagation (NOT USED)
extern double sw; //!< spectral width
extern double proton_freq; //!< proton Larmor frequency (0 if unset)
extern double dphs; //!< receiver phase shift
extern size_t ndims; //!< number of dimensions in data set
extern size_t nacqdims; //!< number of acquisition dimensions
extern size_t nobs; //!< observations per acquisition period (0 if no sync)
extern HamiltonianStore<space_T>* interactions_MFp; //!< NMR interaction set
extern spin_system* sysp; //!< spin system definition (NMRSIM_NULL if unset)
extern MasterObj* masterobjp; //!< master control object (NMRSIM_NULL if unset)
extern setableoperator_spec* sigma0_specp; //!< sigma0 definition (NMRSIM_NULL if unset)
extern setableoperator_spec* detect_specp; //!< detection operator def (NMRSIM_NULL if unset, sign reflects gamma)
extern BlockedOperator sigma0; //!< current sigma0
extern BlockedOperator detect; //!< current detection operator
extern BlockedOperator density; //!< current density matrix
extern Euler global_powder; //!< current orientation (alpha, beta only)
extern double current_gamma; //!< current gamma angle (-1 if unset)
//extern bool ishomogeneous; //!< \c true is Hamiltonian is homogeneous
extern bool haveactions; //!< \c true if \c pulseq is non-trivial (involves RF)
//extern bool haveduration; //!< \c true if \c pulseq involves non-ideal pulses
extern int nzcw; //!< powder orientation (ZCW) factor
extern bool zcwisindex; //!< \c true if ZCW value should be always be interpreted as index
extern bool havepowderquality(); //!< \c true if powder quality can be varied
bool isfrequency(); //!< \c true if data is frequency domain
//extern FILE* transitions_fp; //!< file pointer for transition list (NMRSIM_NULL if none)
bool transitions_log_active(); //!< \c true if transitions log file is active
bool logfile_active(); //!< \c true if log file is active
extern rmatrix vararray; //!< variable values per row of data set
extern rmatrix valerrs; //!< values and errors matrix
extern LIST<const char*> varnames; //!< associated names of parameters
void save_parameters(); //!< update current parameter values to ::vararray
extern CrystalStructure* cstructp; //!< periodic spin system specification
extern double master_time; //!< evolution time
size_t getarrayindex(size_t n); //!< get current index for dimension \a n
extern int global_prop_flags; //!< flags used to create propagator generators
void post_channels(); //!< clean up after \c channels directive
extern bool need_spinsys; //!< \c true if program requires spin_system to be set
extern double tolerance; //!< timing tolerance
extern double phasetolerance; //!< phase tolerance;
extern int gamma_angles; //!< number of gamma integration steps
extern int np; //!< current number of points in direct dimension
extern double detect_freq; //!< current detection Larmor frequency
double get_gamma1H(); //!< 1H magnetogyric ratio
extern bool havescale; //!< \c true if a scale factor been set
extern bool havestore; //!< \c true if at least one pulse sequence has been defined
extern bool powder_array; //!< \c if powder angle is being arrayed
extern LIST<nuclei_spec> blockingnuclei; //!< nuclei with mz blocking
bool have_spinsys(); //!< \c true if spin system has been defined
void verify_powderaverage();

extern double maxdtv; //!< integration timestep
enum tstep_t {
  INT_UNSET =0, //!< time step unset
  INT_EXPLICIT, //!< time step specified explicitly
  INT_ROTORPERIOD //!< time step specified as fraction of rotor period
};
extern tstep_t tstep_type; //!< specification of integration timestep
void set_inttype(tstep_t); //!< ensure integration time step (error if previously set differently)
extern int steps_per_rotation; //!< (minimum) integration timesteps per rotor cycle
double get_maxdt(); //!< return integration timestep
extern double maxjumpdt; //!< maximum timestep jump (experimental)

//! base class for "system variables"
struct SystemVariableBase
{
  SystemVariableBase(const std::string& namev) : name_(namev) {}
  virtual ~SystemVariableBase() {}
  
  //! return variable value (invalid for some types e.g. string)
  virtual double get() const
  { throw Failed("SystemVariableBase: can't get this type of variable"); }

  virtual const char* format() =0; //!< returns pointer to printable value
  virtual void print(std::ostream& =std::cout) const =0; //!< streams variable value
  void printvariablename(std::ostream& ostr) const { ostr << name_; }
  const char* name() const { return name_.c_str(); }
  
  //! variable name
  /** \note use of \c const means SystemVariables can't be copied (generally the Right Thing) */
  const std::string name_;
  mutable char buf[32]; //!< buffer for value (NB fixed size)
};

inline std::ostream& operator<< (std::ostream& ostr, const SystemVariableBase& a)
{
  a.print(ostr);
  return ostr;
}

template<class T> struct SystemVariable;

//! Specialisation of ::SystemVariable for pointers
template<class T> struct SystemVariable<T*> : public SystemVariableBase
{
  T* value; //!< pointer to value

  SystemVariable(const std::string& namev, T* value_)
    : SystemVariableBase(namev), value(value_) {}

  const char* format();

  void print(std::ostream& ostr, subsid_t subsid =S_NONE) const {
    if (subsid==S_NONE)
      ostr << this->name << " = ";
    ostr << (*value);
  }

  const T& operator()() const { return *value; }
};

enum { 
  V_ISFIXED=1, //!< flag no set method
  V_UPDATEINTS=2, //!< flag Hamiltonian needs updating if changes
  V_ISCONST=4, //!< flag evaluates to const (set by parse_systemvariable)
  V_POWDERSUMDEP=8 //!< flag depends on powder orientation / sum array
};

//! Specialisation of ::SystemVariable<T*> for \c double 
template<> class SystemVariable<double*> : 
  public SystemVariableBase, 
  public Setable,
  public RealVariable
{
public:

  SystemVariable(const std::string& namev, double* value_, double scalef_ =1.0, int flagsv =0);

  SystemVariable<double*>& operator=(double v) 
  { *value=v/scalef; return *this; } //!< assign variable

  void set(double, subsid_t =S_ARG1);
  void printvariablename(std::ostream& ostr, subsid_t) const { SystemVariableBase::printvariablename(ostr); }
  virtual void update() {} //don't check const as variables need to set post-construction

  double get() const { return *value; } //!< return (internal) scalar value
  const BaseList<double> get_list() const { validate(); return BaseList<double>(1,value); } //!< return (internal) value as list

  //! format content for external output
  const char* format() {
    snprintf(buf,sizeof(buf),"%g",scalef*get()); //output is scaled
    return buf;
  }
  void print(std::ostream& ostr) const;

  double operator()() const { return *value; }

private:
  double* value; //!< pointer to value
  double scalef; //!< scale factor between internal and external representation
  bool updateints; //!< \c true if changing value changes Hamiltonian

  const char* getname(std::string& dest, subsid_t) const { return SystemVariableBase::name_.c_str(); }
};

//! Specialisation of ::SystemVariable<T*> for \c int
template<> class SystemVariable<int*> : 
  public SystemVariableBase,
  public Setable,
  public RealVariable
{
public:

  SystemVariable(const std::string& namev, int* value_, int flagsv);

  SystemVariable<int*>& operator=(int v) 
  { *value=v; return *this; } //!< assign variable

  void set(double, subsid_t =S_ARG1);
  virtual void update() {} 
  double get() const { return double(*value); } //!< return scalar
  const BaseList<double> get_list() const { 
    validate();
    fvalue=*value;  
    return BaseList<double>(1,&fvalue); } //!< return as list

  const char* format() {
    snprintf(buf,sizeof(buf),"%i",*value);
    return buf;
  }

  //void printvariable(std::ostream& ostr, subsid_t) const;
  void print(std::ostream& ostr) const;
  void printvariablename(std::ostream& ostr, subsid_t) const { ostr << SystemVariableBase::name_; }
  //  void print(std::ostream& ostr, bool =true) const;

  int operator()() const { return *value; }

private:
  int* value; //!< pointer to value
  mutable double fvalue; //!< value stored as \c double for ::get_list
};

//! Specialisation of ::SystemVariable<T*> for \c size_t
template<> class SystemVariable<size_t*> : 
  public SystemVariableBase,
  public Setable,
  public RealVariable
{
public:
  SystemVariable(const std::string& namev, size_t* value_, int flagsv);

  SystemVariable<size_t*>& operator=(size_t v) 
  { *value=v; return *this; } //!< assign variable

  void set(double, subsid_t =S_ARG1);
  virtual void update() {} 
  double get() const { return double(*value); } //!< return scalar
  const BaseList<double> get_list() const { 
    validate();
    fvalue=*value;  
    return BaseList<double>(1,&fvalue); } //!< return as list

  const char* format() {
    snprintf(buf,sizeof(buf),"%lu",(long unsigned)(*value));
    return buf;
  }

  void printvariablename(std::ostream& ostr, subsid_t) const { SystemVariableBase::printvariablename(ostr); }
  void print(std::ostream& ostr) const;

  size_t operator()() const { return *value; }

private:
  size_t* value; //!< pointer to value
  mutable double fvalue; //!< value stored as \c double for ::get_list
};

//! Specialisation of ::SystemVariable<T> for \c string
template<> class SystemVariable<std::string> : public SystemVariableBase
{
public:
  SystemVariable(const std::string& namev, const std::string& value_)
    : SystemVariableBase(namev), value(value_) {}

  const char* format() { return value.c_str(); }
  void print(std::ostream& ostr) const;

  void set(const char* val) { value=val; }
  void set(const std::string& val) { value=val; }

  const std::string& operator()() const { return value; }

private:
  std::string value;
};

typedef double (*FUNC_REAL)(); //!< function returning \c double

//! Specialisation of ::SystemVariable<T> for function
template<> struct SystemVariable<FUNC_REAL> :
  public SystemVariableBase, 
  public RealVariable
{
public:
  SystemVariable(const std::string& namev, FUNC_REAL value_, double scalef_ =1.0)
    : SystemVariableBase(namev), 
      RealVariable(namev.c_str(),false,false),
      value(value_), scalef(scalef_) {}

  double get() const { return (*value)(); }
  const BaseList<double> get_list() const { 
    validate();
    fvalue=get();
    return BaseList<double>(1,&fvalue); }

  const char* format();
  //  void printvariable(std::ostream& ostr, subsid_t) const;
  void print(std::ostream& ostr) const;
private:
  FUNC_REAL value;
  double scalef;
  mutable double fvalue; //!< value stored as \c double for ::get_list
};

typedef MAPTYPE(SystemVariableBase*) systemvarmap_type;
extern systemvarmap_type systemvarmap; //!< map store for \c SystemVariables
void add_systemvarmap(SystemVariableBase&); //!< add \c SystemVariable to map

extern SystemVariable<double*> v_proton_freq,v_sw,v_sw1;
extern SystemVariable<int*> v_np,v_ni; //!< ::SystemVariable for number of data points
//extern SystemVariable<int*> v_evaluation_index;
bool isallowed_nonconst_var(SystemVariableBase*);
//extern SystemVariable<double*> v_synctime;
template<class T> void parse_system_variable(SystemVariable<T*>&, int flags =0); //!< set ::SystemVariable<T*> from input
template<class T> void parse_system_variable_syntax(const char*, size_t, SystemVariable<T*>&, int flags =0); //!< set ::SystemVariable<T*> from input
//void parse_coherencematrix(const char*); //!< parse coherence matrix filter specification
void parse_totalcoherencematrix(const char*); //!< parse coherence matrix filter specification

typedef std::pair<Setable*,subsid_t> varkey_t; //!< ::Setable and subsid indicator specifying a quantity to be varied
inline std::ostream& operator<< (std::ostream& ostr, const varkey_t& key) { return ostr << '(' << key.first << ',' << key.second << ')'; }
typedef std::map<varkey_t,VariableBase*> varmap_t; //!< map variable specification to ::VariableBase containing data specification

//! lightweight class used to relate variable values with associated quantity
struct Variable {
  Setable* ptr; //!< object that variable is associated with
  subsid_t subsid; //!< sub-type - which parameter is being varied
  VariableBase* valuep; //!< object holding variable values
  mutable char buf[20]; //!< buffer for external representation
  mutable std::string namestr; //!< buffer for name (not created by default)

  Variable(subsid_t subsidv =S_NONE, Setable* ptrv =NMRSIM_NULL, VariableBase* valuepv =NMRSIM_NULL);

  const char* name() const; //!< create and return pointer to name

  void update(bool allowwarnings =true); //!< update variable value
  void update(size_t, const BaseList<size_t>&); //!< update variable based on index into array
  void validate() const;  //!< throw ::InternalError exception if pointers not set

  double get() const; //!< get variable value
  void setptr(Setable*); //!< \internal set quantity pointer
  varkey_t key() const; //!< return key suitable for registry

  VariableBase& variable();
  const VariableBase& variable() const; //!< return source values

  const char* format() const;
  void printname(std::ostream& =std::cout) const; //!< stream name of associated parameter
  //! stream variable value
  /** All source values output if \a full is \c true otherwise just current value */
  void print(std::ostream&, bool =true) const;
  void rawupdate(); //!< \internal raw update
};

inline std::ostream& operator<< (std::ostream& ostr, const Variable& a)
{
  a.print(ostr);
  return ostr << '\n';
}

//! fitting variable
struct VarVariable : public Variable {
  VarVariable(const Variable& varv, Parameter& parv, int arraywhichv =-1);

  const char* ensure_named(); //!< create name if non exists
  void error(double); //!< set error value
  double get() const; //!< need to override Variable method
  void set(double); //!< set parameter
  void fix() { parameter.fix(); } //!< fix parameter value
  void release() { parameter.release(); } //!< release parameter
  bool isfixed() const { return parameter.isfixed(); } //!< \c true if parameter is fixed
  bool isconstrained() const { return parameter.isconstrained(); } //!< \c true if parameter value is constrained
  double error() const { return parameter.error(); } //!< \c return error/step on parameters
  double parse_error(); //!< parse error estimate (absolute value or %age)
  void unconstrain(); //!< release bound
  void constrain(const SimpleBoundsState&); //!< apply bound

  int arraywhich; //!< index into arrayed variable
  bool isnamed() const { return parameter.isnamed(); } //!< \c true if variable has name
  void printname(std::ostream& =std::cout) const; //!< stream name
  std::string& rawname() { return parameter.name(); } //!< return parameter name
  Parameter& parameter;
  SimpleBoundsState boundstate; //!< summary of bounds status
};

std::ostream& operator<< (std::ostream& ostr, const VarVariable&);

typedef LIST<Variable> varpars_t;
typedef LIST<VarVariable*> varvarpars_t;
extern varpars_t varpars; //!< list of variable quantities
extern varvarpars_t var_varpars; //!< list of fitting parameters
extern varpars_t expressions; //!< list of expressions
extern rmatrix covar; //!< covariance matrix
void register_fitting_variable(VarVariable*); //!< register new fitting variable

//! User defined variable
struct UserVariable : public RealVariable, public Setable {
public:

UserVariable(const char* namev, VariableBase& valuev, int flagsv =0)
  : RealVariable(namev,valuev), userflags_(flagsv) {}

  enum {
    IGNORE_UNUSED=1,
    IGNORE_REDEFINITION=2
  };
//   void reset(VariableBase& valuev, bool isconstv) //!< reset to new value / const status
//   { isconstant(isconstv); valuep=&valuev; }

  void print(std::ostream&) const; //!< output value (::Setable derivation)
  void printvariablename(std::ostream&, subsid_t) const; //!< output value (::Setable derivation)
  //  void print(std::ostream&, bool =true) const; //!< output value (::RealVariable derivation)
  void set(double, subsid_t =S_ARG1); //!< setter (single value)
  void set(const BaseList<double>&, subsid_t =S_ARG1); //!< setter (list value)

  static UserVariable* create(const char*, double, int =0); //!< create simple constant
  static UserVariable* create(const char*, VariableBase&, int =0); //!< create new variable \a n
  double get() const; //!< return scalar value
  int userflags() const { return userflags_; }
  //  bool isinternal() const { return (userflags_ && INTERNAL); }
  bool ignore_unused() const { return (userflags_ && IGNORE_UNUSED); }
  const BaseList<double> get_list() const { validate(); return value().value(); } //!< return value as list

  //  static ContextWarning<> replacingusedwithvar_warning; //!< replacing (used) variable with non-const
  static ContextWarning<> redefininguservariable_warning; //!< replacing (used) variable without -ignore_refinition
  static ContextWarning<> replacingunused_warning; //!< replacing unused (redudant?) variable
  static ThreadWarning<> unlinkedvar_warning; //!< non-const var unused
  static ThreadWarning<> unlinkedconstvar_warning; //!< const var unused
  static void check_unused(); //!< checked for unused variables

private:
  int userflags_;

  const char* getname(subsid_t) const { return name(); }
};
  
inline std::ostream& operator<< (std::ostream& ostr, const UserVariable& a) {
  a.print(ostr);
  return ostr;
}

//! dump map contents
template<class T> std::ostream& operator<< (std::ostream& ostr,const MAPTYPE(T*)& a)
{
  const typename MAPTYPE(T*)::const_iterator aend(a.end());
  typename MAPTYPE(T*)::const_iterator astart(a.begin());
  while (astart!=aend) {
    ostr << *(astart->second) << '\n';
    ++astart;
  }
  return ostr;
}

//!< dump map contents (::smartptr variant)
template<class T,bool isPoly> std::ostream& operator<< (std::ostream& ostr,const std::map<std::string, smartptr<T,isPoly> >& a)
{
  typedef typename std::map<std::string, smartptr<T,isPoly> >::const_iterator iterator;
  const iterator aend(a.end());
  iterator astart(a.begin());
  while (astart!=aend) {
    ostr << *(astart->second) << '\n';
    ++astart;
  }
  return ostr;
}

typedef void (*Parse_function)(); //!< parse function (no qualifier)
typedef void (*Parse_function_qual)(int); //!< parse function (with qualifier)

//! class for describing directives for parsing
struct par_t {
  Parse_function parsefunc; //!< pointer to parsing function (no qualifier)
  Parse_function_qual parsefunc_qual; //!< pointer to parsing function (with qualifier)
  bool allowmultiple; //!< directive may appear more than once if \c true
  int ident; //!< qualifier for this variant
  bool done; //!< \c true if directive has been used

  //! unqualified directive
  par_t(Parse_function func_ =NMRSIM_NULL, bool allowmultiple_ =false)
    : parsefunc(func_), parsefunc_qual(NMRSIM_NULL), allowmultiple(allowmultiple_), done(false) {}

  //! qualified directive
  par_t(Parse_function_qual func_, int ident_, bool allowmultiple_ =false)
    : parsefunc(NMRSIM_NULL), parsefunc_qual(func_), allowmultiple(allowmultiple_), ident(ident_), done(false) {}

  void operator()(); //!< parse directive
};

typedef std::map<std::string,par_t> command_Factory_t; //!< factory type for directives
extern command_Factory_t& get_par_Factory(); //!< get par block directive factory
extern command_Factory_t& get_Global_Factory(); //!< get global directive factory
bool havekey(const command_Factory_t&, const std::string&); //!< returns \c true if directive has been used

void update_expressions(); //!< update expressions (using current variable values)
extern bool update_interactions; //!< update Hamiltonian

void parse_channels();
void parse_start_operator();
void parse_detect_operator();
void parse_precision();
void parse_spin_rate();
void parse_spin_rate_jitter();
void parse_gamma_angles();
void parse_gamma_zero();
void parse_rotor_angle();
void parse_cache_limit();
void parse_variable();
void parse_function();
void parse_verbose();
void parse_tolerance();
void parse_makefilter();
void parse_matrix();
void parse_log_file();
void parse_transients();
void parse_time_resolution();

int check_sync(double nfloat,double tol =1e-6); //!< convert floating point ratio to integer or 0 if outside tolerance \a tol

productoperator_spec* create_productoperator(char*, int flags =0); //!< parse productoperator from input
setableoperator_spec* create_setableoperator(char*, int flags =0); //!< parse ::setableoperator_spec from input
productoperator_spec* parse_productoperator(int flags =0);
setableoperator_spec* parse_setableoperator(int flags =0);

void make_data_variables(); //!< minimal set of par variables for defining data set
void make_1d_par_variables(); //!< create minimal set of par variables for 1D simulation (includes data_variables)
void make_par_variables(); //!< initialise \c SystemVariables that may be used in \c par (include 1d_par_variables)
void make_pulseq_variables(); //!< initialise \c SystemVariables that may be used in \c pulseq

void write_matrix(FILE* fp, const Matrix<bool>&, const char*, int);
//  fputs("Matrix<bool> write not supported\n",fp);
//}

double get_rotor_phase(double t); //!< return rotor phase IN DEGREES for time \a t (0 if static)
size_t trysync(double& synctime, double period, const char* desc); //!< test if \a period divides into \a synctime (return 0 if not)
void init_parallel(int&, char**&); //!< initialise parallel computation
bool cleanup_parallel(); //!< cleanup parallel computation (after irregular exit prior to exit). Return \c true if no errors
bool isirregular(); //!< returns \c true if data set is irregular
//void ensure_powder(); 
void get_variables(BaseList<double> pars); //!< get current values of all variable parameters
UserVariable* findvariable(const char* name, bool flagused =true); //!< find user defined variable, automatically flag used if \a flagused set

//! stores reference to current file being parsed (needed for ::dimension_set)
struct parse_state {
  std::string fname; //!< name of file
  size_t curline; //!< line number (from 1)
  size_t lines; //!< number of lines (usually 1)

  parse_state() : curline(0) {}
  explicit parse_state(const char* fname_)
    : fname(fname_), curline(0) {}
};

//! virtual dimension set
struct dimension_set
{
  dimension_set() : block_change(false) {}
  dimension_set(size_t n); //!< simple 1D data set

  //! individual dimension
  struct dimension {
    dimension() : dim(0), hard(0) {}
    void set(size_t,bool); //!< update value
    void set(size_t,size_t); //!< set explicit value and skip
    size_t get(size_t& skip, size_t N) const; //!< get value and skip (dimension \a N)
    bool isset() const { return (hard || dim); } //!< return \c true if dimension has been set
    bool ishard() const { return (hard!=0); }
    static void error(size_t N,const char*); //!< incommensurate error for dimension \a N

    size_t dim; //!< current dimension size
    size_t hard; //!< explicitly set size
    size_t skip; //!< skip index
    parse_state state; //!< reference to last change
  };

  dimension dims[MAX_DIMENSIONS+1]; //!< array of dimensions
  bool block_change;
  
  void set(size_t N,size_t n); //!< set dimension \a N to size \a n
  void set(size_t N,size_t n,size_t s); //!< set size \a n, skip \a s of dimension \a N
  size_t get(LIST<size_t>& ns, LIST<size_t>& ss, size_t& tot) const; //!< return sizes \a ns, skips \a ss and total rows \a tot
  void lock() { block_change=true; } //!< block further changes

  struct incommensurate_warning_t : public ContextWarning<> {
    incommensurate_warning_t(const char*);
    void raise(const parse_state&);
  };
  static incommensurate_warning_t incommensurate_warning;
  static incommensurate_warning_t incommensurate_warning2;
};

void update_auxiliary_vars(bool considercurrent = true);
double gammatimeoffset(); //!< express current powder gamma angle as time offset in Hsys
extern dimension_set array_dims; //!< dimension set for array
extern dimension_set sum_dims; //!< dimension set for sum
extern size_t array_n0; //!< total rows in data set
//extern size_t nD_n0; //!< total rows in indirect dimensions;
extern size_t sum_n0; //!< total 'rows' in sum
extern LIST<size_t> array_ns; //!< sizes for array
extern LIST<size_t> sum_ns; //!< sizes for sum
extern LIST<size_t> skips; //!< rows-per-time-increment for array
extern LIST<int> ns; //!< indirect dimension sizes
//! rows-per-time-increment for sum
/** \todo Do >1 skips make sense for summed variables? */
extern LIST<size_t> array_skips; 
extern LIST<size_t> arrayinds; //!< current set of array indices
extern LIST<double> sws; //!< spectral widths
//extern LIST<double> sfrqs; //!< spectrometer frequencies
extern bool process2D; //!< if \c false then 2D processing is impossible
extern simple_counter global_cache; //!< track overall memory usage for cached propagators
extern int global_argc; //!< original \a argc (after MPI init)
extern char** global_argv; //!< original \a argv (after MPI init)
extern bool randomise; //!< if \c true then any random number generators should be randomly seeded prior to use
extern bool ammaster; //!< \c true if master control "thread"

//! base class for \c echo directives
class BaseEcho
{
public:
  BaseEcho(const char* str_ =NMRSIM_NULL) : str(str_ ? str_ : "") {}
  void exec() const;
  void print(std::ostream& ostr) const
  { ostr << "echo " << str.c_str() << '\n'; }
  int update() const; //!< update buffer and return accumulated 'uses' of $variables
private:
  std::string str; //!< output
  mutable char buffer[256]; //!< buffer
};

//! iterator through dimension set
struct array_iter {
  array_iter(bool updatevars =true);
  ~array_iter() { 
    var_index=-1; 
    row_index=-1;
  } //!< clears global row indices to trap misuse

  bool next(const DataStore* =NMRSIM_NULL); //next 'row', returns \c true if variables need updating, optional pass data set for setting current row

  size_t size() const; //!< number of steps in iterator

  bool updatevars_; //!< if \c true, variables/expressions will be updated on change
  bool finished; //!< \c true when iterator finished
  MultiIterator iter; //!< raw iterator
};

//! combine list sizes
class accumulator {
public:
  accumulator() : length_(-1), failed_(false) {}
  size_t operator()() const; //!< return final value
  bool operator!() const 
    { return (length_ < 0); }  //!< return \c true if undefined (May 2022)
  bool add(size_t n); //!< add length \a n returning \c false if incompatible
private:
  int length_; //!< current length
  bool failed_; //!< \c true if incompatible sizes detected
};

//! stream container
template<class T> void dump(const T& a, const char* term ="")
{
  const typename T::const_iterator end(a.end());
  typename T::const_iterator start(a.begin());
  while (start!=end) {
    std::cout << (*start) << term;
    ++start;
  }
}

//! stream container of pointers
template<class T> void dumpptr(const T& a, const char* term ="")
{
  const typename T::const_iterator end(a.end());
  typename T::const_iterator start(a.begin());
  while (start!=end) {
    std::cout << (**start) << term;
    ++start;
  }
}

//! stream list of pointers
template<class T> std::ostream& operator<< (std::ostream& ostr,const BaseList<T*>& a)
{
  const size_t n=a.size();
  if (n) {
    for (size_t i=0;i<n;i++)
      ostr << *(a(i)) << '\n';
  }
  else
    ostr << "(empty)\n";
  return ostr;
}

template<class T> inline std::ostream& operator<< (std::ostream& ostr,const LIST<T*>& a)
{ return ostr << static_cast< const BaseList<T*>& >(a); }

//! add object to container and return reference to stored object
template<class Store,class T> T& push_and_get(Store& store, const T& t) { 
  store.push_back(t); return store.back();
}

//! string representation of floating point number \a val
void snprintf_prec(char* buf, size_t s, double val, bool ignoreformat =false);

//! default partioning (no RF) - don't do anything for now
template<class T> void set_partitioning(T&) {}

//! set partitioning appropriate for sequence \a seq
template<class T> void set_partitioning(T& propgen, const CompSequenceBase&);

//! pretty print start-end interval
std::ostream& prettyprint_interval(double, double, std::ostream& =std::cout);

extern double zerotolerance; //!< tolerance for rounding small values to zero prior to diagonalisation

typedef void (*simple_callback_t)();
void addprepar_callback(simple_callback_t); //!< add callback prior to par block
struct SimpleCallbackStack {
  LIST<simple_callback_t> stack;
  void push(simple_callback_t callback) { stack.push_back(callback); }
  void exec() const;
};

void flush_calculation(); //!< close any open files, tidy up
void add_flushcallback(simple_callback_t); //!< add flush callback for end of calculation 
//bool hasacquisition_rf(); //!< \c true if RF pulse is present during acquisition
//const rmatrix& make_powderarray();
extern double gamma_zero;
extern int orientation; //!< current orientation index
extern int eigenvalue; //!< symmetry eigenvalue (INVALID if unset)
extern int sideband; //!< sideband selection (INVALID if unset)
extern double def_rotor_angle; //!< default rotor angle
void ensure_sigma0detect(); //!< ensure sigma0/detect are up to date
void ensure_operator_matrices(bool constonly =false); //!< ensure operator matrices (including sigma0/detect) are built
extern SystemVariable<double*> v_gamma; //!< system variable for gamma rotor angle

//void prepar_initialise(); //!< setup prior to par block
void read_par_block(bool optional =false); //!< read par block

template<class M> void spyprint(std::ostream&, const BlockedMatrix<M>&, int flags =0);
void spyprint(std::ostream&, const BlockedOperator&, int flags =0); //!< dump operator matrix to stream (structure if precision is 0)

#include "Histogram.h"
#include "Lineshapes.h"

//! scratch object for creation and tidy-up of histogram using current settings
class HistogramMaker {
public:
  HistogramMaker(BaseList<complex>, double);
  ~HistogramMaker();
  BaseHistogram<complex>& operator()() { return *histp_; } //!< return reference to encapsulated Histogram
  const BaseHistogram<complex>& operator()() const { return *histp_; }

  static ThreadWarning<> loss_warning; //!< warning over histogram loss (only if verbose)
  static ThreadWarning<> foldrange_warning; //!< -fold combined with histogram range

private:
  smartptr< BaseHistogram<complex> > histp_; 
};

//! histogram flags
enum { HIST_BINARY=1, //!< binary output
       HIST_INTERP=2, //!< interpolate histogram
       HIST_DOUBLE=4, //!< save in double precision
       HIST_FOLD=8, //!< fold frequencies into spectral width
       HIST_REAL=16, //!< output only real amplitude
       HIST_RANGE=32,  //!< restrict histogram width
       HIST_APPEND=64 //!< append to log_file
};
extern int hist_flags; //!< accumulated histogram flags
extern double hist_threshold; //!< threshold for transition output
extern double hist_min; //!< histogram minimum (if \c HISTO_REAL active)
extern double hist_max; //!< histogram maximum (if \c HISTO_REAL active)
extern double hist_lw;
extern double hist_gfrac;
command_Factory_t& get_histogram_Factory(); //!< return registry for histogram directives
void update_lineshapegen();
extern LIST<size_t> losttrans; //!< total 'lost' transitions
extern LIST<complex> lostsum; //!< total 'lost
extern LineshapeSpec lineshape_spec;
extern LorentzianGaussian* lineshape_genp;

#define NMRSIM_STRINGISE(x) #x

extern ThreadWarning<> commonPAS_warning; //!< common PAS optimisation not implemented
extern ThreadWarning<> mzsymmetry_warning; //!< mzsymmetry probably enabled inappropriately
extern ThreadWarning<> missingacq_warning; //!< no direct acquisition
extern ThreadWarning<> nothingtodo_warning; //!< nothing to simulate!
extern ThreadWarning<> nointeractions_butMAS_warning; //!< no interactions specified but still under MAS
extern ThreadWarning<> generalisedqpole_warning; //!< warnings related to generalised quadrupole treatment
extern ThreadWarning<> generalisedqpole_rf_warning; 
extern ContextWarning<> zerotolerance_warning; //!< zerotolerance looks too large
extern ContextWarning<> chebyshev_iterations_warning; //!< chebyshev_iterations looks too small
extern ContextWarning<> maxdt_warning; //!< maxdt looks too large
extern ContextWarning<> zeroangle_warning; //!< zero rotor angle
extern ContextWarning<> doubleascii_warning; //!< -double flag has no effect on ASCII output
extern ContextWarning<> streamopenfailed_warning; //!< failed to open log file
extern ContextWarning<> excessivecutoff_warning; //!< suspiciously large histogram lineshape cutoff
extern ThreadWarning<> ignoringsideband_warning; //!< ignoring sideband restriction in static simulation

void dirty_stack_push(Variable*); //!< add to stack of ::Variable that refer to an object in course of creation
extern size_t global_np;
extern LIST<size_t> global_nps; //!< global data set size (only fixed after all processing)
extern LIST<Parameter> parameter_list; //!< list of fitting variables
extern bool update_vars; //!< \c true if variables need updating (at change of row)
//extern LIST<BaseWarning*> resetwarning_list; //!< list of serious warnings that should be reset for each calculation
void add_reset_warning(BaseWarning&); //!< add warning to reset list
void reset_warnings();
void set_n(int, int, int =1); //!< set indirect dimension parameters
void raw_set_n(int, int, int =1); //!< set indirect dimension parameters
void ensure_dimension(size_t);
void ensure_array(); //!< ensure looping constructs enabled
void parse_swn(int i); //!< parse sw<i>
void add_channel(size_t nuc); //!< add channel for nucleus nuc

bool checkoptflag(const char* com, bool& havedisable, bool ignoreunknown =false, bool allowtest =true);
bool checkverboseflag(const char*); //!< check for -verbose flag, \a true if found
bool checkflags(bool& fl, const char* com, const char* name);
void setwarnings(BaseWarning::warning_t); //!< global override of warnings
extern ThreadWarning<> syncfailed_warning;
extern ThreadWarning<> library_version; //mismatch between libcmatrix versions
extern Matrix<accumulating_timer<> >* timers_arrayp; //!< array of timers for profiling individual row/summation indices (NMRSIM_NULL if not being used)
void dump_arrayprofile(const Matrix<accumulating_timer<> >&, double =1.0, const char* ="s"); //!< dump profile
void dump_timings(bool abbrev =false); //!< dump timings profile (abbreviated if \a abbrev set)
void check_powderarray(); //!< set up powder array (if any)
void checklogfile(const char* fname, bool otheractive, bool append); //!< crude check for overlapping log files

typedef FASTMAPTYPE(const char*) spell_dictionary_type;
extern spell_dictionary_type& get_spell_dictionary(); //!< return dictionary for possible spelling errors
bool check_unrecognised(const char* unknown, const char* name =NMRSIM_NULL); //!< check unrecognised quantity \a unknown of type \a name (default -> "directive") for mis-spelling
const char* register_and_getenv(const char*); //!< register name of environment variable (to spot later changes) and get value
bool isregistered(const char*); //!< return true if environment variable name already used
bool parse_optimise_block(); //!< return true if read

extern IntervalSampler jitterinfo; //!< MAS instability info
void testoptions(MasterObj&);

void update_array_loop(size_t nr);
HamiltonianStore<space_T>& get_Hamiltonian();

#endif
