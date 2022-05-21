#ifndef NMRSIM_ACTION_H_
#define NMRSIM_ACTION_H_

/*! \file
 \brief  Pulse sequence directives
*/

#include "NMRsim_RF.h"
#include "Parser.h"
#include "timer.h"  

#ifndef DISABLE_EXEC
#define DECLARE_EXEC void exec(MasterObj*, BlockedOperator&, double&) const;
#define DECLARE_FLUSH_CACHE  void flush_cache(dirty_t =DIRTY_ALL);
#define DECLARE_RESET void reset();
#else
#define DECLARE_EXEC
#define DECLARE_FLUSH_CACHE
#define DECLARE_RESET
#endif

//! (abstract) base class for pulse sequence elements
class ActionCommand : public Setable  {
public:
  ActionCommand(bool beverbose_ =false) : beverbose(beverbose_) {}
  virtual ~ActionCommand() {};

  virtual void printvariablename(std::ostream&, subsid_t) const { throw InternalError("Type does not have Setable parameters"); }
  double time() const { return stopwatch(); } //!< return accummulated calculation time
  static ActionCommand* create() { throw Failed("ActionCommand() should not be called"); } //!< parse new directive
  virtual void print(std::ostream&) const =0;
  virtual double duration() const { return 0.0; } //<! duration of element (fallback is zero duration)
  virtual void reset() {} //!< restart calculation
  virtual void flush_cache(dirty_t) {} //!< clear any propagator caches
  virtual void set(double,subsid_t)
  { throw InternalError("This function is not setable"); } //!< set parameter value

  //! not abstract for "translate" programs that don't define exec's
  virtual void exec(MasterObj*, BlockedOperator&, double&) const
  { throw InternalError("This should never be called\n"); }
  virtual usage_t usage() const { return usage_t(); } //!< return memory usage

  virtual void isconstant(bool) {} //!< by default ignore change of const
  friend class actionstack_t;

protected:
  mutable accumulating_timer<> stopwatch; //!< accumulates total calculation time in directive
  bool beverbose; //!< verbose flag
};

inline std::ostream& operator<< (std::ostream& ostr,const ActionCommand& a)
{
  a.print(ostr);
  return ostr;
}

//! stack for pulse sequence actions 
struct actionstack_t : public LIST<ActionCommand*> {
  actionstack_t();
  // void parse(char*); //!< parse new directive
  void reset(); //!< new calculation reset
  void flush_cache(dirty_t); //!< flush propagator caches
#ifndef DISABLE_EXEC
  void exec(MasterObj*, BlockedOperator& rho ,double& t); //!< apply processing to density matrix \a rho at time \a t
  static void execaction(ActionCommand* actionp, MasterObj* objp, BlockedOperator& dens, double& t, size_t& count);
#endif
  void initialise();
  size_t acq_count; //!< count of number of \c acq dimensions
};


//! base class for control classes 
class ActionControl : public ActionCommand
{
 public:
  ActionControl(const char* typev) : type_(typev), items_(0) {}
  virtual ~ActionControl() {}

  const char* type() const { return type_; }
  void items(size_t sizev) { items_=sizev; }
  size_t items() const { return items_; }
  virtual void exec_control(actionstack_t::const_iterator, const actionstack_t::const_iterator&, size_t&, MasterObj* objp, BlockedOperator& dens,double& t) const =0;

 private:
  const char* type_;
  size_t items_;
};

//! do loop
class ActionDo : public ActionControl
{
 public:
  ActionDo(size_t timesv) : ActionControl("do"), times_(timesv) {}
    void printvariablename(std::ostream& ostr, subsid_t) const { ostr << "do"; }
    void print(std::ostream& ostr) const { ostr << "do " << times_; }
    void set(double,subsid_t);
    static ActionCommand* create();
    void exec_control(actionstack_t::const_iterator,const actionstack_t::const_iterator&, size_t& count, MasterObj* objp, BlockedOperator& dens,double& t) const;

    static ContextWarning<> zeroloops_warning; //!< no iterations!
 private:
    size_t times_;
};

//! \c echo directive in \c pulseq block
struct ActionEcho : public ActionCommand, public BaseEcho  {
  ActionEcho(const char* str_ =NMRSIM_NULL) : BaseEcho(str_) {}
  DECLARE_EXEC
  void print(std::ostream& ostr) const { BaseEcho::print(ostr); }
  static ActionCommand* create(); //!< parse input to create new ::ActionEcho object
};

//! \c put \c matrix directive
class ActionPutMatrix : public ActionCommand
{
public:
  ActionPutMatrix(const matrixspec& curm_, const char* name_, int flagsv =0, statistics_t statsselv =statistics_t(0, (const filter_def*)NMRSIM_NULL) )
    : name(name_), curm(curm_), flags_(flagsv), statssel(statsselv), done_(false) {}
  
  DECLARE_EXEC
  
  void print(std::ostream&) const;

  static ActionCommand* create(); //!< parse input to create new ::ActionPutMatrix object

private:
  std::string name; //!< matrix name
  matrixspec curm; //!< reference to matrix
  int flags_;
  statistics_t statssel; //!< argument for statistics output
  mutable bool done_; //!< flag if done
};

//! \c prop directive
  class ActionProp : public ActionCommand
{
public:
  ActionProp(CycledSequence& seq, int n, int flags); //!< \c prop \a seq \a n \a flags
  ActionProp(CycledSequence& seq, double tau, int flags); //!< \c propfor \a seq \a tau \a flags
  ActionProp(double tau, int flags); //!< \c propfor without sequence

  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  void set(double,subsid_t); //!< set \a n / \a propfor

  //  usage_t usage() const { return Ucache.usage(); }
  double duration() const; //!< return duration of one step
  bool iszeroduration() const; //!< return true if exec if for zero time

  static ActionCommand* create() { return create_(false); }  //!< parse input to create new ::ActionProp object (normal \c prop)
  static ActionCommand* createfor() { return create_(true); }  //!< parse input to create \c propfor variant of ::ActionProp

  DECLARE_FLUSH_CACHE
  DECLARE_EXEC
  DECLARE_RESET

  static ContextWarning<> propafteracq_warning; //!< prop after final acq
  static ContextWarning<> cyclingmismatch_warning; //!<  propagator cycling doesn't match ni step
  static ThreadWarning<> empty_propagator_warning; //!< propagator does nothing
  static ThreadWarning<> disablingcombine_warning; //!< not using propagator combination
  static ContextWarning<> invalidreset_warning; //!< reset not valid without sequence
  static ContextWarning<> nocache_warning; //!< warning that cacheing not possible
  static ContextWarning<> inv_warning; //!< warning that -inv is sub-optimal
  bool match(const ActionProp&) const;

private:
  enum { PROP_INV=1, //!< \c -inverse flag
	 PROP_PUTMATRIX=2, //!< \c -putmatrix flag
	 //	 PROP_FULL=4 //!< -full flag
	 PROP_RESET=8 //!< \c -reset flag
  };
  CycledSequence* seqp; //!< pointer to sequence (NMRSIM_NULL if none)
  double propfor; //!< time to propagate for (contents valid for \c propfor only)
  int ntimes; //!< times to apply (-1 indicates \c propfor form)
  int flags_;
  bool checkedsync; //!< flag whether sync has been checked

  void set_ntimes(int n); //!< \internal set \a n
  void set_propfor(double); //!< \internal set \a propfor
  static ActionCommand* create_(bool); //!< \internal shared create
  static const flagsmap_type& get_flags();
};

//! get coherence object
class ActionGet : public ActionCommand
{
public:
  ActionGet(UserVariable& destv, const char* namev, const BlockedOperator& opv, const char* fromnamev, matrixspec fromv);
  ActionGet(UserVariable& destv, const char* namev, const BlockedOperator& opv); //!< get from density matrix
  static ActionCommand* create();
  void print(std::ostream&) const;

  DECLARE_EXEC
  
  static ThreadWarning<> ignoring_imaginary_warning; //!< ignoring imaginary component

private:
  UserVariable& dest_; //!< destination variable
  std::string name_; //!< name of coherence matrix
  const BlockedOperator& op_; //!< associated operator matrix
  std::string fromname_; //!< name of source (NMRSIM_NULL for current density matrix)
  matrixspec from_; //!< reference to source matrix
  bool isdensity() const { return (*fromname_.c_str()=='\0'); } //!< \c true if source is current density matrix

  static ContextWarning<> notimpl_warning;
};

//! set magnetisation object (currently PROBLEMATIC)
// class ActionSet : public ActionCommand
// {
// public:
//   ActionSet(double scalev, const char* namev, const BlockedOperator& opv);
//   ActionSet(double scalev, const char* namev, const BlockedOperator& opv, const char* tonamev, BlockedOperator&);
//   static ActionCommand* create();
//   void print(std::ostream&, subsid_t =S_NONE) const;
//   void set(double,subsid_t);

//   DECLARE_EXEC

// private:
//   double scale_;
//   std::string name_;
//   const BlockedOperator& op_;
//   std::string toname_;
//   BlockedOperator* top_;
// };

struct named_operator {
  std::string name; //!< name of operator
  const BlockedOperator& op; //!< reference to operator matrix
  named_operator(const char* namev, const BlockedOperator& opv)
    : name(namev), op(opv) {}
  named_operator(char*);
};

//! Transfer coherence object
class ActionTransfer : public ActionCommand
{
public:
  ActionTransfer(const named_operator&, const named_operator&);
  static ActionCommand* create();
  void print(std::ostream&) const;

  DECLARE_EXEC

private:
  const named_operator from_; //!< source coherence matrix
  const named_operator to_; //!< destination coherence matrix
};

//! Exchange coherence object
class ActionExchange : public ActionCommand
{
public:
  ActionExchange(const BaseList<named_operator>&, const ExchangeMatrix&);
  static ActionCommand* create();
  void print(std::ostream&) const;

  DECLARE_EXEC

  static ContextWarning<> exchangematrix_warning;
private:
  mutable LIST<complex> amps0_,amps_;
  LIST<named_operator> ops_;
  const ExchangeMatrix& mat_;
};

//! \c scale directive
class ActionScale : public ActionCommand {
public:
  ActionScale(double scale_, const char* name_ =NMRSIM_NULL, const filter_def* filterp_ =NMRSIM_NULL)
    : scale(scale_), name(name_), filterp(filterp_) {}

  void set(double scale_, subsid_t =S_NONE) { scale=scale_; }

  DECLARE_EXEC
  
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;

  static ActionCommand* create();
  static ContextWarning<> coherencefilter_warning; //!< warn that scale is being used with coherence list filter

private:
  double scale; //!< scale factor
  const char* name; //!< name of (optional) filter matrix
  const filter_def* filterp; //!< optional filter matrix
};

//! \c acqpoint directive
class ActionAcqPoint : public ActionCommand {
public:
  ActionAcqPoint(double dphsv =0.0) : dphs(dphsv) {}
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;

  DECLARE_EXEC
  
  void set(double, subsid_t);
  static ActionCommand* create(); //!< parse to create new ::ActionAcqPoint object
  double recphase() const { return dphs*(M_PI/180.0); }
private:
  double dphs; //!< phase shift
};

//! direct dimension acquisition object
class ActionDirectAcq : public ActionCommand
{
public:
  ActionDirectAcq(double, CompSequenceBase* seqp_ =NMRSIM_NULL, double =0.0); 
  //  ActionDirectAcq(size_t, double, CompSequenceBase* seqp_ =NMRSIM_NULL, double =0.0); 
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  usage_t usage() const; //!< return memory usage

  DECLARE_EXEC
  DECLARE_FLUSH_CACHE
  //  DECLARE_RESET

  const CompSequenceBase* sequence() const { return seqp; } //!< return sequence pointer
  void checksync(); //!< check synchronisation of acquisition (called by ::MasterObj::restart)
  void set(double,subsid_t);
  bool ispoint_by_point() const { return point_by_point; } //!< getter for ::point_by_point
  static ActionDirectAcq* create(double, CompSequenceBase* =NMRSIM_NULL, double =0.0); //!< parse to create new ::ActionDirectAcq object

  static ThreadWarning<> zerodurationacq_warning; //!< sw specified with zero duration acquisition sequence
  static ThreadWarning<> frequencydomainRF_warning; //!< RF + frequency domain unstable?
  static ThreadWarning<> ignoringsynctime_warning; //!< ignoring synchronisation time
  static ThreadWarning<> gammasteps_warning; //!< gamma steps not a multiple of observations
  static ThreadWarning<> disabledgamma_warning; //!< gamma compute algorithm being disabled
  static ThreadWarning<> syncsetnotfound_warning; //!< synchronisation time set but not established
  static ThreadWarning<> syncfailed_warning; //!< failed to find synchronisation
  static ThreadWarning<> empty_warning; //!< empty acquisition

  static void reset_acquisition_count() { acquired_=0U; }
  static size_t acquisition_count() { return acquired_; }
  static void add_acquisition_count(size_t n) { acquired_+=n; }

  enum create_status_t { NONE, PARTIAL, FINAL };
  static create_status_t create_status;
  static bool have_acq() { return (create_status!=NONE); } //!< return \c true if direct dimension acquisition has been defined
  double recphase() const { return dphs*(M_PI/180.0); }
  static void add_acq(bool israw =false); //!< update create status

private:
  double dphs; //!< phase shift
  CompSequenceBase* seqp; //!< pointer to acquisition sequence (NMRSIM_NULLif none)
  bool point_by_point; //!< accumulating time-domain point by point?
  double synctime; //!< time for synchronisation between detection and MAS/RF (0 if unset)
  bool checkedsync; //!< \c true if synchronisation has been checked (only required once)
  //int toacquire_; //!< data points to acquire (-1 to complete rest of FID)

  static size_t acquired_; //!< data points accumulated  
  //static void finishacq(); //!< called when FID complete
};

//extern ActionDirectAcq* acqp; //!< pointer to (unique) direct dimension acquisition object (NMRSIM_NULLif none)

class ActionAcqN;

//! \c acq directive (indirect dimensions)
class ActionAcq : public ActionCommand {
public:
  ActionAcq(CompSequenceBase* seqp_, int =0);
  void store_transients(double t) const { ::store_transients(translist,t); } //! store transients after propagator evaluation (currently disabled)
  void restore_transients() const {
    if (!(translist.empty())) 
      ::restore_transients(translist); } //!< restore transient prior to propagator evaluation (currently disabled)

  void print(std::ostream&) const;
  static ActionCommand* create() { return create_(false); } //!< parse to create new ::ActionCommand
  //  static ActionCommand* createraw() { return create_(true); } //!< parse to create new ::ActionCommand

  DECLARE_FLUSH_CACHE
  DECLARE_EXEC
  DECLARE_RESET

  usage_t usage() const; //!< return memory usage
 const CompSequenceBase* sequence() const { return seqp; } //!< return sequence pointer

  static void transientcleanup(BlockedOperator&);
  double spectralwidth() const; //!< return spectral width for relevant dimension
  double dwelltime() const; //!< return dwell time
  size_t dimensionsize() const; //!< return dimension length

  static bool isdirectdimension(); //!< return \c true if direct dimension

  static ThreadWarning<> zerodurationacquisition_warning; //!< sw set with zero-duration acq sequence
  static ThreadWarning<> excessivesync_warning; //!< synchronisation time longer than acquisition duration
  static ThreadWarning<> activetransients_warning; //!< transients before acq
  static ThreadWarning<> trailingtransients_warning; //!< transients after acq
  static ThreadWarning<InternalError> improperreuse_warning;
  static ThreadWarning<> emptypropagator; //!< propagator doesn't do anything

  void increment(MasterObj*, double, int) const;
  void applyincrement(BlockedOperator&, bool =false) const;

  friend class ActionAcqN;
private:
  CompSequenceBase* seqp; //!< pointer to acquisition sequence (NMRSIM_NULLif none)
  //  double dt; //!< dwell time
  size_t nrep; //!< periodic synchronisation factor
  size_t ndim; //!< indirect dimension number (from 0)
  int flags_;

  mutable LIST< BlockedMatrix<complex> > Us; //!< step propagator cache
  mutable LIST<double> phases; //!< rotor phases corresponding to cached step propagators
  mutable int index; //!< current propagator index
  mutable LIST<transient_state> translist; //!< transient states
  mutable BlockedMatrix<complex> Uacc; //!< accumulated propagator

  static ActionCommand* create_(bool israw);
  size_t checksync(double&) const; //!< check synchronisation between RF, MAS, dt
  static const flagsmap_type& get_flags();
};

class ActionAcqN : public ActionCommand
{
public:
  ActionAcqN(size_t, double, CompSequenceBase*, int =0);
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  void set(double, subsid_t);
  static ActionCommand* create();

  DECLARE_FLUSH_CACHE
  DECLARE_EXEC
  DECLARE_RESET

  usage_t usage() const { return acqp->usage(); }

  static ThreadWarning<> fullFID_warning;
  
private:
  size_t toacquire_;
  smartptr<ActionAcq> acqp;
  smartptr<ActionAcqPoint> acqpointp;  
};

//! \c rotor_angle directive
class ActionRotorAngle : public ActionCommand {
public:
  ActionRotorAngle(double rotor_angle_)
    : rotor_angle(rotor_angle_) {}

  void print(std::ostream& ostr) const {
    ostr << "rotor_angle " << rotor_angle;
  }
  void set(double anglev, subsid_t) { rotor_angle=anglev; }
  DECLARE_EXEC

  static ActionCommand* create(); //!< parse to create new object

private:
  double rotor_angle; //!< rotor angle
};

//! \c apply_inversion directive

namespace libcmatrix {
  class InversionGenerator;
}

class ActionInversion : public ActionCommand {
public:
  ActionInversion(size_t nucv, const InversionGenerator& invgenv, double phasev)
    : nuc_(nucv), invgen_(invgenv), phase_(phasev) {}
  
  void print(std::ostream&) const;
  
  void set(double phasev, subsid_t) { phase_=phasev; }
  DECLARE_EXEC

  static ActionCommand* create();

private:
  size_t nuc_;
  const InversionGenerator& invgen_;
  double phase_; //!< phase to use
  mutable BlockedOperator sigmatmp_;
};


//! \c timeadjust directive
class ActionTimeAdjust : public ActionCommand {
public:
  ActionTimeAdjust(double timev, bool isabsv)
    : time_(timev), isabsolute_(isabsv) {}
  
  void print(std::ostream&) const;
  void set(double timev, subsid_t) { time_=timev; }
  DECLARE_EXEC

  static ActionCommand* create();
  static ThreadWarning<> back_warning; //!< moving time backwards warning
  static ContextWarning<> static_warning; //!< Hamiltonian is time independent
private:
  double time_; //!< time 
  bool isabsolute_;
};

//! \c filter directive
class ActionFilter : public ActionCommand {
public:
  ActionFilter(const char* name_,const filter_def& filter_);
  void print(std::ostream&) const;
  DECLARE_EXEC
  static ActionCommand* create(); //!< parse to create new object

  static ContextWarning<> filterafteracq_warning; //!< filter after final acq
  static ContextWarning<> notcoherencefilter_warning; //!< warning that matrix is not a coherence filter

private:
  std::string name; //!< name of filter mask
  const filter_def& filter; //!< reference to matrix
};

typedef ActionCommand* (*ActionCommand_function)();
typedef FASTMAPTYPE(ActionCommand_function) Action_Factory_t;
Action_Factory_t& get_Action_Factory(); //!< get registry for pulse sequence directives
extern LIST<actionstack_t> actionstacks; //!< pulse sequence(s)
//extern size_t acq_count; 
//void output_line(const char* a); //!< write \a a to current ``output''
//extern double dt; //!< current direct dimension dwell time (set by ActionDirectAcq)
double dwelltime();

extern ThreadWarning<> sync_warning; //!< couldn't find sync
extern ThreadWarning<> disablingcombine_warning; //!< disabling propagator combination due to time dependence
extern bool insidecontrol; //!< true if inside control structure

#endif
