#ifndef NMRsim_RF_h_
#define NMRsim_RF_h_

/*! \file
  \brief  Header file defining RF related components */

#include "NMRsim.h"
#include "NMRsim_MasterObj.h"
#include "Sequence.h"
#include "PhaseModulatedPropagator.h"
#include <set>
#include "NMRsim_matrixspec.h"

extern ThreadWarning<> CWoffres_warning; //!< warn about off-resonance irradiation in CW sequences
extern ThreadWarning<> phasemodulation_cache_warning; //!< failed to cache all propagators
extern ThreadWarning<> transients_problem_warning; //!< problem with trailing transients
extern ThreadWarning<> null_sequence_warning; //!< sequence doesn't appear to do anything
extern ThreadWarning<> emptypulseevent_warning; //!< setting parameter for empty event
extern ContextWarning<> short_pulse_warning; //!< extremely short pulse (<1 ns)
extern ContextWarning<> ignoring_offset_warning; //!< ignoring offset (for ideal pulses)
extern ThreadWarning<> chebyshev_warning; //!< chebyshev propagation and static Hamiltonian
extern ContextWarning<> empty_store_warning; //!< store applied to empty sequence
extern ThreadWarning<> phasemodulation_failed_warning; //!< phase modulation could not be used
extern ThreadWarning<> syncfailed_warning;
extern ThreadWarning<> phasemodulation_timedivision_warning;
extern ThreadWarning<> phasemodulation_notused_warning; //!< couldn't use phase modulation
extern ThreadWarning<> syncwithlist_warning; //!< synchronisation condition specified with event list
extern ThreadWarning<> ignoring_transients_warning; //!< ignoring -transients qualifier
extern ThreadWarning<> phasemodulation_syncconflict;
extern ThreadWarning<> phasemodulation_smalltimebase_warning; //!< timebase is small

template<class OpGen> matrixspec parse_coherencematrix(const OpGen&, const BaseList<size_t>&); //!< parse coherence (filter) matrix


struct usage_t;
class CycledSequence;
class MasterObj;
//! Propagator cache for sequential problems
struct Usequentialcache {
  void clear() { Ulist.clear(); } //!< clear cache

  bool empty() const { return Ulist.empty(); } //!< return true if cache not active
  usage_t usage() const; //!< return memory usage

  void initialise(const MasterObj&, double t, const BlockedMatrix<complex>&, double durv); //!< initialise with 1 cycle propagator referenced to time \a t
  const BlockedMatrix<complex>& operator()(const MasterObj&, double t, size_t n, double durv); //!< get propagator for \a n cycles, creating if necessary
  //  const BlockedMatrix<complex>* find(const MasterObj&, double t, size_t n, double durv) const; //!< return pointer to cached propagator for \a n cycles, NMRSIM_NULL if not cached
  const BlockedMatrix<complex>* rawfind(size_t n) const; //!< return pointer to cached propagator for \a n cycles, NMRSIM_NULL if not cached
  bool isvalid(const MasterObj&, double t, double dur) const; //!< \c true if cache is valid for this \a t and state
  LIST< BlockedMatrix<complex> > Ulist;
  size_t lastn,incr;
  hamiltonian_state_t state; //!< hamiltonian state for which cache is valid
  double phase; //!< rotor phase from which propagators calculated
  double duration_;  //!< duration of 1 cycle

  static ThreadWarning<> ordering_warning; //!< propagators not incrementally ordered
  static ThreadWarning<> irregular_warning; //!< propagators not regularly incremented
};

struct Udurcache {
  void clear() { U.clear(); }
  bool empty() const { return U.empty(); }
  usage_t usage() const;
  const BlockedMatrix<complex>* operator()(const MasterObj&, double t, double dur) const;  //!< \c true if cache is valid for this \a t and state
  void store(const BlockedMatrix<complex>&, const MasterObj&, double t, double dur);
  
  BlockedMatrix<complex> U;
  double phase;
  double duration;
  hamiltonian_state_t state;
};
  
struct cacheset_t {
  void flush() { Ucache_.clear(); Udurcache_.clear(); }
  Usequentialcache Ucache_;
  Udurcache Udurcache_;
  usage_t usage() const { return Ucache_.usage()+Udurcache_.usage(); }
};

struct sequencestate_t {
  sequencestate_t();
  void flush();
  size_t curptr;
  mutable double duration_;
  double lastoffset_; //!< offset into current step where last use finished
  void reset(bool resetduration =true);
  cacheset_t* getcache(double);
  cacheset_t& setcache(double);
  //  void lastoffset(double);
  cacheset_t cachezero_;
  cacheset_t cachenonzero_;
  double cacheoffset_;
  usage_t usage() const { return cachezero_.usage()+cachenonzero_.usage(); }
};

//! class for cycled list of sequence fragments
class CycledSequence : protected ListList<PhasedSequence*>
{
public:
  CycledSequence();
  CycledSequence(CompSequenceBase&);  //!< initialise from simple sequence fragment
  CycledSequence(PhasedSequence*);
      
  const char* name() const { return name_.c_str(); }
  void name(const char* namev) { name_=namev; }

  double duration() const; //!< returns current sequence duration

  //  using ListList<PhasedSequence*>::clear;  //!< show explicitly that clear refers to ListList

  void flush(); //!< flush any cached propagators / state

  size_t currentindex() const;
  BaseList<PhasedSequence*> currentlist() const { return (*this)(currentindex()); }  //!< return current sequence list
  PhasedSequence* current() const; //!< pointer to current sequence fragment

    //! return cycle length
  size_t cyclelength() const { return ListList<PhasedSequence*>::size(currentindex()); }

  //  size_t size() const { return isspecial_ ? ListList<PhasedSequence*>::size() : 1; }
    //! call after each use (increments if reqd)
  void next();

    const char* prettyname() const; //!< return sequence name or <anonymous> if none
  void reset(bool resetduration =false); //!< reset current fragment pointer
    
    void checkdirty(); //!< set dirty flag
    void checksync(double dt); //!< checkany synchronisation / PM for characteristic time dt
    bool wasdirty() const { return dirty_; } //!< return \c true if component fragments were rebuilt

    friend class CycledSeqBuilder;
    void print(std::ostream&) const; //!< print contents
    void printshort(std::ostream&) const; //!< print name or contents if none
    void printvariable(std::ostream&, subsid_t) const; //!< print variable name
    
  size_t arraytag() const { return arraytag_; }
  //    size_t sumtag() const { return sumtag_; }
  //    const BaseList<size_t>& arraytags() const { return arraytags_; }
  bool hascycle() const;  //!< return \c true if contains [] list
  bool hasCW() const; //!< return \c true if contains CW fragment
  const ListList<PhasedSequence*>& aslist() const { return *this; }
  sequencestate_t& currentstate() { return states_(currentindex()); }
  const sequencestate_t& currentstate() const { return states_(currentindex()); }
   
    static subsid_t index_to_subsid(size_t); //!< translate internal index into suitable parameter for set
    static size_t subsid_to_index(subsid_t);
    void set(double,size_t);

//   Usequentialcache& cache() { return states_(currentindex()).Ucache_; }
//   Udurcache& durationcache() { return states_(currentindex()).Udurcache_; }
//   double lastoffset() const { return states_(currentindex()).lastoffset_; }
//   void lastoffset(double);
    
  static ThreadWarning<> nonuniform_warning; //!< components have different durations - likely to be a mistake!

  usage_t usage() const;

 private:
    std::string name_; //! sequence name (empty if anonymous)
    size_t curptr; //!< current fragment pointer
  // bool isspecial_; //!< \c true for [] rather than {} list
    bool dirty_; //!< \c true if component fragments were rebuilt
  //  LIST<size_t> arraytags_; //!< {} virtual dimensions
    size_t arraytag_; //!< [] virtual dimension (0 if none)
  LIST<sequencestate_t> states_;
};

inline std::ostream& operator<< (std::ostream& ostr, const CycledSequence& a) {
  a.print(ostr);
  return ostr;
}

//! raw definition of pulse

//! \note ::VariableBase is used because it can store everything from a simple scalar
//! to a complex variable.  Although wasteful for simple variables, the only
//! significant efficiency losses occur once during parsing
struct RawPulseDef {
  RawPulseDef() : rf(NMRSIM_NULL), phase(NMRSIM_NULL), offset(NMRSIM_NULL), master_subsid(S_NONE) {}
  VariableBase* rf; //!< RF nutation rate
  VariableBase* phase; //!< RF phase
  VariableBase* offset; //!< RF offset
  int flags; //!< pulse flags
  subsid_t master_subsid; //!< last subsid flag (identifies when EventID completely updated)
  bool isconst() const; //!< \c true if definition is constant
};

class CompSequenceBase; // forward declaration

typedef LIST< smartptr<RFEvent,false> > eventlist_t;

//! Describes a complete \c pulse directive
struct EventID : public Setable {

  enum id_t { ideal, soft, automatic }; //!< pulse types

  EventID(CompSequenceBase* seqp, double nomdurv)
    : seqp_(seqp), type(soft), nomdur_(nomdurv), durs_(NMRSIM_NULL), chanoff(-1) { } //!< create delay period

  EventID(CompSequenceBase* seqp, size_t nchans, size_t chanoff_,id_t type_, double nomdurv, char sync_)
    : seqp_(seqp), type(type_), nomdur_(nomdurv), durs_(NMRSIM_NULL), sync(sync_) { init(nchans,chanoff_); } //!< create soft pulse

  EventID(CompSequenceBase* seqp, size_t nchans, size_t chanoff_) //!< create hard pulse
    : seqp_(seqp), type(ideal), nomdur_(0.0), durs_(NMRSIM_NULL), sync('|') { init(nchans,chanoff_); } 

  size_t channels() const { return pulses.size(); } //!< number of RF channels (1 if asychronous)
  void duration(double); //!< set duration
  //  double duration() const { return (type==ideal) ? 0.0 : nomdur_; } //!< "ideal" duration of base event (actual duration could vary due to timing errors)
  double totalduration() const; //!< duration of total event

  void set(double, subsid_t i); //!< modify parameter (index \a i encodes the channel number and parameter) 
  void set(const BaseList<double>&, subsid_t); //!< list variant
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;

  CompSequenceBase* seqp_; //!< pointer to ::CompSequenceBase using pulse definition
  id_t type; //!< timing type
  double nomdur_; //!< nominal duration used to calculate tip angle
  VariableBase* durs_; //!< durations (if list)

  //  int size(bool allownull =false) const; //!< return number of pulse elements. If incompatible, return 0 if \a allowmismatch or fail with error. Return <0 if undefined
  int size() const; //!< return number of pulse elements. If incompatible, return 0 if \a allowmismatch or fail with error. Return <0 if undefined

  bool isconst() const; //!< \c true if sequence fragment is constant

  //! RFEvents created from pulse directive
  /** \note The memory inefficient List< List<T> > is used because 
      a list of raw events may be associated with each channel
      and the size of this list could change during execution */
  LIST<eventlist_t> events;

  char sync; //!< timing synchronisation flag (+, - etc.)
  size_t chanoff; //!< channel number for an asynchronous \c EventID (0 for synchronous \c EventIDs)

  bool isdelay() const { return pulses.empty(); } //!< \c true if event is a delay period (no RF)
  LIST<RawPulseDef> pulses; //!< parsed \c pulse directives for each channel

  static EventID* create(CompSequenceBase*, id_t, char sync_); //!< create \c EventID from text input

  double inset() const; //!< duration of pulse \em before timing point (e.g. half duration for | synchronised pulse)
  void init(size_t n,size_t o); //!< initialise EventID for \a n channels starting at \a o
  bool buildevents(); //!< rebuild \c RFEvents
  void ensureevents(); //!< ensure \c RFEvents are built

  bool checknull(size_t chan); //!< warn if channel sequence fragment is empty (and return \c true)
  void set(double, size_t chan, subsid_t); //!< internal set
  void set(const BaseList<double>&, size_t, subsid_t); //!< internal list set
  template<typename T> void set_(const T&, subsid_t);

  static ThreadWarning<> resizelist_warning; //!< warning that list has been resized
};

//! Propagator cache object
/** \note Cache checks that rotor phases match but nothing else i.e.
    cache must be cleared "manually" if sequence, Hamiltonian etc. have changed */
struct Ucache_obj {
  void clear() { Ulist.clear(); } //!< clear cache
  //  void update_time(const MasterObj&, double t, double dur); //!< clear cache if rotor phase is different (and update)
  
  //! return reference to cached propagator or empty slot if new sequence fragment 
  BlockedMatrix<complex>& operator()(const CompSequenceBase&, const MasterObj&, double, double);

  usage_t usage() const; //!< return memory usage

  LIST< BlockedMatrix<complex> > Ulist;
  hamiltonian_state_t state; //!< hamiltonian state for which cache is valid
  double phase; //!< rotor phase for which cache (first point) is valid
  double duration; //!< step duration for which cache is valid
  int nearest(double) const; //!< return closest cache slot (<0 if none)
  static ThreadWarning<> reset_warning; //!< forced reset
};
   
//! (abstract) base class for storing sequence fragments
/** \todo Too much of this class is public; needs a better defined API */
class CompSequenceBase : public Setable {
public:
  CompSequenceBase();
  virtual ~CompSequenceBase() {}

  const char* name() const { return name_.c_str(); }
  void name(const char* namev) { name_=namev; }

  LIST<Sequence> seqlist; //!< corresponding Sequence's for use the SequencePropagator
  //! channel durations
  /** \todo Does this really belong in the class?  Channel durations must be equal for valid fragment */
  LIST<double> times;
  virtual CompSequenceBase* clone() const =0;
  bool isCW; //!< \c true if sequence corresponds to CW irradiation
  double synctime; //!< sychronisation time (between RF and MAS); 0 if unset

  //! return fragment store to empty state
  /** \note <B>This must be called by derived functions</B> */
  virtual void clear();
  void set(double,subsid_t);
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  void reset(dirty_t);  //!< clear propagator caches

  usage_t usage() const { return Ucache.usage(); } //!< memory usage of propagator cache
 
  //LIST<transient_state> transient_states;
  //void store_transients(double t) { ::store_transients(transient_states,t); }
  //void restore_transients() const { ::restore_transients(transient_states); }
    
  bool empty() const
  { return seqlist.empty(); }  //!< \c true if store is empty

  bool isdirty() const { return (status_!=DIRTY_CLEAN); } //!< \c true unless in a clean state (built Sequences match specification)
  bool wasdirty() const; //!< \c true if sequence has been rebuilt

  //! set dirty status
  /** if \a recurse flag is set, propagate recursively to fragments that refer to this one */
  void flagdirty(dirty_t);

  double duration() const;  //!< return sequence duration
  double safe_duration(); //!< return duration, rebuilding if required

  void refresh();   //!< ensure object is clean, rebuilding as required

  void parse_pulse(EventID::id_t); //!< parse input and add to sequence
  size_t verbose_check_sync(double x,const char* nx, double y, const char* ny) const;

  virtual void push_back(EventID*) =0; //!< add individual pulse directive

  //! evaluate propagator starting at time \a t for duration \a dur with RF phase shift \a p
  BlockedMatrix<complex>& evaluate(MasterObj*,double t, double dur, bool disablecache =false, double p =0.0, double offset =0.0) const;
  
  double synchronisation_time() const { return synctime; }
  const smartptr<PhaseModulation,false>& modulation() const { return phasemodp; }
  smartptr<PhaseModulation,false>& modulation() { return phasemodp; }

  smartptr<PhaseModulatedPropagator,false>& generator() { return propgenp; }
  const smartptr<PhaseModulatedPropagator,false>& generator() const { return propgenp; }

  PhaseModulatedPropagator& updatepm(const MasterObj&);

  bool isphasemodulated() const { return (phasemodp.get()!=NMRSIM_NULL); } //!< return \c true if sequence is phase modulated
  //  void buildpm(double dt); //!< rebuild phase modulation (e.g. following change in synchronisation)
  void checksync(double dt); //!< check synchronisation

  int prop_flags() const { return ::prop_flags(partition_index); } //!< return appropriate propagation flags
  const matrix_partition_set& partitioning() const { return partitions(partition_index); } //!< return appropriate ``partitioning'' (Chebyshev propagation)
  size_t partition_index; //!< index to appropriate partitioning
  bool needsresync() const; //!< \c true if PhaseModulation object is incomplete

private:  
  std::string name_; //!< name of fragment

  void rebuild(); //!< rebuild \c Sequences and tidy object
  void checkpm(dirty_t); //!< check whether phase modulation can be used
  virtual double raw_rebuild() =0; //!< \internal raw rebuild

  smartptr<PhaseModulation,false> phasemodp; //!< pointer to ::PhaseModulation object (NMRSIM_NULL if phase modulation not active)
  smartptr<PhaseModulatedPropagator,false> propgenp; //!< pointer to ::PhaseModulatedPropagator (propagator generator needs to be cached to be worthwhile)

  template<class T> void updatepm_(const T&); //!< \internal

  size_t synchint; //!< hint for number of steps per rotor cycle
  //  mutable BlockedMatrix<complex> Ucache; //!< propagator cache for static problems
  mutable bool donesyncwarn; //!< \true if warning about bad synchronisation already given
  mutable Ucache_obj Ucache;
  friend struct Ucache_obj;

 //! "fragment doesn't do anything flag"
  /** \todo Check whether this can't be factored out of the class def */
  mutable bool doneerror;

  //! duration of fragment
  /** Value is only valid if object is clean */
  mutable double duration_;
  mutable dirty_t status_; //!< "dirty" status
  mutable size_t lastusedmult; //!< check whether PhaseModulation being used with conflicting synchronisation
};

std::ostream& operator<< (std::ostream&, const CompSequenceBase&);

//! Concrete CompSequence class for synchronous fragments

//! See ::CompSequenceBase for usage of member of functions
class SyncSequence : public CompSequenceBase {
public:
  CompSequenceBase* clone() const { return new SyncSequence(*this); }
  void clear(); //!< empty object
  void push_back(EventID* evlist) { rawrep.push_back(evlist); } //!< add ::EventID to sequence

private:
  LIST<EventID*> rawrep; //!< ordered list of RF events
  double raw_rebuild(); //!< \internal
  void flush_delay(double); //!< output accumulated delay period
};

//! Concrete CompSequence class for asynchronous fragments
/*!
 * See ::CompSequenceBase for usage of member of functions
 */
class AsyncSequence : public CompSequenceBase {
public:
  AsyncSequence();

  CompSequenceBase* clone() const { return new AsyncSequence(*this); }

  void clear();
  void push_back(EventID* evlist); //!< add event list to sequence fragment

private:
  LIST< LIST<EventID*> > rawrep; //!< separate ordered list of RF events per channel
  double raw_rebuild();
  void push_back(size_t chan, const BaseList<EventID*>& evl, size_t);
  void flush_delay(size_t n, double); //!< output accumulated delay on channel \a n
};

//! Sequence fragment + setable phase shift
class PhasedSequence : public Setable {
public:
  PhasedSequence(CompSequenceBase& seq_, double phase_ =0.0) 
    : seq(seq_), phase(phase_) {}

  void set(double fval,subsid_t =S_PHASE); //!< set phase shift
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;

  bool empty() const { return seq.empty(); }  //!< \c true if sequence fragment is empty
  bool isdirty() const; //!< \c true if sequence fragment needs rebuilding
  double duration() const { return seq.duration(); } //!< sequence duration
  const CompSequenceBase& get() const { return seq; }
  CompSequenceBase& get() { return seq; }

  const BlockedMatrix<complex>& evaluate(MasterObj*,double t, double dur, bool disablecache =false, double offset =0.0) const; //!< evaluate propagator starting at time \a t
private:
  CompSequenceBase& seq; //!< reference to sequence fragment
  double phase; //!< phase shift
};

typedef MAPTYPE(CycledSequence*) cycledseqmap_type;
extern cycledseqmap_type cycledseqmap; //!< map of phase cycled sequences
extern LIST<CycledSequence*> cycledseqlist; //!< list of all cycled sequences

typedef std::set<CompSequenceBase*> dirtylist_t;
extern dirtylist_t dirty_list; //!< set of "dirty" sequences

extern LIST<smartptr<PulseGeneratorBase> > pgenstack; //!< pulse generators for each RF channel
PulseGeneratorBase& pulse_generator(size_t chan =0U); //!< return pulse generator for channel \a chan

void subsid_to_rf(size_t& chan, subsid_t& actsubsid, int subsid); //!< convert hashed \a subsid into channel + true subsid 
void reset_pgens(); //reset pulse generators
void reset_seqlists(); //!< reset cycled sequences

void rebuild_sequences(dirty_t =DIRTY_ALL); //!< rebuild any dirty sequences
void store_transients(LIST<transient_state>&, double t); //!< store current state of RF transients
void restore_transients(const LIST<transient_state>&); //!< restore previously stored transient state
void flush_transients(BlockedMatrix<complex>& U); //!< add any remaining RF transients to propagator \a U
void flush_transients(); //!< discard any remaining RF transients
bool transients_active(); //!< returns \c true if unflushed transients present

#endif
