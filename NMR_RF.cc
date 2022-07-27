// arbitrary separation of "sequence" related stuff to lighten NMR.cc

#include "NMRsim_RF.h"
#include "Parser.h"
#include <sstream>

//! parse pulse directive
/** qualifier distinguishes hard, soft, auto */
void parse_pulse(int);
void parse_reset();
void parse_store();
void parse_channel();
void parse_delay();

#ifdef NMRSIM_INTUITIVE_OFFSET
inline double convertoffset(double x) { return -x; }
#else
inline double convertoffset(double x) { return x; }
#endif

inline bool quickphasemod(const CompSequenceBase& seq) { return seq.modulation().get() && optcache(); }
//inline bool quickphasemod(const Sequence& seq) { return false; }

//! don't use (phase modulation) steps below this fraction of maxdt
#define NMRSIM_MINSTEP 0.6

dirtylist_t dirty_list;
bool havestore=false;

double phasetolerance=5e-4; // match that must be achieved on rotor phase (degrees) NB. doesn't account for "jitter". Can't be too tight as phase at long times will not be calculated accurately
double tolerance=5e-9; //5 ns
option optcombinepropagators("combinepropagators","combine propagators",option::AUTO,option::NOTUSED);
option optsmartprop("smartprop","",option::AUTO,option::NOTUSED);
LIST<smartptr<PulseGeneratorBase> > pgenstack;
LIST<matrix_partition_set> partitions;
LIST<int> prop_flags;
//LIST<double> transient_amps;
//! RF transient model
enum trans_t { TRANS_NONE, //!< no transients
	       TRANS_AUTO, //!< simple "Vega" model global per channel
	       TRANS_MANUAL //!< specified per pulse
 };
trans_t transient_model=TRANS_NONE; //!< RF transient model (default \c TRANS_NONE)
//! "manual" mode does not involve trailing transients
//bool inline hastransients() { return (transient_model==TRANS_SIMPLE); }
size_t usedchannels=0;
size_t curchannel=0; //!< current channel (::AsyncSequence build)

struct SysVar_trans : public SystemVariable<double*> {
  SysVar_trans(const char* name, size_t which)
    : SystemVariable<double*>(name,&value_), which_(which), value_(0.0) {}
  void update();
  size_t which_;
  double value_;
};

LIST<SysVar_trans*> transient_amps;
option optphasemod("phasemodulation","phase modulated RF");

namespace {
  const double rad_to_deg=180.0/M_PI;
  const double deg_to_rad=M_PI/180.0;
  
  struct Proxy_ {
    Proxy_() {
      //      optional_map_t& optional_map(get_optional_map());
      //optional_map["phasemodulation"]=&optphasemod;
      //optional_map["combinepropagators"]=&optcombinepropagators;
      //optional_map["smartprop"]=&optsmartprop;
      add_option(optphasemod);
      add_option(optcombinepropagators);
      add_option(optsmartprop);

      command_Factory_t& par_Factory(get_par_Factory());
      par_Factory["channel"]=par_t(&parse_channel,true);
      par_Factory["delay"]=par_t(&parse_delay,true);
      par_Factory["pulse"]=par_t(&parse_pulse,EventID::soft,true);
      par_Factory["pulseid"]=par_t(&parse_pulse,EventID::ideal,true);
      par_Factory["pulseauto"]=par_t(&parse_pulse,EventID::automatic,true);
      par_Factory["store"]=par_t(&parse_store,true);
    }
  };  

  Proxy_ proxy;
}

void sequencestate_t::reset(bool resetdur)
{
  curptr=0; //!< reset current fragment pointer
  lastoffset_=0.0; //!< no trailing fragment
  if (resetdur)
    duration_=-1.0;
}
  
void CycledSequence::reset(bool resetduration)
{
  if (states_.empty()) {
    states_.create(size());
    return;
  }
  for (size_t i=states_.size();i--;)
    states_(i).reset(resetduration);
}

CycledSequence::CycledSequence()
  : dirty_(false)
{ reset(true); }

CycledSequence::CycledSequence(CompSequenceBase& seq)
  : ListList<PhasedSequence*>(ExplicitList<1,size_t>(1U),new PhasedSequence(seq,0.0)), 
    name_(seq.name())
    { reset(true); }

CycledSequence::CycledSequence(PhasedSequence* pseq)
  : ListList<PhasedSequence*>(ExplicitList<1,size_t>(1U),pseq)
{ reset(true); }

subsid_t CycledSequence::index_to_subsid(size_t ind)
{
  return subsid_t(ind+100);
}

size_t CycledSequence::subsid_to_index(subsid_t subsid)
{
  int ind=int(subsid)-100;
  if (ind<0) {
    std::cerr << "Bad subsid_to_index parameter: " << subsid << '\n';
    error_abort();
  }
  return ind;
}
    
void SysVar_trans::update()
{
  //  if (!hastransients()) {
  //    PulseGenerator_SimpleTransients* pgenp(dynamic_cast<PulseGenerator_SimpleTransients*>(pgenstack(which_).get()));
  //    if (pgenp)
  pgenstack(which_)->transient_amplitude(transient_amps(which_)->get());
      //    else
      //      throw InternalError("transient_amplitude");
    //flagdirty(DIRTY_HAMILTONIAN);
  //  }
  //! invalidate propagators 
  /** A touch inefficient since only sequences involving the relevant channels are affected */
  update_propagators=DIRTY_ALL; 
}

template<class T> void set_partitioning(T& propgen, const CompSequenceBase& seq)
{
  const matrix_partition_set& part(seq.partitioning());
  if (!part.empty())
    propgen.partitioning(part);
}

void parse_transients()
{
  const char* com=parse_string(F_REPLACEDOLLAR);
//   if (strcmp(com,"none")==0) {
//     transient_model=TRANS_NONE;
//     return true;
//   }
  if (strcmp(com,"manual")==0)
    transient_model=TRANS_MANUAL;
  else {
    if (strcmp(com,"automatic")==0)
      transient_model=TRANS_AUTO;
    else
      error_abort("Unknown transient model");
  }

  char name[32];
  bool nonzero=false;
  transient_amps.create(nchannels);
  for (size_t j=0;j<nchannels;j++) {    
    sprintf(name,"transient_amplitude%" LCM_PRI_SIZE_T_MODIFIER "u",j+1);
    transient_amps(j)=new SysVar_trans(name,j);
    parse_system_variable(*(transient_amps(j)));
    if (transient_amps(j)->get())
      nonzero=true;
    //    Variable cvar(S_ARG1,varp);
    //transient_amps(j)=parse_double(&cvar);
  }
  if (nonzero && optphasemod.isauto())
    optphasemod.set(option::OFF);
  //  parse_system_variable(v_trans);
}

size_t CompSequenceBase::verbose_check_sync(double x,const char* nx, double y, const char* ny) const
{
  const size_t sync=check_sync(x/y);
  if (!donesyncwarn && syncfailed_warning.enabled() && (sync==0)) {//!< check enabled as build is expensive
    std::ostringstream str(std::ostringstream::out);
    str << ny << " (";
    prettyprint_time(y,str) << ") does not divide into " << nx << " (";
    prettyprint_time(x,str) << ") for sequence " << name_ << " (future warnings for this sequence suppressed)";
    syncfailed_warning.raise(str.str().c_str());
    donesyncwarn=true;
  }
  return sync;
}

smartptr<ZshiftCache,false> zshift_cachep;

simple_counter global_cache(NMRSIM_CACHE_LIMIT*1024*1024);

Spectrometer spectrometer;
double tres=0.0;

struct SysVar_tres : public SystemVariable<double*> {
  SysVar_tres(const std::string& name_, double* value_)
    : SystemVariable<double*>(name_,value_,1e6) {}
  void update() {
    spectrometer.time_resolution(tres);
    update_propagators=DIRTY_ALL;
    //    flagdirty(DIRTY_ALL); //all sequences need rebuilding
  }
};


SysVar_tres v_tres("time_resolution",&tres);

void parse_time_resolution()
{
  if (are_left())
    parse_system_variable(v_tres);
  else
    std::cout << "time_resolution=" << (tres*1e6) << " us\n";
}

template<class HType> void MasterObj::propagator_(BlockedMatrix<complex>& Udest, const HType& H, const CompSequenceBase& seq, double t1,double t2,double origin)
{
  SequencePropagator propgen(H,get_maxdt(),seq.seqlist,origin,seq.duration(),tolerance,(verbose & VER_GEN) ? verbose_level : 0,seq.prop_flags());
  set_partitioning(propgen,seq);
  propgen.synchronisation_hint(seq.synchronisation_time());
  propgen(Udest,t1,t2);
}

template<class HType> void MasterObj::propagator_(BlockedMatrix<complex>& Udest, const HType& H, const Sequence& seq, double t1,double t2,double origin)
{
  SequencePropagator propgen(H,get_maxdt(),seq,origin,seq.duration(),tolerance,(verbose & VER_GEN) ? verbose_level : 0);
  //  propgen.synchronisation_hint(seq.synchronisation_time());
  propgen(Udest,t1,t2);
}

ThreadWarning<> phasemodulation_cache_warning("phase modulation: failed to cache all propagators.  Increase cache_limit or use verbose -general for more information",&NMRsim_once_warning);

template<class HType> void CompSequenceBase::updatepm_(const HType& H)
{
  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  const double maxdt=get_maxdt();
  if (propgenp.get()) {
    const double timeoff=gammatimeoffset();
    if (lverbose>1)
      std::cout << "PhaseModulation: updating PhaseModulatedPropagator gamma / time offset to " << (timeoff*1e6) << " us\n";
    propgenp->Hsys_offset(timeoff);
  }
  else {
    const matrix_partition_set* partp=partitioning().empty() ? NMRSIM_NULL : &(partitioning());
    propgenp.clear(); // clear old to clear cache
    size_t jumpsteps=0;
    double tcommon=0.0;
    static double lastjumpprint=-1.0;
    if (maxjumpdt) {
      jumpsteps=int(0.98+H.period()/(2.0*maxjumpdt));
      tcommon=(phasemodp->tick()<maxdt) ? phasemodp->tick() : maxdt;
      if (lverbose) {
	const double jumpprint=0.5e6*H.period();
	if (jumpprint!=lastjumpprint) {
	  std::cout << "PhaseModulation maximum jump time: " << (0.5e6*H.period()/jumpsteps) << " us (" << jumpsteps << " step(s) over rotor period)  Cached time: ";
	  prettyprint_time(tcommon) << '\n';
	  lastjumpprint=jumpprint;
	}
      }
    }
    else {
      if (lverbose && (lastjumpprint!=0.0)) {
	std::cout << "PhaseModulation maximum jump time: disabled\n";
	lastjumpprint=0.0;
      }
    }
    if (lverbose>1)
      std::cout << "PhaseModulation: creating new PhaseModulatedPropagator\n";
    PhaseModulatedPropagator* newpropgenp = new PhaseModulatedPropagator(H,partp,maxdt,master_time,*phasemodp,gammatimeoffset(),&global_cache,tolerance,lverbose,prop_flags(),jumpsteps,tcommon);
    if (!(newpropgenp->cacheok()))
      phasemodulation_cache_warning.raise();
    propgenp.reset(newpropgenp);
  }
}

template<class HType> void MasterObj::spin_propagator_(BlockedMatrix<complex>& Udest, const HType& H, const CompSequenceBase& seq, double t1,double t2,double origin)
{
  if (quickphasemod(seq)) {
    PhaseModulatedPropagator& propgen(const_cast<CompSequenceBase&>(seq).updatepm(*this));
    //    set_partitioning(propgen,seq);
    propgen(Udest,t1,t2);
  }
  else
    propagator_(Udest,H,seq,t1,t2,origin);
}

template<class HType> void MasterObj::spin_propagator_(BlockedMatrix<complex>& Udest, const HType& H, const Sequence& seq, double t1,double t2,double origin)
{
  propagator_(Udest,H,seq,t1,t2,origin);
}

PhaseModulatedPropagator& CompSequenceBase::updatepm(const MasterObj& obj)
{
  if (!(obj.isspinningH()))
    throw InternalError("updatepm: not spinning Hamiltonian!");
  if (obj.Hspinp.iscomplex())
    updatepm_(*(obj.Hspinp.get_complex()));
  else
    updatepm_(*(obj.Hspinp.get_real()));
  return *propgenp;
}

template<class SequenceHolder> void MasterObj::propagator(BlockedMatrix<complex>& U,const SequenceHolder& seq, double t1,double t2,double origin)
{
  //  if ((nchannels==0) && ((hamobjtype!=M_SPIN) || !quickphasemod(seq))) {
  if (nchannels==0) {
    propagator(U,t1,t2);
    return;
  }
  switch (hamobjtype) {
  case M_SPIN:
    if (Hspinp.iscomplex())
      spin_propagator_(U,*(Hspinp.get_complex()),seq,t1,t2,origin);
    else
      spin_propagator_(U,*(Hspinp.get_real()),seq,t1,t2,origin);
    break;
  case M_STATIC:
    if (Hstaticp.iscomplex())
      propagator_(U,*(Hstaticp.get_complex()),seq,t1,t2,origin);
    else
      propagator_(U,*(Hstaticp.get_real()),seq,t1,t2,origin);
    break;
  case M_STATICD:
    propagator_(U,*Hstaticdp,seq,t1,t2,origin);
    break;    
  default:
    throw Failed("Unsupported combination of Hamiltonian and RF");
  }
}

template void MasterObj::propagator(BlockedMatrix<complex>&, const CompSequenceBase&, double t1,double t2,double origin);
template void MasterObj::propagator(BlockedMatrix<complex>&, const Sequence&, double t1,double t2,double origin);

namespace {
  std::ostream& printhstate(const hamiltonian_state_t& state, std::ostream& ostr =std::cout)
  {
    return ostr << "Rotor angle=" << (state*rad_to_deg) << " degrees";
  }
}

ThreadWarning<> Ucache_obj::reset_warning("Unable to use previously calculated propagators because requests don't coincide with previous rotor phases / durations. If these are 'near misses' due to rounding errors, try relaxing the rotor phase tolerance (tolerance <value> -phase). Otherwise you may be able to improve efficiency by not re-using the same sequence at different points in the sequence or by changing the gamma angle sampling. Sequence: ",&NMRsim_once_warning);

// NB. don't bother checking whether matrices are used or not
usage_t Ucache_obj::usage() const {
  return Ulist.empty() ? usage_t() : usage_t(Ulist.size()*Ulist.front().size(),Type2Type<complex>());
}

// NB. don't bother checking whether matrices are used or not
usage_t Usequentialcache::usage() const {
  return Ulist.empty() ? usage_t() : usage_t(Ulist.size()*Ulist.front().size(),Type2Type<complex>());
}

usage_t Udurcache::usage() const {
  return U.empty() ? usage_t() : usage_t(U.size(),Type2Type<complex>());
}

// void Ucache_obj::update_time(const MasterObj& obj, double t, double dur)
// {
//   const bool isverb=(verbose & VER_GEN) && (verbose_level>1);
//   const double rphase(get_rotor_phase(t));
//   const hamiltonian_state_t newstate(obj.hamiltonian_state());
//   if (!(Umap.empty())) {
//     if ((fabs(rphase-phase)<1e-8) && (newstate==state) && (fabs(duration-dur)<1e-9)) {
//       if (isverb)
// 	std::cout << "Re-using cached propagator\n";
//       return;
//     }
//     if (isverb) {
//       std::cout << "Clearing MAS propagator cache for rotor phase=" << rphase << " degrees, ";
//       printhstate(obj.hamiltonian_state()) << "  Duration: ";
//       prettyprint_time(dur) << '\n';
//     }
//     clear();
//   }
//   duration=dur;
//   state=newstate;
//   phase=rphase;
// }

int Ucache_obj::nearest(double rphase) const
{
  const size_t n=Ulist.size();
  if (n==0)
    throw InternalError("Ucache_obj::nearest");
  double delta=rphase-phase;
  if (delta<0)
    delta+=360.0;
  size_t ind=0;
  if (fabs(delta)>phasetolerance) {
    if (n>1) {
      const double step=360.0/n; 
      ind=size_t(delta/step+0.5);
      delta-=ind*step;
    }
    else {
      if (delta>180.0)
	delta-=360.0;
    }
  }
  return (fabs(delta)>phasetolerance) ? -1 : ind;
}

void Usequentialcache::initialise(const MasterObj& obj, double t, const BlockedMatrix<complex>& U, double durv)
{
  phase=get_rotor_phase(t);
  state=obj.hamiltonian_state();
  Ulist.create(3U);
  Ulist.front()=U;
  lastn=incr=0;
  duration_=durv;
}

bool Usequentialcache::isvalid(const MasterObj& obj, double t, double durv) const
{
  if (empty())
    return false;
  const bool beverb=(verbose & VER_GEN) && (verbose_level>1);
  const double rphase(get_rotor_phase(t));
  if (fabs(rphase-phase)>phasetolerance) {
    if (beverb)
      std::cout << "Usequential cache: rejecting cache match because current rotor phase (" << rphase << " doesn't match stored (" << phase << ") within tolerance of " << phasetolerance << "\n";
    return false;
  }
  if (obj.hamiltonian_state()!=state) {
    if (beverb)
      std::cout << "Usequential cache: rejecting cache match because Hamiltonian states differ\n";
    return false;
  }
  if (fabs(durv-duration_)>tolerance) {
    if (beverb)
      std::cout << "Usequential cache: rejecting cache match because periods don't match\n";
    return false;
  }
  return true;
}

const BlockedMatrix<complex>* Usequentialcache::rawfind(size_t n) const
{
  if (empty())
    throw Failed("Usequentialcache used before initialisation");
  if (n==0)
    throw InvalidParameter("Usequentialcache::find");
  if (n==1)
    return &(Ulist.front());
  if (n==lastn)
    return &(Ulist(1U));
  if (n==incr)
    return &(Ulist(2U));
  return NMRSIM_NULL;
}

const BlockedMatrix<complex>* Udurcache::operator()(const MasterObj& obj, double t, double dur) const
{
  if (empty())
    return NMRSIM_NULL;
  const bool beverb=(verbose & VER_GEN) && (verbose_level>1);
  if (fabs(dur-duration)>tolerance) {
    if (beverb)
      std::cout << "Uduration cache: rejecting cache match because durations differ\n";
    return NMRSIM_NULL;
  }
  const double rphase(get_rotor_phase(t));
  if (fabs(rphase-phase)>phasetolerance) {
    if (beverb)
      std::cout << "Uduration cache: rejecting cache match because rotor phase doesn't match\n";
    return NMRSIM_NULL;
  }
  if (obj.hamiltonian_state()!=state) {
    if (beverb)
      std::cout << "Uduration cache: rejecting cache match because Hamiltonian states differ\n";
    return NMRSIM_NULL;
  }
  return &U;
}
  
void Udurcache::store(const BlockedMatrix<complex>& Ulast, const MasterObj& obj, double t, double dur)
{
  state=obj.hamiltonian_state();
  phase=get_rotor_phase(t);
  duration=dur;
  U=Ulast;
}

// const BlockedMatrix<complex>* Usequentialcache::find(const MasterObj& obj, double t, size_t n, double dur) const
// {
//   if (!isvalid(obj,t,dur))
//     return NMRSIM_NULL;
//   return rawfind(n);
// }

ThreadWarning<> Usequentialcache::ordering_warning("could not optimise re-use of propagators (smartprop) as not used in increasing order",&NMRsim_once_warning);
ThreadWarning<> Usequentialcache::irregular_warning("could not optimise re-use of propagators (smartprop) as not regularly incremented",&NMRsim_once_warning);

const BlockedMatrix<complex>& Usequentialcache::operator()(const MasterObj& obj, double t, size_t n, double durv)
{
  if (!nochecks && !isvalid(obj,t,durv))
    throw Failed("Usequentialcache used improperly");

  const BlockedMatrix<complex>* Up(rawfind(n));
  if (Up)
    return *Up;
  if ((lastn==0) || (n<lastn)) {
    if (n<lastn)
      ordering_warning.raise();
    pow(Ulist(1U),Ulist.front(),n);
  }
  else {
    BlockedMatrix<complex> U;

    const size_t delt=n-lastn;
    if (delt==1)
      multiply(U,Ulist.front(),Ulist(1U));
    else {
      if (delt==lastn) {
	Ulist(2U)=Ulist(1U);
	incr=delt;
      }
      else {
	if (delt!=incr) {
	  if (incr)
	    irregular_warning.raise();      
	  pow(Ulist(2U),Ulist.front(),delt);
	  incr=delt;
	}
      }
      multiply(U,Ulist(2U),Ulist(1U));
    }
    Ulist(1U).swap(U);
  }
  lastn=n;
  return Ulist(1U);
}
    
BlockedMatrix<complex>& Ucache_obj::operator()(const CompSequenceBase& seq, const MasterObj& obj, double t, double dur)
{ 
  const int lverb=(verbose & VER_GEN) ? verbose_level : 0;
  if ((spin_rate==0) || (dur==0.0)) {
    Ulist.create(1U);
    if ((phase>=0.0) || duration) 
      Ulist.front().clear(); //!< empty any object
    phase=-1.0; //!< flag not spinning
    duration=0.0;
    return Ulist.front();
  }
  const double rphase(get_rotor_phase(t));
  const hamiltonian_state_t newstate(obj.hamiltonian_state());
  if (!(Ulist.empty())) {
    if ((newstate==state) && (fabs(duration-dur)<1e-9)) {
      const int ind=nearest(rphase);
      if (ind>=0) {
	if (lverb>1)
	  std::cout << "Found cache match at index " << (ind+1) << '\n';
	return Ulist(size_t(ind));
      }
    }
    if (lverb) {
      std::cout << "Clearing MAS propagator cache for rotor phase=" << rphase << " degrees, ";
      printhstate(obj.hamiltonian_state()) << "  Duration: ";
      prettyprint_time(dur) << '\n';
    }
    clear();
    reset_warning.raise(seq.name());
  }
  size_t nsteps=1;
  if (seq.synchint)
    nsteps=seq.synchint;
  else {
    const double rotor_period=1.0/spin_rate;
    if (dur>rotor_period)
      seq.verbose_check_sync(dur,"propagator duration",rotor_period,"rotor period");
    else {
      const size_t nsync=seq.verbose_check_sync(1.0/spin_rate,"rotor period",dur,"propagator duration");
      if (nsync)
	nsteps=nsync;
      else {
	if (!(seq.donesyncwarn) && syncfailed_warning.enabled()) {
	  seq.donesyncwarn=true;
	  std::cerr << "If appropriate, supply synchronisation hint in store for sequence: " << seq.name() << '\n';
	}
      }
    }
  }
  Ulist.create(nsteps);

  size_t offset=0U;
  phase=rphase;
  if (nsteps>1) {
    const double step=360.0/nsteps;
    offset=size_t(rphase/step);
    if (offset)
      phase-=offset*step;
  }
  if (lverb)
    std::cout << "Built MAS propagator cache with " << nsteps << " step(s) and initial rotor phase " << phase << " degrees\n";
  duration=dur;
  state=newstate;
  return Ulist(offset);
}

BlockedMatrix<complex> Utmp;

ThreadWarning<> transients_problem_warning("trailing transients and non-zero sequence phase shift detected!",&NMRsim_repeat_warning);      
ThreadWarning<> null_sequence_warning("sequence doesn't seem to do anything: ",&NMRsim_once_warning);

CompSequenceBase::CompSequenceBase()
    : synctime(0.0),
      partition_index(-1),
      synchint(0),
      donesyncwarn(false),
      status_(DIRTY_ALL),
      lastusedmult(0)
{}

BlockedMatrix<complex>& CompSequenceBase::evaluate(MasterObj* objp, double t, double dur, bool disablecache, double phase, double offset) const
{
  if (status_!=DIRTY_CLEAN)
    throw InternalError("CompSequenceBase::evaluate called before sequence rebuilt");

  //if static, use global store, otherwise use object store
  BlockedMatrix<complex>* Up=NMRSIM_NULL;
  bool transproblem=transients_active();
  if (!optcache || transproblem || disablecache || offset)
    Up=&Utmp;
  else {
    Up=&(Ucache(*this,*objp,t,dur));
    if (!!(*Up))
      optcache.setusage(true);
  }

  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  if ((Up==&Utmp) || !(*Up)) {

    if (lverbose) {
      std::cout << "Evaluating propagator for " << name_;
      if (dur) {
	std::cout << " from t=";
	prettyprint_interval(t,t+dur);
	if (offset)
	  std::cout << " (offset into sequence " << (offset*1e6) << " us)";
      }
      else {
	std::cout << " at t=";
	prettyprint_time(t);
      }
      std::cout << '\n';
    }
        
    objp->propagator(*Up,*this,t,t+dur,t-offset);
    if (lverbose>1)
      std::cout << (*Up);

    if (!transproblem && transients_active())
      transproblem=true;
    
    if (transproblem && phase)
      transients_problem_warning.raise();
    
    if (!(*Up) && !doneerror) {
      null_sequence_warning.raise(name(),true);
      doneerror=true; //!< implemented per sequence 
    }
    if (transproblem && (Up!=&Utmp)) { //prevent cacheing
      Utmp=*Up;
      Up->clear();
      Up=&Utmp;
    }
    if (lverbose)
      std::cout << "Caching resulting propagator: " << ((Up==&Utmp) ? "No\n" : "Yes\n");
  }
  else {
    if (lverbose) {
      std::cout << "Using cached propagator for " << name_ << " at t=";
      prettyprint_time(t) << '\n';
    }
  }
  if (phase) {
    Utmp=*Up;
    Up=&Utmp;
    (*zshift_cachep)(Utmp,deg_to_rad*phase); //don't modify cached value
  }
  return *Up;
}

const BlockedMatrix<complex>& PhasedSequence::evaluate(MasterObj* objp, double t, double dur, bool disablecache, double offset) const
{
  const BlockedMatrix<complex>& U(seq.evaluate(objp,t,dur,disablecache,phase,offset));
  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  if (lverbose) {
    std::cout << "Phase shift: " << phase << "\n";
    if (phase && (lverbose>1)) 
      std::cout << U;
  }
  return U;
}

LIST<CycledSequence*> cycledseqlist; //!< list of all cycled sequences
seqmap_type seqmap;

CompSequenceBase* acqseqp=NMRSIM_NULL;

cycledseqmap_type cycledseqmap;
static CompSequenceBase* CurSeqp=NMRSIM_NULL; //!< sequence under construction

//! reset current sequence to empty state
void reset_sequence() 
{
  CurSeqp=NMRSIM_NULL;
  curchannel=0; // revert to synchronous
  usedchannels=0U;
}

//! ensure we have a current sequence
void ensure_sequence()
{
  if (!CurSeqp)
    CurSeqp=new SyncSequence();
}

void
CompSequenceBase::set(double v, subsid_t subsid) 
{
  switch (subsid) {
  case S_ARG1:
    if (v<0.0)
      error_abort("synchronisation time can't be <0");    
    synctime=v*1e-6;
    break;
  case S_ARG2:
    if (v<0.0)
      error_abort("synchronisation hint can't be <0");    
    synchint=round_int(v);
    break;
  default:
    throw InternalError("CompSequenceBase::set");
  }
}

  /** \c Failed exception if object is dirty */
double CompSequenceBase::duration() const 
{
  if (status_>=DIRTY_ALL)
    throw Failed("CompSequence: info requested while sequence is in undefined state");
  return duration_;
}

double CompSequenceBase::safe_duration()
{
  if (status_>=DIRTY_ALL)
    refresh();
  return duration_;
}
  
std::ostream& operator<< (std::ostream& ostr, dirty_t status)
{
  switch (status) {
  case DIRTY_CLEAN:
    return ostr << "update not required";
  case DIRTY_GAMMA:
    return ostr << "gamma angle changed";
  case DIRTY_HAMILTONIANRF:
    return ostr << "Hamiltonian/RF changed";
  case DIRTY_ALL:
    return ostr << "needs full rebuild";
  }
  throw InternalError("dirty_t");
}

void
CompSequenceBase::print(std::ostream& ostr) const
{
  ostr << name_ << " (" << status_ << ")\n";
  if (status_!=DIRTY_ALL) {
    for (size_t chan=0;chan<seqlist.size();chan++) {
      if (!(seqlist(chan).empty())) {
       if (nchannels>1)
	 ostr << "Channel " << (chan+1) << '\n';
       ostr << seqlist(chan);
      }
    }  
    ostr << "Total duration: ";
    prettyprint_time(duration_) << '\n';
    ostr << "CW: " << (isCW ? "Yes\n" : "No\n");
  }
  ostr << "Synchronise over: ";
  if (synctime)
    ostr << (synctime*1e6) << " us\n";
  else
    ostr << "<unspecified>\n";
  if ((verbose & VER_GEN) && (verbose_level>1))
    ostr << "Partitioning type: " << partition_index << '\n';
}

void
CompSequenceBase::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << name_ << "_synctime";
}

void
CompSequenceBase::parse_pulse(EventID::id_t type)
{
  if (nchannels==0)
    error_abort("can't define pulse sequences without active RF channels (add channels directive to spinsys block)\n");

  char *syncp=get_curline();
  const bool issync=isspace(syncp[1]) && strchr("|+-",*syncp);
  char sync='\0';
  if (issync) {
    sync=*syncp;
    (void)get_token(); //swallow argument
  }
  switch (type) {
   case EventID::automatic:
     if (!issync)
       error_abort("pulseauto: synchronisation must be one of | + -\n");
     break;
  case EventID::soft:
    if (issync)
      error_abort("pulse: synchronisation cannot be specified");
    break;
  default:
    break;
  }
  push_back(EventID::create(this,type,sync)); 
}

void EventID::print(std::ostream& ostr) const
{
  std::cout << "Duration: ";
  if (durs_)
    durs_->print(ostr);
  else
    prettyprint_time(nomdur_);

  ostr << ": ";
  for (size_t chan=0;chan<events.size();chan++) {
    const eventlist_t& evl(events(chan));
    if (!evl.empty()) {
      if (nchannels>1)
	ostr << "Channel " << (chan+1) << ": ";
      for (size_t i=0;i<evl.size();i++)
	ostr << (*evl(i));
    }
  }    
}

std::ostream& operator<< (std::ostream& ostr, const EventID& a)
{
  a.print(ostr);
  return ostr;
}

std::ostream& operator<< (std::ostream& ostr, const CompSequenceBase& a)
{
  a.print(ostr);
  return ostr;
}

int rf_to_subsid(size_t chan, subsid_t subsid)
{
  return (subsid<<3)+chan;
}

void subsid_to_rf(size_t& chan, subsid_t& actsubsid, int subsid)
{
  chan = subsid & 7;
  actsubsid=subsid >> 3;
}

bool RawPulseDef::isconst() const
{
  if (rf && !(rf->isconst()))
    return false;
  if (phase && !(phase->isconst()))
    return false;
  if (offset && !(offset->isconst()))
    return false;
  return true;
}

void EventID::ensureevents()
{
  if (events.empty()) {
    events.create(channels());
    buildevents();
  }
}

void EventID::init(size_t nchans, size_t chanoff_)
{
  chanoff=chanoff_;
  if ((nchans>1) && chanoff_)
    throw InternalError("EventID::init");
  pulses.create(nchans);
}

void EventID::set(double fval_int, size_t chan, subsid_t subsid)
{
  if (isconst())
    error_abort("attempt to change const EventID");

  if (events.empty())
    return;

  eventlist_t& curlist(events(chan));

  const bool doset=(!(curlist.empty()) || (subsid==pulses(chan).master_subsid));
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Waiting for complete event: " << (doset ? "No\n" : "Yes\n");
  if (!doset)
    return;
    
  if (curlist.empty()) {
    buildevents();
    if (!checknull(chan))
      seqp_->flagdirty(DIRTY_ALL); //!< if changing list size, flag that timing will have changed
    return;
  }
    
  for (size_t item=curlist.size();item--;) {
    RFEvent* eventp(curlist(item).get());
    switch (subsid) {
    case S_RF:
      dynamic_cast<RFPulseEvent*>(eventp)->rf(fval_int);
      break;
    case S_PHASE:
      dynamic_cast<RFPulseEvent*>(eventp)->phase(fval_int);
      break;
    case S_OFFSET:
      eventp->offset(convertoffset(fval_int));
      break;
    default:
      throw InternalError("EventID::set");
    }
  }
  seqp_->flagdirty(DIRTY_HAMILTONIANRF);
}

ThreadWarning<> emptypulseevent_warning("setting RF parameter for empty pulse event in sequence: ",&NMRsim_once_warning); 

bool EventID::checknull(size_t chan)
{
  const eventlist_t& curlist(events(chan));
  if (curlist.empty()) {
    emptypulseevent_warning.raise(seqp_->name());
    return true;
  }
  return false;
}

ThreadWarning<> EventID::resizelist_warning("Changing number of elements in a multi-step pulse. This may be inefficient compared to keeping a fixed number of steps, particularly for simple phase-modulated sequences.", &NMRsim_once_warning);

void EventID::set(const BaseList<double>& fval_int, size_t chan, subsid_t subsid)
{
  if (isconst())
    error_abort("attempt to change const EventID");

  if (events.empty())
    return;

  eventlist_t& curlist(events(chan));
  if ((curlist.size()!=fval_int.size()) && !(curlist.empty())) {
    resizelist_warning.raise();
    curlist.clear();
  }

  const bool doset=(!(curlist.empty()) || (subsid==pulses(chan).master_subsid));
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Waiting for complete event: " << (doset ? "No\n" : "Yes\n");
  if (!doset)
    return;

  if (curlist.empty()) {
    buildevents(); //!< note we are assuming that the list argument has come from the VariableBase objects used to build the events
    seqp_->flagdirty(DIRTY_ALL); //!< if changing list size, flag that timing will have changed
    return;
  }

  if (checknull(chan))
    return;
  const bool issoft(type!=ideal);
  for (size_t item=curlist.size();item--;) {
    RFEvent* eventp(curlist(item).get());
    double val=fval_int(item);
    switch (subsid) {
//     case S_DUR:
//       val*=1e-6;
//       if (issoft)
// 	eventp->duration(val);
//       else
// 	eventp->nominal_duration(val);
//       break;
    case S_RF:
      dynamic_cast<RFPulseEvent*>(eventp)->rf(val);
      break;
    case S_PHASE:
      dynamic_cast<RFPulseEvent*>(eventp)->phase(val*deg_to_rad);
      break;
    case S_OFFSET:
      eventp->offset(convertoffset(val));
      break;
    default:
      throw InternalError("EventID::set");
    }
  }
  seqp_->flagdirty((issoft && (subsid==S_DUR)) ? DIRTY_ALL : DIRTY_HAMILTONIANRF);
}

namespace {
  template<typename T> struct tmp {};
  template<> struct tmp<double> { typedef double type; };
  template<typename T> struct tmp< BaseList<T> > { typedef LIST<T> type; };
}
    
// template<class T> void EventID::set_(const T& fval,subsid_t subsid)
// {
//   size_t chan;
//   subsid_t actsubsid;
//   subsid_to_rf(chan,actsubsid,subsid);

//   RawPulseDef& pulsedef(pulses(chan));
//   switch (actsubsid) {
//   case S_RF:
//     pulsedef.rf->set(fval);
//     break;
//   case S_PHASE: {
//     typename tmp<T>::type lval(fval);
//     lval*=deg_to_rad;
//     pulsedef.phase->set(lval);
//   }
//     break;
//   case S_OFFSET:
//     pulsedef.offset->set(fval);;
//     break;
//   default:
//     throw InternalError("pulsedef::set");
//   }
//   set(fval,chan,actsubsid);
// }

void EventID::set(double fval,subsid_t subsid)
{
  if (subsid==S_DUR)
    duration(fval);
  else {
    size_t chan;
    subsid_t actsubsid;
    subsid_to_rf(chan,actsubsid,subsid);
    set(fval,chan,actsubsid);
  }
}

void EventID::set(const BaseList<double>& fval,subsid_t subsid)
{
  if (subsid==S_DUR) {
    for (size_t chan=events.size();chan--;)
      events.clear(); //!< just simply clear everything
    seqp_->flagdirty(DIRTY_ALL);
  }
  else {
    subsid_t actsubsid;
    size_t chan;
    subsid_to_rf(chan,actsubsid,subsid);  
    set(fval,chan,actsubsid);
  }
}

bool EventID::isconst() const
{
  for (size_t i=pulses.size();i--;) {
    if (!(pulses(i).isconst()))
      return false;
  }
  return true;
}

ContextWarning<> short_pulse_warning("pulse duration is <1 ns (time is specified in us)",&NMRsim_once_warning);
ContextWarning<> ignoring_offset_warning("transmitter offsets are ignored for hard pulses",&NMRsim_once_warning);

void checkdur(double dur)
{
  if (dur<0.0)
    error_abort("pulse duration cannot be negative");

  if ((dur<1e-9) && (dur>0.0))
    short_pulse_warning.raise();
}

EventID* EventID::create(CompSequenceBase* seqp,id_t type_, char sync_)
{
  subsid_t basesubsid=S_NONE;
  Variable cvar(S_DUR);
  Mark markobj;
  VariableBase* durp=NMRSIM_NULL;

  const double dur=parse_double(&cvar,F_ALLOWLIST)*1e-6; 
  if (cvar.valuep==NMRSIM_NULL)
    checkdur(dur);
  else {
    durp=&(cvar.variable());
    basesubsid=S_DUR;
  }
  
  const size_t argsleft=count_left(false);
  bool incoffset;
  const size_t usechannels= curchannel ? 1 : nchannels;
  const size_t chanoffv= curchannel ? curchannel-1 : 0;
  if (argsleft==2*usechannels)
    incoffset=false;
  else {
    if (argsleft==3*usechannels) { //may be optional flags
      if (type_==ideal)
	ignoring_offset_warning.raise();
      incoffset=true;
    }
    else
      error_abort("expect 2 or 3 RF arguments per (active) RF channel\n");
  }

  EventID* evidp= new EventID(seqp,usechannels,chanoffv,type_,dur,sync_);
  evidp->durs_=durp;

  static flagsmap_type pulseflags;
  if (pulseflags.empty()) {
    pulseflags["coherent"]=RFPulseEvent::COHERENT;
    pulseflags["transients"]=RFPulseEvent::TRANSIENTS;
  }
  for (size_t chan=0;chan<usechannels;chan++) {
    cvar.subsid=rf_to_subsid(chan,S_RF);
    RawPulseDef& pulsedef(evidp->pulses(chan));
    pulsedef.master_subsid=basesubsid;

    pulsedef.rf=parse_double_variable(cvar,F_DENYZERO | F_ALLOWLIST);
    if (cvar.subsid!=S_NONE)
      pulsedef.master_subsid=S_RF;

    if (pulsedef.rf || cvar.get()) { //non-zero RF or rf is variable
      usedchannels|= 1 << (chan+chanoffv); //!< flag non-zero on channel
      cvar.subsid=rf_to_subsid(chan,S_PHASE);      
      pulsedef.phase=parse_double_variable(cvar,F_ISPHASE | F_ALLOWLIST);
      if (cvar.subsid!=S_NONE)
	pulsedef.master_subsid=S_PHASE;

      if (incoffset) {
	cvar.subsid=rf_to_subsid(chan,S_OFFSET);
	pulsedef.offset=parse_double_variable(cvar,F_ALLOWLIST);
	if (cvar.subsid!=S_NONE)
	  pulsedef.master_subsid=S_OFFSET;
      }
      pulsedef.flags=parse_flags(pulseflags);
      if (!incoffset && (pulsedef.flags & RFPulseEvent::COHERENT))
	error_abort("-coherent flag only valid when offset specified");
      if ((pulsedef.flags & RFPulseEvent::TRANSIENTS) && (type_==ideal))	
	error_abort("-transients cannot be specified for ideal pulses");
    }
    else { //0 RF
      pulsedef.rf=NMRSIM_NULL;
      parse_double(); //swallow phase
      if (incoffset)
	parse_double(); //swallow offset
    }
  }
  // if (evidp->isconst()) 
  //  evidp->buildevents(); //!< only build events if EventID is const (otherwise will be constructed later)
  markobj.flush(evidp);
  return evidp;
}

void EventID::printvariablename(std::ostream& ostr, subsid_t subsid) const
{  
  if (subsid==S_DUR) {
    ostr << "duration"; //NB this will not be unique
    return;
  }
  size_t chan;
  subsid_t actsubsid;
  subsid_to_rf(chan,actsubsid,subsid);
  if (nchannels>1)
    ostr << "Channel " << (chan+1) << ' ';
  switch (actsubsid) {
  case S_RF: ostr << "nutation rate"; break;
  case S_OFFSET: ostr << "offset"; break;
  case S_PHASE: ostr << "phase"; break;
  default: throw InternalError("Unknown RF subsid type");
  }
}

void parse_pulse(int which)
{
  ensure_sequence();
  CurSeqp->parse_pulse(EventID::id_t(which));
}
  
//! call after each use (increments if reqd)
void CycledSequence::next()
{
  const size_t curindex=currentindex();
  sequencestate_t& curstate(states_(curindex));
  if (++(curstate.curptr)==size(curindex))
    curstate.curptr=0;
  curstate.lastoffset_=0.0; //!< clear offset
}

// void sequencestate_t::lastoffset(double v)
// {
//   sequencestate_t& curstate(currentstate());
//   if ((v<0.0) || (v>=current()->duration()))
//     throw InvalidParameter("CycledSequence::lastoffset");
//   curstate.lastoffset_=v;
// }

ThreadWarning<> CycledSequence::nonuniform_warning("components of cycled list have unequal durations.  Sure about this?",&NMRsim_once_warning);

double CycledSequence::duration() const
{
  const sequencestate_t& curstate(currentstate());
  if (curstate.duration_<0.0) {
    const BaseList<PhasedSequence*> asrow(currentlist());
    const double durstep=asrow.front()->duration();
    curstate.duration_=durstep;
    for (size_t i=1;i<asrow.size();i++) {
      const double ldurstep=asrow(i)->duration();
      if (!nochecks && (fabs(ldurstep-durstep)>tolerance))
	nonuniform_warning.raise(name());
      curstate.duration_+=ldurstep;
    }
  }
  return curstate.duration_;
}
  
void CycledSequence::set(double phase, size_t ind)
{
  BaseList<PhasedSequence*> asrow(row());
  if (ind>=asrow.size()) {
    std::cerr << "Error: trying to set phase for sequence fragment " << (ind+1) << " when there are only " << asrow.size() << " fragment(s)!\n";
    error_abort();
  }
  asrow(ind)->set(phase);
}
  
  void PhasedSequence::set(double fval,subsid_t subsid) 
  {
    assert(subsid==S_PHASE); //!< assert is harmless
    phase=fval;
  }

void PhasedSequence::print(std::ostream& ostr) const
{
  if (phase)
    ostr << phase << '+';
  ostr << seq.name();
}

void PhasedSequence::printvariablename(std::ostream& ostr, subsid_t) const
{
  ostr << seq.name() << "_phaseshift";
}

void parse_delay()
{
  ensure_sequence();
  Variable cvar(S_DUR,NMRSIM_NULL);
  Mark markobj;
  const double dt=parse_double(&cvar)*1e-6;
  EventID* evidp=new EventID(CurSeqp,dt);
  CurSeqp->push_back(evidp);
  if (cvar.subsid!=S_NONE)
    markobj.flush(evidp);
}

ThreadWarning<> chebyshev_warning("Chebyshev propagation unlikely to work for time-independent Hamiltonians",&NMRsim_once_warning);

void ensure_partition(size_t ind)
{
  matrix_partition_set& partset(partitions(ind));
  if (!partset.empty() || !optpartition)
    return; //!< already built or not using partitions?
  
  for (size_t chan=nchannels;chan--;) { //!< determine partitioning for active channels + Hamiltonian
    if (ind & (1<<chan))
      partset|=pgenstack(chan)->partitioning();
  }
  partset.adddiagonal();
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Created partitioning for type " << ind << ":\n" << partset << '\n';

  bool use_cheby=false;
  switch (prop_method) {
//   case PROP_AUTO:
//     if (verbose & VER_GEN)
//       std::cout << "Density for partition " << ind << ": " << (100.0*partset.density()) << "\%\n";
//     use_cheby=(partset.density()<chebyshev_crossover);    
//     break;
  case PROP_DIAG:
    use_cheby=false;
    break;
  case PROP_CHEBYSHEV:
    use_cheby=true;
    break;
  }

  if (verbose & VER_GEN)
    std::cout << "Use Chebyshev for partition " << ind << ": " << (use_cheby ? "Yes\n" : "No\n");
  
  if (use_cheby) {
    if (spin_rate==0.0)
      chebyshev_warning.raise();
    prop_flags(ind)|=BaseMetaPropagator::usechebyshev;
  }

}

ContextWarning<> empty_store_warning("store applied to empty sequence",&NMRsim_repeat_warning);
  
void parse_store()
{
  const char* seqname(parse_string(F_REPLACEDOLLAR));
  if (cycledseqmap.count(seqname) || seqmap.count(seqname)) {
    parser_printcontext() << "can't redefine sequence/sequence list " << seqname << '\n';
    error_abort();
  }
  havestore=true;

  double synctime=0.0;
  int synchint=0;
  Mark markobj;
  if (are_left()) {
    char* ptr=parse_string(F_IGNOREDOLLAR | F_ALLOWMISSING); //!< don't use REPLACEDOLLAR here, as argument can also be numeric
    CycledSequence* cseqp=parse_cycledsequence(ptr,F_ALLOWMISSING);
    if (cseqp) {
      cseqp->name(seqname);
      cycledseqmap[seqname]=cseqp;
      return;
    }
    Variable cvar(S_ARG1,CurSeqp);
    synctime=parse_double_raw(&cvar,ptr,0); //!< no scaling as this is done in set()

    if (are_left()) {
      cvar.subsid=S_ARG2;
      synchint=parse_unsigned(&cvar);
    }
  }
  if (!CurSeqp) {
    ensure_sequence();
    empty_store_warning.raise();
  }
  else
    CurSeqp->flagdirty(DIRTY_ALL);
  CurSeqp->partition_index=usedchannels;
  CurSeqp->name(seqname);
  CurSeqp->set(synctime,S_ARG1);
  CurSeqp->set(synchint,S_ARG2);
  seqmap[seqname]=CurSeqp;
  markobj.flush(CurSeqp);
  reset_sequence();
}

void CompSequenceBase::clear() 
{
  Ucache.clear();
  doneerror=false;
  seqlist.clear();
  duration_=0.0;
  status_=DIRTY_CLEAN;
  isCW=true;
}

void SyncSequence::clear() 
{
  rawrep.clear();
  CompSequenceBase::clear();
}

void AsyncSequence::clear()
{
  for (size_t i=rawrep.size();i--;)
    rawrep(i).clear();
  CompSequenceBase::clear();
}

void AsyncSequence::push_back(EventID* evlist)
{
  if (evlist->channels()!=1)
    throw Failed("AsyncSequence::push_back: multi-channel sequence element");
  if (evlist->chanoff!=curchannel-1)
    throw Failed("AsyncSequence::push_back: sequence fragments are defined for different nuclei");
  rawrep(curchannel-1).push_back(evlist);
}

void AsyncSequence::push_back(size_t chan, const BaseList<EventID*>& evl, size_t ntimes)
{
  LIST<EventID*> dest(rawrep(chan));
  const size_t len(evl.size());
  dest.reserve(dest.size()+len*ntimes);
  for (;ntimes--;)
    for (size_t i=0;i<len;i++)
      dest.push_back(evl(i));
}

// void AsyncSequence::push_back(CompSequenceBase* propseqp, size_t ntimes)
// {
//   AsyncSequence* seqp=dynamic_cast<AsyncSequence*>(propseqp);
//   if (!seqp)
//     throw Failed("Incompatible sequence type");
//   propseqp->links.insert(this);
//   if (curchannel)
//     push_back(curchannel-1,seqp->rawrep(curchannel-1),ntimes);
//   else {
//     for (size_t chan=nchannels;chan--;)
//       push_back(chan,seqp->rawrep(chan),ntimes);
//   }
// }

// void SyncSequence::push_back(CompSequenceBase* propseqp, size_t ntimes)
// {
//   const SyncSequence* seqp=dynamic_cast<SyncSequence*>(propseqp);
//   if (!seqp)
//     throw Failed("Incompatible sequence type");
//   if (seqp==this)
//     throw InternalError("SyncSequence::push_back");
//   propseqp->links.insert(this);
//   const BaseList<EventID*>& evl(seqp->rawrep);
//   const size_t len(evl.size());
//   rawrep.reserve(rawrep.size()+len*ntimes);
//   for (;ntimes--;)
//     for (size_t i=0;i<len;i++)
//       push_back(evl(i));
// }

AsyncSequence::AsyncSequence()
{
  rawrep.create(nchannels);
}

void AsyncSequence::flush_delay(size_t chan, double dt)
{
  if (dt)
    times(chan)+=seqlist(chan).push_delay(times(chan),dt);
}

void SyncSequence::flush_delay(double dt)
{
  if (dt) {
    for (size_t chan=nchannels;chan--;)
      times(chan)+=seqlist(chan).push_delay(times(chan),dt);
  }
}

void CycledSequence::printshort(std::ostream& ostr) const
{
  const char* namep=name();
  if (namep[0])
    ostr << namep;
  else
    print(ostr);
}

// double CycledSequence::synchronisation_time() const
// {
//   const BaseList<PhasedSequence*> asrow(this->row());
//   double synct=0.0;
//   double synct=asrow.front()->synchronisation_time();
//   for (size_t i=<size();i--;) {
//     const double lsync=asrow(i)->sychronisation_time();
//     if (lsync) {
//       if (synct) {
// 	if (fabs(lsync-synct)>1e-6) {
// 	  std::cerr << "CycledSequence " << name << ": incompatible synchronisation times\n";
	  
//       error_abort(ERR_FAILED);
//     }
//   }
//   return dur;
// }

PhasedSequence* CycledSequence::current() const
{ 
  const size_t curindex=currentindex();
  return (*this)(curindex,states_(curindex).curptr);
}

size_t CycledSequence::currentindex() const
{
  size_t transindex;
  (void)getindex(transindex,arraytag_,size());
  return transindex;
}

// PhasedSequence* CycledSequence::current() const
// { 
//   size_t which=0;
//   if (size()>1) {
//     if (sum_index<0)
//       throw Failed("sum variable used outside valid region");
//     which=(sumtag_ ? suminds(sumtag_-1) : sum_index) % size();
//   }
//   const BaseList<PhasedSequence*> lrow((*this)(which));
//   if (isspecial_)
//     return lrow(curptr);
//   if (lrow.size()==1)
//     return lrow.front();
//   size_t transindex;
//   (void)getindex(transindex,arraytags_(which),lrow.size());
//   //  const size_t i=tag ? varinds(tag-1) : var_index;
//   //  return lrow(i % lrow.size());
//   return lrow(transindex);
// }

std::ostream& operator<< (std::ostream& ostr, const PhasedSequence& a)
{
  a.print(ostr);
  return ostr;
}

ThreadWarning<> CWoffres_warning("Sequence with continuous off-resonance irradiation detected.  This may behave strangely.  Consider whether -coherent flag is appropriate for: ",&NMRsim_once_warning);

void CompSequenceBase::rebuild()
{
  if (status_<DIRTY_ALL)
    return;

  //  masterobjp->ensure_channels();
  ensure_partition(partition_index);

  duration_=raw_rebuild();

  isCW=true;
  bool CWoffres=false;
  for (size_t n=0;n<seqlist.size();n++) {
    if ((verbose & VER_GEN) && (verbose_level>1))
      std::cout << "Overall duration on channel " << (n+1) << ": " << (1e6*times(n)) << " us\n";
    Sequence& curseq(seqlist(n));
    if (curseq.empty() && (times(n)==0.0))
      continue;
    curseq.autosync(tolerance,times(n)-tolerance);
    if (n==0)
      duration_=times(n);
    else {
      if (fabs(times(n)-duration_)>1e-8)
	throw Failed("Sync failure: sequences on different channels have different durations");
    }
    if ((curseq.size()>1) || ((curseq.size()==1) && (curseq.front().duration()!=duration_)))
     //    if (!curseq.empty() && ((curseq.size()>1) || (curseq.front().duration()!=duration_)))
      isCW=false;
    else {
      const RFPulseEvent* cweventp=dynamic_cast< const RFPulseEvent* >(curseq.front().event());
      if (cweventp && cweventp->offset() && !((cweventp->flags()) && RFPulseEvent::COHERENT))
	CWoffres=true;
    }
  }
  //  store_transients(duration_);
  
  reset(DIRTY_ALL);

  status_=DIRTY_CLEAN;
  if (verbose & VER_GEN) {
    std::cout << "Reconstructed sequence: " << name() << '\n' << (*this) << '\n';
    if (verbose_level>1)
      std::cout << "Overall duration: " << (1e6*duration_) << " us\n";
  }
  if (isCW && CWoffres)
    CWoffres_warning.raise(name(),true);
}

/** 
    \post \c status_==DIRTY_CLEAN
*/
void CompSequenceBase::refresh()
{
  switch (status_) {
  case DIRTY_CLEAN:
    return;
  case DIRTY_ALL:
    rebuild();
    break;
  default:
    reset(status_);    
  }
  status_=DIRTY_CLEAN;
}
    
void CompSequenceBase::reset(dirty_t level)
{
  checkpm(level);
  Ucache.clear();
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Clearing propagator cache for " << name() << '\n';
}

bool CompSequenceBase::needsresync() const
{
  return ((!!phasemodp) && (!*phasemodp)) || wasdirty() || !optupdate; //!< \c true if PhaseModulation exists but needs 'build'
}

ThreadWarning<> phasemodulation_failed_warning("NB. Fast MAS+RF propagation could not be applied to sequence (verbose -general gives more information): ",&NMRsim_once_warning,BaseWarning::Inherit,std::cout); //!< should really be sequence-specific

void CompSequenceBase::checkpm(dirty_t level)  
{
  const int lverb=(verbose & VER_GEN) ? verbose_level : 0;
  const bool needrefresh=(level>=DIRTY_HAMILTONIANRF);
  if (lverb>1)
    std::cout << "PhaseModulation: needs refresh: " << (needrefresh ? "Yes\n" : "No\n");
  if (needrefresh) {
    phasemodp.clear();    
    propgenp.clear(); //!< should be ignored anyway, but no harm done to be sure
    if (spin_rate) { //!< only bother if spinning
      double rf_period=duration_; //!< can't call duration() as dirty status may not have been cleared
      if (optphasemod() && rf_period) {
	if (isCW) {
	  double new_rf_period=1.0/spin_rate; //!< set RF period to match rotor period
	  if (new_rf_period<rf_period)
	    rf_period=new_rf_period; //avoid changing to longer than sequence duration, otherwise won't wrap round
	}
	if (lverb>1) {
	  std::cout << "PhaseModulation: checking whether phase modulation is possible for sequence " << name() << '\n';
	  std::cout << "PhaseModulation: effective RF period: " << (rf_period*1e6) << " us\n";
	}
		
	try {
	  //    if (seqlist.size()!=1)
	  //  throw Failed("PhaseModulation: not single channel sequence");
	  phasemodp.reset(new PhaseModulation(seqlist,rf_period,tolerance,lverb));
	  if (lverb>1)
	    std::cout << "PhaseModulation: created new PhaseModulation object\n";	  
	  if (isCW && !(phasemodp->isCW())) {
	    std::cerr << "Internal error: PhaseModulation object created from CW sequence " << name() << " is not itself CW!\n";
	    error_abort();
	  }
	}
	catch (MatrixException& exc) {
	  if (lverb)
	    std::cout << "Sequence " << name() << " could not be reduced to phase modulation: " << exc << '\n';
	  else
	    phasemodulation_failed_warning.raise(name());

	  phasemodp.clear(); //!< failed to create ::PhaseModulation
	}	  
      }
      char buf[256];
      snprintf(buf,sizeof(buf)-1,"phase modulated optimisation for sequence %s",name());
      optphasemod.check( (phasemodp.get()!=NMRSIM_NULL), buf);
    }
  }
  if (level)
    propgenp.clear(); //!< clear propagator generator
}

ThreadWarning<> phasemodulation_timedivision_warning("phase modulation: cannot further sub-divide time base within synchronisation period",&NMRsim_once_warning);

namespace {

  bool addsync(size_t& nsync, size_t n, const char* name, double maxrat)
  {
    switch (n) {
    case 0:
      throw InternalError("addsync");
    case 1:
      return true;
    }
    const size_t newsync=lcm(nsync,n);
    if (newsync>maxrat) {
      if (phasemodulation_timedivision_warning.enabled()) {
	char buf[256];
	snprintf(buf,sizeof(buf)-1," to %lu to accommodate %s sampling. Limit set by maxdt is %i steps.",(long unsigned)newsync,name,(int)(maxrat+0.5));
	phasemodulation_timedivision_warning.raise(buf);
      }
      return false;
    }
    nsync=newsync;
    return true;
  }
}

void CycledSequence::checksync(double dur)
{
  const bool isverb=((verbose & VER_GEN) && (verbose_level>1));
  std::set<CompSequenceBase*> doneset;
  BaseList<PhasedSequence*> rawseqs(this->row());
  for (size_t i=rawseqs.size();i--;) {
    CompSequenceBase& curseq(rawseqs(i)->get());
    const bool rebuild=(doneset.count(&curseq)==0);
    if (isverb)
      std::cout << "Checking synchronisation for " << curseq.name() << " over duration " << (1e6*dur) << " us: " << (rebuild ? "Yes\n" : "No\n");
    if (rebuild) {
      curseq.checksync(dur);
      doneset.insert(&curseq);
    }
  }
}

namespace {
  void ltrysync(double& synctime, const char*& syncnamep, double period, const char* name)
  {
    if (period>synctime) {
      synctime=period;
      syncnamep=name;
    }
  }
}

ThreadWarning<> phasemodulation_notused_warning("phase modulation disabled due to lack of synchronisation",&NMRsim_once_warning);
ThreadWarning<> phasemodulation_smalltimebase_warning("phase modulation base time is much less than integration step",&NMRsim_once_warning);
ThreadWarning<> phasemodulation_syncconflict("sequence seems to be shared with conflicting synchronisation conditions - phase modulation disabled: ",&NMRsim_once_warning);

void CompSequenceBase::checksync(double dt)
{
  if (!phasemodp || (dt==0.0))
    return;

  if (spin_rate==0.0)
    throw InternalError("checksync called for static simulation");
  const double rotor_period=1.0/fabs(spin_rate);
  const double rf_period=duration();
  const char* syncnamep="synchronisation time";
  double usesync=synctime;

  if (synctime==0) {
    ltrysync(usesync,syncnamep,rotor_period,"MAS period");
    if (!isCW)
      ltrysync(usesync,syncnamep,rf_period,"RF period");
    ltrysync(usesync,syncnamep,dt,"observation period");
  }
    
  const size_t nR=verbose_check_sync(usesync,syncnamep,rotor_period,"MAS period");
  const size_t nC=isCW ? 1 : verbose_check_sync(usesync,syncnamep,rf_period,"RF period");
  if ((nR==0) || (nC==0)) {
    phasemodulation_notused_warning.raise();
    phasemodp.clear();
    return;
  } 

  const double maxdt=get_maxdt();
  size_t nmult=1;
  if (maxjumpdt==0.0) {
    const size_t nsteps=isCW ? 1 : nC*phasemodp->timesteps();
    size_t nsync=lcm(nR,nsteps);

    const double maxrat=usesync/(NMRSIM_MINSTEP*maxdt); // allow step to drop below maxdt without becoming excessively small

    bool addsyncok=true;

    if (dt<usesync) {
      const size_t ndt=verbose_check_sync(usesync,syncnamep,dt,"dwell time");
      if (ndt)
	addsyncok &= addsync(nsync,ndt,"dwell time",maxrat); // try to accommodate dwell time first
    }
    if (gamma_angles>0)
      addsyncok &= addsync(nsync,nR*gamma_angles,"gamma angles",maxrat);

    nmult = nsync / nsteps;
    if (!addsyncok && !(optphasemod.isenabled())) {
      phasemodulation_notused_warning.raise();
      phasemodp.clear();
      return;
    }
  } 
      
//   if (ndt)
//     nsync=lcm(nsync,ndt);
//   if (synctime) { //!< synctime is zero for CW i.e. always synced
//     const size_t nrotor=check_sync(synctime/rotor_period);
//     if (nrotor)
//       nsync=lcm(nsync,nrotor);
//     else {
//       if (!nochecks)
// 	std::cerr << "Warning: sequence synchronisation time (" << (synctime*1e6) << " us) is not multiple of rotor period (" << rotor_period << " us) and is being ignored\n";
//     }
//   }

  const int lverb=(verbose & VER_GEN) ? verbose_level : 0;
  const size_t curmult=!(*phasemodp) ? 0 : phasemodp->multiplier();
  if (curmult && (lastusedmult!=nmult)) {
    phasemodulation_syncconflict.raise(name(),true);
    phasemodp.clear();
    return;
  }
  lastusedmult=nmult;
  if (curmult && ((curmult % nmult)==0)) {
    if (lverb)
      std::cout << "No need to rebuild phase modulation with factor of " << nmult << " (divides into existing factor)\n";
  }
  else {    
    phasemodp->build(nmult);
    if (lverb) {
      std::cout << "Built phase modulation subdividing characteristic time into " << nmult << " steps\n";
      if (lverb>1)
	std::cout << (*phasemodp) << '\n';
    }
    if (phasemodp->tick() && (phasemodp->tick()<NMRSIM_MINSTEP*maxdt) && phasemodulation_smalltimebase_warning.enabled()) {
      std::ostringstream str(std::ostringstream::out);
      str << " (";
      prettyprint_time(phasemodp->tick(),str) << " compared to " << (1e6*maxdt) << " us)";
      if (nmult>1)
	str << ".\nCalculation efficiency may be improved by specifying -disable:phasemodulation flag or by allowing the timing to vary using maxjumpdt.";
      phasemodulation_smalltimebase_warning.raise(str.str().c_str());
    }
  }
}
 
bool CompSequenceBase::wasdirty() const
{
  return dirty_list.count(const_cast<CompSequenceBase*>(this));
}

void rebuild_sequences(dirty_t overalldirt)
{
  seqmap_type::iterator curp=seqmap.begin();
  const seqmap_type::iterator endp=seqmap.end();
  dirty_list.clear();
  if (curp==endp)
    return;

  bool setdirt=(overalldirt!=DIRTY_CLEAN);
  while (curp!=endp) {
    CompSequenceBase* curseqp(curp->second);
    if (setdirt)
      curseqp->flagdirty(overalldirt);
    if (curseqp->isdirty()) {
      dirty_list.insert(curseqp);
      curseqp->refresh();
    }
    ++curp;
  }
  
  const int lverb=(verbose & VER_GEN) ? verbose_level : 0;

  if (dirty_list.empty()) {
    if (lverb)
      std::cout << "No sequences required rebuilding\n";
  }
  else {
    if (lverb) {
      std::cout << "Sequence fragments that have been rebuilt:";
      const dirtylist_t::const_iterator end(dirty_list.end());
      dirtylist_t::const_iterator start(dirty_list.begin());
      while (start!=end) {
	std::cout << ' ' << (*start)->name();
	++start;
      }
      std::cout << '\n';
    }
    //scan through cycled lists
    for (size_t i=0;i<cycledseqlist.size();i++) {
      CycledSequence& cseq(*(cycledseqlist(i)));
      cseq.checkdirty();
      if (lverb>1)
	std::cout << "Status of cycled sequence " << cseq.name() << ": " << (cseq.wasdirty() ? "needs flush\n" : "clean\n");
    }
  }
}

const char* CycledSequence::prettyname() const
{
  const char* ptr=name();
  return *ptr ? ptr : "<anonymous>";
}

cacheset_t& sequencestate_t::setcache(double offsetv)
{
  if (offsetv==0.0)
    return cachezero_;
  cacheoffset_=offsetv;
  return cachenonzero_;
}

cacheset_t* sequencestate_t::getcache(double offsetv)
{
  if (offsetv==0.0)
    return &cachezero_;
  if ((cacheoffset_>0.0) && (fabs(offsetv-cacheoffset_)<tolerance))
    return &cachenonzero_;
  return NMRSIM_NULL;
}
  
sequencestate_t::sequencestate_t()
{
  cacheoffset_=-1.0;
  reset(true);
}

void sequencestate_t::flush()
{
  cachezero_.flush();
  cachenonzero_.flush();
  cacheoffset_=-1.0;
  duration_=-1.0; //!< flag duration needs recalculating
}

void CycledSequence::flush()
{
  for (size_t i=states_.size();i--;)
    states_(i).flush();
}

namespace {
  void dumparray(std::ostream& ostr, const BaseList<PhasedSequence*>& a)
  {
    if (a.size()==1) {
      ostr << *(a.front());
      return;
    }
    ostr << '[';
    for (size_t i=0;i<a.size();i++) {
      if (i)
	ostr << ',';
      ostr << *(a(i));
    }
    ostr << ']';
//     if (tag)
//       ostr << ':' << tag;
  }
}

void CycledSequence::print(std::ostream& ostr) const
{
  if (size()==1)
    dumparray(ostr,front());
  //	      arraytags_.empty() ? 0 : arraytags_.front());
  else {
    ostr << '{';
    for (size_t i=0;i<size();i++) {
      if (i)
	ostr << ',';
      dumparray(ostr,(*this)(i));
      //		arraytags_.empty() ? 0 : arraytags_(i));
    }
    ostr << '}';
    if (arraytag_)
      ostr << ':' << arraytag_;
  }
//   ostr << "  Total duration: ";
//   if (duration_<0.0)
//     ostr << "<undetermined>";
//   else
//     prettyprint_time(duration_,ostr);
}

bool CycledSequence::hascycle() const
{
  for (size_t i=size();i--;) {
    if (size(i)>1)
      return true;
  }
  return false;
}

bool CycledSequence::hasCW() const
{
  for (size_t i=size();i--;) {
    if ((size(i)==1) && (*this)(i).front()->get().isCW)      
      return true;
  }
  return false;
}

void CycledSequence::checkdirty()
{
  LIST<PhasedSequence*> aslist(row());
  dirty_=false;
  for (size_t i=aslist.size();i--;) {
    if (aslist(i)->isdirty())
      dirty_=true;
  }
  if ((verbose_level>1) && (verbose & VER_GEN))
    std::cout << "Flagging composite sequence " << prettyname() << " as dirty: " << (dirty_ ? "Yes\n" : "No\n");
  if (dirty_)
    flush();
}

bool PhasedSequence::isdirty() const
{
  return (dirty_list.count(&seq)!=0);
}

// void flagdirty(dirty_t status)
// {
//   seqmap_type::iterator curp=seqmap.begin();
//   seqmap_type::iterator endp=seqmap.end();
//   while (curp!=endp) {
//     curp->second->flagdirty(status,false);
//     ++curp;
//   }
// }

ThreadWarning<> ignoring_transients_warning("-transients specified when transients not set (ignored)",&NMRsim_once_warning);

//int EventID::size(bool allownull) const
int EventID::size() const
{
  if (isdelay())
	return 1; // simple delays have effectively one pulse
  accumulator acclength;
  if (durs_)
    acclength.add(durs_->value().size());

  for (size_t chan=pulses.size();chan--;) {
    const RawPulseDef& pulsedef(pulses(chan));
    if (!pulsedef.rf) //no RF
      continue;
    const size_t amplen=pulsedef.rf->value().size();
    if (!amplen && !(pulsedef.rf->isconst()))
      return -1; 
    const size_t phaselen=pulsedef.phase->value().size();
    if (!phaselen && !(pulsedef.phase->isconst()))
      return -1; 
//     if ((!amplen || !phaselen) && allownull)
//       return 0;
    bool ok = acclength.add(amplen) && acclength.add(phaselen);
    if (pulsedef.offset) {
      const size_t offsetlen=pulsedef.offset->value().size();
    if (!offsetlen && !(pulsedef.offset->isconst()))
      return -1; 
    //      if (!offsetlen && allownull)
    //	return 0;
      ok &= acclength.add(offsetlen);
    }
    if (!ok)
      error_abort("Different parameters of a pulse definition have incompatible list sizes");
  }
  if (!acclength) // allowed to return -1 to indicate size not known
    return -1;
  return acclength();
}

bool EventID::buildevents()
{ 
  const int common_length=size();
//   if (allowfail && (common_length==0))
//     return false;
  const bool timevar=(durs_ && (durs_->value().size()>1));

  const size_t usechannels=pulses.size();
  events.create(usechannels);

  if (common_length<=0) {
    for (size_t i=usechannels;i--;)
      events(i).clear();
    return false;
  }

  //  syncs.create(nchannels);
  bool dirty=false;
  if ((common_length>1) && sync && (sync!='|')) {
    sync='\0';
    syncwithlist_warning.raise();
  }
  for (size_t chan=usechannels;chan--;) {
    const RawPulseDef& pulsedef(pulses(chan));
    if (!pulsedef.rf) //no RF
      continue;
    eventlist_t& dest(events(chan));  
    const BaseList<double>& rflist(pulsedef.rf->value());
    const BaseList<double>& phaselist(pulsedef.phase->value());
    bool simple=(rflist.size()==1) && (phaselist.size()==1);
    BaseList<double> offsetlist;
    if (pulsedef.offset) {
      offsetlist.create(pulsedef.offset->value());
      simple&=(offsetlist.size()==1);
    }
    //    bool buildtransients=(pulsedef.flags & RFPulseEvent::TRANSIENTS) || (transient_model==TRANS_SIMPLE);
    if (common_length==dest.size()) //no change
      continue;
    const size_t actchan=chanoff+chan;    
    int flags=pulsedef.flags;
    switch (transient_model) {
    case TRANS_NONE:
      if (flags & RFPulseEvent::TRANSIENTS)
	ignoring_transients_warning.raise();
      break;
    case TRANS_MANUAL:
      break;
    case TRANS_AUTO:
      if (!transient_amps(actchan)->isconstant() || (transient_amps(actchan)->get()))
	flags|=RFPulseEvent::TRANSIENTS;
      break;
    }
    dirty=true;
    dest.create(common_length);
    const bool rfvar(rflist.size()>1);
    const bool phasevar(phaselist.size()>1);
    const bool offsetvar(offsetlist.size()>1);
    double rf=rflist.front();
    double phase=phaselist.front()*deg_to_rad;
    double offset=offsetlist.empty() ? 0.0 : offsetlist.front();
    PulseGeneratorBase& pgen=*(pgenstack(actchan));

    RFEvent* newpulse=NMRSIM_NULL;
    for (size_t i=0;i<common_length;i++) {
      if (!simple || (i==0)) {
	const double dur=timevar ? durs_->value()(i)*1e-6 : nomdur_;
	if (rfvar)
	  rf=rflist(i);
	if (phasevar) 
	  phase=phaselist(i)*deg_to_rad;
	if (offsetvar)
	  offset=offsetlist(i);
	
	if (type==ideal)
	  newpulse=new HardPulse(pgen,dur,rf,phase);
	else
	  newpulse=new CWPulse(pgen,dur,rf,phase,convertoffset(offset),flags);
      }
      dest(i).reset(newpulse);
    }
  }
  if (dirty) { //flag that sequences will need rebuilding from new events
//     if ((verbose & VER_GEN) && (verbose_level>1))
//       std::cout << "Flag rebuild for sequence: " << seqp_->name << '\n';
    seqp_->flagdirty(DIRTY_ALL);
  }
  return true;
}
 
void EventID::duration(double dur)
{
  dur*=1e-6;
  if (nomdur_==dur)
    return;

  const bool issoft(type!=ideal);
  for (size_t chan=events.size();chan--;) {
    eventlist_t& curlist(events(chan));
    for (size_t item=curlist.size();item--;) {
      RFEvent* evp(curlist(item).get());
      if (issoft)
	evp->duration(dur);
      else
	evp->nominal_duration(dur);
    }
  }
  nomdur_=dur;
  seqp_->flagdirty(issoft ? DIRTY_ALL : DIRTY_HAMILTONIANRF);
}

namespace {
  char cleansync(char sync, size_t n, size_t nels)
  {
    if (sync!='x')
      return sync;
    if (n==0)
      return '+';
    if (n==nels-1)      
      return '-';
    return '|';
  }
}

double EventID::inset() const
{
  const double dur=totalduration();
  if ((dur<0) && (sync!='+'))
    error_abort("Can't use - or | synchronisation when pulse fragment duration is undetermined");

  switch (sync) {
  case '+':
    return 0.0;
  case '-':
    return dur;
    break;
  case '|':
    return dur/2;
  }
  throw InternalError("EventID::inset");
}

ThreadWarning<> syncwithlist_warning("synchronisation specified with list of RF events",&NMRsim_repeat_warning);

ThreadWarning<> nullsequencelist_warning("RF sequence will do nothing: ",&NMRsim_once_warning);

namespace {
  char getsync(const EventID& evid, size_t chan =0)
  {
    char sync=evid.sync;
    switch (sync) {
    case '|': break;
    case '\0':
      sync='|';
      break;
    default:
      if (evid.events(chan).size()>1)
	syncwithlist_warning.raise();
    }
    return sync;
  }
}

double SyncSequence::raw_rebuild()
{
  seqlist.clear();
  if (nchannels) { //allow nchannels to be zero
    seqlist.create(nchannels,Sequence(spectrometer));
    times.create(nchannels);
    times=0;
  }
  const bool isverb=(verbose & VER_GEN) && (verbose_level>1);
  double totdelay=0.0;
  double delay=0.0;
  double needdelay=0.0;
  bool hasduration=false;
  bool foundevents=false;
  for (size_t n=0;n<rawrep.size();n++) {
    EventID* evidp(rawrep(n));
    double dur=evidp->totalduration();
    if (dur>tolerance)
      hasduration=true;
    if (dur<0) { //!< if undetermined, flag that sequence will exist, but use nominal duration of 0
      dur=0;
      hasduration=true;
    }
    if (evidp->isdelay()) {
      flush_delay(delay);
      delay=dur-needdelay;
      totdelay+=delay;
      if (delay<0)
	throw Failed("Insufficient delay period after +/| synchronised auto pulse");
      needdelay=0.0;
    }
    else {      
      evidp->ensureevents();
      if (evidp->type==EventID::automatic) {
	const double cor=evidp->inset();
	if (delay<cor)
	  throw Failed("auto event overlaps with previous event / start of sequence");
	needdelay=dur-cor;
	delay-=cor;
      }
      flush_delay(delay);
      delay=0.0;
      for (size_t chan=0;chan<nchannels;chan++) {
	if (isverb && (nchannels>1))
	  std::cout << "Channel " << (chan+1) << '\n';
	eventlist_t& curlist(evidp->events(chan));
	Sequence& curseq(seqlist(chan));
	//	const BaseList<char> cursync(evidp->syncs(chan));
	double& curtime(times(chan));
	const char sync=getsync(*evidp,chan);
	for (size_t item=0;item<curlist.size();item++) {
	  RFEvent* p(curlist(item).get());
	  const TimedEvent& tev(curseq.push_back(p,curtime,sync));
	  if (isverb)
	    std::cout << "Added: " << tev;
	  foundevents=true;
	  curtime+=tev.duration();
	}
      }
    }
  }
  if (needdelay)
    throw Failed("+/| synchronised auto pulse at end of sequence fragment");
  flush_delay(delay);
  if (!nochecks && (verbose & VER_GEN) && !foundevents && (delay==0.0) && !hasduration) //!< made conditional on verbose as this warning is often bogus
    nullsequencelist_warning.raise(name(),true);
  return totdelay;
}

// ought to be cached really
double EventID::totalduration() const
{
  if (type==ideal)
    return 0.0;

  if (durs_==NMRSIM_NULL) {
    const int npulses=size();
    if (npulses<0)
      return -1.0; //!< flag duration unknown
    return npulses*nomdur_;
  }

  const BaseList<double> times(durs_->value());
  double total=0.0;
  for (size_t i=times.size();i--;)
    total+=times(i);
  return total;
}

double AsyncSequence::raw_rebuild()
{
  seqlist.clear();
  seqlist.create(nchannels,Sequence(spectrometer));
  times.create(nchannels);
  times=0;
  bool foundevents=false;
  double delay;
  bool hasduration=false;
  for (size_t chan=nchannels;chan--;) {
    const LIST<EventID*>& rawrepchan(rawrep(chan));
    if (rawrepchan.empty())
      continue;
    double& timeschan(times(chan));
    Sequence& curseq(seqlist(chan));

    const size_t nels(rawrepchan.size());
    delay=0.0;
    double needdelay=0.0;      
    for (size_t n=0;n<nels;n++) {
      EventID* evidp(rawrepchan(n));
      evidp->ensureevents();
      double dur=evidp->totalduration();
      if (dur>tolerance)
	hasduration=true;
      if (dur<0) { //!< if undetermined, flag that sequence will exist, but use nominal duration of 0
	dur=0;
	hasduration=true;
      }
      switch (evidp->channels()) {
      case 0:
	flush_delay(chan,delay);
	delay=dur-needdelay;
	if (delay<0)
	  throw Failed("Insufficient delay period after +/| synchronised auto pulse");
	needdelay=0.0;
	break;
      case 1: {
	if (evidp->type==EventID::automatic) {
	  const double cor=evidp->inset();
	  if (delay<cor)
	    throw Failed("auto event overlaps with previous event / start of sequence");
	  needdelay=dur-cor;
	  delay-=cor;
	}
	flush_delay(chan,delay);
	delay=0.0;	
	const eventlist_t& curlist(evidp->events.front());
	const char sync=getsync(*evidp);
	for (size_t item=0;item<curlist.size();item++) {
	  RFEvent* p(curlist(item).get());
	  foundevents=true;
	  const TimedEvent& tev(curseq.push_back(p,timeschan,sync));
	  timeschan+=tev.duration();
	}
      }
	break;
      default:
	throw InternalError("AsyncSequence::raw_rebuild");
      }
    }
    if (needdelay)
      throw Failed("+/| synchronised auto pulse at end of sequence fragment");
    flush_delay(chan,delay);
  }
  if (!foundevents && (delay==0.0) && !hasduration)
    nullsequencelist_warning.raise(name(),true);
  return 0.0; //can't have delay period with async
}

void CompSequenceBase::flagdirty(dirty_t statusv) 
{ 
  if (statusv>status_) {
    status_=statusv;
    if ((verbose & VER_GEN) && (verbose_level>1)) 
      std::cout << "Increased 'dirty' status of " << name() << " to: " << statusv << '\n';
  }    
}

void parse_tolerance()
{
  if (!parser_isnormal()) {
    std::cout << "tolerances:  " << (tolerance*1e6) << " us (timing)  " << (phasetolerance) << " degrees (phase)\n";
    return;
  }

  const double val=parse_double();
  if (val<=0.0)
    error_abort("tolerance must be >0");

  static flagsmap_type flags;
  if (flags.empty()) {
    flags["timing"]=0;
    flags["phase"]=1;
  }
  const int flagsval=parse_flags(flags,0,true);
  if (flagsval)
    phasetolerance=val;
  else
    tolerance=val*1e-6;
}

void parse_cache_limit()
{
  static const size_t M_to_bytes=1024*1024;
  if (are_left()) {
    const double value=parse_double();
    if (value>NMRSIM_HUGE_CACHE_LIMIT)
      error_abort("Unfeasibly large cache_limit - value is specified in M");

    const long bytes=long(value*M_to_bytes);
    if (bytes<0)
      error_abort("cache_limit must be >=0");
    global_cache.limit(bytes);    
  }
  else
    std::cout << "cache_limit=" << global_cache.limit()/double(M_to_bytes) << " M (Used " << ((100.0*global_cache())/global_cache.limit()) << "%)\n";
}

void parse_channel()
{
  if (CurSeqp && (curchannel==0))
    error_abort("can't specify channel during synchronous sequence");

  const size_t chan=parse_int();
  if (chan<1 || chan>nchannels)
    error_abort("channel outside range");

  curchannel=chan;
  ensure_sequence();
}

template<class OpGen> PulseGeneratorBase* create_pgen(const OpGen& opgen, size_t j)
{
  const size_t nucid(nucids(j));
  static const int loc_verbose=(verbose & VER_GEN) ? verbose_level : 0;
  PulseGeneratorBase* pgenp= new PulseGenerator(opgen,nucid,loc_verbose,&global_cache);
  if (!transient_amps.empty())
    pgenp->transient_amplitude(transient_amps(j)->get());
  return pgenp;
//   switch (transient_model) {
//   case TRANS_NONE: case TRANS_MANUAL:
//     return new PulseGenerator(opgen,nucid,loc_verbose);
//   case TRANS_SIMPLE: {
//     PulseGenerator_SimpleTransients* pgenp= new PulseGenerator_SimpleTransients(opgen,nucid,loc_verbose);
//     pgenp->transient_amplitude(transient_amps(j));
//     return pgenp;
//   }
//   default:
//     throw Failed("Unknown transient model");
//   }
}

template<class OpGen> void ensure_pgens(const OpGen& opgen)
{
  static bool reentrylock=false;
  if (reentrylock)
    throw InternalError("re-entering ensure_pgens");
  reentrylock=true;
  pgenstack.create(nchannels);
  ListList<double> globalFz;
  for (size_t j=nchannels;j--;) {
    pgenstack(j).reset(create_pgen(opgen,j));
    globalFz+=opgen.diag_Fz(nucids(j));    
  }
  const int verb=(verbose & VER_GEN) && (verbose_level>1) ? 1 : 0;
  zshift_cachep.reset(new ZshiftCache(globalFz,verb));
}

void MasterObj::ensure_channels()
{
  if (!partitions.empty())
    return;
  partitions.create(1<<nchannels);

  const int defflags=( optcombinepropagators() ? BaseMetaPropagator::combinepropagators : 0) | (!optsync ? BaseMetaPropagator::nosynchronisation : 0);
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Default propagator flags: " << defflags << '\n';

  prop_flags.create(1<<nchannels,defflags);
  if (nchannels==0)
    return;
#ifndef NOPERIODIC
  if (!simple_opgenp)
    ensure_pgens(*crystal_opgenp);
  else
#endif
    ensure_pgens(*simple_opgenp);
}

void reset_pgens()
{
  for (size_t i=pgenstack.size();i--;)
    pgenstack(i)->reset();
}

void store_transients(LIST<transient_state>& states, double t)
{
  if (transient_model==TRANS_NONE)
    return;
  states.create(nchannels);
  if (verbose & VER_GEN) {
    std::cout << "Storing transient state at t=";
    prettyprint_time(t) << '\n';
  }
  for (size_t i=nchannels;i--;)
    pgenstack(i)->store_transient(states(i));
}

void restore_transients(const LIST<transient_state>&)
{
//   if (!hastransients())
//     return;
//   if (states.empty())
//     throw Failed("restore_transients: called before sequence transients have been determined");
//   if (verbose & VER_GEN)
//     std::cout << "Restoring transient state\n";
//   for (size_t i=nchannels;i--;)
//     pgenstack(i)->restore_transient(states(i));
}

//flush any remaining transients
void flush_transients(BlockedMatrix<complex>&)
{
//   if (hastransients()) {
//     if (verbose & VER_GEN)
//       std::cout << "Flushing transients\n";
//     for (size_t i=nchannels;i--;)
//       pgenstack(i)->clear_transient(U);
//   }
}

void flush_transients()
{
  for (size_t i=nchannels;i--;)
    pgenstack(i)->clear_transient();
}

    
bool transients_active()
{
//   if (hastransients()) {
//     for (size_t i=nchannels;i--;) {
//       if (!(pgenstack(i)->isfinished()))
// 	return true;
//     }
//   }
  return false;
}

void reset_seqlists()
{
  for (size_t i=cycledseqlist.size();i--;)
    cycledseqlist(i)->reset();
}

template<class T> void MasterObj::set_Hsystem(const BlockedMatrix<T>& H)
{
  //! slightly wasteful as each pgen has its own copy...
  reset_pgens(); //reset pulse generators (otherwise Hsystem may fail)
  for (size_t j=nchannels;j--;)
      pgenstack(j)->Hsystem(H);      
}

template void MasterObj::set_Hsystem(const BlockedMatrix<complex>&);
template void MasterObj::set_Hsystem(const BlockedMatrix<double>&);

PulseGeneratorBase& pulse_generator(size_t chan)
{
  if (chan>=nchannels)
    throw BadIndex("pulse_generator: bad channel index",chan,nchannels);
  if (pgenstack.empty())
    throw Failed("pulse_generator: called before pulse generators created");
  return *(pgenstack(chan));
}

usage_t CycledSequence::usage() const
{
  usage_t tot;
  for (size_t i=states_.size();i--;)
    tot+=states_(i).usage();
  return tot;
}
