/*! \file 
  \brief  Pulse sequence "actions" */

#include "Action.h"
#include "NMRsim_RF.h"
#include "expression_definition.hpp"
#include <set>
#include "rational.h"
#include "NMRsim_logfile.h"
#include "InversionGenerator.h"
#include <sstream>

//static LIST<size_t> control_stack;
bool insidecontrol=false; //!< true if inside control structure

using namespace libcmatrix;

const double rational::tolerance=1e-8;  //!< tolerance within which two ratios are regarded as equal

// struct inversionhandle {
//   explicit inversionhandle(nuclei_spec nucv =H_NUCLEUS)
//     : nuc_(nucv), invgenp(NMRSIM_NULL) {}
//   const InversionGenerator& operator()() const;
//   nuclei_spec nuc_;
//   mutable const InversionGenerator* invgenp;
// };

const InversionGenerator* create_inversiongenerator(size_t nuc)
{
  if (!masterobjp || !cstructp)
    throw InternalError("create_inversiongenerator");
#ifndef NOPERIODIC
  if (!(masterobjp->simple_opgenp))
    return new InversionGenerator(*(masterobjp->crystal_opgenp),*(masterobjp->Hstructp), *cstructp, nuc);
#endif
  return new InversionGenerator(*(masterobjp->simple_opgenp),*(masterobjp->Hstructp), *cstructp, nuc);
}

// const InversionGenerator& inversionhandle::operator()() const
// {
//   if (!invgenp) {
//     if (!masterobjp || !cstructp)
//       throw InternalError("inversionhandle");
// #ifndef NOPERIODIC
//     if (!(masterobjp->simple_opgenp))
//       invgenp=new InversionGenerator(*(masterobjp->crystal_opgenp),masterobjp->Hstruct, *cstructp, nuc_);
//     else
// #endif
//       invgenp=new InversionGenerator(*(masterobjp->simple_opgenp),masterobjp->Hstruct, *cstructp, nuc_);
//   }
//   return *invgenp;
// }

typedef std::map<size_t,const InversionGenerator*> inversionmap_type;
static inversionmap_type inversionmap;

//void parse_makeinversion();

bool istimedependent() { return (spin_rate!=0.0); }

static size_t local_acq_count=0;
static double dt=0.0;
double dwelltime() { return dt; }

Action_Factory_t& get_Action_Factory()
{
  static Action_Factory_t Action_Factory;
  return Action_Factory;
}

// Action_Factory_t& get_initialise_Factory()
// {
//   static Action_Factory_t Action_Factory;
//   return Action_Factory;
// }

LIST<actionstack_t> actionstacks;
//actionstack_t initialisestack;

//bool haveduration=false;
bool haveactions=false;

static const double deg_to_rad=M_PI/180.0;
static const double rad_to_deg=180.0/M_PI;

EXPR_DECLARE_FUNCTION(FunctionSyncRatio,"sync_ratio")

const BlockedOperator& validate_operator(const BlockedOperator& a)
{
  if (!a)
    throw InternalError("Operator matrix undefined");
  return a;
}

struct Action_Proxy_ {
  Action_Proxy_() {

    command_Factory_t& par_Factory(get_par_Factory());
    par_Factory["filter"]=par_t(&parse_makefilter,true);
    //    par_Factory["inversion"]=par_t(&parse_makeinversion,true);

    Action_Factory_t& Action_Factory(get_Action_Factory());
    Action_Factory["acq"]=&ActionAcq::create;
    Action_Factory["acqn"]=&ActionAcqN::create;
    Action_Factory["acqpoint"]=&ActionAcqPoint::create;
    Action_Factory["do"]=&ActionDo::create;
    Action_Factory["pulse180"]=&ActionInversion::create;
    Action_Factory["prop"]=&ActionProp::create;
    Action_Factory["propfor"]=&ActionProp::createfor;
    Action_Factory["rotor_angle"]=&ActionRotorAngle::create;
    Action_Factory["filter"]=&ActionFilter::create;
    Action_Factory["putmatrix"]=&ActionPutMatrix::create;
    Action_Factory["echo"]=&ActionEcho::create;
    Action_Factory["get"]=&ActionGet::create;
    //    Action_Factory["set"]=&ActionSet::create;
    Action_Factory["scale"]=&ActionScale::create;
    Action_Factory["timeadjust"]=&ActionTimeAdjust::create;
    Action_Factory["transfer"]=&ActionTransfer::create;
    Action_Factory["exchange"]=&ActionExchange::create;

    const static FunctionSyncRatio syncratio2obj(2U);
    const static FunctionSyncRatio syncratio3obj(3U);
    declare_function("sync_ratio",2U,syncratio2obj);
    declare_function("sync_ratio",3U,syncratio3obj);
  }
};

//declare ActionCommands
static Action_Proxy_ action_proxy_;

void ActionProp::set_ntimes(int ntimesv)
{
  if (ntimes<0)
    throw InternalError("ActionProp::set_ntimes");
  ntimes=ntimesv;
  if (ntimes<0)
    throw InvalidParameter("ActionProp::set_ntimes");
  reset();
}

void ActionProp::set_propfor(double propforv)
{
  if (ntimes>=0)
    throw InternalError("set_propfor");
  propfor=propforv;
  if (propfor<0.0)
    throw InvalidParameter("ActionProp::set_propfor");
  reset();
}

ActionAcqN::ActionAcqN(size_t npts, double phasev, CompSequenceBase* seqp, int flagsv)
  : toacquire_(npts),
    acqp(new ActionAcq(seqp,flagsv)),
    acqpointp(new ActionAcqPoint(phasev)) {}

usage_t ActionAcq::usage() const
{
  const usage_t use= Us.empty() ? usage_t() : usage_t(Us.size()*Us.front().size(),Type2Type<complex>());
  return use+usage_t(Uacc);
}

// ActionAcqBase::ActionAcqBase(CompSequenceBase* seqp_, bool beverbosev)
//   : ActionCommand(beverbosev), seqp(seqp_) {}

ThreadWarning<> ActionAcq::zerodurationacquisition_warning("sw set with zero-duration acquisition sequence",&NMRsim_repeat_warning);

ActionAcq::ActionAcq(CompSequenceBase* seqp_, int flagsv)
  : seqp(seqp_), nrep(0), flags_(flagsv)
{
  ndim=local_acq_count;
  const double sw=spectralwidth();

  if (seqp && (seqp->duration()==0.0) && sw)
    zerodurationacquisition_warning.raise();

  // dt = sw ? 1.0/sw : 0.0;

  bool checksw=true;
  if (seqp) {
    haveactions=true; //flag that we have pulses etc.
//     if (seqp->duration())
//       haveduration=true;
//     else {
    if (seqp->duration()==0.0) {
      process2D=false; //can't do processing
      checksw=false;
    }
  }
  if (checksw && (sw==0.0)) {
    parser_printcontext() << "Spectral width for dimension " << (ndim+1) << " has not been set\n";
    error_abort();
  }
  reset();
}

size_t ActionAcq::dimensionsize() const { return (ndim==ndims-1) ? np : ns(ndim); }
double ActionAcq::spectralwidth() const { return (ndim==ndims-1) ? sw : sws(ndim); }
double ActionAcq::dwelltime() const
{
  const double sw=spectralwidth();
  return sw ? 1.0/sw : 0.0;
}

ActionProp::ActionProp(CycledSequence& seq_, int ntimesv, int flags)
  : seqp(&seq_), propfor(0.0), ntimes(0), flags_(flags), checkedsync(false)
    // invert(flags & PROP_INV)
{
  haveactions=true; //flag that we have pulses etc.
//   if (seqp->duration())
//     haveduration=true;
  set_ntimes(ntimesv);
}

ActionProp::ActionProp(CycledSequence& seq_, double propforv, int flags)
  : seqp(&seq_), ntimes(-1), flags_(flags), checkedsync(false)
    // invert(flags & PROP_INV)
{
  haveactions=true; //flag that we have pulses etc.
//   if (seqp->duration())
//     haveduration=true;
  set_propfor(propforv);
}

ActionProp::ActionProp(double propforv, int flags)
  : seqp(NMRSIM_NULL), ntimes(-1), flags_(flags)
    //invert(flags & PROP_INV)
{
  set_propfor(propforv);
}

ActionCommand* ActionEcho::create()
{
  ActionCommand* newac = new ActionEcho(get_curline());
  set_curline(NMRSIM_NULL);
  return newac;
}

ActionCommand* ActionInversion::create()
{
  size_t nuc=NULL_NUCLEUS;
  if (count_left()>1)
    nuc=parse_nucleusname(parse_string(F_REPLACEDOLLAR));
  else {
    if (sysp->ishomonuclear())
      nuc=(*sysp)(0).nucleus();
    else
      error_abort("can only omit nucleus type if system is homonuclear");
  }
  const InversionGenerator* invgenp=NMRSIM_NULL;
  const inversionmap_type::const_iterator iter(inversionmap.find(nuc));
  if (iter==inversionmap.end()) {
    invgenp=create_inversiongenerator(nuc);
    inversionmap[nuc]=invgenp;
  }
  else
    invgenp=iter->second;
  Variable cvar(S_PHASE);
  const double phase=parse_double(&cvar,F_ISPHASE);
  return new ActionInversion(nuc,*invgenp,phase);
}

void ActionInversion::print(std::ostream& ostr) const
{
  ostr << "pulse180 " << nuctolabel(nuc_) << ' ' << phase_;
}

void ActionInversion::exec(MasterObj*, BlockedOperator& sigma, double&) const
{
  invgen_(sigmatmp_,sigma,phase_*deg_to_rad);
  sigma.swap(sigmatmp_);
}

size_t ActionDirectAcq::acquired_=0U;

ActionDirectAcq::ActionDirectAcq(double dphsv, CompSequenceBase* seqp_, double synctime_)
  : dphs(dphsv), seqp(seqp_), synctime(synctime_), checkedsync(false)
				//, toacquire_(-1)
{}
//   if (np<1)
//     error_abort("np not set (or ridiculous!)");
// }

// ActionDirectAcq::ActionDirectAcq(size_t toacquirev, CompSequenceBase* seqp_, double synctime_)
//   : seqp(seqp_), synctime(synctime_), checkedsync(false), toacquire_(toacquirev)
// {}

void ActionDirectAcq::add_acq(bool israw)
{
  if (create_status==FINAL) 
    error_abort("can't have more than one direct dimension acq");
  if (create_status==NONE)
    local_acq_count++;
  create_status=israw ? PARTIAL : FINAL;
}

ActionDirectAcq* ActionDirectAcq::create(double phasev, CompSequenceBase* acqseqp, double synctimev)
{
//   bool israw=(np>=0);  
  add_acq(false);
  if (acqseqp)
    haveactions=true;
  return new ActionDirectAcq(phasev,acqseqp,synctimev);
//   bool israw=(np>=0);  
//   add_acq(israw);
//   ActionDirectAcq* acqp=israw ? new ActionDirectAcq(np,acqseqp,synctimev) : new ActionDirectAcq(acqseqp,synctimev);
//   if (acqseqp)
//     haveactions=true;
//   return acqp;
}

// bool hasacquisition_rf()
// {
//   return acqp ? (acqp->sequence()!=NMRSIM_NULL) : false;
// }

ActionCommand* ActionAcqPoint::create()
{
  if (!ActionAcq::isdirectdimension())
    error_abort("acqpoint can only be used in direct dimension");
  if (isfrequency())
    error_abort("acqpoint can not be used in frequency (histogram) mode");
  Variable cvar(S_PHASE);
  const double phasesh=parse_double(&cvar,F_ISPHASE);
  ActionCommand* comp=new ActionAcqPoint(phasesh);
  ActionDirectAcq::add_acq(true); //!< do afterwards (updates acq_count)
  return comp;
}

bool ActionAcq::isdirectdimension()
{
  bool isfinal=(local_acq_count==nacqdims-1);
  if ((local_acq_count==nacqdims) && (ActionDirectAcq::create_status==ActionDirectAcq::PARTIAL))
    isfinal=true;
  return isfinal;
}

ActionCommand* ActionAcqN::create()
{
  ActionCommand* comp=ActionAcq::create_(true);
  ActionDirectAcq::add_acq(true); //!< do afterwards (updates acq_count)
  return comp;
}

const flagsmap_type& ActionAcq::get_flags()
{
  static flagsmap_type flags;
  if (flags.empty())
    flags["putmatrix"]=1;
  return flags;
}

#define BASEACTIONSYNTAX "[<acquistion sequence>|-]#[<sync time / us>]"
const char* actionfinalsyntax="<N pts>#<acq phase>#" BASEACTIONSYNTAX;
const char* actionsyntax="<N pts>#" BASEACTIONSYNTAX;

ActionCommand* ActionAcq::create_(bool israw)
{
  const bool isfinal=isdirectdimension();
  if (!isfinal) {
    if (local_acq_count>=nacqdims)
      error_abort("can't specify more acq's than dimensions!");
    if (israw)
      error_abort("can't use acqn for indirect dimensions");
  }
  const char* syntax = isfinal ? actionfinalsyntax : actionsyntax;
  size_t offset = isfinal ? 1 : 0;
  if (!israw)
    offset++;

  int nacq=-1; 
  double phasesh=0.0;
  if (isfinal) {
    if (israw) {
      nacq=parse_unsigned_syntax(actionfinalsyntax,1);   
      if (nacq==0)
	ActionDirectAcq::empty_warning.raise();
    }
    if (are_left()) {
      Variable cvar(S_PHASE);
      phasesh=parse_double_syntax(actionfinalsyntax,offset,&cvar,F_ISPHASE);
    }
    //      parse_system_variable(v_dphs,F_ISPHASE);
  }

  CompSequenceBase* acqseqp=NMRSIM_NULL;
  if (parser_isnormal()) {
    char* seqname(parse_string_syntax(syntax,offset+1,F_REPLACEDOLLAR));
    if (strcmp(seqname,"-")!=0)
      acqseqp=findmap(seqmap,seqname,"sequence");
  }

  if (isfinal && !israw) {
    double synctime=0.0;
    if (are_left()) {
      Variable cvar(S_ARG1);
      synctime=1e-6*parse_double_syntax(syntax,offset+2,&cvar);
    }
    return ActionDirectAcq::create(phasesh,acqseqp,synctime);
  }
      
//     error_abort("Additional argument to acq: note that synchronisation times are now specified when sequence is store'd");

  const int lflags=parse_flags(get_flags());

  if (israw)
    return new ActionAcqN(nacq,phasesh,acqseqp,lflags);

  ActionAcq* acqp=new ActionAcq(acqseqp,lflags);
  local_acq_count++;
  return acqp;
}

static BlockedOperator* verifyspec(const char* name, const matrixspec& op, bool allowH =false)
{
  switch (op.type) {
  case matrixspec::OPERATOR:
    return &(op.asoperator().op);
    break;
  case matrixspec::SPECIAL:
    if (allowH)
      return NMRSIM_NULL;
  default:
    break;
  }
  parser_printcontext() << name << " is not a valid operator matrix\n";
  error_abort();
  return NMRSIM_NULL; //just to avoid warning
}

// void ActionSet::print(std::ostream& ostr, subsid_t subsid) const
// {
//   ostr << "set " << name_;
//   ostr << ' ' << scale_;
//   if (top_)
//     ostr << ' ' << toname_;
//   if (subsid==S_NONE)
//     ostr << '\n';
// }

// void ActionSet::set(double val,subsid_t)
// {
//   scale_=val;
// }

// void ActionSet::exec(MasterObj*, BlockedOperator& sigma, double&) const
// {
//   BlockedOperator& dest(top_ ? *top_ : sigma);
//   dest=op_;
//   dest*=scale_;
// }

ActionTransfer::ActionTransfer(const named_operator& fromv, const named_operator& tov)
  : from_(fromv), to_(tov)
{}

void ActionTransfer::print(std::ostream& ostr) const
{
  ostr << "transfer " << from_.name << ' ' << to_.name;
}

void ActionTransfer::exec(MasterObj*, BlockedOperator& sigma, double&) const
{
  const complex scale=NMR_trace(validate_operator(from_.op),sigma);
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Transfering with scale factor: " << scale << '\n';
  sigma=validate_operator(to_.op);
  sigma.row()*=scale;
}

ActionExchange::ActionExchange(const BaseList<named_operator>& opsv, const ExchangeMatrix& matv)
  : ops_(opsv),
    mat_(matv) {}

void ActionExchange::print(std::ostream& ostr) const
{
  ostr << "exchange ";
  for (size_t i=0;i<ops_.size();i++)
    ostr << ops_(i).name << ' ';
  mat_.printraw(ostr);
}

void ActionExchange::exec(MasterObj*, BlockedOperator& sigma, double& t) const
{
  const size_t nsites=mat_.sites();
  if (ops_.size()!=nsites) {
    parser_printthread(std::cerr) << "exchange: mismatch between number of operators (" << ops_.size() << ") and size of exchange matrix (" << nsites << ")\n";    
    error_abort();
  }
  amps0_.create(nsites);
  for (size_t i=nsites;i--;) {
    const BlockedOperator& curop(validate_operator(ops_(i).op));
    const complex trace=NMR_trace(curop,sigma);
    const double opnorm=norm(curop);
    amps0_(i)=trace/opnorm;
  }
  multiply(amps_,mat_(),amps0_);
  if (verbose & VER_GEN) {
    if (verbose_level>1)
      std::cout << "Applied exchange matrix\n" << mat_;
    std::cout << "Initial coherence amplitudes: " << amps0_ << '\n';
    std::cout << "Final coherence amplitudes: " << amps_ << '\n';
  }
  for (size_t i=nsites;i--;)
    mla(sigma,amps_(i)-amps0_(i),ops_(i).op);
}

ThreadWarning<> ActionGet::ignoring_imaginary_warning("get: ignoring imaginary component",&NMRsim_once_warning);

void ActionGet::exec(MasterObj*, BlockedOperator& sigma, double& t) const
{
  const BlockedOperator& which(isdensity() ? sigma : from_.asoperator().op);
  const complex val=NMR_trace(validate_operator(op_),which);
  if (imag(val)>tolerance) {
    std::ostringstream str(std::ostringstream::out);
    str << " of " << val << " at time t=";
    prettyprint_time(t,str);
    ignoring_imaginary_warning.raise(str.str().c_str());
  }
  dest_.set(real(val));
}

ContextWarning<> ActionExchange::exchangematrix_warning("exchange matrix ought not to be used directly",&NMRsim_once_warning);

ActionCommand* ActionExchange::create()
{
  const size_t n=count_left();
  if (n<2)
    error_abort("syntax: exchange <op1>+ <exchange matrix name>");
  const size_t nsites=n-1;
  LIST<named_operator> ops;
  ops.reserve(nsites);
  for (size_t i=0;i<nsites;i++)
    ops.push_back(named_operator(parse_string(F_REPLACEDOLLAR)));
  const char* name=parse_string(F_REPLACEDOLLAR);
  const matrixspec& matref=findmap(matrixmap,name,"matrix");
  if (matref.type!=matrixspec::EXCHANGE) {
    parser_printcontext() << name << " is not an exchange matrix (created by matrix set exchange in par)";
    error_abort();
  }
  const ExchangeMatrix& exchref(matref.asexchange());
  if (exchref.isexchange())
    exchangematrix_warning.raise();
  return new ActionExchange(ops,exchref);
}
  
static BlockedOperator& getop(char* name)
{
  const matrixmap_type::iterator curp=matrixmap.find(name);
  if (curp==matrixmap.end()) {
    productoperator_spec* newopp=create_productoperator(name);
    BlockedOperator* opexprp=new BlockedOperator();
    *opexprp=masterobjp->make(*newopp);
    delete newopp;
    return *opexprp;
  }
  return *verifyspec(name,curp->second);
}
 
named_operator::named_operator(char* namev)
: name(namev),
op(getop(namev)) {}

ActionCommand* ActionTransfer::create()
{
  char* fromname=parse_string(F_REPLACEDOLLAR);
  const BlockedOperator& fromop=getop(fromname);
  char* toname=parse_string(F_REPLACEDOLLAR);
  const BlockedOperator& toop=getop(toname);  
  return new ActionTransfer(named_operator(fromname,fromop),named_operator(toname,toop));
} 
 
ActionGet::ActionGet(UserVariable& destv, const char* namev, const BlockedOperator& opv, const char* fromnamev, matrixspec fromv)
  : dest_(destv), name_(namev), op_(opv), fromname_(fromnamev), from_(fromv)
{}

ActionGet::ActionGet(UserVariable& destv, const char* namev, const BlockedOperator& opv)
  : dest_(destv), name_(namev), op_(opv)
{}

void ActionGet::print(std::ostream& ostr) const
{
  ostr << "get " << dest_.name() << ' ' << name_;
  if (!isdensity())
    ostr << ' ' << fromname_;
}

// ActionSet::ActionSet(double scalev, const char* namev, const BlockedOperator& opv)
//   : scale_(scalev), name_(namev), op_(opv), top_(NMRSIM_NULL) {}

// ActionSet::ActionSet(double scalev, const char* namev, const BlockedOperator& opv, const char* tonamev, BlockedOperator& tov)
//   : scale_(scalev), name_(namev), op_(opv),
//     toname_(tonamev), top_(&tov) {}

// ActionCommand* ActionSet::create()
// {
//   Variable cvar(S_ARG1);
//   const double scale=parse_double(&cvar);
//   if ((cvar.subsid==S_NONE) && allowwarnings())
//     parser_printcontext() << "warning: set with constant scale factor doesn't make much sense!";
//   const char* name=parse_string();
//   const BlockedOperator& op=getop(name);
//   if (are_left()) {
//     const char* toname=parse_string();
//     BlockedOperator& to=getop(toname);
//     return new ActionSet(scale,name,op,toname,to);
//   }
//   return new ActionSet(scale,name,op);
// }

ContextWarning<> ActionGet::notimpl_warning("get from non-operator matrix not implemented (ignored)",&NMRsim_repeat_warning);

ActionCommand* ActionGet::create()
{
  UserVariable& var=*findnewvariable(parse_variable_name());
  char* name=parse_string(F_REPLACEDOLLAR);
  const BlockedOperator& op=getop(name);    
  if (are_left()) {
    const char* fromname=parse_string(F_REPLACEDOLLAR);
    matrixspec from=findmap(matrixmap,name,"matrix");
    if (from.type!=matrixspec::OPERATOR) {
      notimpl_warning.raise();
      return NMRSIM_NULL;
    }
    return new ActionGet(var,name,op,fromname,from);
  }
  return new ActionGet(var,name,op);
}

void ActionAcqN::print(std::ostream& ostr) const
{
  ostr << "acqn " << toacquire_;
  if (acqp->sequence())
    ostr << ' ' << acqp->sequence()->name();
}

void ActionAcq::print(std::ostream& ostr) const
{
  ostr << "acq (dimension " << (ndim+1) << ')';
  if (nrep)
    ostr << " (x" << nrep << ") ";
  else
    ostr << " (aperiodic) ";
  if (seqp)
    ostr << seqp->name();
  if (flags_) {
    ostr <<  ' ';
    dumpflags(ostr,get_flags(),flags_);
  }
}

void ActionDirectAcq::print(std::ostream& ostr) const
{
  ostr << "acq";
//   if (toacquire_>=0)
//     ostr << "n " << toacquire_ << ' ';
//   else
    ostr << " (direct) ";
  ostr << dphs;
  if (seqp)
    ostr << ' ' << seqp->name();
  if (synctime)
    ostr << ' ' << (synctime*1e6);
}

void ActionAcqPoint::print(std::ostream& ostr) const
{
  ostr << "acqpoint " << dphs;
}


void ActionDirectAcq::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << ((subsid==S_PHASE) ? "acq_phase" : "acq_synctime");
}

void ActionAcqN::printvariablename(std::ostream& ostr, subsid_t) const
{
  ostr << "acq_phase";
}

void ActionAcqN::set(double v, subsid_t subsid)
{
  acqpointp->set(v,subsid);
}

void ActionAcqPoint::printvariablename(std::ostream& ostr, subsid_t) const
{
  ostr << "acq_phase";
}

void ActionProp::print(std::ostream& ostr) const
{
  ostr << ((ntimes<0) ? "propfor " : "prop ");
  if (ntimes<0) {
    ostr << propfor << ' ';
    if (seqp)
      seqp->printshort(ostr);
  }
  else {
    if (seqp==NMRSIM_NULL)
      throw InternalError("ActionProp::print");
    seqp->printshort(ostr);
    //  ostr << seqp->name;
    if (ntimes>1)
      ostr << ' ' << ntimes;
  }
  if (flags_) {
    ostr << ' ';
    dumpflags(ostr,get_flags(),flags_);
  }
//   if (invert)
//     ostr << " -inv";
}
  
void ActionProp::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << ((ntimes<0) ? "propfor_" : "prop_");
  if (seqp) {
    seqp->printshort(ostr);
    ostr << '_';
  }
  if (subsid==S_ARG1) 
    ostr << ((ntimes<0) ? "duration" : "repeat_count");
  else
    ostr << "element_" << (1+CycledSequence::subsid_to_index(subsid)) << "_phase";
}
  
ActionFilter::ActionFilter(const char* name_,const filter_def& filter_)
  : name(name_), filter(filter_) {}
//  if (!filter_) //!< undefined filter should be caught later
//    throw Undefined("ActionFilter");
//}

void ActionFilter::print(std::ostream& ostr) const
{
  ostr << "filter " << name << ' ' << filter.spec << '\n';
  if ((verbose & VER_GEN) && !!filter)
    ostr << filter();
}

ContextWarning<> ActionProp::propafteracq_warning("prop after final acq (ignored)",&NMRsim_once_warning);
ContextWarning<> ActionProp::cyclingmismatch_warning("propagator cycling doesn't match ni step",&NMRsim_repeat_warning);
ContextWarning<> ActionProp::nocache_warning("prop used without sequence does not cache propagators - consider introducing named delay periods to indicate linked fragments that would benefit from caching",&NMRsim_once_warning);
ContextWarning<> ActionProp::inv_warning("-inv is not valid in many cases - consider using pulse180 or alternatives",&NMRsim_once_warning);
ContextWarning<> ActionProp::invalidreset_warning("-reset is not applicable without sequence (ignored)",&NMRsim_repeat_warning);

const flagsmap_type& ActionProp::get_flags()
{
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["inv"]=PROP_INV;
    flags["putmatrix"]=PROP_PUTMATRIX;
    flags["reset"]=PROP_RESET;
    //    flags["full"]=PROP_FULL;
  }
  return flags;
}

ActionCommand* ActionProp::create_(bool ispropfor)
{
  if (local_acq_count>nacqdims) {
    propafteracq_warning.raise();
    return NMRSIM_NULL;
  }

  Variable cvar(S_ARG1);
  double lpropfor=-1.0;
  int lntimes=1;
  if (ispropfor) {
    lntimes=-1;
    lpropfor=parse_double(&cvar);
    if (lpropfor<0.0)
      error_abort("duration cannot be negative");
  }

  CycledSequence* curseqp=NMRSIM_NULL;
  if (!ispropfor || parser_isnormal()) {
    curseqp=parse_cycledsequence();
    const size_t len=curseqp->cyclelength();
    if (len>1) {
//       if (!active2D)
// 	error_abort("can only cycle propagators in nD simulations");
      if (active2D && (len!=skips(local_acq_count)))
	cyclingmismatch_warning.raise();
      if (ispropfor) {
	if (curseqp->hascycle())
	  error_abort("sequence list containing [] cycle cannot be used with propfor");
      }
      else {
	if (curseqp->hasCW())
	  error_abort("CW fragments cannot be sensibly used with prop, use propfor instead");
      }
    }  
    if (!ispropfor && parser_isnormal())
      lntimes=parse_unsigned(&cvar);
  }
  if (!curseqp && !nochecks)
    nocache_warning.raise();

  int lflags=parse_flags(get_flags());
  if ((lflags & PROP_RESET) && (curseqp==NMRSIM_NULL)) {
    invalidreset_warning.raise();
    lflags-=PROP_RESET;
  }

  if (lflags & PROP_INV) {
    if (!(sysp->ishomonuclear()))
      error_abort("-inv cannot be used with heteronuclear spin system");
    inv_warning.raise();
  }

  if (ispropfor) {
    return curseqp
      ? new ActionProp(*curseqp,lpropfor,lflags)
      : new ActionProp(lpropfor,lflags);
  }
  return new ActionProp(*curseqp,lntimes,lflags);
}

void ActionProp::set(double val,subsid_t subsid)
{
  if (subsid==S_ARG1) {
    if (ntimes<0)
      set_propfor(val);
    else
      set_ntimes(round_int(val));
    return;
  }
  if (!seqp)
    throw InternalError("ActionProp::set");
  const size_t ind=CycledSequence::subsid_to_index(subsid);
  seqp->set(val,ind);
}

// void ActionAcqBase::set(double synctime_, subsid_t)
// {
//   synctime=synctime_*1e-6;
//   if (synctime<0.0)
//     error_abort("acq: synchronisation time can't be <0");
// }

// ActionCommand* ActionMaxdt::create()
// {
//   if (auto_vars & A_MAXDT)
//     std::cerr << "Warning: changing maxdt during sequence overrides maxdt auto_opt\n";
//   return new ActionMaxdt(parse_double()*1e-6);
// }

ActionCommand* ActionRotorAngle::create()
{
  Variable cvar(S_ARG1);
  return new ActionRotorAngle(parse_double(&cvar));
}

ContextWarning<> ActionTimeAdjust::static_warning("timeadjust has no effect when Hamiltonian is time independent (zero spin rate) - don't confuse with propagation with a delay period!",&NMRsim_once_warning);

ActionCommand* ActionTimeAdjust::create()
{
  if (!istimedependent())
    static_warning.raise();
  Variable cvar(S_ARG1);
  const double dt=parse_double(&cvar);
  static flagsmap_type flags;
  if (flags.empty())
    flags["-absolute"]=1;
  const bool isabs=parse_flags(flags);
  return new ActionTimeAdjust(dt,isabs);
}

void ActionTimeAdjust::print(std::ostream& ostr) const 
{
  ostr << "timeadjust " << time_;
  if (isabsolute_)
    ostr << " -absolute";  
}

ThreadWarning<> ActionTimeAdjust::back_warning("timeadjust: moving time backwards!",&NMRsim_once_warning);

void ActionTimeAdjust::exec(MasterObj*, BlockedOperator&, double& t) const
{
  const double usetime=time_*1e-6;
  if (isabsolute_) {
    if (usetime<t)
      back_warning.raise();
    t=usetime;
  }
  else {
    if (time_<0.0)
      back_warning.raise();
    t+=usetime;
  }
  if (verbose & VER_GEN) {
    std::cout << "Adjusted time to ";
    prettyprint_time(t,std::cout) << '\n';
  }
}

// void actionstack_t::parse(char* lbuf)
// {
//   ActionCommand* newac=parse_factory(lbuf,get_Action_Factory(),Type2Type<ActionCommand>());
//   if (newac)
//     this->push_back(newac);
//   char* keyname=parse_string();
//   Action_Factory_t& Action_Factory(get_Action_Factory());
//   const Action_Factory_t::const_iterator curcom=Action_Factory.find(keyname);
//   if (curcom==Action_Factory.end()) {
//     parser_printcontext() << "directive not recognised " << keyname << std::endl;
//     error_abort();
//   }
//   ActionCommand* newac;
//   Mark markobj;
//   try {
//     newac=(curcom->second)();
//     parser_checkfinished();
//   } catch (MatrixException& exc) {
//     parser_printcontext() << exc << '\n';
//   }
//   if (newac) {
//     this->push_back(newac);
//     markobj.flush(newac);
//   }
//   else 
//     error_abort();
//}

void actionstack_t::reset()
{
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Resetting action stack\n";
  const actionstack_t::const_iterator finish(end());
  actionstack_t::const_iterator start(begin());
  while (start!=finish)
    (*start++)->reset();
}

void actionstack_t::flush_cache(dirty_t level)
{
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Flushing action stack: " << level << '\n';
  const actionstack_t::const_iterator finish(end());
  actionstack_t::const_iterator start(begin());
  while (start!=finish) {    
    (*start)->flush_cache(level);
    ++start;
  }
}

std::ostream& operator<< (std::ostream& ostr, const accumulating_timer<>& stopwatch)
{
  ostr << "Called: " << stopwatch.entries() << " times";
  if (stopwatch.entries()) {
    ostr << "  Total: " << stopwatch() << " s  Mean: ";
    const double res=stopwatch.resolution();
    if (res && (stopwatch()<stopwatch.resolution()*stopwatch.entries()))
      ostr << "(less than timing resolution)";
    else
      ostr << (1e6*stopwatch()/stopwatch.entries()) << " us";
  }
  return ostr;
}

void get_variables(BaseList<double> pars)
{
  if (pars.size()<varpars.size())
    throw Mismatch("get_variables");
  const varpars_t::iterator end(varpars.end());
  varpars_t::iterator start(varpars.begin());
  size_t ind=0;
  while (start!=end) {
    pars(ind++)=start->get();
    ++start;
  }
}

 ThreadWarning<> sync_warning("synchronisation not found",&NMRsim_repeat_warning);

namespace {

  double checksynctol(double t)
  {
    if ((t<=0.0) || (t>=1.0))
      error_abort("invalid tolerance (must be >0 and <1)");
    return t;
  }
  
  rational calcsync(double ratio, int cycles, double tol)
  {
    if (ratio>1.0) {
      const rational best=calcsync(1.0/ratio,cycles,1e6).reciprocal();
      return (fabs(best-ratio)<tol) ? best : rational(ratio);
    }
        
    if (cycles<0)
      throw InvalidParameter("sync");
    if (cycles==0)
      return rational(ratio);

    if (ratio>1.02*cycles) {
      char buf[256];
      snprintf(buf,sizeof(buf),": desired ratio (%g:1) is outside range that can be accommodated (%i:1).  Increase number of synchronisation cycles?",ratio,cycles);
      sync_warning.raise(buf);
    }
      
    //static std::set<rational> ratios;
    //    static size_t ratios_cycles=0;
    static std::map<size_t, std::set<rational> > ratiosmap;

    std::set<rational>& ratios(ratiosmap[cycles]);
    if (ratios.empty()) {
      ratios.insert(rational(1,1));
      for (size_t N=2;N<=cycles;N++)
	for (size_t M=N-1;M>0;M--)
	  ratios.insert(rational(M,N));
     
      if (verbose & VER_GEN) {
	std::cout << ratios.size() << " synchronisation ratios found over " << cycles << " cycles: ";
	if ((verbose_level>1) || (ratios.size()<=50))
	  dump(ratios," ");
	std::cout << '\n';
      }

      if (ratios.empty())
	error_abort("No synchronisation factors found: bad argument to calcsync?\n");

      //      ratios_cycles=cycles;
    }

    typedef std::set<rational>::const_iterator iterator;

    const iterator after=std::lower_bound(ratios.begin(),ratios.end(),ratio);
    if (after==ratios.end())
      throw InternalError("calcsync");
    const rational* bestp=NMRSIM_NULL;
    if (after==ratios.begin())
      bestp=&(*after);
    else {
      iterator before=after;
      --before;      
      bestp=(fabs(ratio-*before)<fabs(ratio-*after)) ? &(*before) : &(*after);
    }
    return (fabs(*bestp-ratio)<tol) ? *bestp : rational(ratio);
  }
}

// void dosyncratio(LIST<double>& dest, const BaseList< BaseList<double> >& args, double tol)
// {
//   assert(args.front().size()==1);
//   assert(args(1U).size()==1);
// }

void FunctionSyncRatio::operator()(LIST<double>& dest, const BaseList<double>& args) const
//void FunctionSyncRatio2::operator()(LIST<double>& dest, const Expression& expr) const//const BaseList< BaseList<double> >& args) const
{
  const double tol= (args.size()==3) ? checksynctol(args(2U)) : 1e6;
  //  dosyncratio(dest,args,tol);
  const int cycles=round_int(args(1U));
  const double inratio=args.front();
  if (cycles) {
    const rational actratio=calcsync(inratio,cycles,tol);
    if (verbose & VER_GEN) {
      if (!actratio)
	std::cout << "Ratio " << inratio << " not synchronised\n";
      else
	std::cout << "Synchronises " << inratio << " to " << actratio << '\n';
    }
    if (!!actratio) {
      dest.push_back(actratio.numerator());
      dest.push_back(actratio.denominator());
      return;
    }
  }
  dest.push_back(0.0);
  dest.push_back(0.0);
}

// void FunctionSyncRatio3::operator()(LIST<double>& dest, const BaseList< BaseList<double> >& args) const
// {
//   assert(args.size()==3);
//   assert(args(2U).size()==1);
//   dosyncratio(dest,args,checksynctol(args(2U).front()));
// }

// double FunctionSyncTime::operator()(const BaseList<double>& args) const
// {
//   const double tol=(args.size()==4) ? checksynctol(args(3U)) : 1e6;
//   const int cycles=round_int(args(2U));
//   if (cycles==0)
//     return 0.0;

//   const double x=args.front();
//   const double y=args(1U);
//   const rational actratio=calcsync(x/y,cycles,tol);
//   return actratio.denominator()*y;
// }

#ifndef DISABLE_EXEC

namespace {
  void log_propagator(const BlockedMatrix<complex>& U, const char* rawtitle ="propagator")
  {
    logfile_controller* logfilep(get_logfile());
    if (!logfilep)
      std::cout << rawtitle << '\n' << U;
    else {
      char title[128];
      logfile_controller::maketitle(title,sizeof(title),rawtitle);
      logfilep->write(U,title);
    }
  }

  void accumulate_(BlockedMatrix<complex>& Utotal, const BlockedMatrix<complex>& U, BlockedMatrix<complex>& Utmp)
  {
    if (!Utotal)
      Utotal=U;
    else {
      multiply(Utmp,U,Utotal);
      Utotal.swap(Utmp);
    }
  }
}

double ActionProp::duration() const
{
  return (ntimes<0) ? propfor*1e-6 : seqp->duration();
}

ThreadWarning<> smartprop_warning("smart prop caching can't be used if combinepropagators optimisation is disabled",&NMRsim_once_warning);

ThreadWarning<> ActionProp::empty_propagator_warning("empty propagator in prop",&NMRsim_once_warning);

ThreadWarning<> weirdsync_warning("synchronisation time is not a multiple of sequence duration - ignoring",&NMRsim_repeat_warning);

ThreadWarning<> ActionProp::disablingcombine_warning("disabling propagator combination as duration is not multiple of rotation period. Sequence: ",&NMRsim_once_warning);

double redmod(double newoff, double step)
{
  double ntimes=floor((newoff+tolerance)/step);
  double resid=newoff-ntimes*step;
  return fabs(resid>tolerance) ? resid : 0.0;
}

bool ActionProp::iszeroduration() const
{
  return (ntimes<0) ? (fabs(propfor)<tolerance) : (ntimes==0);
}
  
void ActionProp::exec(MasterObj* objp, BlockedOperator& sigma, double& mastert) const
{
  if (flags_ & PROP_RESET) {
    if (seqp==NMRSIM_NULL)
      throw InternalError("ActionProp::exec");
    seqp->currentstate().reset(false);
  }
  if (iszeroduration())
    return;

  const double initt=mastert;
  const double lastoff=seqp ? seqp->currentstate().lastoffset_ : 0.0;

  const option::optional_t nocachestatus=optcache.get(); //!< store current nocache
  const bool invert=flags_ & PROP_INV;
  if (invert && interactions_MFp) {//!don't cache Us if Hamiltonian inverted
    interactions_MFp->invert_linear();
    optcache.set(option::OFF);
  }
  const bool disablecache=!optcache;

  bool finished=false;
  BlockedMatrix<complex> Utotal,Utmp;
    
  double synctime=0.0;
  size_t nsynctimes=0;
  size_t lntimes=ntimes;
  if (seqp) {

    sequencestate_t& cstate(seqp->currentstate());
    const BaseList<PhasedSequence*> curlist(seqp->currentlist());
    const CompSequenceBase& frontseq(curlist.front()->get());
    const size_t cyclelength=curlist.size();
    const bool isCW=(cyclelength==1) && (frontseq.isCW);
    const double seqdur=seqp->duration();
    if (lastoff) {
      if (ntimes>=0) {
	parser_printthread(std::cerr) << "Can't use prop with partially used sequence (" << seqp->name() << ") i.e. don't mix prop and propfor using same sequence.\n";
	error_abort();
      }
      if (isCW)
	throw InternalError("non-zero offset with CW sequence");
      if (lastoff>seqdur-tolerance)
	throw InternalError("prop: offset exceeds fragment duration");
    }
    
    if (ntimes>=0) {
      synctime=seqdur;	
      nsynctimes=ntimes / cyclelength;
    }
    else {	
      if (seqdur<=0.0) {
	parser_printthread(std::cerr)<< "trying to use propfor with sequence (" << seqp->name() << ") with zero duration!\n";
	error_abort();
      }
      double propforleft=propfor*1e-6;
      if (cyclelength==1) {
	if (isCW) {
	  if (spin_rate)
	    synctime=fabs(1.0/spin_rate);
	}
	else {
	  synctime=frontseq.synchronisation_time();
	  if ((synctime!=0.0) && (check_sync(synctime/seqdur)==0)) {
	    weirdsync_warning.raise();
	    synctime=0.0;
	  }
	}
      }
      if (synctime==0.0)
	synctime=seqdur;

      nsynctimes=size_t(0.5+floor((tolerance+propforleft)/synctime)); //!< find complete repetitions of sequence length
    }
    
    if ((verbose & VER_GEN) && (verbose_level>1))
      std::cout << "Synchronisation time: " << (synctime*1e6) << " us   Number: " << nsynctimes << "\n";
    
    if (nsynctimes) {
      const size_t repeats=nsynctimes;
      cacheset_t* cachesp=disablecache ? NMRSIM_NULL : cstate.getcache(lastoff);
      //      const BlockedMatrix<complex>* Ucachep = cachesp ? (cachesp->Ucache_).find(*objp,mastert,repeats,synctime) : NMRSIM_NULL; //!< see if cached
      const bool canuse=cachesp ? (cachesp->Ucache_).isvalid(*objp,mastert,synctime) : false;
      const BlockedMatrix<complex>* Ucachep = canuse ? &((cachesp->Ucache_)(*objp,mastert,repeats,synctime)) : NMRSIM_NULL; //!< see if cached
      if (verbose & VER_GEN)
	std::cout << "Can use already cached propagators (smartprop) for " << repeats << " repeat(s) (" << seqp->name() << "): " << (Ucachep ? "Yes\n" : "No\n");
      if (Ucachep) {
	optcache.setusage(true);
	Utotal&=*Ucachep; //!< if found, copy out
	if (lntimes>0)
	  lntimes-=repeats*cyclelength;
	mastert+=repeats*synctime;
      }
      else {
	bool usesmart=false;
	if (optsmartprop()) {
	  if (optcombinepropagators())
	    usesmart=true;
	  else
	    smartprop_warning.raise();
	}
	bool docombine=(optcombinepropagators() && (nsynctimes>=2)) || usesmart;
	optcombinepropagators.setusage(docombine);
	if (docombine && (spin_rate!=0.0) && (seqdur!=0.0)) {
	  const double rotor_period=1.0/fabs(spin_rate);
	  if ((synctime==0.0) || (check_sync(synctime/rotor_period)==0)) {
	    disablingcombine_warning.raise(seqp->name());
	    docombine=false;
	    usesmart=false;
	  }
	}
	if (docombine || usesmart) { 
	  BlockedMatrix<complex> Ubase;
	  if (ntimes<0) {
	    if (cyclelength!=1)
	      throw InternalError("ActionProp::exec - cyclelength");
	    Ubase=(curlist.front())->evaluate(objp,mastert,synctime,disablecache,lastoff);
	  }
	  else {
	    if (lastoff!=0.0)
	      throw InternalError("ActionProp::exec - lastoff");
	    double totdur=0.0;
	    for (size_t i=cyclelength;i--;) {
	      const PhasedSequence* pseqp(seqp->current());
	      const double dur=pseqp->duration();
	      if (!(pseqp->empty()) || dur) {
		accumulate_(Ubase,pseqp->evaluate(objp,mastert+totdur,dur,disablecache),Utmp);
		totdur+=dur;
	      }
	      seqp->next();
	    }
	    if (fabs(seqdur-totdur)>tolerance)
	      throw InternalError("Mismatch between cache and calculated cycle durations");
	  }
	  if (!Ubase)
	    empty_propagator_warning.raise();
    else {
      if (usesmart && !disablecache) {
        if (verbose & VER_GEN) 
          std::cout << "smartprop: initialising (" << seqp->name() << ")\n";
        optsmartprop.setusage(true);
        Usequentialcache& Useqcache(cstate.setcache(lastoff).Ucache_);
        Useqcache.initialise(*objp,mastert,Ubase,synctime);
        Utotal&=Useqcache(*objp,mastert,repeats,synctime);
      }
      else {
        pow(Utmp,Ubase,repeats);
        Utotal&=Utmp;
      }
	  }
	  mastert+=synctime*repeats;
	  lntimes-=repeats*cyclelength;
	}
	optcombinepropagators.check(docombine);
      }
    }

    if (ntimes<0) {
      const double dur=(propfor*1e-6)-(mastert-initt);
      finished=(fabs(dur)<tolerance);
      if (!finished) {
	cacheset_t* cachesp=disablecache ? NMRSIM_NULL : cstate.getcache(lastoff);
	PhasedSequence* pseqp=curlist.front();
	const BlockedMatrix<complex>* Udurp=cachesp ? (cachesp->Udurcache_)(*objp,mastert,dur) : NMRSIM_NULL;
	if (verbose & VER_GEN)
	  std::cout << "Found cache match for duration " << (dur*1e6) << " us: " << (Udurp ? "Yes\n" : "No\n");
	const double seqdur=pseqp->duration();
	if (seqdur==0.0)
	  throw InternalError("ActionProp::exec - seqdur");
	if (Udurp) {
	  optcache.setusage(true);
	  accumulate_(Utotal,*Udurp,Utmp);
	  cstate.lastoffset_=dur;
	}
	else {
	  Utmp=pseqp->evaluate(objp,mastert,dur,true,lastoff);
	  if (!disablecache)
	    cstate.setcache(lastoff).Udurcache_.store(Utmp,*objp,initt,dur);
	  Utotal&=Utmp;
	}
	finished=true;
	mastert+=dur;
	if (isCW)
	  cstate.lastoffset_=0.0;
	else {
	  const double newoff=redmod(lastoff+dur,seqdur); //!< note that "nsynctimes" component is assumed to be multiple of sequence duration and so ignored
	  if ((verbose & VER_GEN) && (verbose_level>1))
	    std::cout << "Setting sequence offset to " << (newoff*1e6) << " us\n";
	  cstate.lastoffset_=newoff;
	}
      }
    }
  }    
  if (!finished) {
    if (seqp) {
      if (ntimes<0)
	throw InternalError("ActionProp::exec - ntimes");
      for (size_t n=0;n<lntimes;n++) {
	const PhasedSequence* pseqp(seqp->current());
	const double dur=pseqp->duration();
	if (!(pseqp->empty()) || dur) { //need to check as exec may be called before H is created for null sequences
	  accumulate_(Utotal,pseqp->evaluate(objp,mastert,dur,disablecache),Utmp);
	  mastert+=dur;
	}
	seqp->next();
      }
    }
    else {
      if (ntimes>=0)
	throw InternalError("ActionProp::exec - ntimes");
      const double dur=(propfor*1e-6)-(mastert-initt); //!< mastert is not expected to have incremented, but safe
      if (fabs(dur)>tolerance) {
	objp->propagator(Utmp,mastert,mastert+dur);
	Utotal&=Utmp;
	mastert+=dur;
      }
    }      
  }

  if (flags_ & PROP_PUTMATRIX)
    log_propagator(Utotal);
  
  if (!!Utotal) {
    if ((verbose & VER_GEN) && (verbose_level>1)) {
      std::cout << "Applying overall propagator for duration ";
      prettyprint_time(mastert-initt) << '\n' << Utotal << '\n';
    }
    sigma.unitary_simtrans(Utotal);
  }

  if (invert && interactions_MFp) {//restore 
    interactions_MFp->invert_linear();  
    optcache.set(nocachestatus);
  }
}

// void ActionDirectAcq::reset()
// {
//   const int lverb=(verbose & VER_GEN) ? verbose_level : 0;
//   if (lverb>1)
//     std::cout << "Resetting direct dimension\n";
//   checksync();
// }

 ThreadWarning<> ActionAcq::excessivesync_warning("synchronisation period exceeds acquisition duration (indirect dimension)",&NMRsim_once_warning);

void ActionAcqN::reset()
{  
  //  acqp->reset(); //!< don't actually need to reset - done at each exec
}

void ActionAcq::reset() 
{
  const int lverb=(verbose & VER_GEN) ? verbose_level : 0;
  if (lverb>1)
    std::cout << "Resetting acq in dimension " << (ndim+1) << '\n';
  index=-1;
  translist.clear();
  Uacc.clear();
  double sync_period;
  nrep=checksync(sync_period);
  char buf[256];

  if (nrep>dimensionsize()) {//!< don't de-synchronise in case we benefit from other caching
    snprintf(buf,sizeof(buf),": synchronisation period in dimension %lu (%lu dwell times) exceeds length of the dimension",(unsigned long)(ndim+1),(unsigned long)nrep);
    excessivesync_warning.raise(buf);
  }
  if (lverb) {
    if (sync_period && (verbose_level>1))
      std::cout << "Acq synchronisation time: " << (sync_period*1e6) << " us\n";
    std::cout << "Acquisition RF / MAS synchronised in dimension " << (ndim+1) << ": ";
    if (nrep)
      std::cout << "Yes (x" << nrep << ")\n";
    else
      std::cout << "No\n";
  }
}

 ThreadWarning<> ActionDirectAcq::zerodurationacq_warning("sw specified with zero duration acquisition sequence",&NMRsim_repeat_warning);

void ActionDirectAcq::flush_cache(dirty_t level)
{
  dt=0.0;

  //! check key acquisition parameters in case not caught earlier
  if (sw<0.0)
    error_abort("sw cannot be <0.0");
  if (np<1)
    error_abort("np cannot be <1");
  
  if (seqp && (seqp->duration()==0.0)) {
    if (sw)
      zerodurationacq_warning.raise();
  }
  else {
    if ((sw==0.0) && (np>1))
      error_abort("sw unset (or zero) for direct acquisition with >1 data points");
  }
  
  if ((sw==0.0) && (acc_mode!=ACC_TIME))
    error_abort("frequency domain calculation requires sw to be set");
  
  dt=sw ? 1.0/sw : 0.0;

  bool needcheck=!checkedsync;
  if (seqp && ((level>DIRTY_GAMMA) || seqp->needsresync() || !optupdate)) {
    seqp->checksync(dt);
    needcheck=true;
  }
  if (needcheck) {
    checksync();
    checkedsync=true;
  }
}

// void ActionAcqBase::flush_cache_(flush_t level)
// {
//   if (seqp && phasemodp.get() && ((level>FLUSH_GAMMAONLY) || dirty_list.count(seqp))) {
//     if ((verbose & VER_GEN) && (verbose_level>1)) {
//       std::cout << "Resetting phase modulation propagator cache for ";
//       print(std::cout);
//     }
//     checkpm(seqp->duration()); //!< potential problem with acq_period defn
//     propgenp.clear(); //!< clear propagator generator
//   }
// }

void ActionAcqN::flush_cache(dirty_t level)
{
  acqp->flush_cache(level);
}

void ActionAcq::flush_cache(dirty_t level)
{
  //  flush_cache_(level);
  if (level || (seqp && seqp->wasdirty())) {
    if ((verbose & VER_GEN) && (verbose_level>1)) {
      std::cout << "Flushing propagator cache for ";
      print(std::cout);
      std::cout << '\n';
    }
    Us.clear();
  }
  //  if (seqp)
  //  checkpm(seqp->duration());
}

// void ActionDirectAcq::flush_cache(dirty_t level)
// {
//   flush_cache_(level);
//   //checksync();
// }

usage_t ActionDirectAcq::usage() const { return usage_t(); }

void ActionDirectAcq::set(double v, subsid_t subsid)
{
  switch (subsid) {
  case S_ARG1:
    if (v<0)
      error_abort("acq: synctime cannot be negative");
    synctime=1e-6*v;
    break;
  case S_PHASE:
    dphs=v;
    break;
  default:
    throw InternalError("ActionDirectAcq::set");
  }
}

void ActionAcqPoint::set(double v, subsid_t)
{
  dphs=v;
}

void ActionProp::flush_cache(dirty_t level)
{ 
  if (!seqp)
    return;

  const bool needflush=(seqp->wasdirty() || level);
  const bool needresync=!iszeroduration() && (!checkedsync || needflush);
    
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Need resync for ";
    print(std::cout);
    std::cout << ": " << (needresync ? "Yes\n" : "No\n");
  }
  
  if (needresync) {
    seqp->checksync(duration());    
    checkedsync=false;
  }
}

ThreadWarning<> ActionAcq::activetransients_warning("transients before acq",&NMRsim_repeat_warning);
ThreadWarning<InternalError> ActionAcq::improperreuse_warning("caught improper propagator re-use",&NMRsim_repeat_warning);
ThreadWarning<> ActionAcq::trailingtransients_warning("trailing transients created during acq",&NMRsim_once_warning);
ThreadWarning<> ActionAcq::emptypropagator("indirect dimension propagator doesn't seem to do anything!",&NMRsim_once_warning);
ThreadWarning<> ActionAcqN::fullFID_warning("acquiring complete FID using acqn is inefficient - use normal acq if possible",&NMRsim_once_warning);

void ActionAcqPoint::exec(MasterObj* objp, BlockedOperator& sigma, double& t) const
{
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Acquiring single data point at t=";
    prettyprint_time(t) << '\n';
  }
  objp->raw_add(*this,sigma);
}

void ActionAcqN::exec(MasterObj* objp, BlockedOperator& sigma, double& t) const
{
  const double dt=acqp->dwelltime();
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Acquiring " << toacquire_ << " points at t=";
    prettyprint_time(t) << " with dwell time of " << (dt*1e6) << " us\n";
  }
  acqp->reset(); //!< need to reset Acq each time
  acqp->transientcleanup(sigma);
  BlockedOperator lsigma(sigma);
  if (toacquire_==np)
    fullFID_warning.raise();
  for (size_t i=0;i<toacquire_;i++) {
    acqpointp->exec(objp,sigma,t);
    acqp->increment(objp,t,i);
    sigma=lsigma;
    acqp->applyincrement(sigma,beverbose && (i==0));
  }  
  t+=toacquire_*dt;
}

void ActionAcq::applyincrement(BlockedOperator& sigma, bool verb) const
{
  if (verb)
    log_propagator(Uacc);
  if (!Uacc)
    emptypropagator.raise();
  else
    sigma.unitary_simtrans(Uacc);
}

void ActionAcq::transientcleanup(BlockedOperator& sigma)
{
  if (transients_active()) {
    activetransients_warning.raise();
    
    BlockedMatrix<complex> Utmp;
    
    flush_transients(Utmp);
    if (!!Utmp) {
      log_propagator(Utmp,"transient_propagator");
      sigma.unitary_simtrans(Utmp);
    }
  }
}

void ActionAcq::increment(MasterObj* objp, double t, int oldindex) const
{
  const double dt=dwelltime();
  BlockedMatrix<complex> Utmp;
  
  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  
  const double starttime=t+oldindex*dt;
  const double endtime=starttime+dt;
  BlockedMatrix<complex>* Up=NMRSIM_NULL;
  
  double junk;
  if (checksync(junk)!=nrep)
    error_abort("synchronisation of indirect dimension has changed between increments");

  if (nrep && optcache()) {
    
    if (Us.empty()) {
      Us.create(nrep);
      phases.create(nrep);
    }
    
    const size_t cacheindex=oldindex % nrep;
    Up=&(Us(cacheindex));
    
    const double rphase=get_rotor_phase(starttime);
    if (!(Up->empty()) && (fabs(rphase-phases(cacheindex))>1e-5)) {
      Up->clear();
      if (lverbose)
	improperreuse_warning.raise();
    }
    if (lverbose) {
      std::cout << (Up->empty() ? "Creating" : "Using") << " cached propagator from t=";
      prettyprint_interval(starttime,endtime) << '\n';
    }
    
    if (Up->empty()) {
      if (seqp) {
	objp->propagator(*Up,*seqp,starttime,endtime,starttime);
	if (transients_active()) {
	  trailingtransients_warning.raise();
	  flush_transients(*Up); //incorporate transient
	}
      }
      else
	objp->propagator(*Up,starttime,endtime);
      phases(cacheindex)=rphase;
    }
    else
      optcache.setusage(true);
  }
  else {
    Up=&Utmp;
    if (lverbose) {
      std::cout << "Determining acquisition propagator from ";
      prettyprint_interval(starttime,endtime) << '\n';
    }
    if (seqp) {
      restore_transients();
      objp->propagator(*Up,*seqp,starttime,endtime,starttime);
      store_transients(endtime);
    }
    else
      objp->propagator(*Up,starttime,endtime);
  }
  if (lverbose>1)
    std::cout << (*Up) << '\n';
  
  //propagate
  if (!(Up->empty()))
    Uacc&=*Up;
}

void ActionAcq::exec(MasterObj* objp, BlockedOperator& sigma, double& t) const
{
  if (insidecontrol)
    error_abort("acq cannot be inside control / loop structure - did you mean to use acqn or acqpoint?");

  const size_t newindex=getarrayindex(ndim);
  const double endtime=t+(index+1)*dwelltime();

  const bool isverb=((verbose & VER_GEN) && (verbose_level>1));
  if (isverb)
    std::cout << "Dimension: " << ndim << "  Current index: " << index << "  New index: " << newindex << std::endl;

  if (newindex!=index) {
    if (newindex!=index+1)
      error_abort("ActionAcq: out of order increment");
    index=newindex;
    
    if (isverb)
      std::cout << "Incrementing dimension " << (ndim+1) << '\n';
          
    if (index==0)
      return; //first point in t1 (do nothing)

    transientcleanup(sigma);
    increment(objp,t,index-1);
  }

  if (index==0)
    return; //first point in t1 (do nothing)

  applyincrement(sigma,beverbose);

  t=endtime;
}

void ActionPutMatrix::print(std::ostream& ostr) const
{
  ostr << "putmatrix " << name;
  if (flags_) {
    ostr << ' ';
    dumpflags(ostr,putmatrix_flags(),flags_);
    if (statssel.first)
      ostr << " <spin index selection>";     //!< can't dump stats arg as compressed!
    if (statssel.second)
      ostr << " <filter matrix>"; // can't output filter matrix as name unknown!
  }  
}

void ActionPutMatrix::exec(MasterObj* Objp, BlockedOperator&, double&) const
{
  if (done_ && (flags_ & PM_ONCE))
    return;
  dump_matrix(curm,name.c_str(),Objp,flags_,statssel);
  done_=true;
}

void ActionFilter::exec(MasterObj*, BlockedOperator& sigma, double&) const
{
  sigma.apply(filter());
}

void ActionRotorAngle::exec(MasterObj* objp, BlockedOperator&, double&) const {
  objp->reset_angle(rotor_angle*deg_to_rad);
}

actionstack_t::actionstack_t()
{
  local_acq_count=0; //!< reset count of acquisitions with new action stack
  acq_count=-1; //!< flag uninitialised
}

void actionstack_t::initialise()
{
  acq_count=local_acq_count; //!< store overall acq_count
}

void actionstack_t::execaction(ActionCommand* actionp, MasterObj* objp, BlockedOperator& dens, double& t, size_t& count)
{
  timer_guard_t guard(actionp->stopwatch);
  actionp->exec(objp,dens,t);
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Density matrix after step " << count << " (t=";
    prettyprint_time(t) << ")\n" << dens;
  }
  ++count;
}

void actionstack_t::exec(MasterObj* objp, BlockedOperator& dens,double& t)
{
  const actionstack_t::const_iterator finish(end());
  actionstack_t::const_iterator start(begin());
  const bool dumpsigma=(verbose & VER_GEN) && (verbose_level>1);
  if (dumpsigma) {
    std::cout << "Density matrix at t=";
    prettyprint_time(t) << " (start)\n" << dens;
  }
  if (!dens)
    throw InternalError("Starting pulseq block with empty density matrix!");
  ActionDirectAcq::reset_acquisition_count();
  size_t count=1;
  while (start!=finish) {
    ActionControl* controlp=dynamic_cast<ActionControl*>(*start);
    if (controlp) {
      insidecontrol=true;
      ++start;
      controlp->exec_control(start,finish,count,objp,dens,t);
      start+=controlp->items();
    }
    else {
      insidecontrol=false;
      execaction(*start,objp,dens,t,count);
      ++start;
    }
  }
}

ContextWarning<> ActionDo::zeroloops_warning("zero iterations of loop!",&NMRsim_once_warning);

ActionCommand* ActionDo::create()
{
  Variable cvar(S_ARG1);
  const size_t ntimes=parse_unsigned(&cvar);  
  if (!nochecks && (cvar.subsid==S_NONE) && (ntimes==0))
    zeroloops_warning.raise();
  return new ActionDo(ntimes);
}

void ActionDo::set(double timesv,subsid_t)
{
  const int timesi=round_int(timesv);
  if (timesi<0) {
    parser_printcontext() << "do loop iterations cannot be negative (" << timesi << ")\n";    
    error_abort();
  }
  if (!nochecks && (timesi==0))
    zeroloops_warning.raise();
  times_=timesi;
}

void ActionDo::exec_control(actionstack_t::const_iterator start, const actionstack_t::const_iterator& end, size_t& count, MasterObj* objp, BlockedOperator& dens,double& t) const
{  
  for (size_t n=times_;n--;) {    
    actionstack_t::const_iterator lstart(start);
    for (size_t i=items();i--;) {
      if (lstart==end) 
	throw InternalError("ActionDo");
      actionstack_t::execaction(*lstart,objp,dens,t,count);
      ++lstart;
    }
  }
}
    
void ActionEcho::exec(MasterObj*, BlockedOperator&, double&) const 
{ 
  BaseEcho::exec(); 
}

ContextWarning<> ActionFilter::filterafteracq_warning("filter after final acq",&NMRsim_repeat_warning);

const filter_def& getfiltermatrix(const char* name)
{
  const matrixspec& curm(findmap(matrixmap,name,"matrix"));
  if (curm.type!=matrixspec::BOOL) {
    parser_printcontext() << name << " is not a coherence mask\n";
    error_abort();
  }
  return curm.asfilter();    
}

ContextWarning<> ActionFilter::notcoherencefilter_warning("filter used with matrix that is not a coherence matrix filter which is not physical (disable warning with -nochecks)",&NMRsim_once_warning);

ActionCommand* ActionFilter::create()
{
  if (local_acq_count>nacqdims)
    filterafteracq_warning.raise();
  const char* name=parse_string(F_REPLACEDOLLAR);
  const filter_def& filterspec(getfiltermatrix(name));
  if (!nochecks && filterspec.isspinorder())
    notcoherencefilter_warning.raise();
  return new ActionFilter(name,filterspec);
}

ActionCommand* ActionPutMatrix::create()
{
  const char* thisname(parse_string(F_REPLACEDOLLAR));
  const matrixspec& curm(findmap(matrixmap,thisname,"matrix"));
  const std::pair<int,statistics_t> flaginfo(get_putmatrix_flags());
  return new ActionPutMatrix(curm,thisname,flaginfo.first, flaginfo.second);
}

void ActionProp::reset()
{
  if (seqp) {
    seqp->reset();
    if (ntimes<0)
      checkedsync=false;
    seqp->flush(); //!< can't assume always flushed by dirty flag (e.g. complete action stack flush)
  }
}

void ActionScale::print(std::ostream& ostr) const 
{
  ostr << "scale " << scale;
  if (filterp)
    ostr << ' ' << name;
}

void ActionScale::printvariablename(std::ostream& ostr, subsid_t) const 
{
  ostr << "scale";
}

ContextWarning<> ActionScale::coherencefilter_warning("scale used with coherence filter matrix which is not very physical (disable warning with -nochecks)",&NMRsim_once_warning);

ActionCommand* ActionScale::create()
{
  Variable cvar(S_ARG1);
  const double scale=parse_double(&cvar,F_DENYZERO);
  if (cvar.subsid)
    havescale=true;      

  if (are_left()) {
    const char* name=parse_string(F_REPLACEDOLLAR);
    const filter_def& filterspec(getfiltermatrix(name));
    if (!nochecks && !(filterspec.isspinorder()))
      coherencefilter_warning.raise();
    return new ActionScale(scale,name,&(filterspec));
  }
  return new ActionScale(scale);
}

void ActionScale::exec(MasterObj*, BlockedOperator& sigma, double&) const
{
  if (filterp)
    sigma.scale(scale,(*filterp)());
  else
    sigma*=scale;
}

#endif
