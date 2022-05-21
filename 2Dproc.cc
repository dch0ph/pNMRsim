/*! \file
  \brief  Processing directives for 2D processing
*/

#define NMRSIM_AUTOPHASE_DEFAULT_MA 10

#include "NMRsim_common.h"
#include "NMRsim_Process.h"
static const double deg_to_rad=M_PI/180.0;

class ProcessSplitSinCos : public ProcessCommand {
public:
  ProcessSplitSinCos()
    : ProcessCommand(PROC_HAS2D | PROC_CHANGES_ROWS | PROC_CHANGES_COLUMNS) {}
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create() { return new ProcessSplitSinCos(); }
  void print(std::ostream& ostr) const {  ostr << "splitsincos"; }
  static ThreadWarning<> fd_warning; //!< applying split to frequency domain data
};

class ProcessFTIndirect : public ProcessCommand {
public:
  ProcessFTIndirect(size_t flagsv)
    : ProcessCommand(PROC_HAS2D), flags_(flagsv) {}
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  void print(std::ostream&) const;
private:
  size_t flags_;
};

#ifdef USE_MINUIT
class ProcessAutoPhase : public ProcessCommand {
public:
  ProcessAutoPhase(int flagsv =0);
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  void print(std::ostream&) const;
  static ThreadWarning<> tddata_warning;

  enum { ZEROONLY=1, IS2D =2};

private:
  bool zeroonly; //!< only optimise zero-order phase
};
#else
class ProcessAutoPhase {
public:
  static ProcessCommand* create() {
    error_abort("autophase unsupported (requires MINUIT)");
    return NMRSIM_NULL;
  }
};
#endif

class ProcessBaselineCorrect : public ProcessCommand {
public:
  ProcessBaselineCorrect(const BaseList<size_t>& indsv)
    : ProcessCommand(PROC_HAS1D), inds_(indsv) {
    if (inds_.empty())
      throw InvalidParameter("ProcessBaselineCorrect()");
  }
  static ProcessCommand* create();
  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream& ostr) const { ostr << "baselinecorrect " << inds_; }
private:
  LIST<size_t> inds_;
};

class ProcessSum : public ProcessCommand {
public:
  ProcessSum()
    : ProcessCommand(PROC_HAS2D | PROC_CHANGES_ROWS) {}
  ProcessSum(const BaseList<size_t>& rowselv)
    : ProcessCommand(PROC_HAS2D | PROC_CHANGES_ROWS), rowsel_(rowselv) {}

  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  void print(std::ostream& ostr) const { ostr << "sum"; }
private:
  LIST<size_t> rowsel_;
};

class ProcessAdd : public ProcessCommand {
public:
  ProcessAdd(const char*, double =1.0);
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  void print(std::ostream&) const;
private:
  std::string name_;
  cmatrix data_;
  filestruct pars_;
  double scale_;
};

namespace {
  struct Proxy_ {
    Proxy_() {
      Process_Factory_t& factory(get_Process_Factory());

      factory["autophase"]=&ProcessAutoPhase::create;
      factory["splitsincos"]=&ProcessSplitSinCos::create;
      factory["ftindirect"]=&ProcessFTIndirect::create;
      factory["add"]=&ProcessAdd::create;
      factory["baselinecorrect"]=&ProcessBaselineCorrect::create;
      factory["integrate"]=&ProcessIntegrate::create;
      factory["sum"]=&ProcessSum::create;
      factory["system"]=&ProcessSystem::create;
    }
  };
  
  Proxy_ proxy;
}

ProcessAdd::ProcessAdd(const char* namev, double scalev)
  : ProcessCommand(PROC_HAS2D),
    name_(namev), scale_(scalev)
{
  raw_read_file(data_,pars_,namev);
}

void ProcessAdd::print(std::ostream& ostr) const
{
  ostr << "add " << name_;
  if (scale_!=1.0)
    ostr << ' ' << scale_;
}

ThreadWarning<> combine_sw_warning("apparently mismatched spectral widths",&NMRsim_repeat_warning);
ThreadWarning<> combine_domain_warning("current and added data set apparently in different domains - dimension: ",&NMRsim_repeat_warning);

void checkdomainsw(const processing_state& pflag, double sw, domain_t dom, const char* domname)
{
  if (sw && pflag.sw) {
    if (fabs(sw-pflag.sw)>1e-5*sw) {
      char buf[256];
      snprintf(buf,sizeof(buf),": sw of current data set=%g kHz, sw of added data set=%g kHz (%s dimension)",pflag.sw*1e-3,sw*1e-3,domname);
      combine_sw_warning.raise(buf);
    }
  }
  if (dom && (pflag.istimedomain!=(dom==D_TIME)))
    combine_domain_warning.raise(domname);
}

ProcessCommand* ProcessSum::create()
{
  if (are_left()) {
    static const char synstr[]="[<row selection>]";
    const LIST<size_t> row_sel(parse_unsignedintarray_syntax(synstr,0,1));
    if (row_sel.empty())
      error_abort("sum: row selection is empty!");
    return new ProcessSum(row_sel);
  }
  return new ProcessSum();
}

void ProcessSum::exec(cmatrix& a, LIST<processing_state>& pflags) const
{
  const size_t nr=a.rows();
  if (!(rowsel_.empty()))
    verify_index_list(rowsel_,nr,"row");

  if (nr<2)
    return;
  
  cmatrix d(size_t(1),a.cols(),complex(0.0,0.0));
  BaseList<complex> drow(d.row());

  if (rowsel_.empty()) {
    for (size_t i=nr;i--;)
      drow+=a.row(i);
  }
  else {
    for (size_t i=rowsel_.size();i--;) {
      const size_t r=rowsel_(i);
      drow+=a.row(r-1); //!< already verified that indices are OK
    }
  }
  a.swap(d);

  //  for (size_t ndim=1;ndim<pflags.size();ndim++) {
  //  array_dims.set(dim,n,ni_skip);
  //  raw_set_n(ndim+1,1,1);
  //}
  array_ns.clear();
  array_n0=1;

  const processing_state colstate(pflags.back()); //!< explicit copy for safety
  pflags.create(size_t(1),colstate);
}

void ProcessAdd::exec(cmatrix& a, LIST<processing_state>& pflags) const
{
  if (!arematching(a,data_)) {
    std::cerr << "Can't add data set " << name_ << " (" << data_.rows() << 'x' << data_.cols() << ") to current data set (" << a.rows() << 'x' << a.cols() <<")\n";
    error_abort();
  }
  if (scale_==1.0)
    a+=data_;
  else
    mla(a,scale_,data_);

  checkdomainsw(pflags.back(),pars_.sw,pars_.domain,"direct");
  if (a.rows()>1)
    checkdomainsw(pflags.front(),pars_.sw1,pars_.domain1,"indirect");
}

ProcessCommand* ProcessAdd::create()
{
  const char* fname=parse_string(F_REPLACEDOLLAR); //!< read in immediately - consider delaying $ parsing
  const double scale=are_left() ? parse_double() : 1.0;
  return new ProcessAdd(fname,scale);
}

void ProcessFTIndirect::print(std::ostream& ostr) const
{
  ostr << "ftindirect ";
  print_ftflags(ostr,flags_);
}

ProcessCommand* ProcessFTIndirect::create()
{
  return new ProcessFTIndirect(get_ftflags());
}

ProcessCommand* ProcessBaselineCorrect::create()
{
  static const char selsyn[]="baselinecorrect " NMRSIM_RANGESTR;
  const LIST<size_t> sel(parse_unsignedintarray_syntax(selsyn,0,1));
  if (sel.empty())
    error_abort("baselinecorrect: column selection cannot be empty");
  return new ProcessBaselineCorrect(sel);
}

void ProcessBaselineCorrect::exec(BaseList<complex> data, processing_state&) const
{
  LIST<complex> sel;
  get_selection(sel,data,inds_);  
  const complex offset=sum(sel)/sel.size();
  data-=offset;
}

ThreadWarning<> ProcessSplitSinCos::fd_warning("splitsincos applied to frequency domain data",&NMRsim_once_warning);

void ProcessSplitSinCos::exec(cmatrix& a, LIST<processing_state>& pflags) const
{
  const size_t nrows=a.rows();
  const size_t ncols=a.cols();
  if (ncols & 1)
    error_abort("splitsincos applied to data with odd number of columns");

  if (!(pflags.back().istimedomain))
    fd_warning.raise();

  const size_t newcols=ncols/2;
  cmatrix d(nrows*2,newcols);
  const range realrange(0,newcols-1);
  const range imagrange(newcols,ncols-1);
  for (size_t r=nrows;r--;) {
    BaseList<complex> destrowr(d.row(2*r));
    const BaseList<complex> sourcerow(a.row(r));
    destrowr=sourcerow(realrange);
    BaseList<complex> destrowi(d.row(2*r+1));
    destrowi=sourcerow(imagrange);
  }
  a.swap(d);
  skips.front()=2;
}

void ProcessFTIndirect::exec(cmatrix& a, LIST<processing_state>& pflags) const
{
  const size_t nrows=a.rows();
  const size_t nrows2=nrows/2;
  LIST<complex> tmpr;
  LIST<complex> tmpi(nrows2);
  FTObject obj(nrows,flags_);

  processing_state& curflags=pflags.front();

  if (verbose & VER_GEN)
    std::cout << "Performing ftindirect over " << nrows << " data rows assuming " << ((skips.front()==2) ? "amplitude (States)" : "phase") << " modulated data.\n";

  for (size_t c=a.cols();c--;) {
    switch (skips.front()) {
    case 1: 
      tmpr.create(nrows);
      for (size_t r=nrows;r--;)
	tmpr(r)=a(r,c);
      obj.exec(tmpr,curflags);
      for (size_t r=nrows;r--;)
	a(r,c)=tmpr(r);
      break;

    case 2: {
      tmpr.create(nrows2);
      size_t r,hr;
      for (r=hr=0;hr<nrows;r++,hr+=2) {
	tmpr(r)=complex(a(hr,c).real(),a(hr+1,c).real());
	tmpi(r)=complex(a(hr,c).imag(),a(hr+1,c).imag());
      }
      obj.exec(tmpr,curflags);
      obj.exec(tmpi,curflags);
      for (r=hr=0;hr<nrows;r++,hr+=2) {
	a(hr,c)=complex(tmpr(r).real(),tmpi(r).real());
	a(hr+1,c)=complex(tmpr(r).imag(),tmpi(r).imag());
      }
    }
      break;

    default:
      error_abort("ftindirect - skip can only be 1 or 2");
    }
  }
  curflags.istimedomain=!curflags.istimedomain;
}

#ifdef USE_MINUIT
#if MINUITVER==2
#include "Minuit2/MnStrategy.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#else
#include "Minuit/MnStrategy.h"
#include "Minuit/FunctionMinimum.h"
#include "Minuit/VariableMetricMinimizer.h"
#include "Minuit/SimplexMinimizer.h"
#endif

class MinDataObj : public BaseMinFunction {
 public:
  MinDataObj(const BaseList<complex>& datav)
    : data(datav) { reset(); }
  double operator()(const BaseList<double>& pars) const;
  void reset() { fcncount=0; }

 private:
  mutable size_t fcncount;
  const BaseList<complex>& data;
};

double MinDataObj::operator()(const BaseList<double>& pars) const
{
  const int lverbose = (verbose & VER_OPTIM) ? verbose_level : 0;

  fcncount++;

  const double p0=deg_to_rad*pars.front();
  const double p1=0.0;//deg_to_rad*pars(size_t(1));

  double sumreal=0.0;

  const double scale=1.0/data.size();
  for (size_t i=data.size();i--;) {
    const double phasec=phase_correction(i*scale,p0,p1);
     sumreal+=real(data(i)*expi(phasec));
  }
  
  if (lverbose)
     std::cout << fcncount << ": " << pars << "\t " << sumreal << '\n';
  return -sumreal;  
}

ThreadWarning<> ProcessAutoPhase::tddata_warning("autophase applied to time-domain data",&NMRsim_once_warning);

ProcessAutoPhase::ProcessAutoPhase(int flagsv)
  : ProcessCommand( (flagsv & IS2D) ? PROC_HAS2D : PROC_HAS1D),
  zeroonly(flagsv & ZEROONLY) {}

ProcessCommand* ProcessAutoPhase::create()
{
  static flagsmap_type autophase_flags;
  if (autophase_flags.empty()) {
    autophase_flags["zeroonly"]=ZEROONLY;
    autophase_flags["twod"]=IS2D;
  }
  const int flags=parse_flags(autophase_flags);
  return new ProcessAutoPhase(flags);
}

void ProcessAutoPhase::print(std::ostream& ostr) const
{ 
  ostr << "autophase";
  if (zeroonly)
    ostr << " -zeroonly";
}

void ProcessAutoPhase::exec(cmatrix& data, LIST<processing_state>& pflags) const
{
  exec(data.row(),pflags.back());
}

void ProcessAutoPhase::exec(BaseList<complex> data, processing_state& pflag) const
{
  if (pflag.istimedomain)
    tddata_warning.raise();

  LCM_MINUITNAMESPACE::MnUserParameters mnparas;
  static bool doneinit=false;

  static smartptr<LCM_MINUITNAMESPACE::ModularFunctionMinimizer,false> minimizerp;

  if (!doneinit) {
    const double err=10.0;
    
#if MINUITVER==2
    mnparas.Add("p0",0.0,err);
    if (!zeroonly)
      mnparas.Add("p1",0.0,err);
#else
    mnparas.add("p0",0.0,err);
    if (!zeroonly)
      mnparas.add("p1",0.0,err);
#endif

    minimizerp.reset(new LCM_MINUITNAMESPACE::SimplexMinimizer); //!< Simplex minimiser is more stable in the presence of experimental noise
  }

  MinDataObj minobj(data);
  MinuitAdaptor<MinDataObj> theFCN(minobj);

  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;

  //do optim
  const LCM_MINUITNAMESPACE::MnStrategy mnstrategy(1);
  const int maxevals=200;
  const double tol=1e-1; //convert 0..1 tolerance into scale that Minuit works with
#if MINUITVER>1
  LCM_MINUITNAMESPACE::FunctionMinimum min(minimizerp->Minimize(theFCN,mnparas,mnstrategy,maxevals,tol));
  const int iters=min.NFcn();
  const bool isvalid=min.IsValid();
#else
  LCM_MINUITNAMESPACE::FunctionMinimum min(minimizerp->minimize(theFCN,mnparas,mnstrategy,maxevals,tol));
  const int iters=min.nfcn();
  const bool isvalid=min.isValid();
#endif
  if (isvalid || (iters>maxevals)) {
    if (lverbose) {
      if (isvalid)
	std::cout << "Optimisation converged after " << iters << " evaluations\n";
      else
	std::cout << "Optimisation unfinished after maximum evaluations\n";
    }

    const LCM_MINUITNAMESPACE::FunctionMinimum minvals(min);
    double p1=0;

#if MINUITVER==2
    const LCM_MINUITNAMESPACE::MnUserParameterState& parstate(minvals.UserState());  
    const double p0=parstate.Value( (unsigned int)0);
    if (!zeroonly)
      p1=parstate.Value(1);
#else
    const LCM_MINUITNAMESPACE::MnUserParameterState& parstate(minvals.userState());  
    const double p0=parstate.value( (unsigned int)0);
    if (!zeroonly)
      p1=parstate.value(1);
#endif
    if (!silent) {
      std::cout << "Autophase successful with zero-order phase: " << p0 << " deg";
      if (zeroonly)
	std::cout << '\n';
      else
	std::cout << ", first-order phase: " << p1 << " deg\n";
    }
    phase_correct_ip(data, deg_to_rad*p0, deg_to_rad*p1);
  }
  else
    std::cerr << "Optimisation did not converge or failed (autophase ignored)\n";
}

#endif
// end of MINUIT only functions


ContextWarning<> ProcessIntegrate::overlap_warning("integration regions overlap - this is probably an error (suppress this warning with -nochecks)",&NMRsim_once_warning);

namespace {
  bool isoverlap(size_t x1, size_t x2, size_t y1, size_t y2) {
    return (x1 <= y2) && (y1 <= x2);
  }
}

ProcessIntegrate::ProcessIntegrate(UserVariable& destv, const BaseList<size_t>& rangesv) :
    ProcessCommand(PROC_HAS1D),
    dest_(destv), ranges_(rangesv) 
{
  const size_t rlen=ranges_.size();
  if (rlen & 1)
    error_abort("odd number of arguments to integrate - must be <start> <end> pairs");
  maxindex_=0;
  for (size_t i=0;i<rlen;i+=2) {
    if (ranges_(i+1)<ranges_(i))
      std::swap(ranges_(i),ranges_(i+1));
    if (ranges_(i+1)>maxindex_)
      maxindex_=ranges_(i+1);
  }
  if (!nochecks) {
    bool isproblem=false;
    for (size_t i=0;i<rlen;i+=2) {
      const size_t starti=ranges_(i);
      const size_t endi=ranges_(i+1);
      for (size_t j=i+2;j<rlen;j+=2) {
	if (isoverlap(starti,endi,ranges_(j),ranges_(j+1))) {
	  isproblem=true;
	  break;
	}
      }
    }
    if (isproblem)
      overlap_warning.raise();
  }
}

ProcessCommand* ProcessIntegrate::create()
{
  UserVariable& var=*findnewvariable(parse_variable_name());
  const size_t nleft=count_left();
  LIST<size_t> rangevs;
  static const char selsyn[]="integrate " NMRSIM_RANGESTR;
  if (nleft==1)
    rangevs=parse_unsignedintarray_syntax(selsyn,0,1);
  else {
    rangevs.reserve(nleft);
    for (size_t i=0;i<nleft;i++)
      rangevs.push_back(parse_unsigned_offset(1));
  }
  return new ProcessIntegrate(var,rangevs);
}

void ProcessIntegrate::print(std::ostream& ostr) const
{
  ostr << "integrate " << dest_.name();
  for (size_t j=0;j<ranges_.size();j++)
    ostr << ' ' << ranges_(j);
}

void ProcessIntegrate::exec(BaseList<complex> dest, processing_state&) const
{
  if (!(dest_.hasbase()))
    throw InternalError("integrate: destination variable has no data store");

  if (dest.size()<=maxindex_) {
    parser_printthread(std::cerr) << "Integration involves indices up to " << (maxindex_+1) << " while data set only has " << dest.size() << " points\n";
    error_abort();
  }

  VariableBase& varbase(dest_.value());

  const size_t nranges=ranges_.size()/2;

  if (!(varbase.isarray())) {
    if (array_n0<1)
      throw InternalError("integrate: array_n0 unset");
    varbase.initialise(array_n0,nranges);
  }

  if (row_index>=array_n0)
    throw InternalError("integrate: row index out of range");

  tmp_.create(nranges);
  size_t j=0;
  for (size_t i=0;i<nranges;i++) {
    double v=0.0;
    const size_t start=ranges_(j++);
    const size_t end=ranges_(j++);
    for (size_t k=start;k<end;k++)
      v+=real(dest(k));
    tmp_(i)=v;
  }
  varbase.set_current_row(tmp_);
}
