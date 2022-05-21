// "Module" for fitting

#include <sstream>
#include "NMRsim_MasterObj.h"
#include "fit.h"
#include "NMRsim_Process.h"
#include "Interaction.h"
#include "Parser.h"
#include "optim.h"
#include "ttyio.h"
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
#endif

DataStore residuals,fit_set;

void getpars(LIST<VarVariable*>&);
void ensure_named();
void ensure_dataset(const DataStore&);
void dump_parameters(bool onlyvars =false, bool ignoreconstraints =false);
void calc_parameters(rmatrix&);
//bool calc_statistics(rmatrix&, bool ignoreconstraints, bool lsilent);
void display_statistics(bool ignoreconstraints =false);
size_t prepare_paras(bool lsilent, bool nomods =false);

bool used_reverse=false; //!< true if -reversefrequency flag has ever been used

static double minuit_precision=0.0; //!< explicit eps (0 if using default)

LIST<Parameter> parameter_list;
VariableBase para_store,error_store;

command_Factory_t optimiseblock_Factory;

void doupdate(VariableBase& var, BaseList<double> tmp, size_t col)
{
  if (tmp.size()!=valerrs.rows())
    throw InternalError("doupdate");
  tmp=valerrs(range(),col);
  var=tmp;
}
    
//! copy out parameter values and errors
void updatevarerrors()
{
  ScratchList<double> tmp(valerrs.rows());
  doupdate(para_store,tmp,0U);
  doupdate(error_store,tmp,1U);

  static bool donecreate=false;
  if (!donecreate) {
    (void)UserVariable::create("final_parameters",para_store,UserVariable::IGNORE_UNUSED);
    (void)UserVariable::create("final_errors",error_store,UserVariable::IGNORE_UNUSED);
    donecreate=true;
  }
}

class FitObj : public BaseFitFunction {
 public:
  FitObj(MasterObj& obj_)
    : obj(obj_) {}
  void operator()(BaseList<double> dest, const BaseList<double>& pars) const;

 private:
  MasterObj& obj;
  mutable DataStore tmp;
};

static LIST<cmatrix> fit_sets;

namespace {
  enum { STAT_RELEASE=1, STAT_SILENT=2, STAT_IFCONS=4 };
  double global_chisqr=0;
  FitObj* fitobjp=NMRSIM_NULL;
#ifdef USE_MINUIT
  MasterObj* objp;
  smartptr<LCM_MINUITNAMESPACE::FunctionMinimum,false> globalminp;
  bool min_funcset=false;
  OptimiseFunction min_externalp=NMRSIM_NULL;
#endif
  int dirtystart=-1; //!< pointer to un-named variables in fitting vars stack (-1 nothing to do)
  bool haverun=false;
  bool haveresults=false;

  //! Could be more intelligent and cache number of actives, but time saving minimal
  size_t countactive() {
    size_t active=0;
    for (size_t i=var_varpars.size();i--;) {
      if (!(var_varpars(i)->isfixed()))
	active++;
    }
    return active;
  }

  enum {
    OPT_GRADIENT =1,
    OPT_SIMPLEX =2,
    OPT_RANDOMISE =4,
    OPT_COMPLEX=8,
    OPT_REAL=16,
    OPT_IMAG=32,
    OPT_NOCONSTRAINTS=64
  };

  //! store variable names
  void getvarnames()
  {
    varnames.create(var_varpars.size());
    for (size_t i=varnames.size();i--;) 
      varnames(i)=var_varpars(i)->name();
  }

  const int opt_types(OPT_GRADIENT | OPT_SIMPLEX);
  const int complexity_types(OPT_COMPLEX | OPT_REAL | OPT_IMAG);

  int fit_verbose_level()
  {
    if (silent)
      return 0;
    const int fitverb = (verbose & VER_OPTIM) ? verbose_level : 0;
    return (fitverb>2) ? 2 : fitverb; //level 3 gives too much
  }

  LIST<double> noisevals;
  LIST<double> noisevector;
  LIST<size_t> mask;
}

void build_fitobjs(MasterObj&);
void ensure_noisevector(size_t);

//! optimisation mode
enum optim_t { OPTIM_NONE=0, //!< none active
	       OPTIM_FIT, //!< fitting
	       OPTIM_OPT  //!< max/minisation
};

struct rawpara {
  double value;
  double step;
  bool isconst;
};

class FitNormalise : public FitCommand
{
public:
  FitNormalise() {}
  void exec() const {} //!< don't actually do anything here.  Object created earlier and normalisation applied explicitly
  static FitCommand* create();
  static void apply(DataStore&); //!< apply normalisation
  static bool isactive() { return !!normaliseptr; }
private:
  static smartptr<ProcessNormalise,false> normaliseptr;
};

class FitSystem : public FitCommand
{
public:
  FitSystem(const char* comv)
    : com_(comv) {}
  void exec() const;
  static FitCommand* create();
private:
  std::string com_;
};

class FitMask : public FitCommand
{
public:
  FitMask(const BaseList<size_t>& indv =BaseList<size_t>())
    : inds_(indv) {}
  static FitCommand* create();
  void exec() const { mask=inds_; }

private:
  LIST<size_t> inds_;
};

void validatesel(LIST<size_t>& sel, LIST<bool>& ind, const char* name, size_t which, size_t limit)
{
  static const char* synstr=NMRSIM_ROWCOLSTR;
  try {
    sel=parse_unsignedintarray_syntax(synstr,which,1);
    ind.create(limit,true);
    for (size_t i=sel.size();i--;) {
      const size_t curind=sel(i);
      if (curind>=limit) {
	parser_printcontext() << "index (" << (curind+1) << ") exceeds data range: " << limit << '\n';
	error_abort();
      }
      bool& cur(ind(curind));
      if (cur)
	cur=false;
      else {
	parser_printcontext() << "index selection contains repeated elements: " << sel << '\n';
	error_abort();
      }
    }
  } catch (MatrixException& exc) {
    parser_printcontext() << exc << '\n';
    error_abort();
  }
}

void inversemask(LIST<size_t>& sel, const BaseList<bool>& cin)
{
  const size_t cols=cin.size();
  LIST<size_t> isel1(cols-sel.size());
  isel1.create(size_t(0));
  for (size_t i=0;i<cols;i++) {
    if (cin(i))
      isel1.push_back(i);
  }
  sel.swap(isel1);
}

void addmask(LIST<size_t>& lmask, const BaseList<size_t>& inds, size_t cols, size_t base, bool rev)
{
  if (!rev && (base==0)) {
    lmask.push_back(inds);
    return;
  }
  const size_t ninds=inds.size();
  lmask.reserve(lmask.size()+ninds);
  if (rev) {
    for (size_t i=ninds;i--;) { //!< reverse stacking order so indices are nicely sequential
      const size_t useind= cols-1-inds(i);
      lmask.push_back(useind+base);
    }
  }
  else {
    for (size_t i=0;i<ninds;i++)
      lmask.push_back(inds(i)+base);
  }
}

void addmask(LIST<size_t>& lmask, size_t cols, size_t base)
{
  lmask.reserve(lmask.size()+cols);
  for (size_t i=0;i<cols;i++)
    lmask.push_back(i+base);
}

static Warning<> fitmasknotrev_warning("-reversefrequency used with at least one data set but not with fit mask. This may mask out the wrong region!",&NMRsim_once_warning);
static Warning<> fitmaskrev_warning("-reversefrequency used with fit mask but not with any data set. This may mask out the wrong region!",&NMRsim_once_warning);

void parse_maskflags(bool& inv, bool& rev)
{
  enum { EXCLUDE=1, REVFREQ=2 };

  static flagsmap_type mask_flags;
  if (mask_flags.empty()) {
    mask_flags["exclude"]=EXCLUDE;
    mask_flags["reversefrequency"]=REVFREQ;
  }
  const int flags=parse_flags(mask_flags);
  inv = flags & EXCLUDE;
  rev = flags & REVFREQ;
  if (!nochecks) {
    if (used_reverse) {
      if (!rev)
	fitmasknotrev_warning.raise();
    }
    else {
      if (rev)
	fitmaskrev_warning.raise();
    }
  }
}

//#define FITMASKFLAGS "[-exclude|-reversefrequency]"
//const char* fitmasksyntax2D="<column mask>#[<row mask>]#" FITMASKFLAGS;
//const char* fitmasksyntaxrow="<mask set 1> [<mask set 2>...]#" FITMASKFLAGS;

FitCommand* FitMask::create()
{
  if (fit_sets.empty())
    error_abort("fit mask can only be used after data sets are defined");

  if (!are_left())
    return new FitMask();


  LIST<size_t> lmask;
  size_t base=0;
  bool inv,rev;

  if ((fit_sets.size()==1) && (fit_sets.front().rows()>1)) {
    const cmatrix& fit_set(fit_sets.front());
    LIST<size_t> sel1,rsel;
    LIST<bool> cind,rind;
    validatesel(sel1,cind,"Column selection (first argument)",1,fit_set.cols());
    if (parser_isnormal())
      validatesel(rsel,rind,"Row selection (second argument)",2,fit_set.rows());

    parse_maskflags(inv,rev);
    if (inv)
      inversemask(sel1,cind);
    const size_t cols=fit_set.cols();
    for (size_t i=0;i<fit_set.rows();i++) {
      if (rind(i))
	addmask(lmask,sel1,cols,base,rev);
      else {
	if (inv)
	  addmask(lmask,cols,base);
      }
      base+=cols;
    }
  }
  else { //1D
    LIST< LIST<size_t> > sels(fit_sets.size());
    LIST< LIST<bool> > inds(fit_sets.size());
    char scr[60];
    for (size_t i=0;i<fit_sets.size();i++) {
      sprintf(scr,"Selection for data set %" LCM_PRI_SIZE_T_MODIFIER "u",i+1);
      validatesel(sels(i),inds(i),scr,2,fit_sets(i).cols());
    }
    parse_maskflags(inv,rev);
    for (size_t i=0;i<fit_sets.size();i++) {
      LIST<size_t>& csel(sels(i));
      const size_t csize=fit_sets(i).size();
      if (inv)
	inversemask(csel,inds(i));
      addmask(lmask,csel,csize,base,rev);
      base+=csize;
    }
  }
  return new FitMask(lmask);
}
     
// void dumprange(bool& done, size_t start, size_t end, std::ostream& ostr)
// {
//   if (done)
//     ostr << ',';
//   else
//     done=true;
//   ostr << start+1;
//   if (start!=end)
//     ostr << '-' << end;
// }

// void outputcompressed(const BaseList<size_t>& a, std::ostream& ostr =std::cout)
// {
//   bool on=false;
//   bool done=false;
//   for (size_t i=0;i<a.size();i++) {
//     if (a(i)) {
//       if (!on) {
// 	on=true;
// 	start=i;
//       }
//     }
//     else {
//       if (on) {
// 	on=false;
// 	dumprange(done,start,i-1,ostr);
//       }
//     }
//   }
//   if (on)
//     dumprange(done,start,a.size()-1,ostr);
// }

smartptr<ProcessNormalise,false> FitNormalise::normaliseptr;

void FitNormalise::apply(DataStore& FID)
{
  if (!normaliseptr)
    return; //!< nothing to do
  normaliseptr->rawexec(FID);
//   if (FID.isnD())
//     normaliseptr->exec(FID.matrix().row(),states());
//   else
//     normaliseptr->exec(FID.listlist().row());
}

FitCommand* FitSystem::create()
{
  static char comline[MAXLINE];
  const char* restcom=get_curline();
  if (restcom && *restcom) {
    substitute_string(comline,sizeof(comline),restcom,SUB_NUMERIC);
    FitCommand* newcom=new FitSystem(comline);
    set_curline(NMRSIM_NULL);
    return newcom;
  }
  else
    throw Failed("missing command line");
}

void FitSystem::exec() const
{
  system(com_.c_str());
}

FitCommand* FitNormalise::create()
{
  if (isinteractive)
    std::cerr << "normalisation mode must be fixed before starting and cannot be changed interactively\n";
  else
    normaliseptr.reset(static_cast<ProcessNormalise*>(ProcessNormalise::create()));
  return NMRSIM_NULL;  //!< don't return an object
}

class FitSet : public FitCommand
{
public:
  FitSet() {}
  explicit FitSet(const BaseList<VarVariable*>& varsv)
    : vars_(varsv), haveval_(false) {}
  FitSet(const BaseList<VarVariable*>& varsv, double valv)
    : vars_(varsv), haveval_(true), val_(valv) {}
  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> vars_;
  bool haveval_;
  double val_;
};

class FitError : public FitCommand
{
public:
  explicit FitError(const BaseList<VarVariable*>& varsv, double errv =0.0) : vars_(varsv), err_(errv) {}
  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> vars_;
  double err_;
};

class FitFix : public FitCommand
{
public:
  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> varvs;  
};

class FitConstrain : public FitCommand
{
public:
  FitConstrain(const BaseList<VarVariable*>& varsv, const SimpleBoundsState& statev)
    : vars(varsv), state(statev) {}

  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> vars;
  SimpleBoundsState state;
};

class FitUnconstrain : public FitCommand
{
public:
  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> varvs;
};

class FitRelease : public FitCommand
{
public:
  void exec() const;
  static FitCommand* create();
private:
  LIST<VarVariable*> varvs;
};

#ifdef USE_MINUIT
struct FitMinuitRun : public FitCommand
{
  FitMinuitRun() { haverun=true; }
  void exec() const;
  static FitCommand* create() { return new FitMinuitRun(); }
  static FitCommand* create_interactive();
};
#endif

struct FitHelp : public FitCommand
{
  static FitCommand* create() { return NMRSIM_NULL; }
};

struct FitAbort : public FitCommand
{
  static FitCommand* create() { return NMRSIM_NULL; }
};

struct FitFinish : public FitCommand
{
  static FitCommand* create() { return NMRSIM_NULL; }
};

class FitTolerance : public FitCommand
{
public:
  explicit FitTolerance(double);
  FitTolerance() : tolerance_(-1.0) {}
  void exec() const;
  static FitCommand* create();
private:
  double tolerance_;
};

class FitIterations : public FitCommand
{
public:
  explicit FitIterations(size_t iters) : iters_(iters) {}
  FitIterations() : iters_(-1) {}
  void exec() const;
  static FitCommand* create();
private:
  int iters_;
};

class FitMethod : public FitCommand
{
public:
  explicit FitMethod(size_t =0);
  void exec() const;
  static FitCommand* create();
private:
  size_t flags_;
  static flagsmap_type flags;
};

class FitPrecision : public FitCommand
{
public:
  explicit FitPrecision(double epsv);
  FitPrecision() : eps_(0.0) {}
  void exec() const;
  static FitCommand* create();
 private:
  double eps_;
};

class FitStatistics : public FitCommand
{
public:
  explicit FitStatistics(int flagsv) : flags_(flagsv) {}
  void exec() const;
  static FitCommand* create();
private:
  int flags_;
};

flagsmap_type FitMethod::flags;

class FitRun : public FitCommand
{
public:
  FitRun(int flagsv, bool isinterv =false) : flags(flagsv), isinter(isinterv) {
    haverun=true;
  }
  enum { allowinteractive=1 };
  void exec() const;
  static FitCommand* create();
  static FitCommand* create_interactive();
private:
  int flags;
  bool isinter;
};

const FitCommand_t finish_desc(&FitFinish::create,"finish\t - terminate interactive session");
const FitCommand_t abort_desc(&FitAbort::create, "abort\t - abort optimisation");
const FitCommand_t help_desc(&FitHelp::create, "help\t - display this help message");
const FitCommand_t mask_desc(&FitMask::create, "mask [<indices set 1> <indices set 2>...|<column indices> [<row indices>]] [-exclude] - set data points to include in fitting (-exclude inverts selection).  Without arguments, remove mask"); 
const FitCommand_t release_desc(&FitRelease::create, "release [<parameter>*]\t - release previously fixed parameter");
const FitCommand_t constrain_desc(&FitConstrain::create, "constrain <parameter> [<min> -minimum|<max> -maximum|<min> <max>] - constrain minimum and/or maximum value of a parameter");
const FitCommand_t unconstrain_desc(&FitUnconstrain::create, "unconstrain [<parameter>*]\t - release constraints on all or specified parameters");
const FitCommand_t method_desc(&FitMethod::create, "method [-simplex|-gradient] [-randomise] [-real|-imag|-complex] [-no_eta_constraints] set optimisation method and/or flags");
const FitCommand_t fix_desc(&FitFix::create, "fix <parameter>\t - fix parameter value");
const FitCommand_t set_desc(&FitSet::create, "set [<parameter> [<value>]] - set (or display) parameter value ('set' displays all parameters)");
const FitCommand_t normalise_desc(&FitNormalise::create, "normalise|normalize [<to>] [-integral|-minmax|-abs]");
const FitCommand_t error_desc(&FitError::create, "error <parameter> [<value>]]\t - set (or display) estimated error on parameter value");
const FitCommand_t system_desc(&FitSystem::create, "system <command> \t - execute system command");
const FitCommand_t precision_desc(&FitPrecision::create, "precision [<machine eps>]\t - set (or display) machine precision");
const FitCommand_t tolerance_desc(&FitTolerance::create, "tolerance [<tolerance>]\t - set (or display) convergence tolerance");
const FitCommand_t iterations_desc(&FitIterations::create, "iterations [<iteractions/evaluations>] - set (or display) maximum number of interations/function evaluations");
const FitCommand_t statistics_desc(&FitStatistics::create, "statistics [-release_constraints|-silent] - calculate 'fitting' statistics");
const FitCommand_t run_desc(&FitRun::create, "run\t - run fitting");
const FitCommand_t interactive_desc(&FitRun::create_interactive);

varvarpars_t var_varpars;

static int opt_flags=OPT_COMPLEX;
static size_t opt_method=OPT_GRADIENT;
static int maxiter=-1;
static double opt_tolerance=0.0;
static optim_t opttype=OPTIM_NONE;
static double fit_noise=0.0;

//optim_t optimisetype(); //!< return optimisation mode
bool perform_optim(MasterObj&, DataStore&);
bool perform_fit(MasterObj&, DataStore&);

// void VarVariable::error(double err)
// {
//   if (err<=0.0)
//     throw InvalidParameter("VarVariable::error");
//   step=err;
// }

void VarVariable::printname(std::ostream& ostr) const
{
  if (isnamed())
    ostr << "- " << name();
  else {
    ostr << '(';
    Variable::printname(ostr);
    ostr << ')';
  }
}

typedef MAPTYPE(VarVariable*) namedvarmap_t;
namedvarmap_t namedvarmap; //!< list of named fitting variables

//const double default_step_scale=0.1; //default error is 10%

rmatrix covar;
static LIST<processing_state> fit_procs;
static LIST<size_t> row_sel,col_sel;

ThreadWarning<> impossiblesave_warning("-statistics/-original save requested but no fitting data available (fitting not active or save erroneously placed in proc rather than finalise?)",&NMRsim_once_warning);

size_t sizefactor()
{
  return ((opt_flags & complexity_types)==OPT_COMPLEX) ? 2U : 1U;
}

namespace {

  inline size_t fitnpts(size_t start) { return sizefactor()*(mask.empty() ? start : mask.size()); }

  bool fitsaveok() {
    if (!covar && residuals.empty()) {
      impossiblesave_warning.raise();
      return false;
    }
    return true;
  }

  void extractdata(BaseList<double> dest, const BaseList<complex>& source)
  {
    static LIST<complex> tsource;
    const BaseList<complex>* sourcep=&source;
    if (!(mask.empty())) {
      tsource=source(mask);
      sourcep=&tsource;
    }
    switch (opt_flags & complexity_types) {
    case OPT_COMPLEX:
      dest=asdoubles(*sourcep);
      break;
    case OPT_REAL:
      dest=reals(*sourcep);
      break;
    case OPT_IMAG:
      dest=imags(*sourcep);
      break;
    default:
      throw InternalError("extractdata");
    }
  }

  void create_residuals(LIST<double>& dest)
  {
    dest.create(fitnpts(residuals.size()));
    extractdata(dest,residuals.row());
  }
}

Fit_Factory_t& get_fit_Factory()
{
  static Fit_Factory_t fit_Factory;
  return fit_Factory;
}

Fit_Factory_t& get_minmax_Factory()
{
  static Fit_Factory_t minmax_Factory;
  return minmax_Factory;
}

GlobalFit_Factory_t& get_setup_fit_Factory()
{
  static GlobalFit_Factory_t setup_fit_Factory;
  return setup_fit_Factory;
}

void dump_factory(const Fit_Factory_t& factory)
{
  const Fit_Factory_t::const_iterator end(factory.end());
  Fit_Factory_t::const_iterator start(factory.begin());
  while (start!=end) {
    const char* desc= (start->second).description;
    if (desc) //!< NMRSIM_NULL indicates command that should not be used interactively
      std::cout << desc << '\n';
    ++start;
  }
}

void save_statistics(const ProcessSave& saveobj, const ProcessSave_state& cstate, int)
{
  if (fitsaveok()) {
    if (!covar.empty())
      saveobj.save(covar,cstate,"covariance");
    if (!residuals.empty()) {
      LIST<double> aresiduals;
      create_residuals(aresiduals);
      saveobj.save(aresiduals,cstate,"residuals");
    }
  }
}

void applymask(DataStore& dest, const DataStore& source, const BaseList<size_t>& mask, bool flip)
{
  if (source.empty())
    throw InternalError("applymask");
  dest.duplicate_structure_processing(source);
  BaseList<complex> destrow(dest.row());
  const BaseList<complex>& sourcerow(source.row());
  destrow=complex(0.0);
  const size_t Nm1=source.size()-1;
  for (size_t i=mask.size();i--;) {
    const size_t ind= flip ? Nm1-mask(i) : mask(i);
    destrow(ind)=sourcerow(ind);
  }
}

void save_original(const ProcessSave& saveobj, const ProcessSave_state& cstate, int)
{
  if (fitsaveok()) {
    // ought to copy out useful information on sw etc. 
    if (mask.empty())
      saveobj.save(fit_set,cstate,"original");
    else {
      const bool isflipped=saveobj.original_isflipped();
      const bool needflip=saveobj.flags() & ProcessSave::REVERSEFREQ;
      if (isflipped!=needflip)
	throw InternalError("original reserve state doesn't match reversefreq flag");
      DataStore tmpfit_set;
      applymask(tmpfit_set,fit_set,mask,isflipped);
      saveobj.save(tmpfit_set,cstate,"original");
    }
  }
}

double final_optimisation=0.0;
double final_chisquared=0.0;
static SystemVariable<double*> v_final_optimisation("final_optimisation",&final_optimisation,1.0,V_ISFIXED);
static SystemVariable<double*> v_final_chisquared("final_chisquared",&final_chisquared,1.0,V_ISFIXED);

prefinalise_callback_t old_prefinalise_callback=NMRSIM_NULL;

bool parse_optimise_block()
{
  return read_block("optimise",optimiseblock_Factory,true);
}

void create_finalisevars()
{
  if (old_prefinalise_callback)
    (*old_prefinalise_callback)();

  add_systemvarmap(v_final_optimisation);
  add_systemvarmap(v_final_chisquared);
}

ContextWarning<> callback_override_warning("over-riding previously modified calculation callback",&NMRsim_repeat_warning);

void replace_callback()
{  
  if (calculation_callback)
    callback_override_warning.raise();
  else {
    if ((verbose & VER_OPTIM) && (verbose_level>1))
      std::cout << "Setting optimisation callback\n";
  }
  switch (opttype) {
  case OPTIM_FIT:
    calculation_callback=perform_fit; break;
  case OPTIM_OPT:
    calculation_callback=perform_optim; break;
  default:
    throw InternalError("replace_callback");
  }

  old_prefinalise_callback=prefinalise_callback;
  prefinalise_callback=create_finalisevars;
}

class ProcessStatistics : public ProcessCommand {
public:
  explicit ProcessStatistics(int flagsv)
    : ProcessCommand(PROC_HAS2D), flags_(flagsv) {}
  void print(std::ostream& ostr) const;
  void rawexec(DataStore&) const;
  static ProcessCommand* create();
private:
  int flags_;
};

void parse_fit();
void parse_minmax(int ismax);

struct Fit_Proxy_ {
  Fit_Proxy_() {
    command_Factory_t& par_Factory(get_par_Factory());
    par_Factory["fit"]=par_t(&parse_fit,true);
    par_Factory["minimise"]=par_t(&parse_minmax,0,true);
    par_Factory["maximise"]=par_t(&parse_minmax,1,true);

    optimiseblock_Factory["fit"]=par_t(&parse_fit,true);
    optimiseblock_Factory["minimise"]=par_t(&parse_minmax,0,true);
    optimiseblock_Factory["maximise"]=par_t(&parse_minmax,1,true);

    Process_Factory_t& finalise_Factory(get_Finalise_Factory());
    finalise_Factory["fit"]=&ProcessStatistics::create;

    register_save_function("statistics",save_statistics);
    register_save_function("original",save_original);
  }
};

//declare Fit additions
static const Fit_Proxy_ fit_proxy_;

static LIST<FitCommand*> fitstack;

//ContextWarning<> fitdirective_ignored_warning("optimisation/fitting directive ignored: ",&NMRsim_repeat_warning);

void make_optim_command(const char* keyname, const Fit_Factory_t& factory)
{
  const Fit_Factory_t::mapped_type funcdesc(factory_parse(keyname,factory));

  FitCommand* fitcom=(*(funcdesc.funcp))();
  if (fitcom) {
    parser_checkfinished();
    fitstack.push_back(fitcom);
  }
  // else
  //  fitdirective_ignored_warning.raise(keyname);
}

void runinteractive()
{
  const bool isfit=(opttype==OPTIM_FIT);
  const char* intro = isfit ? "fit> " : "optimise> ";
  const Fit_Factory_t& factory(isfit ? get_fit_Factory() : get_minmax_Factory());
  for (;;) {
    parser_getinteractive(intro);
    const char* keyname=parse_string(F_REPLACEDOLLAR);
    const Fit_Factory_t::const_iterator curdesc=factory.find(keyname);
    if (curdesc==factory.end()) {
      std::cerr << "Unknown or misplaced directive: " << keyname << " (try 'help'?)\n";
      continue;
    }
    const Fit_Factory_t::mapped_type funcdesc(curdesc->second);
    FitCommand* fitcom=NMRSIM_NULL;
    bool parsefailed=false;
    try {
      fitcom=(*(funcdesc.funcp))();
      if (fitcom && are_left())
	parsefailed=true;
    } catch (const std::exception&) {
      parsefailed=true;
    }
    if (parsefailed) {
      std::cerr << "Syntax: " << funcdesc.description << '\n';
      if (fitcom)
	delete fitcom;
      continue;
    }
    if (fitcom) {
      try {
	fitcom->exec();
      } catch (const std::exception& exc) {
	std::cerr << "Optimisation directive failed: " << exc.what() << '\n';
      }
      delete fitcom;      
    }
    else {
      if (funcdesc.funcp==FitAbort::create) {
	haveresults=false;
	return;
      }
      if (funcdesc.funcp==FitFinish::create)
	return;
      if (funcdesc.funcp==FitHelp::create) {
	std::cout << "Available commands:\n";
	dump_factory(factory);
	std::cout << '\n';
      }
      else
	throw InternalError("Unrecognised optimisation directive");
    }
  }
}

Optimise_Factory_t& get_Optimise_Factory()
{
  static Optimise_Factory_t optimise_Factory;
  return optimise_Factory;
}

//static LIST<size_t> actord; //!< actual parameter index

double get_tolerance()
{
  if (opt_tolerance)
    return opt_tolerance;
  return (opttype==OPTIM_FIT) ? NMRSIM_DEFAULT_FIT_TOLERANCE : NMRSIM_DEFAULT_OPTIM_TOLERANCE;
}

size_t get_maxiter()
{
  if (maxiter>=0)
    return maxiter;

  if (opttype==OPTIM_FIT)
    return NMRSIM_DEFAULT_ITERATIONS;

  const size_t npar=countactive();
  return 200+100*npar+5*npar*npar; //!< Minuit's default MAXFCN
}
 
ContextWarning<> nonterminalstar_warning("* can only be used for parameter matching at end of name e.g. t1* not *t1",&NMRsim_once_warning);

template<typename T> void findpartialmap(LIST<T>& dest, const MAPTYPE(T)& map, const char* name, const char* mapname, bool dumpoptions =true)
{
  const size_t nchars=strlen(name)-1;
  const typename MAPTYPE(T)::const_iterator end(map.end());
  const size_t oldlen=dest.size();
  
  if (name[nchars]!='*') {
    const typename MAPTYPE(T)::const_iterator curp(map.find(name));
    if (curp!=end) {
      dest.push_back(curp->second);
      return;
    }
  }
  else {
    typename MAPTYPE(T)::const_iterator start(map.begin());
    while (start!=end) {
      const std::string& key(start->first);
      if (key.compare(0,nchars,name,nchars)==0)
	dest.push_back(start->second);
      ++start;
    }
  }
  if (dest.size()==oldlen) {
    dumpmap(map,name,mapname,dumpoptions);
    if (!nochecks) {
      const char* starptr=strchr(name,'*');
      if (starptr && ((starptr-name)!=nchars))
	nonterminalstar_warning.raise();
    }
    error_abort();
  }
}

void getpar(LIST<VarVariable*>& destlist)
{
  ensure_named();
  const char* dest(parse_string(F_REPLACEDOLLAR));
  char* tail;  
  long n=strtol(dest,&tail,10);
  if (*tail=='\0') {
    if ((n<1) || (n>var_varpars.size())) {
      parser_printcontext() << "parameter number must be between 1 and " << var_varpars.size() << std::endl;
      if (isinteractive)
	throw BadIndex("optimisation parameter",n,var_varpars.size()); //!< need to throw so can be caught in interactive mode
      error_abort();
    }
    destlist.push_back((var_varpars(n-1)));
  }
  else
    findpartialmap(destlist,namedvarmap,dest,"fitting parameter");
}

void getpars(LIST<VarVariable*>& dest)
{
  dest.create(0U);
  while (are_left())
    getpar(dest);
}

FitCommand* FitRelease::create()
{
  FitRelease* newp=new FitRelease();
  getpars(newp->varvs);
  return newp;
}

void FitRelease::exec() const 
{
  if (varvs.empty()) {
    for (size_t i=var_varpars.size();i--;)
      var_varpars(i)->release();
  }
  else {
    for (size_t i=varvs.size();i--;)
      varvs(i)->release();
  }
}

void VarVariable::constrain(const SimpleBoundsState& state)
{
  parameter.constrain(*(state.function()));
  boundstate=state;
}

void VarVariable::unconstrain()
{
  if (boundstate.unconstrain())
    parameter.unconstrain();
}

FitCommand* FitUnconstrain::create()
{
  FitUnconstrain* newp=new FitUnconstrain();
  getpars(newp->varvs);
  return newp;
}

void FitUnconstrain::exec() const 
{
  if (varvs.empty()) {
    for (size_t i=var_varpars.size();i--;)
      var_varpars(i)->unconstrain();
  }
  else {
    for (size_t i=varvs.size();i--;)
      varvs(i)->unconstrain();
  }
}

ThreadWarning<> fit_allfixed_warning("all parameters are fixed.  One or more parameters must be released before running optimisation",&NMRsim_repeat_warning);

void FitFix::exec() const
{
  for (size_t i=varvs.size();i--;)
    varvs(i)->fix();
  if (countactive()==0)
    fit_allfixed_warning.raise();
}

ContextWarning<> fit_noparameters_warning("fix with no parameters specified (does nothing)",&NMRsim_once_warning);

FitCommand* FitFix::create()
{
  FitFix* newp=new FitFix();
  getpars(newp->varvs);
  if (newp->varvs.empty())
    fit_noparameters_warning.raise();
  return newp;
}

namespace {

  flagsmap_type stat_flags;

  int parse_stat_flags() {
    if (stat_flags.empty()) {
      stat_flags["release_constraints"]=STAT_RELEASE;
      stat_flags["silent"]=STAT_SILENT;
      stat_flags["if_constrained"]=STAT_IFCONS;
    }
    return parse_flags(stat_flags);
  }
}

void ProcessStatistics::print(std::ostream& ostr) const
{
  ostr << "fit statistics ";
  dumpflags(ostr,stat_flags,flags_);
}

FitCommand* FitStatistics::create()
{
  return new FitStatistics(parse_stat_flags());
}

ProcessCommand* ProcessStatistics::create()
{
  const char* keyname=parse_string(F_REPLACEDOLLAR);
  if (strcmp(keyname,"statistics")!=0)
    error_abort("only fit statistics is valid in finalise block");
  return new ProcessStatistics(parse_stat_flags());
}

ThreadWarning<> fit_nonoiselevel_warning("Used default noise level of 1.0 i.e. errors should be scaled by true noise level",&NMRsim_repeat_warning);
ThreadWarning<> fit_randomassumed_warning("Errors are estimated assuming residual is random noise i.e. overestimated if fit is poor.  Use fit noiselevel to specify noise level explicitly to avoid this.",&NMRsim_repeat_warning);

bool havenoisespec() { return !(noisevals.empty()) || (fit_noise!=0.0); }

void ensure_noisevector(size_t npuse)
{
  if (!(noisevector.empty()))
    return;
  if (npuse==0)
    throw InternalError("ensure_noisevector");
  const size_t sfactor=sizefactor();
  noisevector.create(fitnpts(npuse));
  //  bool ok=true;
  const int fitverb=fit_verbose_level();
  if (noisevals.empty()) {
    //if (!fit_noise)
    //  ok=false;
    noisevector=fit_noise ? fit_noise : 1.0;
  }
  else {
    if (noisevals.size()!=global_nps.size())
      error_abort("Noise specification doesn't match data set");
    noisevector.create(0U);     
    for (size_t i=0;i<noisevals.size();i++) {
      if (noisevals(i)==0.0)
	error_abort("Zero noise level detected");
      noisevector.push_back(global_nps(i)*sfactor,noisevals(i));
    }
  }
  if (fitverb>1)
    std::cout << "Noise vector: " << noisevector << '\n';
}

void docalc_statistics(int flags)
{
  if (flags & STAT_IFCONS) {
    bool havecons=false;
    for (size_t i=var_varpars.size();i--;) {
      const VarVariable& cvar(*(var_varpars(i)));
      if (!(cvar.isfixed()) && cvar.isconstrained()) {
	havecons=true;
	break;
      }
    }
    if (!havecons) //!< no active constraints so do nothing
      return;
  }

  const bool ignoreconstraints=flags & (STAT_RELEASE | STAT_IFCONS);
  prepare_paras(true,true);
  if (fitobjp==NMRSIM_NULL)
    throw InternalError("docalc_statistics");
  ensure_noisevector(global_np);
  const int fitverb=fit_verbose_level();
  covariance(covar,*fitobjp,parameter_list,noisevector,fitverb,ignoreconstraints);
  if (!havenoisespec() && global_chisqr)
    covar*=global_chisqr;
  calc_parameters(valerrs);
  if (covar.rows()==1)
    covar.clear(); //!< don't bother with 1x1

  updatevarerrors();
  getvarnames();
  if (!silent && !(flags & STAT_SILENT))
    display_statistics(ignoreconstraints);
}

void FitStatistics::exec() const
{
  docalc_statistics(flags_);
}

void ProcessStatistics::rawexec(DataStore& FID) const
{
  //  build_noisevector();
  if (masterobjp==NMRSIM_NULL)
    throw InternalError("ProcessStatistics::rawexec");
  build_fitobjs(*masterobjp);
  ensure_named();
  docalc_statistics(flags_);
}

void FitConstrain::exec() const
{
  for (size_t i=vars.size();i--;)
    vars(i)->constrain(state);
}

FitCommand* FitConstrain::create()
{
  LIST<VarVariable*> vars;
  getpar(vars);
  //  VarVariable& var(getpar());
  const double val=parse_double();
  double val2=0.0;
  int type=-1;
  static flagsmap_type flags;
  if (parser_isnormal())
    val2=parse_double();
  else {
    if (flags.empty()) {
      flags["maximum"]=SimpleBoundsState::MAX;
      flags["minimum"]=SimpleBoundsState::MIN;
    }
    type=parse_flags(flags);
  }
  SimpleBoundsState state;
  switch (type) {
  case -1:
    state.lowerupper(val,val2);
    break;
  case SimpleBoundsState::MIN:
    state.lower(val);
    break;
  case SimpleBoundsState::MAX:
    state.upper(val);
    break;
  default:
    error_abort("Syntax: constrain <parameter> <min> <max>|<min> -minimum|<max> -maximum");
  }
  return new FitConstrain(vars,state);
}
    
double VarVariable::parse_error()
{
  const char* str(parse_string(F_REPLACEDOLLAR));
  char* tail;
  double err=strtod(str,&tail);  
  switch (*tail) {
  case '\0':
    break;
  case '%':
    if (get()==0.0) {
      parser_printcontext() << "can't use % scale factor if parameter is zero\n";
      throw Failed("parse_error"); //!< don't error abort (so interactive can catch)
    }
    err*=get()*0.01;
    break;
  default:
    parser_printcontext() << "failed to parse as error estimate: " << str << '\n';
    throw Failed("parse_error"); //!< don't error abort (so interactive can catch)
  }
  return err;    
}

FitCommand* FitError::create()
{ 
  LIST<VarVariable*> vars;
  getpar(vars);
  const double error=are_left() ? vars.front()->parse_error() : 0.0;
  return new FitError(vars,error);
}

void FitError::exec() const
{
  for (size_t i=0;i<vars_.size();i++) {
    VarVariable& var_(*(vars_(i)));
    if (err_)
      var_.error(err_);
    else
      std::cout << "Estimated error on parameter " 
		<< var_.name()
		<< ": " << var_.error()
		<< " (" << 100.0*(var_.error()/var_.get()) << "%)\n";
  }
}

FitCommand* FitSet::create()
{ 
  if (!are_left())
    return new FitSet();
  LIST<VarVariable*> varsv;
  getpar(varsv);
  //  VarVariable& varv(getpar());
  return are_left() ? new FitSet(varsv,parse_double()) : new FitSet(varsv);
}

void FitSet::exec() const
{
  if (vars_.empty()) {
    prepare_paras(true); //!< might need to update bounds
    dump_parameters();
    return;
  }
  for (size_t i=0;i<vars_.size();i++) {
    VarVariable* varp_(vars_(i));
    if (haveval_)
      varp_->set(val_);
    else
      std::cout << *varp_ << '\n';
  }
}

FitTolerance::FitTolerance(double v)
{
  if (v<=0.0 || v>=1.0)
    error_abort("tolerance criterion must be <1 and >0");
  tolerance_=v;
}

FitPrecision::FitPrecision(double v)
{
  if (v<=0.0 || v>=1e-3)
    error_abort("machine precision must be <1e-3 and >0");
  eps_=v;
}

void FitTolerance::exec() const
{
  if (tolerance_<0.0) {
    if (!silent) {
      std::cout << "Tolerance criterion: " << opt_tolerance;
      if (opt_tolerance==0.0)
	std::cout << " (" << get_tolerance() << ')';
      std::cout << '\n';
    }
  }
  else {
    opt_tolerance=tolerance_;
    if (verbose & VER_OPTIM)
      std::cout << "Tolerance criterion set to " << get_tolerance() << '\n';
  }
}

void FitPrecision::exec() const
{
  if (eps_==0.0) {
    if (!silent) {
      std::cout << "machine eps: ";
      if (minuit_precision)
	std::cout << minuit_precision << '\n';
      else
	std::cout << "<Minuit default>\n";
    }
  }
  else {
#if MINUITVER>1
    minuit_precision=eps_;
    if (verbose & VER_OPTIM)
      std::cout << "Minuit machine precision set to " << eps_ << '\n';
#endif
  }
}

void FitIterations::exec() const 
{
  if (iters_<0) {
    if (!silent) {
      std::cout << "Maximum iterations: evaluations: " << maxiter;
      if (maxiter<0)
	std::cout << " (" << get_maxiter() << ')';
      std::cout << '\n';
    }
  }
  else {
    maxiter=iters_; 
    if (verbose & VER_OPTIM)
      std::cout << "Maximum iterations/evaluations set to " << get_maxiter() << '\n';
  }
}

FitCommand* FitRun::create()
{
  static flagsmap_type flags;
  if (flags.empty())
    flags["allowinteractive"]=allowinteractive;
  return new FitRun(parse_flags(flags));
}

FitCommand* FitRun::create_interactive()
{
  if (silent)
    error_abort("can't optimise interactively with -silent");
  if (global_workers)
    error_abort("can't optimise interactively in multi-processor mode");
  isinteractive=true;  
  return new FitRun(0,true);
}

FitMethod::FitMethod(size_t flagsv)
{
  if (flagsv) {
    switch (flagsv & opt_types) {
    case OPT_SIMPLEX: case OPT_GRADIENT:
      break;
    case 0:
      flagsv|=NMRSIM_DEFAULT_OPTIM_METHOD;
      break;
    default:
      error_abort("can only set one of -simplex / -gradient");
    }
    switch (flagsv & complexity_types) {
    case OPT_COMPLEX: case OPT_IMAG: case OPT_REAL:
      break;
    case 0:
      flagsv|=OPT_COMPLEX;
      break;
    default:
      error_abort("can only set one of -real / -imag / -complex");
    }
  }
  flags_=flagsv;
}

void FitMethod::exec() const
{
  if (flags_) {
    opt_method=(flags_ & opt_types);
    opt_flags=flags_ & ~opt_types;
  }
  else {
    std::cout << "fit method ";
    dumpflags(std::cout,flags,(opt_flags | opt_method)) << '\n';
  }
}

void parse_fit_exp();
void parse_fit_add();
void parse_fit_noiselevel();

ThreadWarning<> fit_directiveafterrun_warning("fitting/optimisation instructions after last 'run' will have no effect",&NMRsim_once_warning);

void execfitstack()
{
  ensure_named();
  if (!haverun) {
#ifdef USE_MINUIT
    if (opttype!=OPTIM_FIT)
      fitstack.push_back(new FitMinuitRun());
    else
#endif
      fitstack.push_back(new FitRun(FitRun::allowinteractive));
  }
  else {
    if (dynamic_cast<FitRun*>(fitstack.back())==NMRSIM_NULL)
      fit_directiveafterrun_warning.raise();
  }

  for (size_t i=0;i<fitstack.size();i++)
    fitstack(i)->exec();
}
 
ThreadWarning<> fit_sw_warning("apparently mismatched spectral parameters",&NMRsim_repeat_warning);
   
namespace {
  //! set variables from raw values
  void set_variables(const BaseList<double>& pars)
  {
    if (pars.size()!=parameter_list.size())
      throw Mismatch("set_variables");
    for (size_t i=pars.size();i--;) {
      if (!(parameter_list(i).isfixed()))
	var_varpars(i)->set(pars(i));
    }
    update_expressions();
  }

  //! set variables from parameter list
  void set_variables()
  {
    for (size_t i=parameter_list.size();i--;) {
      const Parameter& lpar(parameter_list(i));
      if (!lpar.isfixed()) {
	VarVariable& lvar(*(var_varpars(i)));
	lvar.set(lpar.get());
      }
    }
    update_expressions();
  }

  void rawcheck(processing_state& datastate, const processing_state& calcstate, const char* name, double processing_state::* ptr, bool fillempty =false)
  {
      const double v(calcstate.*ptr);
      double& data_v(datastate.*ptr);
      if (fillempty && (data_v==0.0))
	data_v=v;
      else {
	const double maxv=std::max(v,data_v);
	if (fabs(v-data_v)>1e-5*maxv) {
	  char buf[256];
	  snprintf(buf,sizeof(buf),": %s of fitting data=%g kHz, %s of calculated data=%g kHz",name,data_v*1e-3,name,v*1e-3);
	  fit_sw_warning.raise(buf);
	}
      }
  }

  void rawcheck(BaseList<processing_state> datastates, const BaseList<processing_state>& calcstates, const char* name, double processing_state::* ptr, bool fillempty =false)
  {
    for (size_t i=0;i<datastates.size();i++)
      rawcheck(datastates(i),calcstates(i),name,ptr,fillempty);
  }

  void checksw(DataStore& data, const DataStore& calc)
  {
    BaseList<processing_state> datastates(data.states());
    const BaseList<processing_state>& calcstates(calc.states());

    if (active2D) {
      if (datastates.size()!=calcstates.size()) {
	std::cerr << "Number of dimensions in fitting and calculated data sets are different (" << datastates.size() << " vs. " << calcstates.size() << ")\n";
	error_abort();
      }
      
      rawcheck(datastates,calcstates,"sw",&processing_state::sw,true);
      rawcheck(datastates,calcstates,"ref",&processing_state::ref);
      //!< don't bother comparing sfrq
    }
    else {
      // if row-oriented, just compare last
      rawcheck(datastates.back(),calcstates.back(),"sw",&processing_state::sw,true);
      rawcheck(datastates.back(),calcstates.back(),"ref",&processing_state::ref);
    }
  }

#ifdef USE_MINUIT
  void set_variables(const LCM_MINUITNAMESPACE::MnUserParameterState& pars)
  {
    for (size_t i=var_varpars.size();i--;)
      var_varpars(i)->set(
#if MINUITVER==2
			  pars.Value(i)
#else
			  pars.value(i)
#endif
);
    update_expressions();
  }
#endif
}

const char* VarVariable::ensure_named()
{
  std::ostringstream str(std::ostringstream::out);
  ptr->printvariablename(str,subsid); 
  if (variable().isarray())
    str << '_' << (arraywhich+1);
  namestr=rawname()=str.str(); //!< use same name in both RealVariable and Parameter
  return name();
}
  
void ensure_named()
{
  if (dirtystart<0)
    return;
  const bool isverb=(verbose & VER_OPTIM) && (verbose_level>1);
  for (int i=dirtystart;i<var_varpars.size();i++) {
    VarVariable& varv(*(var_varpars(i)));
    varv.validate();
    const char* name(varv.ensure_named());
    if (*name=='\0') {
      if (isverb)
	std::cout << "Failed to name fitting/optimisation variable " << (i+1) << '\n';
      continue;
    }
    namedvarmap[name]=&varv;
    
    if (isverb)
      std::cout << "Naming fitting/optimisation variable " << (i+1) << ": " << varv.name() << '\n';
  }
  dirtystart=-1;
}
 
void register_fitting_variable(VarVariable* vp)
{
  if (dirtystart<0)
    dirtystart=var_varpars.size();
  if ((verbose & VER_OPTIM) && (verbose_level>1))
    std::cout << "Creating (array) fitting/optimisation variable " << (var_varpars.size()+1) << '\n';
  var_varpars.push_back(vp);
  dirty_stack_push(vp);
}

void FitObj::operator()(BaseList<double> dest, const BaseList<double>& pars) const
{
  if ((verbose & VER_OPTIM) && (verbose_level>1))
    std::cout << "Fitting variables: " << pars << '\n';
  set_variables(pars);
  obj.calc(tmp,0);
  FitNormalise::apply(tmp); //!< apply any normalisation
  if (!(fit_set.empty()))
    ensure_dataset(tmp);
  extractdata(dest,tmp.row());
}

SystemVariable<double*> v_fit_noise("noiselevel",&fit_noise);

FitCommand* FitMethod::create()
{
  if (flags.empty()) {
    flags["randomise"]=OPT_RANDOMISE;
    flags["simplex"]=OPT_SIMPLEX;
    flags["gradient"]=OPT_GRADIENT;
    flags["imag"]=OPT_IMAG;
    flags["real"]=OPT_REAL;
    flags["complex"]=OPT_COMPLEX;
    flags["noconstraints"]=OPT_NOCONSTRAINTS;
    flags["no_eta_constraints"]=OPT_NOCONSTRAINTS;
  }
  return new FitMethod(parse_flags(flags));
}

FitCommand* FitIterations::create()
{
  return are_left() ? new FitIterations(parse_unsigned()) : new FitIterations();
}
  
FitCommand* FitTolerance::create()
{
  return are_left() ? new FitTolerance(parse_double()) : new FitTolerance();
}

ContextWarning<> precwarning("SetPrecision not supported by this version of Minuit (directive ignored)",&NMRsim_once_warning);

FitCommand* FitPrecision::create()
{
#if MINUITVER>1
  return are_left() ? new FitPrecision(parse_double()) : new FitPrecision();
#else
  precwarning.raise();
  set_curline(NMRSIM_NULL);
  return NMRSIM_NULL;
#endif
}

void checknotbuilt()
{
  static ThreadWarning<> alreadybuilt_warning("fit data set has already been referenced.  Further changes to the data set are probably an error.",&NMRsim_once_warning);
  if (!fit_set.empty())
    alreadybuilt_warning.raise();
}

void read_file(cmatrix& tmp_set, filestruct& info, const char* fname, LIST<std::string>* varnamesp =NMRSIM_NULL)
{
  checknotbuilt();
  raw_read_file(tmp_set,info,fname,varnamesp);
}

size_t add_data_set(cmatrix& tmp_set, const processing_state& tmp_proc)
{
  fit_sets.push_back(cmatrix());
  fit_sets.back().swap(tmp_set);
  fit_procs.push_back(tmp_proc);  
  return fit_procs.size();
}

static void ensure_fit()
{
  switch (opttype) {
  case OPTIM_NONE:
    opttype=OPTIM_FIT;
    replace_callback();
  case OPTIM_FIT:
    break;
  default:
    throw InternalError("ensure_fit");
  }
}

std::pair<domain_t,bool> getdomain()
{
  static flagsmap_type flags;
  enum { LD_FREQUENCY=1, LD_TIME=2, LD_REVFREQ=4 };
  if (flags.empty()) {
    flags["time"]=LD_TIME;
    flags["frequency"]=LD_FREQUENCY;
    flags["reversefrequency"]=LD_REVFREQ;
  }

  const int dflags=parse_flags(flags,0,true);
  domain_t domain=D_UNKNOWN;
  bool rev=false;
  switch (dflags) {
  case LD_FREQUENCY:
    domain=D_FREQUENCY;
    break;
  case LD_TIME:
    domain=D_TIME;
    break;
  case LD_REVFREQ:
    domain=D_FREQUENCY;
    rev=true;
    break;
  }
  return std::pair<domain_t,bool>(domain,rev);
}

static LIST<size_t> npstore;
static LIST<double> swstore;
//static LIST<size_t> nistore;

void makevars(size_t setno, size_t np, size_t ni, const filestruct& info)
{
  char scratch[64];
  sprintf(scratch,"sw_exp%" LCM_PRI_SIZE_T_MODIFIER "u",setno);
  (void)UserVariable::create(scratch,info.sw,UserVariable::IGNORE_UNUSED);
  sprintf(scratch,"np_exp%" LCM_PRI_SIZE_T_MODIFIER "u",setno);
  (void)UserVariable::create(scratch,double(np),UserVariable::IGNORE_UNUSED);
  sprintf(scratch,"ni_exp%" LCM_PRI_SIZE_T_MODIFIER "u",setno);
  (void)UserVariable::create(scratch,double(ni),UserVariable::IGNORE_UNUSED);
  sprintf(scratch,"ref_exp%" LCM_PRI_SIZE_T_MODIFIER "u",setno);
  (void)UserVariable::create(scratch,info.ref,UserVariable::IGNORE_UNUSED);

  const bool inseq=((setno-1)==swstore.size());
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "sw/np in sequence: " << (inseq ? "Yes\n" : "No\n");
  
  if (inseq) {
    swstore.push_back(info.sw);
    npstore.push_back(np);
    //    nistore.push_back(ni);
  }
  else {
    swstore.clear();
    npstore.clear();
    //    nistore.clear();
  }
}

ThreadWarning<> autoset_notimp_warning("Setting np/sw from experimental data sets not implemented for muitple data sets",&NMRsim_once_warning);

bool try_find_sw()
{
  switch (swstore.size()) {
  case 1:
    if (verbose & VER_GEN)
      std::cout << "Setting sw from experimental data\n";
    v_sw=swstore.front();
    return true;
  case 0:
    break;
  default:
    autoset_notimp_warning.raise();
  }
  return false;
}
    
bool try_find_np()
{
  switch (npstore.size()) {
  case 1:
    if (verbose & VER_GEN)
      std::cout << "Setting np from experimental data\n";
    v_np=npstore.front();
    return true;
  case 0:
    break;
  default:
    autoset_notimp_warning.raise();
  }
  return false;
}
    
ContextWarning<> fit_expafterprevious_warning("fit exp used after fit add/exp.  Previous data sets / data masks will be ignored",&NMRsim_repeat_warning);
ContextWarning<> fit_nofilenames_warning("no filenames specified in fit exp; fitting disabled.",&NMRsim_repeat_warning);
ContextWarning<> fit_overwritingdomain_warning("time/frequency domain specification overwriting that of file read",&NMRsim_repeat_warning);
ContextWarning<> fit_domainunclear_warning("unsure how data point reduction should affect spectral width so no scaling applied. Use explicit -frequency flag if spectral width needs scaling",&NMRsim_repeat_warning);

void parse_fit_exp()
{
  ensure_fit();
  if (!fit_sets.empty() || !mask.empty()) {
    fit_sets.clear();
    mask.clear();
    fit_expafterprevious_warning.raise();
  }
  const char* fname;
  cmatrix tmp_set;
  filestruct info;
  while (parser_isnormal() && (fname=parse_string(F_REPLACEDOLLAR | F_ALLOWMISSING))) {    
    read_file(tmp_set,info,fname);
    const processing_state pstate(info.sw,info.domain==D_TIME,info.sfrq,info.ref);
    const size_t lnp=tmp_set.cols(); //!< need to preserve as tmp_set destroyed!
    const size_t lni=tmp_set.rows();
    const size_t setno=add_data_set(tmp_set,pstate);
    makevars(setno,lnp,lni,info);
  }
  if (fit_sets.empty())
    fit_nofilenames_warning.raise();
}

void parse_fit_add()
{
  if (!(mask.empty()))
    error_abort("Can't add data sets after fit mask");
  ensure_fit();
  const char* fname=parse_string(F_REPLACEDOLLAR);

  LIST<std::string> readvarnames; //!< read list of items to read in
  while (parser_isvariable_name())
    readvarnames.push_back(clean_variable_name(parse_string()));

  cmatrix tmp_set;
  filestruct info;
  read_file(tmp_set,info,fname,&readvarnames);

  bool doreverse=false;
  if (are_left()) {
    double scalef=1.0;
    if (parser_isnormal()) {
      const size_t oldnp=tmp_set.cols();
      check_extract(tmp_set);
      scalef=double(tmp_set.cols())/oldnp; //!< scale spectral width
    }
    const std::pair<domain_t,bool> domflags=getdomain();
    const domain_t sdom=domflags.first;
    doreverse=domflags.second;
    if (doreverse) {
      used_reverse=true;
      for (size_t r=tmp_set.rows();r--;)
	rawreverse(tmp_set.row(r));
    }
    if (sdom!=D_UNKNOWN) {
      if ((info.domain!=D_UNKNOWN) && (info.domain!=sdom))
	fit_overwritingdomain_warning.raise();
      info.domain=sdom;
    }
    if (scalef!=1.0) {
      switch (info.domain) {
      case D_TIME: break;
      case D_UNKNOWN:
	//	(info.sw)*=scalef;
	fit_domainunclear_warning.raise();
	break;
      case D_FREQUENCY:
	(info.sw)*=scalef;
	break;
      }
    }
  }    
  if (tmp_set.cols()==1)
    tmp_set.transpose();
  const size_t lnp=tmp_set.cols(); //!< store value as tmp_set is destroyed
  const size_t lni=tmp_set.rows();
  const bool istd=(info.domain==D_TIME);
  const size_t setno=add_data_set(tmp_set,processing_state(info.sw,istd,info.sfrq,info.ref));
  makevars(setno,lnp,lni,info);
  if (verbose & VER_OPTIM)
    std::cout << "Read in file " << fname << "  SW: " << (info.sw*1e-3) << " kHz  Time domain: " << (istd ? "Yes" : "No") << "  Reversed: " << (doreverse ? "Yes" : "No") << '\n';
}

void parse_fit()
{
  if (isinteractive)
    error_abort("Can't add new fit directives after interactive");
  
  if ((opttype!=OPTIM_NONE) && (opttype!=OPTIM_FIT))
    error_abort("Can't change fitting / optimisation type");

  Fit_Factory_t& fit_Factory(get_fit_Factory());
  GlobalFit_Factory_t& setup_fit_Factory(get_setup_fit_Factory());
  static bool doneadd=false;
  if (!doneadd) {
    doneadd=true;
    setup_fit_Factory["exp"]=&parse_fit_exp;
    setup_fit_Factory["add"]=&parse_fit_add;
    setup_fit_Factory["noiselevel"]=&parse_fit_noiselevel;

    fit_Factory["tolerance"]=tolerance_desc;
    fit_Factory["iterations"]=iterations_desc;
    fit_Factory["method"]=method_desc;
    fit_Factory["release"]=release_desc;  
    fit_Factory["statistics"]=statistics_desc;
    fit_Factory["normalize"]=normalise_desc;
    fit_Factory["normalise"]=normalise_desc;
    fit_Factory["constrain"]=constrain_desc;
    fit_Factory["unconstrain"]=unconstrain_desc;
    fit_Factory["fix"]=fix_desc;
    fit_Factory["run"]=run_desc;
    fit_Factory["set"]=set_desc;
    fit_Factory["system"]=system_desc;
    fit_Factory["error"]=error_desc;
    fit_Factory["abort"]=abort_desc;
    fit_Factory["help"]=help_desc;
    fit_Factory["mask"]=mask_desc;
    fit_Factory["finish"]=finish_desc;
    fit_Factory["interactive"]=interactive_desc;

    const bool fitverb = (verbose_level>1) && (verbose & VER_OPTIM);
    if (!silent || !nochecks)
      optimisation_base_warning.type(fitverb ? BaseWarning::Always : BaseWarning::FirstOnly);
  }

  const char* keyname=parse_string(F_REPLACEDOLLAR);
  const GlobalFit_Factory_t::const_iterator iter(setup_fit_Factory.find(keyname));

  if (iter!=setup_fit_Factory.end()) {
    if (haverun)
      error_abort("can't modify data set parameters after run directive");
    (iter->second)();
    return;
  }
  
  make_optim_command(keyname,fit_Factory);
}
  
ContextWarning<> fit_ignoringnoiselevel_warning("ignoring previous setting of data set noiselevel",&NMRsim_once_warning);
ThreadWarning<> notinteraction_warning("fitted eta parameter is not part of an NMR interaction",&NMRsim_repeat_warning);

void parse_fit_noiselevel()
{
  parse_system_variable(v_fit_noise,F_DENYSUM | F_DENYVAR);
  if (fit_noise<0.0)
    error_abort("fit noise must be >0 (or zero if not specified)");
  static bool done=false;  
  if (done)
    fit_ignoringnoiselevel_warning.raise();
  else
    done=true;
}

void adjustconstraint(VarVariable& cvar, bool constraineta)
{
  Parameter& cpar(cvar.parameter);
  if (constraineta) {
    if (cvar.boundstate.lowerupper(0.0,1.0))
      cpar.constrain(*(cvar.boundstate.function()));
  }
  else {
    if (cvar.boundstate.unconstrain())
      cpar.unconstrain();
  }
}

size_t prepare_paras(bool lsilent, bool nomods)
{    
  const size_t npars=parameter_list.size();  
  
  set_seed(); //randomise random number generator
  
  //   const varvarpars_t::iterator end(var_varpars.end());
  //   varvarpars_t::iterator start(var_varpars.begin());
  //   const size_t npars=var_varpars.size();
  //   dest.create(npars);
  
  size_t varcount=0;
  const bool constraineta=!nomods && !(opt_flags & OPT_NOCONSTRAINTS);

  for (size_t i=0;i<npars;i++) {
    VarVariable& cvar(*(var_varpars(i)));
    Parameter& cpar(cvar.parameter);
    switch (cvar.subsid) {
      case S_ASYM: {
	const Interaction* interp=dynamic_cast<Interaction*>(cvar.ptr);
	if (interp==NMRSIM_NULL)	
	  notinteraction_warning.raise();
	else {
	  if (!(interp->isxy()))
	    adjustconstraint(cvar,constraineta);
	}
      }
	break;
    case S_GFRAC: case S_GFRAC1:
      adjustconstraint(cvar,constraineta);
      break;
    }
    cpar.set(cvar.get()); //!< ensure parameter contains same value
    // 	interp->convert_anisotropy(cvar);
    //    rawpara& cdest(dest(i));
    
//     if (converteta && (cvar.subsid==S_ETA)) { //turn ETAs into XYs
//       const Interaction* interp=dynamic_cast<Interaction*>(cvar.ptr);
//       if (interp==NMRSIM_NULL)
// 	notinteraction_warning.raise();
//       else {
// 	interp->convert_anisotropy(cvar);
//       }
//     }
//     cdest.value=cvar.get();
//     cdest.step=cvar.step;
//     cdest.isconst=cvar.isfixed;
    if (!(cvar.isfixed())) {
//       if (step==0.0) {
// 	if (val==0.0)
// 	  error_abort("Step unset with zero value parameter");
// 	step=val*default_step_scale;
//       }
//       if (cvar.step==0.0)
// 	throw InternalError("Fitting variable with zero error");
      if (!nomods && (opt_flags & OPT_RANDOMISE))
	cvar.set(cvar.get()+cvar.error()*(random(2.0)-1.0));
      //      actord.push_back(varcount);
      varcount++;
   }
    
    //    cdest.value=val;
    //    vals.push_back(val);
    //    cvar.set(cdest.value);
    //    errs.push_back(step);
    //    cdest.step=step;
    //    ++start;
  }
  if (varcount==0)
    throw Failed("Can't fit/optimise with no active variables!");
  
//   if ((verbose & VER_OPTIM) && (verbose_level>1)) {
//     std::cout << "Parameter steps:";
//     for (size_t i=0;i<npars;i++) {
//       if (!(var_varpars(i).isfixed()))
// 	std::cout << ' ' << var_varpars(i).step;
//     }
//     std::cout << '\n';
//   }
  
  if (!lsilent) {
    std::cout << "\nFitting/optimisation parameters: " << npars << '\n';
    //    dump(var_varpars);
    dump_parameters();
    std::cout << '\n';
  }
  return varcount;
}

#ifdef USE_MINUIT 

static bool min_invert=false;
 
class FitExternal : public FitCommand
{
public: 
  explicit FitExternal(OptimiseFunction);
  void exec() const { min_externalp=min_externalp_; }
  static FitCommand* create();
private:
  OptimiseFunction min_externalp_;
};

FitExternal::FitExternal(OptimiseFunction min_externalpv)
  : min_externalp_(min_externalpv) 
{
  if (!min_externalpv)
    throw Failed("NMRSIM_NULL optimisation function!");    
  min_funcset=true;
}

FitCommand* FitExternal::create()
{
  const OptimiseParse parsefuncp=findmap(get_Optimise_Factory(),parse_string(F_REPLACEDOLLAR),"optimisation function");
  OptimiseFunction optfuncp((*parsefuncp)());
  return new FitExternal(optfuncp);
}

class FitSum : public FitCommand
{
public:
  FitSum() { min_funcset=true; }
  explicit FitSum(const BaseList<size_t>& cselv, const BaseList<size_t>& rselv =LIST<size_t>()) : col_sel_(cselv), row_sel_(rselv) { min_funcset=true; }
  void exec() const;
  static FitCommand* create();
private:
  LIST<size_t> col_sel_,row_sel_;
};

FitCommand* FitSum::create()
{
  if (!are_left())
    return new FitSum();
  const LIST<size_t> col_sel(parse_unsignedintarray(1));
  LIST<size_t> row_sel;
  if (are_left())
    row_sel=parse_unsignedintarray(1);
  return new FitSum(col_sel,row_sel);
}

void FitSum::exec() const
{
  col_sel=col_sel_;
  row_sel=row_sel_;
  min_externalp=NMRSIM_NULL;
}

 void parse_minmax(int invert)
{
//   if (global_workers)
//     error_abort("Sorry! optimisation/fitting is not currently compatible with multi-processor operation!");

  if (isinteractive)
    error_abort("Can't add new optimisation directives after interactive");

  switch (opttype) {
  case OPTIM_FIT:
    error_abort("can't combine fitting and optimisation");
  case OPTIM_OPT:
    if (min_invert!=invert)
      error_abort("can't switch between maximise and minimise");
    break;
  case OPTIM_NONE:
    min_invert=invert;
    opttype=OPTIM_OPT;
    replace_callback();
    break;
  }

  Fit_Factory_t& minmax_Factory(get_minmax_Factory());
  static bool doneadd=false;
  if (!doneadd) {
    doneadd=true;
    minmax_Factory["tolerance"]=tolerance_desc;
    minmax_Factory["precision"]=precision_desc;
    minmax_Factory["iterations"]=iterations_desc;
    minmax_Factory["method"]=method_desc;
    minmax_Factory["release"]=release_desc;    
    minmax_Factory["fix"]=fix_desc;
    minmax_Factory["set"]=set_desc;
    minmax_Factory["system"]=system_desc;
    minmax_Factory["constrain"]=constrain_desc;
    minmax_Factory["unconstrain"]=unconstrain_desc;
    minmax_Factory["statistics"]=statistics_desc;
    minmax_Factory["error"]=error_desc;
    minmax_Factory["abort"]=abort_desc;
    minmax_Factory["help"]=help_desc;
    minmax_Factory["finish"]=finish_desc;
    minmax_Factory["run"]=FitCommand_t(&FitMinuitRun::create,"run\t - run optimisation");
    minmax_Factory["sum"]=FitCommand_t(&FitSum::create,"sum [<column range> [<row range>]] - use sum of data set (or selected region) as optimisation function");
    minmax_Factory["external"]=FitCommand_t(&FitExternal::create,"external <expression>\t - specify (externally parsed function) for optimisation");
    minmax_Factory["interactive"]=interactive_desc;
  }

  make_optim_command(parse_string(F_REPLACEDOLLAR),minmax_Factory);
}

class MinObj : public BaseMinFunction {
 public:
  MinObj(MasterObj& obj_)
    : obj(obj_) { reset(); }
  double operator()(const BaseList<double>& pars) const;
  void reset() { fcncount=0; }

 private:
  MasterObj& obj;
  mutable size_t fcncount;
  mutable DataStore tmp;
};

double MinObj::operator()(const BaseList<double>& pars) const
{
  int lverbose = (verbose & VER_OPTIM) ? verbose_level : 1;
  if (silent)
    lverbose=0;
  fcncount++;
  if (lverbose) {
    std::cout << fcncount << ": " << pars;
    std::cout.flush();
  }
  set_variables(pars);
  obj.calc(tmp,0);
  if (tmp.empty())
    error_abort("Data set is empty - no actual acquisition?");
  double val;
  if (min_externalp)
    val=(*min_externalp)(tmp);
  else {
    if (tmp.rows()>1) {
      if (!tmp.isnD())
	error_abort("can't apply 2D sum to irregular data set");
      const cmatrix& tmpmat(tmp.matrix());
      val= real(row_sel.empty()
		? sum(tmpmat(range(),col_sel))
		: sum(tmpmat(row_sel,col_sel)));
    }
    else {
      if (!row_sel.empty())
	error_abort("minimise: 1D data set but 2D sum");
      
      const BaseList<complex>& asrow(tmp.row());
      val=real(col_sel.empty()
	       ? sum(asrow)
	       : sum(asrow(col_sel)));
    }
  }
  if (lverbose)
    std::cout << "\t " << val << '\n';
  return min_invert ? -val : val;  
}

void FitMinuitRun::exec() const
{
  (void)prepare_paras(silent);

  LCM_MINUITNAMESPACE::MnUserParameters mnparas;
#if MINUITVER>1
  if (minuit_precision)
    mnparas.SetPrecision(minuit_precision);
#endif

  for (size_t i=0;i<var_varpars.size();i++) {
    const VarVariable& cvar(*(var_varpars(i)));
    const char* cname=cvar.name();
    if (cvar.isfixed())
#if MINUITVER==2
      mnparas.Add(cname,cvar.get());
    else {
      mnparas.Add(cname,cvar.get(),cvar.error());
#else
      mnparas.add(cname,cvar.get());
    else {
      mnparas.add(cname,cvar.get(),cvar.error());
#endif
      cvar.boundstate.apply(mnparas,i);
    }
  }

    if (objp==NMRSIM_NULL)
      throw InternalError("FitMinuitRun::exec");

  static MinObj minobj(*objp);
  static MinuitAdaptor<MinObj> theFCN(minobj);

  static int lastopt_method=-1;
  static smartptr<LCM_MINUITNAMESPACE::ModularFunctionMinimizer,false> minimizerp;

  if (lastopt_method!=opt_method) {
    switch (opt_method) {
    case OPT_GRADIENT:
      minimizerp.reset(new LCM_MINUITNAMESPACE::VariableMetricMinimizer);
      break;
    case OPT_SIMPLEX:
      minimizerp.reset(new LCM_MINUITNAMESPACE::SimplexMinimizer);
      break;
    default:
      throw InternalError("FitMinuitRun");
    }
    lastopt_method=opt_method;
  }
  minobj.reset(); //! reset function count
  //do optim
  const LCM_MINUITNAMESPACE::MnStrategy mnstrategy(NMRSIM_OPTIM_STRATEGY);
  const int maxevals=get_maxiter();
  const double tol=get_tolerance()*1e3; //convert 0..1 tolerance into scale that Minuit works with
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
    if (!silent) {
      if (isvalid)
	std::cout << "Optimisation converged after " << iters << " evaluations\n";
      else
	std::cout << "Optimisation unfinished after maximum evaluations\n";
    }
    globalminp.reset(new LCM_MINUITNAMESPACE::FunctionMinimum(min));
    haveresults=true;
  }
  else
    std::cerr << "Optimisation did not converge or failed (optimisation ignored)\n";
}

ContextWarning<> fit_convergefailed_warning("optimisation did not properly converge",&NMRsim_repeat_warning);

bool perform_optim(MasterObj& obj, DataStore& Spec)
{
  objp=&obj; //!< copy out MasterObj pointer
  execfitstack();

  if (!haveresults || !globalminp)
    return false;
#if MINUITVER==2
  const LCM_MINUITNAMESPACE::MnUserParameterState& parstate(globalminp->UserState());  
  const bool isvalid=globalminp->IsValid();
  final_optimisation=globalminp->Fval();
#else
  const LCM_MINUITNAMESPACE::MnUserParameterState& parstate(globalminp->userState());  
  const bool isvalid=globalminp->isValid();
  final_optimisation=globalminp->fval();
#endif
  if (!isvalid)
    fit_convergefailed_warning.raise();
  if (min_invert)
    final_optimisation=-final_optimisation; //reverse sign change

  if (!silent) {
    std::cout << "Optimised parameters:\n";
    for (size_t i=0;i<var_varpars.size();i++) {
      std::cout << (i+1);
      const VarVariable& cvar(*(var_varpars(i)));
      if (cvar.isnamed())
	std::cout << " (" << cvar.name() << ')';
      std::cout << ": " << 
#if MINUITVER==2
	parstate.Value(i)
#else
	parstate.value(i)
#endif
		<< '\n';
    }
    // << " +/- " << parstate.error(i) << '\n'; //! error surely is not meaningful?
    std::cout << "\nOptimised value: " << final_optimisation << '\n';
  }

  set_variables(parstate);
  
  if (!finalisestacks.empty())
    (void)obj.calc(Spec,verbose,false); //may need save, so calc optimised data set (no need to share result)

  return true;
}
#else
 void parse_minmax(int)
{
  error_abort("Optimisation support (Minuit) not enabled");
}

bool perform_optim(MasterObj&, DataStore&) //shouldn't be called
 { throw InternalError("perform_optim"); }
#endif

void dump_parameters(bool onlyvars, bool ignoreconstraints)
{
  //  static double asymtol=1e-8;
  const bool verb=(verbose & VER_OPTIM);

  size_t active=0;
  for (size_t i=0;i<var_varpars.size();i++) {
    VarVariable cvar(*(var_varpars(i)));
    size_t ind=i;
    if (onlyvars) {
      if (cvar.isfixed())
	continue;
      ind=active++;
    }
    std::cout << (ind+1) << ' ';
    cvar.printname();
    const double val=cvar.get();
    std::cout << ":\t " << val;
//     Interaction* interp=NMRSIM_NULL;
//     if (cvar.subsid==S_ASYM) { //report XY as ETA
//       interp=dynamic_cast<Interaction*>(cvar.ptr);
//       if (interp) {
// 	if (interp->isxy() || ((val>-asymtol) && (val<1+asymtol)))
// 	  interp=NMRSIM_NULL; //!< no problem
//       }
//       else
//  	notinteraction_warning.raise();
//       const double asym=interp->asymmetry();
// //       std::cout << "[Asymmetry: " << asym;
// //       if (cvar.step)
// // 	std::cout << " +/- " << ((cvar.step)/interp->anisotropy()) << " (approx)";
// //       std::cout << ']';
// //       if (
// // 	interp=NMRSIM_NULL; //! everything OK	
//     }
    if (cvar.isfixed())
      std::cout << "  (fixed)\n";
    else {
      std::cout << "\t +/- " << cvar.error();
      if (!ignoreconstraints) {
	std::cout << "  (" << cvar.boundstate << ')';
	if (verb && cvar.isconstrained())
	  std::cout << "  Internal value: " << cvar.parameter.get_internal();
      }
      std::cout << '\n';
    }
//     if (interp) { //bad asymmetry
//       ScratchList<double,3> carts(3U);
//       anisotropy_to_cartesian(carts(0U),carts(1U),carts(2U),0.0,interp->anisotropy(),interp->asymmetry());
//       sort(carts);
//       double junk,laniso,lasym;
//       cartesian_to_anisotropy(junk,laniso,lasym,carts(1U),carts(0U),carts(2U));
//       std::cout << "[ Invalid asymmetry - reordering principal values gives anisotropy=" << laniso << ", asymmetry=" << lasym << " (Euler angles will need changing)]\n";
//     }
  }
}
 
void build_fitset()
{
  if (fit_set.empty()) {
    raw_build_set(fit_set,fit_sets,fit_procs,(verbose & VER_OPTIM) ? verbose_level : 0);
    FitNormalise::apply(fit_set); //!< apply any normalisation
    
    if (fit_verbose_level()) {
      std::cout << "Fit selection: ";
      if (mask.empty())
	std::cout << "None\n";
      else {
	if (verbose_level>1)
	  std::cout << mask << '\n';
	else
	  std::cout << mask.size() << " points\n" ;
      }
    }
  }
}

void calc_parameters(rmatrix& lvalerrs)
{
  const bool have_covar=!!covar;
  lvalerrs.create(var_varpars.size(),2);
  size_t varind=0;
  if (have_covar && (fit_verbose_level()>1))
    std::cout << "Updating standard errors using covariance matrix\n" << covar << '\n';

  for (size_t ind=0;ind<var_varpars.size();ind++) {
    VarVariable& cvar(*(var_varpars(ind)));
    lvalerrs(ind,0U)=cvar.get();
    const bool haveerr=have_covar && !(cvar.isfixed());
    const double lerr=haveerr ? sqrt(covar(varind,varind)) : 0.0;
    lvalerrs(ind,1U)=lerr;
    cvar.error(lerr);
    if (haveerr)
      varind++;
  }
}

ThreadWarning<> emptydata_warning("Signal-to-noise ratio appears to be <2. Is the noiselevel parameter too large or is there no signal to fit?",&NMRsim_once_warning);
ThreadWarning<> badmask_warning("Signal-to-noise ratio of data selection appears to be <2. Is the noiselevel parameter too large or is the wrong section being selected by fit mask?",&NMRsim_once_warning);

void FitRun::exec() const
{
  if (isinter) {
    runinteractive();
   return;
  }

  prepare_paras(silent);
  build_fitset();

  const BaseList<complex> source(fit_set.row());
  List<double> rawdata(fitnpts(source.size()));
  extractdata(rawdata,source);

  if (!nochecks && (fit_noise!=0)) {
    double datamax=max(rawdata);
    const double datamin=std::abs(min(rawdata));
    if (datamin>datamax)
      datamax=datamin;
    if (datamax<2*fit_noise) {
      if (mask.empty())
	emptydata_warning.raise();
      else {
	badmask_warning.raise();
	if (verbose & VER_OPTIM)
	  std::cout << "Raw masked fitting data: " << rawdata << '\n';
      }
    }
  }
      
  const int fitverb=fit_verbose_level();
  rmatrix newcovar;

  const size_t actiter=get_maxiter();
  const double acttol=get_tolerance();
  if (fitobjp==NMRSIM_NULL)
    throw InternalError("FitRun::exec");
  ensure_noisevector(source.size());
  const bool havenoise=havenoisespec();
  const bool isinter=(flags & allowinteractive);

  try {
    optimisation_info info;

    switch (opt_method) {
    case OPT_GRADIENT:
      info=fitdata(newcovar,parameter_list,*fitobjp,rawdata,noisevector,fitverb,acttol,actiter);
      break;
    case OPT_SIMPLEX: 
      do {
	covar.clear();
	info=simplex_fitdata(newcovar,parameter_list,*fitobjp,rawdata,noisevector,fitverb,acttol,actiter);
      } while (isinter && getlogical("Continue? ",1));
      break;
    default:
      throw InternalError("perform_fit");
    }
    if (!silent && (actiter!=0)) {
      if (info.converged) {
	std::cout << "Fitting converged in " << info.iterations << " iteration";
	std::cout << ((info.iterations>1) ? "s\n" : "\n");
      }
      else
	std::cout << "Fitting did not converge\n";
    }
    global_chisqr=info.optimum;
    set_variables();
    covar.swap(newcovar); //!< only swap in if successful
    if (!havenoise)
      covar*=global_chisqr; //adjust covariance matrix

    if (info.converged && !silent) {
      if (havenoise)
	std::cout << "Normalised chi^2 (should be ~1): " << global_chisqr << " (using specified noiselevel of " << fit_noise << ")\n";
      else
	std::cout << "N.B. cannot use chi^2 to assess quality of fit (noise level unspecified)\n";
    }
    haveresults=true;
  }
  catch (SingularMatrix<double>& exc) {
    const rmatrix& Hessian(exc.matrix());
    bool doneintro=false;
    size_t active=0;
    for (size_t i=0;i<var_varpars.size();i++) {
      VarVariable cvar(*(var_varpars(i)));    
      if (cvar.isfixed())
	continue;
      const double Hii=Hessian(active,active);
      if ((Hii==0.0) || NMRSIM_ISNAN(Hii)) {
	if (!doneintro) {
	  std::cerr << "Hessian matrix was singular.  Simulation appears to be invariant to the following parameter(s):\n";
	  doneintro=true;
	}
	std::cerr << (i+1) << ' ';
	cvar.printname(std::cerr);
	std::cerr << '\n';
      }
      active++;
    }
    if (!doneintro) {
      std::cerr << "Inversion of Hessian failed without obviously invariant parameters: contents in HessianDump\n";
      spy(std::cerr,Hessian);
      try {
	write_matrix("HessianDump",Hessian); //trap errors
      } catch (...) {}
    }
  }
  catch (std::exception& exc) {
    std::cerr << "Fitting failed (" << exc.what() << ") - results discarded\n";
  }
}

ThreadWarning<> fit_constraintsactive_warning("error bounds are more doubtful in presence of constraints (particularly if parameter values are close to their limits).  Consider re-calculating from converged parameters without constraints or using fit statistics -release_constraints in finalise.",&NMRsim_once_warning);

void ensure_dataset(const DataStore& Spec)
{
  static bool doneensure=false;  
  if (doneensure)
    return;
  doneensure=true;

  if (!arematching(Spec,fit_set)) {
    std::cerr << "Supplied and calculated data sets differ in size:\n";
    std::cerr << "Supplied data set: ";
    fit_set.print_structure(std::cerr);
    std::cerr << "\nCalculated data set: ";
    Spec.print_structure(std::cerr);      
    std::cerr << '\n';
    error_abort();
  }
  checksw(fit_set,Spec);
}

void display_statistics(bool ignoreconstraints)
{
  std::cout << "\nFitting parameters and errors\n";
  dump_parameters(true,ignoreconstraints);
  std::cout << '\n';
  bool conwarning=false;
  LIST<size_t> actord(var_varpars.size());
  actord.create(0);
  for (size_t i=0;i<var_varpars.size();i++) {
    const VarVariable& cvar(*(var_varpars(i)));
    if (!(cvar.isfixed())) {
      actord.push_back(i);
      if (!ignoreconstraints && cvar.isconstrained())
	conwarning=true;
    }
  }
  if (!havenoisespec()) {
    if (global_chisqr)
      fit_randomassumed_warning.raise();
    else
      fit_nonoiselevel_warning.raise();
  }
  if (conwarning)
    fit_constraintsactive_warning.raise();

  const size_t varind=covar.rows();
  if (varind==0)
    return;

  if (varind!=actord.size())
    throw InternalError("display_statistics");
  if (varind>1) {
    std::cout << "\nCorrelation matrix\n";
    for (size_t j=0;j<varind-1;j++)
      std::cout << '\t' << (actord(j)+1);
    std::cout << '\n';
    for (size_t i=1;i<varind;i++) {
      std::cout << (actord(i)+1) << ": ";
      for (size_t j=0;j<i;j++)
	std::cout << covar(i,j)/sqrt(covar(i,i)*covar(j,j)) << " \t";
      std::cout << '\n';      
    }
  }
  std::cout << '\n';
}

ThreadWarning<> fit_noscale_warning("no variable scaling or addsignals parameter in fitting.  Suggest adding scaling command into proc",&NMRsim_once_warning);
ThreadWarning<> fit_scale_warning("Variable scaling in processing when data is being normalised.  Fitting will probably fail due to this dependent parameter.",&NMRsim_once_warning);
ThreadWarning<> fit_zerospectrum_warning("Calculated spectrum has zero maximum!",&NMRsim_once_warning);

bool perform_fit(MasterObj& obj, DataStore& Spec)
{  
  //const int verb = (verbose & VER_OPTIM) ? verbose_level : 0;

//   if ((fit_sets.size()==1) && fit_sw && (fit_sw!=sw))
//     std::cerr << "Warning: spectral width of fitted and calculated sets seem to differ\n";

//  build_noisevector();
  build_fitset();
  build_fitobjs(obj);

  if (FitNormalise::isactive()) {
    if (havescale)
      fit_scale_warning.raise();
  }
  else {
    if (have_spinsys() && !havescale && fit_noscale_warning.enabled() && (maxiter!=0)) {
      obj.calc(Spec,0);
      ensure_dataset(Spec);
      const double maxexp=max(real(fit_set.row()));
      const double maxcalc=max(real(Spec.row()));
      if ((verbose & VER_OPTIM) && (verbose_level>1))
	std::cout << "Max experimental: " << maxexp << "\tMax calculated: " << maxcalc << '\n';
      if (maxcalc==0)
	fit_zerospectrum_warning.raise();
      else {
	char buf[256];
	snprintf(buf,sizeof(buf),"\n\tscale %gV",maxexp/maxcalc);
	fit_noscale_warning.raise(buf);
      }
    }
  }

  execfitstack();
  if (!haveresults)
    return false;
  
  final_chisquared=global_chisqr;

  calc_parameters(valerrs);
  if (!silent && (maxiter!=0))
    display_statistics();

  obj.calc(Spec,verbose,false); //!< no need to share
  if (ammaster) {
    updatevarerrors();
    getvarnames();
    residuals=fit_set;
    residuals-=Spec;
    if (!silent) {
      if (residuals.row().size()>1) {//!< little paranoid - hard to imagine fitting with <2 points!
	LIST<double> extractresid;
	create_residuals(extractresid);
	std::cout << "\nStandard deviation of residuals (noiselevel if random): " << stddev(extractresid) << '\n';
      }
      std::cout << '\n';
    }
  }
  return true;
}

std::ostream& operator<< (std::ostream& ostr, const VarVariable& a)
{
  if (a.variable().isarray())
    ostr << "Element " << a.arraywhich << " of ";
  a.printname(ostr);
  ostr << '=' << a.get();
  if (a.isfixed())
    return ostr << " (fixed)\n";

  ostr << " +/- ";
  if (a.error())
    ostr << a.error();
  else
    ostr << "<default>";
  return ostr << '\n';
}

namespace {
  const BaseList<double> validate_array(const VariableBase& var)
  {
    if (var.array_item_length()!=1)
      error_abort("Can't perform fitting on list of values");
    return var.get_array().row();
  }
}

double VarVariable::get() const
{
  if (variable().isarray()) {
    const BaseList<double> array(validate_array(variable()));
    return array(size_t(arraywhich));
  }
  return Variable::get();
}

void VarVariable::set(double fval)
{
  if (variable().isarray()) {
    BaseList<double> array(validate_array(variable()));
    array(size_t(arraywhich))=fval;
    return;
  }
  if (!ptr)
    throw InternalError("Variable: ptr unset");
  if (!valuep->issimple())
    throw Failed("Can only set simple variables");
  BaseList<double>& value(valuep->value());
  if (value.size()!=1)
    throw InternalError("Variable: can't set list variable");

  if (fval!=value.front()) {
    ptr->set(fval,subsid);
    value.front()=fval;
  }
}

void build_fitobjs(MasterObj& obj)
{
  if (!fitobjp) {
    noisevals.create(size_t(0));
    if (v_fit_noise.isdefined() && !v_fit_noise.isconstant()) {
      if (v_fit_noise.uses()!=A_ARRAY)
	error_abort("fit_noise is not a simple number or {} array");
      array_iter iter;
      while (iter.next())
	//      noisevals.create(noisevals.size()+np,fit_noise);
	noisevals.push_back(fit_noise);
    }
    fitobjp=new FitObj(obj);
  }
}

void VarVariable::error(double err)
{
  if (err)
    parameter.error(err); //!< update both Parameter and VariableBase (for 'Errorof'), but error can't be zero for parameter
  variable().set_error(err,arraywhich);
}

