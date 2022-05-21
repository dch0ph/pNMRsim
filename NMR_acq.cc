/* routines particular to complex simulations */

// first load config for USEMPI
// Set up MPI stuff here to avoid conflicts after using libcmatrix 
#include "config.h"
#if defined(HAVE_LIBMKL) && defined(LCM_USE_EXTERNAL)
#include "mkl.h"
#define USEMKL 1
#else
#undef USEMKL
#endif

#include "NMRsim.h"

#ifndef DISABLE_EXEC
#ifdef USEMPI
#include "cmatrix_MPI.h"
#define PARTYPE libcmatrix::MPI_controller
#else
#ifdef HAVE_FORK_CONTROLLER
#include "Fork_controller.h"
#define USE_FORK_CONTROLLER
#define PARTYPE libcmatrix::Fork_controller
#endif
#endif
#endif

#ifdef HAVE_PARSYS
#include "smartptr.h"
libcmatrix::smartptr<PARTYPE,false> parconp;
bool assumeparallel=true; //!< assume parallel computation unless impossible

size_t get_thread_num() { 
  return !parconp ? 0 : parconp->get_thread_num();
}
#else
size_t get_thread_num() { 
  return 0;
}
#endif

#include "NMRsim_Process.h"
#include "NMRsim_spinsys.h"
#include "NMRsim_logfile.h"
#include "Action.h"
#include "powder.h"
#include <sstream>

FILE* transitions_fp=NMRSIM_NULL;
const char* lastlogfilename=NMRSIM_NULL;

size_t global_workers=0; // number of *workers* (0 if unset)
Matrix<accumulating_timer<> >* timers_arrayp;
LIST<BaseWarning*> resetwarning_list;

size_t global_np=0; //!< finalise data set size not fixed
LIST<size_t> global_nps;
static bool powder_print=false; //!< if set to true, force print number of powder orientations each time (re-)created

//! cutoff warning trigger (%)
#define CUTOFF_TRIGGER 1.0

std::pair<double,double> actual_swref()
{
  return (hist_min==hist_max) 
    ? std::pair<double,double>(sw,0.0)
    : std::pair<double,double>(hist_max-hist_min,0.5*(hist_min+hist_max));
}

static DataStore ResultSpec; //!< temporary object used by distributed sums
static DataStore OrientSpec; //!< temporary object used for individual orientation

namespace {

void init_row(DataStore& a)
{
  const bool istd=!isfrequency();
  if (!isirregular()) {
    if (a.empty())
      a.create(array_n0,np,complex(0.0));
  }
  else {
    if (a.empty())
      a.create(array_n0);
    if (a.rows()<=row_index) {
      const std::pair<double,double> swref(actual_swref());
      a.createrow(np,complex(0.0),processing_state(swref.first,istd,detect_freq,swref.second));
    }
    else {
      if (a.state(row_index).sw!=sw) {
	char buf[256];
	snprintf(buf,sizeof(buf),": row %lu",(unsigned long)(row_index+1));
	changedsw_warning.raise(buf);
      }
    }
    processing_state& cstate(a.state(row_index));
    cstate.istimedomain=istd;    
  }
}

  void create_empty(DataStore& a)
  {
    if (!isirregular()) {
      if (np==0)
	a.clear();
      else
	a.create(array_n0,np,complex(0.0));
    }
    else {
      a.clear();
      array_iter variter;
      while (variter.next())
	init_row(a);
    }
  }  

  void make_zero(DataStore& Spec)
  {
    if (Spec.empty())
      create_empty(Spec);
    else
      Spec=complex(0.0);
  }

}

enum powder_t { POW_NONE, POW_ZCW, POW_BETA, POW_FILE, POW_ZCW3, POW_ALPHABETA, POW_SPHERICALZCW };
powder_t powder_type=POW_NONE;
char* crystal_filep=NMRSIM_NULL;
void parse_powderquality();

void parse_histogram(int);
void parse_crystal_file();

const BlockedOperator empty_op;
const BlockedOperator* use_detectp=NMRSIM_NULL;
bool detectED=false;

option optgamma("gammacompute");
option optparallel("parallel");
option optED("EDmatching","matching initial density matrix and detect operator");
option optforceeigenbasis("forceeigenbasis","",option::AUTO,option::NOTUSED);
option optforcepointbypoint("forcepointbypoint","",option::AUTO,option::NOTUSED);
//optional_t optlongdtsync;

static const double rad_to_deg=180.0/M_PI;
static const double deg_to_rad=M_PI/180.0;

static const double dotlimit=2.5; //show progress counter if > this time (s)
static const double siglimit=10.0; //show ETA if > this limit (s)
static size_t use_gamma_angles=0; //!< gamma integration steps (0 if not gamma-COMPUTE)
static bool update_gamma=false;
static double detect_period=0.0;
static double restrict_gamma=2*M_PI;
static double globalscale=0.0;
dirty_t update_propagators=DIRTY_ALL; //flag that cached propagators shouldn't be used
size_t orientations=1; //!< if not set, implicitly a single orientation
int orientation=-1;
int nzcw=0;
bool zcwisindex=true; //!< by default always treat as index
range_t rangequal=sphere;
void ensure_basepowder();
bool powder_array=false;
bool powder_reverse=false;
static rmatrix powder_weights;
static int curnzcw=-1;
static int nalpha=1;
//smartptr<PowderMethod> pregammapowdm; //!< powder averaging method before any gamma_angles
smartptr<PowderMethod> powdm; //!< powder averaging method
bool have_orientation=false;
bool process2D=true;
static double gammaX=0.0;

rmatrix vararray;
rmatrix valerrs;
LIST<const char*> varnames;

ThreadWarning<> shortCW("Continuous wave sequence with short duration may lead to inefficiences - consider using effectively infinite duration e.g. 1e6 for sequence: ",&NMRsim_once_warning);

double get_rf_period(const CompSequenceBase* seqp)
{
  double rf_period=0.0;
  if (seqp) {
    rf_period=seqp->duration();
    if (seqp->isCW) {
      if (detect_period<rf_period+tolerance)
	rf_period=detect_period;
      else
	shortCW.raise(seqp->name(),true);
    }
  }
  return rf_period;
}

void flush_transitionlog()
{
  if (transitions_fp && (transitions_fp!=stdout))
    fclose(transitions_fp);    
}

void ensure_vararray()
{
  if (!!vararray)
    return;

  ensure_array(); //!< can be called twice
  
  if (array_n0>1) {
    size_t usevars=varpars.size();
    if (powder_array)
      usevars+=3;
    if (usevars) {//NB array parameters within sum are not included
      vararray.create(array_n0,usevars,0.0);
      varnames.create(usevars);
      size_t i=0;
      for (;i<varpars.size();i++)
	varnames(i)=varpars(i).name();
      if (powder_array) {
	varnames(i++)="alpha";
	varnames(i++)="beta";
	varnames(i)="gamma";
      }
    }
  }
}

void save_parameters()
{
  ensure_vararray();
  if (!(vararray.empty())) {
    if (row_index>=vararray.rows())
      throw BadIndex("save_parameters",row_index,vararray.rows());
    get_variables(vararray.row(row_index));      
  }
}

void parse_n(int);

struct SysVar_histogram : public SystemVariable<double*> {
  SysVar_histogram(const std::string& name_, double* value_)
    : SystemVariable<double*>(name_,value_) {}
  void update() { update_lineshapegen(); }
};

static double curscale=0.0; //!< current powder weighting

SystemVariable<double*> v_alpha("alpha",&(global_powder.alpha),rad_to_deg,V_UPDATEINTS | V_POWDERSUMDEP);
SystemVariable<double*> v_beta("beta",&(global_powder.beta),rad_to_deg,V_UPDATEINTS | V_POWDERSUMDEP);
SystemVariable<double*> v_gamma("gamma",&(global_powder.gamma),rad_to_deg,V_UPDATEINTS | V_POWDERSUMDEP);
SystemVariable<double*> v_weight("weight_orientation",&curscale,1.0,V_ISFIXED | V_POWDERSUMDEP);
SysVar_histogram v_histlw("histogram_linewidth",&hist_lw);
SysVar_histogram v_histgfrac("histogram_gaussianfraction",&hist_gfrac);
SystemVariable<int*> v_powderquality("powderquality",&nzcw,V_ISFIXED | V_ISCONST);
SystemVariable<int*> v_orientation("i_orientation",&orientation,V_ISFIXED | V_POWDERSUMDEP);
//SystemVariable<int*> v_var_index("i_array",&var_index,true,false);
//SystemVariable<int*> v_row_index("i_row",&row_index,true,false);

complex single_point(const BlockedOperator& sigma)
{
  return (!detect || detectED) ? NMR_trace(sigma) : NMR_trace(detect,sigma);
}
    
static PowderMethod::include_t offsetqual=PowderMethod::middle;
bool havepowderquality() { return (nzcw>=2); }

enum { CRYS_SPHERE =1, //!< complete sphere
       CRYS_HEMISPHERE =2, //!< hemisphere
       CRYS_OCTANT =4, //!< octant
       CRYS_PRINT =128, //!< verbose
       CRYS_START =256, //!< include pole
       CRYS_END =512, //!< include max beta angle
       CRYS_MIDDLE =1024, //!< between limits
       CRYS_BOTH =2048, //!< include both limits
       CRYS_REVERSE =8192 //!< reverse order
};

bool ammaster=true;

ThreadWarning<> parfailed_warning("Not using parallelisation as fewer summation steps than processors",&NMRsim_once_warning); 
ThreadWarning<> parnone_warning("Ignoring -enable:parallel as no parallelisation mechanism enabled at compile time",&NMRsim_once_warning); 
ThreadWarning<> npfailed_warning("Cannot determine number of processor cores - use NMRSIM_NUM_CORES environment to set explicitly",&NMRsim_once_warning); 
ThreadWarning<> numcoreslarge_warning("Number of threads set using NMRSIM_NUM_CORES exceeds number of (online) processor cores",&NMRsim_once_warning); 

bool cleanup_parallel()
{
#ifdef USE_FORK_CONTROLLER
  if (!!parconp) {
    try {
      parconp->join_kill();
    } catch (...) { //!< don't allow exceptions to leak
      return false;
    }
  }
#endif
  return true;
}

void init_parallel(int& argc, char**& argv)
{
  if (!optparallel)
    return;
#if (defined(USEMPI) || defined(HAVE_FORK_CONTROLLER)) && defined(USEMKL)
  mkl_set_num_threads(1); //!< disable MKL threading to prevent oversubscription
#endif
  const int lverbose=(verbose & (VER_POWDER | VER_PARALLEL)) ? verbose_level : 0;
#ifdef USEMPI
  parconp.reset(new MPI_controller(argc,argv,lverbose));
  global_workers=parconp->get_workers();
  optparallel.setusage(true);
#else
#ifdef USE_FORK_CONTROLLER
  size_t nproc=0;
  const char* envval=register_and_getenv("NMRSIM_NUM_CORES");
  if (envval) {
    nproc=parse_counting(envval);
    if (nproc==0)
      error_abort("NMRSIM_NUM_CORES does not evaluate to a valid number of processors");
  }
  else {
    if (optparallel.isauto())
      assumeparallel=false; //!< don't assume parallel mode
    else {
      if (optparallel.isenabled() && lverbose)
	std::cout << "NMRSIM_NUM_CORES unset and so defaulting to using all available processor cores\n";
    }
  }
#ifdef _SC_NPROCESSORS_ONLN
  const size_t nonln=sysconf(_SC_NPROCESSORS_ONLN);
  if (nproc) {
    if (!nochecks && (nproc>nonln))
      numcoreslarge_warning.raise();
  }
  else {
    if (assumeparallel)
      nproc=nonln;
  }
  if (nproc)
    parconp.reset(new Fork_controller(nproc,lverbose));
#endif
  global_workers=nproc;
  optparallel.setusage(global_workers!=0);
  if (global_workers==0) {
    if (optparallel.isenabled())
      npfailed_warning.raise();
    return;
  }
#else
  optparallel.setusage(false);
  if (optparallel.isenabled())
    parnone_warning.raise();
  return;
#endif
#endif
}

void ensure_parallel()
{
  static bool done_parallel_init=false;
  if (done_parallel_init)
    return;
  done_parallel_init=true;
#ifdef HAVE_PARSYS
  if (!parconp)
    throw InternalError("ensure_parallel");
  ammaster=parconp->ammaster();
  //      disablelogging=!ammaster;
  if (!silent && ammaster)
    std::cout << "Starting parallel system with " << global_workers << " workers." << (optparallel.isauto() ? " Disable parallel execution with -disable:parallel\n" : "\n");
//   else {
//     parfailed_warning.raise();
//     parconp.clear();
//   }
#else
  throw InternalError("ensure_parallel");
#endif
}

std::ostream& parser_printthread(std::ostream& ostr)
{
#ifdef HAVE_PARSYS
#ifdef USE_FORK_CONTROLLER
  static const int offset=1;
#else
  static const int offset=0; //!< don't shift as master process is 0 and workers are numbered from 1
#endif
  if (!!parconp)
    ostr << 'T' << (offset+parconp->get_thread_num()) << ": ";
#endif  
  return ostr;
}

static void make_nD_par_variables()
{
  add_systemvarmap(v_ni);
  add_systemvarmap(v_sw1);

  command_Factory_t& par_Factory(get_par_Factory());
  char vname[6];
  for (size_t i=1;i<=MAX_DIMENSIONS;i++) {
    sprintf(vname,"n%" LCM_PRI_SIZE_T_MODIFIER "u",i);
    par_Factory[vname]=par_t(&parse_n,i,true);
    sprintf(vname,"sw%" LCM_PRI_SIZE_T_MODIFIER "u",i);
    par_Factory[vname]=par_t(&parse_swn,i,true);
  }
  par_Factory["ni"]=par_t(&parse_n,1,true);
}

void add_reset_warning(BaseWarning& warning)
{
  warning.type(BaseWarning::FirstOnly);
  resetwarning_list.push_back(&warning);
}

void reset_warnings()
{
  for (size_t i=resetwarning_list.size();i--;)
    resetwarning_list(i)->reset();
}

static void prepar_callback()
{
  add_systemvarmap(v_powderquality);  
  add_systemvarmap(v_alpha);
  add_systemvarmap(v_beta);
  add_systemvarmap(v_gamma);
  add_systemvarmap(v_weight);
  add_systemvarmap(v_orientation);

  make_nD_par_variables();

  command_Factory_t& par_Factory(get_par_Factory());
  par_Factory["powderquality"]=&parse_powderquality;
  par_Factory["histogram"]=par_t(&parse_histogram,ACC_FREQ,true);
  par_Factory["pseudohistogram"]=par_t(&parse_histogram,ACC_PSEUDO,true);
  par_Factory["method"]=par_t(&parse_ignored,true);
  par_Factory["crystal_file"]=&parse_crystal_file;

  add_flushcallback(flush_transitionlog); //!< flush transition log when calculation terminates

  add_reset_warning(chebyshev_convergence_warning);
  add_reset_warning(chebyshev_nonunitary_warning);
  add_reset_warning(cmatrix_eigensystem_controller.nonorthogonal_warning);
}

namespace {
  struct Proxy_ {
    Proxy_() {
      //      optional_map_t& optional_map(get_optional_map());
      add_option(optgamma);
      add_option(optED);
      //      optional_map["parallel"]=&optparallel;
      add_option(optforceeigenbasis);
      //      optional_map["longdtsync"]=&optlongdtsync;
      add_option(optforcepointbypoint);

      addprepar_callback(prepar_callback);
    }
  };  
  static Proxy_ proxy;
}

#define NMRSIM_ENABLE_ROMAN

class DotMaker {
public:
  DotMaker(double dottimev =dotlimit)
    : makedots(dottimev>=0.0), donedots(false), dottime(dottimev), unhandled(0) {}
  void flush();  
  void print();

private:
  bool makedots,donedots;
  timer<> stopwatch;
  double dottime;
  size_t unhandled;
  size_t basestep;
  char basechar;
};

void DotMaker::flush()
{
  if (!donedots)
    return;  
#ifdef NMRSIM_ENABLE_ROMAN  
  if (unhandled>=500) {
    std::cout << 'D';
    unhandled-=500;
  }
  while (unhandled>=100) {
    std::cout << 'C';
    unhandled-=100;
  }  
  if (unhandled>=50) {
    std::cout << 'L';
    unhandled-=50;
  }
  while (unhandled>=10) {
    std::cout << 'X';
    unhandled-=10;
  }
  if (unhandled>5) {
    std::cout << 'V';
    unhandled-=5;
  }
  for (size_t i=unhandled;i--;)
    std::cout << 'I';
#endif
  std::cout << '\n';
}

void DotMaker::print() 
{
  if (!makedots)
    return;       
  if (!donedots) {
    const double t=stopwatch();
    if (t>dottime) {
      donedots=true;      
#ifdef NMRSIM_ENABLE_ROMAN
      const double rate=unhandled/t;      
      if (rate>=1000.0) {
	basestep=1000;
	basechar='M';
      }
      else {
	if (rate>=100.0) {
	  basestep=100;
	  basechar='C';
	}
	else {
	  if (rate>=10.0) {
	    basestep=10;
	    basechar='X';
	  }
	  else {
	    basestep=1;
	    basechar='I';
	  }
	}
      }
#else
      basechar='.';
      basestep=1;
#endif
      while (unhandled>=basestep) {
	std::cout << basechar;
	unhandled-=basestep;
      }
      std::cout.flush();
    }
    else
      unhandled++;
  }
  else {
    unhandled++;
    if (unhandled==basestep) {	
      std::cout << basechar;
      std::cout.flush();
      unhandled-=basestep;
    }
  }
}    

static void set_orientation(Euler& powder, double& scale, size_t orient, int locverbose)
{
  powder.gamma=gamma_zero;
  if (!!powdm) {
    if (powder_reverse)
      orient = powdm->orientations()-orient-1;    
    powdm->orientation(powder,scale,orient);
    if (scale<0.0) {
      parser_printthread(std::cerr) << "Powder orientation " << orient << " returned -ve scale factor: " << powder << ' ' << scale << '\n';
      error_abort();
    }
    if (powdm->angles()==3)
      powder.gamma+=gamma_zero; //for 3 angle sets, use gamma_zero as offset
  }
  orientation=orient; //set global orientation index
  if (locverbose & VER_POWDER) {
    std::cout << "Orientation: " << powder;
    if (verbose_level>1)
      std::cout << " (" << scale << ')';
    std::cout << std::endl;
  }
  global_powder=powder;
}

size_t nobs=0;
double hist_include_threshold=0.0;
double hist_log_threshold=0.0;

void parse_histogram_range()
{
  //      if (hist_flags & HIST_RANGE) {
  hist_min=parse_double();
  hist_max=parse_double();
  if (hist_min==hist_max)
    error_abort("Bad histogram range: histogram -range <min> <max>");
  if (hist_min>hist_max)
    ::std::swap(hist_min,hist_max);
}

ContextWarning<> doubleascii_warning("-double has no effect on ASCII output",&NMRsim_once_warning);
ContextWarning<> streamopenfailed_warning("failed to open transitions stream",&NMRsim_repeat_warning);
ContextWarning<> logfileMPI_warning("histogram log_file is not compatible with MPI operation - any file will almost certainly be garbled",&NMRsim_repeat_warning);

void parse_histogram_log_file()
{
#ifdef USEMPI
  if (global_workers)
    logfileMPI_warning.raise();
#endif
  const char* fname=parse_string();
  if (parser_isnormal())
    hist_log_threshold=parse_double();

  static flagsmap_type histflags;

  if (histflags.empty()) {
    histflags["binary"]=HIST_BINARY;
    histflags["double"]=HIST_DOUBLE;
    histflags["real"]=HIST_REAL;
    histflags["append"]=HIST_APPEND;
  }
  static const int mask=HIST_BINARY | HIST_DOUBLE | HIST_REAL | HIST_APPEND;
  hist_flags=(hist_flags & ~mask) | parse_flags(histflags);

  //      if ((hist_flags & HIST_BINARY) && !fname)
  //	parser_printcontext() << "Warning: -binary specified without file\n";
  if ((hist_flags & HIST_DOUBLE) && !(hist_flags & HIST_BINARY))
    doubleascii_warning.raise();

  if (strcmp(fname,"-")==0) {
    transitions_fp=stdout;
    hist_flags &= ~HIST_BINARY; //ensure ASCII;
  }
  else {
    if (!nochecks)
      checklogfile(fname,logfile_active(),hist_flags & HIST_APPEND);
    static char openflags[3]={'\0'};
    openflags[0]=(hist_flags & HIST_APPEND) ? 'a' : 'w';
    openflags[1]=(hist_flags & HIST_BINARY) ? 'b' : '\0';
    transitions_fp=fopen(fname,openflags);
    if (!transitions_fp)
      streamopenfailed_warning.raise();
  }
  //else {
  //  if (hist_flags & HIST_REAL)
  //    parser_printcontext() << "Warning: real flag has no effect without output file\n";
  //}
}

bool transitions_log_active() { return (transitions_fp!=NMRSIM_NULL); }

ContextWarning<> excessivecutoff_warning("histogram lineshape cutoff is more than " NMRSIM_STRINGISE(CUTOFF_TRIGGER) "%",&NMRsim_once_warning);

int makeshapeflags(bool dofold)
{
  int flags=LineshapeSpec::variable;
  if (!dofold)
    flags|=LineshapeSpec::nofold;
  if ((verbose & VER_GEN) && (verbose_level>1))
    flags|=LineshapeSpec::verbose;  
  return flags;
}

void parse_histogram_lineshape()
{
  if (acc_mode!=ACC_FREQ)
    error_abort("can only specify lineshape in histogram mode (use addlb for time-domain data)");

  const char histogram_lineshape_syntax[]="histogram lineshape <linewidth / Hz>#<Guassian fn>#[<lineshape steps>#<cutoff>]";
  
  lineshape_spec.cutoff=NMRSIM_DEFAULT_LINESHAPE_CUTOFF;
  lineshape_spec.resolution_steps=NMRSIM_DEFAULT_LINESHAPE_STEPS;
  lineshape_spec.cache_maximum= !optcache ? 0 : NMRSIM_DEFAULT_LINESHAPES_CACHE;
  parse_system_variable_syntax(histogram_lineshape_syntax,1,v_histlw);
  if (are_left()) {
    parse_system_variable_syntax(histogram_lineshape_syntax,2,v_histgfrac);
    if (are_left()) {
      lineshape_spec.resolution_steps=parse_unsigned_syntax(histogram_lineshape_syntax,3);
      if (lineshape_spec.resolution_steps==0)
	error_abort("histogram steps must be >0");
    
      lineshape_spec.cutoff=parse_double_syntax(histogram_lineshape_syntax,4);
      if ((lineshape_spec.cutoff>=1.0) || (lineshape_spec.cutoff<0.0))
	error_abort("histogram lineshape cutoff must be >=0 and <1");
    
      if (lineshape_spec.cutoff>=0.01*CUTOFF_TRIGGER)
	excessivecutoff_warning.raise();
    }
  }
  lineshape_spec.flags=makeshapeflags(hist_flags & HIST_FOLD);
  update_lineshapegen();
}

void parse_histogram_flags()
{
  if (parser_isnormal())
    hist_include_threshold=parse_double();
    
  static flagsmap_type histflags;
  if (histflags.empty()) {
    //  histflags["binary"]=HIST_BINARY;
    //histflags["interpolate"]=HIST_INTERP;
    //histflags["double"]=HIST_DOUBLE;
    histflags["fold"]=HIST_FOLD;
    //histflags["real"]=HIST_REAL;
    //histflags["range"]=HIST_RANGE;
  }  
  static const int mask=HIST_FOLD;
  if (are_left()) {
    if (lineshape_genp)
      error_abort("can't change flags after histogram lineshape"); //!< fold flag is needed when generator created
    hist_flags=(hist_flags & ~mask) | parse_flags(histflags);
  }
          
//   if ((hist_flags & HIST_BINARY) && !fname)
//     parser_printcontext() << "Warning: -binary specified without file\n";
//   if ((hist_flags & HIST_DOUBLE) && !(hist_flags & HIST_BINARY))
//     parser_printcontext() << "Warning: -double has no effect on ASCII output\n";
  
//       if (fname) {
// 	if (strcmp(fname,"-")==0) {
// 	  transitions_fp=stdout;
// 	  hist_flags &= ~HIST_BINARY; //ensure ASCII;
// 	}
// 	else {
// 	  transitions_fp=fopen(fname,(hist_flags & HIST_BINARY) ? "wb" : "w");
// 	  if (!transitions_fp)
// 	    parser_printcontext() << "Warning: failed to open transitions stream\n";
// 	}
//       }
//       else {
// 	if (hist_flags & HIST_REAL)
// 	  parser_printcontext() << "Warning: real flag has no effect without output file\n";
//       }
//     }
//   }
}

//command_Factory_t& get_histogram_Factory()
//{
//  return histogram_Factory;
//}

void parse_histogram(int newm)
{  
  const acc_t newaccm=static_cast<acc_t>(newm);
  if (acc_mode==ACC_TIME)
    acc_mode=newaccm;
  else {
    if (acc_mode!=newaccm)
      error_abort("can't switch between acquisition modes");
  }

  //command_Factory_t& histogram_Factory(get_histogram_Factory());
  //static bool doneadd=false;

  static command_Factory_t histogram_Factory;
  if (histogram_Factory.empty()) {
    histogram_Factory["range"]=&parse_histogram_range;
    histogram_Factory["log_file"]=&parse_histogram_log_file;
    histogram_Factory["lineshape"]=&parse_histogram_lineshape;
  }

  if (are_left()) {
    const char* keyname=parse_string(F_IGNORESYNTAX);
    //    const command_Factory_t::iterator iter(histogram_Factory.find(keyname));
    // if (iter!=histogram_Factory.end()) {
    //   (iter->second)();
    //   return;
    // }
    (factory_parse(keyname,histogram_Factory))();
    return;
  }
  parse_histogram_flags();
}

template<class T> void write_binary(const complex& amp, double freq)
{
  T binvals[3];
  if (hist_flags & HIST_REAL) {
    binvals[0]=real(amp);
    binvals[1]=freq;
    fwrite(binvals,sizeof(T),2,transitions_fp); 
  }
  else {
    binvals[0]=real(amp);
    binvals[1]=imag(amp);
    binvals[2]=freq;
    fwrite(binvals,sizeof(T),3,transitions_fp); 
  }
}

void raw_write_transitions(const complex& amp, double freq)
{
  if (hist_flags & HIST_BINARY) {
    if (hist_flags & HIST_DOUBLE)
      write_binary<double>(amp,freq);
    else
      write_binary<float>(amp,freq);
  }
  else {
    if (hist_flags & HIST_REAL)
      fprintf(transitions_fp,"%g %g\n",real(amp),freq);
    else
      fprintf(transitions_fp,"%g %g %g\n",real(amp),imag(amp),freq);
  }
}

void add_transition(BaseHistogram<complex>& hist,const complex& amp, double freq)
{
//   if (verbose_level>1)
//     std::cout << "Adding " << amp << " of " << freq << '\n';
  static const double threshold2=hist_include_threshold*hist_include_threshold;
  if (norm(amp)<threshold2)
    return;

  if (acc_mode==ACC_PSEUDO) {
    if (sw==0.0)
      throw InternalError("add_transition");
    const BaseList<complex> data(hist.row());
    add_FID(data,amp,-freq/sw); //!< flip sign for consistency between pseudohist & TD
  }
  else
    hist.add(amp,freq);

  static const double logthreshold2=hist_log_threshold*hist_log_threshold;
  if (transitions_fp && (norm(amp)>=logthreshold2)) {
    if (hist_flags & HIST_FOLD) {
      const double minf=hist.minimum();
      freq=minf+mod(freq-minf,hist.range());
    }
    if ((hist_max==hist_min) || ((freq<=hist_max) && (freq>=hist_min))) {
#ifdef USE_FORK_CONTROLLER
      if (global_workers) {
	static Mutex<> lock;
	MutexLock<> guard(lock);
	raw_write_transitions(amp,freq);
      }
      else
#endif
	raw_write_transitions(amp,freq);
    }
  }
}

template<typename T> void add_spectrum_(BaseHistogram<complex>& hist, const complex& scale, SpectrumIterator<T>& Obj, double fold_sw =0.0)
{
  const double fold_sw2=fold_sw/2.0;
  const bool check_fold(fold_sw>=0.0);

  T amp;
  double freq;
  while (Obj(amp,freq)) {
    if (check_fold) {
      if (freq>=fold_sw2) //fold frequencies into +/- fold_sw/2
	freq-=fold_sw;
      else {
	if (freq<-fold_sw2)
	  freq+=fold_sw;
      }
    }
    //if (isreal)	
    //  add_transition(hist,amp*rscale,freq);
    //else
      add_transition(hist,amp*scale,freq);
  }
}

template<typename F> void add_spectrum(BaseHistogram<complex>& hist, const complex& scale, F& Obj, double fold_sw =0.0)
{
  if (eigenvalue!=NMRSIM_INVALID)
    Obj.eigenvalue(eigenvalue);
  add_spectrum_(hist,scale,Obj,fold_sw);
}

template<class Obj,class Ugen> void add_metaspectrum(BaseHistogram<complex>& hist, const complex& scale, double fold_sw,Obj& obj, const Ugen& propgen, size_t npts, double t1, double t2, const BlockedOperator& ldensity, const BlockedOperator& ldetect, int lverbose)
{
  if (sideband!=NMRSIM_INVALID)
    obj.sideband(sideband);
  RawMetaSpectrum<Ugen,Obj,BlockedOperator> mobj(obj,propgen,npts,t1,t2,ldensity,ldetect,lverbose);
  add_spectrum(hist,scale,mobj,fold_sw);
}

ThreadWarning<> ignoringsideband_warning("sideband ignored in static/time-domain simulation",&NMRsim_once_warning);

template<class Obj,class Ugen> void add_metaspectrum(BaseHistogram<complex>& hist, const complex& scale, Obj& obj, const Ugen& propgen, double deltat, const BlockedOperator& ldensity, const BlockedOperator& ldetect, int lverbose)
{
  if (sideband!=NMRSIM_INVALID)
    ignoringsideband_warning.raise();

  RawMetaSpectrum<Ugen,Obj,BlockedOperator> mobj(obj,propgen,1,master_time,master_time+deltat,ldensity,ldetect,lverbose);
  add_spectrum(hist,scale,mobj);
}

namespace {
  template<class Holder,class T> void ensureobj_(Holder& holder,size_t,Type2Type<T>,Bool2Type<false>)
  {
    holder.set(Type2Type< smartptr<T,false> >()).reset(new T((verbose & VER_GEN) ? verbose_level : 0));
  }

  template<class Holder,class T> void ensureobj_(Holder& holder,size_t lnobs,Type2Type<T>,Bool2Type<true>)
  {
    holder.set(Type2Type< smartptr<T,false> >()).reset(new T(lnobs,(verbose & VER_GEN) ? verbose_level : 0));
  }

  template<class T> bool matchobs(const T& obj,Bool2Type<true>) { return (obj.observations()==nobs); }
  template<class T> bool matchobs(const T&, Bool2Type<false>) { return true; }
}

void checkED()
{
  const bool ison=optED();
  if (ison) {
    density.makeidentitycache();
    if (use_detectp)
      use_detectp->makeidentitycache();
  }
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Checking for sigma0/detect matching identity operator: " << (ison ? "Yes\n" : "No\n");
}

template<class T> T& MasterObj::ensuretobj()
{
  typedef Type2Type< smartptr<T,false> > switch_type;
  if (!(objholder.istype(switch_type())) || !matchobs(*(objholder.get(switch_type())),Bool2Type<SignalGenerator_traits<T>::gamma>()) )
    ensureobj_(objholder,nobs,Type2Type<T>(),Bool2Type<SignalGenerator_traits<T>::gamma>());
  use_detectp = (detectED & SignalGenerator_traits<T>::allowED) ? &empty_op : &detect;
  checkED();
  return *(objholder.get(switch_type()));
}

template<class PropGen> void MasterObj::add_freq(BaseHistogram<complex>& hist, complex scale, PropGen& propgen)
{
  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  const double fold_sw=nobs/detect_period;

  if (use_gamma_angles) {
    if (detectED)
      add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaPeriodicSpectrumED>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,empty_op,lverbose);
    else
      add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaPeriodicSpectrum>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,detect,lverbose);
  }
  else
    add_metaspectrum(hist,scale,fold_sw,ensurefobj<PeriodicSpectrum>(),propgen,nobs,master_time,master_time+detect_period,density,detect,lverbose);
}

template<class StaticObj,class SpinObj> void MasterObj::add_freq(BaseHistogram<complex>& hist, complex scale, const StaticObj& Hstaticp, const SpinObj& Hspinp, const ActionDirectAcq& acq)
{
  //  assert(acqp!=NMRSIM_NULL);
  const CompSequenceBase* acqseqp = acq.sequence();

  if (Hempty() && !acqseqp) {
    const complex val=scale*single_point(density);
    add_transition(hist,val,0.0);
    return;
  }
  const bool point_by_point=acq.ispoint_by_point();
  if (point_by_point || !detect_period || !nobs)
    throw InternalError("add_freq");

  const double dt=dwelltime();
  const int lverbose = (verbose & VER_GEN) ? verbose_level : 0;
  const double fold_sw=nobs/detect_period;
  //  const double act_acq_period=acqseqp ? acqseqp->duration() : 0.0;
  const int prop_flags= acqseqp ? acqseqp->prop_flags() : global_prop_flags;
  const double rf_period=get_rf_period(acqseqp);
  const double maxdt=get_maxdt();

  switch (hamobjtype) {
  case M_SPIND: {
    MetaPropagator propgen(*Hspindp,0.0,lverbose,prop_flags);
    set_partitioning(propgen);
    if (use_gamma_angles) {
      if (detectED)
	add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaInhomogeneousSpectrumED>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,empty_op,lverbose);
      else
	add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaInhomogeneousSpectrum>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,detect,lverbose);
    }
    else
      add_metaspectrum(hist,scale,fold_sw,ensurefobj<InhomogeneousSpectrum>(),propgen,nobs,master_time,master_time+detect_period,density,detect,lverbose);
    break;
  }
  case M_SPIN:
    if (acqseqp) {
      if (acqseqp->modulation().get()) {
	PhaseModulatedPropagator& propgen(const_cast<CompSequenceBase*>(acqseqp)->updatepm(*this));
	add_freq(hist,scale,propgen);
	//	add_metaspectrum(hist,scale,fold_sw,ensurefobj<PeriodicSpectrum>(),propgen,nobs,master_time,master_time+detect_period,density,detect,lverbose);
      }
      else {
	SequencePropagator propgen(*Hspinp,maxdt,acqseqp->seqlist,master_time,rf_period,tolerance,lverbose,prop_flags);
	set_partitioning(propgen,*acqseqp);
	propgen.synchronisation_hint(acqseqp->synchronisation_time());
	add_freq(hist,scale,propgen);
      }
    }
    else {
      MetaPropagator propgen(*Hspinp,maxdt,lverbose,prop_flags);
      set_partitioning(propgen);
      if (use_gamma_angles) {
	if (detectED)
	  add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaPeriodicSpectrumED>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,empty_op,lverbose);
	else
	  add_metaspectrum(hist,scale,fold_sw,ensurefobj<GammaPeriodicSpectrum>(),propgen,use_gamma_angles,master_time,master_time+detect_period,density,detect,lverbose);
      }
      else
	add_metaspectrum(hist,scale,fold_sw,ensurefobj<PeriodicSpectrum>(),propgen,nobs,master_time,master_time+detect_period,density,detect,lverbose);
    }
    break;
    
  case M_STATIC:
    if (acqseqp) {
      SequencePropagator propgen(*Hstaticp,acqseqp->seqlist,master_time,rf_period,tolerance,lverbose,prop_flags);
      set_partitioning(propgen,*acqseqp);
      if (nobs>1)
	add_metaspectrum(hist,scale,fold_sw,ensurefobj<PeriodicSpectrum>(),propgen,nobs,master_time,master_time+detect_period,density,detect,lverbose);
      else {
	if (detectED)
	  add_metaspectrum(hist,scale,ensurefobj<StaticSpectrumED>(),propgen,dt,density,empty_op,lverbose);
	else
	  add_metaspectrum(hist,scale,ensurefobj<StaticSpectrum>(),propgen,dt,density,detect,lverbose);
      }
    }
    else {
      if (detectED)
	add_metaspectrum(hist,scale,ensurefobj<StaticSpectrumED>(),*Hstaticp,dt,density,empty_op,lverbose);
      else
	add_metaspectrum(hist,scale,ensurefobj<StaticSpectrum>(),*Hstaticp,dt,density,detect,lverbose);
    }
    break;
    
  case M_STATICD:
    if (detectED)
      add_metaspectrum(hist,scale,ensurefobj<StaticSpectrumED>(),*Hstaticdp,dt,density,empty_op,lverbose);
    else
      add_metaspectrum(hist,scale,ensurefobj<StaticSpectrum>(),*Hstaticdp,dt,density,detect,lverbose);
    break;
    
  default:
    throw InternalError("Unknown Hamiltonian type");
  } //hamobjtype 
}

template<class T> T& MasterObj::ensurefobj()
{
  typedef Type2Type< smartptr<T,false> > switch_type;
  if (!(fobjholder.istype(switch_type())) || !matchobs(*(fobjholder.get(switch_type())),Bool2Type<SignalGenerator_traits<T>::gamma>()) )
    ensureobj_(fobjholder,nobs,Type2Type<T>(),Bool2Type<SignalGenerator_traits<T>::gamma>());
  checkED();
  return *(fobjholder.get(switch_type()));
}


ThreadWarning<> suggestfold_warning("frequency domain MAS propagation but without folding of spectral histogram\nConsider adding histogram method -fold to par",&NMRsim_once_warning);

void MasterObj::add_freq(BaseList<complex> FIDSpec, complex cscale, bool liscomplex, const ActionDirectAcq& acq)
{
  HistogramMaker hist(FIDSpec,sw);
  if (!(hist_flags & HIST_FOLD) && spin_rate)
    suggestfold_warning.raise();
  add_freq(hist(),cscale,liscomplex,acq);
}

void MasterObj::add_freq(BaseHistogram<complex>& Spec, complex cscale, bool iscomplex, const ActionDirectAcq& acq)
{
  if (iscomplex)
    add_freq(Spec,cscale,Hstaticp.get_complex(),Hspinp.get_complex(),acq);
  else
    add_freq(Spec,cscale,Hstaticp.get_real(),Hspinp.get_real(),acq);
}

ThreadWarning<> changedsw_warning("spectral width for row appears to have changed",&NMRsim_repeat_warning);

ActionDirectAcq::create_status_t ActionDirectAcq::create_status=ActionDirectAcq::NONE;

ThreadWarning<> ActionDirectAcq::empty_warning("acquisition with zero data points to acquire",&NMRsim_once_warning);

void MasterObj::raw_add(const ActionAcqPoint& acq, BlockedOperator& sigma)
{
  if (transients_active()) {
    ActionAcq::activetransients_warning.raise();
    flush_transients();
  }    
  const complex cscale=(curscale*globalscale)*expi(acq.recphase());

  DataStore& Specall(*curstorep);
  init_row(Specall);

  const size_t nstart=ActionDirectAcq::acquisition_count();
  if (nstart>=np)
    error_abort("acqpoint failed: FID already fully acquired");

  BaseList<complex> fullFIDSpec(Specall.row(row_index));
  const complex res=cscale*single_point(sigma);
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Adding to data point " << (nstart+1) << ": " << res << '\n';
  fullFIDSpec(nstart)+=res;
  ActionDirectAcq::add_acquisition_count(1U);
  if (nstart==np-1)
    sigma.clear();
}

ThreadWarning<> freq_warning("frequency domain propagation has not been actively maintained - use with caution! Disable warning with -nochecks.",&NMRsim_once_warning);

void MasterObj::raw_add(const ActionDirectAcq& acq, BlockedOperator& sigma, double& t, int toacq)
{
  const double dt=dwelltime();
  if (verbose & VER_GEN) {
    std::cout << "Acquisition starts at t=";
    prettyprint_time(master_time) << '\n';
  }
  if (toacq==0)
    return;

  if (insidecontrol)
    error_abort("acq cannot be inside control / loop structure - did you mean to use acqn or acqpoint?");

  if (transients_active()) {
    ActionAcq::activetransients_warning.raise();
    flush_transients();
  }    
  
  const complex cscale=(curscale*globalscale)*expi(acq.recphase());

  DataStore& Specall(*curstorep);
  init_row(Specall);

  const bool liscomplex=iscomplex();

  BaseList<complex> fullFIDSpec(Specall.row(row_index));
  const size_t nstart=ActionDirectAcq::acquisition_count();
  BaseList<complex> FIDSpec;
  size_t nend;
  if (toacq<0)
    nend=np-1;
  else {
    nend=nstart+toacq-1;    
    if (nend>=np) {
      parser_printthread(std::cerr) << "acquiring " << toacq << " points would overflow data data - too many points specified in acqn\n";
      error_abort();
    }
  } 
  const size_t npts=nend-nstart+1;
  if (npts==0)
    ActionDirectAcq::empty_warning.raise();
  else {
    if (nend<nstart)
      throw InternalError("raw_add");
    BaseList<complex> FIDSpec(fullFIDSpec(range(nstart,nend)));
    if (acc_mode!=ACC_TIME) {
      if (isfrequency() && (fullFIDSpec.size()!=FIDSpec.size()))
	throw InternalError("frequency domain propagation combined with partial acquisition");
      if (!nochecks)
	freq_warning.raise();
      add_freq(FIDSpec,cscale,liscomplex,acq);
    }
    else {
      if (liscomplex)
	add_td(FIDSpec,cscale,Hstaticp.get_complex(),Hspinp.get_complex(),acq);
      else
	add_td(FIDSpec,cscale,Hstaticp.get_real(),Hspinp.get_real(),acq);
    }
    
    t+=npts*dt; //!< increment time
    ActionDirectAcq::add_acquisition_count(npts);
  }
  if (nend==np-1)
    sigma.clear(); //!< clear density matrix as no longer updated
}

void ActionDirectAcq::exec(MasterObj* objp, BlockedOperator& sigma, double& t) const
{
  objp->raw_add(*this,sigma,t,-1);
}

template<class Obj> void add_RawMetaFID_(Obj& FIDgen,BaseList<complex> FID, complex scale)
{
  if (eigenvalue!=NMRSIM_INVALID)
    FIDgen.eigenvalue(eigenvalue);
  if (sideband!=NMRSIM_INVALID)
    ignoringsideband_warning.raise();
  //  if (rratio)
  //  FIDgen.reduction_factor(rratio);
  FIDgen.add_FID(FID,scale);
}

template<class Obj,class UGen> void add_RawMetaFID(Obj& obj, BaseList<complex> FID, complex scale, const UGen& propgen, size_t npts, double t1,double t2,const BlockedOperator& ldensity, const BlockedOperator& ldetect,int lverbose)
{
  RawMetaFID<UGen,Obj,BlockedOperator> FIDgen(obj,propgen,npts,t1,t2,ldensity,ldetect,lverbose);  
  add_RawMetaFID_(FIDgen,FID,scale);
}

template<class Obj,class UGen> void add_RawMetaFID(Obj& obj, BaseList<complex> FID, complex scale, const UGen& propgen, double dt,const BlockedOperator& ldensity, const BlockedOperator& ldetect,int lverbose)
{
  RawMetaFID<UGen,Obj,BlockedOperator> FIDgen(obj,propgen,1,master_time,master_time+dt,ldensity,ldetect,lverbose);
  add_RawMetaFID_(FIDgen,FID,scale);
}

void MasterObj::getgammapars(GammaPeriodicFID& obj, double& use_period, size_t& steps, double period) const
{
  use_period=detect_period;
  steps=use_gamma_angles;
  if (!optsync)
    return;
  const size_t reduction_ratio=check_sync(detect_period/period);
  const bool usereduction=(reduction_ratio>1);
  //  check_optional(optlongdtsync,"exploiting dwell time being multiple of synchronisation period",usereduction,"longdtsync");
  if (reduction_ratio==0)
    return;
  //    throw InternalError("add_td: failed to find sync");
  if (usereduction) {
    if (verbose & VER_GEN)
      std::cout << "FID sampling reduction factor: " << reduction_ratio << '\n';
    obj.reduction_factor(reduction_ratio);
    use_period/=reduction_ratio;
  }
  else {
    if (verbose & VER_GEN)
      std::cout << "FID sampling reduction factor: None\n";
  }
  steps=use_gamma_angles/reduction_ratio;
  if (steps<1)
    throw InternalError("add_td: gamma angles < reduction factor");
}

// template<class PropGen> double getsyncperiod(const PropGen& propgen, Bool2Type<false>)
// {
//   return propgen.period();
// }

// template<class PropGen> double getsyncperiod(const PropGen& propgen, Bool2Type<true>)
// {
//   double period=propgen.synchronisation_hint();
//   if (period==0.0)
//     period=propgen.period();
//   return period;
// }

template<class PropGen> void MasterObj::add_td(BaseList<complex> FID, complex scale, PropGen& propgen, double prop_period)
{
  const int prop_verbose = (verbose & VER_GEN) ? verbose_level : 0;

  if (use_gamma_angles) {
    GammaPeriodicFID& obj(ensuretobj<GammaPeriodicFID>());
    size_t steps;
    double use_period;
    if (prop_period==0.0)
      prop_period=propgen.period();
    if (prop_verbose)
      std::cout << "Propagator generator period: " << (1e6*prop_period) << " us\n";
    getgammapars(obj,use_period,steps,prop_period);
    add_RawMetaFID(obj,FID,scale,propgen,steps,master_time,master_time+use_period,density,*use_detectp,prop_verbose);
  }
  else {
    //    check_optional(optlongdtsync,"exploiting dwell time being multiple of synchronisation period",usereduction,"longdtsync");
    PeriodicFID& obj(ensuretobj<PeriodicFID>());
    add_RawMetaFID(obj,FID,scale,propgen,nobs,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
  }
}

//Warning<> pmpbp_warning("Experimental code combining PhaseModulation and point-by-point (unsynchronised) acquisition triggered. Confirm same result obtained with -disable:phasemodulation.",&NMRsim_once_warning);
ThreadWarning<> pmpbp_warning("Unimplemented potential optimisation - PhaseModulation + point-by-point (unsynchronised)",&NMRsim_once_warning);

template<class StaticObj,class SpinObj> void MasterObj::add_td(BaseList<complex> FID, complex scale, const StaticObj& Hstaticp, const SpinObj& Hspinp, const ActionDirectAcq& acq)
{
  const double dt=dwelltime();
  const CompSequenceBase* acqseqp = acq.sequence();
  const size_t npts(FID.size());
  const int prop_verbose = (verbose & VER_GEN) ? verbose_level : 0;
  if ((Hempty() && !acqseqp) || (npts==1)) {
    const complex val=scale*single_point(density);
    if (prop_verbose>1)
      std::cout << "Adding " << val << " to FID\n";
    FID+=val;
    return;
  }

  //  const double act_acq_period=acqseqp ? acqseqp->duration() : 0.0;
  const int prop_flags= acqseqp ? acqseqp->prop_flags() : global_prop_flags;
  const double rf_period=get_rf_period(acqseqp);
  const double maxdt=get_maxdt();

  if (acq.ispoint_by_point()) {
    typedef Type2Type< smartptr<AsynchronousFID,false> > switcher_t;
    objholder.set(switcher_t()).reset(new AsynchronousFID(master_time,dt,(verbose & VER_GEN) ? verbose_level : 0));    
    AsynchronousFID& obj(*(objholder.get(switcher_t())));
    const double end_time=master_time+(npts-1)/sw;
    switch (hamobjtype) {
    case M_SPIN:
      if (acqseqp) {
	if (acqseqp->modulation().get()) {
	  if (!nochecks) //!< Need to merge MinimalPropagator into PhaseModulatedPropagator
	    pmpbp_warning.raise();
	  //PhaseModulatedPropagator& propgen(const_cast<CompSequenceBase*>(acqseqp)->updatepm(*this));
	  //add_RawMetaFID(obj,FID,scale,propgen,npts-1,master_time,master_time+detect_period,density,detect,prop_verbose);
	}
	//else {
	SequencePropagator propgen(*Hspinp,maxdt,acqseqp->seqlist,master_time,rf_period,tolerance,prop_verbose,prop_flags);
	set_partitioning(propgen,*acqseqp);
	propgen.synchronisation_hint(acqseqp->synchronisation_time());
	add_RawMetaFID(obj,FID,scale,propgen,npts-1,master_time,end_time,density,detect,prop_verbose);
	//}
      }
      else { //!< same as synchronised case
	MetaPropagator propgen(*Hspinp,maxdt,prop_verbose,prop_flags);
	set_partitioning(propgen);
	add_RawMetaFID(obj,FID,scale,propgen,npts-1,master_time,end_time,density,detect,prop_verbose);
      }
      return;
    case M_STATIC:
      if (acqseqp) {
	SequencePropagator propgen(*Hstaticp,acqseqp->seqlist,master_time,rf_period,tolerance,prop_verbose,prop_flags);
	set_partitioning(propgen,*acqseqp);
	add_RawMetaFID(obj,FID,scale,propgen,npts-1,master_time,end_time,density,detect,prop_verbose);
      }
      else {
	StaticFID_H& obj(ensuretobj<StaticFID_H>()); //!< actually just the same as synchronised
	add_RawMetaFID(obj,FID,scale,*Hstaticp,dt,density,*use_detectp,prop_verbose);
      }
      return;
    default:
      break;
    }
    error_abort("Cannot combine diagonal Hamiltonian with unsychronised acquisition. Unexpected combination of -enable/-disable flags or lack of synchronisation? Suppress diagonal Hamiltonian using -disable:realhamiltonian.");    
  }

  if (!detect_period && sw)
    throw InternalError("add_td: detect_period unset");
 
  //not point_by_point
  switch (hamobjtype) {
  case M_SPIND: {
    //const double maxt=1.0/fabs(spin_rate);
    // 	DiagonalInhomogeneousPropagator propgen((*Hspindp)());
    // 	if (use_gamma_angles) {
    // 	  propagators(dUs,gamma_angles,propgen,t,t+maxt);
    // 	  GammaInhomogeneousFID& obj(*objholder(Type2Type<GIFID_t>()));
    // 	  obj.set_Us(dUs);
    // 	  obj.add_FID(FID,scale,density(),detect());
    // 	}
    // 	else {
    // 	  propagators(dUs,nobs,propgen,t,t+maxt);
    // 	  InhomogeneousFID& obj(*objholder(Type2Type<IFID_t>()));
    // 	  obj.set_Us(dUs);	  
    // 	  obj.add_FID(FID,scale,density(),detect());
    // 	}
    if (acqseqp!=NMRSIM_NULL)
      throw InternalError("add_td");
    MetaPropagator propgen(*Hspindp,0.0,prop_verbose,prop_flags);
    set_partitioning(propgen);
    if (use_gamma_angles) {
      GammaInhomogeneousFID& obj(ensuretobj<GammaInhomogeneousFID>());
      add_RawMetaFID(obj,FID,scale,propgen,use_gamma_angles,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
    }
    else {
      InhomogeneousFID& obj(ensuretobj<InhomogeneousFID>());
      add_RawMetaFID(obj,FID,scale,propgen,nobs,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
    }
  }
    break;
    
  case M_SPIN:
    if (acqseqp) {
      if (acqseqp->modulation().get()) {
	PhaseModulatedPropagator& propgen(const_cast<CompSequenceBase*>(acqseqp)->updatepm(*this));
	add_td(FID,scale,propgen,acqseqp->synchronisation_time());
	//PeriodicFID& obj(ensuretobj<PeriodicFID>());	
	//	add_RawMetaFID(obj,FID,scale,propgen,nobs,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
      }
      else {
	SequencePropagator propgen(*Hspinp,maxdt,acqseqp->seqlist,master_time,rf_period,tolerance,prop_verbose,prop_flags);
	set_partitioning(propgen,*acqseqp);
	const double synctime=acqseqp->synchronisation_time();
	propgen.synchronisation_hint(synctime);
	//	if ((verbose & VER_GEN) && (verbose_level>1))
	//  std::cout << propgen;
	add_td(FID,scale,propgen,synctime);
      }
    }
    else {
      MetaPropagator propgen(*Hspinp,maxdt,prop_verbose,prop_flags);
      set_partitioning(propgen);
      add_td(FID,scale,propgen);
//       if (use_gamma_angles) {
// 	GammaPeriodicFID& obj(ensuretobj<GammaPeriodicFID>());
// 	add_RawMetaFID(obj,FID,scale,propgen,use_gamma_angles,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
//       }
//       else {
// 	PeriodicFID& obj(ensuretobj<PeriodicFID>());
// 	add_RawMetaFID(obj,FID,scale,propgen,nobs,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
//       }
    }
    break;
    
  case M_STATIC:
    if (acqseqp) {
      SequencePropagator propgen(*Hstaticp,acqseqp->seqlist,master_time,rf_period,tolerance,prop_verbose,prop_flags);
      set_partitioning(propgen,*acqseqp);
      if (nobs>1) {
	PeriodicFID& obj(ensuretobj<PeriodicFID>());
	add_RawMetaFID(obj,FID,scale,propgen,nobs,master_time,master_time+detect_period,density,*use_detectp,prop_verbose);
      }
      else {
	const bool forcebasis=optforceeigenbasis.isenabled();
	if (verbose & VER_GEN)
	  std::cout << "Force eigenbasis propagation: " << (forcebasis ? "Yes\n" : "No\n");
	if (forcebasis) {
	  StaticFID_H& obj(ensuretobj<StaticFID_H>());
	  add_RawMetaFID(obj,FID,scale,propgen,1,master_time,master_time+dt,density,*use_detectp,prop_verbose);
	}
	else {
	  StaticFID_U& obj(ensuretobj<StaticFID_U>());
	  add_RawMetaFID(obj,FID,scale,propgen,1,master_time,master_time+dt,density,*use_detectp,prop_verbose);
	}
      }
    }
    else {
      StaticFID_H& obj(ensuretobj<StaticFID_H>());
      add_RawMetaFID(obj,FID,scale,*Hstaticp,dt,density,*use_detectp,prop_verbose);
    }
    break;
    
  case M_STATICD: {
    StaticFID_H& obj(ensuretobj<StaticFID_H>());
    add_RawMetaFID(obj,FID,scale,*Hstaticdp,dt,density,*use_detectp,prop_verbose);
  }
    break;
    
  default:
    throw InternalError("Unknown Hamiltonian type");
  } //hamobjtype
}

ThreadWarning<> ActionDirectAcq::frequencydomainRF_warning("frequency domain propagation (histogram) can be numerically unstable in the presence of non-trival RF irradiation.  Check results against equivalent time domain simulation.",&NMRsim_once_warning);
ThreadWarning<> ActionDirectAcq::ignoringsynctime_warning("synchronisation time irrelevant for zero duration sequence (ignored)",&NMRsim_once_warning);
ThreadWarning<> ActionDirectAcq::gammasteps_warning("gamma steps not a multiple of observation steps",&NMRsim_once_warning);
ThreadWarning<> ActionDirectAcq::syncsetnotfound_warning("synchronisation time set, but no sync found",&NMRsim_once_warning);
ThreadWarning<> ActionDirectAcq::syncfailed_warning("could not synchronise acquisition, so calculation may be very slow. Turn on verbose -general for more details or set explicit synchronisation time in acq. Use -nochecks to disable warning.",&NMRsim_once_warning);
ThreadWarning<> ActionDirectAcq::disabledgamma_warning("not enabling gammacompute algorithm due to non-trivial pulseq. Gamma compute can be enabled explicity (if appropriate) using -enable:gammacompute.",&NMRsim_once_warning);

void ActionDirectAcq::checksync()
{
  //determine nature of acquistion
  const double rotor_period=spin_rate ? 1.0/fabs(spin_rate) : 0.0;
  double acq_period=0.0;
  bool gammacompatible=true;
  double use_synctime=synctime;
  if ((use_synctime==0.0) && seqp) //!< if no explicit syncronisation time, use sequence's
    use_synctime=seqp->synchronisation_time();

  if (seqp) {
    if (verbose & VER_GEN)
      std::cout << "CW during acquistion: " << (seqp->isCW ? "Yes\n" : "No\n");

    if (!(seqp->isCW)) {
      gammacompatible=false;
      acq_period=seqp->safe_duration(); //!< may need to rebuild
      if (acq_period) {
	if (acc_mode!=ACC_TIME)
	  frequencydomainRF_warning.raise();
      }
      else {
	if (seqp->synchronisation_time())
	  ignoringsynctime_warning.raise();
	if (spin_rate)
	  error_abort("MAS incompatible with zero duration acquisition sequences");
      }
    }
  }

  const double dt=sw ? 1.0/sw : 0.0;
  if ((np==1) || optforcepointbypoint.isenabled()) {
    detect_period=0.0;
    point_by_point=true;
    optforcepointbypoint.setusage(true);
  }
  else {
    detect_period=use_synctime;
    //if (!detect_period) { // unless set explicitly, use longest time
      if (acq_period>detect_period)
	detect_period=acq_period;
      if (rotor_period>detect_period)
	detect_period=rotor_period;
      if (dt>detect_period)
	detect_period=dt;
      //}    
    point_by_point=false;
  }

  //  size_t nsync_acqseq=0;
  //size_t nsync_rotor=0;
  if (detect_period) {
    if (spin_rate)
      (void)trysync(detect_period,rotor_period,"MAS period");
    
    if (acq_period)
      (void)trysync(detect_period,acq_period,"RF period");
    
    if ((acc_mode==ACC_TIME) && (dt!=0.0))
      trysync(detect_period,dt,"dwell time");

    point_by_point=(detect_period==0.0);  
  }

  //Check whether we can use gamma symmetry
  //Disallow if (non CW) RF during acquisition or non-trivial pulse sequence
  bool allowgamma=(gamma_angles>1) && (spin_rate!=0.0);
  if (!optgamma)
    allowgamma=false; //!< turn off without warning
  else {
    if (allowgamma && optgamma.isauto() && (!gammacompatible || haveactions)) {
      allowgamma=false;
      disabledgamma_warning.raise();
    }
  }
  //  bool allowgamma=(gamma_angles>1) && (spin_rate!=0.0) && ((gammacompatible && !haveactions && optgamma) || (optgamma==OPTIONAL_ON));
  use_gamma_angles=0;

  nobs=0;
  if (!point_by_point) {
    if (acc_mode!=ACC_TIME) {
      nobs=int(0.1+ceil(detect_period/dt));
      if (allowgamma) //round up so also multiple of gamma_steps
	nobs=gamma_angles*int(1+(nobs-1)/gamma_angles);
    }
    else
      nobs=check_sync(detect_period/dt);
    if (nobs==0)
      throw InternalError("checksync");
    optsync.setusage(true);

    if (allowgamma) {
//       if (phasemodp.get()) {
// 	error_abort("Not implemented");
//       }
//       else {
      //! check against number of observations
      if (gamma_angles % nobs) {
	char buf[256];
	snprintf(buf,sizeof(buf),": gamma angles=%i, observations=%" LCM_PRI_SIZE_T_MODIFIER "u",gamma_angles,nobs);
	gammasteps_warning.raise(buf);
      }
      else {
	use_gamma_angles= check_sync(detect_period/rotor_period)*gamma_angles;
	if (use_gamma_angles==0)
	  throw InternalError("checksync");
      }
    }
  }
  else { //!point_by_point
    if (np>1) {
      if (!(optforcepointbypoint.isenabled())) { //!< don't warn if explicitly using forcepointbypoint
	if (use_synctime)
	  syncsetnotfound_warning.raise();
	else {      
	  if (verbose & VER_GEN)
	    parser_printthread(std::cerr) << "Warning: no synchronisation condition found\n";
	  else
	    syncfailed_warning.raise();
	}
      }
      else
	optforcepointbypoint.setusage(true);
      detect_period=np*dt;
    }
    if (acc_mode!=ACC_TIME)
      error_abort("frequency domain calculation is incompatible with asynchronous calculation");
    optsync.setusage(false);
  }
       
  if (verbose & VER_GEN) {
    if (np!=1)
      std::cout << "Point by point propagation: " << (point_by_point ? "Yes\n" : "No\n");
    std::cout << "Overall repeat period: ";
    prettyprint_time(detect_period) << '\n';
    if (nobs)
      std::cout << "Observations per repeat period: " << nobs << '\n';
    if (seqp) {
      std::cout << "Acquisition RF period: ";
      if (acq_period)
	std::cout << "0 - CW?\n";
      else
	prettyprint_time(acq_period) << '\n';
    }
    if (spin_rate && (gamma_angles>1))
      std::cout << "Using gamma compute (gammacompute): " << (use_gamma_angles ? "Yes\n" : "No\n");
  }
  optgamma.setusage(use_gamma_angles!=0);
}

size_t trysync(double& synctime, double period, const char* desc)
{
  if (synctime && optsync()) {
    size_t n=check_sync(synctime/period);
    if (n) {
      if ((verbose & VER_GEN) && (verbose_level>1))
	std::cout << desc << " divides (x" << n << ") into synchronisation time (" << (synctime*1e6) << " us)\n";    
      return n;
    }    
    //! sync failed
    if (verbose & VER_GEN) {
      std::cout << desc << " (";
      prettyprint_time(period) << ") does not divide into sync time (" << (synctime*1e6) << " us)\n";
    }
  }
  synctime=0.0;
  return 0;
}

size_t ActionAcq::checksync(double& sync_period) const
{    
  double rf_period = 0.0;
  const double synctime=seqp ? seqp->synchronisation_time() : 0.0;
  if (seqp && !seqp->isCW) {
    rf_period=seqp->safe_duration();
    if (spin_rate && (rf_period==0.0))
      throw Failed("Zero duration sequences incompatible with MAS");
  }

  const double rotor_period=spin_rate ? 1.0/fabs(spin_rate) : 0.0;
  const double dt=dwelltime();
  sync_period=synctime;
  if (!sync_period) { // unless set explicitly, use longest time
    if (rf_period>sync_period)
      sync_period=rf_period;
    if (rotor_period>sync_period)
      sync_period=rotor_period;
    if (dt>sync_period)
      sync_period=dt;
  }

  size_t lnrep=0;
  if (sync_period) {

    if (spin_rate)
      trysync(sync_period,rotor_period,"MAS period");
    
    if (rf_period)
      trysync(sync_period,rf_period,"RF period");
    
    if (dt)
      lnrep=trysync(sync_period,dt,"indirect dwell time");
  }

  if (seqp)
    seqp->checksync(dt);

  return lnrep;
}

void checkgamma(double gamma)
{
  if (spin_rate) {
    if (verbose & VER_GEN) {
      std::cout << "Updating gamma angle: ";
      if (update_gamma)
	std::cout << (gamma*rad_to_deg) << " degrees\n";
      else
	std::cout << "No\n";
    }      
    if (update_gamma) {
      update_propagators=DIRTY_GAMMA;
      masterobjp->rotor_phase(gamma);
    }
  }
}

ThreadWarning<> gammacompute_warning("Gamma angle is being integrated explicitly.  Consider using 3-angle sets, such as the builtin 3zcw sets, in place of setting gamma_angles.",&NMRsim_once_warning);

//! add single alpha/beta orientation
void MasterObj::add_FID_orientation(DataStore& Specall, size_t lorient)
{ 
  curstorep=&Specall;
  double local_curscale=1.0;
  Euler powder(global_powder);

  update_gamma=false;
  static bool have_orient=false;
  static Euler last_powder;

  if (powder_array)
    update_interactions=update_gamma=true;
  else {
    if (!!powdm) {
      set_orientation(powder,local_curscale,lorient,verbose);
      current_gamma=powder.gamma;
      if (have_orient) {
	if ((powder.alpha!=last_powder.alpha) || (powder.beta!=last_powder.beta))
	  update_interactions=true;
	if (powder.gamma!=last_powder.gamma)
	  update_gamma=true;
      }
      else {
	update_interactions=update_gamma=true;
	have_orient=true;
      }
      last_powder=powder;
    }  
  }	

  if (verbose & VER_GEN)
    std::cout << "Powder alpha: " << (powder.alpha*rad_to_deg) << "  beta: " << (powder.beta*rad_to_deg) << '\n';

  if (local_curscale==0.0)
    return;
   
  //  reset();
  const size_t nactionstacks=actionstacks.size();
  for (size_t i=nactionstacks;i--;)
    actionstacks(i).reset(); //clear any stored propagators
  
  array_iter variter;
  while (variter.next()) {
    accumulating_timer<>* timerp= timers_arrayp ? &( (*timers_arrayp)(sum_index,row_index)) : NMRSIM_NULL;
    if (timerp)
      timerp->enter();

    setup_cp->enter();
    
    reset_logcount();

    if (verbose & VER_GEN)
      reset_warnings(); //!< if being verbose, allow serious warnings to repeat
        
    if (powder_array) {
      double dummy;
      set_orientation(powder,dummy,row_index,verbose);
      size_t base=varpars.size();	
      vararray(row_index,base++)=powder.alpha*rad_to_deg;
      vararray(row_index,base++)=powder.beta*rad_to_deg;
      vararray(row_index,base)=powder.gamma*rad_to_deg;
      update_interactions=true;
    }
    
    if (!optupdate) {
      update_gamma=true;
      update_interactions=true;
      update_propagators=DIRTY_ALL;
    }

    if (active2D && (verbose & VER_GEN))
      std::cout << "2D position: " << arrayinds << '\n';

    ensure_operator_matrices(); //!< always check operator matrices (quick) since these might need update independently of variables

    if (verbose & VER_GEN)
      std::cout << "Updating Hamiltonian: " << (update_interactions ? "Yes\n" : "No\n");
    
    if (update_vars) {
      save_parameters();
      restart(); //!< need restart to ensure MAS angle etc. set
    }

    if (update_interactions) {
      refresh_interactions();
      recreate_Hs();
      update_interactions=false;
      update_propagators=DIRTY_ALL; //cached propagators will be invalid
    }
    else
      checkgamma(powder.gamma);

    jitterinfo.reset(); //!< reset any MAS jitter to ensure consistent starting point
    
    if (verbose & VER_GEN)
      std::cout << "Global propagator flush level: " << update_propagators << '\n';
//     dirty_t overalldirt=DIRTY_CLEAN;
//     switch (update_propagators) {
//     case DIRTY_NONE:
//       break;
//     case DIRTY_GAMMAONLY:
//       overalldirt=DIRTY_GAMMA;
//       //      flagdirty(DIRTY_GAMMA);
//       break;
//     case DIRTY_ALL:
//       overalldirt=DIRTY_HAMILTONIAN;
//       //      flagdirty(DIRTY_HAMILTONIAN);
//       break;
//     default:
//       throw InternalError("update_propagators");
//     }

    rebuild_sequences(update_propagators);
    actionstack_t& curactionstack((nactionstacks>1) ? actionstacks(row_index) : actionstacks.front());
    curactionstack.flush_cache(update_propagators); //!< always called even if propagators are valid

    const size_t gamma_loop_steps=(use_gamma_angles || (gamma_angles<1)) ? 1 : gamma_angles;
    curscale=local_curscale/gamma_loop_steps;
    if (gamma_loop_steps>1)
      gammacompute_warning.raise();

    update_gamma=false;
    update_propagators=DIRTY_CLEAN;

    setup_cp->leave(); //!< actionstack counted separately

    for (size_t gammastep=0;gammastep<gamma_loop_steps;gammastep++) {
      if (gammastep) {
	powder.gamma=global_powder.gamma+(restrict_gamma*gammastep)/gamma_loop_steps;
	update_gamma=true;
	checkgamma(powder.gamma);
	curactionstack.flush_cache(DIRTY_GAMMA);
      }
      current_gamma=powder.gamma;
      if (verbose & VER_GEN) {
	std::cout << "Gamma angle: " << (powder.gamma*rad_to_deg) << " degrees\n";
	if (spin_rate)
	  std::cout << "Rotor phase at t=0: " << (phase(0.0)*rad_to_deg) << " degrees\n";
      }

      reset_seqlists(); //Reset "special" CycledSequences
      master_time=0.0;
      
      density=sigma0;
            
      //exec_action();    
      curactionstack.exec(this,density,master_time);

      if (transients_active()) {
	ActionAcq::trailingtransients_warning.raise();
	for (size_t i=nchannels;i--;)
	  pgenstack(i)->clear_transient();
      }
    }
    logfile_controller* logfilep(get_logfile());
    if (logfilep)
      logfilep->flush();
    if (timerp)
      timerp->leave();
  } //loop over arrayed parameters / 2D index
  if ((verbose & VER_PROFILE) && (verbose_level>1)) {
    std::cout << "Resource profile after orientation " << (lorient+1) << '\n';
    dump_timings(true);
  }
}

 
void MasterObj::operator()(size_t st,size_t end,size_t nthr) const
{
  const size_t orients = !powdm ? 1 : powdm->orientations();
  if (end<st)
    throw InternalError("MasterObj: end orientation before start orientation!");
  //   if (!powdm)
//     throw InternalError("Threaded add_FID without powder average!");
  const size_t sum_steps=sum_n0 ? sum_n0 : 1;
  const size_t chunk_steps=(end-st)*sum_steps;
  MasterObj* mutable_this=const_cast<MasterObj*>(this);
  MultiIterator sumiter(sum_ns);
  if (verbose & VER_POWDER) {
    const size_t totalsteps=orients*sum_steps;
    std::cout << "Process " << nthr << " computing summation steps " << (st+1) << " to " << end << " (out of " << totalsteps << ")\n"; //!< adjusted to totalsteps 4/8/15
  }
  bool doneETA=silent;
  timer<WallTimer> wstopwatch;    
#ifdef USE_FORK_CONTROLLER
  bool donecheck=nochecks;
#endif
  const bool processorient=!(procstacks.empty());

  for (size_t step=st;step<end;step++) {
    const size_t orient = step % orients;
    if ((orient==0) || (step==st)) {
      sum_index = step / orients;
      sumiter.reverse(suminds,sum_index);
      update_interactions=true;
      if (verbose & VER_GEN)
	std::cout << "Process " << nthr << " working on sum array position: " << suminds << '\n';
    }
    DataStore* destp=&ResultSpec;
    if (processorient) {
      make_zero(OrientSpec);
      destp=&OrientSpec;
    }
    mutable_this->add_FID_orientation(*destp,orient);
#ifdef USE_FORK_CONTROLLER
    if (!donecheck) {
      if (!(parconp->verify(destp->row()))) {
	parser_printthread() << "Pre-allocation of shared memory for data failed (try reducing data set size - currently " << destp->row().size() << " points - and/or number of processors)";
	error_abort();
      }
      donecheck=true;
    }
#endif
    if (processorient) {
      apply_procstacks(procstacks,OrientSpec);
      ResultSpec+=OrientSpec;
    }
    if (!doneETA) {
      const double curt=wstopwatch();
      if (curt>dotlimit) {
	const size_t orientsdone=step+1-st;
	const double eta=chunk_steps*curt/orientsdone;
	if (eta>siglimit)  {
	  std::cout << orientsdone << " summation step(s) (out of " << chunk_steps << ") " << (orientsdone ? "have" : "has") << " taken " << curt << " seconds.  Total estimated CPU time: ";
	  prettyprint_time(eta);
	  std::cout << std::endl;
	  doneETA=true;
	}
      }	      	    
    }
  }
}

void check_powderarray()
{
  if (!powder_array)
    return;

  ensure_basepowder();
  const size_t orients=powdm->orientations();
  //  if (orients!=array_n0) {
  //  assert(powder_type!=POW_NONE);
  // if (array_n0!=1)
  //  error_abort("{} does not match number of powder orientations");
  if (orients>1) {//!< if variable orientation, set by {} on variables
    array_dims.set(1,orients,true);
    array_dims.get(array_ns,array_skips,array_n0); //extract dimensionality
    // assert(array_n0==orients);
  }
}

//Things to do for new simulation
// void MasterObj::initrun()
// {
//   ensure_channels();
//   global_prop_flags=prop_flags.front();

//   update_interactions=true;

//   //static size_t curgamma_steps=0;

//   //  if ((curnzcw==nzcw) && (curgamma_steps==gamma_loop_steps))
//   //   return;

//   ensure_basepowder();

//   //  if (gamma_loop_steps<1)
//   //  throw InternalError("ensure_powder: called before gamma_loop_steps initialised");

//   //  curgamma_steps=gamma_loop_steps;
// //   if (gamma_loop_steps==1)
// //     powdm.reset(pregammapowdm.get(),mxflag::nondynamic); //reference
// //   else {
// //     if (pregammapowdm->angles()==3)
// //       error_abort("Can't add gamma integration to powder method already including gamma angle");
// //     powdm.reset(new WithGamma(*pregammapowdm,gamma_loop_steps),mxflag::normal);
// //     //    delete pregammapowdm;
// //   }

//   orientations=powdm->orientations();
//   if (verbose & VER_POWDER)
//     std::cout << "Powder orientations: " << orientations << '\n';

// //   if (powder_array) {
// //     if (nacqdims!=1)
// //       error_abort("powder array cannot be combined with >1 dimensions");
// //     const size_t orients=powdm->orientations();
// //     if (orients!=array_n0) {
// //       assert(powder_type!=POW_NONE);
// //       if (array_n0!=1)
// // 	error_abort("{} does not match number of powder orientations");
// //       array_dims.set(1,orients,true);
// //       array_dims.get(array_ns,array_skips,array_n0); //extract dimensionality
// //       assert(array_n0==orients);
// //     }
// //   }

//   ensure_array(); //!< can be called twice

//   if (array_n0>1) {
//     size_t usevars=varpars.size();
//     if (powder_array)
//       usevars+=3;
//     if (usevars) {//NB array parameters within sum are not included
//       vararray.create(array_n0,usevars);
//       varnames.create(usevars);
//       size_t i=0;
//       for (;i<varpars.size();i++)
// 	varnames(i)=varpars(i).name();
//       if (powder_array) {
// 	varnames(i++)="alpha";
// 	varnames(i++)="beta";
// 	varnames(i)="gamma";
//       }
//     }
//   }
//   evaluation_state=CONTEXT_MAINLOOP; //!< flag now in main loop
// }

//Things to do before sequence run
// void MasterObj::reset()
// {
//   //reset_angle(def_rotor_angle); // done in restart at lower level
//   //  maxdt=def_maxdt;
//   //  reset_pgens(); //ensure that PulseGenerators start without RF
// }


ThreadWarning<> ignoringrev_warning("-reverse flag being ignored e.g. during parallel powder averaging",&NMRsim_once_warning);
ThreadWarning<> badlychunked_warning("Number of parallel workers does not divide evenly into number of integration steps",&NMRsim_repeat_warning);

void MasterObj::calc(DataStore& Sumspec, int locverbose, bool share)
{
  evaluation_state=CONTEXT_PRE;
  //  bool isfirst=true;
#ifdef HAVE_PARSYS
  bool dothread=false;
#endif
  //  DataStore Specall;

  evaluation_index++; //incr evaluation index
  if (verbose & VER_GEN)
    std::cout << "Incrementing evaluation index to " << evaluation_index << '\n';

  const bool processorient=!(procstacks.empty());
  const size_t sumsteps=sum_n0 ? sum_n0 : 1;

  //  Sumspec.reset(!isfrequency());
  if (!have_spinsys()) {
    create_empty(Sumspec);
    if (processorient) {

      ensure_vararray();
    
      evaluation_state=CONTEXT_MAINLOOP; //!< flag now in main loop
      MultiIterator sumiter(sum_ns);

      for (sum_index=0;sum_index<sumsteps;sum_index++) {

	NMRSIM_EXPECT(sumiter.next(suminds)); //get indices
	
	if ((sumsteps>1) && (verbose & VER_GEN))
	  std::cout << "Sum array position: " << suminds << '\n';

	make_zero(OrientSpec);
	apply_procstacks(procstacks,OrientSpec);
	Sumspec+=OrientSpec;
      }
    }
  }
  else {
    if (verbose & VER_GEN)
      reset_warnings();

    ensure_channels();
    global_prop_flags=prop_flags.front();
    ensure_basepowder();
    orientations=powdm->orientations();
    if (verbose & VER_POWDER)
      std::cout << "Powder orientations: " << orientations << '\n';

    ensure_vararray();
    
    evaluation_state=CONTEXT_MAINLOOP; //!< flag now in main loop
        
    if ((verbose & VER_PROFILE) && !timers_arrayp)
      timers_arrayp=new Matrix<accumulating_timer<> >(sumsteps,array_n0);
    
    MultiIterator sumiter(sum_ns);
    
    const size_t orients= (powder_array || !powdm) ? 1 : powdm->orientations();
    
    const size_t totalsteps=sumsteps*orients;

    create_empty(Sumspec);    
//     if (Sumspec.empty())
//       create_empty(Sumspec);
//     else
//       Sumspec=complex(0.0);

    timer<> stopwatch;
    timer<WallTimer> wstopwatch;    
    
    losttrans.create(0U);
    lostsum.create(0U);

#ifdef HAVE_PARSYS
    dothread=(global_workers>0) && !powder_array && (totalsteps>=global_workers);    
    if ((totalsteps<6*global_workers) && !assumeparallel)
      dothread=false; //!< turn off parallel computation if it doesn't seem worth it

    if (dothread) {
      ensure_parallel();
      if (!parconp)
	dothread=false;
    }
    if (dothread) {
      if (powder_reverse)
	ignoringrev_warning.raise();

      size_t chunk=totalsteps/(global_workers*3);
      if (chunk==0) {
	char errmes[100];
	
	if (totalsteps<global_workers) {
	  snprintf(errmes,sizeof(errmes),"Number of integration steps (%" LCM_PRI_SIZE_T_MODIFIER "u) is less than number of workers (%" LCM_PRI_SIZE_T_MODIFIER "u)",totalsteps,global_workers);    
	  error_abort(errmes);
	}		  
	if (totalsteps % global_workers) {
	  snprintf(errmes,sizeof(errmes),"Integration steps: %" LCM_PRI_SIZE_T_MODIFIER "u, workers: %" LCM_PRI_SIZE_T_MODIFIER "u",totalsteps,global_workers);
	  badlychunked_warning.raise(errmes);
	}
	chunk=totalsteps/global_workers; //!< must be >=1
	if (verbose & VER_POWDER)
	  std::cout << "Using chunks of " << chunk << " orientation(s) over " << global_workers << " workers to integrate " << totalsteps << " powder orientations\n";
      }
      make_zero(ResultSpec);

#if defined(HAVE_FORK_CONTROLLER) && !defined(USEMPI)
      parconp->start();
      ammaster=parconp->ammaster();
      if ((verbose_level<2) && !ammaster) {
	verbose_level=0; //!< disable output from slave threads
	silent=true;
	setwarnings(BaseWarning::Silent);
      }
      share=false; //!< disable broadcast
#endif
      parconp->run(*this,totalsteps,chunk);
    }
    else {
      if (optparallel.isenabled())
	parfailed_warning.raise();
#endif    
      // end HAVE_PARSYS
      const bool dodots=!silent && !(verbose & VER_GEN) && !((verbose & VER_PROFILE) && (verbose_level>1)); //!< suppress dots if verbose or profiling
      DotMaker dotmaker(dodots ? dotlimit : -1.0);
      bool doneETA=silent || (totalsteps<2);
      bool firstsumstep=true;

      for (sum_index=0;sum_index<sumsteps;sum_index++) {

	NMRSIM_EXPECT(sumiter.next(suminds)); //get indices
	
	if ((sumsteps>1) && (verbose & VER_GEN))
	  std::cout << "Sum array position: " << suminds << '\n';
	
	update_interactions=true;
	
	//Specall.reset(!isfrequency());
      //    if (!Specall.empty())
      //  Specall=0.0;
		
	for (size_t lorient=0;lorient<orients;lorient++) {
	  dotmaker.print();
	  //	  add_FID_orientation(Specall,lorient);
	  if (processorient) {
	    make_zero(OrientSpec);
	    add_FID_orientation(OrientSpec,lorient);
	    apply_procstacks(procstacks,OrientSpec);
	    if (firstsumstep) {
	      Sumspec=OrientSpec; //!< copies across processing etc.
	      firstsumstep=false;
	    }
	    else
	      Sumspec+=OrientSpec;
	  }
	  else
	    add_FID_orientation(Sumspec,lorient);

	  const double curt=wstopwatch();
	  if (!doneETA && (curt>dotlimit)) {
	    const double eta=totalsteps*curt/(sum_index*orients+lorient+1);
	    if (eta>siglimit)  {
	      std::cout << "First ";
	      if (lorient)
		std::cout << (lorient+1) << " orientations";
	      else
		std::cout << "orientation";
	      std::cout << " (out of " << totalsteps << " total summation steps) " << (lorient ? "have" : "has") << " taken " << curt << " seconds.  Total estimated calculation time: ";
	      prettyprint_time(eta);
	      std::cout << std::endl;
	    }	      	    
	    doneETA=true;
	  }
	}
      }
      dotmaker.flush();
      //      } //have_spinsys
      //     else {
      //       if (Specall.empty()) //create empty
      // 	create_empty(Specall);
      //     }

      // 	if (isfirst) {
      // 	  Sumspec=Specall;
      // 	  isfirst=false;
      // 	}
      // 	else
      // 	  Sumspec+=Specall;
#ifdef HAVE_PARSYS
    }
    if (dothread) {
#ifdef USE_FORK_CONTROLLER
      parconp->sum_join(Sumspec,ResultSpec);
#else
      parconp->sum(Sumspec,ResultSpec);
#endif
    }
#endif
    if (ammaster) {
      if ((totalsteps>1) && (locverbose & VER_GEN)) {
	std::cout << "Time elapsed: " << stopwatch() << " s\n";
#ifdef HAVE_PARSYS
	if (dothread)
	  std::cout << "Wall clock time: " << wstopwatch() << " s\n";
#endif
      }
      //	  apply_procstacks(procstacks,Specall);
      if (!losttrans.empty() && HistogramMaker::loss_warning.enabled()) {
	std::ostringstream str(std::ostringstream::out);
	str << ": number lost: " << losttrans << "\t sum: " << lostsum;
	HistogramMaker::loss_warning.raise(str.str().c_str());
      }
    }
    //#ifdef HAVE_PARSYS
    //if (share && dothread)
    //  parconp->broadcast(Sumspec.row());
    //#endif      
  }
  sum_index=-1; //indicate sum_index is no longer valid
  evaluation_state=CONTEXT_POSTCALC;
  if (ammaster || share) {
    ensure_array(); //!< will not have been called if no spinsys
    if (ammaster)      
      apply_procstacks(postprocstacks,Sumspec);
#ifdef HAVE_PARSYS
    if (share && dothread)
      parconp->broadcast(Sumspec.row());
#endif
    if (Sumspec.empty()) {
      global_np=0;
      global_nps.clear();
    }
    else {
      global_np=Sumspec.size();
      Sumspec.copy_structure(global_nps);
    }
  }    
}

const rmatrix& make_powderarray()
{
  ensure_array(); //!< need to create array structure in order to create powder method
  static rmatrix pvalues(array_n0,4);
  array_iter iter;
  size_t r=0;
  while (iter.next()) {
    BaseList<double> crow(pvalues.row(r++));
    crow(0U)=global_powder.alpha;
    crow(1U)=global_powder.beta;
    crow(2U)=global_powder.gamma;
    crow(3U)=1.0;
  }
  if ((verbose & VER_POWDER) && (verbose_level>1))
    std::cout << "Explicit powder matrix\n" << pvalues << '\n';
  return pvalues;
}

ThreadWarning<InvalidParameter> invalidweighting_warning("found invalid weighting factor(s) (<=0.0) in .cry file",&NMRsim_repeat_warning);

void checkangles(rmatrix& angles)
{
  const size_t weightcol=angles.cols(); 
  bool foundbad=false;
  for (size_t i=angles.rows();i--;) {
    double* arow(angles.vector(i));
    if (arow[weightcol]<=0.0)
      foundbad=true;
    for (size_t j=weightcol;j--;)
      arow[j]*=deg_to_rad;
  }
  if (foundbad)
    invalidweighting_warning.raise();
}

void read_simpsoncry(rmatrix& angles, const char *fname)
{
  char errmess[128];
  FILE *fp=pathopen(fname,"r");
  if (fp==NMRSIM_NULL) {
    snprintf(errmess,sizeof(errmess),"read_simpsoncry: failed to open %s",fname);
    throw Failed(errmess);
  }
  try {
    int N;

    if (fscanf(fp,"%i",&N)!=1) {
      snprintf(errmess,sizeof(errmess),"read_simpsoncry: failed to parse %s",fname);
      throw Failed(errmess);
    }

    angles.create(N,3);

    for (int i=0;i<N;i++) {
      double *arow = angles.vector(i);
      if (fscanf(fp,"%lg%lg%lg",arow,arow+1,arow+2)!=3)
	throw Failed("Corrupt SIMPSON crystal file?");
    }
  }
  catch (...) {
    fclose(fp);
    throw;
  }
  fclose(fp);
  checkangles(angles);
}

ContextWarning<> ignoringrange_warning("range qualifier flags have no effect on explicit sample (ignored)",&NMRsim_repeat_warning);

void checkzcw(char*& namep)
{
  if (*namep==':') {
    namep++;
    zcwisindex=true;
  }
}

ContextWarning<> unlikelyzcw_warning("ZCW index is very high (>100).  Use zcwXXX rather than zcw:XXX to set explicit number of orientations",&NMRsim_once_warning);

void parse_crystal_file()
{
  static flagsmap_type crystal_flags;
  if (crystal_flags.empty()) {
    crystal_flags["reverse"]=CRYS_REVERSE;
    crystal_flags["sphere"]=CRYS_SPHERE;
    crystal_flags["hemisphere"]=CRYS_HEMISPHERE;
    crystal_flags["octant"]=CRYS_OCTANT;
    crystal_flags["start"]=CRYS_START;
    crystal_flags["end"]=CRYS_END;
    crystal_flags["middle"]=CRYS_MIDDLE;
    crystal_flags["both"]=CRYS_BOTH;
    crystal_flags["print"]=CRYS_PRINT;
  }
  const char* syntaxstr="\nSyntax:\tcrystal_file single <alpha> <beta> <gamma>\nor\tcrystal_file [3zcw<N>|zcw<N>|sphericalzcw<N>|beta<N>|alphabeta<Nalpha>,<Nbeta>] [-sphere|-hemisphere|-octant] [-start|-middle|-end|-both] [-reverse|-print]\nor\tcrystal_file <crystal orientations file> [-reverse|-print]\n";

  zcwisindex=false; //!< by default, use "old" switching behaviour
  if (!are_left())
    error_abort(syntaxstr);

  crystal_filep=parse_string();
  if (strcmp(crystal_filep,"?")==0)
    error_abort(syntaxstr);

  char* last=crystal_filep+strlen(crystal_filep)-1;
  if ((*crystal_filep=='{') ^ (*last=='}'))
    error_abort("Mismatched {}");

  if (*last=='}') {
    powder_array=true;
    crystal_filep++;
    *last--='\0';
  }
  if (strcmp(crystal_filep,"single")==0) {
    parse_system_variable(v_alpha,F_DENYSUM);
    parse_system_variable(v_beta,F_DENYSUM);
    parse_system_variable(v_gamma,F_DENYSUM);
    if (!v_alpha.isconstant() || !v_beta.isconstant() || !v_gamma.isconstant())
      powder_array=true;
    have_orientation=true;
    return;
  }
  size_t nargs=1;
  powder_type=POW_NONE;
  if (STRNCMP(crystal_filep,"3zcw")==0) {
    powder_type=POW_ZCW3;
    crystal_filep+=4;
    checkzcw(crystal_filep);
  }
  else {
    if (STRNCMP(crystal_filep,"zcw")==0) {
      powder_type=POW_ZCW;
      crystal_filep+=3;
      checkzcw(crystal_filep);
    }
    else {
      if (STRNCMP(crystal_filep,"sphericalzcw")==0) {
	powder_type=POW_SPHERICALZCW;
	crystal_filep+=12;
	checkzcw(crystal_filep);
      }
      else {
	if (STRNCMP(crystal_filep,"beta")==0) {
	  powder_type=POW_BETA;
	  crystal_filep+=4;
	}
	else {
	  if (STRNCMP(crystal_filep,"alphabeta")==0) {
	    powder_type=POW_ALPHABETA;
	    crystal_filep+=9;
	    nargs=2;
	  }
	  else {
	    //try to read crystal file
	    try {
	      FILE* fp=pathopen(crystal_filep,"r");
	      if (!fp)
		throw Failed("X");
	      read_matrix(powder_weights,fp);	// don't worry about resource leak on exception (called rarely)
	      fclose(fp);
	      nargs=0;
	      switch (powder_weights.cols()) {
	      case 3: case 4:
		break;
	      default:
		parser_printcontext() << "Crystal orientations file must have 3 or 4 columns: " << crystal_filep << '\n';
		error_abort();
	      }
	      checkangles(powder_weights);
	    }
	    catch (...) {
	      try {
		char fname[128];
		snprintf(fname,sizeof(fname),"%s.cry",crystal_filep);
		read_simpsoncry(powder_weights,fname);
	      } catch(...) {	  
		parser_printcontext() << "Failed to open " << crystal_filep << " for crystal orientations or unknown sampling specification (known: zcw, sphericalzcw, 3zcw, beta, alphabeta, single)\n";
		error_abort();
	      }
	    }
	    nzcw=1;
	    powder_type=POW_FILE;
	  }
	}
      }
    }
  }
    
  switch (nargs) {
  case 0:
    break;
  case 2:
    if (sscanf(crystal_filep,"%i,%i",&nalpha,&nzcw)!=2)
      error_abort("failed to parse numerical qualifier <alpha angles>,<beta angles>");
    break;
  case 1:
    if (sscanf(crystal_filep,"%i",&nzcw)!=1)
      error_abort("failed to parse numerical qualifier (zcw<N>|beta<N>)");
    if ((powder_type!=POW_BETA) && (nzcw>=100) && zcwisindex)
      unlikelyzcw_warning.raise();
    break;
  default:
    throw InternalError("crystal_file");
  }
//     const char* colonpos=strrchr(crystal_filep,':');
//     if (colonpos) {
//       colonpos++;
//       try {
// 	rangequal=PowderMethod::range_type(colonpos);
//       }
//       catch (...) {
// 	parser_printcontext() << "Unrecognised angle range qualifier: " << colonpos << '\n';
// 	return false;
//       } 
//     }
//  }

  if (parser_isnormal() && (powder_type==POW_ZCW3))
    restrict_gamma=deg_to_rad*parse_double();

  int flags=parse_flags(crystal_flags);
  if (flags & CRYS_REVERSE) {
    powder_reverse=true;
    flags-=CRYS_REVERSE;
  }
  if (flags & CRYS_PRINT) {
    powder_print=true;
    flags-=CRYS_PRINT;
  }
  if (powder_type==POW_FILE) {
    if (flags)
      ignoringrange_warning.raise();
    return;
  }
  static const int rangeflags(CRYS_SPHERE | CRYS_HEMISPHERE | CRYS_OCTANT);
  switch (flags & rangeflags) {
  case 0: case CRYS_SPHERE:
    rangequal=sphere;
    break;
  case CRYS_HEMISPHERE:
    rangequal=hemisphere;
    break;
  case CRYS_OCTANT:
    rangequal=octant;
    break;
  default:
    error_abort("cannot specify more than one range qualifier");
  }
  static const int offsetflags(CRYS_START | CRYS_END | CRYS_MIDDLE | CRYS_BOTH);
  switch (flags & offsetflags) {
  case 0: case CRYS_MIDDLE:
    offsetqual=PowderMethod::middle;
    break;
  case CRYS_END:
    offsetqual=PowderMethod::end;
    break;
  case CRYS_START:
    offsetqual=PowderMethod::start;
    break;
  case CRYS_BOTH:
    offsetqual=PowderMethod::both;
    break;
  default:
    error_abort("cannot specify more than one offset qualifier");
  }      
}

template<class PowdM> int getorients(int orients, const char* name, size_t trigger =0)
{
  int n=PowdM::orientations_to_N(orients);
  if (n<0) {
    parser_printthread(std::cerr) << "Unrecognised number of orientations.  ";
    n=1;
    while (PowdM::N_to_orientations(n)<orients)
      n++;
    if (n==1)
      std::cerr << "Next allowed number is ";
    else
      std::cerr << "Nearest allowed numbers are " << PowdM::N_to_orientations(n-1) << " or ";
    std::cerr << PowdM::N_to_orientations(n) << '\n';
    if (n<trigger)
      std::cerr << "If an orientation set index is intended, please use the syntax <type>:<set no.> e.g. zcw:8\n";
    error_abort();
  }
  return n;
}

ThreadWarning<> orientationsfixed_warning("number of orientations is fixed for this powder sampling",&NMRsim_once_warning);
ThreadWarning<> gamma3angle_warning("explicit gamma angle integration (gamma_angles) makes little sense combined with a 3-angle powder method",&NMRsim_once_warning);

void ensure_basepowder()
{
  if (nzcw==curnzcw)
    return;

  switch (powder_type) {
  case POW_NONE:
    if (powder_array) {
      const rmatrix& pvalues(make_powderarray());
      powdm.reset(new ExplicitSampling(pvalues));
    }
    else
      //      powdm.reset(new PowderSingle(Euler(global_powder.alpha,global_powder.beta,gamma_zero)));
      powdm.reset(new PowderSingle(Euler(global_powder.alpha,global_powder.beta,global_powder.gamma)));
    break;
  case POW_ZCW3:
    if (!zcwisindex)
      nzcw=getorients<ZCW3>(nzcw,"3zcw",40);
    powdm.reset(new ZCW3(nzcw,rangequal,offsetqual,restrict_gamma));
    break;
  case POW_ZCW:
    if (!zcwisindex)
      nzcw=getorients<PlanarZCW>(nzcw,"zcw",20);
    powdm.reset(new PlanarZCW(nzcw,rangequal,offsetqual));
    break;
  case POW_SPHERICALZCW:
    if (!zcwisindex)
      nzcw=getorients<SphericalZCW>(nzcw,"sphericalzcw",20);
    powdm.reset(new SphericalZCW(nzcw,rangequal,offsetqual));
    break;    
  case POW_BETA: case POW_ALPHABETA:
    powdm.reset(new PlanarGrid(nalpha,nzcw,rangequal,offsetqual));
    break;
  case POW_FILE:
    if (nzcw!=1)
      orientationsfixed_warning.raise();
    powdm.reset(new ExplicitSampling(powder_weights));
    break;
  default:
    throw InternalError("Unknown powder type");
  }
  curnzcw=nzcw;

  if ((verbose & VER_POWDER) || powder_print)
    std::cout << "Set up powder averaging with quality factor " << nzcw << " corresponding to " << powdm->orientations() << " orientation(s)\n";

  if ((gamma_angles>1) && (powdm->angles()==3))
    gamma3angle_warning.raise();
}

void dump_operator(const char* name, const productoperator_spec& spec, const BlockedOperator& op, std::ostream& ostr =std::cout)
{
  ostr << name << ": ";
  spec.print(ostr);
  ostr << '\n';
  if (verbose_level>1)
    spyprint(ostr,op);
}

void update_auxiliary_vars()
{
  if (detect_specp) {
    const productoperator_spec& detect_spec(detect_specp->current());
    gammaX=gamma(nuclei_spec(detect_spec.nucleus()));
    detect_freq=proton_freq*gammaX/get_gamma1H();
  }
}

void ensure_sigma0detect()
{
  if (!sigma0_specp || !detect_specp)
    throw Failed("start_operator and/or detect_operator not set e.g. no default values in heteronuclear systems");

  const bool needupdate=!sigma0 || !detect || !(sigma0_specp->isconstant()) || (!detect_specp->isconstant());

  if (verbose & VER_GEN)
    std::cout << "Updating start/detect operators: " << (needupdate ? "Yes\n" : "No\n");
  
  if (!needupdate)
    return;

  productoperator_spec sigma0_spec(sigma0_specp->current());
  const productoperator_spec& detect_spec(detect_specp->current());

  detectED=false;
  bool changesigma0=!sigma0_specp->isconstant();
  globalscale=1.0;
  if (optED() && ((!haveactions && !active2D) | optED.isenabled()) && arematching(sigma0_spec,detect_spec)) {
    if (sigma0_spec==detect_spec)
      detectED=true;
    else {
      if (masterobjp->Hstructp->isblocked(sigma0_spec.nucleus())) { //If truncating from x -> +/-, then should multiply by 0.5
	detectED=true;
	globalscale=0.5;
	sigma0_spec&=detect_spec;
	changesigma0=true;
      }
    }
//     if (detectED && enablespecial)
//       //check for spin-1/2 special
//       detectspecial= (detect_spec.unique().op=='+') && (phasep->structure().mzblocks()==2);
  }    
  optED.check(detectED);
//     if (detectED)
//       std::cout << "Heteronuclear detect special: " << (detectspecial ? "Yes\n" : "No\n");
  if (!sigma0 || changesigma0) {
    sigma0=masterobjp->make(sigma0_spec);
    if (verbose & VER_GEN)
      dump_operator("Initial density matrix (start)",sigma0_spec,sigma0);
  }
  if (!detect || !detect_specp->isconstant()) {
    conj_transpose(detect,masterobjp->make(detect_spec)); //can't always use detectED, so create anyway
    //if detect nucleus gamma<0, conjugate signal before further processing
    if (verbose & VER_GEN) {
      //      if (!detectED)
      dump_operator("Detection operator (detect)",detect_spec,detect);
      std::cout << "conjugate FID / reverse spectrum: " << ((detect_freq<0) ? "Yes\n" : "No\n");    
    }
  }  

  //! make spin order filter matrices (depends on sigma0 block structure)
  //! postpone to point of use
//   const bool isverb=((verbose & VER_GEN) && (verbose_level>1));
//   const matrixmap_type::iterator end(matrixmap.end());
//   matrixmap_type::iterator start(matrixmap.begin());
//   while (start!=end) {
//     matrixspec& mspec(start->second);
//     if (mspec.type==matrixspec::BOOL) {
//       filter_def& def( const_cast<filter_def& >(mspec.asfilter()));
//       if (def.isspinorder() && (changesigma0 || !(def.filterp))) {
// 	const spinorder_spec& spec(def.spec(Type2Type<spinorder_spec>()));
// 	if (isverb)
// 	  std::cout << "Making spin order matrix: " << spec << '\n';	
// 	def.filterp=masterobjp->create_spinorder(spec,sigma0.blockstructure());
//       }
//     }
//     ++start;
//   }
}

void parse_powderquality()
{
  if (!are_left()) {
    std::cout << "powderquality: ";
    if (havepowderquality())
      std::cout << "<undefined>\n";
    else
      std::cout << nzcw << '\n';
  }
  else {
    if (powder_array)
      error_abort("Cannot combine powderquality with array over powder orientations");
    if (havepowderquality()) {
      parse_system_variable(v_powderquality);
      //      nzcw=parse_unsigned();
      if (!havepowderquality()) //!< not redundant - may have changed!
	error_abort("Invalid powder quality (>=2)");
    }
    else
      error_abort("Cannot modify powderquality for this sampling");
  }		    
}

void parse_n(int dim)
{
  int n=parse_int();
  int ni_skip= are_left() ? parse_int() : 1;
  set_n(dim,n,ni_skip);
  if (n==1)
    v_ni.isconstant(true); //!< flag ni now set (and fixed)
}

bool MasterObj::check_powder(const HamiltonianStore<space_T>& interactions, const Euler& testangle) const
{
  size_t skip=4;
  if (powder_type!=POW_BETA) {
    switch (rangequal) {
    case sphere: return true;
    case hemisphere: skip=4; break;
    case octant: skip=1; break;
    default: throw InternalError("Unknown angle range qualifier");
    }
  }
  BlockedMatrix<complex> H;
  {
    const HamiltonianStore<space_T> interactions_LF(interactions,testangle);    
    add_Hamiltonian(H,interactions_LF);
  }
  if (!H) {
    if (verbose & VER_GEN)
      std::cout << "System Hamiltonian: empty (no active interactions)\n";
    return true;
  }
  if (verbose & VER_GEN)
    std::cout << "System Hamiltonian (no RF):\n" << H << '\n';
  const double normH=norm(H);
  const double talpha=testangle.alpha;
  const double tbeta=testangle.beta;
  
  for (size_t whichq=skip;whichq<8;whichq+=skip) {
    Euler testangle2;
    if (powder_type==POW_BETA)
      testangle2=Euler(talpha+M_PI/3,tbeta,0); //check invariance w.r.t. alpha
    else {
      switch (whichq) {
      case 1: testangle2=Euler(talpha+M_PI/2,tbeta,0); break;
      case 2: testangle2=Euler(talpha+M_PI,tbeta,0); break;
      case 3: testangle2=Euler(talpha+3*M_PI/2,tbeta,0); break;
      case 4: testangle2=Euler(talpha+M_PI,M_PI-tbeta,0); break;
      case 5: testangle2=Euler(talpha+3*M_PI/2,M_PI-tbeta,0); break;
      case 6: testangle2=Euler(talpha,M_PI-tbeta,0); break;
      case 7: testangle2=Euler(talpha+M_PI/2,M_PI-tbeta,0); break;
      default: throw InternalError("Shouldn't happen");
      }
    }
    BlockedMatrix<complex> Htmp;
    const HamiltonianStore<space_T> interactions_LF(interactions,testangle2);    
    add_Hamiltonian(Htmp,interactions_LF);
    Htmp-=H;
    const double fnorm=norm(Htmp)/normH;
    const bool beverb=(verbose & VER_GEN) && (verbose_level>1);
    if (beverb)
      std::cout << "Quadrant " << whichq << ": " << fnorm << '\n';
    if (fnorm>1e-6) {
      if (powder_type==POW_BETA)
	parser_printthread(std::cerr) << "Hamiltonian is not invariant with alpha.  Adjust powder averaging method in crystal_file\n";
      else {
	if (beverb)
	  std::cout << "Angle check failed. Difference matrix\n" << Htmp << '\n';
	else
	  std::cerr << "Hamiltonian does not have required spatial symmetry.  Remove any -octant/-hemisphere in crystal_file or set -debug for more verbose output.\n";
      }
      return false;
    }
  }
  return true;
}

ThreadWarning<> nogammaspec_warning("neither gamma_angles, gamma_zero, nor 3 angle integration specified in spinning simulation.  Defaulting to single gamma angle of zero.",&NMRsim_once_warning);
ThreadWarning<> gammasetinstatic_warning("gamma_angles set in static simulation (ignored)",&NMRsim_once_warning);

void verify_powderaverage()
{
  if (gamma_angles) {
    if (gamma_angles<2)
      error_abort("gamma_angles must be >1 (set gamma_zero for gamma_angles=1)");
    if (spin_rate==0.0) {
      gammasetinstatic_warning.raise();
      gamma_angles=0;
    }
  }
  ensure_basepowder();
  command_Factory_t& par_Factory(get_par_Factory());
  static const bool have_gamma_zero=havekey(par_Factory,"gamma_zero");
  static const bool have_gamma_angles=havekey(par_Factory,"gamma_angles");
  const bool is3angle=(powdm->angles()==3) || have_orientation;
  if (is3angle && have_gamma_zero)
    error_abort("Can't set both gamma_zero and single crystal orientation / 3 angle integration");

  if (spin_rate && !have_gamma_angles && !have_gamma_zero && !is3angle) {
    nogammaspec_warning.raise();
    global_powder.gamma=0.0;
  }
}

void DataStore::create(size_t r, size_t c, const complex& v)
{
  //  sws_.create(r,swv);
  if (c==0)
    throw InternalError("DataStore::create: zero columns");
  update_auxiliary_vars(); //!< bodge to ensure detect_freq created
  procstates_.reserve(ndims);
  procstates_.create(size_t(0));
  for (size_t i=0;i<ndims-1;i++)
    procstates_.push_back(processing_state(sws(i),true));//,sfrqs(i)));
  //  sfrqs(ndims-1)=detect_freq;
  const std::pair<double,double> swref(actual_swref());
  procstates_.push_back(processing_state(swref.first,!isfrequency(),detect_freq,swref.second));
  set(Int2Type<2>()).create(r,c,v);  
}

void dump_arrayprofile(const Matrix<accumulating_timer<> >& timers_array, double mult, const char* units)
{
  const size_t sumrows=timers_array.rows();
  const size_t rowrows=timers_array.cols();
  std::cout << "Profile by summation index (columns) and row index (rows).  Units: " << units << "\n\n";
  
  if (sumrows>1) {
    for (size_t i=0;i<sumrows;i++) 
      std::cout << '\t' << (i+1);
    std::cout << "\tTotal\n";
  }
  LIST<double> rsum(sumrows,0.0);
  for (size_t r=0;r<rowrows;r++) {
    if (rowrows>1)
      std::cout << (r+1) << '\t';
    double tots=0.0;
    for (size_t i=0;i<sumrows;i++) {
      const double val=mult*timers_array(i,r)();
      std::cout << val << '\t';
      tots+=val;
      rsum(i)+=val;
    }
    if (sumrows>1)
      std::cout << tots;
    std::cout << '\n';
  }
  if (rowrows>1) {
    std::cout << "Total";
    for (size_t i=0;i<sumrows;i++)
      std::cout << '\t' << rsum(i);
    std::cout << '\n';
  }
}
