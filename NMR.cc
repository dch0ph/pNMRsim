#include "NMRsim.h"
#include "NMRsim_logfile.h"
#include "Parser.h"
#include "NMRsim_matrixspec.h"
#include "NMRsim_MasterObj.h"

matrixmap_type matrixmap;
prop_t prop_method=PROP_DIAG;
//bool have_spinsys=false;

#define MAXDT_TRIGGER 10
#define CHEBYSHEV_TRIGGER 5
#define ZEROTOLERANCE_TRIGGER 0.1

void make_par_variables();
static size_t lastrank=0;

int global_prop_flags=0;
double chebyshev_crossover=DEFAULT_CHEBYSHEV_CROSSOVER;
double eigtolerance=1e-6;
double zerotolerance=0.0; //don't round close-to-zero values by default

void parse_operator();

static const double rad_to_deg=180.0/M_PI;
static const double deg_to_rad=M_PI/180.0;

#ifndef DISABLE_EXEC
//#include "NMRsim_Process.h"
#include "NMRsim_spinsys.h"
#include "Propagation.h"
#include "AsynchronousPropagation.h"
#include "timer.h"

void parse_sideband();
void parse_eigenvalue();
void parse_prop_method();
void parse_chebyshev_iterations();
void parse_ordermatrix(const char*, bool, bool =false);
void parse_zerotolerance();
void parse_maxdt();
void parse_steps_per_rotation();
void parse_maxjumpdt();
void parse_makeoperator();

acc_t acc_mode=ACC_TIME;
bool isfrequency() { return (acc_mode==ACC_FREQ); }

option optsync("synchronise");
option optrealH("realhamiltonian","real Hamiltonian");

struct statsobj {

  void create(size_t);

  List<count_t> count; //!< used 32 bit integer as Matlab format is not well defined for 64 bit integers
  List<count_t> nzcount;
  List<complex> sum;
  List<double> normsum;
  mutable List<double> sqrtnormsum; //!< need to store in non-volatile location

  void add(size_t n,const complex&);
  void add(size_t n) { count(n)++; }
};

void statsobj::create(size_t n)
{
  count.create(n,count_t(0));
  nzcount.create(n,count_t(0));
  sum.create(n,complex(0.0,0.0));
  normsum.create(n,0.0);
}

void statsobj::add(size_t n,const complex& v) 
{
  nzcount(n)++;
  count(n)++;
  sum(n)+=v;
  normsum(n)+=norm(v);
}

class operator_stats {
public:
  operator_stats(const BlockedOperator&, const statistics_t&, double tolv);
  void print(std::ostream&) const;
  void log(logfile_controller&, const char*) const;
private:
  size_t nstates;
  size_t nspins_all;
  size_t nspins_sel;
  statsobj coherence_stats;
  statsobj spinorder_stats;

  List<size_t> countnz; //!< internal buffer of number of set bits in each state
  static void printstats(std::ostream&, const statsobj&, int);
  static void logstats(logfile_controller&, const statsobj&, const char*, const char*);
};

void operator_stats::printstats(std::ostream& ostr, const statsobj& obj, int offset)
{
  for (size_t i=0;i<obj.nzcount.size();i++) {
    const size_t nz=obj.nzcount(i);    
    ostr << (int(i)+offset) << ": " << nz << " non-zero states out of " << obj.count(i);
    if (nz)
      ostr << " with sum of " << obj.sum(i) << " and norm of " << std::sqrt(obj.normsum(i));
    ostr << '\n';
  }
}

void operator_stats::print(std::ostream& ostr) const
{
  ostr << "Decomposition of " << nstates << " x " << nstates << " operator matrix over " << nspins_sel << " active spins:\n";
  ostr << "By coherence order:\n";
  printstats(ostr,coherence_stats,-int(nspins_sel));
  ostr << "By spin order:\n";
  printstats(ostr,spinorder_stats,0);
}

void operator_stats::logstats(logfile_controller& logcon, const statsobj& stats, const char* base, const char* qualtitle)
{
  char superbase[256];
  if (qualtitle)
    snprintf(superbase,sizeof(superbase),"%s_%s",qualtitle,base);
  else
    strncpy(superbase,base,sizeof(superbase));

  char title[280];
  snprintf(title,sizeof(title),"%s_nonzerocounts",superbase);
  logcon.write(stats.nzcount,title);
  snprintf(title,sizeof(title),"%s_counts",superbase);
  logcon.write(stats.count,title);
  snprintf(title,sizeof(title),"%s_sums",superbase);
  logcon.write(stats.sum,title);

  snprintf(title,sizeof(title),"%s_normsum",superbase);
  //  List<double> sqrtnorm(stats.normsum.size());
  stats.sqrtnormsum.create(stats.normsum.size());
  for (size_t i=stats.normsum.size();i--;)
    stats.sqrtnormsum(i)=std::sqrt(stats.normsum(i));
  logcon.write(stats.sqrtnormsum,title);
}

void operator_stats::log(logfile_controller& logcon, const char* title) const
{
  smartptr<logfile_controller::composite_guard,false> guardp;

  const char* usetitle=logcon.supports_composite() ? NMRSIM_NULL : title;
  if (usetitle==NMRSIM_NULL)
    guardp.reset(new logfile_controller::composite_guard(logcon,title));
  
  logstats(logcon,coherence_stats,"coherence",usetitle);
  logstats(logcon,spinorder_stats,"spinorder",usetitle);
}

size_t log2(size_t n)
{
  size_t logval=0;
  while (n) {
    if (n==1)
      return logval;
    logval++;
    n>>=1;
  }
  throw InvalidParameter("log2: input is not a power of 2");
}

ThreadWarning<> nomatchingstates_warning("No states survive both spin selection and coherence filter - statistics will all be zero!",&NMRsim_once_warning);

operator_stats::operator_stats(const BlockedOperator& op, const statistics_t& statssel, double tol)
{
  const cmatrix a(op.full());

  const int lverbose=(verbose & VER_GEN) ? verbose_level : 0;
  state_t sel=statssel.first;

  const Matrix<bool>* filterp = statssel.second ? &((statssel.second)->ensurefull()) : NMRSIM_NULL;

  nstates=a.rows();
  nspins_all=log2(nstates); //!< deduce number of spins from matrix size (emergency catch of non-spin-1/2-system)
  if (sel==0) {
    sel=(1<<nspins_all)-1;
    nspins_sel=nspins_all;
  }
  else {
    nspins_sel=countnonzero(sel);
    if ((nspins_sel<1) || (nspins_sel>nspins_all))
      throw Failed("operator_stats: unexpected number of selected spins");
  }

  coherence_stats.create(2*nspins_sel+1);
  spinorder_stats.create(nspins_sel+1);

  if (lverbose>1)
    std::cout << "Selection mask: " << std::hex << sel << '\n';

  countnz.create(nstates);
  for (state_t r=nstates;r--;) {
    const state_t state = r & sel;
    countnz(r)=countnonzero(state);
  }

  bool foundstate=false;
  const double tol2=tol*tol;
  for (state_t r=nstates;r--;) {
    const state_t bra= r & sel;
    const size_t bracount=countnz(bra);
    for (state_t c=nstates;c--;) {
      if (filterp && !(*filterp)(r,c))
	continue; //!< uncounted state
      foundstate=true;

      const state_t ket=c & sel;
      const size_t ketcount=countnz(ket);
      const size_t spinorder=countnz(bra ^ ket);
      const int coherence=int(ketcount)-int(bracount);
      const size_t coherence_i=nspins_sel+coherence;
      const complex& val(a(r,c));
      if (norm(val)>tol2) {
	coherence_stats.add(coherence_i,val);
	spinorder_stats.add(spinorder,val);
      }
      else {
	coherence_stats.add(coherence_i);
	spinorder_stats.add(spinorder);
      }
    }
  }
  if (!foundstate)
    nomatchingstates_warning.raise();
} 

namespace {
  struct Proxy_ {
    Proxy_() {
      optional_map_t& optional_map(get_optional_map());
      optional_map["minimalupdating"]=&optupdate;
      optional_map["synchronise"]=&optsync;
      optional_map["realhamiltonian"]=&optrealH;
      
      addprepar_callback(make_par_variables);
    }
  };  
  static Proxy_ proxy;
}

LIST<InternalCommand*> internal_timers;

InternalCommand* create_timer(const char* name)
{
  InternalCommand* timerp=new InternalCommand(name); // will give leak in valgrind
  internal_timers.push_back(timerp);
  return timerp;
}

InternalCommand* setup_cp=create_timer("Simulation set up");

BlockedOperator density; //!< current density matrix
int evaluation_index=0;
size_t quadrupole_order=1;
//size_t gamma_loop_steps=0;
double detect_freq=0.0;
double proton_freq=0.0;
double spin_rate=0.0;
double spin_rate_jitter=0.0;
double spin_rate_jitter_cortime=0.0;
int gamma_angles=0;
double def_rotor_angle=MAGIC_ANGLE; 
double rotor_orient=0.0;
CrystalStructure* cstructp=NMRSIM_NULL;
double maxdtv=1.0e-6;
double maxjumpdt=0.0;
//bool forcegcompute=false;

IntervalSampler jitterinfo;

int eigenvalue=NMRSIM_INVALID;
int eigenvalues=NMRSIM_INVALID;
int sideband=NMRSIM_INVALID;
double current_gamma=NMRSIM_INVALID;
spin_system* sysp=NMRSIM_NULL;
Euler global_powder(0.0,0.0,0.0);
double gamma_zero=0.0;

double gammatimeoffset() {
  if ((spin_rate==0.0) || (current_gamma==NMRSIM_INVALID))
    throw InternalError("gammatimeoffset");
  return current_gamma/(TWO_PI*spin_rate);
}

static size_t getnspins()
{
  if (!sysp)
    error_abort("cannot refer to spin indices before spin system defined");
  return sysp->nspins();
}

 void adjust_and_verify_spin_indices(BaseList<size_t> arg1)
 {
   const size_t n=getnspins();
   for (size_t i=0;i<arg1.size();i++) {
     size_t& ind(arg1(i));
     if ((ind<1) || (ind>n)) {
       parser_printcontext() << "spin index (" << ind << ") out of range\n";
       error_abort();
     }
     ind--;
   }
 }

//! NB look for exact zeros as these are likely to correspond to empty blocks
template<class M> void dump_density(const BlockedMatrix<M>& a, std::ostream& ostr =std::cout)
{
  if ((verbose & VER_PROFILE)==0)
    return;
  const BaseList<M> asrow(a.row());
  const size_t items=asrow.size();
  size_t nz=0;
  const M zero(0);
  for (size_t i=items;i--;) {
    if (asrow(i)!=zero)
      nz++;
  }
  ostr << "Density: " << (100.0*nz/items) << "%\n";
}

template<class M> void print_structure(std::ostream& ostr,const BlockedMatrix<M>& a)
{
  for (size_t i=0;i<a.size();i++) {
    if (i)
      ostr << ", ";
    ostr << a.rows(i) << 'x' << a.cols(i);
  }
}

//! output blocked matrix ('spy' if precision 0)
template<class M> void spyprint(std::ostream& ostr,const BlockedMatrix<M>& a, int structureflags)
{
  switch (structureflags) {
  case PM_STRUCTURE:
    print_structure(ostr,a);
    break;
  case 0:
    //  if (precision)
    ostr << a;
    //  else
    //spy(ostr,a);
    dump_density(a,ostr);
    break;
  default:
    error_abort("putmatrix: only -structure flag is valid for a simple matrix");
  }
}

//! output operator matrix ('spy' if precision 0)
void spyprint(std::ostream& ostr,const BlockedOperator& a, int structureflags)
{
  if (structureflags) {
    const size_t eigblks=a.eigblocks();    
   if (eigblks!=1)
     error_abort("putmatrix: structure output not implemented for operators with eigenvalue structure e.g. periodic bases");
    
    size_t r,c,mzeig;
    const block_pattern& blkstr(a.blockstructure());
    if (structureflags==PM_STRUCTURE)
      ostr << blkstr;
    block_pattern::iterator mziter(blkstr);
    bool first=true;
    while (mziter.next(r,c,mzeig)) {
      if (!first) 
	ostr << ", ";
      switch (structureflags) {
      case PM_STRUCTURE: {
	const cmatrix& ablk(a(mzeig,size_t(0)));
	ostr << ablk.rows() << 'x' << ablk.cols();
      }
	break;
      case PM_EIGENBASIS:
	ostr << (r+1) << ',' << (c+1);
	break;
      default:
	throw InternalError("spyprint: unhandled structure output for BlockedOperator");
      }
      first=false;
    }
    if (structureflags==PM_STRUCTURE)
      ostr << " (sizes of matrix blocks)\n";
    else
      ostr << " (row, column eigenbasis block of associated Hamiltonian)\n";
    return;
  }
  
  //  if (precision)
  a.print(ostr,verbose & VER_GEN);
  //   else
  //     spy(ostr,a.row());
  dump_density(a.row(),ostr);
}

//! default - 0 precision not special 
template<class M> void spyprint(std::ostream& ostr,const M& a) {
  ostr << a;
}

MasterObj* masterobjp=NMRSIM_NULL;
bool istimedomain;

double get_rotor_phase(double t)
{
  if (!masterobjp)
    throw InternalError("get_rotor_phase: called before Hamiltonian created");
  const double lim360=1.0-1e-3*phasetolerance;
  
  try {
    const double x=masterobjp->phase(t)/TWO_PI;
    const double frac=x-floor(x);
    return (frac>lim360) ? 0.0 : 360.0*frac;
  } catch (MatrixException&) {
    return 0.0;
  }
}

inline double get_rotor_phase() { return get_rotor_phase(master_time); }

static timer<> global_stopwatch;
inline double get_cputime() { return global_stopwatch(); }

ContextWarning<> phasemodjitter_warning("phase modulation optimisation is likely to interact poorly with spin rate jitter and has been disabled (use -enable:phasemodulation to override",&NMRsim_once_warning);
ContextWarning<> jitteropts_warning("The following optimisations may be inappropriate in the presence of spinning rate jitter: combinepropagators, smartprop, gammacompute. Consider explicitly disabling them and/or explicitly enabling forcepointbypoint.",&NMRsim_once_warning);
ContextWarning<> excessivejitter_warning("Spinning rate jitter parameter exceeds 200 ns. This is unlikely to be simulated well.",&NMRsim_once_warning);
ContextWarning<> cortimeset_warning("Doesn't make sense to set correlation time for spin rate jitter when amplitude of jitter is zero",&NMRsim_once_warning);
void check_jitter()
{
  if (spin_rate_jitter<0)
    error_abort("spinning rate jitter cannot be <0");
  if (spin_rate_jitter_cortime<0)
    error_abort("spinning rate jitter correlation time cannot be <0");
  if (!nochecks && (spin_rate_jitter_cortime!=0.0) && (spin_rate_jitter==0.0))
    cortimeset_warning.raise();

  static bool donechecks=nochecks;
  if (donechecks || (spin_rate_jitter==0.0))
    return;

  if (optphasemod.isauto()) {
    phasemodjitter_warning.raise();
    optphasemod.set(option::OFF);
  }
    
  if ( optsmartprop.isauto() || optcombinepropagators.isauto() || optgamma.isauto() || !(optforcepointbypoint.isenabled()))
    jitteropts_warning.raise();

  if (spin_rate_jitter>0.2)
    excessivejitter_warning.raise();
}

bool build_jitter()
{
  static double last_jitter=-1e10;
  static double last_cortime=-1e10;
  if ((spin_rate_jitter!=last_jitter) || (spin_rate_jitter_cortime!=last_cortime)) {

    last_jitter=spin_rate_jitter;
    last_cortime=spin_rate_jitter_cortime;
    check_jitter();
    
    jitterinfo=IntervalSampler(spin_rate_jitter*1e-6,spin_rate_jitter_cortime*1e-6,randomise);
    if (spin_rate_jitter && (verbose & VER_GEN)) {
      std::cout << "Building spin rate jitter info using gaussian jitter with standard deviation of " << (spin_rate_jitter*1e3) << " ns";
      if (spin_rate_jitter_cortime)
	std::cout << " and correlation time of " << spin_rate_jitter_cortime << "us\n";
      else
	std::cout << '\n';
    }
  }

  return (spin_rate_jitter!=0.0);
}

// struct SysVar_chebyshev : public SystemVariable<size_t*> {
//   SysVar_chebyshev()
//     : SystemVariable<size_t*>("chebyshev_iterations",&(cmatrix_eigensystem_controller.chebyshev_iterations),false,true) {}
  
//   void update();
// };
  
// void SysVar_chebyshev::update() 
// {
//   if ((*this)()>eigensystem_controller::max_chebyshev_iterations) {
//     std::cerr << "Maximum Chebyshev iterations is " << eigensystem_controller::max_chebyshev_iterations <<'\n';
//     error_abort();
//   }
// }

SystemVariable<double*> v_spin_rate("spin_rate",&spin_rate,1.0,V_UPDATEINTS | V_ISCONST); //change to spin rate forces Hamiltonian re-build
SystemVariable<double*> v_spin_rate_jitter("spin_rate_jitter",&spin_rate_jitter,1.0,V_UPDATEINTS | V_ISCONST); //change to spin rate jitter forces Hamiltonian re-build
SystemVariable<double*> v_spin_rate_jitter_cortime("spin_rate_jitter_correlationtime",&spin_rate_jitter_cortime,1.0,V_UPDATEINTS | V_ISCONST); //change to spin rate jitter forces Hamiltonian re-build
SystemVariable<double*> v_rotor_angle("rotor_angle",&def_rotor_angle,rad_to_deg,V_UPDATEINTS | V_ISCONST);
SystemVariable<int*> v_gamma_angles("gamma_angles",&gamma_angles,V_ISCONST);
SystemVariable<double*> v_gammazero("gamma_zero",&gamma_zero,rad_to_deg,V_UPDATEINTS | V_ISCONST);
SystemVariable<double*> v_proton_freq("proton_frequency",&proton_freq,1.0,V_ISCONST | V_UPDATEINTS);
SystemVariable<double*> v_observe_freq("observe_frequency",&detect_freq,1.0,V_ISFIXED);
//SystemVariable<int*> v_sum_index("i_sum",&sum_index,true,false); //fixed, but not constant
SystemVariable<int*> v_evaluation_index("i_evaluation",&evaluation_index,V_ISFIXED);
SystemVariable<int*> v_sideband("sideband",&sideband,0);
SystemVariable<int*> v_eigenvalue("eigenvalue",&eigenvalue,0);
SystemVariable<int*> v_eigenvalues("eigenvalues",&eigenvalues,V_ISFIXED | V_ISCONST);
SystemVariable<size_t*> v_chebyshev("chebyshev_iterations",&(cmatrix_eigensystem_controller.chebyshev_iterations),V_ISCONST);
//SysVar_chebyshev v_chebyshev;
//SystemVariable<int*> v_ntimes("acqseq_ntimes",&acqseq_ntimes,false,true);
SystemVariable<double*> v_time("time",&master_time,1e6,V_ISFIXED); //output time in us
//SysVar_trans v_trans;
SystemVariable<FUNC_REAL> v_rotor_phase("rotor_phase",get_rotor_phase);
SystemVariable<FUNC_REAL> v_cputime("cputime",get_cputime);

bool isallowed_nonconst_var(SystemVariableBase* varp)
{
  return (varp==&v_evaluation_index);
}

bool proton_freq_isconstant()
{
  return v_proton_freq.isconstant();
}

BlockedOperator sigma0,detect;
LIST<nuclei_spec> blockingnuclei;

ThreadWarning<> syncfailed_warning("failed to find synchronisation: ",&NMRsim_repeat_warning);

void make_par_variables()
{
  if (!optsync)
    syncfailed_warning.type(BaseWarning::Ignore); //!< disable warnings if synchronisation explicitly disabled
  
  make_1d_par_variables();

  add_systemvarmap(v_spin_rate);
  add_systemvarmap(v_rotor_angle);
  add_systemvarmap(v_gamma);
  add_systemvarmap(v_gamma_angles);
  add_systemvarmap(v_observe_freq);
  add_systemvarmap(v_eigenvalue);
  add_systemvarmap(v_eigenvalues);
  add_systemvarmap(v_sideband);  
  //  add_systemvarmap(v_sum_index);
  add_systemvarmap(v_evaluation_index);

  command_Factory_t& par_Factory(get_par_Factory());
  par_Factory["maxdt"]=par_t(&parse_maxdt,true);
  par_Factory["maxjumpdt"]=par_t(&parse_maxjumpdt,true);
  par_Factory["steps_per_rotation"]=&parse_steps_per_rotation;
  par_Factory["spin_rate_jitter"]=&parse_spin_rate_jitter;
  //  par_Factory["operator"]=par_t(&parse_makeoperator,true);
  par_Factory["eigenvalue"]=&parse_eigenvalue;
  par_Factory["sideband"]=&parse_sideband;  
  par_Factory["propagation_method"]=&parse_prop_method;
  par_Factory["chebyshev_iterations"]=&parse_chebyshev_iterations;
  par_Factory["zerotolerance"]=&parse_zerotolerance;
}

template<class T> inline bool hasnonzero(const T& a) {
  const typename T::const_iterator aend=a.end();
  const typename T::value_type zero(0);
  return (std::find_if(a.begin(),aend,std::bind1st(std::not_equal_to<typename T::value_type>(),zero))!=aend);
}

ThreadWarning<> matrixspec::writefailed_warning("write to log file failed: object not defined at this point",&NMRsim_repeat_warning);

void matrixspec::log(logfile_controller& ctrl, const char* name, int flags) const
{
  const bool full = flags & PM_FULL;
  const int sflags= flags & PM_STRUCTUREFLAGS;
  try {
    switch (type) {
    case OPERATOR: {
	const BlockedOperator& op(asoperator().op); 
	if (full)
	  ctrl.write(op.full(),name);
	else {
	  if (sflags) {	    
	    const size_t eigblks=op.eigblocks();    
	    if (eigblks!=1)
	      error_abort("structure output not implemented for operators with eigenvalue structure e.g. periodic bases");

	    const block_pattern blkstr(op.blockstructure());
	    const size_t blks=blkstr.blocks;
	    List<size_t> rs(blks);
	    List<size_t> cs(blks);
	    
	    size_t r,c,mzeig;
	    size_t count=0;
	    block_pattern::iterator mziter(blkstr);
	    while (mziter.next(r,c,mzeig)) {
	      switch (sflags) {
	      case PM_STRUCTURE: {
		const cmatrix& ablk(op(mzeig,size_t(0)));
		rs(count)=ablk.rows();
		cs(count)=ablk.cols();
	      }
		break;
	      case PM_EIGENBASIS:
		rs(count)=r+1;
		cs(count)=c+1;
		break;
	      default:
		throw InternalError("matrixspec::log: unhandled structure output for BlockedOperator");
	      }
	      count++;
	    }
	    if (count!=blks)
	      throw Mismatch("matrixspec::log: number of coherence blocks found does not match expected",count,blks);

	    char tmpname[128];
	    const char* type=(sflags==PM_STRUCTURE) ? "sizes" : "blockindices";
	    snprintf(tmpname,sizeof(tmpname),"%s_row_%s",name,type);
	    ctrl.write(rs,tmpname);
	    snprintf(tmpname,sizeof(tmpname),"%s_column_%s",name,type);
	    ctrl.write(cs,tmpname);
	  }
	  else
	    ctrl.write(op.row(),name);
	}
    }
      break;
//     case DOUBLE:
//       ctrl.write(asdouble(),name);
//       break;
    case BOOL:
      if (sflags)
	error_abort("-structure / -eigenbasis only implemented for Hamiltonian / operator matrices");
      ctrl.write(asfilter()().row(),name);
      break;
    default:
      throw InternalError("matrixspec: can only log operator or filter matrices");
    }
  } catch (Undefined&) {
    writefailed_warning.raise();
  }
}

std::ostream& operator<< (std::ostream& ostr, const matrixspec& curm)
{
  switch (curm.type) {
  case matrixspec::SPECIAL:
    return ostr << "Can't output here (local)\n";
  case matrixspec::BOOL: {
    const BlockedFilter* filterp(curm.asfilter().filterp);
    if (filterp)
      return ostr << *filterp;
    else
      return ostr << "<undefined here>\n";
  }
  case matrixspec::OPERATOR:    
    spyprint(ostr,curm.asoperator().op);
    return ostr;
//   case matrixspec::DOUBLE:
//     spyprint(ostr,curm.asdouble());
//     return ostr;
  case matrixspec::EXCHANGE:
    return ostr << curm.asexchange();
  case matrixspec::NONE:
    return ostr << "<undefined>\n";
  }
  throw InternalError("Unknown matrix type");
}

usage_t matrixspec::usage() const
{
  switch (type) {
  case BOOL:
    return usage_(asfilter()());
  case OPERATOR:
    return usage_(asoperator().op.row());
    //  case DOUBLE:
    //return usage_(asdouble());
  case EXCHANGE:
    return usage_(asexchange()());
  default:
    break;
  }
  return usage_t(); //!< 0 for special types    
}

matrixspec& matrixspec::operator= (const matrixspec& from)
{
  switch (from.type) {
  case BOOL:
    set_filter(from.asfilter());
    return *this;
  case OPERATOR:
    set_operator(from.asoperator());
    return *this;
//   case DOUBLE:
//     set_double(from.asdouble());
//     return *this;
  case NONE:
    type=NONE;
    ptr=NMRSIM_NULL;
    return *this;
  case EXCHANGE:
    set_exchange(from.asexchange());
    return *this;
  case SPECIAL:
    type=SPECIAL;
    ptr=from.ptr;
    return *this;
  default:
    throw Undefined("matrixspec=");
  }
}      

bool matrixspec::operator!() const 
{
  switch (type) {
  case OPERATOR:
    return !(asoperator().op);
  case BOOL:
    return !asfilter();
  case EXCHANGE:
    return !(asexchange()());
    //  case DOUBLE:
    // return static_cast< const BlockedMatrix<double>* >(ptr)->empty();
  case SPECIAL:
    return false;
  default:
    break;
  }
  return true;
}

void matrixspec::set_operator(const operator_def& a)
{
  type=OPERATOR;
  ptr=static_cast<const void*>(&a);
}

// void matrixspec::set_double(const BlockedMatrix<double>& a)
// {
//   type=DOUBLE;
//   ptr=static_cast<const void*>(&a);
// }

void matrixspec::set_exchange(const ExchangeMatrix& a)
{
  type=EXCHANGE;
  ptr=static_cast<const void*>(&a);
}

void matrixspec::set_filter(const filter_def& a)
{
  type=BOOL;
  ptr=static_cast<const void*>(&a);
}

const operator_def& matrixspec::asoperator() const
{
  if (type!=OPERATOR)
    throw Failed("matrixspec::asoperator: not an operator matrix");
  return *static_cast< const operator_def* >(ptr);
}

const ExchangeMatrix& matrixspec::asexchange() const
{
  if (type!=EXCHANGE)
    throw Failed("matrixspec::asexchange: not an exchange matrix");
  return *static_cast< const ExchangeMatrix* >(ptr);
}

// const BlockedMatrix<double>& matrixspec::asdouble() const 
// {
//   if (type!=DOUBLE)
//     throw Failed("matrixspec::ascomplex: not an rmatrix");
//   return *static_cast< const BlockedMatrix<double>* >(ptr);
// }

const filter_def& matrixspec::asfilter() const
{
  if (type!=BOOL)
    throw Failed("matrixspec::ascomplex: not a filter matrix");
  return *static_cast<const filter_def*>(ptr);
}

// template<class R,class C> matrixspec deref(const RealComplexHolder<R,C>& a)
// {
//   typedef RealComplexHolder<R,C> holder_t;
//   switch (a.type()) {
//   case holder_t::REAL:
//     return matrixspec((*(a.get_real()))());
//   case holder_t::COMPLEX:
//     return matrixspec((*(a.get_complex()))());
//   }
//   throw Undefined("dereference Hamiltonian");
// }

template<class R,class C> bool isempty(const RealComplexHolder<R,C>& a)
{
  typedef RealComplexHolder<R,C> holder_t;
  switch (a.type()) {
  case holder_t::REAL:
    return !(*(a.get_real()));
  case holder_t::COMPLEX:
    return !(*(a.get_complex()));
  }
  return true;
}

bool MasterObj::Hempty() const
{
  switch (hamobjtype) {
  case M_STATIC:
    return isempty(Hstaticp);
  case M_STATICD:
    return !(*Hstaticdp);
  case M_SPIND:
    return !(*Hspindp);
  case M_SPIN:
    return isempty(Hspinp);
  default:
    throw InternalError("Hempty");
  }
}

usage_t MasterObj::usage() const
{
  return Husage();
}

ThreadWarning<> nullhamiltonian_warning("Usage statistics requested for undefined Hamiltonian",&NMRsim_once_warning);

usage_t MasterObj::Husage() const
{
  switch (hamobjtype) {
  case M_STATIC:
    switch (Hstaticp.type()) {
    case Hstatic_t::REAL:
      return usage_t((*(Hstaticp.get_real()))().row().size(),Type2Type<double>());
    case Hstatic_t::COMPLEX:
      return usage_t((*(Hstaticp.get_complex()))().row().size(),Type2Type<complex>());
    }
    break;
  case M_STATICD:
    return usage_t(Hstaticdp->usage(),Type2Type<double>());
  case M_SPIND:
    return usage_t(Hspindp->usage(),Type2Type<double>());
  case M_SPIN:
    switch (Hspinp.type()) {
    case Hspin_t::REAL:
      return usage_t((Hspinp.get_real())->usage(),Type2Type<double>());
    case Hspin_t::COMPLEX:
      return usage_t((Hspinp.get_complex())->usage(),Type2Type<complex>());
    }
  case M_NULL:
    nullhamiltonian_warning.raise();
    return usage_t();

  default:
    break;
  }
  throw InternalError("Unknown hamobjtype");
}


void MasterObj::printH(std::ostream& ostr, int structureflags) const
{
  if (structureflags) {
    const SpinOpGeneratorBase& opgen(spinop_generator());
    const size_t eigblks=opgen.eigblocks();
    switch (structureflags) {
    case PM_EIGENBASIS:
      if (eigblks!=1)
	error_abort("cannot currently output Hamiltonian eigenbasis in this case");    
      for (size_t i=0;i<opgen.actual_mzblocks();i++) {
	const BaseList<size_t> curinds(opgen.blockindices(i));
	ostr << "State order for block " << (i+1) << ": " << curinds << '\n';
      }
      break;
    case PM_STRUCTURE: {
      const BaseList<size_t>& diagsizes(opgen.diagonal_structure());
      size_t norder=0;
      for (size_t i=0;i<diagsizes.size();i++) {
	if (i)
	  ostr << ", ";
	const size_t n=diagsizes(i);
	ostr << n << 'x' << n;
	norder+=n*n*n;
      }
      ostr << " (block sizes)  sum n^3: " << norder << '\n';
    }
      break;
    default:
      throw InternalError("printH: unhandled structure type");
    }
    return;
  }

  switch (hamobjtype) {
  case M_STATIC:
    switch (Hstaticp.type()) {
    case Hstatic_t::REAL:
      spyprint(ostr,*(Hstaticp.get_real()));
      break;
    case Hstatic_t::COMPLEX:
      spyprint(ostr,*(Hstaticp.get_complex()));
      break;
    }
    break;
  case M_STATICD:
    ostr << (*Hstaticdp);
    break;
  case M_SPIND:
    spyprint(ostr,(*Hspindp)(master_time));
    break;
  case M_SPIN:
    switch (Hspinp.type()) {
    case Hspin_t::REAL:
      spyprint(ostr,(*(Hspinp.get_real()))(master_time));
      break;
    case Hspin_t::COMPLEX:
      spyprint(ostr,(*(Hspinp.get_complex()))(master_time));
      break;
    default:
      ostr << "<unset>\n";
    }
    break;
  default:
    throw Failed("Unknown hamobjtype");
  }
  if (verbose & VER_PROFILE)
    ostr << Husage();
}

void MasterObj::printfullH(std::ostream& ostr) const
{
  if (!spin_rate) {
    printH(ostr);
    return;
  }
  switch (hamobjtype) {
  case M_SPIND:
    ostr << *Hspindp;
    break;
  case M_SPIN:
    switch (Hspinp.type()) {
    case Hspin_t::REAL:
      ostr << *(Hspinp.get_real());
      break;
    case Hspin_t::COMPLEX:
      ostr << *(Hspinp.get_complex());
      break;
    default:
      ostr << "<unset>\n";
    }
    break;
  default:
    throw Failed("Unknown hamobjtype");
  }
  if (verbose & VER_PROFILE)
    ostr << Husage();
}

namespace {
  void write_diag(logfile_controller& ctrl, const ListList<double>& Hdiag, const char* name)
  {
    BlockedMatrix<double> Hfull;
    full(Hfull,Hdiag);
    ctrl.write(Hfull,name);
  }
}

void MasterObj::logH(logfile_controller& ctrl, const char* name, int structureflags) const
{
  if (structureflags) {
    const SpinOpGeneratorBase& opgen(spinop_generator());
    const size_t eigblks=opgen.eigblocks();
    switch (structureflags) {
    case PM_EIGENBASIS: {
      if (eigblks!=1)
	error_abort("cannot currently output Hamiltonian eigenbasis in this case"); 
      ListList<size_t> mzblks; //!< OK, slightly inefficient, but output slow anyway
      const size_t nblks=opgen.actual_mzblocks();
      for (size_t i=0;i<nblks;i++)
	mzblks.push_back(opgen.blockindices(i));
      ctrl.write(mzblks,name);
    }
      break;
    case PM_STRUCTURE:
      ctrl.write(opgen.diagonal_structure(),name);
      break;
    default:
      throw InternalError("logH: unhandled structure type");
    }
    return;
  }

  switch (hamobjtype) {
  case M_STATIC:
    switch (Hstaticp.type()) {
    case Hstatic_t::REAL:
      ctrl.write((*(Hstaticp.get_real()))(),name);
      break;
    case Hstatic_t::COMPLEX:
      ctrl.write((*(Hstaticp.get_complex()))(),name);
      break;
    }
    return;
  case M_STATICD: 
    write_diag(ctrl,(*Hstaticdp)(),name);
    return;
  case M_SPIND: {
    const ListList<double> Hdiag((*Hspindp)(master_time));
    write_diag(ctrl,Hdiag,name);
  }
    return;
  case M_SPIN:
    switch (Hspinp.type()) {
    case Hspin_t::REAL:
      ctrl.write((*(Hspinp.get_real()))(master_time),name);
      return;
    case Hspin_t::COMPLEX:
      ctrl.write((*(Hspinp.get_complex()))(master_time),name);
      return;
    }
  default:
    throw Failed("Unknown/undefined Hamiltonian");
  }
}

template<class OpGen> void MasterObj::create_Hstatic(const OpGen& opgen, Type2Type<complex>)
{
  (Hstaticp.set_complex()).reset(new BlockedStaticHamiltonian<complex>(opgen));
}

template<class OpGen> void MasterObj::create_Hstatic(const OpGen& opgen, Type2Type<double>)
{
  (Hstaticp.set_real()).reset(new BlockedStaticHamiltonian<double>(opgen));
}

template<class OpGen> void MasterObj::create_Hspin(const OpGen& opgen, double gamma,Type2Type<complex>)
{
  if (build_jitter())
    (Hspinp.set_complex()).reset(new BlockedSpinningHamiltonian<complex>(opgen,spin_rate,gamma,rinfo,jitterinfo));//,CWdesc));
  else
    (Hspinp.set_complex()).reset(new BlockedSpinningHamiltonian<complex>(opgen,spin_rate,gamma,rinfo));//,CWdesc));
}

template<class OpGen> void MasterObj::create_Hspin(const OpGen& opgen, double gamma,Type2Type<double>)
{
  if (build_jitter())
    (Hspinp.set_real()).reset(new BlockedSpinningHamiltonian<double>(opgen,spin_rate,gamma,rinfo,jitterinfo));
  else
    (Hspinp.set_real()).reset(new BlockedSpinningHamiltonian<double>(opgen,spin_rate,gamma,rinfo));
}

bool MasterObj::reset_angle(double angle)
{
  if (spin_rate) {
    const size_t rank= Hstructp->isclassicsecondorder() ? 4 : 2;
    if ((rank!=lastrank) || (angle!=rotor_angle)) {
      if ((verbose & VER_GEN) && (verbose_level>1))
	std::cout << "Changing/setting rotor angle to " << (angle*rad_to_deg) << " degrees (rank " << rank << ")\n";
      rinfo=RotorInfo(rank,angle,rotor_orient);
      rotor_angle=angle;
      lastrank=rank;
      return true;
    }
  }
  else
    rotor_angle=0.0;
  return false;
}
 
bool MasterObj::change_angle(double angle)
{
  const bool changed_angle=reset_angle(angle);
  if (changed_angle)
    recreate_Hs();
  return changed_angle;
}

template<class HType,class StoreType> void set_interactions(HType& H, const Euler& powder, Type2Type<StoreType>)
{
  HamiltonianStore<StoreType> Hstore(*interactions_MFp,powder);
  if (verbose & VER_GEN)
    std::cout << Hstore;
  H.interactions(Hstore);
}    

template<class StaticObj,class SpinObj> void MasterObj::set_interactions(const Euler& powder, StaticObj& Hstaticp, SpinObj& Hspinp)
{
  switch (hamobjtype) {
  case M_SPIND:
    // Need to trap empty Hamiltonian otherwise propagation will fail
    if (interactions_MFp->empty())
      error_abort("Inhomogeneous Hamiltonian but no interactions - what is there to simulate?");
    ::set_interactions(*Hspindp,powder,Type2Type<space_T>());
    break;
  case M_SPIN:
    ::set_interactions(*Hspinp,powder,Type2Type<space_T>());
    break;
  case M_STATICD:
    if (quadrupole_order!=1)
      ::set_interactions(*Hstaticdp,powder,Type2Type<space_T>());
    else 
      ::set_interactions(*Hstaticdp,powder,Type2Type<double>());    
    break;
  case M_STATIC:
    if (quadrupole_order!=1)
      ::set_interactions(*Hstaticp,powder,Type2Type<space_T>());
    else {
      ::set_interactions(*Hstaticp,powder,Type2Type<double>());
      set_Hsystem((*Hstaticp)());
    }
    break;
  default:
    throw Failed("Unknown hamobjtype");
  }

  if (verbose & VER_GEN) {
    std::cout << "Full Hamiltonian\n";
    printfullH(std::cout);
  }
}

/**
   \pre Hamiltonian can not be static (otherwise \c Failed exception)
*/
void MasterObj::rotor_phase(double phase0)
{
  switch (hamobjtype) {
  case M_SPIN:
    if (Hspinp.iscomplex())
      (Hspinp.get_complex())->rotor_phase(phase0);
    else
      (Hspinp.get_real())->rotor_phase(phase0);
    break;
  case M_SPIND:
    Hspindp->rotor_phase(phase0);
    break;
  default:
    throw Failed("rotor_phase called for static Hamiltonian!");
  }
}

/**
   \pre Hamiltonian can not be static (otherwise \c Failed exception)
*/
double MasterObj::phase(double t) const
{
  switch (hamobjtype) {
  case M_SPIN:
    return (Hspinp.iscomplex()) ?
      (Hspinp.get_complex())->phase(t) : (Hspinp.get_real())->phase(t);
  case M_SPIND:
    return Hspindp->phase(t);

  default:
    throw Failed("phase called for static Hamiltonian!");
  }
}

int check_sync(double nfloat,double tol)
{
  if (nfloat<=0.0) {
    char buf[256];
    snprintf(buf,sizeof(buf),"check_sync: asked to round to an integer a number (%g) which is <=0",nfloat);
    throw InvalidParameter(buf);
  }
  if (!optsync)
    return 0;
  int n=int(nfloat+0.5);
  return (fabs(nfloat-n)>tol) ? 0 : n;
}

bool MasterObj::canbuildoperators(const BaseList<productoperator_spec*>& specs) const
{
  for (size_t i=specs.size();i--;) {
    const productoperator_spec& curspec(*(specs(i)));
    try {
#ifndef NOPERIODIC
      if (!simple_opgenp)
	block_pattern bp(*crystal_opgenp,curspec);
      else
#endif
	block_pattern bp(*simple_opgenp,curspec);

    } catch (Failed&) {
      if (verbose & VER_GEN)
	parser_printthread(std::cerr) << "Failed to build spin operator: " << curspec << '\n';
      return false;
    }
  }
  return true;
}

MasterObj::MasterObj()
{
}

ThreadWarning<> commonPAS_warning("Optimisation for special case of common PAS is not implemented",&NMRsim_once_warning);

MasterObj::MasterObj(const basespin_system& sys,const HamiltonianStore<space_T>& Hstorev,const BaseList<nuclei_spec>& blockingnuc,int flags, bool use_crystal)
    : 
  Hstructp(new HamiltonianStructure(sys,Hstorev,blockingnuc,weakints,flags)), 
  simple_opgenp(use_crystal ? NMRSIM_NULL : new SpinOpGenerator(*Hstructp,*cstructp,(verbose & VER_GEN) ? verbose_level : 0))
#ifndef NOPERIODIC
  ,crystal_opgenp(use_crystal ? new CrystalOpGenerator(*Hstructp,*cstructp,(verbose & VER_GEN) ? verbose_level : 0) : NMRSIM_NULL)
#endif
  ,hamtype(H_NULL), hamobjtype(M_NULL)
{
  if (!canbuildoperators(opspecstack))
    throw Failed("Can't build spin operators");

  update();

  if (!have_spinsys()) 
    return; 
  //    restart();
  //  return;
  //}
  
  if (verbose & VER_GEN) {
    std::cout << "Hamiltonian structure:\n" << *Hstructp;
#ifndef NOPERIODIC
    if (!simple_opgenp)
      std::cout << (*crystal_opgenp);
    else
#endif
      std::cout << (*simple_opgenp);
  }

//   if (Hstructp->isdiagonal())
//     ishomogeneous=false;
  
//   if (verbose & VER_GEN) {
//     std::cout << "Inhomogeneous: "; 
//     if (ishomogeneous)
//       std::cout << "No\n";
//       else {
// 	if (iscommonPAS)
// 	  std::cout << "Yes (common PAS)\n";
// 	else
// 	  std::cout << "Yes\n";
//       }
//   }
  
  if (sw && (verbose & VER_GEN)) {
    std::cout << "Acquisition time: ";
    prettyprint_time(np/sw) << '\n';
  }
  
  quadrupole_order=Hstructp->quadrupole_order();

  if (verbose & VER_GEN)
    std::cout << "Overall quadrupole order: " << quadrupole_order << '\n';
  
  ensure_channels();
  //  restart();
}

//! these warnings output on cout rather than cerr
ThreadWarning<> generalisedqpole_warning("The generalised quadrupole treatment has not been fully evaluated.  Use at your own risk!",&NMRsim_once_warning,BaseWarning::Inherit,std::cout);
ThreadWarning<> qpole_rf_warning("RF sequences have been defined with non-first-order quadrupoles present.  Treatment is only approximate and has not been fully evaluated.  Use at your own risk!",&NMRsim_once_warning,BaseWarning::Inherit,std::cout);

void MasterObj::update()
{
  if (!have_spinsys())
    return;

  //!< can't do this in restart since havestore is not set until par block (after MasterObj created)
  if (quadrupole_order!=1) {
    if (havestore)
      qpole_rf_warning.raise();
    if (!isclassicQ())
      generalisedqpole_warning.raise();
  }

  if (proton_freq && !!simple_opgenp)
    simple_opgenp->proton_frequency(proton_freq);
  eigenvalues=spinop_generator().eigblocks();
}

BlockedOperator MasterObj::make(const productoperator_spec& spec) const
{
#ifndef NOPERIODIC
  if (!simple_opgenp)
    return BlockedOperator(*crystal_opgenp,spec);
#endif
  return BlockedOperator(*simple_opgenp,spec);
}


void MasterObj::restart()
{
  timer_guard_t guard(*setup_cp,true); //!< ignore if already timing

// //   if (phasep)
// //     phasep->update();
  update();

// //   if ((verbose & VER_GEN) && acqseqp)
// //     std::cout << "Acquisition sequence\n" << (*acqseqp) << '\n';

  if (have_spinsys())
    reset_angle(def_rotor_angle);
//   //  if (acqp)
//   //  acqp->checksync();
}

#ifndef NOPERIODIC
void MasterObj::recreate_Hs(const CrystalOpGenerator& opgen, double gamma)
{
  switch (hamobjtype) {
  case M_STATIC:
    create_Hstatic(opgen,Type2Type<complex>());
    break;
  case M_SPIN:
    create_Hspin(opgen,gamma,Type2Type<complex>());
    break;
  default:
    throw Failed("Unknown/unsupported hamobjtype");
  }
}
#endif

void MasterObj::recreate_Hs(const SpinOpGenerator& opgen, double gamma)
{
  switch (hamobjtype) {
  case M_STATIC:
    if (hamtype==H_COMPLEX)
      create_Hstatic(opgen,Type2Type<complex>());
    else
      create_Hstatic(opgen,Type2Type<double>());
    break;
  case M_SPIN:
    if (hamtype==H_COMPLEX)
      create_Hspin(opgen,gamma,Type2Type<complex>());
    else
      create_Hspin(opgen,gamma,Type2Type<double>());
    break;
  case M_SPIND:
    if (spin_rate_jitter)
      throw Failed("Spinning rate jitter not implemented for diagonal Hamiltonians"); //!< should have been caught earlier
    Hspindp.reset(new BlockedDiagonalSpinningHamiltonian(opgen,spin_rate,gamma,rinfo));
    break;
  case M_STATICD:
    Hstaticdp.reset(new BlockedDiagonalStaticHamiltonian(opgen));
    break;
  default:
    throw Failed("Unknown/unset hamobjtype");
  }
}

ThreadWarning<> diagonaljitter_warning("Spinning rate jitter not currently implemented for diagonal Hamiltonians - reverting to full Hamiltonian",&NMRsim_once_warning);

/**
   \post ::MasterObj::Hempty is \c false
*/
void MasterObj::recreate_Hs()
{
  static hamobj_t lasthamobjtype=M_NULL;
  static ham_t lasthamtype=H_NULL;

  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Recreating Hamiltonian\n";

  //const CompSequenceBase* acqseqp = acqp ? acqp->sequence() : NMRSIM_NULL
  hamtype=H_COMPLEX;
  hamobjtype=M_NULL;

  //  if (!acqseqp && !noshortcuts && !haveactions) { //can't optimise if RF present
  const bool allowreal=!(Hstructp->iscomplex()) && optrealH();
  optrealH.check(allowreal);
  if (allowreal) {
    //    const bool allowdiagonal=!acqseqp && !haveactions && !!(phasep->simple_opgenp); //CrystalOpGen doesn't support diagonal 
    //    if (Hstruct.isdiagonal() && !hasacquisition_rf() && !haveactions && (quadrupole_order!=0))
    if (Hstructp->isdiagonal() && !haveactions && ((quadrupole_order==1) || ((quadrupole_order!=0) && !isclassicQ()))) {
      if (spin_rate_jitter) {
	diagonaljitter_warning.raise();
	hamtype=H_REAL; // consolation prize
      }
      else
	hamtype=H_DIAGONAL;
    }
    else {
      if ((quadrupole_order==1) || isclassicQ())
	hamtype=H_REAL;
    }
  }

  if (hamtype!=lasthamtype) {
    if (verbose & VER_GEN) {
      std::cout << "Hamiltonian type: ";
      switch (hamtype) {
      case H_DIAGONAL:
	std::cout << "diagonal\n";
      break;
      case H_COMPLEX:
	std::cout << "complex\n";
	break;
      case H_REAL:
	std::cout << "real\n";
	break;
      default:
	throw InternalError("Hamiltonian type");
      }
    }
    lasthamtype=hamtype;
  }

  if (spin_rate) {
    if (lastrank==0)
      throw InternalError("recreate_Hs called before RotorInfo established");

    switch (hamtype) {
    case H_DIAGONAL:
      hamobjtype=M_SPIND;
      break;
    default:
      if (iscommonPAS && !haveactions)
	commonPAS_warning.raise();
      hamobjtype=M_SPIN;
      break;
    }
  }
  else {
    switch (hamtype) {
    case H_DIAGONAL:
      hamobjtype=M_STATICD;
      break;
    case H_COMPLEX: case H_REAL:
      hamobjtype=M_STATIC;
      break;
    default:
      throw InternalError("Hamiltonian type");
    }
  }

#ifndef NOPERIODIC
  if (!simple_opgenp)
    recreate_Hs(*crystal_opgenp,global_powder.gamma);
  else
#endif
    recreate_Hs(*simple_opgenp,global_powder.gamma);
  
  const Euler powder_no_gamma(global_powder.alpha,global_powder.beta,0.0); //!< gamma angle is set in rotor phase
  if (iscomplex())
    set_interactions(powder_no_gamma,Hstaticp.get_complex(),Hspinp.get_complex());
  else
    set_interactions(powder_no_gamma,Hstaticp.get_real(),Hspinp.get_real());  

  if (lasthamobjtype!=hamobjtype) {
    verify_powderaverage(); //only verify powder average if changing type of Hamiltonian
    lasthamobjtype=hamobjtype;
  }
}

void parse_eigenvalue()
{
  parse_system_variable(v_eigenvalue);
}

void parse_sideband()
{
  parse_system_variable(v_sideband);
}

ContextWarning<> zerotolerance_warning("zerotolerance should be a small value",&NMRsim_once_warning);

void parse_zerotolerance()
{
  if (are_left()) {
    zerotolerance=parse_double();
    if (zerotolerance<0.0)
      error_abort("zerotolerance must be >=0.0");
    if (zerotolerance>ZEROTOLERANCE_TRIGGER)
      zerotolerance_warning.raise();
  }
  else
    std::cout << "zerotolerance=" << zerotolerance << '\n';
}
    
void parse_prop_method()
{
  const char* meth=parse_string(F_REPLACEDOLLAR);
//   if (strcmp(meth,"automatic")==0) {
//     prop_method=PROP_AUTO;
//     if (are_left()) {
//       chebyshev_crossover=0.01*parse_double();
//       if ((chebyshev_crossover<0.0) || (chebyshev_crossover>1.0))
// 	error_abort("chebyshev crossover density should be between 0 and 100 (%)\n");
//       if (allowwarnings() && (chebyshev_crossover!=0.0) && (chebyshev_crossover<0.01))
// 	std::cerr << "Warning: chebyshev crossover less than 1%  Value specified as fraction rather than percentage? (-nochecks disables this warning)\n";
//     }
//   }
//   else {
    if (strcmp(meth,"diagonalise")==0)
      prop_method=PROP_DIAG;
    else {
    if (strcmp(meth,"chebyshev")==0)
      prop_method=PROP_CHEBYSHEV;
    else
      error_abort("Unknown propagation method.  Known: diagonalise, chebyshev");
    }
    //  }
}

ContextWarning<> chebyshev_iterations_warning("chebyshev_iterations set to less than " NMRSIM_STRINGISE(CHEBYSHEV_TRIGGER) "; this is unlikely to work properly",&NMRsim_once_warning);

void parse_chebyshev_iterations()
{
  if (!are_left())
    std::cout << "chebyshev_iterations=" << v_chebyshev() << '\n';
  else {
    parse_system_variable(v_chebyshev);
    if (v_chebyshev()<CHEBYSHEV_TRIGGER)
      chebyshev_iterations_warning.raise();
  }
}

void parse_spin_rate()
{
  parse_system_variable(v_spin_rate,F_DENYZERO);
}

void parse_spin_rate_jitter()
{
  static const char* synstr="<jitter amplitude (us)>#[<correlation time (us)>]";

  parse_system_variable_syntax(synstr,1,v_spin_rate_jitter,F_DENYSUM);

  if (are_left())
    parse_system_variable_syntax(synstr,2,v_spin_rate_jitter_cortime,F_DENYSUM);
  else
    spin_rate_jitter_cortime=0.0;

  check_jitter();
}

void parse_gamma_zero()
{
  parse_system_variable(v_gammazero);
}

tstep_t tstep_type=INT_UNSET;
int steps_per_rotation=0;

double get_maxdt()
{
  if (tstep_type==INT_ROTORPERIOD) {
    if (spin_rate==0)
      error_abort("integration timestep cannot be specified by steps_per_rotation in static simulation");
    return (1.0/spin_rate)/steps_per_rotation;
  }
  return maxdtv;
}

void set_inttype(tstep_t newtype)
{
  if (tstep_type && (tstep_type!=newtype))
    error_abort("cannot mix steps_per_rotation and maxdt when specifying integration time step");
  tstep_type=newtype;
}

ContextWarning<> maxdt_warning("maxdt set to >" NMRSIM_STRINGISE(MAXDT_TRIGGER) " us - are you sure about this?",&NMRsim_once_warning);
ContextWarning<> parameteroverwrite_warning("parameter (maxdt/maxjumpdt) set here overrides earlier setting",&NMRsim_repeat_warning); 
ContextWarning<> jumpdt_warning("Use maxjumpdt with care e.g. check convergence of results with respect to parameter (-nochecks disables this warning)",&NMRsim_once_warning);

double parsedt(const char* parname, bool& doneset)
{
  if (doneset)
    parameteroverwrite_warning.raise();
  else
    doneset=true;
  const double dt=parse_double()*1e-6;
  if (dt<=0.0) {
    parser_printcontext() << parname << " must be >0\n";
    error_abort();
  }
  return dt;
}

void parse_maxdt()
{
  static bool doneset=false;
  if (!are_left())
    std::cout << "maxdt=" << (maxdtv*1e6) << " us\n";
  else {
    set_inttype(INT_EXPLICIT);
    maxdtv=parsedt("maxdt",doneset);
    if (maxdtv>MAXDT_TRIGGER)
      maxdt_warning.raise();
  }
}

void parse_steps_per_rotation()
{
  set_inttype(INT_ROTORPERIOD);
  steps_per_rotation=parse_int();
  if (steps_per_rotation<1)
    error_abort("steps_per_rotation must be >0");
}

void parse_maxjumpdt()
{
  static bool doneset=false;
  if (!are_left())
    std::cout << "maxjumpdt=" << (maxjumpdt*1e6) << " us\n";
  else {
    maxjumpdt=parsedt("maxjumpdt",doneset);
    if (!nochecks && maxjumpdt)
      jumpdt_warning.raise();
  }
}

void parse_gamma_angles()
{
  //  static flagsmap_type flags;
//   if (flags.empty())
//     flags["forcegcompute"]=1;
  parse_system_variable(v_gamma_angles);
  //  forcegcompute=(parse_flags(flags)!=0);
}

ContextWarning<> zeroangle_warning("rotor_angle of zero makes little sense!",&NMRsim_once_warning);

void parse_rotor_angle()
{
  //  if (are_left()) {
  parse_system_variable(v_rotor_angle);
  if (def_rotor_angle==0.0)
    zeroangle_warning.raise();
  //}
  //else
  //  std::cout << "rotor_angle=" << (def_rotor_angle*rad_to_deg) << " degrees\n";
}


void parse_generalmatrix(const char* name, bool isgeneral =true)
{
  //  const size_t nargs=count_left();
  //if (nargs<2)
  //  error_abort("Too few arguments to matrix set <name> exchange <real matrix elements> [<imaginary elements>]\n");
  //const size_t sites=(nargs>2) ? nargs-1 : 2;
  Mark markobj;
  Variable cvar(S_ARG1);
  VariableBase* elsp=parse_double_variable(cvar,F_ALLOWLIST);
  bool isconst=elsp->isconst();
  //  ExchangeMatrix* matrixp=NMRSIM_NULL;
//   if (are_left()) {
//     if (!isgeneral)
//       error_abort("cannot specify imaginary component to exchange matrix (use matrix set general to create a general complex (square) matrix)");
//     cvar.subsid=S_ARG2;
//     VariableBase* imagp=parse_double_variable(cvar,F_ALLOWLIST);
//     if (!(imagp->isconst()))
//       isconst=false;
//     matrixp=new ExchangeMatrix(name,realp->value(),imagp->value(),isconst);
//   }
//   else
  ExchangeMatrix* matrixp=new ExchangeMatrix(name,elsp->value(),isconst,isgeneral);
  markobj.flush(matrixp);
  matrixmap[name]=matrixspec(*matrixp);
}

bool matchtag(const char* name)
{
  const char* ptr=get_curline();
  const size_t len=strlen(name);
  if (strncmp(ptr,name,len)!=0)
    return false;
  const char term=ptr[len];
  const bool found=(term=='\0') || isspace(term);
  if (found)
    (void)get_token(); //!< swallow token
  return found;
}
  
void parse_matrix()
{
  if (strcmp(parse_string(F_REPLACEDOLLAR),"set")!=0)
    error_abort("syntax: matrix set <to> <from>");

  const char* matname= parse_string(F_REPLACEDOLLAR);

  if (strcmp(matname,"detect")==0)
    error_abort("use detect_operator to set detection operator");
  if (strcmp(matname,"start")==0)
    error_abort("use start_operator to set initial density matrix");

  if (matrixmap.count(matname))
    error_abort("can't redefine matrices");

  if (matchtag("totalcoherence")) {
    parse_ordermatrix(matname,true,true);
    return;
  }

  if (matchtag("coherenceorder")) {
    parse_ordermatrix(matname,true);
    return;
  }
  
  if (matchtag("spinorder")) {
    parse_ordermatrix(matname,false);
    return;
  }

//   if (matchtag("exchange")) {
//     parse_generalmatrix(matname,false);
//     return;
//   }
  if (matchtag("general")) {
    parse_generalmatrix(matname);
    return;
  }
  
  setableoperator_spec* specp=create_setableoperator(get_curline());
  BlockedOperator* opexprp=new BlockedOperator();
//   if (specp->isconstant()) {    
//     *opexprp=masterobjp->make(specp->current());
//     delete specp; //!< can now delete definition
//     specp=NMRSIM_NULL;
//   }
  operator_def* defp=new operator_def(specp,*opexprp);
  matrixmap[matname]=matrixspec(*defp);
  set_curline(NMRSIM_NULL);
}
   
//! ensure all operator matrices are up-to-date
/** Not very clever for non-const operators - these are remade each time
    rather than being cached for one powder orienantion to the next.
    Unlikely that this will be time limiting **/

void ensure_operator_matrices(bool constonly)
{
  if (!constonly)
    ensure_sigma0detect(); //!< sigma0/detect are special case
  const matrixmap_type::iterator end(matrixmap.end());
  matrixmap_type::iterator start(matrixmap.begin());
  while (start!=end) {
    matrixspec& mspec(start->second);
    if (!constonly || mspec.isconstant())
      mspec.update();
    ++start;
  }
}

const BlockedFilter& filter_def::operator()() const
{
  filter_def& asmut=const_cast< filter_def& >(*this);
  if (!(asmut.ensure()))
    throw Undefined("filter_def");
  return *filterp;
} 

ThreadWarning<> zerofilter_warning("Filter matrix created that is all zero - probably not intended!",&NMRsim_once_warning);
ThreadWarning<> zerofilterhermitian_warning("All zero filter matrix created. Filter matrices may not work as expected for blocked problems - try using -disable:mzblocking",&NMRsim_repeat_warning);

bool filter_def::ensure(bool needed)
{
  const char* sourcename=NMRSIM_NULL;
  const BlockedOperator* sourcep=NMRSIM_NULL;
  if (!!density) {
    sourcename="density";
    sourcep=&density;
  }
  else {
    if (!!sigma0) {
      sourcename="sigma0";
      sourcep=&sigma0;
    }
    else {
      if (needed)
	error_abort("Can't construct filter matrix as neither density matrix nor start_operator is defined");
      return false;
    }
  }
  const block_pattern& sourcepattern=sourcep->blockstructure();
  if ((filterp!=NMRSIM_NULL) && (filterp->blockstructure()==sourcepattern))
    return true;

  if (verbose & VER_GEN)
    std::cout << ((filterp==NMRSIM_NULL) ? "" : "Re-") << "Building spin-order/coherence filter " << spec << " using " << sourcename << " structure\n";
      
  //  blkspec=sourcepattern;
  
  //if (isspinorder()) {
    //    const spinorder_spec& lspec(spec(Type2Type<spinorder_spec>()));
  filterp=masterobjp->create_spinorder(spec,sourcepattern);
  if (!nochecks && !(filterp->isnonzero())) {
    if (sourcepattern.ishermitian())
      zerofilterhermitian_warning.raise();
    else
      zerofilter_warning.raise();
  }
    //}
    // else {
    //  throw Failed("Withdrawn functionality");
    //    const filter_spec& lspec(spec(Type2Type<filter_spec>())); //!< coherence list only
    // filterp=masterobjp->create_filter(lspec);
    // }
  return true;
}


bool matrixspec::isconstant() const
{
  switch (type) {
  case OPERATOR: {
    const operator_def& def(asoperator());
    return def.specp ? (def.specp)->isconstant() : false;
  }
  case BOOL:
    //   if (asfilter().isspinorder()) {
    if (!sigma0_specp)
      return false;
      //	throw Failed("can't determine const status of spin-order filter");
    return sigma0_specp->isconstant();
    //}
    return true;
  case EXCHANGE:
    return asexchange().isconstant();

  default:
    break;
  }
  return true;
}

void matrixspec::update()
{
  const bool isverb=((verbose & VER_GEN) && (verbose_level>1));

  switch (type) {
  case OPERATOR: {
    const operator_def& def(asoperator());
    if ( (&(def.op)==&sigma0) || (&(def.op)==&detect) || (def.specp==NMRSIM_NULL))
      return; //!< ignore built-ins
    if ( !(def.specp->isconstant()) || !(def.op)) {
      const productoperator_spec spec(def.specp->current());
      if (isverb)
	std::cout << "Making operator: " << spec << '\n';
      def.op=masterobjp->make(spec); //!< if empty or non-constant, remake
    }    
  }
    break;
  case EXCHANGE: {
    ExchangeMatrix& mat(const_cast<ExchangeMatrix& >(asexchange()));
    if (!mat || (!(mat.isconstant())))
      mat.update();
  }
    break;

    //!< build as required, since potentially depends on current density matrix structure 
//   case BOOL: {
//     filter_def& def( const_cast<filter_def& >(asfilter()));
//     def.ensure();
//   }
//     break;

  default:
    break; //!< do nothing otherwise
  }
}
  
void post_channels()
{
  if (blockingnuclei.empty())
    error_abort("channels called before spin system");

  for (size_t i=nucids.size();i--;) {
    const size_t nuc=nucids(i);
    nuclei_spec* iter=::std::find(blockingnuclei.begin(),blockingnuclei.end(),nucids(i));
    if (iter==blockingnuclei.end()) {
      parser_printcontext() << "RF channel " << nuctolabel(nuc) << " not present in spin system (or used twice)\n";
      error_abort(ERR_FAILED);
    }
    blockingnuclei.erase(iter);
  }
}

const basespin_system* get_spin_system()
{
  if (!sysp)
    error_abort("No spin system defined: missing spinsys block?");
  return sysp;
}

void make_pulseq_variables()
{
  add_systemvarmap(v_time);  
  if (spin_rate)
    add_systemvarmap(v_rotor_phase);
}

#endif

//double acqseq_synctime=0.0; //!< common synchronisation time for acqusition (0 if unset)

bool update_interactions=true;

size_t nchannels=0; //!< number of RF channels (0 before channels)
//double dphs=0.0; //!< receiver phase shift
double master_time=0.0; //!< current sequence time

setableoperator_spec* sigma0_specp=NMRSIM_NULL; //!< initial density matrix specification (\c NMRSIM_NULLif unset)
setableoperator_spec* detect_specp=NMRSIM_NULL; //!< detect operator specification (\c NMRSIM_NULLif unset)

//SystemVariable<double*> v_dphs("dphs",&dphs,rad_to_deg);
//SystemVariable<double*> v_synctime("acq_synctime",&acqseq_synctime,1e6);

LIST<size_t> nucids;

void add_channel(size_t nuc)
{
  nucids.push_back(nuc);
  nchannels=nucids.size();
}

void parse_channels()
{
  char* cptr;
  while ((cptr=parse_string(F_REPLACEDOLLAR | F_ALLOWMISSING)))
    add_channel(parse_nucleusname(cptr));

  post_channels();
}

ThreadWarning<> fullnotimp_warning("-full currently only supported for operator matrices.  Ignoring flag for: ",&NMRsim_once_warning);

// void MasterObj::expandoperator(cmatrix& dest, const BlockedOperator& op) const
// {
// #ifndef NOPERIODIC
//   if (!simple_opgenp)
//     op.full(dest,*crystal_opgenp);
//   else
// #endif
//     op.full(dest,*simple_opgenp);
// }

const Matrix<bool>& filter_def::ensurefull() const
{
  if (!filterp)
    filterfull.clear();
  else {
    if (filterfull.empty() || (blkspec!=filterp->blockstructure())) {
      filterp->full(filterfull);
      blkspec=filterp->blockstructure();
    }
  }
  return filterfull;
}
			       
ThreadWarning<> undefinedop_warning("Operator undefined at this point: ",&NMRsim_repeat_warning);

void dump_matrix(const matrixspec& curm, const char* name, MasterObj* Objp, int flags, statistics_t statssel)
{
  const bool full= flags & PM_FULL;
  const int sflags= flags & PM_STRUCTUREFLAGS;

  logfile_controller* logfilep(get_logfile());
  if (curm.type==matrixspec::SPECIAL) {
    if (!Objp)
      error_abort("Hamiltonian is not defined at this point");
  }
  else {
    if (curm.type!=matrixspec::OPERATOR) {
      if (full)
	fullnotimp_warning.raise(name);
      if (sflags)
	error_abort("-structure / -eigenbasis only implemented for Hamiltonians and operators");
    }
    if (curm.type==matrixspec::BOOL) {
      filter_def& asmut = const_cast< filter_def& >(curm.asfilter());
      asmut.ensure(false); //!< try to create
    }
  }

  static double tol=0.0;
  if (tol==0.0) {
    const char* tolstr=register_and_getenv("NMRSIM_COHERENCE_TOLERANCE");
    if (tolstr) {
      char* endptr;
      tol=strtod(tolstr,&endptr);
      if (*endptr!='\0') {
	parser_printcontext() << "Couldn't convert '" << tolstr << "' to valid floating point number\n";
	error_abort();
      }
      if (tol<=0.0) {
	parser_printcontext() << "Negative or zero tolerance specified in NMRSIM_COHERENCE_TOLERANCE\n";
	error_abort();
      }
    }
    else
      tol=NMRSIM_DEFAULT_COHERENCE_TOLERANCE;
  }

  if (!logfilep) {
    std::cout << name;
    switch (sflags) {
    case PM_STRUCTURE:
      std::cout << " structure:\n";
      break;
    case PM_EIGENBASIS:
      std::cout << " eigenbasis:\n";
      break;
    case PM_STATISTICS:
      std::cout << " statistics:\n";
      break;
    default:
      std::cout << '\n';
    }
    if (!full && (curm.type==matrixspec::OPERATOR) && (sflags!=PM_STATISTICS) && (curm.asoperator().op.blockstructure().ishermitian()))
      std::cout << "[Information for blocks above diagonal only; use -full for full matrix]\n";

    switch (curm.type) {
    case matrixspec::SPECIAL:
      Objp->printH(std::cout,sflags);
      break;
    case matrixspec::OPERATOR: {
      const BlockedOperator& op(curm.asoperator().op);
      if (sflags==PM_STATISTICS) {
	operator_stats stats(op,statssel,tol);
	stats.print(std::cout);
      }
      else {
	if (full) {
	  const cmatrix tmp(op.full());
	  std::cout << tmp;	
	}
	else {
	  //std::cout << op.blockstructure() << '\n';
	  spyprint(std::cout,op,sflags);
	  if ((verbose & VER_PROFILE) && (sflags==0))
	    std::cout << curm.usage();
	}
      }
    }
      break;
    default:
      std::cout << curm;
      if (verbose & VER_PROFILE)
	std::cout << curm.usage();
    }
  }
  else {
    char tmptitle[128];
    char title[128];
    switch (sflags) {
    case PM_STRUCTURE:      
      snprintf(tmptitle,sizeof(tmptitle),"%s_structure",name);
      break;
    case PM_EIGENBASIS:
      snprintf(tmptitle,sizeof(tmptitle),"%s_eigenbasis",name);
      break;
    case PM_STATISTICS:
      snprintf(tmptitle,sizeof(tmptitle),"%s_statistics",name);
      break;      
    default:
      strncpy(tmptitle,name,sizeof(tmptitle));
    }
    logfile_controller::maketitle(title,sizeof(title),tmptitle);

    if (curm.type==matrixspec::SPECIAL)
      Objp->logH(*logfilep,title,sflags);
    else {
      if (sflags==PM_STATISTICS) {
	const BlockedOperator& op(curm.asoperator().op);
	if (!op)
	  undefinedop_warning.raise(name);
	else {
	  operator_stats stats(curm.asoperator().op,statssel,tol);
	  stats.log(*logfilep,title);
	}
      }
      else	
	curm.log(*logfilep,title,flags);
    }
  }
}

const flagsmap_type& putmatrix_flags()
{
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["full"]=PM_FULL;
    flags["structure"]=PM_STRUCTURE;
    flags["once"]=PM_ONCE;
    flags["eigenbasis"]=PM_EIGENBASIS;
    flags["statistics"]=PM_STATISTICS;
  }
  return flags;
}

std::pair<int,statistics_t> get_putmatrix_flags()
{
  const int flags=parse_flags(putmatrix_flags());
  const int sflags=flags & PM_STRUCTUREFLAGS;
  statistics_t statssel(0, (const filter_def*)NMRSIM_NULL );
  if (sflags) {
    if ( (flags & PM_FULL) && sflags)
      error_abort("-full and -structure/-eigenbasis flags cannot be sensibly combined");
    switch (sflags) {
    case PM_STRUCTURE: case PM_EIGENBASIS:
      break;
    case PM_STATISTICS:
      if (are_left()) {
	List<size_t> spinsel(parse_unsignedintarray());
	adjust_and_verify_spin_indices(spinsel);
	const size_t nspins=sysp->nspins(); //!< don't need to re-verify sysp
	statssel.first=spinorder_spec::compress_spinindices(nspins,spinsel);
	if (are_left())
	  statssel.second=&(getfiltermatrix(parse_string()));
      }
      break;
    default:
      error_abort("can only use one of -structure, -eigenbasis, -statistics");
    }
  }
  return std::pair<int,statistics_t>(flags,statssel);
}

ContextWarning<> nonflags_warning("-once flag is meaningless in par block. More appropriate for use in pulseq block once operators / matrices have been defined",&NMRsim_repeat_warning);

void parse_par_putmatrix()
{
  const char* name(parse_string(F_REPLACEDOLLAR));
  const matrixspec curm(findmap(matrixmap,name,"matrix"));
  const std::pair<int,statistics_t> flaginfo=get_putmatrix_flags();
  const int flags=flaginfo.first;
  if (flags & PM_ONCE)
    nonflags_warning.raise();
  if (flags & PM_STATISTICS) {
    if (curm.type!=matrixspec::OPERATOR)
      error_abort("statistics output only defined for operator matrices");
  }
  dump_matrix(curm,name, NMRSIM_NULL, flags, flaginfo.second);
}

void make_1d_par_variables()
{
  make_data_variables();
  command_Factory_t& par_Factory(get_par_Factory());
  par_Factory["putmatrix"]=par_t(&parse_par_putmatrix,true);
  add_systemvarmap(v_proton_freq);
  add_systemvarmap(v_cputime);
}

context_t evaluation_state=CONTEXT_PRE; //!< pre evaluation

void initialise_simulation_environment()
{
  //adjust diagonalisation environment
  cmatrix_eigensystem_controller.tolerance = nochecks ? 0.0 : eigtolerance;
  cmatrix_eigensystem_controller.effectivezero = zerotolerance;
  int everbose = ((verbose & VER_GEN) && (verbose_level>1)) ? verbose_level : 0;
  if ((everbose==0) && abortonwarning) //!< ensures that matrices will be dumped
    everbose=1;
  cmatrix_eigensystem_controller.verbose = everbose;
  cmatrix_eigensystem_controller.throwexception = abortonwarning;
  evaluation_state=CONTEXT_PRE;
}

void MasterObj::initialise_simple_simulation()
{
  restart();
  ::initialise_simulation_environment();
  recreate_Hs(); //!< set up Hamiltonians
}

void parse_detect_operator()
{  
  detect_specp=parse_setableoperator();
}

void parse_start_operator()
{
  sigma0_specp=parse_setableoperator();
}

void multioperator_spec::swap(multioperator_spec& a)
{
	LIST<productoperator_spec>::swap(a);
	std::swap(arraytag_, a.arraytag_);
	std::swap(nuc_, a.nuc_);
}

void multioperator_spec::setnucleus()
{
	if (empty())
		throw InternalError("setnucleus() called before any operators defined");
	nuc_ = front().nucleus();
	for (size_t j=1; j<size(); j++) {
		if ((*this)(j).nucleus() != nuc_) {
			nuc_ = NULL_NUCLEUS;
			return;
		}
	}
}
			
size_t setableoperator_spec::currentindex() const
{
	if (size()==1)
		return 0;
	size_t transindex;
	(void)getindex(transindex, arraytag_, size()); //!< note ignores update status (potential missed optimisation with virtual dimension)
	return transindex;
}

ListList<int> parse_coherencelist()
{
  if (count_left()!=nchannels) {
    parser_printcontext() << "number of remaining arguments (" << count_left() << ") doesn't match number of RF channels (" << nchannels << ") - check directive syntax\n";
    error_abort();
  }
  ListList<int> cohers;
  for (size_t n=0;n<nchannels;n++)
    cohers.push_back(parse_intarray());
  if (verbose & VER_PARSE)
    std::cout << "Coherence list: " << cohers << '\n';
  return cohers;
}

// BlockedFilter* MasterObj::create_filter(const filter_spec& cohers) const
// { 
// #ifndef NOPERIODIC
//   if (!simple_opgenp)
//     return new BlockedFilter(*crystal_opgenp,nucids,cohers);
//   else
// #endif
//     return new BlockedFilter(*simple_opgenp,nucids,cohers);
// }

BlockedFilter* MasterObj::create_spinorder(const spinorder_spec& spec, const block_pattern& blkp) const
{
#ifndef NOPERIODIC
  if (!simple_opgenp)
    error_abort("filter not implemented with periodic systems. Use -disable:periodic");
#endif
  try {
    return new BlockedFilter(*simple_opgenp,blkp,spec);
  }
  catch (MatrixException& exc) {
    std::cerr << "Couldn't create filter matrix - " << exc;
  }
  error_abort("Re-try with -disable:mzblocking flag?");
  return NMRSIM_NULL; //!< suppresses warning
}

// void parse_coherencematrix(const char* matname)
// {
//   if (nchannels==0)
//     error_abort("coherence filters can only be applied to RF channels (no RF channels active!)");
  
//   filter_def* fdefp=new filter_def(parse_coherencelist());
//   //  BlockedFilter* newf=masterobjp->create_filter(cohers);

//   //  if (verbose & VER_GEN)
//   //  std::cout << "Filter matrix " << matname << '\n' << (*newf) << '\n';
//   matrixmap[matname]= matrixspec(*fdefp);
// }

ContextWarning<> overlappingranges_warning("overlapping sets of spin indices in spin order specification - this may be an error (disable warning with -nochecks)",&NMRsim_once_warning);

ContextWarning<> coherunset_warning("filter matrix defined without specifying coherences of all spins - may give unexpected results",&NMRsim_once_warning);

template<bool iscoher =false> struct order_traits {
  typedef size_t base_t;
  static inline LIST<size_t> parse() { return parse_unsignedintarray(); }
};
template<> struct order_traits<true> {
  typedef int base_t;
  static inline LIST<int> parse() { return parse_intarray(); }
};

template<bool iscoher> filter_def* parse_ordermatrix_(bool onlytotal)
{
  typedef typename order_traits<iscoher>::base_t base_t;

  const size_t n=getnspins();

  try {

    int nleft=onlytotal ? 1 : count_left();

    if (--nleft==0) { // forced if onlytotal set
      LIST<base_t> arg1( order_traits<iscoher>::parse() );
      return new filter_def(spinorder_spec(n,arg1));
    }
    
    if ((nleft & 1)==0)
      error_abort("spin- / coherence- order filter must be specified as either <order list> (all spins) or as repeating pairs of <spin list> <order list>");

    LIST<size_t> arg1;    
    spinorder_spec spec(n,iscoher);
    while (nleft>0) {
      arg1=parse_unsignedintarray();
      adjust_and_verify_spin_indices(arg1);
      spec.add(arg1, order_traits<iscoher>::parse() );
      nleft-=2;
    }
    if (!nochecks) {
      if (spec.overlappingranges())
	overlappingranges_warning.raise();

      if (iscoher && !(spec.allspinsused()))
	coherunset_warning.raise();
    }
    return new filter_def(spec);
  }
  catch (MatrixException& exc) {
    parser_printcontext() << "Couldn't create spin- / coherence- order specification: " << exc.what() << '\n';
    error_abort();
  }
  return NMRSIM_NULL; //!< shouldn't get this far
}

void parse_ordermatrix(const char* matname, bool iscoher, bool onlytotal)
{  
  filter_def* fdefp= iscoher ? parse_ordermatrix_<true>(onlytotal) : parse_ordermatrix_<false>(onlytotal);
  if (verbose & VER_PARSE)
    parser_printcontext() << "parsed spin- / coherence- order specification to give matrix " << matname << ": " << fdefp->spec << '\n';
  matrixmap[matname]= matrixspec(*fdefp);
}

void parse_makefilter()
{
  error_abort("filter deprecated - use matrix set <name> coherence etc.");
  //  parse_coherencematrix(parse_string(F_REPLACEDOLLAR));
}
  
std::ostream& operator<< (std::ostream& ostr, const usage_t& usage)
{
  ostr << (usage.bytes/1024.0);
  if (usage.bytes)
    ostr << " K (" << usage.items << " items)";
  return ostr;
}

template<class T> void makeU(BlockedMatrix<complex>& U, double t1, double t2, const T& H, double maxdt =0.0)
{
  if (!H) {
    U.clear();
    return;
  }
  const int verb = (verbose & VER_GEN) ? verbose_level : 0;
  MetaPropagator propgen(H,maxdt,verb,global_prop_flags);
  set_partitioning(propgen);
  propgen(U,t1,t2);
}

void MasterObj::propagator(BlockedMatrix<complex>& U, double t1, double t2)
{
  switch (hamobjtype) {
  case M_SPIN:
    if (Hspinp.iscomplex())
      makeU(U,t1,t2,*(Hspinp.get_complex()),get_maxdt());
    else
      makeU(U,t1,t2,*(Hspinp.get_real()),get_maxdt());
    break;
  case M_STATICD:
    makeU(U,t1,t2,*Hstaticdp);
    break;
  case M_SPIND:
    makeU(U,t1,t2,*Hspindp);
    break;
  case M_STATIC:
    if (Hstaticp.iscomplex())
      makeU(U,t1,t2,*(Hstaticp.get_complex()));
    else
      makeU(U,t1,t2,*(Hstaticp.get_real()));
    break;
  default:
    throw Failed("Unknown hamobjtype");
  }
//   BlockedMatrix<complex> Utrans;
//   flush_transients(Utrans); //need to incorporate any unflushed transients
//   if (!!Utrans)
//     U*=Utrans;//U = U * Utrans (i.e. Utrans applied first)
}

bool have_spinsys() { return (sysp!=NMRSIM_NULL); }

ExchangeMatrix::ExchangeMatrix(const char* namev, const BaseList<double>& rawvalsv, bool isconstv, bool isgeneralv) :
  name_(namev),
  sites_(0U),
  rawvals_(rawvalsv),
  isconst_(isconstv),
  isgeneral_(isgeneralv)
{
  update();
}

// ExchangeMatrix::ExchangeMatrix(const char* namev, const BaseList<double>& rawvalsv, const BaseList<double>& rawvalsiv, bool isconstv) :
//   name_(namev),
//   sites_(0U),
//   rawvalsr_(rawvalsv),
//   rawvalsi_(rawvalsiv),
//   isconst_(isconstv),
//   isgeneral_(true)
// {
//   update();
// }

ThreadWarning<> ExchangeMatrix::unphysical_warning("exchange matrix elements are expected to lie between 0 and 1: ",&NMRsim_once_warning);

void ExchangeMatrix::checkset(size_t r, size_t c, double v)
{
  char buf[256];
  if (v<0.0 || v>1.0) {
    snprintf(buf,sizeof(buf),"found %g when filling element %" LCM_PRI_SIZE_T_MODIFIER "u,%" LCM_PRI_SIZE_T_MODIFIER "u",v,r+1,c+1);
    unphysical_warning.raise(buf);
  }
  matrix_(r,c)=v;
}

//!< quick table-based integer sqrt
size_t intsqrt(size_t n)
{
  const size_t defmin=25U; //!< by default find roots up to 25
  static LIST<size_t> cache;
  if (cache.size()<=n) {
    size_t todo=(n<defmin) ? defmin : n;
    cache.create(todo+1,size_t(0));
    for (size_t i=0;;i++) {
      const size_t s=i*i;
      if (s>todo)
	break;
      cache(s)=i;
    }
  }
  return cache(n);
}

void ExchangeMatrix::update()
{
  size_t nels=rawvals_.size();
//   if (nels==0) { //!< special case of purely imaginary
//     nels=rawvalsi_.size();
//     if (nels==0) {
//       matrix_.clear();
//       return;
//     }
//     rawvalsr_.create(nels,0.0);
//   }
//   const bool isr=isreal();
//   if (!isr && (rawvalsi_.size()!=nels))
//     throw Mismatch("ExchangeMatrix::update");
 
  size_t n=0U;
  bool isr=true;
  if (isgeneral_) {
    n=intsqrt(nels);
    if (n==0) {
      if (((nels & 1)==0)) {
	n=intsqrt(nels/2);
	isr=false;
      }
      if (n==0) {
	parser_printcontext() << "number of supplied elements (" << nels << ") can't be number of elements in a square real or complex matrix\n";
	error_abort();
      }
    }
  }
  else {
    const size_t offdiagels=(nels & 1) ? 2*nels : nels;
    n=(1+intsqrt(4*offdiagels+1))/2;
    const size_t nhalf=n*(n-1)/2;
    if ((nels!=nhalf) && (nels!=2*nhalf)) {
      parser_printcontext() << "number of supplied matrix elements (" << nels << ") doesn't make sense as either the upper triangle or all off-diagonal elements of an exchange matrix\n";
      error_abort();
    }
//     if (!isr)
//       throw InternalError("exchange-only matrix must be purely real");  
  }
    
  matrix_.create(n,n);
  isreal_=isr;
  if (isgeneral_) {
    BaseList<complex> dest(matrix_.row());
    if (isr)
      dest=rawvals_;
    else {
      size_t i,j;
      for (i=j=0;j<nels;i++,j+=2)
	dest(i)=complex(rawvals_(j),rawvals_(j+1));
    }
    return;
  }
  const bool full=((nels & 1)==0);
  size_t r=0;
  size_t c=1;
  for (size_t i=0;i<nels;i++) {
    const double v=rawvals_(i);
    checkset(r,c,v);
    matrix_(r,c)=v;
    if (!full)
      matrix_(c,r)=v;
    if (++c==n) {
      r++;
      c= full ? 0U : r+1;
    }
    else {
      if (c==r)
	c++; //skip diagonal
    }
  }
  for (size_t c=n;c--;) {
    double diagv=1.0;
    for (r=n;r--;) {
      if (r!=c)
	diagv-=real(matrix_(r,c));
    }
    checkset(c,c,diagv);
  }
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Updated exchange matrix " << name_ << " from elements: " << rawvals_ << "\nto give matrix\n" << matrix_;
} 

void ExchangeMatrix::set(const BaseList<double>& vals, subsid_t)
{
  rawvals_=vals;
  update();
}

void ExchangeMatrix::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << name_;// << ((subsid==S_ARG1) ? "_real" : "_imag");
}
  
std::ostream& operator<< (std::ostream& ostr, const ExchangeMatrix& a)
{
  if (!a)
    return ostr << "[undefined in this context]\n";
  if (a.isreal())
    return ostr << real(a());
  else
    return ostr << a();
}

void ExchangeMatrix::printraw(std::ostream& ostr) const
{
  if (!(*this)) {
    ostr << "[undefined in this context]\n";
    return;
  }
  ostr << rawvals_;
  //if (!isreal())
  //  ostr << ' ' << rawvalsi_;
}
