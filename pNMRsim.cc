#include "Action.h"
#include "NMRsim_Process.h"
#include "Parser.h"
#include "NMRsim_spinsys.h"
#include "cmatrix_external.h"
#include "NMR.h"
#include "MAS.h"
#include "ttyio.h"
#include "Propagation.h"
#include "InversionGenerator.h"

#ifdef HAVE_SYS_RESOURCE
#include <sys/resource.h>
#endif

using namespace libcmatrix;
using namespace std;

template<typename First, typename Second> std::pair<First,Second> operator+ (const std::pair<First,Second>& a, const std::pair<First,Second>& b)
{
  return std::pair<First,Second>(a.first+b.first,a.second+b.second);
}

#ifdef NDEBUG
#define VERBOSE_LEVEL 2
#else
#define VERBOSE_LEVEL 3
#endif

int F_defaultdataset=F_DENYSUM;

bool debug=false;

bool default_calculation(MasterObj& masterobj, DataStore& Spec)
{
  masterobj.calc(Spec,verbose,false);
  return true;
}

precalculation_callback_t precalculation_callback=NMRSIM_NULL;
calculation_callback_t calculation_callback=NMRSIM_NULL;
prefinalise_callback_t prefinalise_callback=NMRSIM_NULL;

bool need_spinsys=true;
option optperiodic("periodic");
option opteigsym("eigsymmetry","eigenvalue symmetry",option::AUTO,option::NOTUSED);
option optblocking("mzblocking","",option::AUTO,option::NOTUSED);
option optmzsym("mzsymmetry","mz eigenvalue symmetry");
option optmergeproc("mergeprocessing","merge processing blocks",option::AUTO,option::NOTUSED);
bool randomise=false;

LIST<procstack_t> procstacks,postprocstacks,finalisestacks;
InternalCommand* startup_cp=NMRSIM_NULL;

// template<typename T> class Deref {
// public:
//   template<typename T2> Deref(const std::map<T2,T*>&, const std::pair<T2,T*>& v)
//     : value_(v.second) {}
//   Deref(const std::list<T*>&, const T& v)
//     : value_(v) {}
//   Deref(const LIST<T*>&, const T& v)
//     : value_(v) {}
//   const T& operator()() const { return value_; }

// private:
//   const T& value_;
// };

template<class T> void ldumpptr(const BaseList<T>& a, const char* name)
{
  const size_t nstacks=a.size();
  if (nstacks==0)
    return;
  for (size_t i=0;i<nstacks;i++) {
    std::cout << name;
    if (nstacks>1)
      std::cout << ' ' << (i+1);
    std::cout << ":\n";
    dumpptr(a(i),"\n");
  }
  std::cout << '\n';
}

typedef std::pair<double,usage_t> tusage_t;

template<class T> std::ostream& dumpname(const T& a, std::ostream& ostr)
{
  return ostr << a.name();
}

template<> std::ostream& dumpname(const ActionCommand& a, std::ostream& ostr)
{
  a.print(ostr);
  return ostr;
}
template<> std::ostream& dumpname(const InternalCommand& a, std::ostream& ostr)
{
  a.print(ostr);
  return ostr;
}
template<> std::ostream& dumpname(const ProcessCommand& a, std::ostream& ostr)
{
  a.print(ostr);
  return ostr;
}

void dumpfeatures(std::ostream& ostr, const LIST<const option*>& flist, const char* prefix)
{
  if (flist.empty())
    return;
  ostr << prefix << ':';
  for (size_t i=0;i<flist.size();i++)
    ostr << ' ' << flist(i)->name();
  ostr << '\n';
}

void dump_optimisationusage(std::ostream& ostr =std::cout)
{
  const optional_map_t& optional_map(get_optional_map());  
  const optional_map_t::const_iterator end(optional_map.end());
  optional_map_t::const_iterator start(optional_map.begin());
  if (start==end)
    return;
  LIST<const option*> usedlist,nousedlist,noinfolist;
  while (start!=end) {
    const option* curop(start->second);
    ++start;
    if (curop->isdisabled()) {
      if (curop->getusage()==option::USED)
	throw InternalError("Disabled option apparently used"); //!< should have been caught earlier, but why not...      
      nousedlist.push_back(curop); //!< if disabled, automatically add to not used list
    }      
    else {
      switch (curop->getusage()) {
      case option::USED:
	usedlist.push_back(curop);
	break;
      case option::NOTUSED:
	nousedlist.push_back(curop);
	break;
      case option::NOTKNOWN:
	noinfolist.push_back(curop);
	break;
      }
    }
  }  
  ostr << '\n';
  dumpfeatures(ostr,usedlist,"Used features / optimisations");
  dumpfeatures(ostr,nousedlist,"Unused features / optimisations");
  dumpfeatures(ostr,noinfolist,"Unknown status");
}

template<class T> usage_t dump_timer_(const T& v, bool abbrev, std::ostream& ostr, double dur =0.0)
{
  //  const double dur=v.time();
  const usage_t use(v.usage());
  if (!abbrev) {
    dumpname(v,ostr) << ": ";
    if (dur)
      prettyprint_time(dur,ostr);
    else
      ostr << "<negligible>";
    
    if (use)
      ostr << "  Memory: " << use;
    ostr << '\n';
  }
  return use;
}
      

template<class T> tusage_t dump_timer(const MAPTYPE(T*)& a, const char* name, bool abbrev =false, std::ostream& ostr =std::cout)
{
  usage_t tot;
  const typename MAPTYPE(T*)::const_iterator end(a.end());
  typename MAPTYPE(T*)::const_iterator start(a.begin());
  if (start!=end) {    
    if (!abbrev)
      ostr << name << '\n';
    while (start!=end) {
      tot=tot+dump_timer_(*(start->second),abbrev,ostr);
      ++start;
    }
    if (abbrev) {
      ostr << name << ": ";
      //      prettyprint_time(tot.first,ostr);
      ostr << " Memory: " << tot << '\n';
    }
  }
  return tusage_t(0.0,tot);
}

template<class T> tusage_t dump_timer(const T& a, const char* name, bool abbrev =false, std::ostream& ostr =std::cout)
{
  double tott=0.0;
  usage_t tot;
  const typename T::const_iterator end(a.end());
  typename T::const_iterator start(a.begin());
  if (start!=end) {    
    if (!abbrev)
      ostr << name << '\n';
    while (start!=end) {
      const double dur=(*start)->time();
      tott+=dur;
      tot=tot+dump_timer_(**start,abbrev,ostr,dur);
      ++start;
    }
    if (abbrev && (tott || tot.bytes) ) {
      ostr << name << ": ";
      prettyprint_time(tott,ostr);
      ostr << "  memory: " << tot << '\n';
    }
  }
  return tusage_t(tott,tot);
}

template<class T> tusage_t dump_timers(const BaseList<T>& a, const char* name,  bool abbrev =false, std::ostream& ostr =std::cout)
{
  tusage_t tot;
  const size_t nstacks=a.size();
  if (nstacks==0)
    return tusage_t(0.0,usage_t());
  char buf[256];
  for (size_t i=0;i<nstacks;i++) {	  
    const T& curl(a(i));
    if (curl.empty())
      continue;
    const char* usename=name;
    if (nstacks>1) {
      snprintf(buf,sizeof(buf)-1,"%s %" LCM_PRI_SIZE_T_MODIFIER "u",name,i+1);
      usename=buf;
    }
    tot=tot+dump_timer(curl,usename,abbrev,ostr);
  }
  if (!abbrev)
    ostr << '\n';
  return tot;
}

void dump_timings(bool abbrev)
{
  if (!abbrev)
    std::cout << "\nInternal functions\n";
  
  usage_t tott(masterobjp->usage());
  {
    const seqmap_type::const_iterator end(seqmap.end());
    seqmap_type::const_iterator start(seqmap.begin());      
    while (start!=end) {
      tott+=start->second->usage();
      ++start;
    }
  }
  {
    const matrixmap_type::const_iterator end(matrixmap.end());
    matrixmap_type::const_iterator start(matrixmap.begin());      
    while (start!=end) {
      tott+=start->second.usage();
      ++start;
    }
  }
  startup_cp->usage(tott);
  tusage_t total=dump_timer(internal_timers,"Internal timers",abbrev);
  const tusage_t mainusage=dump_timers(actionstacks,"Pulse sequence",abbrev);
  total=total+mainusage;
  total=total+dump_timer(seqmap,"Sequence fragments",abbrev);
  total=total+dump_timer(cycledseqmap,"Cycled sequences",abbrev);
  total=total+dump_timers(procstacks,"Per-row processing",abbrev);
  total=total+dump_timers(postprocstacks,"Per-calculation processing",abbrev);
  total=total+dump_timers(finalisestacks,"Finalise processing",abbrev);
  
  if (timers_arrayp && !abbrev) {
    std::cout << '\n';
    dump_arrayprofile(*timers_arrayp,1000.0,"ms");
  }

  const usage_t& mtotal(total.second);
  static const double bytes_to_M=1.0/(1024*1024);
  if (!abbrev)
    std::cout << '\n';
  std::cout << "Global propagator cache: " << (global_cache()*bytes_to_M) << " M\n";
  std::cout << "Total accounted-for memory: " << ((mtotal.bytes+global_cache())*bytes_to_M) << " M\n\n";

  std::cout << " ---- Low level allocation ---\n";
  static const double words_to_M=4.0/(1024*1024);
  if (verbose_level>1)
    matrix_traits<complex>::allocator::print(std::cout);
  std::cout << "Memory allocated to complex matrices: " << (cmatrix::total_allocated_words()*words_to_M) << " M\n";    
  if (verbose_level>1) {
    std::cout << '\n';
    matrix_traits<double>::allocator::print(std::cout);
  }
  std::cout << "Memory allocated to real matrices: " << (rmatrix::total_allocated_words()*words_to_M) << " M\n";

#ifdef HAVE_SYS_RESOURCE
  struct rusage useobj;
  if (getrusage(RUSAGE_SELF, &useobj)!=0)
    std::cerr << "rusage failed - can't report system statistics\n";
  else
    std::cout << "\nMaximum memory allocated by system to process: " << (useobj.ru_maxrss*(1.0/1024)) << " M\n";
#endif  
}

const char* NEWSstr=
#include "NEWS"

void printnews()
{
  std::cout << "Recent changes:\n";
  std::cout << NEWSstr;
  //#ifdef __DATE__
  //   std::cout << "\n(I was compiled on " << __DATE__ << ")\n";
  //#endif
}

#define XLocalStr(x) #x
#define LocalStr(x) XLocalStr(x)

void printversion()
{
  unsigned short day,month,year;
  if (sscanf(NEWSstr,"%hu/%hu/%hu",&day,&month,&year)!=3)
    throw Failed("Failed to parse data from NEWS");
  printf("%02u.%02u.%02u\n",year,month,day);
  printf("Minuit support: %s",
#ifdef USE_MINUIT
	 "Yes"
#else
	 "No"
#endif
	 );
  printf("\t  Periodic basis sets: %s\n",
#ifdef NOPERIODIC 
	 "No" 
#else
	 "Yes"
#endif
	 );
  printf("Optimised matrix ops: %s",
#ifdef LCM_USE_EXTERNAL
	 LocalStr(LCM_EXTERNAL_NAME) 
#else
	 "No"
#endif
	 );
  printf("\t  Parallel system: %s\n",
#ifdef USEMPI
	 "MPI" 
#else
#ifdef HAVE_FORK_CONTROLLER
	 "Threading"
#else
	 "None"
#endif
#endif
	 );
  printf("Optimised complex type: %s",
#ifdef LCM_USE_SSECOMPLEX
	 "Yes" 
#else
	 "No"
#endif
	 );
  printf("\t  Natural offset sign: %s\n",
#ifdef NMRSIM_INTUITIVE_OFFSET
	 "Yes"
#else
	 "No"
#endif
	 );
  puts("Nuclear spin properties: " NMRSIM_NUCLEUS_PROPERTIES);
}

int global_argc=0;
char** global_argv;

//Warning<> offdiagonal_warning("using nondiagonal product operator - any mz blocking is being suppressed",&NMRsim_once_warning);
ThreadWarning<> parprofile_warning("verbose profile of limited use with parallel execution; stats returned from control thread only",&NMRsim_repeat_warning);
ThreadWarning<> buildfailed_warning("failed to build spin operators using expected blocking and so disabling mz blocking. You may be able to use the channels directive to selectively disable blocking e.g. for non-secular quadrupoles.",&NMRsim_repeat_warning);
ThreadWarning<> mzsymmetry_warning("mzsymmetry optimisation explicitly enabled in case where it is unlikely to be appropriate",&NMRsim_repeat_warning);
ThreadWarning<> missingacq_warning("no direct acquisition (missing acq in pulseq?) and so processing blocks will be skipped",&NMRsim_once_warning);
ThreadWarning<> nothingtodo_warning("nothing to do - terminating!",&NMRsim_repeat_warning);
ThreadWarning<> nointeractions_butMAS_warning("no interactions specified but spin_rate is non-zero - static simulation would be much more efficient and MAS simulation may fail",&NMRsim_repeat_warning);
ThreadWarning<> defaultsigma0_warning("no start_operator specified, defaulting to ",&NMRsim_repeat_warning);

int main(int argc_,char **argv_)
{
  ensure_library_version("3.13.0");
  set_nucleus_properties(NMRSIM_NUCLEUS_PROPERTIES);
  
  startup_cp=create_timer("Initialise/global");
  startup_cp->enter();

  bool havedisable=false;

  //  optional_map_t& optional_map=get_optional_map();
  //  optional_map["parallel"]=&optparallel;
  add_option(optparallel);
  setdefaultoptions();

  //!< bit of a kludge - need to scan for parallel and verbose:powder flags before calling init_parallel which may in principle strip out arguments
  for (size_t i=1;(i<argc_) && (argv_[i][0]=='-');i++) {
    if (!checkoptflag(argv_[i],havedisable,true))
      checkverboseflag(argv_[i]);
  }

  init_parallel(argc_,argv_);

  global_argc=argc_;
  global_argv=argv_;

  try {

  int count=1;

  add_option(optperiodic);
  add_option(optcache);
  add_option(optmzsym);
  add_option(optblocking);
  add_option(opteigsym);
  add_option(optpartition);
  add_option(optmergeproc);
  setdefaultoptions();

  bool noshortcuts=false;  
  const char* qualstr=NMRSIM_NULL;

  while (count<argc_) {
    const char* com=argv_[count];
    if ((com[0]!='-') || (com[1]=='\0'))
      break;
    count++;
    
    if (strcmp(com,"-news")==0) {
      printnews();
      return 0;
    }
    if (strcmp(com,"-version")==0) {
      printversion();
      return 0;
    }
    if (!checkflags(nochecks,com,"nochecks") && !checkflags(noshortcuts,com,"noshortcuts") && !checkflags(debug,com,"debug") && !checkflags(noexecute,com,"noexecute") && !checkflags(abortonwarning,com,"abort") && !checkflags(silent,com,"silent") && !checkflags(randomise,com,"randomise")) {
      if (!checkoptflag(com,havedisable) && !checkverboseflag(com)) {
	if (strcmp(com,"-qualify")==0) {
	  if (count==argc_) {
	    cerr << "Missing argument: -qualify <string>\n";
	    return ERR_INVALID_INPUT;
	  }
	  qualstr=argv_[count++];
	}
	else {
	  cerr << "Unknown flag: " << com << '\n';
	  return ERR_INVALID_INPUT;
	}
      }
    }
  }
  if (randomise)
    set_seed(); //!< reset random number generator

  if (noshortcuts) { //!< turn off optimisations than have not been specifically enabled
    if (!nochecks && havedisable)
      std::cerr << "Warning: -disable flags redundant in presence of -noshortcuts\n";
    optional_map_t& optional_map(get_optional_map());  
    const optional_map_t::iterator end(optional_map.end());
    optional_map_t::iterator start(optional_map.begin());
    while (start!=end) {
      (start->second)->setdefault(option::OFF);
      ++start;
    }
  }
  ProcessSave::defaultstopoverwrite=optsaveprotect.isenabled(); //!< by default allow overwrite

  if (count==argc_) {
    cerr << "Syntax: pNMRsim [-news|-version]|[-noshortcuts|-nochecks|-noexecute|-debug|-enable:<option>|-disable:<option>|-qualify <string>|-silent|-randomise|-abort] <.in file>|-\n";
    return 1;
  }
  if (nochecks && abortonwarning)
    std::cerr << "Warning: -abort unlikely to have an effect with -nochecks enable\n";  
  verbose_level = debug ? VERBOSE_LEVEL : 1;
  if (debug && noexecute)
    verbose=VER_PARSE; //!< if debugging file read, turn on VER_PARSE
  if (!optcache)
    global_cache.limit(0);

  if (debug)
    setwarnings(BaseWarning::Always); //!< if debugging leave warnings turned on
  else {
    inconsistent_const_warning.type(BaseWarning::Ignore);
    inconsistent_nouses_warning.type(BaseWarning::Ignore);
    if (nochecks)
      setwarnings(BaseWarning::Ignore);
    else {
      if (silent)
	setwarnings(BaseWarning::Silent);
      else {
	NMR_asymmetry_warning.type(BaseWarning::FirstOnly);	
	propagation_closetodiagonal_warning.type(BaseWarning::FirstOnly);
	PhaseModulation::base_warning.type(BaseWarning::FirstOnly);
	lcm_sequence_warning.type(BaseWarning::FirstOnly);
      }
    }      
  }
  
  char* fname=argv_[count];
  char fnamebase[256]="";
  if (strcmp(fname,"-")!=0) {    
    strncat(fnamebase,fname,sizeof(fnamebase)-1); //!< strncat ensures proper termination (although filename is truncated)
    stripleaf(fnamebase,".in");
    if (qualstr)
      strncat(fnamebase,qualstr,sizeof(fnamebase)-strlen(fnamebase)-1);
    systemvarmap["name"]=new SystemVariable<std::string>("name",fnamebase);
  }
  systemvarmap["qualifier"]= new SystemVariable<std::string>("qualifier",qualstr ? qualstr : "");

  size_t argc=argc_-count-1;
  char** argv=argv_+count+1;

  if (argc && (argv[0][0]=='-'))
    std::cerr << "Warning: $1 begins with - (" << argv[0] << ").  Misplaced flag?\n";

  declare_builtin_block("spinsys");
  declare_builtin_block("par");
  //  declare_builtin_block("initialise");
  declare_builtin_block("pulseq");
  declare_builtin_block("initialproc");
  declare_builtin_block("proc");
  declare_builtin_block("optimise");
  declare_builtin_block("finalise");

  parser_init(fname,argc,argv);

  if (noexecute) { //scan through blocks
    scan_block("spinsys",true);
    scan_block("par");
    //    scan_block("initialise",true);
    scan_block("pulseq",true,true);
    scan_block("initialproc",true,true);
    scan_block("proc",true,true);
    scan_block("optimise",true);
    scan_block("finalise",true);
    return 0;
  }

  command_Factory_t& spinsys_Factory(initialise_spinsys_Factory());
  spinsys_Factory["channels"]=&parse_channels;
  //spinsys_Factory["variable"]=par_t(&parse_variable,true);
  //spinsys_Factory["function"]=par_t(&parse_function,true);
  //spinsys_Factory["verbose"]=&parse_verbose;
  spinsys_Factory["time_resolution"]=&parse_time_resolution;
  spinsys_Factory["transients"]=&parse_transients;
  
  (void)read_block("spinsys",spinsys_Factory,true);

  command_Factory_t& par_Factory(get_par_Factory());  
  par_Factory["start_operator"]=&parse_start_operator;
  par_Factory["detect_operator"]=&parse_detect_operator;
  par_Factory["spin_rate"]=&parse_spin_rate;
  par_Factory["gamma_angles"]=&parse_gamma_angles;
  par_Factory["gamma_zero"]=&parse_gamma_zero;
  //par_Factory["np"]=&parse_np;
  //par_Factory["sw"]=&parse_sw;
  par_Factory["rotor_angle"]=&parse_rotor_angle;
  par_Factory["log_file"]=par_t(&parse_log_file,true);
  //par_Factory["function"]=par_t(&parse_function,true);
  //par_Factory["puts"]=par_t(&parse_par_echo,0,true);
  //par_Factory["delay"]=par_t(&parse_delay,true);
  //par_Factory["variable"]=par_t(&parse_variable,true);
  //par_Factory["verbose"]=&parse_verbose;
  par_Factory["matrix"]=par_t(&parse_matrix,true);
  par_Factory["tolerance"]=par_t(&parse_tolerance,true);
  par_Factory["cache_limit"]=&parse_cache_limit;

  if (have_spinsys()) {
    if (!interactions_MFp)
      error_abort("spinsys block failed to defined nuclei / Hamiltonian");

    if (verbose & VER_GEN)
      dump_interactions();
    
    if (!nochecks && !(interactions_MFp->verify(std::cout,*cstructp,1e-1)))
      error_abort("Boundary conditions of periodic system incorrect");

    matrixmap["start"]=matrixspec(*(new operator_def(sigma0_specp,sigma0)));
    matrixmap["density"]=matrixspec(*(new operator_def(NMRSIM_NULL,density)));
    matrixmap["detect"]=matrixspec(*(new operator_def(detect_specp,detect)));
    matrixmap["hamiltonian"]=matrixspec(NMRSIM_NULL,matrixspec::SPECIAL);
  }
  
  make_par_variables();
  read_par_block(); //!< Always need some par set up e.g. np, so not optional
  NMRSIM_EXPECT(sum_dims.get(sum_ns,array_skips,sum_n0)==0); //array_skips is ignored
  process_array_ns();
  check_powderarray();
  //if (powder_array && (nacqdims!=1))
  //  error_abort("powder array cannot be combined with >1 dimensions");
  if (sum_ns.empty() && (sum_n0>1))
    sum_ns.create(1,sum_n0);

  sum_dims.lock(); //!< prevent subsequent changes
  array_dims.lock();
  ensure_array();

  //factory_read(initialisestack,"initialise",get_initialise_Factory(),Type2Type<ActionCommand>());
  size_t nactionstacks=0;
  
  if (have_spinsys()) {
//     if (!blockingnuclei.empty()) {
//       if (optblocking==0) 
// 	blockingnuclei.clear();
//       else {
// 	if ((optblocking==OPTIONAL_AUTO) && nondiagonal_opspec) {
// 	  offdiagonal_warning.raise();
// 	  blockingnuclei.clear();
// 	}
//       }	
//     }
    if (!optblocking) {
      blockingnuclei.clear();
      InversionGenerator::nucleusnotfound_warning.type(BaseWarning::Ignore); //!< suppress warning about nucleus not being blocked
    }
    if (verbose & VER_GEN) {
      cout << "Blocked nuclei: ";
      if (blockingnuclei.empty())
	cout << "none\n";
      else
	cout << blockingnuclei << '\n';
    }

    bool use_crystal=(havekey(spinsys_Factory,"cells") && optperiodic());
#ifdef NOPERIODIC
    if (use_crystal) {
      cout << "Program compiled without support for periodic systems - reverting to conventional calculation\n";
      use_crystal=false;
    }
#endif

    int flags=0;
    if (use_crystal && opteigsym())
      flags|=MetaFlags::UseEigSymmetry;
    if (optmzsym.isenabled()) {
      if (interactions_MFp->haslinear() || !(sysp->ishomonuclear()))
	mzsymmetry_warning.raise();
      flags|=MetaFlags::UseMzSymmetry;
    }
    if (optpartition())
      flags|=MetaFlags::UsePartitioning;

    if (isclassicQ())
      flags|=MetaFlags::ClassicSecondOrder;

    optperiodic.check(use_crystal);
    optmzsym.check(flags & MetaFlags::UseMzSymmetry);
    optpartition.check(flags & MetaFlags::UsePartitioning);
    if (use_crystal) 
      opteigsym.check(flags & MetaFlags::UseEigSymmetry);

    try {
      masterobjp=new MasterObj(*sysp,*interactions_MFp,blockingnuclei,flags,use_crystal);    
    }
    catch (Failed&) {
      if (blockingnuclei.empty())
	throw; //!< not expecting this!
      blockingnuclei.clear();
      buildfailed_warning.raise();
      masterobjp=new MasterObj(*sysp,*interactions_MFp,blockingnuclei,flags,use_crystal);    
    }
    optblocking.setusage(!(blockingnuclei.empty()));
    initialise_simulation_environment();
    if (sysp->ishomonuclear()) {
      const nuclei_spec nuc((*sysp)(0).nucleus());
      if (!detect_specp)
	detect_specp=new setableoperator_spec(operator_spec(nuc,'+'),*sysp);
      if (!sigma0_specp) {
	defaultsigma0_warning.raise(haveactions ? "Fz" : "Fx");
	sigma0_specp=new setableoperator_spec(operator_spec(nuc,haveactions ? 'z' : 'x'),*sysp);
      }
    }
    
    make_pulseq_variables();
    ensure_operator_matrices(true);

    rebuild_sequences(); //necessary before sequences can be used

    nactionstacks=multiple_factory_read(actionstacks,"pulseq",get_Action_Factory(),Type2Type<ActionCommand>(),Type2Type<ActionControl>());
    if (nactionstacks==0) {
      if (active2D)
	error_abort("nD sequences must define a pulseq block");
      
      actionstacks.create(1U);
      actionstacks.back().push_back(ActionDirectAcq::create(0.0)); //!< by default x phase detect
      if (verbose & VER_GEN)
	std::cout << "Creating default direct acquisition\n";
      actionstacks.back().initialise(); //!< tidy up
      nactionstacks++;
    }
    else {
      if (nactionstacks>1) {
	if (active2D)
	  error_abort("Can't define multiple pulseq blocks for nD acquisition");
	if (array_n0!=nactionstacks) {
	  std::cerr << "Number of pulse sequence blocks (" << nactionstacks << ") does not match rows in data set (" << array_n0 << ")\n";
	  error_abort();
	}
      }
    }
    const size_t acq_count=actionstacks.front().acq_count;
    for (size_t i=1;i<nactionstacks;i++) {
      if (actionstacks(i).acq_count!=acq_count)
	error_abort("Number of acquisition dimensions in different pulseq blocks must be the same");
    }
    if (acq_count!=nacqdims) {
      if (acq_count!=nacqdims-1) {
	std::cerr << "Number of acquisition dimensions (" << nacqdims << ") does not match number of acq's (" << acq_count << ")\n";
	error_abort();
      }
    }    
    if (verbose & VER_GEN) {
      for (size_t i=0;i<nactionstacks;i++) {	  
	cout << "Pulse sequence";
	if (nactionstacks>1)
	  cout << ' ' << (i+1);
	cout << ":\n" << actionstacks(i) << '\n';
      }
    }
  }
  else {
    if (verbose & VER_PARSE)
      std::cout << "Skipping parsing of pulseq (no spin system)\n";
    masterobjp=new MasterObj();
  }
  //  masterobjp->restart();

  if (verbose & VER_PARSE)
    parser_printcontext() << "Trying to parse proc blocks\n";

  bool haveproc=false;
  size_t nprocstacks=0;
  if (ActionDirectAcq::have_acq() || !have_spinsys()) {
    Process_Factory_t& Process_Factory(get_Process_Factory());
    //    if (ActionDirectAcq::have_acq() && read_proc_blocks(procstacks,nprocstacks,"initialproc",Process_Factory,true)) //!< initialproc only valid when used with DirectAcq
    if (read_proc_blocks(procstacks,nprocstacks,"initialproc",Process_Factory,true)) 
      haveproc=true;
    if (read_proc_blocks(postprocstacks,nprocstacks,"proc",Process_Factory,false))
      haveproc=true;
    if (!procstacks.empty() && optmergeproc()) {      
      bool donemerge=try_merge_proc_blocks(procstacks,postprocstacks);
      optmergeproc.check(donemerge);
    }
    if (have_spinsys()) {
      LIST<procstack_t>& usestacks(procstacks.empty() ? postprocstacks : procstacks);
      if (usestacks.empty())
	usestacks.create(1U);
      for (size_t i=usestacks.size();i--;)
	usestacks(i).push_front(new ProcessSignReverse(true));
    }

    if (verbose & VER_PARSE)
      parser_printcontext() << "Trying to parse optim block\n";

    parse_optimise_block(); //!< parse optional optimise block

    if (prefinalise_callback)
      (*prefinalise_callback)();
    Process_Factory_t& Finalise_Factory(get_Finalise_Factory());
    if (read_proc_blocks(finalisestacks,nprocstacks,"finalise",Finalise_Factory,false))
      haveproc=true;

    if (!haveproc && have_spinsys()) { //no explicit processing?
      sprintf(fname,isfrequency() ? "%s.spe" : "%s.fid",fnamebase);
      finalisestacks.create(1U);
      finalisestacks.front().push_back(new ProcessSave(fname));
      haveproc=true;
    }
    if (verbose & VER_GEN) {
      ldumpptr(procstacks,"Per-spectrum processing (initialproc)");
      ldumpptr(postprocstacks,"Processing (proc)");
      ldumpptr(finalisestacks,"Finalisation (finalise)");
    }
  }
  else {
    if (verbose & VER_PARSE)
      std::cout << "Skipping parsing of proc/finalise\n";
    if (have_spinsys())
      missingacq_warning.raise();
  }

  parser_flush();

  if (verbose & VER_GEN) {
    cout << "\nStandard variables:\n" << systemvarmap << '\n';
    if (!(vardefmap.empty()))
      cout << "User variable definitions:\n" << vardefmap << '\n';

    if (verbose_level>1) {
      if (!varpars.empty()) {
	cout << "\nArrayed/summed parameters: " << varpars.size() << '\n';
	dump(varpars);
      }

      if (!expressions.empty()) {
	cout << "\nExpressions:\n";
	dump(expressions,"\n");
      }
    }
  }

  if (!haveproc && !have_spinsys()) {
    nothingtodo_warning.raise();
    return 0;
  }

  //! sanity checks before starting calculation
  if (!nochecks && have_spinsys() && interactions_MFp && havekey(par_Factory,"crystal_file")) {
    const Euler testangle(0.3,0.6,0); //non-special angle
    if (!masterobjp->check_powder(*interactions_MFp,testangle))
      error_abort(ERR_FAILED);
    
    if (spin_rate && (!have_spinsys() || interactions_MFp->empty()))
      nointeractions_butMAS_warning.raise();
  }
  
  startup_cp->leave(); //<! finished initialisation

  if (masterobjp==NMRSIM_NULL)
    throw InternalError("masterobjp is NULL");

  if (precalculation_callback)
    (*precalculation_callback)(*masterobjp);
  
  if (!testedoptions.empty()) {
    testoptions(*masterobjp);    
    exit(0);
  }

  DataStore Spec;

  if (!calculation_callback)
    calculation_callback=default_calculation;
  const bool ok=(*calculation_callback)(*masterobjp,Spec);
  if (!ok)
    return ERR_FAILED;

  if (ammaster) {
    evaluation_state=CONTEXT_FINALISE;
    apply_procstacks(finalisestacks,Spec);
    evaluation_state=CONTEXT_TERMINATE;
  }

  flush_calculation();

  if (verbose & VER_PROFILE) {
    if (global_workers)
      parprofile_warning.raise();

    if (ammaster) {
      std::cout.precision(4);
      dump_timings();
    }
  }
  if (ammaster) {
    if (verbose & (VER_GEN | VER_PROFILE))
      dump_optimisationusage();

    if (!nochecks)
      UserVariable::check_unused();

    if (verbose & VER_GEN) {
      const size_t totwarn=lcm_base_warning.count()+lcm_serious_warning.count()+NMRsim_repeat_warning.count()+NMRsim_once_warning.count();
      std::cout << "\nTotal warnings: " << totwarn << " \tMajor: " << lcm_serious_warning.count() << '\n';
    }
  }

  } catch (const MatrixException& exc) {
    cerr << exc << '\n';
    if (!cleanup_parallel())
      std::cerr << "Warning: failed to cleanup parallel system - this may leave stuck processes\n";
    return ERR_FAILED;
  }
  return 0;
}
