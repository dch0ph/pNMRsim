#include "NMRsim_Process.h"
#include "Parser.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

int F_defaultdataset=F_SIMPLE;
double detect_freq=0.0; //!< don't know detection frequency

bool havesave=false;
bool debug=false;

#ifdef NDEBUG
#define VERBOSE_LEVEL 2
#else
#define VERBOSE_LEVEL 3
#endif

const char* NEWSstr=
#include "NEWS"

void printversion()
{
  unsigned short day,month,year;
  if (sscanf(NEWSstr,"%hu/%hu/%hu",&day,&month,&year)!=3)
    throw Failed("Failed to parse data from NEWS");
  printf("%02u.%02u.%02u\n",year,month,day);
}

ProcessCommand* HookedSave_create() {
  havesave=true;
  return ProcessSave::create(); //!< not a true hook!
}

void checkdomain(domain_t& dom, const char* name)
{
  static ThreadWarning<> domainunclear_warning("time/frequency domain could not be determined from data set - defaulting to time domain (use setdomain to set explicitly).  Domain: \n",&NMRsim_repeat_warning);

  if (dom==D_UNKNOWN) {
    std::cerr << "For " << name << " dimension: ";
    domainunclear_warning.raise(name);
    dom=D_TIME;
  }
}

void checksw(double& sw, double datasw,const char* name)
{
  static ThreadWarning<> swmismatch_warning("overriding spectral width specified in data set with a different explicitly specified value.",&NMRsim_repeat_warning);

  if (sw) {
    if (datasw && (fabs(datasw-sw)>1e-4*sw)) {
      char errm[256];
      snprintf(errm,sizeof(errm)," Dimension: %s  Data set: %g kHz, explicit sw: %g kHz",name,datasw*1e-3,sw*1e-3);
      swmismatch_warning.raise(errm);
      datasw=sw;
    }
  }
  else
    sw=datasw;
}

int global_argc=0;
char** global_argv;

int main(int argc_,char **argv_)
{
  global_argc=argc_;
  global_argv=argv_;

  try {

  int count=1;
  bool randomise=false;
  bool havedisable=false;
  bool noexecute=false;
  setdefaultoptions();

  while (count<argc_) {
    const char* com=argv_[count];
    if ((com[0]!='-') || (com[1]=='\0'))
      break;
    count++;

    if (strcmp(com,"-version")==0) {
      printversion();
      return 0;
    }
    
    if (!checkflags(randomise,com,"randomise") && !checkflags(debug,com,"debug") && !checkflags(abortonwarning,com,"abort") && !checkflags(noexecute,com,"noexecute") && !checkflags(silent,com,"silent")) {
      if (!checkoptflag(com,havedisable) && !checkverboseflag(com)) {
	cerr << "Unknown flag: " << com << '\n';
	return ERR_INVALID_INPUT;
      }
    } 
  }

  ProcessSave::defaultstopoverwrite=!(optsaveprotect.isdisabled()); //!< by default allow overwrite

  const char syntaxstring[]="Syntax: pNMRproc [-randomise|-debug|-abort|-silent|-enable:<option>|-disable:<option>|-noexecute] <.in file> <data file> [<arguments>]";

  if (randomise)
    set_seed(); //!< reset random number generator
  if (count==argc_)
    error_abort(syntaxstring);

  verbose_level = debug ? VERBOSE_LEVEL : 1;
  if (debug && noexecute)
    verbose=VER_PARSE; //!< if debugging file read, turn on VER_PARSE
  if (debug)
    setwarnings(BaseWarning::Always); //!< if debugging leave warnings turned on
  else {
    if (silent)
      setwarnings(BaseWarning::Silent);
  }
  
  char* fname=argv_[count++];
  if (count>=argc_)
    error_abort(syntaxstring);

  const char* pname=argv_[count];
  char pnamebase[256]="";
  strncat(pnamebase,pname,sizeof(pnamebase)-1); //!< strncat ensures proper termination (although filename is truncated)
  char* dot=strrchr(pnamebase,'.');
  if (dot && (dot>strrchr(pnamebase,'/')))
    *dot='\0'; // strip off terminator
  systemvarmap["name"]=new SystemVariable<std::string>("name",pnamebase);

  size_t argc=argc_-count-1;
  char** argv=argv_+count+1;

  if (argc && (argv[0][0]=='-'))
    std::cerr << "Warning: $1 begins with - (" << argv[0] << ").  Misplaced flag?\n";

  declare_builtin_block("par");
  declare_builtin_block("proc");

  List<cmatrix> datasets(1U);
  filestruct data;
  raw_read_file(datasets.front(),data,pname);
  const bool is2D=(datasets.front().rows()>1);
  if (is2D)
    set_n(1,datasets.front().rows());

  checkdomain(data.domain,"direct");
  sw=data.sw; //!< set if possible
  if (is2D) {
    checkdomain(data.domain1,"indirect");
    sws.front()=data.sw1;
  }

  parser_init(fname,argc,argv);

  if (noexecute) { //scan through blocks
    scan_block("par");
    scan_block("proc",true,true);
    return 0;
  }

  make_data_variables();

  const size_t lnp=datasets.front().cols();
  np=lnp;
  v_np.isconstant(true);
  
  command_Factory_t& par_Factory(get_par_Factory());
  par_Factory["sw1"]=par_t(&parse_swn,1,true); //!< just add sw1

  read_par_block(true);
  if (np!=lnp)
    error_abort("Cannot change np - this is determined from data set");
  
  checksw(sw,data.sw,"direct");
  if (sw==0.0)
    error_abort("Spectral width not determined from data set.  Must be set explicitly in par block with sw");
  if (is2D)
    checksw(sws.front(),data.sw1,"indirect");

  NMRSIM_EXPECT(sum_dims.get(sum_ns,array_skips,sum_n0)==0); //array_skips is ignored
  if (sum_n0>1)
    error_abort("Sum || arrays cannot be used");
  process_array_ns();
  
  LIST<processing_state> dataprocs;
  if (is2D && (ndims==2))
    dataprocs.push_back(processing_state(sws.front(),data.domain1==D_TIME));
  dataprocs.push_back(processing_state(sw,data.domain==D_TIME, data.sfrq, data.ref));

  ThreadWarning<> nDunhandled_warning(">2D data sets not fully supported - indirect dimension information is being ignored",&NMRsim_once_warning);
  if (ndims>3)
    nDunhandled_warning.raise();
      
//   if (sum_ns.empty() && (sum_n0>1))
//     sum_ns.create(1,sum_n0);

//   sum_dims.lock(); //!< prevent subsequent changes
  array_dims.lock();
  ensure_array();
      
  if (verbose & VER_GEN) {
    cout << "\nStandard variables:\n" << systemvarmap << '\n';
    if (!(vardefmap.empty()))
      cout << "User variable definitions:\n" << vardefmap << '\n';
  }

  DataStore Spec;
  raw_build_set(Spec,datasets,dataprocs,(verbose & VER_GEN) ? verbose_level : 0);

  Process_Factory_t Process_Factory(get_Process_Factory());
  Process_Factory["save"]=HookedSave_create; //!< insert hook so we can detect save
  procstack_t procstack;
  if (!read_proc_block(procstack,"proc",Process_Factory,false))
    error_abort("Failed to find proc block");

  if (verbose & VER_GEN) {
    std::cout << "Processing to be applied:\n";
    dumpptr(procstack,"\n");
  }

  parser_flush();

  ThreadWarning<> nosave_warning("no save command in processing. This is probably a mistake!",&NMRsim_repeat_warning);
  if (!havesave)
    nosave_warning.raise();

  procstack.exec(Spec);

  } catch (const MatrixException& exc) {
    cerr << exc << '\n';
    return ERR_FAILED;
  }
  return 0;
}
