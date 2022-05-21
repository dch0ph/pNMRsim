/*! \file
  \brief  components shared between pNMRsim-base programs 
*/

#include "NMRsim.h"
#include "Parser.h"
#include "NMRsim_logfile.h"
#include "Lineshapes.h"
#include <sstream>

#define NMRSIM_DEFAULT_LOGFILE_FORMAT ASCII

option optcache("cache","cache propagators",option::AUTO,option::NOTUSED);
option optpartition("mzpartitioning","mz partitioning");
option optupdate("minimalupdating");

void option::setdefault(optional_t def)
{
  if (value_==AUTO)
    value_=def;
}

void add_option(option& opt)
{
  optional_map_t& optional_map(get_optional_map());  
  optional_map[opt.name()]=&opt;
}

double sw=0.0; //!< spectral width (0 if unset)
int np=0; //!< number of acquisition points
int verbose=0; //!< verbose factor

LIST<const char*> testedoptions;

SystemVariable<int*> v_np("np",&np,V_ISCONST);
SystemVariable<double*> v_sw("sw",&sw,1.0,V_ISCONST);
static bool donenp=false;
static bool donesw=false;

bool logfile_active() { return (get_logfile()!=NMRSIM_NULL); }
inline bool allowoutput() { return (!silent || (get_logfile()!=NMRSIM_NULL)); }

flagsmap_type verbose_flags;

spell_dictionary_type& get_spell_dictionary()
{
  static spell_dictionary_type spell_dictionary;
  return spell_dictionary;
}

bool check_unrecognised(const char* com, const char* name)
{
  parser_printcontext() << (name ? name : "directive") << " not recognised: " << com << std::endl;
  spell_dictionary_type& spell_dict(get_spell_dictionary());
  if (spell_dict.empty()) {
    spell_dict["logfile"]="log_file";
    spell_dict["auto_opt"]="autoopt";
    spell_dict["minimize"]="minimise";
    spell_dict["maximize"]="maximise";
    spell_dict["addsignal"]="addsignals";
  }
  const spell_dictionary_type::const_iterator iter(spell_dict.find(com));
  if (iter!=spell_dict.end()) {
    parser_printthread(std::cerr) << "Correct to " << (iter->second) << '?' << std::endl;
    return true;
  }
  return false;
}

#ifndef NDEBUG
int def_verbose=VER_GEN | VER_POWDER | VER_OPTIM | VER_PARSE | VER_PARALLEL;
#else
int def_verbose=VER_GEN | VER_POWDER | VER_OPTIM | VER_PARALLEL;
#endif

const flagsmap_type& getverboseflags()
{
  static bool doneaddmap=false;
  if (!doneaddmap) {
    verbose_flags["general"]=VER_GEN;
    verbose_flags["powder"]=VER_POWDER;
    verbose_flags["optimise"]=VER_OPTIM;
    verbose_flags["parse"]=VER_PARSE;
    verbose_flags["profile"]=VER_PROFILE;
    verbose_flags["parallel"]=VER_PARALLEL;
    doneaddmap=true;
  }
  return verbose_flags;
}

ContextWarning<> logfileactiveclash_warning("opening log file while another log file with the same name is active - this is almost certainly an error. Close one before re-opening.",&NMRsim_repeat_warning);
ContextWarning<> logfileoverwrite_warning("opening log file with same name as one used previously, which is expected to overwrite the first. Add -append flag?",&NMRsim_repeat_warning);

void checklogfile(const char* fname, bool otheractive, bool append)
{
  static std::string lastfname;
  if (strcmp(fname,lastfname.c_str())==0) {
    if (otheractive)
      logfileactiveclash_warning.raise();
    else {
      if (!append)
	logfileoverwrite_warning.raise();
    }
  }
  lastfname=fname;
}

static void output_line(const char* buffer)
{
  logfile_controller* logfilep(get_logfile());
#ifdef HAVE_FORK_CONTROLLER
  if ((evaluation_state==CONTEXT_MAINLOOP) && global_workers) {
    static Mutex<> lock; //!< should be safe using static - Mutex only created when function called (only within multi-threaded region)
    MutexLock<> lockguard(lock);
    const size_t thrnum=get_thread_num();
    if (!logfilep)
      std::cout << 'T' << (thrnum+1) << ": " << buffer << '\n';
    else {
      char qualecho[16];
      snprintf(qualecho,sizeof(qualecho),"echoT%" LCM_PRI_SIZE_T_MODIFIER "u",thrnum+1);
      char title[32];
      logfile_controller::maketitle(title,sizeof(title),qualecho);
      logfilep->write(buffer,title);
    }
    return;
  }
#endif
  if (!logfilep)
    std::cout << buffer << '\n';
  else {
    char title[16];
    logfile_controller::maketitle(title,sizeof(title),"echo");
    logfilep->write(buffer,title);
  }
}
 
void parse_verbose()
{
  verbose=parse_flags(getverboseflags());
  if (!verbose)
    verbose=def_verbose;
}

#if NMRSIM_USE_HASH==0
bool FunctionCompare::operator()(const function_spec& f1, const function_spec& f2) const
{
  if (f1.second!=f2.second)//!< sort first on argument number since this is quicker
    return (f1.second<f2.second);
  return (strcmp(f1.first,f2.first)<0);
}
#endif

std::ostream& prettyprint_interval(double t1, double t2, std::ostream& ostr)
{
  prettyprint_time(t1,ostr) << " to ";
  return prettyprint_time(t2,ostr);
}

void parse_precision()
{
  const int prec=parse_int();
  const int matrixprec=parser_isnormal() ? parse_int() : prec;
  if ((prec<-1) || (matrixprec<-1))
    error_abort("precision must be >=-1");
  if (prec>=0) {
    std::cout.precision(prec);
    std::cerr.precision(prec);
  }    
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["complexpair"]=ostream_controller::pair;
    flags["complexi"]=ostream_controller::withi;
    flags["complexcompact"]=ostream_controller::compact;
  }
  ostream_controller& ctrl(cmatrix_ostream_controller()); //!< only need to modify one (default)
  if (matrixprec>=0)
    ctrl.matrixprecision=matrixprec;
  if (are_left())
    ctrl.complexview=ostream_controller::complexview_t(parse_flags(flags,ostream_controller::withi,true));
}

void parse_par_echo(int isecho)
{
  if (!allowoutput()) {
    set_curline(NMRSIM_NULL); //!< need to swallow rest of line
    return;
  }
  char buffer[MAXLINE];
  substitute_string(buffer,sizeof(buffer),isecho ? get_curline() : get_token(),SUB_ESCAPE);
  if (isecho)
    set_curline(NMRSIM_NULL);
  output_line(buffer);
}

void parse_np()
{
  parse_system_variable(v_np,F_defaultdataset);
  if ((np<1) && v_np.isconstant())
    error_abort("np cannot be <1");
  donenp=true;
}

void parse_sw()
{
  parse_system_variable(v_sw,F_defaultdataset); //can't allow sw to differ between summed spectra
  if ((sw<=0.0) && v_sw.isconstant())
    error_abort("sw cannot be <=0!");
  donesw=true;
}

void make_data_variables()
{
  command_Factory_t& par_Factory(get_par_Factory());
  par_Factory["echo"]=par_t(&parse_par_echo,1,true);
  par_Factory["precision"]=par_t(&parse_precision,true);
  par_Factory["np"]=&parse_np;
  par_Factory["sw"]=&parse_sw;
  add_systemvarmap(v_np);
  add_systemvarmap(v_sw);
}

systemvarmap_type systemvarmap;

SystemVariable<int*>::SystemVariable(const std::string& namev, int* value_, int flagsv)
  : SystemVariableBase(namev), 
    RealVariable(namev.c_str(),!(flagsv & V_ISFIXED),(flagsv & V_ISCONST), (flagsv & V_POWDERSUMDEP) ? A_SUM : 0),
    value(value_)
{
  if (flagsv & ~(V_ISFIXED | V_ISCONST | V_POWDERSUMDEP))
    error_abort("Unprocessed flags in SystemVariable<int*>");  
}

SystemVariable<size_t*>::SystemVariable(const std::string& namev, size_t* value_, int flagsv)
  : SystemVariableBase(namev),
    RealVariable(namev.c_str(),!(flagsv & V_ISFIXED),(flagsv & V_ISCONST), (flagsv & V_POWDERSUMDEP) ? A_SUM : 0),
    value(value_)
{
  if (flagsv & ~(V_ISFIXED | V_ISCONST | V_POWDERSUMDEP))
    error_abort("Unprocessed flags in SystemVariable<size_t*>");  
}

SystemVariable<double*>::SystemVariable(const std::string& namev, double* value_, double scalef_, int flagsv)
  //variables are const unless altered by parse_system_variable
  : SystemVariableBase(namev),
    RealVariable(namev.c_str(),!(flagsv & V_ISFIXED),(flagsv & V_ISCONST), (flagsv & V_POWDERSUMDEP) ? A_SUM : 0), 
    value(value_), scalef(scalef_),
    updateints(flagsv & V_UPDATEINTS)
{
  if ((flagsv & V_ISFIXED) && updateints)
    throw InvalidParameter("SystemVariable: fixed and update are incompatible");
  if (flagsv & ~(V_ISFIXED | V_ISCONST | V_UPDATEINTS | V_POWDERSUMDEP))
    error_abort("Unprocessed flags in SystemVariable<double*>");

  if (scalef_==0.0)
    throw InvalidParameter("SystemVariable: zero scale factor!");
}
  
void SystemVariable<double*>::set(double fval, subsid_t subsid)
{
  if (!issetable())
    throw InternalError("Attempt to change fixed variable");
  assert(subsid==S_ARG1); //!< assert is harmless
  fval/=scalef;
  if (*value!=fval) {
    *value=fval;
    update();
    if (updateints)
      update_interactions=true;
  }
}

void SystemVariable<int*>::set(double fval, subsid_t subsid)
{
  if (!issetable())
    throw InternalError("Attempt to change fixed variable");
  assert(subsid==S_ARG1); //!< assert is harmless
  const int newval=round_int(fval);
  if (*value!=newval) {
    *value=newval;
    update();
  }
}

void SystemVariable<size_t*>::set(double fval, subsid_t subsid)
{
  if (!issetable())
    throw InternalError("Attempt to change fixed variable");
  assert(subsid==S_ARG1); //!< assert is harmless
  const int newval=round_int(fval);
  if (newval<0) {
    parser_printthread(std::cerr) << "Attempt to set " << SystemVariableBase::name() << " to negative value\n";
    error_abort();
  }
  if (*value!=newval) {
    *value=newval;
    update();
  }
}

void printattributes(std::ostream& ostr, int uses)
{
  ostr << '[';
  bool done=false;
  if (uses & A_SUM) {
    ostr << "sum";
    done=true;
  }
  if (uses & A_ARRAY) {
	ostr << (done ? "," : "") << "array";
	done=true;
  }
  if (uses & A_VAR) {
    ostr << (done ? "," : "") << "fit-variable";
    done=true;
  }
  if (!done)
    ostr << "<none>";
  ostr << ']';
}

namespace {
  void dumptype(std::ostream& ostr, const RealVariable& var, const char* type)
  {
    //    ostr << '(' << (var.isconstant() ? "const " : "") << type << (var.isused() ? ")" : ", unused)");
    ostr << '(' << (var.isconstant() ? "const " : "") << type << ')';
    const int uses=var.uses();
    if (uses) {
      ostr << ' ';
      printattributes(ostr,uses);
    }      
  }
}

void SystemVariable<double*>::print(std::ostream& ostr) const
{
  ostr << SystemVariableBase::name() << '=' << (scalef*get()) << ' ';
  dumptype(ostr,*this,"real");
}

void SystemVariable<int*>::print(std::ostream& ostr) const
{
  ostr << SystemVariableBase::name() << '=' << (*value) << ' ';
  dumptype(ostr,*this,"integer");
}

void SystemVariable<size_t*>::print(std::ostream& ostr) const
{
  ostr << SystemVariableBase::name() << '=' << (*value) << ' ';
  dumptype(ostr,*this,"positive integer");
}

void SystemVariable<std::string>::print(std::ostream& ostr) const 
{
  ostr << name_ << '=' << value << " (string)";
}

void add_systemvarmap(SystemVariableBase& var)
{
  const char* name=var.name_.c_str();
  if (vardefmap.find(name)!=vardefmap.end()) {
    parser_printcontext() << "System variable clashes with previously defined user variable: " << name << '\n';
    error_abort();
  }
  systemvarmap[var.name_.c_str()]=&var;
}

bool isirregular()
{
  if (!v_np.isconstant() || !v_sw.isconstant()) {
    if (active2D)
      throw InternalError("isirregular"); //!< sanity check
    return true;
  }
  return false;      
}


SimpleCallbackStack& get_prepar_stack()
{
  static SimpleCallbackStack prepar_stack;
  return prepar_stack;
}

void addprepar_callback(simple_callback_t callback)
{
  get_prepar_stack().push(callback);
}

void SimpleCallbackStack::exec() const
{
  for (size_t i=stack.size();i--;)
    (*stack(i))();
}

void read_par_block(bool opt)
{
  if (verbose & VER_PARSE)
    parser_printcontext() << "Trying to parse par block\n";

  get_prepar_stack().exec();
  (void)read_block("par",get_par_Factory(),opt);
  if (!donenp)
    try_find_np();
  if (!donesw)
    try_find_sw();

  if (verbose & VER_PARSE)
    parser_printcontext() << "Finished parsing of par block\n";
}

static LIST< smartptr<logfile_controller,false> > logfile_stack;

logfile_controller* get_logfile()
{
  return logfile_stack.empty() ? NMRSIM_NULL : logfile_stack.back().get();
}

logfile_controller::composite_guard::~composite_guard()
{
  logcon.close_composite();
}

ThreadWarning<> logfile_controller::unexpectedcompositeclose_warning("logfile_controller: attempt to close non-existent composite object",&NMRsim_repeat_warning);

void logfile_controller::close_composite()
{
  if (filep_composite==NMRSIM_NULL)
    unexpectedcompositeclose_warning.raise();
  if (!(filep_composite->ok_to_close())) {
    parser_printthread(std::cerr) << "Can't close structured object - no items written? File will be corrupt\n";
    error_abort();
  }
  delete filep_composite;
  filep_composite=NMRSIM_NULL;
}

void logfile_controller::open_composite(const char* name)
{
  if (filep_matlab==NMRSIM_NULL)
    throw Failed("logfile_controller::open_composite can only be used for Matlab output format");

  if (filep_composite)
    throw Failed("logfile_controller: can't recursively create structured objects");

  filep_composite=new matlab_controller::composite(*filep_matlab,name,matlab_controller::STRUCT);
}

logfile_controller::composite_guard::composite_guard(logfile_controller& logconv, const char* name)
  : logcon(logconv)
{
  logcon.open_composite(name);
}

logfile_controller::logfile_controller(const char* name, int flagsv) 
  : flags(flagsv), filep_composite(NMRSIM_NULL)
{
  if (flags & MATLAB) {
    int mcflags=0;
    if (!(flags & DOUBLE))
      mcflags|=matlab_controller::singleprecision;
    if (flags & APPEND)
      mcflags|=matlab_controller::append;

    filep_matlab=new matlab_controller(name,5,mcflags);
    filep_ascii=NMRSIM_NULL;
  }
  else {
    filep_matlab=NMRSIM_NULL;
    filep_ascii=fopen(name,(flags & APPEND) ? "a" : "w");
    if (!filep_ascii)
      throw Failed("logfile_controller: couldn't open file");
  }
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Opened log file: " << name << '\n';
}

static int logcount=0;
static int logsupercount=0;

void logfile_controller::maketitle(char* dest, int n, const char* source)
{
  if (logsupercount)
    snprintf(dest,n,"%s_%i_%i",source,logsupercount,++logcount);
  else
    snprintf(dest,n,"%s_%i",source,++logcount);
}

void reset_logcount()
{
  logcount=0;
  logsupercount++;
}

logfile_controller::~logfile_controller()
{
  if (filep_matlab) {
    if (filep_composite)
      delete filep_composite;
    delete filep_matlab;
  }
  else
    fclose(filep_ascii);
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Closed log file\n";
}

void logfile_controller::write(const char* str, const char* name) const
{
  if (!try_write_matlab(str,name)) {
    fputs(str,filep_ascii);
    fputc('\n',filep_ascii);
  }
}

void logfile_controller::flush() const
{
  fflush(filep_matlab ? filep_matlab->file_pointer() : filep_ascii);
}

template<class T> void logfile_controller::write(const BlockedMatrix<T>& a, const char* name) const
{
  if (!try_write_matlab(a,name)) {
    int outflags=mxflag::block;
    if (flags & DOUBLE)
      outflags|=mxflag::doublep;
    for (size_t i=0;i<a.size();i++)
      write_matrix(filep_ascii,a(i),name,outflags);
  }
}

template<typename T> void write_ascii_list(FILE* fp, const BaseList<T>& curlist, const char* formatstr)
{
  for (size_t c=0;c<curlist.size();c++) {
    if (c)
      fputc(' ',fp);
    fprintf(fp,formatstr,curlist(c));
  }
  fputc('\n',fp);
}

template<> void logfile_controller::write(const ListList<size_t>& a, const char* name) const
{
  if (!try_write_matlab(a,name)) {
    fprintf(filep_ascii,"%%%s\n",name);
    for (size_t r=0;r<a.size();r++)
      write_ascii_list(filep_ascii,a(r),"%" LCM_PRI_SIZE_T_MODIFIER "u");
  }
}

template<typename T> struct output_traits {};
template<> struct output_traits<size_t> { static const char* formatstr; };
const char* output_traits<size_t>::formatstr="%" LCM_PRI_SIZE_T_MODIFIER "u";
template<> struct output_traits<double> { static const char* formatstr; };
const char* output_traits<double>::formatstr="%g";

#ifdef NEED_SEPARATE_COUNT_T
template<> struct output_traits<unsigned int> { static const char* formatstr; };
const char* output_traits<unsigned int>::formatstr="%u";
#endif

template<typename T> void logfile_controller::write(const BaseList<T>& curlist, const char* name) const
{
  if (!try_write_matlab(curlist,name)) {
    fprintf(filep_ascii,"%%%s\n",name);
    write_ascii_list(filep_ascii,curlist,output_traits<T>().formatstr);
  }
}

template<> void logfile_controller::write(const BaseList<complex>& curlist, const char* name) const
{
  if (!try_write_matlab(curlist,name)) {
    fprintf(filep_ascii,"%%%s\n",name);
    write_vector(filep_ascii,curlist,mxflag::rmat); //!< output R C
  }
}
  
template<class T> void logfile_controller::write(const Matrix<T>& a, const char* name) const
{
  if (!try_write_matlab(a,name)) {
    int outflags=mxflag::block;
    if (flags & DOUBLE)
      outflags|=mxflag::doublep;
    write_matrix(filep_ascii,a,name,outflags);
  }
}

template void logfile_controller::write(const BlockedMatrix<double>&, const char*) const;
template void logfile_controller::write(const BlockedMatrix<complex>&, const char*) const;
template void logfile_controller::write(const BlockedMatrix<bool>&, const char*) const;
template void logfile_controller::write(const ListList<size_t>&, const char*) const;
template void logfile_controller::write(const BaseList<size_t>&, const char*) const; //!< slightly problematic for Matlab
#ifdef NEED_SEPARATE_COUNT_T
template void logfile_controller::write(const BaseList<count_t>&, const char*) const;
#endif
template void logfile_controller::write(const BaseList<double>&, const char*) const;
template void logfile_controller::write(const BaseList<complex>&, const char*) const;

template void logfile_controller::write(const Matrix<complex>&, const char*) const;

SimpleCallbackStack& getflushstack() 
{
  static SimpleCallbackStack flushstack;
  return flushstack;
}

void add_flushcallback(simple_callback_t proc)
{
  getflushstack().push(proc);
}

void flush_calculation()
{
  logfile_stack.clear();
  getflushstack().exec();
}

void logfile_controller::parse_raw(char*& name, int& log_flags)
{
  static flagsmap_type llog_flags;
  if (llog_flags.empty()) {
    llog_flags["double"]=DOUBLE;
    llog_flags["matlab"]=MATLAB;
    llog_flags["ascii"]=ASCII;
    llog_flags["append"]=APPEND;
  }
  if (are_left()) {
    name=parse_string(F_REPLACEDOLLAR);
    log_flags=parse_flags(llog_flags);
    switch (log_flags & (MATLAB | ASCII)) {
    case 0:
      log_flags|=NMRSIM_DEFAULT_LOGFILE_FORMAT;
      break;
    case MATLAB: case ASCII:
      break;
    default:
      error_abort("can't combine -matlab and -ascii flags");
    }
  }
  else
    name=NMRSIM_NULL;
}

ContextWarning<> logfile_controller::ignoring_warning("logging commands will be ignored during powder averaging",&NMRsim_once_warning);
ContextWarning<> logfile_controller::badclose_warning("log_file close requested with no active file",&NMRsim_repeat_warning);

void update_log_file(const char* name, int flags)
{
  if (!ammaster) {//!< silently ignore if not master
    logfile_controller::ignoring_warning.raise();//!< NB logging outside of powder averaging will raise this warning even though harmless (and not actually ignored by master)
    return; 
  }
  if (name && name[0]) {
    //    log_filep.clear(); //close any existing log_file
    //log_filep.reset(new logfile_controller(name,flags));
    if (!nochecks)
      checklogfile(name,transitions_log_active() || logfile_active(),flags & logfile_controller::APPEND); //!< crude check for overlapping log_file usage
    logfile_stack.push_back();
    logfile_stack.back().reset(new logfile_controller(name,flags));
  }
  else {
    if (logfile_stack.empty())
      logfile_controller::badclose_warning.raise();
    else
      logfile_stack.pop_back();
  }      
}

logfile_controller::logfile_controller(const logfile_controller&)
{
  throw InternalError("logfile_controller object cannot be copied");
}

//bool disablelogging=false;

void parse_log_file()
{
  char* name;
  int log_flags;
  logfile_controller::parse_raw(name,log_flags);
  update_log_file(name,log_flags);
}

LIST<size_t> losttrans; //!< total 'lost' transitions
LIST<complex> lostsum; //!< total 'lost'
LineshapeSpec lineshape_spec;
double hist_lw=0.0;
double hist_gfrac=0.0;
int hist_flags=0;
double hist_min=0.0;
double hist_max=0.0;
LorentzianGaussian* lineshape_genp=NMRSIM_NULL;
ThreadWarning<> HistogramMaker::loss_warning("lost histogram intensity: ",&NMRsim_repeat_warning);

namespace {
  template<class T> void checkloss(const T& Spec)
  {
    const size_t lost=Spec.lost();
    if (lost) {
      if (losttrans.empty()) {
	losttrans.create(array_n0,size_t(0));
	lostsum.create(array_n0,complex(0.0));
      }
      losttrans+=lost;
      lostsum+=Spec.lostsum();
      if ((verbose & VER_GEN) && HistogramMaker::loss_warning.enabled()) {
	std::ostringstream str(std::ostringstream::out);
	str << "Transitions: " << lost << "  Total intensity: " << Spec.lostsum();
	HistogramMaker::loss_warning.raise(str.str().c_str());
      }
    }
  }
}

HistogramMaker::~HistogramMaker()
{
  if (!!histp_)
    checkloss(*histp_);
}

ThreadWarning<> HistogramMaker::foldrange_warning("-fold combined with specified histogram range",&NMRsim_once_warning);

HistogramMaker::HistogramMaker(BaseList<complex> dest, double lsw)
{
  const size_t lnp=dest.size();
  if (lnp==0)
    throw InvalidParameter("HistogramMaker");
  //create histogram
  double hist_range,hist_left;
  if (hist_min!=hist_max) { //range specified
    hist_range=hist_max-hist_min;
    hist_left=hist_min;
    if (hist_flags & HIST_FOLD)
      foldrange_warning.raise();
  }
  else {
    if (lsw==0.0)
      throw Failed("Histogram (unspecified range) requires specification of spectral width");
    hist_range=lsw;
    hist_left=(-lsw-lsw/lnp)/2.0;
  }
  if (lineshape_genp)
    histp_.reset(new LineshapeHistogram<complex>(dest,*lineshape_genp,hist_range,hist_left));
  else { //! no lineshape
    switch (hist_flags & (HIST_INTERP | HIST_FOLD)) {
      //   case (HIST_INTERP | HIST_FOLD): {
      //     FoldingInterpHistogram<complex> Spec(FIDSpec,hist_range,hist_left);
      //     add_freq(Spec,cscale,liscomplex);
      //   }
      //     break;
      //   case HIST_INTERP: {
      //     InterpHistogram<complex> Spec(FIDSpec,hist_range,hist_left);
      //     add_freq(Spec,cscale,liscomplex);
      //     checkloss(Spec);
      //   }
      //     break;
    case HIST_FOLD:
      histp_.reset(new FoldingHistogram<complex>(dest,hist_range,hist_left));
      break;
    case 0:
      histp_.reset(new Histogram<complex>(dest,hist_range,hist_left));
      break;
    default:
      throw InternalError("HistogramMaker: invalid histogram mode");
    }
  }
}

void update_lineshapegen()
{  
  if (hist_lw) {
    if (lineshape_genp) {
      lineshape_genp->lw(hist_lw);
      lineshape_genp->fraction(hist_gfrac);
    }
    else
      lineshape_genp=new LorentzianGaussian(hist_lw,hist_gfrac,lineshape_spec);
  }
  else {
    if (lineshape_genp) { //!< kill lineshape generator
      delete lineshape_genp;
      lineshape_genp=NMRSIM_NULL;
    }
  }
}

const char* SystemVariable<FUNC_REAL>::format()
{ 
  try {
    snprintf_prec(buf,sizeof(buf),scalef*get());
  } catch (...) {
    return "<undefined>";
  }
  return buf;
}

void SystemVariable<FUNC_REAL>::print(std::ostream& ostr) const 
{
  ostr << SystemVariableBase::name() << '=';
  try {
    ostr << (scalef*get()) << " (function)";
  } catch (...) {
    ostr << "<undefined in this context>";
  }
}

// void SystemVariable<double*>::print(std::ostream& ostr, bool full) const
// {
//   if (full)
//     this->print(ostr);
//   else {
//     this->printvariablename(ostr);
//     ostr << '=' (get()*scalef);
//   }
// }

// void SystemVariable<int*>::print(std::ostream& ostr, bool full) const
// {
//   this->print(ostr,full ? S_NONE : S_ARG1);
// }

// void SystemVariable<size_t*>::print(std::ostream& ostr, bool full) const
// {
//   this->print(ostr,full ? S_NONE : S_ARG1);
// }

// void SystemVariable<FUNC_REAL>::print(std::ostream& ostr, bool full) const
// {
//   this->print(ostr,full ? S_NONE : S_ARG1);
// }

command_Factory_t& get_par_Factory() { 
  static command_Factory_t par_Factory;
  return par_Factory;
}

void update_expressions()
{
  const varpars_t::iterator end=expressions.end();
  varpars_t::iterator start=expressions.begin();
  while (start!=end) {
    start->update(false); //!< suppress warnings since updating everything
    ++start;
  } 
}   

/**
   \pre quantity pointer must be set (otherwise \c Undefined exception) */
VariableBase& Variable::variable()
{
  if (!valuep)
    throw Undefined("Variable");
  return *valuep;
}

/**
   \pre quantity pointer must be set (otherwise \c Undefined exception) */
const VariableBase& Variable::variable() const 
{
  if (!valuep)
    throw Undefined("Variable");
  return *valuep;
}

void Variable::rawupdate()
{
  const BaseList<double>& val(variable().value());
  if (verbose & VER_GEN) {
    std::cout << "Setting ";
    ptr->printvariablename(std::cout,subsid);
    std::cout << " to ";
    if (val.size()==1)
      std::cout << val.front();
    else
      std::cout << val;
    std::cout << '\n';
  }
  if (val.size()==1) 
    ptr->set(val.front(),subsid);
  else
    ptr->set(val,subsid);
}

void Variable::update(bool showwarnings)
{
  VariableBase& var(variable());
  if (!var.isconst()) {//may get called for constant expressions (ignore)
    var.update(showwarnings);
    rawupdate();
  }
}

varpars_t expressions;
varpars_t varpars;

int BaseEcho::update() const
{
  int uses=0;
  substitute_string(buffer,sizeof(buffer),str.c_str(),SUB_ESCAPE,uses);
  return uses;
}

void BaseEcho::exec() const
{ 
  if (allowoutput()) {
    update();
    output_line(buffer);
  }
}

bool checkflags(bool& fl, const char* com, const char* name)
{
  if (com[0]!='-' || strcmp(com+1,name))
    return false;
  fl=true;
  return true;
}

bool option::isstate(optional_t v) const
{
  if (value_==TEST)
    throw InternalError("Caught flag in TEST state");
  return (value_==v);
}

bool option::isnotstate(optional_t v) const
{
  if (value_==TEST)
    throw InternalError("Caught flag in TEST state");
  return (value_!=v);
}

void checkoptflags(const std::string com, option::optional_t v, bool ignoreunknown, bool allowtest)
{
  optional_map_t& optional_map(get_optional_map());
  const optional_map_t::iterator end(optional_map.end());
  const optional_map_t::iterator curp(optional_map.find(com));
  if (curp==end) {
    if (ignoreunknown)
      return;
    parser_printthread(std::cerr) << "Unrecognised option flag: " << com << "\nAllowed values are:";
    optional_map_t::iterator start(optional_map.begin());
    while (start!=end) {
      std::cerr << ' ' << (start->first);
      ++start;
    }
    std::cerr << '\n';
    exit(ERR_INVALID_INPUT);
  }
  (curp->second)->set(v);
  if (v==option::TEST) {
    if (allowtest) {
      const char* str=(curp->first).c_str();
      if (std::find(testedoptions.begin(),testedoptions.end(),str)==testedoptions.end()) //!< have to check as options may be checked twice
	testedoptions.push_back(str);
    }
    else
      error_abort("-test:<flag> feature not available");
  }
}

void option::setusage(bool v)
{
  if (v && (value_==OFF))
    throw InvalidParameter("Can't use disabled feature");
  used_=v ? USED : NOTUSED;
}

 
option::option(const char* labelv, const char* descv, optional_t def, used_t defusage)
  : label_(labelv), desc_(descv), 
    value_(def), used_(defusage) {
  if (defusage==USED)
    throw InvalidParameter("option: cannot be initialised to used");
}
 
void option::check(bool val, const char* optdesc)
{
  const char* desc=optdesc;
  if (!desc)
    desc=desc_.c_str();
  if (*desc=='\0')
    desc=label_.c_str();

  //  assert(get_optional_map().find(name)->second==&request); // if debugging check that name matches map
  if (value_==TEST) {
    parser_printthread(std::cerr) << "Flag " << label_ << " can't (currently) be tested.\n";
    error_abort();
  }
  setusage(val);
  if (verbose & VER_GEN)
    std::cout << "Using " << desc_ << " (" << label_ << "): " << (val ? "Yes\n" : "No\n");  
  if ((value_==AUTO) || ((value_==ON)==val))
    return;
  parser_printthread(std::cerr) << "Warning: " << ((value_==OFF) ? "dis" : "en") << "able:" << label_ << " request was not fulfilled\n";
}

bool checkverboseflag(const char* com)
{
  if (strcmp(com,"-verbose")==0) {
    verbose=def_verbose;
    return true;
  }
  if (STRNCMP(com,"-verbose:")!=0)
    return false;
  verbose|=parse_flag(getverboseflags(),com+9);
  return true;
}
      
bool checkoptflag(const char* com, bool& havedisable, bool ignoreunknown, bool allowtest)
{
  if (*com++!='-')
    return false;
  if (STRNCMP(com,"enable:")==0) {
    checkoptflags(com+7,option::ON, ignoreunknown, allowtest);
    return true;
  }
  if (STRNCMP(com,"disable:")==0) {
    checkoptflags(com+8,option::OFF, ignoreunknown, allowtest);
    havedisable=true;
    return true;
  }
  if (STRNCMP(com,"test:")==0) {
    checkoptflags(com+5,option::TEST, ignoreunknown, allowtest);
    return true;
  }
  return false;
}

void setdefaultoptions(option::optional_t state)
{
  optional_map_t& optional_map(get_optional_map());
  const optional_map_t::iterator end(optional_map.end());
  optional_map_t::iterator start(optional_map.begin());
  while (start!=end) {
    (start->second)->set(state);
    ++start;
  }
}

//!< global override of warnings
void setwarnings(BaseWarning::warning_t type)
{
  lcm_base_warning.type(type);
  NMRsim_repeat_warning.type(type);
  NMRsim_once_warning.type(type);
}

ThreadWarning<> library_version("Potential incompatibility between libcmatrix versions",&NMRsim_repeat_warning);

void ensure_library_version(const char* needversion)
{
  if (strcmp(cmatrix_abi_version,needversion)<0) {
    char buf[256];
    snprintf(buf,sizeof(buf)," Found V%s, but expecting V%s",cmatrix_abi_version,needversion);
    library_version.raise(buf);
  }
}

//! greatest common divisor (Euclid's algorithm)
int gcd(int a, int b)
{
  return b ? gcd(b, a % b) : a;
}

int lcm(int a, int b)
{
  const int top=a*b;
  if (top==0)
    throw InvalidParameter("lcm: zero input");
  return top/gcd(a,b);
}
