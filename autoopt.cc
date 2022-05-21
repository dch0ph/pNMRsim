// "Module" for autoopt

#include "NMRsim_MasterObj.h"
#include "NMRsim_Process.h"
#include "Parser.h"

void parse_autoopt();

//! base for auto-optimisation classes
template<class T> class auto_opt_base {
protected:
  T& var; //!< reference to variable optimised
public:
  auto_opt_base(T& var_) : var(var_) {}
  virtual ~auto_opt_base() {}
  T get() const { return var; } //!< return current value
  void set(T val) { var=val; } //!< set value
  virtual T& operator++ () =0; //!< 'next' value
};

//! step by increment
template<class T> class auto_opt_incr : public auto_opt_base<T> {
private:
  const T incr; //!< increment size
public:
  using auto_opt_base<T>::var;
  auto_opt_incr(T& var_, T incr_);
  T& operator++ () { var+=incr; return var; } //!< apply increment
};

//! step by multiplication factor
template<class T> class auto_opt_scale : public auto_opt_base<T> {
private:
  const T scale; //!< scale factor
public:
  using auto_opt_base<T>::var;
  auto_opt_scale(T& var_, T scale_);
  T& operator++ () { var*=scale; return var; } //!< apply scale
};

//! auto_opt enable flags
enum autoopt_t {
  A_NULL=0, //!< dummy type
  A_POWDER=1, //!< auto_opt number of powder orientations
  A_GAMMA_ANGLES=2, //!< auto_opt number of gamma angle steps
  A_MAXDT=4, //!< auto_opt integration timestep
  A_CHEBYSHEV=8, //!< auto_opt Chebyshev iterations
  A_MAXJUMPDT=16, //!< auto_opt maxjumpdt
  A_STEPS_PER_ROTATION=32 //!< auto_opt steps_per_rotation
};

struct autoopt_item {
  autoopt_item(autoopt_t whichv, double chistopv, int flagsv)
    : which(whichv), chistop(chistopv), flags(flagsv) {}
  
  bool isreset() const { return (flags!=0); }
  
  autoopt_t which;
  double chistop;
  int flags;    
};

//! auto-optimise using supplied function object  
template<class T> T auto_opt(MasterObj& Obj, auto_opt_base<T>& optobj, T startval, const autoopt_item&, const char* name, T maxval =T(0));

namespace {
  precalculation_callback_t old_precalculation_callback=NMRSIM_NULL;
  VariableBase* autostorep=NMRSIM_NULL;
  LIST<autoopt_item> autoopt_stack;
}

struct Autoopt_Proxy_ {
  Autoopt_Proxy_() {
    command_Factory_t& par_Factory(get_par_Factory());
    par_Factory["autoopt"]=par_t(&parse_autoopt,true);
  }
};

//declare Autoopt additions
static const Autoopt_Proxy_ autooptproxy_;

template<class T> auto_opt_scale<T>::auto_opt_scale(T& var_, T scale_)
  : auto_opt_base<T>(var_), scale(scale_) {
  if (scale_<=T(0) || scale_==T(1))
    throw InvalidParameter("auto_opt_scale: scale factor can't be <0 or 1!");
}

template<class T> auto_opt_incr<T>::auto_opt_incr(T& var_, T incr_)
  : auto_opt_base<T>(var_), incr(incr_) {
  if (incr==T(0))
    throw InvalidParameter("auto_opt_incr: increment can't be zero!");
}

static DataStore buffers[2];

template<class T> T auto_opt(MasterObj& Obj, auto_opt_base<T>& optobj, T startval, const autoopt_item& curopt, const char* name, T maxval)
{
  DataStore* lastbufp=buffers+0;
  DataStore* curbufp=buffers+1;
  LIST<complex> diff;
  double normfac=0.0;
  //  const double chi2stop=(curopt.chistop)*(curopt.chistop);
  const T initval=optobj.get();
  if (maxval && (initval>maxval))
    throw InvalidParameter("auto_opt: initial value exceeds maximum!");
  if (initval!=startval)
    optobj.set(startval);
  T prevval;
  size_t iter=0;
  double chi2;
  const bool doreset=curopt.isreset();
  bool abort=false;
  for (;!abort;iter++) {
    ::std::swap(curbufp,lastbufp);
    if (!(curbufp->empty()))
      *curbufp=complex(0.0);
    Obj.calc(*curbufp,0);
    if (iter) {
      diff=curbufp->row();
      diff-=lastbufp->row();
      if ((verbose & VER_OPTIM) && (verbose_level>1))
	std::cout << lastbufp->row().front() << "  " << curbufp->row().front() << "  " << diff.back() << '\n';
      chi2=norm(diff)/normfac;
      const double sqrtchi2=sqrt(chi2);
      if (!silent)
	std::cout << "Fractional difference (compared to previous iteration) with " << name << "=" << optobj.get() << ": " << sqrtchi2 << std::endl;
      if (NMRSIM_ISNAN(sqrtchi2))
	error_abort("Calculation of chi2 failed - optimised parameter is probably out of range");
      if (sqrtchi2<curopt.chistop)
	break;
    }
    else {
      normfac=norm(curbufp->row());
      if (normfac==0)
	error_abort("Can't optimise empty spectrum!");
    }
    prevval=optobj.get();
    ++optobj;
    if (maxval && (optobj.get()>maxval)) {
      if (!silent)
	std::cout << "Optimisation of " << name << " stopped at " << prevval << " (next increment would exceed limit)\n";
      abort=true;
    }
  }

  if (!abort && !silent) {
    std::cout << name << " optimised at " << prevval << '\n';
    if (iter==1 && (chi2<1e-8))
      std::cerr << "Warning: are you sure spectrum depends on " << name << " or is starting point already at or beyond optimum?\n";
  }
  optobj.set(doreset ? initval : prevval);
  return prevval;
}

//instantiate auto_opt templates
template class auto_opt_incr<int>;
template class auto_opt_incr<size_t>;
template class auto_opt_scale<double>;
template int auto_opt(MasterObj&, auto_opt_base<int>&, int, const autoopt_item&, const char* name, int);
template double auto_opt(MasterObj&, auto_opt_base<double>&, double, const autoopt_item&, const char* name, double);
template size_t auto_opt(MasterObj&, auto_opt_base<size_t>&, size_t, const autoopt_item&, const char* name, size_t);

namespace {
  const char autooptsyntax[]="autoopt <stop fraction> powderquality|gamma_angles|maxdt|maxjumpdt|chebyshev_iterations|steps_per_rotation [-reset]";

  autoopt_t name_to_autoopt(const char* name)
  {
    if (strcmp(name,"powderquality")==0)
      return A_POWDER;
    if (strcmp(name,"gamma_angles")==0)
      return A_GAMMA_ANGLES;
    if (strcmp(name,"maxdt")==0)
      return A_MAXDT;
    if (strcmp(name,"steps_per_rotation")==0)
      return A_STEPS_PER_ROTATION;
    if (strcmp(name,"maxjumpdt")==0)
      return A_MAXJUMPDT;
    if (strcmp(name,"chebyshev_iterations")==0)
      return A_CHEBYSHEV;
    parser_printcontext() << "not an autopt variable: " << name << '\n';
    error_abort(autooptsyntax);
    return A_NULL; //!< dummy to suppress compilation warning
  }
}

ThreadWarning<> silent_warning("-silent flag makes no sense with autoopt - disabled",&NMRsim_repeat_warning);

void testoptions(MasterObj& Obj)
{
  const size_t nflags = testedoptions.size();
  const size_t ntests=1 << nflags;
  DataStore& usebuf(buffers[0]);
  cmatrix alldata;
  optional_map_t& optmap(get_optional_map());

  LIST<option*> optvals(nflags);
  for (size_t i=nflags;i--;) {
    const optional_map_t::iterator iter(optmap.find(testedoptions(i)));
    if (iter==optmap.end())
      throw InternalError("testoptions: flag has disappeared!");
    optvals(i)=iter->second;
  }

  std::string origname;
  SystemVariable<std::string>* strpos=NMRSIM_NULL;
  const systemvarmap_type::iterator namepos(systemvarmap.find("name"));
  if (namepos!=systemvarmap.end()) {
    strpos=dynamic_cast< SystemVariable<std::string>* >(namepos->second);
    if (strpos)
      origname=(*strpos)();
  }
  const bool havename=!(origname.empty());
  char scratch[512];
  const bool silentstatus=silent;
  if ((verbose & VER_OPTIM)==0) {
    silent=true; //!< suppress output
    setwarnings(BaseWarning::Silent);
  }

  double normfac=0.0;  

  std::cout << '\n';
  for (size_t test=0; test<ntests; test++) {
    std::cout << (test+1) << ':';
    if (test==0)
      std::cout << " <none>";
    size_t mask=1;
    char offset=0;
    if (havename)
      offset=snprintf(scratch,sizeof(scratch),"%s",origname.c_str());

    for (size_t i=0;i<nflags;i++) {
      const bool ison=(test & mask);
      const option::optional_t val = ison ? option::ON : option::OFF;
      if (ison) {
	std::cout << ' ' << testedoptions(i);
	if (havename)
	  offset+=snprintf(scratch+offset,sizeof(scratch)-offset,"_%s",testedoptions(i));
      }
      optvals(i)->set(val);
      mask <<=1;
    }
    if (havename)
      strpos->set(scratch); //!< over-ride $name
    if (!usebuf.empty())
      usebuf=0.0;
    timer<> stopwatch;
    Obj.calc(usebuf,0);
    const double timetaken=stopwatch();
    std::cout << " \tTime: " << (timetaken*1e3) << " ms\n";    
    if (normfac==0.0)
      normfac=norm(usebuf.row());
    if (alldata.empty())
      alldata.create(ntests,usebuf.row().size());
    BaseList<complex> crow(alldata.row(test));
    crow=usebuf.row();
  }
  if (normfac==0.0)
    error_abort("Can't compare options on empty spectrum!");

  const double threshold=1e-4;
  std::cout << "\nNormalised differences (should be of order of numerical errors):\n";
  for (size_t c=0;c<ntests-1;c++) 
    std::cout << " \t \t" << (c+1);
  std::cout << '\n';
  LIST<complex> diff;
  bool problem=false;
  for (size_t r=1;r<ntests;r++) {
    std::cout << (r+1) << ": ";
    for (size_t c=0;c<r;c++) {
      diff=alldata.row(r);
      diff-=alldata.row(c);
      const double sqrtchi2=sqrt(norm(diff)/normfac);
      std::cout << " \t" << sqrtchi2;
      if (sqrtchi2>threshold)
	problem=true;
    }
    std::cout << '\n';
  }
  if (problem)
    std::cout << "\nOh dear! At least one pair of calculations differ by more than tolerance of " << threshold << ".\n";
  else
    std::cout << "\nSuccess! All calculations match within tolerance of " << threshold << ".\n";

  if (havename)
    strpos->set(origname);
  silent=silentstatus;
}    

void autoopt_callback(MasterObj& masterobj)
{
  if (old_precalculation_callback)
    (*old_precalculation_callback)(masterobj);

  if (autoopt_stack.empty())
    return;

  //  if (!testedoptions.empty())
  //  error_abort("Can't sensibily combine -test:<flags> with autoopt. Best to test flags first, then optimise simulation parameters.");

  LIST<double> finalvals;
  
  for (size_t opt_i=0;opt_i<autoopt_stack.size();opt_i++) {
    const autoopt_item& curopt(autoopt_stack(opt_i));
    
    double finalval=0;
    switch (curopt.which) {
    case A_POWDER: {
      if (!havepowderquality())
	error_abort("Can't optimise this powder type");      
      auto_opt_incr<int> optobj(nzcw,1);
      zcwisindex=true; //!< ensure nzcw is always treated as index
      finalval=auto_opt(masterobj,optobj,nzcw,curopt,"powderquality");
    }
      break;
    case A_GAMMA_ANGLES: {
      const int incr=(nobs<2) ? 2 : nobs;
      auto_opt_incr<int> optobj(gamma_angles,incr);	
      finalval=auto_opt(masterobj,optobj,gamma_angles ? gamma_angles : incr,curopt,"gamma_angles");
    }
      break;
    case A_STEPS_PER_ROTATION: {
      set_inttype(INT_ROTORPERIOD);
      auto_opt_scale<int> optobj(steps_per_rotation,2);
      finalval=auto_opt(masterobj,optobj,steps_per_rotation,curopt,"steps_per_rotation");
    }
      break;
    case A_MAXDT: {
      set_inttype(INT_EXPLICIT);
      auto_opt_scale<double> optobj(maxdtv,0.5);
      finalval=auto_opt(masterobj,optobj,maxdtv,curopt,"maxdt");
    }
      break;
    case A_CHEBYSHEV: {
      auto_opt_incr<size_t> optobj(cmatrix_eigensystem_controller.chebyshev_iterations,1);
      finalval=auto_opt(masterobj,optobj,cmatrix_eigensystem_controller.chebyshev_iterations,curopt,"chebyshev_iterations",eigensystem_controller::max_chebyshev_iterations);
    }
      break;
    case A_MAXJUMPDT: {
      if (maxjumpdt==0)
	error_abort("can't optimise maxjumpdt if is zero (disabled!)");
      auto_opt_scale<double> optobj(maxjumpdt,0.5);
      finalval=auto_opt(masterobj,optobj,maxjumpdt,curopt,"maxjumpdt");
    }
      break;
    default:
      throw InternalError("unknown autoopt");
    }
    finalvals.push_back(finalval);
  }
  if (autostorep)
    autostorep->set(finalvals);
  //VariableBase& store(*(new VariableBase(finalvals)));
  //UserVariable::create("autoopt_results",store,true);
}

void parse_autoopt()
{
  if (!are_left())
    error_abort(autooptsyntax);

  static flagsmap_type auto_flags;

  if (auto_flags.empty()) {
    auto_flags["reset"]=1;

    if ((verbose & VER_OPTIM) && (verbose_level>1))
      std::cout << "Inserting autoopt callback\n";
    old_precalculation_callback=precalculation_callback;
    precalculation_callback=autoopt_callback;

    VariableBase* autostorep=new VariableBase();
    (void)UserVariable::create("autoopt_results",*autostorep,UserVariable::IGNORE_UNUSED);
  }

  const double chistop=parse_double();
  if (chistop<=0.0 || chistop>=1.0)
    error_abort("autoopt stop parameter must be >0 and <1");

  static LIST<autoopt_t> optvars;
  optvars.create(0U);
  while (parser_isnormal())
    optvars.push_back(name_to_autoopt(parse_string(F_REPLACEDOLLAR)));
  const int flags=parse_flags(auto_flags);
  static int autoopt_resetlist=0;
  static ContextWarning<> repeatautoopt_warning("request for autoopt on previously reset variable",&NMRsim_repeat_warning);
  for (size_t i=0;i<optvars.size();i++) {
    const autoopt_t which(optvars(i));
    if (autoopt_resetlist & which)
      repeatautoopt_warning.raise();
    const autoopt_item item(which,chistop,flags);
    if (item.isreset())
      autoopt_resetlist|=which;
    autoopt_stack.push_back(item);
  }
}
