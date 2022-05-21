#include "expression_definition.hpp"
#include <cmath>
#include <sstream>
#include "NMRsim.h"
#include "Parser.h"

static const double deg_to_rad=M_PI/180.0;
static const double rad_to_deg=180.0/M_PI;

namespace {

  //! absolutely ensure close-to-zero is zero
  inline double ensure_zero(double x)
  { return (fabs(x)<1e-8) ? 0.0 : x; }

  size_t getindex(double v, size_t max) {
    const int ind=round_int(v);
    if ((ind<1) || (ind>max))
      throw BadIndex("Index out of range",ind,max);    
    return ind-1;
  }
}

Function_Factory_t& get_Function_Factory()
{
  static Function_Factory_t Function_Factory;
  return Function_Factory;
}

double ExpressionFunctionBase::checknan(double v) const
{
  if (NMRSIM_ISNAN(v))
    throw InvalidParameter(name_);
  return v;
}

double ExpressionFunctionBase::checkhuge(double v) const
{
  if (v==HUGE_VAL)
    throw Failed(name_);
  return v;
}

double ExpressionFunctionBase::checkneghuge(double v) const
{
  if (v==-HUGE_VAL)
    throw Failed(name_);
  return v;
}

namespace {
  void getmn(int& m, int& n, const BaseList<double>& in, const char* name)    
  {
    m=round_int(in.front());
    n=round_int(in(size_t(1)));
    if ((m!=0) && (n!=0))
      return;
    char buf[256];
    snprintf(buf,sizeof(buf),"arguments to %s must be non-zero integers",name);
    throw InvalidParameter(buf);
  }
}

EXPR_DECLARE_FUNCTION(FunctionGcd,"gcd");
void FunctionGcd::operator()(LIST<double>& dest, const BaseList<double>& in) const {  
  int m,n;
  getmn(m,n,in,"gcd");
  dest.push_back(gcd(m,n));
}

EXPR_DECLARE_FUNCTION(FunctionLcm,"lcm");
void FunctionLcm::operator()(LIST<double>& dest, const BaseList<double>& in) const {  
  int m,n;
  getmn(m,n,in,"lcm");
  dest.push_back(lcm(m,n));
}

size_t check_index(double x)
{
  const int n=round_int(x);
  if (n<1)
    error_abort("index into data set cannot be <1");
  return n-1;
}

void getrealimag(LIST<double>& dest, const BaseList<double>& in, bool isreal)
{
  if (current_data_row.empty())
    error_abort("Can't use getreal/imag outside processing context (no current data set)");

  size_t start=check_index(in.front());
  size_t end=start;
  switch (in.size()) {
  case 1:
    break;
  case 2:
    end=check_index(in(size_t(1)));
    if (start>end)
      ::std::swap(start,end);
    break;
  default:
    throw InternalError("getrealimag");
  }
  if (end>=current_data_row.size()) {
    parser_printthread() << "getreal/imag index (" << end << ") outside range of data set (" << current_data_row.size() << " points)\n";
    error_abort();
  }
  const size_t npts=end-start+1;
  dest.reserve(dest.size()+npts);
  if (isreal) {
    for (size_t i=start;i<=end;i++)
      dest.push_back(real(current_data_row(i)));
  }
  else {
    for (size_t i=start;i<=end;i++)
      dest.push_back(imag(current_data_row(i)));
  }
}
  
EXPR_DECLARE_FUNCTION_USEARRAY(FunctionGetReal,"getreal");
EXPR_DECLARE_FUNCTION_USEARRAY(FunctionGetImag,"getimag");

void FunctionGetReal::operator()(LIST<double>& dest, const BaseList<double>& in) const
{
  getrealimag(dest,in,true);
}  
void FunctionGetImag::operator()(LIST<double>& dest, const BaseList<double>& in) const
{
  getrealimag(dest,in,false);
}  

EXPR_DECLARE_FUNCTION(FunctionSin,"sin");
void FunctionSin::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(sin(in.front()*deg_to_rad)); }

EXPR_DECLARE_FUNCTION(FunctionCos,"cos");
void FunctionCos::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(cos(in.front()*deg_to_rad)); }

EXPR_DECLARE_FUNCTION(FunctionTan,"tan");
void FunctionTan::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(tan(in.front()*deg_to_rad)); }

EXPR_DECLARE_FUNCTION(FunctionAsin,"asin");
void FunctionAsin::operator()(LIST<double>& dest, const BaseList<double>& in) const
{ dest.push_back(checknan(asin(in.front()))*rad_to_deg); }

EXPR_DECLARE_FUNCTION(FunctionAcos,"acos");
void FunctionAcos::operator()(LIST<double>& dest, const BaseList<double>& in) const
{ dest.push_back(checknan(acos(in.front()))*rad_to_deg); }

EXPR_DECLARE_FUNCTION(FunctionAtan,"atan");
void FunctionAtan::operator()(LIST<double>& dest, const BaseList<double>& in) const
{ 
  double res;
  switch (in.size()) {
  case 1:
    res=atan(in.front());
    break;
  case 2:
    res=atan2(in.front(),in(size_t(1)));
    break;
  default:
    throw InternalError("atan");
  }
  dest.push_back(res*rad_to_deg); 
}

EXPR_DECLARE_STRING_FUNCTION(FunctionError,"Error");
bool FunctionError::rawisconstant() const { return false; } //!< always false to prevent premature optimisation
void FunctionError::operator()(LIST<double>&) const
{
  error_abort(args_.front().c_str());
}

struct FunctionWarn : public ExpressionStringFunction {	
  explicit FunctionWarn(int nargsv, bool oncev)
    : ExpressionStringFunction(oncev ? "WarnOnce" : "Warn",nargsv), once_(oncev) {}

  ExpressionBase* clone() const { return new FunctionWarn(*this); }
  bool rawisconstant() const { return false; }
  void operator()(LIST<double>& dest) const { 
    dest.create(size_t(0)); //!< note string functions need to set result (not accumulated)
    if (args_.size()>1)
      warnp->raise(args_(size_t(1)).c_str());
    else
      warnp->raise();
  }
  void create(BaseList<std::string>& args) {
    if (args.front().size()==0)
      error_abort("Warn: warning message is empty!");
    ExpressionStringFunction::create(args);
    BaseWarning* warnbasep(once_ ? &NMRsim_once_warning : &NMRsim_repeat_warning);
    warnp.reset(new ThreadWarning<>(args_.front().c_str(),warnbasep));
  }
  bool once_;
  smartptr<Warning<>, false> warnp;
};

struct FunctionIndicesof : public ExpressionStringFunction {	
  explicit FunctionIndicesof(int nargsv)
    : ExpressionStringFunction("Indicesof",nargsv), nuc(NULL_NUCLEUS) {}

  ExpressionBase* clone() const { return new FunctionIndicesof(*this); }
  bool rawisconstant() const { return (sysp!=NMRSIM_NULL); } //!< will block evaluation until spin system defined

  void operator()(LIST<double>& dest) const;
  void create(BaseList<std::string>& args) {
      nuc=args.empty()
	? NULL_NUCLEUS 
	: parse_nucleusname(args.front().c_str());
  }
  size_t nuc;
};

  void FunctionIndicesof::operator()(LIST<double>& dest) const
  { 
    if (!sysp)
      throw Failed("Cannot evaluate Indicesof before spin system is defined");
    const size_t n=sysp->nspins();
    if (nuc==NULL_NUCLEUS) {
      dest.create(n);
      for (size_t i=n;i--;)
	dest(i)=i+1;
    }
    else {
      dest.reserve(n);
      dest.create(size_t(0));
      for (size_t i=0;i<n;i++) {
	if ((*sysp)(i).nucleus()==nuc)
	  dest.push_back(double(i+1));
      }
    }
  }

// EXPR_DECLARE_FUNCTION(FunctionSinh,"sinh");
// void FunctionSinh::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(sinh(in.front())); }

// EXPR_DECLARE_FUNCTION(FunctionCosh,"cosh");
// void FunctionCosh::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(cosh(in.front())); }

// EXPR_DECLARE_FUNCTION(FunctionTanh,"tanh");
// void FunctionTanh::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(tanh(in.front())); }

namespace {

  typedef double (*P_RANDFUNC)(double,bool);

  void dorandom(LIST<double>& dest, const BaseList<double>& in, P_RANDFUNC func)
  {
    const double sigma=in.front();
    const size_t n=(in.size()>1) ? round_int(in.back()) : 1;
    dest.reserve(dest.size()+n);
    for (size_t i=n;i--;)
      dest.push_back(func(sigma,false));
  }
}

EXPR_DECLARE_FUNCTION(FunctionRandom,"random");
void FunctionRandom::operator()(LIST<double>& dest, const BaseList<double>& in) const
{ dorandom(dest,in,libcmatrix::random); }

EXPR_DECLARE_FUNCTION(FunctionNoise,"noise");
void FunctionNoise::operator()(LIST<double>& dest, const BaseList<double>& in) const
{ dorandom(dest,in,libcmatrix::gauss); }

EXPR_DECLARE_FUNCTION(FunctionExp,"exp");
void FunctionExp::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(checkhuge(exp(in.front()))); }

EXPR_DECLARE_FUNCTION(FunctionLn,"ln");
void FunctionLn::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(checkneghuge(log(in.front()))); }

// EXPR_DECLARE_FUNCTION(FunctionLog10,"log10");
// void FunctionLog10::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(checkneghuge(log10(in.front()))); }

EXPR_DECLARE_FUNCTION(FunctionAbs,"abs");
void FunctionAbs::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(fabs(in.front())); }

EXPR_DECLARE_FUNCTION(FunctionSign,"sign");
void FunctionSign::operator()(LIST<double>& dest, const BaseList<double>& in) const {
  const double v=in.front();
  double res=0.0;
  if (ensure_zero(v))
    res=(v<0.0) ? -1.0 : 1.0;  
  dest.push_back(res);
}

EXPR_DECLARE_FUNCTION(FunctionRound,"round");
void FunctionRound::operator()(LIST<double>& dest, const BaseList<double>& in) const 
{
  double res;
  switch (in.size()) {
  case 1:
    res=ensure_zero(floor(in.front()+0.5));
    break;
  case 2: {
    const double y=in(size_t(1));
    if (y)
      res=y*ensure_zero(floor(0.5+in.front()/y));
    else
      throw Failed("round: divide by zero");
  }
    break;
  default:
    throw InternalError("FunctionRound");
  }
  dest.push_back(res);
}

EXPR_DECLARE_FUNCTION(FunctionCeil,"ceil");
void FunctionCeil::operator()(LIST<double>& dest, const BaseList<double>& in) const 
{
  double res;
  switch (in.size()) {
  case 1:
    res=ensure_zero(ceil(in.front()));
    break;
  case 2: {
    const double y=in(size_t(1));
    if (y)
      res=y*ensure_zero(ceil(in.front()/y));
    else
      throw Failed("ceil: divide by zero");
  }
    break;
  default:
    throw InternalError("FunctionCeil");
  }
  dest.push_back(res);
}

EXPR_DECLARE_FUNCTION(FunctionFloor,"floor");
void FunctionFloor::operator()(LIST<double>& dest, const BaseList<double>& in) const 
{
  double res;
  switch (in.size()) {
  case 1:
    res=ensure_zero(floor(in.front()));
    break;
  case 2: {
    const double y=in(size_t(1));
    if (y)
      res=y*ensure_zero(floor(in.front()/y));
    else
      throw Failed("floor: divide by zero");
    break;
  }
  default:
    throw InternalError("FunctionFloor");
  }
  dest.push_back(res);
}

EXPR_DECLARE_FUNCTION(FunctionSqrt,"sqrt");
void FunctionSqrt::operator()(LIST<double>& dest, const BaseList<double>& in) const { dest.push_back(checknan(sqrt(in.front()))); }

const argdef_t LSSdef(ExpressionGeneralFunction::makeislist(ScratchList<bool>(true,false,false)));
const argdef_t LLdef(ExpressionGeneralFunction::makeislist(ScratchList<bool>(true,true)));
const argdef_t Ldef(ExpressionGeneralFunction::makeislist(ExplicitList<1,bool>(true))); //!< ScratchList no good here as single argument form works differently...
static const argdef_t LSdef(ExpressionGeneralFunction::makeislist(ScratchList<bool>(true,false)));
static const argdef_t SLLdef(ExpressionGeneralFunction::makeislist(ScratchList<bool>(false,true,true)));
static const argdef_t LLLdef(ExpressionGeneralFunction::makeislist(ScratchList<bool>(true,true,true)));

//EXPR_DECLARE_GENERAL_FUNCTION(FunctionNorm,"norm",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionSize,"size",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionRepeat,"repeat",LSdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionHead,"head",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionHead2,"head",LSdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionTail,"tail",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionTail2,"tail",LSdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionLast,"last",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionLast2,"last",LSdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionRev,"rev",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionIf,"if",SLLdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionSortIndex,"sortindex",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionReplace,"replace",LLLdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionSum,"sum",Ldef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionStdDev,"stddev",Ldef);

EXPR_DECLARE_MODIFIER_FUNCTION(FunctionMutable,"mutable",false);
EXPR_DECLARE_MODIFIER_FUNCTION(FunctionConst,"const",true);

EXPR_DECLARE_MULTIVARIATE_FUNCTION(FunctionSwitch,"switch");
EXPR_DECLARE_MULTIVARIATE_FUNCTION(FunctionSwitchDefault,"switchdefault");

EXPR_DECLARE_GENERAL_FUNCTION(FunctionEcho,"echo",Ldef);
void FunctionEcho::operator()(LIST<double>& dest) const
{
  dest=getarg(size_t(0));
  dumparray(std::cout,dest) << '\n';
}

// EXPR_DECLARE_FUNCTION(FunctionFailed,"failed");
// void FunctionFailed::operator()(LIST<double>&, const BaseList<double>&) const { throw Failed("expression cannot be evaluated"); }

EXPR_DECLARE_GENERAL_FUNCTION(FunctionBadParameter,"badparameter",Ldef);
void FunctionBadParameter::operator()(LIST<double>&) const 
{ 
  std::ostringstream ostr(std::ostringstream::out);
  dumparray(ostr,getarg(size_t(0)));
  throw InvalidParameter(ostr.str().c_str());
}

void FunctionIf::operator()(LIST<double>& dest) const
{
  const double cond=getarg(size_t(0)).front();
  if (ensure_zero(cond))
    dest=getarg(size_t(1));
  else
    dest=getarg(size_t(2));
}

void FunctionSwitch::operator()(LIST<double>& dest) const
{
  const BaseList<double> switchvals(getarg(size_t(0)));
  const int max=nargs()-1;
  if (max<2)
    error_abort("switch() has less than three arguments!");
  for (size_t i=0;i<switchvals.size();i++) {
    const size_t switchval(getindex(switchvals(i),max));
    dest.push_back(getarg(1+switchval));
  }
}

void FunctionSwitchDefault::operator()(LIST<double>& dest) const
{
  const BaseList<double> switchvals(getarg(size_t(0)));
  const int max=nargs()-2;
  if (max<1)
    error_abort("switchdefault() has less than three arguments!");
  for (size_t i=0;i<switchvals.size();i++) {
    const int rind=round_int(switchvals(i));
    const size_t switchval=((rind<1) || (rind>max)) ? max+1 : rind-1;
    dest.push_back(getarg(1+switchval));
  }
}

// void FunctionNorm::operator()(LIST<double>& dest) const
// {
//   double s=0.0;
//   const BaseList<double>& in(getarg(size_t(0)));
//   for (size_t i=in.size();i--;)
//     s+=in(i)*in(i);
//   dest.create(1,s);
// }

void FunctionSum::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  const double v= a.empty() ? 0.0 : sum(a);
  dest.create(1,v);
}
  
void FunctionStdDev::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  const size_t npts=a.size();
  if (npts<2)
    error_abort("stddev called on list with less than 2 items");
  const double mean=sum(a)/npts;
  double dev2=0.0;
  for (size_t i=npts;i--;) {
    const double d=a(i)-mean;
    dev2+=d*d;
  }
  dest.create(1,std::sqrt(dev2/npts));
}
  
void FunctionHead::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  if (a.empty())
    throw Failed("head applied to empty list");
  dest.create(1,a.front());
}

void FunctionHead2::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  const int n=round_int(getarg(size_t(1)).front());
  if (n<0)
    throw InvalidParameter("head: number of elements cannot be negative!");
  if (n>a.size())
    throw Failed("head: insufficient elements in list");
  dest=a.truncate(n);
}

void FunctionRev::operator()(LIST<double>& dest) const 
{
  const BaseList<double> a(getarg(size_t(0)));
  const size_t n=a.size();
  dest.create(n);
  for (size_t i=n;i--;)
    dest(i)=a(n-i-1);
}

void FunctionTail::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  if (a.empty())
    error_abort("tail applied to empty list");
  dest=a(range(1,a.size()-1));
}

void FunctionTail2::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  const int n=round_int(getarg(size_t(1)).front());
  if (n<0)
    throw InvalidParameter("tail: list index cannot be less than zero");
  if (n>a.size())
    throw Failed("tail: insufficient elements in list");
  if (n==a.size()) {
    dest.create(size_t(0));
    return;
  }
  dest=a(range(n,a.size()-1));
}

void FunctionReplace::operator()(LIST<double>& dest) const
{
  const BaseList<double> source(getarg(size_t(0)));
  const BaseList<double> inds(getarg(size_t(1)));
  const BaseList<double> values(getarg(size_t(2)));
  const size_t ninds=inds.size();
  const bool matching=(values.size()==ninds);
  if (!matching && (values.size()!=1))
    throw Mismatch("replace: number of values doesn't match number of indices",values.size(),ninds);
  const size_t max=source.size();
  const size_t baseindex=dest.size();
  dest.push_back(source); //!< copy original values and then modify in place
  if (ninds) {
    const double v=values.front();
    for (size_t i=ninds;i--;) {
      const size_t ind=baseindex+getindex(inds(i),max);
      dest(ind)= matching ? values(i) : v;
    }
  }
}
      
void FunctionLast::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  if (a.empty())
    error_abort("last applied to empty list");
  dest.create(size_t(1),a.back());  
}

void FunctionLast2::operator()(LIST<double>& dest) const
{
  const BaseList<double> a(getarg(size_t(0)));
  const int n=round_int(getarg(size_t(1)).front());
  if (n<0)
    throw InvalidParameter("last: number of elements cannot be negative!");
  if (n==0) {
    dest.create(size_t(0));
    return;
  }
  if (n>a.size())
    throw Failed("last: insufficient elements in list");
  dest=a(range(a.size()-n,a.size()-1));
}

void FunctionSize::operator()(LIST<double>& dest) const
{
  dest.create(1,double(getarg(size_t(0)).size()));
}

void FunctionRepeat::operator()(LIST<double>& dest) const
{
  const BaseList<double>& arg2(getarg(size_t(1)));
  if (arg2.size()!=1)
    throw InternalError("assertion failed: FunctionRepeat");
  const int ntimes=round_int(arg2.front());
  if (ntimes<0)
    throw InvalidParameter("repeat: count must be >=0");
  bool donefirst=false;
  for (size_t i=ntimes;i--;) {
    const BaseList<double> vals(getarg(size_t(0)));
    if (!donefirst) {
      dest.reserve(dest.size()+ntimes*vals.size());
      donefirst=true;
    }
    dest.push_back(vals); //!< don't cache getarg(0) in case non-const 
  }
}


void FunctionSortIndex::operator()(LIST<double>& dest) const
{
  const BaseList<double>& slist(getarg(size_t(0)));
  const size_t n=slist.size();
  dest.create(n);
  if (n==0)
    return;

  const LIST<size_t> indices(indexed_sort(slist));

  for (size_t i=n;i--;)
    dest(i)=indices(i)+1;
}

void FunctionExtract::operator()(LIST<double>& dest) const
{
  const BaseList<double>& slist(getarg(size_t(0)));
  const BaseList<double>& indices(getarg(size_t(1)));
  dest.create(indices.size());
  const size_t max=slist.size();
  for (size_t i=indices.size();i--;)
    dest(i)=slist(getindex(indices(i),max));
}

void FunctionSplitRow::operator()(LIST<double>& dest) const
{
  const BaseList<double>& slist(getarg(size_t(0)));
  const int num=round_int(getarg(size_t(1)).front());
  if (num<1)
    throw InvalidParameter("splitrow argument 2");
  const int which=round_int(getarg(size_t(2)).front());
  if ((which<1) || (which>num))
    throw InvalidParameter("splitrow argument 3");
  if (slist.size() % num) {
    char scratch[256];
    snprintf(scratch,sizeof(scratch),"can't split list of size %" LCM_PRI_SIZE_T_MODIFIER "u evenly between %i destination variables",slist.size(),num);
    throw Failed(scratch);
  }
  const size_t need=slist.size()/num;
  dest.reserve(dest.size()+need);
  const size_t start=(which-1)*need;
  const size_t end=which*need;
  for (size_t i=start;i<end;i++)
    dest.push_back(slist(i));
}

void FunctionSplitColumn::operator()(LIST<double>& dest) const
{
  const BaseList<double>& slist(getarg(size_t(0)));
  const int num=round_int(getarg(size_t(1)).front());
  if (num<1)
    throw InvalidParameter("splitcolumn argument 2");
  const int which=round_int(getarg(size_t(2)).front());
  if ((which<1) || (which>num))
    throw InvalidParameter("splitcolumn argument 3");
  const size_t need=(slist.size()+num-1)/num;
  dest.reserve(dest.size()+need);
  for (size_t i=which-1;i<slist.size();i+=num)
    dest.push_back(slist(i));
}
  
void declare_function(const char* name, size_t args, const ExpressionNamedBase& obj)
{
  static Function_Factory_t& Function_Factory(get_Function_Factory());
  const std::string* namep=new std::string(name); //!< create string so name is not volatile - Will give leak in valgrind
  const char* usename=namep->c_str();
  const function_spec testkey(usename,-1);
  const Function_Factory_t::const_iterator iter(Function_Factory.find(testkey));
  if (iter!=Function_Factory.end()) {
    const function_spec key(iter->first);
    parser_printcontext() << "Cannot declare new function (" << name << ") with fixed arguments to overwrite existing multi-variate function\n";
    error_abort();
  }
  Function_Factory[function_spec(usename,int(args))]=function_def_t(&obj);
}

void declare_function(const char* name, const ExpressionNamedBase& obj)
{
  static Function_Factory_t& Function_Factory(get_Function_Factory());
  //! N.B. doesn't bother to search for existing functions! (might be a "feature" that a new multi-variate function overwrites existing fixed argument functions ...)
  std::string* namep=new std::string(name); //!< create string so name is not volatile
  Function_Factory[function_spec(namep->c_str(),-1)]=function_def_t(&obj);
}
 
#if NMRSIM_USE_HASH
size_t FunctionHash::operator()(const function_spec& a) const 
{ 
  const size_t val= stringhasher_(a.first) ^ a.second;
  //  if (verbose & VER_PARSE) 
  //std::cout << "Hashed " << a.first << " with " << a.second << " arguments to " << std::hex << val << std::dec << '\n';
  return val;
}
#endif
 
struct proxy_ {
  proxy_() {
    const static FunctionSin sinobj(size_t(1));
    declare_function("sin",size_t(1),sinobj);
    const static FunctionCos cosobj(size_t(1));
    declare_function("cos",size_t(1),cosobj);
    const static FunctionTan tanobj(size_t(1));
    declare_function("tan",size_t(1),tanobj);
    //    const static FunctionFailed failedobj(size_t(0));
    //declare_function("failed",size_t(0),failedobj);
    const static FunctionBadParameter badparamobj(size_t(1));
    declare_function("badparameter",size_t(1),badparamobj);
    //declare_function("sinh",size_t(1),new FunctionSinh(size_t(1)));
    //declare_function("cosh",size_t(1),new FunctionCosh(size_t(1)));
    //declare_function("tanh",size_t(1),new FunctionTanh(size_t(1)));
    const static FunctionRandom rand1obj(size_t(1));
    const static FunctionRandom rand2obj(size_t(2));
    declare_function("random",size_t(1),rand1obj);
    declare_function("random",size_t(2),rand2obj);
    const static FunctionNoise noise1obj(size_t(1));
    const static FunctionNoise noise2obj(size_t(2));
    declare_function("noise",size_t(1),noise1obj);
    declare_function("noise",size_t(2),noise2obj);
    const static FunctionAsin asinobj(size_t(1));
    declare_function("asin",size_t(1),asinobj);
    const static FunctionAcos acosobj(size_t(1));
    declare_function("acos",size_t(1),acosobj);
    const static FunctionAtan atan1obj(size_t(1));
    const static FunctionAtan atan2obj(size_t(2));
    declare_function("atan",size_t(1),atan1obj);    
    declare_function("atan",size_t(2),atan2obj);
    const static FunctionIf ifobj(3U);
    declare_function("if",3U,ifobj);
    const static FunctionSwitch switchobj;
    declare_function("switch",switchobj);
   const static FunctionSwitchDefault switchdefaultobj;
    declare_function("switchdefault",switchdefaultobj);
    const static FunctionEcho echoobj(size_t(1));
    declare_function("echo",size_t(1),echoobj);
    const static FunctionLn lnobj(size_t(1));
    declare_function("ln",size_t(1),lnobj);
    //    declare_function("log10",size_t(1),new FunctionLog10(size_t(1)));
    const static FunctionExp expobj(size_t(1));
    declare_function("exp",size_t(1),expobj);
    const static FunctionSign signobj(size_t(1));
    declare_function("sign",size_t(1),signobj);
    const static FunctionFloor floor1obj(size_t(1));
    const static FunctionFloor floor2obj(size_t(2));
    declare_function("floor",size_t(1),floor1obj);
    declare_function("floor",size_t(2),floor2obj);
    const static FunctionCeil ceil1obj(size_t(1));
    const static FunctionCeil ceil2obj(size_t(2));
    declare_function("ceil",size_t(1),ceil1obj);
    declare_function("ceil",size_t(2),ceil2obj);
    const static FunctionRound round1obj(size_t(1));
    const static FunctionRound round2obj(size_t(2)); 
    declare_function("round",size_t(1),round1obj);
    declare_function("round",size_t(2),round2obj);
    const static FunctionLcm lcmobj(size_t(2));
    declare_function("lcm",size_t(2),lcmobj);
    const static FunctionLcm gcdobj(size_t(2));
    declare_function("gcd",size_t(2),gcdobj);
    const static FunctionSqrt sqrtobj(size_t(1));
    declare_function("sqrt",size_t(1),sqrtobj);
    const static FunctionAbs absobj(size_t(1));
    declare_function("abs",size_t(1),absobj);
    const static FunctionExtract extractobj(size_t(2));
    declare_function("extract",size_t(2),extractobj);
    const static FunctionRepeat repeatobj(size_t(2));
    declare_function("repeat",size_t(2),repeatobj);
    //    declare_function("norm",size_t(1),new FunctionNorm(size_t(1)));
    const static FunctionSize sizeobj(size_t(1));
    declare_function("size",size_t(1),sizeobj);
    const static FunctionGetReal getreal1obj(size_t(1));
    const static FunctionGetReal getreal2obj(size_t(2)); 
    declare_function("getreal",size_t(1),getreal1obj);
    declare_function("getreal",size_t(2),getreal2obj);
    const static FunctionGetImag getimag1obj(size_t(1));
    const static FunctionGetImag getimag2obj(size_t(2)); 
    declare_function("getimag",size_t(1),getimag1obj);
    declare_function("getimag",size_t(2),getimag2obj);
    const static FunctionHead head1obj(size_t(1));
    const static FunctionHead2 head2obj(size_t(2)); 
    declare_function("head",size_t(1),head1obj);
    declare_function("head",size_t(2),head2obj);
    const static FunctionTail tail1obj(size_t(1));
    const static FunctionTail2 tail2obj(size_t(2)); 
    declare_function("tail",size_t(1),tail1obj);
    declare_function("tail",size_t(2),tail2obj);
    const static FunctionLast last1obj(size_t(1));
    const static FunctionLast2 last2obj(size_t(2)); 
    declare_function("last",size_t(1),last1obj);
    declare_function("last",size_t(2),last2obj);
    const static FunctionRev revobj(size_t(1));
    declare_function("rev",size_t(1),revobj);
    const static FunctionReplace replaceobj(size_t(3));
    declare_function("replace",size_t(3),replaceobj);
    const static FunctionMutable mutableobj;
    declare_function("mutable",size_t(1),mutableobj);
    const static FunctionConst constobj;
    declare_function("const",size_t(1),constobj);  
    const static FunctionSortIndex sortindexobj(size_t(1));
    declare_function("sortindex",size_t(1),sortindexobj);

    const static FunctionSum sumobj(size_t(1));
    declare_function("sum",size_t(1),sumobj);

    const static FunctionError errorobj(size_t(1));
    declare_function("Error",size_t(1),errorobj);

    const static FunctionStdDev stddevobj(size_t(1));
    declare_function("stddev",size_t(1),stddevobj);

    const static FunctionWarn warnobj(size_t(1),false);
    declare_function("Warn",size_t(1),warnobj);
    const static FunctionWarn warnonceobj(size_t(1),true);
    declare_function("WarnOnce",size_t(1),warnonceobj); //!< 2 argument forms unlikely to be useful, so not declared
    const static FunctionIndicesof indicesobj0(size_t(0));
    const static FunctionIndicesof indicesobj1(size_t(1));
    declare_function("Indicesof",size_t(0),indicesobj0);
    declare_function("Indicesof",size_t(1),indicesobj1);
  }
};

static const proxy_ proxy;

