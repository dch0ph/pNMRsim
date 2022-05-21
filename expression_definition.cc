#include "expression_definition.hpp"
#include "rmatrix.h"

//bool allowwarnings();

static bool evaluatingisconstant=false;

template<typename T> class guardset {
public:
  guardset(T& varv, const T& valv)
    : var_(varv) {
    store_=var_;
    var_=valv;
  }
  ~guardset() { var_=store_; }
private:
  T& var_;
  T store_;
};
  
struct ExpressionUserFunction::stack {
  LIST< LIST< LIST<double> >* > stack_;

  size_t size() const { return stack_.size(); }
  const BaseList<double>& get(size_t) const;
  void pop() { stack_.pop_back(); }
  void push(LIST< LIST<double> >& listv) { stack_.push_back(&listv); }
  size_t nargs() const { return stack_.back()->size(); }
  void clear() { stack_.clear(); }
};

ExpressionUserFunction::ExpressionUserFunction(const function_spec& fspec, argdef_t islistv, const Expression*& exprpv)
  : ExpressionFunctionBase(fspec.first,fspec.second), 
    exprp_(exprpv),   //!< NMRSIM_NULL expression is valid for parsing
    islist_(islistv)
{}

const BaseList<double>& ExpressionUserFunction::getslot(size_t n) { return argstack_.get(n); }
void ExpressionUserFunction::clearstack() { argstack_.clear(); } //!< ensure argument stack is empty
void ExpressionGeneralFunction::clearstack() { argstack_.clear(); } //!< ensure argument stack is empty

ContextWarning<> recursion_warning("Highly recursive function call - infinite loop?",&NMRsim_once_warning);


void clearevalstacks()
{
  ExpressionUserFunction::clearstack();
  ExpressionGeneralFunction::clearstack();
  recursion_warning.reset();
}

ExpressionFunctionBase* ExpressionUserFunction::create(size_t nargsv) const
{
  if (nargsv!=nargs())
    throw InternalError("ExpressionUserFunction::create");
  return clone();
}

double ExpressionBase::get(const Expression& expr) const 
{
  LIST<double> scr;
  this->get(scr,expr);
  if (scr.size()!=1)
    throw Failed("Scalar result expected in this context");
  return scr.front();
}

void checkstack(int recurse_level)
{
  if (!nochecks)
    return;
  if (NMRSIM_RECURSE_ABORT && (recurse_level==NMRSIM_RECURSE_ABORT))
    recursion_warning.raiseas(BaseWarning::RaiseException);
  if (NMRSIM_RECURSE_WARN && (recurse_level==NMRSIM_RECURSE_WARN))
    recursion_warning.raise();
}

void Expression::print(std::ostream& ostr) const
{
  if (empty())
    ostr << "<undefined>";
  else
    back()->print(ostr,*this);
}

void ExpressionSlot::print(std::ostream& ostr, const Expression&) const
{
  ostr << '#';
  //  if (isstar_)
  //  ostr << '*';
  const size_t nargs=argnamestack.size();
  if (nargs==0)
    ostr << slot_;
  else {
    if (slot_ > nargs)
      ostr << slot_ << "<OUT OF RANGE>"; //!< something not right!
    else
      ostr << argnamestack(slot_-1);
  }    
}

size_t ExpressionSlot::findslot(const char* slotname)
{
  if (argnamestack.empty())
    error_abort("argument names used in definition, but none defined in declaration");

  for (size_t slotm1=argnamestack.size();slotm1--;) {
    if (strcmp(slotname,argnamestack(slotm1))==0)
      return slotm1+1;
  }
  parser_printcontext() << "failed to find argument name: " << slotname << '\n';
  error_abort();
  return -1; //!< to stop warning
}

ExpressionSlot::ExpressionSlot(size_t slotv)
  : slot_(slotv)//, isstar_(isstarv)
{
  if (slotv<1)
    throw InvalidParameter("#<n> function arguments are numbered from 1");
  if (!argnamestack.empty()) {
    if (slot_ > argnamestack.size())
      parser_printcontext() << "Argument " << slotv << " is out of range - function has only " << argnamestack.size() << " arguments)\n";
    if (*(argnamestack(slot_-1))=='\0')
      parser_printcontext() << "Argument " << slotv << " has been declared as unused\n";
  }
}

void ExpressionBase::setchildren(size_t start, size_t end)
{
  children=slice(start,end-start+1);
}

Expression::Expression(const ExpressionBase& obj, const BaseList< LIST<double> >& args)
{
  const size_t nargs=args.size();
  for (size_t n=0;n<nargs;n++)
    push_original(new ExpressionReference(const_cast< LIST<double>& >(args(n))));
  ExpressionBase* headp=obj.clone();
  push_original(headp);
  headp->setchildren(0U,nargs-1); //!< set pointers
}
  
void ExpressionConstant::print(std::ostream& ostr,const Expression&) const 
{
  if (value_.empty())
    ostr << "<empty>";
  else 
    dumparray(std::cout,value_);

  if (isinteger())
    ostr << " (integer)";
}

void ExpressionReference::print(std::ostream& ostr,const Expression&) const 
{
  if (value_.empty())
    ostr << "<empty>";
  else
    dumparray(std::cout,value_);

  if (isinteger())
    ostr << " (integer)";
}

double ExpressionDivide::operator()(double a,double b) const 
{
  if (b==0.0)
    throw Failed("divide failed: divide by zero");
  return a/b;
}

bool ExpressionShift::isconstant(int&, const Expression&) const
{
  static const bool isconst=proton_freq_isconstant();
  return isconst;
}

ExpressionShift::ExpressionShift(double valuev, size_t qualifier) 
  : grat_(getgrat(qualifier)), value_(1,valuev)
{
  if ((grat_==0.0) && proton_freq_isconstant())
    error_abort("proton_frequency must be set for ppm qualifiers to be used");
}

void ExpressionShift::print(std::ostream& ostr,const Expression&) const 
{
  if (value_.size()>1)
    ostr << value_;
  else
    ostr << value_.front();
  ostr << " ppm";
}

void ExpressionShift::get(LIST<double>& res, const Expression&) const
{
  if (!proton_freq)
    error_abort("proton_frequency not set");
  multiply(res,value_,grat_*proton_freq*1e-6);  
}

void ExpressionMinus::print(std::ostream& ostr, const Expression& expr) const
{
  ostr << '-';
  expr(children.front())->print(ostr,expr);
}

double ExpressionPower::operator()(double a, double b) const
{
  return pow(a,b);
}

double ExpressionPower::operator()(double a, int n) const
{
  switch (n) {
  case 0:
    return 1.0;
  case 1:
    return a;
  case 2:
    return a*a;
  case 3:
    return a*a*a;
  }
  return (*this)(a,(double)n);
}

double ExpressionModulus::operator()(double a,double b) const 
{
  if (b==0.0)
    throw Failed("modulus failed: divide by zero");
  return fmod(a,b);
    //a-b*floor(a/b);
}

void ExpressionRange::print(std::ostream& ostr, const Expression& expr) const
{
  expr(children.front())->print(ostr,expr);
  ostr << ':';
  expr(children(1U))->print(ostr,expr);
  if (children.size()==3) {
    ostr << ':';
    expr(children(2U))->print(ostr,expr);
  }
}

void ExpressionRange::get(LIST<double>& res, const Expression& expr) const 
{
  const double step=(children.size()==3) ? expr(children(1U))->get(expr) : 1.0;
  const double start=expr(children.front())->get(expr);
  const double end=expr(children.back())->get(expr);
  const double fsteps=(end-start)/step;
  const int steps=int(1.5+floor(fsteps+NMRSIM_ROUNDTOL));
  if (steps<1) {
    if (children.size()==3)
      throw Failed("<start>:<step>:<end> gives empty/invalid range");
    res.clear();
    return;
  }
  res.create(steps);
  for (size_t i=0;i<steps;i++)
    res(i)=start+i*step;
}

ExpressionUserFunction::stack ExpressionUserFunction::argstack_;
ExpressionGeneralFunction::stack ExpressionGeneralFunction::argstack_;

const BaseList<double>& ExpressionUserFunction::stack::get(size_t n) const
{
  if (stack_.empty())
    throw InternalError("stack::get");
  LIST< LIST<double> >& args(*(stack_.back()));
  if (n>args.size())
    throw BadIndex("ExpressionUserFunction::stack",n,args.size());
  return args(n-1);
}

void ExpressionSlot::get(LIST<double>& res, const Expression&) const
{
  res=ExpressionUserFunction::getslot(slot_);
}

bool ExpressionSlot::isconstant(int&, const Expression&) const
{
  return evaluatingisconstant; //!< Slot in general is non-const, but when evaluating isconstant, it is the constness of the argument themselves that is relevant
}

void ExpressionUnary::get(LIST<double>& res, const Expression& expr) const 
{
  expr(children.front())->get(res,expr);
  for (size_t i=res.size();i--;)
    res(i)=(*this)(res(i));
}

void ExpressionBinary::get(LIST<double>& res, const Expression& expr) const 
{
  expr(children.front())->get(res,expr);
  const ExpressionBase& arg2(*(expr(children(1U))));
  LIST<double> scr;
  arg2.get(scr,expr);
  const bool isint(arg2.isinteger());
  if (res.size()==scr.size()) {
    for (size_t i=res.size();i--;)
      res(i)=docalc(res(i),scr(i),isint);
    return;
  }  
  static char errmess[256];

  switch (scr.size()) {
  case 0:
    if (res.size()==1) {
      res.clear();
      return;
    }
    snprintf(errmess,sizeof(errmess),"Empty (right) argument to operator %c",op_);
    throw Failed(errmess);
  case 1: {
    const double val=scr.front();
    for (size_t i=res.size();i--;)
      res(i)=docalc(res(i),val,isint);
  }
    return;
  }

  switch (res.size()) {
  case 0:
    if (scr.size()==1)
      return;
    snprintf(errmess,sizeof(errmess),"Empty (left) argument to operator %c",op_);
    throw Failed(errmess);
  case 1: {
    const double val=res.front();
    res.create(scr.size());
    for (size_t i=scr.size();i--;)
      res(i)=docalc(val,scr(i),isint);
  }
    return;
  }
  snprintf(errmess,sizeof(errmess),"Arguments to operator %c have incompatible lengths",op_);
  throw Mismatch(errmess,res.size(),scr.size());
}

void ExpressionBinary::print(std::ostream& ostr, const Expression& expr) const
{
  ostr << '(';
  expr(children.front())->print(ostr,expr);
  ostr << op_;
  expr(children(1U))->print(ostr,expr);
  ostr << ')';
}

bool ExpressionBase::isconstant(int& root, const Expression& expr) const
{
  const bool isverb=(verbose & VER_PARSE) && (verbose_level>1);
  if (isverb) {
    std::cout << "Evaluating const status of expression " << this << ": ";
    print(std::cout, expr);
    std::cout << ": ";
  }
  bool isconst=true;
  for (size_t i=children.size();i--;) {      
    const size_t off=children(i);
    if (!(expr(off)->isconstant(root,expr))) {
      isconst=false;
      break;
    }
    if ((root<0) || (off<root))
      root=off;
  }
  if (isverb)
    std::cout << (isconst ? "Yes\n" : "No\n");
  return isconst;
}

int ExpressionBase::uses(const Expression& expr) const
{
  int res=uses();
  for (size_t i=children.size();i--;) {
    const size_t off=children(i);
    res |= expr(off)->uses(expr);
  }
  return res;
}

void ExpressionSimpleFunction::get(LIST<double>& res, const Expression& expr) const
{
  const size_t N=children.size();
  if (isfunction()) {
    LIST<double> tmpvals(N);
    LIST< LIST<double> > tmpargs(N);
    size_t maxrep=1;
    for (size_t i=0;i<N;i++) {
      LIST<double>& curarg(tmpargs(i));
      expr(children(i))->get(curarg,expr);  
      const size_t n=curarg.size();
      if ((n!=1) && (n!=maxrep)) {
	if (maxrep==1)
	  maxrep=n;
	else
	  throw Mismatch("ExpressionFunction: list arguments have incompatible lengths",maxrep,n);
      }
    }    
    res.create(0U);
    for (size_t i=0;i<maxrep;i++) {
      for (size_t M=0;M<N;M++) {
	const BaseList<double>& curarg(tmpargs(M));
	tmpvals(M)=(curarg.size()>1) ? curarg(i) : curarg.front();
      }
      (*this)(res,tmpvals);
    }
  }
  else {
    if (children.size()==1)
      expr(children.front())->get(res,expr);
    else {
      res.clear();
      //      if ((&scrlist)==(&res))
      //	throw ArgumentClash("ExpressionSimpleFunction");
      LIST<double> scr;
      for (size_t i=0;i<children.size();i++) {
	expr(children(i))->get(scr,expr);
	res.push_back(scr);
      }
    }
  }  
}

void ExpressionFunctionBase::print(std::ostream& ostr, const Expression& expr) const
{
  if (isfunction())
    ostr << name_ << '(';
  else
    ostr << '[';
  for (size_t i=0;i<children.size();i++) {
    if (i)
      ostr << ',';
    expr(children(i))->print(ostr,expr);
  }
  ostr << (isfunction() ? ')' : ']');
}

void ExpressionStringFunction::rawprint(std::ostream& ostr) const
{
  ostr << name_ << '(';
  for (size_t i=0;i<args_.size();i++) {
    if (i)
      ostr << ',';
    ostr << '\'' << args_(i) << '\'';
  }
  ostr << ')';
}

ContextWarning<> ExpressionVariableFunction::appliedtoexpr_warning("variable function applied to expression rather than value set",&NMRsim_repeat_warning);

ExpressionGeneralFunction::ExpressionGeneralFunction(const char* name, const BaseList<size_t>& offs)
  : ExpressionFunctionBase(name,offs) {}

ExpressionGeneralFunction::ExpressionGeneralFunction(const char* name, argdef_t islistv, const BaseList<size_t>& offs)
  : ExpressionFunctionBase(name,offs),
    islist_(islistv) {}

ExpressionGeneralFunction::ExpressionGeneralFunction(const char* name, argdef_t islistv, size_t nargs)
  : ExpressionFunctionBase(name,nargs),
    islist_(islistv) {}

ExpressionGeneralFunction::ExpressionGeneralFunction(const char* name)
  : ExpressionFunctionBase(name,-1) {}

bool ExpressionGeneralFunction::islist(unsigned short M) const 
{
  if (ismultivariate())
    return true;
  return (islist_ >> M) & 1;
}

struct accargs {
  accargs(const char* namev =NMRSIM_NULL)
    : maxrep(1U),name_(namev) {}

  void add(size_t);
  size_t operator()() const { return maxrep; }
  size_t maxrep;
  const char* name_;
};

ContextWarning<> zerolengthscalar_warning("zero length scalar argument - argument not specified as list(L) - ",&NMRsim_repeat_warning);

void accargs::add(size_t n)
{
  if (n==0)
    zerolengthscalar_warning.raise(name_);
  if (maxrep!=n) {
    if (maxrep==1)
      maxrep=n;
    else {
      if (n!=1)
	throw Mismatch("scalar arguments to function have incompatible lengths",maxrep,n);
    }
  }
}

void ExpressionUserFunction::setexpression(const Expression& expr)
{
  exprp_=&expr; //!< this will set expr in all copies of function
  //  setisconst(expr.isconstant()); //!< this will not set const in recursively used copies, but not a problem since only the "master" version's const status counts
}

bool ExpressionUserFunction::isconstant(int& root, const Expression& expr) const
{
  if (!(ExpressionBase::isconstant(root,expr)))
    return false;
  if (!exprp_) //!< function must be recursively defined, so return false
    return false; 
    //   throw Failed("ExpressionUserFunction::isconstant: called before expression set");
  guardset<bool> guard(evaluatingisconstant,true); //!< reinstated Feb 13 as recursively defined functions in fulltrans.inc don't play well otherwise
  return exprp_->isconstant();
}

void ExpressionUserFunction::get(LIST<double>& res, const Expression& expr) const
{
  if (!exprp_)
    throw Failed("ExpressionUserFunction::get: called before expression set");
  accargs acc(name_);
  const size_t N=children.size();
  LIST< LIST<double> > scalarargs(N);
  LIST< LIST<double> > args(N);
  for (size_t i=0;i<N;i++) {
    if (islist(i))
      expr(children(i))->get(args(i),expr);  
    else {
      expr(children(i))->get(scalarargs(i),expr);  
      acc.add(scalarargs(i).size());
      args(i).create(1,-1e6);
    }
  }  
  argstack_.push(args); //!< add arguments to stack
  checkstack(argstack_.size());
  res.create(0U);
  LIST<double> tmp;
  const bool doverb=(verbose & VER_PARSE) && (verbose_level>1);
  if (doverb) {
    std::cout << "Entering " << name_ << " (" << N << " argument(s)):";
    if (acc()!=1)
      std::cout << " repeat: " << acc();
    std::cout << '\n';
  }
  for (size_t which=0;which<acc();which++) {
    for (size_t M=0;M<N;M++) {
      if (!islist(M)) {
	const BaseList<double> curarg(scalarargs(M));
	const size_t actwhich=(curarg.size()>1) ? which : 0U;
	args(M).front()=curarg(actwhich);
      }
    }
    exprp_->back()->get(tmp,*exprp_); //!< don't call expr_->get() directly (resets evaluation)
    res.push_back(tmp);
    if (doverb) {
      std::cout << name_ << '(';
      for (size_t M=0;M<N;M++) {
	if (M)
	  std::cout << ',';
	dumparray(std::cout,args(M));
      }
      std::cout << "): ";
      dumparray(std::cout,tmp) << '\n';
    }     
  }
  argstack_.pop(); //!< pop argument stack (don't bother with exception safety since evaluation stack will be aborted and then reset)
}

class ExpressionGeneralFunction::Evaluator {
public:
  Evaluator(const ExpressionGeneralFunction& objv, const Expression& exprv);
  ~Evaluator() { ExpressionGeneralFunction::argstack_.pop_back(); }

  const BaseList<double>& getarg(size_t n) const;
  size_t nargs() const { return tmpargs.size(); } 
  void operator()(LIST<double>&);

private:
  const ExpressionGeneralFunction& obj_;
  const Expression& expr_;
  mutable LIST< LIST<double> > tmpargs;
  mutable LIST< BaseList<double> > actargs;
  mutable ScratchList<bool> done;
  accargs acc;
};

ExpressionGeneralFunction::Evaluator::Evaluator(const ExpressionGeneralFunction& objv, const Expression& exprv)
  : obj_(objv), expr_(exprv),
    tmpargs(objv.nargs()),
    actargs(objv.nargs()),    
    done(objv.nargs(),false),
    acc(objv.name_)
{
  ExpressionGeneralFunction::argstack_.push_back( (Evaluator*)NMRSIM_NULL); //new stack frame, but flag not complete
  checkstack(argstack_.size());
  for (size_t i=0;i<nargs();i++) {
    if (!(obj_.islist(i))) {
      LIST<double>& curarg(tmpargs(i));
      expr_(obj_.children(i))->get(curarg,expr_);
      acc.add(curarg.size());
      done(i)=true;
    }
  }
  ExpressionGeneralFunction::argstack_.back()=this;
}

void ExpressionGeneralFunction::Evaluator::operator()(LIST<double>& res)
{
  const bool doverb=(verbose & VER_PARSE) && (verbose_level>1);
  res.create(0U);
  const size_t N=nargs();
  for (size_t which=0;which<acc();which++) {
    for (size_t M=0;M<N;M++) {
      if (!obj_.islist(M)) {
	const LIST<double>& curarg(tmpargs(M));
	BaseList<double>& actarg(actargs(M));
	if (curarg.size()>1)
	  actarg.create(1,curarg.vector()+which);
	else
	  actarg.create(1,curarg.vector());
      }
    }
    obj_(res);
    if (doverb) {
      std::cout << obj_.name_ << '(';
      for (size_t M=0;M<N;M++) {
	if (M)
	  std::cout << ',';
	dumparray(std::cout,actargs(M));
      }
      std::cout << ") -> ";
      dumparray(std::cout,res) << '\n';
    }     
  }
}

const BaseList<double>& ExpressionGeneralFunction::getarg(size_t n) const
{
  if (argstack_.empty())
    throw InternalError("Empty GeneralFunction eval stack");
  Evaluator* evalp=argstack_.back();
  if (!evalp)
    throw InternalError("GeneralFunction::getarg");
  return evalp->getarg(n);
}

const BaseList<double>& ExpressionGeneralFunction::Evaluator::getarg(size_t n) const  
{
  if (n>=nargs())
    throw BadIndex("Out of range argument",n,nargs());
  BaseList<double>& curarg(actargs(n));
  if (!done(n)) { //!< lazy evaluation of list arguments
    if (!obj_.islist(n))
      throw InternalError("ExpressionGeneralFunction::Evaluator");
    const ExpressionBase& expr(*(expr_(obj_.children(n))));
    expr.get(tmpargs(n),expr_);  
    curarg.create(tmpargs(n));
    if (expr.isconstant(expr_)) //!< if object is non-const don't set done flag
      done(n)=true;
  }
  return curarg;
}

void ExpressionGeneralFunction::get(LIST<double>& res, const Expression& expr) const 
{
  Evaluator evalobj(*this,expr);
  evalobj(res);
}

ExpressionVariable::ExpressionVariable(RealVariable& varv)
  : var_(varv) {}

void ExpressionVariable::get(LIST<double>& res, const Expression&) const
{
  res=var_.get_list();
}

void ExpressionVariable::print(std::ostream& ostr,const Expression&) const
{ 
  ostr << '$' << var_.name();
}

// void ExpressionExpression::get(LIST<double>& res, const Expression&) const
// {
//   expr_.get(res);
// }

// void ExpressionExpression::print(std::ostream& ostr,const Expression&) const
// { 
//   expr_.print(ostr);
// }

argdef_t ExpressionGeneralFunction::makeislist(const BaseList<bool>& islistv)
{
  argdef_t islist=0;
  for (size_t M=islistv.size();M--;) {
    islist<<=1U;
    if (islistv(M))
      islist|=1U;
  }
  return islist;
}
  
void error_abort(const char* str, int lerrno)
{
  if (isinteractive)
    throw Failed(str);
  parser_printcontext(std::cerr) << str << std::endl;
  error_abort(lerrno);
}

void error_abort(int lerrno)
{  
  std::cerr.flush();
  std::cout.flush(); //make sure output is flushed
  if (isinteractive)
    throw Failed("Operation failed");
  exit(lerrno);
}

ContextWarning<> include_nosubstitution_warning("`` expression with no variable substitution",&NMRsim_repeat_warning);

void ExpressionInclude::validate()
{
  if (strchr(filename_.c_str(),'$')) {
    if (!isdynamic_)
      error_abort("Can't use $ variables in static includes");
  }
  else {
    if (isdynamic_)
      include_nosubstitution_warning.raise();
  }
}

void ExpressionInclude::get(LIST<double>& res, const Expression&) const
{
  if (!isdynamic_) {
    FILE* fp=pathopen(filename_.c_str(),"r");
    int fail=(fp==NMRSIM_NULL);
    if (fail)
      parser_printcontext(std::cerr) << "failed to read: " << filename_ << '\n';
    else {
      try {	
	read_vector(res,fp);
      }
      catch (std::exception& exc) {
	fail=true;
	parser_printcontext(std::cerr) << "failed to parse " << filename_ << ": " << exc.what() << '\n';
      }
      fclose(fp);
    }
    if (fail)
      throw Failed("File open failed");
    return;
  }
  substitute_string(buffer,sizeof(buffer),filename_.c_str(),SUB_ESCAPE); //substitute variables
  if ((verbose && VER_GEN) && (verbose_level>1))
    std::cout << "Executing dynamic include: " << buffer << '\n';
  FILE* pfp=popen(buffer,"r");
  if (pfp==NMRSIM_NULL) {
    parser_printcontext(std::cerr) << "`` expression failed: " << buffer << '\n';
    error_abort();
  }
  bool failed=false;
  try {
    read_vector_ascii(res,pfp);    
  }
  catch (...) {
    failed=true;
  }
  int fcode=pclose(pfp);
  if (fcode<0) 
    perror("pclose"); //report system error
  if (failed || fcode) {
    parser_printcontext(std::cerr) << "`` expression (" << buffer << ") failed or didn't return simple number\n";
    error_abort();
  }
}
