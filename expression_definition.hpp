#ifndef expression_definition_hpp_
#define expression_definition_hpp_

/*! \file
 \brief  Definition of expressions
*/

#include "NMRsim_common.h"
#include "ScratchList.h"

using namespace libcmatrix;

//! node for $ variable reference
class ExpressionVariable : public ExpressionBase {
public:
  explicit ExpressionVariable(RealVariable&);
  ExpressionBase* clone() const { return new ExpressionVariable(*this); }
  void get(LIST<double>&, const Expression&) const;
  void print(std::ostream&, const Expression&) const;
  bool isconstant(int&, const Expression&) const { return var_.isconstant(); }
  int uses() const { return var_.uses(); }
  struct Add;
  struct ParseFailed;
private:
  RealVariable& var_; //!< reference to $ variable
};

//! reference to n-th argument (user-defined function)
class ExpressionSlot : public ExpressionBase {
public:
  explicit ExpressionSlot(size_t);
  ExpressionBase* clone() const { return new ExpressionSlot(*this); }
  void get(LIST<double>&, const Expression&) const;
  void print(std::ostream&, const Expression&) const;
  bool isconstant(int&, const Expression&) const;
  static size_t findslot(const char*);
private:
  size_t slot_; //!< argument number
};

//! node referencing another expression  
// class ExpressionExpression : public ExpressionBase {
// public:
//   explicit ExpressionExpression(const Expression& exprv) : ExpressionBase(exprv.isconstant()), expr_(exprv) {}
//   ExpressionBase* clone() const { return new ExpressionExpression(*this); }
//   void get(LIST<double>&, const Expression&) const;
//   void print(std::ostream&, const Expression&) const;
// private:
//   const Expression& expr_; //!< reference to another expression
// };

//! node for file include
class ExpressionInclude : public ExpressionBase {
public:
  ExpressionInclude(const char* filenamev, bool isdynamic)
    : filename_(filenamev), isdynamic_(isdynamic)
  { validate(); }

  ExpressionBase* clone() const { return new ExpressionInclude(*this); }

  void get(LIST<double>& res, const Expression&) const;
  bool isconstant(int&, const Expression&) const { return !isdynamic_; }

  template<bool> struct Add;

  void print(std::ostream& ostr,const Expression&) const {
    const char paren(isdynamic_ ? '`' : '"');
    ostr << paren << filename_ << paren;
  } 
private:
  std::string filename_; //!< name of file
  bool isdynamic_;
  mutable char buffer[256]; //!< temporary buffer
  void validate();
};

//! node for start:step:end
struct ExpressionRange : public ExpressionBase {
public:
  ExpressionRange(size_t st, size_t end)
    : ExpressionBase(st,end) {} //!< constructor without step

  ExpressionRange(size_t st, size_t step, size_t end)
    : ExpressionBase(ExplicitList<3,size_t>(st,step,end)) {} //!< constructor with step

  ExpressionBase* clone() const { return new ExpressionRange(*this); }

  void print(std::ostream&, const Expression&) const;
  void get(LIST<double>&, const Expression&) const;

  template<bool HasMiddle> struct Add;
};

//! node for constant
class ExpressionConstant : public ExpressionBase {
public:
  explicit ExpressionConstant(double valuev, bool isintv =false) 
    : value_(1,valuev), isint_(isintv) {} //!< scalar constant
  
  explicit ExpressionConstant(const BaseList<double>& valuev, bool isintv =false)
    : value_(valuev), isint_(isintv) {} //!< list constant

  ExpressionConstant(const ExpressionBase& item, const Expression& expr, bool isintv =false)
    : isint_(isintv) 
  { item.get(value_,expr); } //!< create constant from constant sub-expression

  ExpressionBase* clone() const { return new ExpressionConstant(*this); }
  void get(LIST<double>& res, const Expression&) const { res=value_; }
  void print(std::ostream&, const Expression&) const;
  bool isinteger() const { return isint_; } //!< \c true if constant is integral
  const BaseList<double>& operator()() const { return value_; }

  struct Add;
  template<int> struct AddPhase;

  friend class VariableBuilder;
private:
  LIST<double> value_; //!< value
  bool isint_; //!< \c true if constant is integral
};

//! node for reference
class ExpressionReference : public ExpressionBase {
public:
  explicit ExpressionReference(LIST<double>& valuev, bool isintv =false)
    : value_(valuev), isint_(isintv) {}

  ExpressionBase* clone() const { return new ExpressionReference(*this); }
  void get(LIST<double>& res, const Expression&) const { res=value_; }
  void print(std::ostream&, const Expression&) const;
  bool isinteger() const { return isint_; } //!< \c true if constant is integral
  const BaseList<double>& operator()() const { return value_; }
private:
  LIST<double>& value_; //!< value
  bool isint_; //!< \c true if constant is integral
};

//! node for shift (as ppm value)
class ExpressionShift : public ExpressionBase {
public:
  explicit ExpressionShift(double valuev, size_t =0); 
  
  ExpressionBase* clone() const { return new ExpressionShift(*this); }
  void get(LIST<double>&, const Expression&) const;
  void print(std::ostream&, const Expression&) const;
  bool isconstant(int&, const Expression&) const; //!< override isconstant
  struct Add;
private:
  double grat_; //!< gammaX/gamma1H ratio
  LIST<double> value_; //!< value(s)
};

//! (abstract) base class for unary functions
struct ExpressionUnary : public ExpressionBase {
  explicit ExpressionUnary(size_t sourcev)
    : ExpressionBase(sourcev) {}

  virtual double operator()(double) const =0; //!< apply function
  void get(LIST<double>&, const Expression&) const;

  template<typename T> struct Add;
};

//! unary minus node
struct ExpressionMinus : public ExpressionUnary {
  explicit ExpressionMinus(size_t sourcev) : ExpressionUnary(sourcev) {}
  double operator()(double v) const { return -v; }
  void print(std::ostream&, const Expression&) const;
  ExpressionBase* clone() const { return new ExpressionMinus(*this); }
};

//! (abstract) base class for binary functions
class ExpressionBinary : public ExpressionBase {
public:
  ExpressionBinary(size_t leftv, size_t rightv, char opv)
    : ExpressionBase(leftv,rightv), op_(opv) {}

  virtual double operator()(double,double) const =0; //!< apply function
  virtual double operator()(double a,int b) const { return (*this)(a,(double)b); } //!< fallback for integer specialisation
  virtual void print(std::ostream&, const Expression&) const;
  void get(LIST<double>&, const Expression&) const;
  template<typename T> struct Add;
private:
  char op_; //!< symbol for function
  double docalc(double a, double b, bool isint) const {
    return isint ? (*this)(a,round_int(b)) : (*this)(a,b);
  }
};

//! multiply expression node
struct ExpressionMultiply : public ExpressionBinary {
  ExpressionMultiply(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'*') {}
  double operator()(double a,double b) const { return a*b; }
  ExpressionBase* clone() const { return new ExpressionMultiply(*this); }
};  

//! divide expression node
struct ExpressionDivide : public ExpressionBinary {
  ExpressionDivide(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'/') {}
  double operator()(double,double) const;
  ExpressionBase* clone() const { return new ExpressionDivide(*this); }
};  

//! % expression node
struct ExpressionModulus : public ExpressionBinary {
  ExpressionModulus(size_t leftv,size_t rightv)
    : ExpressionBinary(leftv,rightv,'%') {}
  double operator()(double,double) const;
  ExpressionBase* clone() const { return new ExpressionModulus(*this); }
};  

//! ^ expression node
struct ExpressionPower : public ExpressionBinary {
  ExpressionPower(size_t leftv,size_t rightv)
    : ExpressionBinary(leftv,rightv,'^') {}
  double operator()(double,double) const;
  double operator()(double,int) const;
  ExpressionBase* clone() const { return new ExpressionPower(*this); }
};  

//! add expression node
struct ExpressionAdd : public ExpressionBinary {
  ExpressionAdd(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'+') {}
  double operator()(double a,double b) const { return a+b; }
  ExpressionBase* clone() const { return new ExpressionAdd(*this); }
};  

struct ExpressionGT : public ExpressionBinary {
  ExpressionGT(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'>') {}
  double operator()(double a,double b) const { return (a > b) ? 1.0 : 0.0; }
  ExpressionBase* clone() const { return new ExpressionGT(*this); }
};  

struct ExpressionLT : public ExpressionBinary {
  ExpressionLT(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'<') {}
  double operator()(double a,double b) const { return (a < b) ? 1.0 : 0.0; }
  ExpressionBase* clone() const { return new ExpressionLT(*this); }
};  

//! subtract expression node
struct ExpressionSubtract : public ExpressionBinary {
  ExpressionSubtract(size_t leftv,size_t rightv) : ExpressionBinary(leftv,rightv,'-') {}
  double operator()(double a,double b) const { return a-b; }
  ExpressionBase* clone() const { return new ExpressionSubtract(*this); }
};  

//! base class for nodes taking string arguments
class ExpressionStringFunction : public ExpressionNamedBase
{
public:
  // ExpressionStringFunction(const char* name, BaseList<std::string>& args)
  //   : name_(name), nargs_(args.size()) { create(args); }

  explicit ExpressionStringFunction(const char* name, int nargsv =1)
    : ExpressionNamedBase(name), nargs_(nargsv) {}
  
  bool ismultivariate() const { return (nargs_<0); } //!< \c true if multivariate
  //  bool isfunction() const { return (name_!=NMRSIM_NULL); } //!< returns \c true if function (rather than list)
  void print(std::ostream& ostr, const Expression&) const { rawprint(ostr); }
  virtual void rawprint(std::ostream&) const;
  size_t nargs() const { return (nargs_<0) ? args_.size() : nargs_; } //!< return number of arguments
  virtual void create(BaseList<std::string>&);
  bool isconstant(int&, const Expression&) const { return rawisconstant(); } //!< \c true if value is constant
  virtual bool rawisconstant() const { return true; } //!< simplified isconstant interface
  void get(LIST<double>& dest, const Expression&) const { (*this)(dest); }
  virtual void operator()(LIST<double>&) const =0;

  struct AddString;

protected:
  LIST<std::string> args_; //!< arguments
  int nargs_; //!< -1 for multi-variate
};

class UserVariable;

//!< base class for functions operating on variable names
class ExpressionVariableFunction : public ExpressionStringFunction
{
public:
  explicit ExpressionVariableFunction(const char* name, int nargsv =0)
    : ExpressionStringFunction(name,nargsv), varp(NMRSIM_NULL) {}

  void create(BaseList<std::string>&);
  void rawprint(std::ostream&) const;
  static ContextWarning<> appliedtoexpr_warning; //!< can't apply variable function to expression
  //  static ContextWarning<> nonarray_warning; //!< Values/Errorsof type function applied to variable that is not an array
  static ContextWarning<> outofmainloop_warning; //!< evaluated outside of main loop
  static ContextWarning<> outofmainloop_expression_warning; 
protected:
  UserVariable* varp; //!< pointer to user-defined variable
  void verify_initialised() const;
};

void declare_function(const char*, size_t n, const ExpressionNamedBase&); //!< declare function taking \a n arguments
void declare_function(const char*, const ExpressionNamedBase&); //!< declare multi-variate function

//! base class for function node
class ExpressionFunctionBase : public ExpressionNamedBase
{
public:
  ExpressionFunctionBase(const char* name, const BaseList<size_t>& offs)
    : ExpressionNamedBase(name,offs), nargs_(offs.size()) {}

  explicit ExpressionFunctionBase(const char* name, int nargsv =0)
    : ExpressionNamedBase(name), nargs_(nargsv) {}
  
  bool ismultivariate() const { return (nargs_<0); }
  bool isfunction() const { return (name_!=NMRSIM_NULL); } //!< returns \c true if function (rather than list)
  void print(std::ostream&, const Expression&) const;
  size_t nargs() const { return (nargs_<0) ? children.size() : nargs_; }

  struct Start;
  struct End;
  static ExpressionFunctionBase* create_function(const char*, size_t); //!< parse and create new object

protected:
  int nargs_; //!< -1 for multi-variate
  double checknan(double) const; //!< filter throwing \c InvalidParameter if NaN  
  double checkhuge(double) const; //!< filter throwing \c Failed if \c HUGE_VAL
  double checkneghuge(double) const; //!< filter throwing \c Failed if \c -HUGE_VAL
};

//! base class for simple (non-list) functions
class ExpressionSimpleFunction : public ExpressionFunctionBase
{
public:
  explicit ExpressionSimpleFunction()
    : ExpressionFunctionBase(NMRSIM_NULL) {}

  ExpressionSimpleFunction(const char* name, const BaseList<size_t>& offs)
    : ExpressionFunctionBase(name,offs),
      nargs_(offs.size()) {}

  explicit ExpressionSimpleFunction(const char* name, size_t nargsv =0)
    : ExpressionFunctionBase(name,nargsv),
      nargs_(nargsv) {}

  virtual ExpressionBase* clone() const { return new ExpressionSimpleFunction(*this); }
  virtual void operator()(LIST<double>&, const BaseList<double>&) const
  { throw InternalError("ExpressionSimpleFunction"); }

  void get(LIST<double>&, const Expression&) const;

private:
  size_t nargs_;
};

typedef unsigned long argdef_t; //!< signature for pattern of arguments

//! (abstract) base class for general functions
class ExpressionGeneralFunction : public ExpressionFunctionBase {
public:
  ExpressionGeneralFunction(const char* name, argdef_t islistv, const BaseList<size_t>& offs);
  ExpressionGeneralFunction(const char* name, argdef_t islistv, size_t nargs);
  ExpressionGeneralFunction(const char* name); //!< declare multi-variate
  ExpressionGeneralFunction(const char* name, const BaseList<size_t>& offs);

  virtual ~ExpressionGeneralFunction() {}
  virtual void operator()(LIST<double>&) const =0;
  void get(LIST<double>&, const Expression&) const;

  static argdef_t makeislist(const BaseList<bool>&); //!< create signature for pattern of arguments
  static void clearstack(); //!< \internal
  bool islist(unsigned short M) const; //!< \c true if argument \a M is a list (rather than scalar)

protected:
  const BaseList<double>& getarg(size_t) const;
  class Evaluator;
  typedef LIST<Evaluator*> stack;
private:  
  argdef_t islist_; //!< argument pattern signature
  static stack argstack_;
};

#define EXPR_DECLARE_FUNCTION(CLASSNAME,NAME)		\
struct CLASSNAME : public ExpressionSimpleFunction {\
  explicit CLASSNAME(size_t nargs) : ExpressionSimpleFunction(NAME,nargs) {} \
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionSimpleFunction(NAME,childv) {} \
ExpressionBase* clone() const { return new CLASSNAME(*this); }\
  void operator()(LIST<double>&, const BaseList<double>&) const;\
};

#define EXPR_DECLARE_FUNCTION_USEARRAY(CLASSNAME,NAME)		\
struct CLASSNAME : public ExpressionSimpleFunction {\
  explicit CLASSNAME(size_t nargs) : ExpressionSimpleFunction(NAME,nargs) {} \
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionSimpleFunction(NAME,childv) {} \
  ExpressionBase* clone() const { return new CLASSNAME(*this); }	\
  bool isconstant(int&, const Expression&) const { return false; }\
  int uses() const { return A_ARRAY; }\
  void operator()(LIST<double>&, const BaseList<double>&) const;\
};

#define EXPR_DECLARE_GENERAL_FUNCTION(CLASSNAME,NAME,ARGL)	\
struct CLASSNAME : public ExpressionGeneralFunction {\
  explicit CLASSNAME(size_t nargs) : ExpressionGeneralFunction(NAME,ARGL,nargs) {} \
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionGeneralFunction(NAME,ARGL,childv) {} \
  ExpressionBase* clone() const { return new CLASSNAME(*this); }\
  void operator()(LIST<double>&) const;\
};

#define EXPR_DECLARE_GENERAL_FUNCTION_MODIFIER(CLASSNAME,NAME,ARGL,CONST)	\
struct CLASSNAME : public ExpressionGeneralFunction {\
  explicit CLASSNAME(size_t nargs) : ExpressionGeneralFunction(NAME,ARGL,nargs) {} \
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionGeneralFunction(NAME,ARGL,childv) {} \
  ExpressionBase* clone() const { return new CLASSNAME(*this); }\
  bool isconstant(int&, const Expression&) const { return CONST; }\
  void operator()(LIST<double>&) const;\
};

extern const argdef_t Ldef;

#define EXPR_DECLARE_MODIFIER_FUNCTION(CLASSNAME,NAME,CONST)	\
struct CLASSNAME : public ExpressionGeneralFunction {\
  explicit CLASSNAME() : ExpressionGeneralFunction(NAME,Ldef,1U) {}	\
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionGeneralFunction(NAME,Ldef,childv) {} \
  ExpressionBase* clone() const { return new CLASSNAME(*this); }\
  bool isconstant(int&, const Expression&) const { return CONST; }\
  void operator()(LIST<double>& dest) const { dest=getarg(0U); }\
}

#define EXPR_DECLARE_MULTIVARIATE_FUNCTION(CLASSNAME,NAME)	\
struct CLASSNAME : public ExpressionGeneralFunction {\
  explicit CLASSNAME() : ExpressionGeneralFunction(NAME) {}		\
  explicit CLASSNAME(const BaseList<size_t>& childv) : ExpressionGeneralFunction(NAME,childv) {} \
  ExpressionBase* clone() const { return new CLASSNAME(*this); }	\
  void operator()(LIST<double>&) const;\
};

#define EXPR_DECLARE_VARIABLE_FUNCTION(CLASSNAME,NAME)\
  struct CLASSNAME : public ExpressionVariableFunction {\
    explicit CLASSNAME() : ExpressionVariableFunction(NAME) {}\
    ExpressionBase* clone() const { return new CLASSNAME(*this); }\
    bool rawisconstant() const;\
    void operator()(LIST<double>&) const;\
  };

#define EXPR_DECLARE_STRING_FUNCTION(CLASSNAME,NAME)\
  struct CLASSNAME : public ExpressionStringFunction {\
    explicit CLASSNAME(size_t nargsv =1U) : ExpressionStringFunction(NAME, nargsv) {} \
    ExpressionBase* clone() const { return new CLASSNAME(*this); }\
    bool rawisconstant() const;\
    void operator()(LIST<double>&) const;\
  };

extern const argdef_t LLdef;
EXPR_DECLARE_GENERAL_FUNCTION(FunctionExtract,"extract",LLdef);
extern const argdef_t LSSdef;
EXPR_DECLARE_GENERAL_FUNCTION(FunctionSplitRow,"splitrow",LSSdef);
EXPR_DECLARE_GENERAL_FUNCTION(FunctionSplitColumn,"splitcolumn",LSSdef);

//! unlike above, there is a separate ExpressionUserFunction for each <name,args> pair
class ExpressionUserFunction : public ExpressionFunctionBase {
public:
  ExpressionUserFunction(const function_spec&, argdef_t, const Expression*&);
  ExpressionFunctionBase* clone() const { return new ExpressionUserFunction(*this); }
  ExpressionFunctionBase* create(size_t) const;
  void get(LIST<double>&, const Expression&) const;
  void setexpression(const Expression&);
  static const BaseList<double>& getslot(size_t);
  static void clearstack();
  bool isconstant(int&, const Expression&) const;
  bool islist(unsigned short M) const { return ((islist_ >> M) & 1); } //!< \c true if argument \a M is a list (rather than scalar)
  
private:
  struct stack;
  static stack argstack_;
  const Expression*& exprp_;
  argdef_t islist_; //!< argument pattern signature
};
  
#endif
