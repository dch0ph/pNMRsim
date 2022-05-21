#include "Parser.h"
#include "parser_common.hpp"
#include "ScratchList.h"
#include "ttyio.h"
#include "expression_definition.hpp"
#if NMRSIM_USE_HASH
#include <unordered_set>
typedef std::unordered_set<const char*> hashed_set_t;
#else
#include "cmatrix_hash_set.h"
typedef hashed_set<> hashed_set_t;
#endif
#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <list>
#include <stack>
#include <sstream>
#include <set>

bool have_virtualdimensions=false;

typedef std::set<std::string> envstore_t;
envstore_t envstore; //!< store unique environment variables used - we assume that pointers are unique

namespace {
  const int NMRSIM_DEFAULT_PRECISION=5;
  const int NMRSIM_DEFAULT_MATRIXPRECISION=5;
  const int NMRSIM_DEFAULT_TIMEPRECISION=5;

  const char* getstoredenv(const char* name)
  {
    const char* env=getenv(name);
    if (env)
      envstore.insert(name);
    return env;
  }

}

void dumpenvironment(std::ostream& ostr)
{
  const envstore_t::const_iterator end(envstore.end());
  envstore_t::const_iterator start(envstore.begin());
  while (start!=end) {
    const char* name=(*start++).c_str();
    ostr << name << '=' << getenv(name) << ' ';
  }
}
        
//! line spacer in output .in files
#define OUTPUT_SPACER '\t'

bool isinteractive=false;

BaseWarning NMRsim_repeat_warning(BaseWarning::Always); //!< base for warnings that normally always repeat
BaseWarning NMRsim_once_warning(BaseWarning::FirstOnly); //!< base for warnings that normally are shown just once

ContextWarning<> emptyloop_warning("empty loop/control object",&NMRsim_repeat_warning);
ContextWarning<> postcontinuation_whitespace_warning("whitespace after \\ continuation",&NMRsim_repeat_warning);
ContextWarning<> extradollar_warning("variable names should supplied without leading $ (ignored)",&NMRsim_once_warning);
ContextWarning<> namequote_warning("variable names and other string literals ought to be quoted using \'",&NMRsim_once_warning);
ContextWarning<> filetail_ignored_warning("rest of file is being ignored",&NMRsim_repeat_warning);
ContextWarning<> function_redefinition_warning("redefining function: ",&NMRsim_repeat_warning);
//ContextWarning<> misplaced_sum_warning("sum || parameter not meaningful in this context (e.g. in proc rather than initialproc?)",&NMRsim_repeat_warning);
//ContextWarning<> misplaced_sum_warning("sum || parameter not meaningful in this context",&NMRsim_repeat_warning);
ContextWarning<> includeonce_withlist_warning("includeonce combined with list arguments.  Doesn't make much sense!",&NMRsim_repeat_warning);
//ContextWarning<> nofunctionarguments_warning("function with no arguments (#<n>): ",&NMRsim_once_warning);
//ContextWarning<> constfunction_warning("function evaluates to a constant! Function: ",&NMRsim_once_warning);
ContextWarning<> xy_variable_warning("using variable names $x or $y may cause confusion with x or y used as phases",&NMRsim_once_warning);
ContextWarning<> trailingcharacters_warning("trailing characters after } are being ignored",&NMRsim_once_warning);
ContextWarning<> ignoreddirective_warning("instruction ignored",&NMRsim_repeat_warning);

dimension_set::incommensurate_warning_t dimension_set::incommensurate_warning("data set size incommensurate with number of previously set rows");
dimension_set::incommensurate_warning_t dimension_set::incommensurate_warning2("list length is incommensurate with previous established length");

#define NMRSIM_OLD_FUNCTION '_'
static char linebuf[MAXLINE];
typedef LIST<std::string> pathlist_t;
pathlist_t pathstore;
varmap_t varmap;
size_t substitute_maxchars=MAXLINE/3;

bool nochecks=false;
bool abortonwarning=false;
bool noexecute=false;
int verbose_level=1;
bool silent=false;

static const double pi=M_PI;
SystemVariable<double*> v_pi("pi",const_cast<double*>(&pi),1.0,V_ISFIXED | V_ISCONST); //const_cast is harmless

static LIST<Variable*> dirty_stack; //!< stack of ::Variable with unset quantity pointers

ContextWarning<> sharedtarget_warning("dirty links have same apparent target",&NMRsim_repeat_warning);

void dirty_stack_push(Variable* varp)
{
  if (varp==NMRSIM_NULL)
    throw InternalError("dirty_stack_push");
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Pushing variable link onto dirty links stack (type=" << (varp->subsid) << ")\n";
  if ((verbose_level>1) && sharedtarget_warning.enabled()) { //!< only check if debugging since arrayed variables may have multiple links to the same target
    subsid_t subsid=varp->subsid;
    for (size_t i=dirty_stack.size();i--;) {
      if (subsid==dirty_stack(i)->subsid) {
	sharedtarget_warning.raise();
	break;
      }
    }
  }
  dirty_stack.push_back(varp);
}

hashed_set_t onceonlyset;
uservarmap_type vardefmap;

optional_map_t& get_optional_map()
{
  static optional_map_t optional_map;
  return optional_map;
}

command_Factory_t& get_Global_Factory() {
  static command_Factory_t Global_Factory;
  return Global_Factory; 
}

void parse_include(int);
void parse_setenv();

void initialise_precision()
{
  ostream_controller& ctrl(cmatrix_ostream_controller(std::cout));
  std::cout.precision(NMRSIM_DEFAULT_PRECISION);
  std::cerr.precision(NMRSIM_DEFAULT_PRECISION);
  ctrl.matrixprecision=NMRSIM_DEFAULT_MATRIXPRECISION;
  ctrl.timeprecision=NMRSIM_DEFAULT_TIMEPRECISION;
  ctrl.indexbase=1;
}

EXPR_DECLARE_VARIABLE_FUNCTION(FunctionValueof,"Valueof");
EXPR_DECLARE_VARIABLE_FUNCTION(FunctionValuesof,"Valuesof");
EXPR_DECLARE_VARIABLE_FUNCTION(FunctionErrorof,"Errorof");
EXPR_DECLARE_VARIABLE_FUNCTION(FunctionErrorsof,"Errorsof");

namespace {

  struct Parser_Proxy_ {
    Parser_Proxy_() {
      command_Factory_t& Global_Factory(get_Global_Factory());     
      Global_Factory["includeonce"]=par_t(&parse_include,1,true);
      Global_Factory["include"]=par_t(&parse_include,0,true);
      Global_Factory["variable"]=par_t(&parse_variable,true);
      Global_Factory["function"]=par_t(&parse_function,true);
      Global_Factory["setenv"]=par_t(&parse_setenv,true);
      Global_Factory["verbose"]=par_t(&parse_verbose,true);
      const static FunctionValueof Valueofobj;
      declare_function("Valueof",1U,Valueofobj);
      const static FunctionValuesof Valuesofobj;
      declare_function("Valuesof",1U,Valuesofobj);
      const static FunctionErrorof Errorofobj;
      declare_function("Errorof",1U,Errorofobj);
      const static FunctionErrorsof Errorsofobj;
      declare_function("Errorsof",1U,Errorsofobj);
    }
  };
    
  Parser_Proxy_ parser_proxy_; //!< declare globals

  //! return \c true if \a p is a (readable) directory
  bool isdirectory(const char* p)
  {
    struct stat statbuf;
    
    return stat(p,&statbuf) ? false : S_ISDIR(statbuf.st_mode);
  }

  //! return pointer to continuation point or NMRSIM_NULL if none
  char* continuation(char* buf)
  {
    char* p=buf+strlen(buf)-1;
    bool doneinc=false;
    for (;;) {
      if (p<buf)
	return NMRSIM_NULL;
      if (!isspace(*p)) {
	if (*p=='\\') {
	  if (doneinc)
	    postcontinuation_whitespace_warning.raise();
	  return p;
	}
	return NMRSIM_NULL;
      }
      if ((*p!='\n') && (*p!='\r'))
	doneinc=true; //!< spot any non-CR/LF whitespace
      p--;
    }
  }

	ThreadWarning<> unixpath_warning("File name matches Unix-style absolute path, when expecting Windows-style filenames: ",&NMRsim_once_warning);
	ThreadWarning<> windowspath_warning("File name matches Windows-style absolute path, when expecting Unix-style filenames: ",&NMRsim_once_warning);

  //! return \c true if \a p is an absolute pathname (won't catch DOS style X:\)
  bool isabsolute(const char* p)
    {
	const bool iswindowsabs=(*p=='\\') || ((p[1]==':') && (p[2]=='\\'));
	const bool isunixabs=(*p==DIR_SEPARATOR);
#ifdef WINDOWS_PATHS
	if (isunixabs)
		unixpath_warning.raise(p);
	const bool result=iswindowsabs;
#else
	if (iswindowsabs)
		windowspath_warning.raise(p);
	const bool result=isunixabs; 
#endif
	if ((verbose & VER_PARSE) && (verbose_level>1))
		parser_printthread() << p << " is absolute path: " << (result ? "Yes\n" : "No\n");
	return result;	
    }
}

ThreadWarning<> path_unreadable_warning("not a readable directory (skipped): ",&NMRsim_repeat_warning);
ThreadWarning<> path_emptyelement_warning("empty path element (skipped) in ",&NMRsim_repeat_warning);

//! parse path specification
void checkpathname(pathlist_t& pstore, const char* path)
{
  if (!path)
    path=".";
  std::string tmp(path);
  char* p=const_cast<char*>(tmp.c_str()); // a bit naughty - but we are not changing the string size
  char* token;
  while ((token=strtok(p,":"))) {
    p=NMRSIM_NULL;
    // don't strip white space?
    if (*token) {
      if (isdirectory(token)) {
	pstore.push_back(token);
	if ((verbose & VER_PARSE) && (verbose_level>1))
	  std::cout << "Adding to path: " << token << '\n';
      }
      else
	path_unreadable_warning.raise(token);      
    }
    else
      path_emptyelement_warning.raise(path);
  }
}

FILE* pathopen(const char* fname, const char* mode)
{
  if (isabsolute(fname))
    return fopen(fname,mode);

  char scr[MAXPATH];
  
  const pathlist_t::const_iterator end(pathstore.end());
  pathlist_t::const_iterator start(pathstore.begin());
  FILE* fp;
  while (start!=end) {
    if (snprintf(scr,MAXPATH,"%s%c%s",start->c_str(),DIR_SEPARATOR,fname)>=MAXPATH)
      throw Failed("pathopen: filename too long");
    fp=fopen(scr,mode);
    if (fp)
      return fp;
    ++start;
  }
  return NMRSIM_NULL;
}

void overflowerror()
{
  error_abort("overflow when constructing input line");
}

void pushback(char*& dest, const char* end, char s)
{
  if (dest+1>=end)
    overflowerror();
  *dest++=s;
  *dest='\0';
}

//! return \c true if \a ptr corresponds to a flag argument
/** Flag consists of - followed by only (but at least one) alphanumeric characters 
 */
bool parser_isflag(const char* ptr)
{
  if (!ptr || (*ptr!='-') || isspace(ptr[1]))
    return false;
  if (strcmp(ptr+1,"?")==0)
    return true;
  
  while (*(++ptr)) {
    if (isspace(*ptr))
      return true;
    if (!isalpha(*ptr) && (*ptr!='_'))
      return false;
  }
  return true;
}  

//Note "-" doesn't count as flag
bool parser_isflag()
{
  return parser_isflag(get_curline());
}

bool parser_isvariable_name(const char* ptr)
{
  if (!ptr)
    return false;
  return (*ptr=='\'');
}

bool parser_isvariable_name()
{
  return parser_isvariable_name(get_curline());
}

void pushback(char*& dest, const char* end, const char* source)
{
  const size_t n=strlen(source);
  if (dest+n>=end)
    overflowerror();
  //  const bool overlap=(dest<source) ? (dest+n>=source) : (source+n>=dest);     //  strcpy(dest,source);
  memmove(dest,source,n+1);
  dest+=n;
}

void pushback(char*& dest, const char* end, const char* source, const char* send)
{
  const size_t n=send-source;
  if (dest+n>=end)
    overflowerror();
  memmove(dest,source,n);
  dest+=n;
}

Mark::Mark()
  : stackpos(dirty_stack.size()) {}

bool Mark::needsflush() const
{
  return (dirty_stack.size()>stackpos);
}

ContextWarning<> Mark::unflushed_warning("mark destroyed before flush (ignore if triggered by exception)",&NMRsim_repeat_warning);

Mark::~Mark()
{
  if (needsflush())
    unflushed_warning.raise();
}

void Mark::flush(Setable* ptr)
{
  if ((verbose & VER_PARSE) && (verbose_level>1))
    std::cout << "Flushing dirty links stack\n";
  if (dirty_stack.size()<stackpos)
    throw InternalError("Mark::flush");
  int uses=0;
  while (needsflush()) {
    if (ptr) {
      Variable* curp=dirty_stack.back();
      curp->setptr(ptr);
      uses |= curp->variable().uses();
    }      
    dirty_stack.pop_back();
  }
  if (ptr)
    ptr->setable_uses(uses);
  //  doneflush=true;
}

void checkdirty()
{
  if (!dirty_stack.empty()) {
    dirty_stack.back()->print(std::cerr,false);
    std::cerr << std::endl;
    error_abort("Line finished with dirty objects uncleared");
  }
}

template<> int parse_raw(char*, Variable*, int);
template<> size_t parse_raw(char*, Variable*, int);
void parse_array(LIST<double>&, char*, int flags =F_DENYEXPR);
void parse_array_syntax(const char*, size_t, LIST<double>&, char*, int flags =F_DENYEXPR);

void Setable::set(const BaseList<double>&, subsid_t subsid)
{
  parser_printthread(std::cerr) << "Can't set this setable quantity (";
  printvariablename(std::cerr,subsid);
  std::cerr << ") to a list\n";
  error_abort();
}

inline std::ostream& operator<< (std::ostream& ostr, const parse_state& a)
{
  if (a.curline) {
    if (a.lines>1)
      ostr << "lines " << (a.curline-a.lines+1) << '-';
    else
      ostr << "line ";
    ostr << a.curline;
  }
  else
    ostr << "start"; //!< line no. is zero before first line is read
  return ostr << " of " << a.fname;
}

dimension_set::incommensurate_warning_t::incommensurate_warning_t(const char* namev)
  : ContextWarning<>(namev,&NMRsim_repeat_warning) {}
void dimension_set::incommensurate_warning_t::raise(const parse_state& state)
{
  std::ostringstream str(std::ostringstream::out);
  str << " at " << state;
  Warning<>::raise(str.str().c_str());
}

dimension_set array_dims,sum_dims;

struct funcdef
{  
  parse_state state_;
  FILE* fp;
  long fstart;
  size_t startline;
  size_t endline;
  bool active;

  funcdef() : fp(NMRSIM_NULL) {}
  funcdef(const char* fname_)
    : state_(fname_), startline(0), endline(0) { open(fname_); }
  
  funcdef(const char* fname_, long fstart_, size_t startline_, size_t endline_) 
    : state_(fname_), fstart(fstart_), startline(startline_), endline(endline_)
  { open(fname_); }

  ~funcdef();
  
  funcdef* define(const char*);
  void open(const char*);
  bool getline(char* linebuf, size_t max);
  const char* getfname() const { return state_.fname.c_str(); }
  void restart();

  const parse_state& state() const { return state_; }
};

std::ostream& operator<< (std::ostream& ostr, const funcdef& a)
{
  if (a.startline)
    return ostr << a.state_.fname << ':' << a.startline << '-' << a.endline;
  else
    return ostr << a.state_.fname;
}

class frame {
public:
  smartptr<funcdef,false> sourcep;
  LIST<std::string> args;
  std::list<std::string> coms;

  void open(funcdef*, const BaseList<char*>&,bool =false);
  bool skipcomments();
  bool isactive() { return !!sourcep; }
  const char* arg(size_t) const;
  const size_t nargs() const { return getargs().size(); }
  bool isargdefined(size_t n) const { return (n<=nargs()); }
  const parse_state& state() const { return sourcep->state(); }
  void insertline(const char* str) { coms.push_back(str); }
  const char* getfname() const { return sourcep->getfname(); }
  const BaseList<std::string> getargs() { return args; }

private:
  const LIST<std::string>& getargs() const { return args; }
  frame& operator= (const frame&); //disable
};
  

const char* frame::arg(size_t n) const
{
  if (!isargdefined(n)) {
    std::cerr << state() << ": invalid argument number: " << n << '\n';
    error_abort();
  }
  return getargs()(n-1).c_str();
}

frame* curframep(NMRSIM_NULL);
frame* topframep(NMRSIM_NULL);
std::string inname;

const char* parser_getfname()
{
  return inname.c_str();
}

double curgrat=0.0;
double defgrat=0.0; //!< set to unique nucleus grat for default

static const double rad_to_deg=180.0/M_PI;

template<> double parse_raw(char*, Variable*, int);

class Strpbrk {
public:
  Strpbrk(const char*);
  char* operator()(char*) const;
private:
  ScratchList<bool,256> tableok_;
};

Strpbrk::Strpbrk(const char* matchchars)
  : tableok_( sizeof(char) << 8, true)
{
  while (*matchchars) {
    tableok_(size_t(*matchchars))=false;
    matchchars++;
  }
  tableok_(size_t(0))=false; //!< convenient if also stop on '\0';
}
 
char* Strpbrk::operator()(char* p) const
{
  for (;tableok_(*p);p++);
  return *p ? p : NMRSIM_NULL;
}

char matchingbracket(char start)
{
  static char matchtable[sizeof(char)<<8]={'\0'};
  static bool havetable=false;
  if (!havetable) {
    matchtable[size_t('{')]='}';
    matchtable[size_t('}')]='X';
    matchtable[size_t('(')]=')';
    matchtable[size_t(')')]='X';
    matchtable[size_t('[')]=']';
    matchtable[size_t(']')]='X';
    matchtable[size_t('\'')]='\'';
    matchtable[size_t('"')]='"';
    matchtable[size_t('`')]='`';
    matchtable[size_t('|')]='|';
  }
//   switch (start) {
//   case '{':
//     return '}';
//   case '(':
//     return ')';
//     break;
//   case '[':
//     return ']';
//   case '"':
//     return '"';
//   case '\'':
//     return '\'';
//   }
  const char res(matchtable[size_t(start)]);
  switch (res) {
  case 'X':
    return '\0';
  case '\0':
    throw InternalError("matchingbracket");
  }
  return res;
}

bool checkbracket(char*& cptr, char start)
{
  if (*cptr!=start)
    return false;
  const char end=matchingbracket(start);
  char* last=cptr+strlen(cptr+1);
  if (*last!=end) {
    parser_printcontext() << "Failed to find matching " << end << '\n';
    error_abort();
  }
  cptr++;
  *last='\0';
  return true;
}

VariableBase* parse_double_variable_syntax(Variable& cvar, const char* syn, size_t narg, int flags)
{
  char* tok=get_token(0,"floating point");
  verify_syntax(tok,syn,narg);
  (void)parse_double_raw(&cvar,tok,flags | F_ALWAYSNEW);
  return cvar.valuep;
}

VariableBase* parse_double_variable(Variable& cvar, int flags)
{
  return parse_double_variable_syntax(cvar,NULL,0,flags);
}

const char* clean_variable_name(char* name)
{
  if (*name=='$') { 
    extradollar_warning.raise();
    name++;
  }
  else {
    if (!checkbracket(name,'\''))
      namequote_warning.raise();
  }
  return name;
}

ContextWarning<> inconsistent_const_warning("Inconsistency between constant status but non-zero attributes flag: ",&NMRsim_repeat_warning);
ContextWarning<> inconsistent_nouses_warning("Inconsistency between lack of const status and zero attributes flag: ",&NMRsim_repeat_warning);

RealVariable::RealVariable(const char* namev, VariableBase& varv)
  : name_(namev), varp_(&varv), issetable_(true), isconst_(varv.isconst()), isused_(false)
{
  if (!namev)
    throw InvalidParameter("RealVariable");
  uses_=varv.uses();
  if (inconsistent_nouses_warning.enabled() && !(varv.empty())) {
    const bool isinconsis= isconst_ ? (uses_!=0) : (uses_==0);
    if (isinconsis) {
      inconsistent_const_warning.raise(namev);
      varv.print(std::cerr,true);
      std::cerr << '\n';
    }
  }
}  

RealVariable& parse_variable_name(const char* cptr, bool flagused)
{
  const systemvarmap_type::const_iterator curp(systemvarmap.find(cptr));
  if (curp!=systemvarmap.end()) {
    RealVariable* gptr=dynamic_cast<RealVariable*>(curp->second);
    if (!gptr) {
      parser_printcontext() << '$' << cptr << " does not resolve to a numeric variable" << std::endl;
      throw Failed("parse_variable_name");
    }
    if (flagused)
      gptr->flagused();
    return *gptr;
  }
  return *findvariable(cptr,flagused);
}

void parse_array(LIST<double>& tmparray, char* cptr, int flags)
{
  chewwhite(cptr);
  if (*cptr=='\0') {
    tmparray.clear();
    return;
  }
//   if (!checkbracket(cptr,'{') && !checkbracket(cptr,'[')) {
//     parser_printcontext() << "Expected {} or [] delimited list while parsing: " << cptr << '\n';
//     error_abort();
//   }
    
  expr_grammar& expr_parser(get_expr_parser());
  expr_parser.reset(flags | F_ALLOWLIST);
  const parse_info<> info(parse(cptr,
			  expr_parser.use_parser<expr_grammar::rawlist_def>(),
			  space_p));
  if (info.full) {
    VariableBase& var(expr_parser.variable());
    if (!var.isconst()) {
      parser_printcontext() << "expecting real constant (list), but expression is not constant: " << cptr << '\n';
      error_abort();
    }
    tmparray.swap(var.value());
    return;
  }
  parser_printcontext() << "list parsing failed at: \"" << info.stop << "\"\n";
  error_abort();
}

void fudge_shift(double& val, size_t qualifier)
{
  val*=get_nmrfreq(qualifier)*1e-6;
}

void register_expression(const Variable& var)
{
  Variable& link=push_and_get(expressions,var);
  if (verbose & VER_PARSE)
    std::cout << "Created expression: " << *(var.valuep) << '\n';
  if (!link.ptr)
    dirty_stack_push(&link);
}

VarVariable::VarVariable(const Variable& varv, Parameter& parv, int arraywhichv)
  :
  Variable(varv), 
  arraywhich((arraywhichv==-1) ? 0U : arraywhichv),
  parameter(parv)
{
  const double step=variable().get_error(arraywhich);
  if (step==0.0)
    throw InternalError("VarVariable created from fixed parameter");
  parameter.error(step);
  //  updatebound();
}

// void VarVariable::updatebound()
// {  
//   const BaseBoundFunction* funcp=boundstate.function();
//   if (funcp==NMRSIM_NULL) {
//     parameter.unconstrain();
//     return;
//   }
//   const double val=get();
//   if (!(funcp->isvalid(val))) {
//     parser_printcontext() << "Current value (" << val << ") of parameter " << parameter.name() << " is outside bound range: " << boundstate << '\n';
//     error_abort();
//   }
//   parameter.constrain(*funcp);
// }

void register_variable(Variable& vardef, VariableBase& newvar)
{
  const BaseList<double> steps(newvar.get_errors());
  if (newvar.isexpr()) {
    register_expression(vardef);
    if (!steps.empty())
      throw InternalError("register_variable");
  }
  else {
    if (newvar.isarray()) { //'register' array sizes      
      const ListList<double> arrayvals(newvar.get_array());
      //   if ((arrayvals.size()>1) && (flags & F_DENYSUM))
      //	misplaced_sum_warning.raise();
      const size_t item_length=newvar.array_item_length();
      verify_arraysizes(arrayvals,newvar.sumtag(),newvar.arraytags(),item_length);
      Variable& link=push_and_get(varpars,vardef);
      if (!(link.ptr))
	dirty_stack_push(&link);
      //check for fitted elements
      for (size_t i=0;i<steps.size();i++) {
	if (steps(i)) {
	  if (item_length!=1)
	    error_abort("Can't have errors on list values");
	  //	  dirty_stack_push(create_fitting_variable(VarVariable(link,steps(i),i)));
	  parameter_list.push_back(Parameter(arrayvals.row()(i),steps(i)));
	  register_fitting_variable(new VarVariable(link,parameter_list.back(),i));
	}
      }
    }
    else {
      double err=0.0;
      switch (steps.size()) {
      case 0:
	break;
      case 1: 
	err=steps.front();
	break;
      default:
	throw InternalError("register_variable");
      }
      if (err) {
	parameter_list.push_back(Parameter(newvar.value().front(),steps.front()));
	register_fitting_variable(new VarVariable(vardef,parameter_list.back()));
      }
      else {
	if (newvar.isconst())
	  vardef.subsid=S_NONE; // flag not variable      
	//! NB does not remove any previous Variable from map
      }
    }
  }
}

double handle_variable(Variable* vardefp, int flags)
{
  VariableBase* newvarp=NMRSIM_NULL;
  if (vardefp && vardefp->ptr) {
    varmap_t::iterator iter=varmap.find(vardefp->key());
    if (iter!=varmap.end()) {
      if (!(flags & F_ALWAYSNEW))
	throw InternalError("Redefine demands ALWAYSNEW");
      newvarp=iter->second;
      if (!(newvarp->isconst()))
	error_abort("Can only redefine const variable");
    }
  }
  expr_grammar& expr_parser(get_expr_parser());
  VariableBase& cvar(expr_parser.variable());
  const BaseList<double> steps(cvar.get_errors());
  const double value=cvar.value().size() ? cvar.value().front() : 0.0; //0 if empty array

  if (cvar.issimple() && (cvar.value().size()<2) && steps.empty() && !(flags & F_ALWAYSNEW)) {
    //    if (cvar.value().size()>2) //!< allow 0 or 1 element
    //  throw InternalError("handle_variable:2");
    if (vardefp) {
      vardefp->subsid=S_NONE;
      vardefp->valuep=NMRSIM_NULL;
    }
    expr_parser.builder_.reset(); //reset parser so another variable can be parsed
    return value;
  }
  if (!vardefp || (vardefp->subsid==S_NONE))
    error_abort("Can't use non-constant value here");

  if (newvarp) {
    if (!(cvar.isconst()))
      error_abort("Can only overwrite existing variable with constant expression");
  }
  else
    newvarp=new VariableBase();

  newvarp->swapin(cvar); //swap in contents
  vardefp->valuep=newvarp;

  const int uses=newvarp->uses();

  if ((uses & A_SUM) && (flags & F_DENYSUM))
    error_abort("Dependency on sum || parameter not permitted in this context");
    //    misplaced_sum_warning.raise();
  
  if ((uses & A_VAR) && (flags & F_DENYVAR))
    error_abort("Dependency on fitting variable not permitted in this context");

  if ((uses & A_ARRAY) && (flags & F_DENYARRAY))
    error_abort("Dependency on arrayed variable not permitted in this context");

  if (newvarp->isexpr() && (uses==A_VAR))
    newvarp->update(); //!< if expression only depends on fitting variables, force evaluation - added on experimental basis 7/2/14

  register_variable(*vardefp,*newvarp);

  expr_parser.builder_.reset(); //reset parser so another variable can be parsed

  return value;
}

double handle_variable(int flags, size_t subsid)
{
  Variable cvar(subsid_t(subsid),NMRSIM_NULL);
  return handle_variable(&cvar,flags);
}

double parse_double_raw(Variable* vardefp, char* cptr, int flags)
{
  if (!vardefp && (flags & F_ALWAYSNEW))
    throw InternalError("parse_double_raw");

  expr_grammar& expr_parser(get_expr_parser());
  expr_parser.reset(flags);
  const parse_info<> info(parse(cptr,
			  expr_parser.use_parser<expr_grammar::variable_def>(),
			  space_p));
  if (!info.full) {
    parser_printcontext() << "Failed to parse numerical / expression value at: \"" << info.stop << "\"\n";
    error_abort();
  }  
  return handle_variable(vardefp,flags);
}
  
double parse_shift(Variable* vardef, double grat, int flags)
{
  curgrat=grat;
  //  double retval=parse_double(vardef,flags | F_ISSHIFT);
  double retval=parse_double(vardef,flags);
  curgrat=defgrat; //unset to trap attempts to parse 'p' without defining nucleus
  return retval;
}

template<> double parse_raw(char* cptr, Variable* cvarp, int flags)
{
  return parse_double_raw(cvarp,cptr,flags);
}

double parse_double(Variable* vardef, int flags)
{
  return parse_double_syntax(NULL,0,vardef,flags);
}

double parse_double_syntax(const char* syn, size_t narg, Variable* vardef, int flags)
{
  char* tok=get_token(0,"floating point");
  verify_syntax(tok,syn,narg);
  return parse_raw<double>(tok,vardef,flags);
}

size_t parser_verbose_level()
{
  return (verbose & VER_PARSE) ? verbose_level : 0;
}

funcdef::~funcdef()
{
  if (parser_verbose_level()>1)
    std::cout << "Closing frame: " << state_.fname << '\n';
  if ((fp!=stdin) && fp)
    fclose(fp);
}

void funcdef::open(const char* fname_)
{
  static const char stdinname[]="<stdin>";
  if (strcmp(fname_,"-")==0) {
    state_.fname=stdinname;
    fp=stdin;
  }
  else {//postpone opening file until required
    state_.fname=fname_;
    fp=NMRSIM_NULL;
  }
  active=false;
}

void parser_insertline(const char* buf)
{
  if (verbose & VER_PARSE)
    std::cout << "Pushing line onto input stack: " << buf << '\n';
  curframep->insertline(buf);
}

bool funcdef::getline(char* linebuf, size_t max)
{
  if (!active)
    return false;

  const char* cptr=NMRSIM_NULL;
  const bool usestack=!(curframep->coms.empty());
  size_t lines=0;
  if (usestack) {
    cptr=curframep->coms.front().c_str();
    lines=1;
  }
  else {
    cptr=linebuf;
    char* dest=linebuf;
    int left=max;   
    while (active) {
      if (fgets(dest,left,fp)) {
	lines++;
	state_.curline++;
	char* contp=continuation(dest);
	if (state_.curline+1==endline) {//reached end of parse lines?
	  if (contp)
	    error_abort("continuation \\ at end of include block");
	  active=false;
	  if (verbose & VER_PARSE)
	    parser_printthread() << "End of frame: " << (*this) << " (reached stop point at line " << endline << ")\n";	    
	}
	else {
	  if (contp) {
	    left-=contp-dest;
	    if (left<0)
	      throw InternalError("funcdef::getline");
	    dest=contp;
	  }
	  else
	    break;
	} 
      }
      else {
	if (verbose & VER_PARSE)
	  parser_printthread() << "End of frame: " << (*this) << " (end-of-file / error state)\n";
	active=false;
      }
    }
  }

  if (lines==0) { // Failed to find input line
    active=false;
    return false;
  }
  state_.lines=lines;

  if (verbose & VER_PARSE)
    parser_printcontext() << "Read line: " << cptr << '\n';
  const size_t n=strlen(cptr);
  if (n>=max-1) {
    parser_printcontext() << "input line too long (max of " << max << " chars)\n";
    error_abort();
  }
  if (usestack) {
    strcpy(linebuf,cptr);
    curframep->coms.pop_front();
  }

  if (n) { // if non-zero string, strip terminating white space
    char* endchar=linebuf+n-1;
    while (isspace(*endchar))
      endchar--;
    endchar[1]='\0';
  }
  return true;
}

funcdef* funcdef::define(const char* blockname)
{
  char linebuf[MAXLINE];
  const long newfstart=ftell(fp);
  const size_t newstartline=state_.curline+1;
  size_t newendline;
  for (;;) {
    if (!fgets(linebuf,sizeof(linebuf),fp)) {
      parser_printcontext() << "Exhausted file when looking for end of " << blockname << '\n';
      error_abort();
    }
    if (continuation(linebuf))
      error_abort("\\ continuation not allowed in block definition");
    state_.curline++;
    char* cptr=linebuf;
    chewwhite(cptr);
    if (*cptr=='}') {
      newendline=state_.curline;
      break;
    }
  }
  funcdef* newfuncp= new funcdef(state_.fname.c_str(),newfstart,newstartline,newendline);
  if (verbose & VER_PARSE)
    std::cout << blockname << " defined: " << (*newfuncp) << " at file pointer: " << newfstart << '\n';
  return newfuncp;
}

void funcdef::restart()
{
  if (active)
    throw Failed("Recursive file open");
  if (!fp) {
    fp=pathopen(state_.fname.c_str(),"r");
    if (!fp) {
      parser_printcontext() << "failed to open " << state_.fname << std::endl;
      error_abort(ERR_FAILED);
    }
  }
  if (startline)
    fseek(fp,fstart,SEEK_SET);
  active=true;
  state_.curline=startline ? startline-1 : 0;
}

typedef MAPTYPE(funcdef*) funcdefmap_type;
funcdefmap_type funcdefmap;

void declare_builtin_block(const char* name)
{
  funcdefmap[name]=NMRSIM_NULL;
}

void par_t::operator()()
{
  if (!allowmultiple) {
    if (done)
      error_abort("Can't reuse this directive");

    done=true;
  }
  if (parsefunc)
    (*parsefunc)();
  else {
    if (parsefunc_qual)
      (*parsefunc_qual)(ident);
    else
      throw InternalError("par_t");
  }
}

Splitter::Splitter(char* source, const char* charsv, bool uniquev)
  : curptr(source), chars(charsv), unique(uniquev)
{ 
  sepchar='\0';
  chewwhite(curptr);
}
  
char* Splitter::next(char& retsepchar)
{
  if (atend())
    return NMRSIM_NULL;

  retsepchar=sepchar;

  s=NMRSIM_NULL;
  for (const char* trychar=chars;*trychar;trychar++) {
    char *cur(strchr((*trychar=='-') ? curptr+1 : curptr,*trychar));
    if ( cur && ((cur<s) || !s))
      s=cur;
  }
  char* end;
  if (s) {
    sepchar=*s;
    end=s-1;
    if (unique && retsepchar && (sepchar!=retsepchar)) {
      parser_printcontext() << "cannot mix item separators: changed from " << retsepchar << " to " << sepchar  << '\n';
      error_abort();
    }
  }
  else {
    end=curptr+strlen(curptr)-1;
    sepchar='\0';
  }

  while (isspace(*end))
    end--;
  end[1]='\0'; //remove trailing spaces

  char* oldptr=curptr;
  if (!*(oldptr)) //separator at start of string?
    throw Failed("Couldn't parse expression");

  if (sepchar) {
    curptr=s+1;
    chewwhite(curptr);
  }
  else
    curptr=NMRSIM_NULL;
  return oldptr;
}

//Parser state
const char* blockname=NMRSIM_NULL;
bool insection=false;
bool needfresh=false;
bool lastblockwasdefinition=false;
bool found_block=false;
bool have_block=false;
const char** keynames;
size_t nkeys;
bool isopt;
std::stack< frame, LIST<frame> > fps;

static bool popinclude();

bool parser_foundlast() { return found_block; }
bool parser_inblock() { return insection; }

int global_flags=0;

void parser_init(const char* fname, int argc, char* argv[],int global_flagsv)
{
  BaseList<char*> args(argc,argv);
  checkpathname(pathstore,register_and_getenv("NMRSIM_PATH"));
  fps.push(frame());
  inname=fname;
  topframep=curframep=&(fps.top());  
  curframep->open(new funcdef(fname),args);
  global_flags=global_flagsv;

  add_systemvarmap(v_pi);
  initialise_precision();
}

bool parser_isopen()
{
  return curframep ? curframep->isactive() : false;
}

const char delims[3]={ ' ','\t','\0' };
static const double deg_to_rad=M_PI/180.0;

static char* placeholder=NMRSIM_NULL;
static char* lasttoken=NMRSIM_NULL;

#ifndef HAVE_STRTOK_R
// replacement for strtok_r if missing
char* strtok_r( char *str, const char *delim, char **save_ptr ) 
{
  if (!str)
    str=*save_ptr;

  str+=strspn(str,delim); //! skip leading delimiters

  if (!*str) {
    *save_ptr=str;
    return NMRSIM_NULL;
  }

  char* ptr = str;

  while (*ptr && !strchr(delim,*ptr))
    ++ptr;

  *save_ptr = *ptr ? ptr + 1 : ptr;
  *ptr = '\0';

  return str;
}
#endif

ContextWarning<> redundantquote_warning("Quoting arguments with \" is no longer necessary",&NMRsim_once_warning);

char* token_search(char*& cptr, char endpar ='\0', bool optional =false)
{
  static const char paras[]="\"\'`(){}[]|";
  static const Strpbrk strpbrk_paraset(paras); //!< use cached strpbrk for larger number of characters
  static LIST<char> parstack(size_t(10));  
  static char pdelims[2]={'\0'};

  const char* ldelims=delims;
  //  bool ispar=false;
  if (endpar) {
    if (strchr(paras,endpar)) {
      if (optional)
	throw InternalError("token_search: optional token separator can't be parenthesis/quote");
      //  ispar=true;
      parstack.create(1,endpar);
    }
    pdelims[0]=endpar;
    ldelims=pdelims;
  }
  else
    parstack.create(size_t(0)); //reset paren stack

  char* lptr=cptr;
  char* nextdelimp=NMRSIM_NULL;
  bool stripquote=false;
  char insideliteral='\0';
  for (;;) {
    if (lptr>nextdelimp) //!< need to search for end delimiter?
      nextdelimp=strpbrk(lptr,ldelims); 
    if (!nextdelimp) {
      if (endpar) {
	if (optional)
	  return NMRSIM_NULL;
	parser_printcontext() << "Failed to find expected " << endpar << '\n';
	error_abort();
      }
      nextdelimp=lptr+strlen(lptr);
    }
    char* nextpar=NMRSIM_NULL;
    if (insideliteral) {
      for (;;) {
	nextpar=strchr(lptr,insideliteral);
	if ((nextpar<=lptr) || (nextpar[-1]!='\\')) //not quoted?
	  break;
	lptr=nextpar+1;
      }
    }
    else {
      const bool shortcircuit=parstack.empty();
      char delimchar='\0';
      if (shortcircuit) {
	delimchar=*nextdelimp;
	*nextdelimp='\0';
      }
      nextpar=strpbrk_paraset(lptr);
      if (shortcircuit)
	*nextdelimp=delimchar; //!< restore
    }
    
    if (nextpar==NMRSIM_NULL)
      break;
    if (insideliteral)
      insideliteral='\0';
    else {
      bool islit=false;
      switch (*nextpar) {
      case '\'': case '`':
	islit=true;
	break;
      case '\"':
	if ((!endpar) && (nextpar==cptr)) {
	  stripquote=true;
	  redundantquote_warning.raise();
	}
	else {
	  if (!(stripquote && (nextpar==nextdelimp-1)))
	    islit=true;
	}
	break;
      }
      if (islit)
	insideliteral=*nextpar;
      else {
	if (!(parstack.empty()) && (*nextpar==parstack.back())) {
	  parstack.pop_back();
	  if (*nextpar==endpar && (parstack.empty()))
	    return nextpar;
	}
	else {
	  const char match=matchingbracket(*nextpar);
	  if (!match) {
	    parser_printcontext() << "Unexpected closing parenthesis (" << *nextpar << ')';
	    if (!(parstack.empty()))
	      std::cerr << ". Expecting " << parstack.back();
	    std::cerr << " when parsing " << cptr << '\n';
	    error_abort();
	  }
	  parstack.push_back(match);
	}
      }
    }
    lptr=nextpar+1;
  }
  if (parstack.size()) {
    parser_printcontext() << "unterminated quote: missing trailing " << parstack.back() << '\n';
    error_abort();
  }
  if (stripquote) {
    if (nextdelimp[-1]!='\"')
      error_abort("argument quoted with \" is not isolated phrase");
    nextdelimp[-1]='\0';
    cptr++;
  }
  return nextdelimp;
}

char* get_token_raw(char*& holder, int flags =0, const char* argtype =NMRSIM_NULL)
{
  char* cptr=holder;
  if (cptr) {
    chewwhite(cptr);
    char* nextdelimp=token_search(cptr);
    if (*nextdelimp) {
      *nextdelimp='\0';
      holder=nextdelimp+1;
    }
    else 
      holder=NMRSIM_NULL;
    
//     if (*cptr=='\"') {
//       holder=strchr(++cptr,'\"');
//       if (!holder)
// 	error_abort("Unbalanced \"");
//       *holder++='\0';
//     }
//     else
//       cptr=strtok_r(NMRSIM_NULL,delims,&holder);
  }
  if ((cptr==NMRSIM_NULL) && !(flags & F_ALLOWMISSING)) {
    if (!argtype)
      argtype="";
    parser_printcontext(std::cerr) << "Missing " << (*argtype) << " argument\n";
    error_abort();
  }
  return cptr;
}

char* get_token(int flags, const char* type)
{
  lasttoken=get_token_raw(placeholder,flags,type);
  return lasttoken;
}

bool are_left()
{
  if (!placeholder)
    return false;
  chewwhite(placeholder);
  return (*placeholder!='\0');
}

size_t count_left(bool countflags)
{
  char* inpargs=get_curline();
  if (!inpargs)
    return 0;
  ScratchList<char> copyargs( BaseList<char>( 1+strlen(inpargs) ,inpargs) );
  //need to copy arguments as strtok will break up...
  char* localholder=copyargs.vector();
  size_t count=0;
  char* ptr;
  while ((ptr=get_token_raw(localholder,F_ALLOWMISSING))) {
    if (countflags || !parser_isflag(ptr))
      count++;
  }
  return count;
}

void parser_newblock(const char* blockname_, bool isopt_, const char* keynames_[],size_t nkeys_)
{
  if (!parser_isopen()) {
    if (isopt)
      return;
    parser_printcontext() << "Expecting to find " << blockname_ << "block, but reached end of input file!\n";
    error_abort();
  }
  if (insection) {
    char scratch[256];
    snprintf(scratch,sizeof(scratch),"parser_newblock: can't start new block (%s) inside active block (%s)",blockname_,blockname); 
    throw Failed(scratch);
  }
  blockname=blockname_;
  keynames=keynames_;
  nkeys=nkeys_;
  isopt=isopt_;
  found_block=false;
}

bool popinclude()
{
  fps.pop();
  if (fps.empty()) {
    curframep=topframep=NMRSIM_NULL; //finished all files?
    if (insection) {
      parser_printcontext() << "Reached end of input while still parsing " << blockname << " block!\n";
      error_abort();
    }
    return false;
  }
  curframep=&(fps.top());
  return true;
}

void frame::open(funcdef* sourcep_, const BaseList<char*>& largs, bool persist)
{
  if (!largs.empty())
    args=largs; //store arguments

  sourcep.reset(sourcep_,persist ? mxflag::nondynamic : mxflag::normal);
  if (verbose & VER_PARSE) {
    std::cout << "Opening " << (*sourcep);
    if (args.size()) {
      std::cout << " with arg list: ";
      for (size_t i=0;i<args.size();i++) {
	const char* arg(args(i).c_str());
	if (*arg)
	  std::cout << arg << ' ';
	else
	  std::cout << "'' ";
      }
    }
    std::cout << '\n';
  }
  if (sourcep->active) {
    parser_printcontext() << "Recursive call to " << (*sourcep) << '\n';
    error_abort();
  }
  sourcep->restart();
}

bool frame::skipcomments()
{
  lasttoken=NMRSIM_NULL; //!< clear references to previous lines
  char rawlinebuf[MAXLINE];
  while (sourcep->getline(rawlinebuf,MAXLINE)) {
    substitute_string(linebuf,sizeof(linebuf),rawlinebuf,SUB_NUMERIC | SUB_ABORTHASH);
    if (linebuf[0]!='#')
      return true;
  }
  return false;
}

void parser_getinteractive(const char* intro, std::istream& istr)
{
  char rawlinebuf[MAXLINE];
  for (;;) {
    std::cout << intro;
    getline(rawlinebuf,MAXLINE);
    substitute_string(linebuf,sizeof(linebuf),rawlinebuf,SUB_NUMERIC | SUB_ABORTHASH);
    if (*linebuf) {
      set_curline(linebuf);
      return;
    }
  }
}

ContextWarning<> lastblockdefinition_warning("final block was parsed as an internal include - mislabelled block?  Parser was expecting: ",&NMRsim_repeat_warning); 

void parser_flush()
{
  if (insection || !needfresh)
    error_abort("closing file before block parsed (use -verbose:parse to diagnose)");
  if (lastblockwasdefinition)
    lastblockdefinition_warning.raise(blockname);
  for (;;) {
    if (!curframep)
      return;
    if (curframep->skipcomments())
      break; //!< oh dear, found more text
    if (!popinclude())
      return; //!< OK, nothing left
  }
  filetail_ignored_warning.raise();
}


namespace {
  struct tmparg_t {
    char str[10];
    char* operator()() { return str; }
    const char* operator()() const { return str; }
    tmparg_t() { *str='\0'; }
  };
}
    
void parse_include(int onceonly)
{
  const char* fname(parse_string());

  LIST<char*> args;
  LIST< BaseList<double> > varargs;
  LIST<tmparg_t> tmpargs;
  size_t maxrep=1;
  while (are_left()) {
    char* tokptr=get_token();
    const bool isquoted=checkbracket(tokptr,'\'');   
    BaseList<double> argptr;
    if (!isquoted) {
      Variable cvar(S_ARG1);
      const double v=parse_double_raw(&cvar,tokptr,F_DENYVAR | F_DENYARRAY | F_ALLOWLIST | F_ALWAYSNEW);
      if (cvar.subsid!=S_NONE)
	throw InternalError("parse_include");
      tmpargs.push_back(tmparg_t()); //create new buffer space
      tokptr=tmpargs.back()();
      argptr.create(cvar.valuep->value());
      const size_t n=argptr.size();
      if (n==1)
	snprintf_prec(tokptr,sizeof(tmparg_t),v,true);
      else {
	if ((n!=1) && (n!=maxrep)) {
	  if (maxrep==1)
	    maxrep=n;
	  else
	    throw Failed("include: arguments have incompatible lengths");
	}
      }
    }      
    args.push_back(tokptr);
    varargs.push_back(argptr);
  }

  if (onceonly) {
    if (maxrep!=1)
      includeonce_withlist_warning.raise();

    const bool already=(onceonlyset.count(fname)!=0);
    if (verbose & VER_PARSE)
      std::cout << "File " << fname << " already included: " << (already ? "Yes\n" : "No\n");
    if (already)
      return; //already included
    onceonlyset.insert(fname);
  }
    
  //! empty list?
  if (!maxrep)
    return;

  static char tlinebuf[MAXLINE]="include XX"; // XXX is fence-post check
  if (maxrep>1) {
    char* const endptr=tlinebuf+sizeof(tlinebuf);
    for (size_t i=0;i<maxrep;i++) {
      char* curptr=tlinebuf+sizeof("include ")-1;
      pushback(curptr,endptr,fname);
      for (size_t j=0;j<varargs.size();j++) {
	const BaseList<double>& largs(varargs(j));
	char* tokptr=args(j);
	if (largs.size()>1) //!< need to change argument?
	  snprintf_prec(tokptr,sizeof(tmparg_t),largs(i),true); //!< re-format numerical argument
	pushback(curptr,endptr," \'"); //!< inter-argument space
	pushback(curptr,endptr,tokptr); //!< add argument
	pushback(curptr,endptr,'\''); //!< need to quote to prevent further evaluation
      }
      parser_insertline(tlinebuf);
    }
    return;
  }
    
  fps.push(frame());
  frame& newframe=fps.top();
  const funcdefmap_type::const_iterator curp(funcdefmap.find(fname));
  const bool newfile(curp==funcdefmap.end());
  funcdef* sourcep=newfile ? new funcdef(fname) : curp->second;
  if (!sourcep) {
    parser_printcontext() << "Attempt to include built-in block " << fname << '\n';
    error_abort();
  }
  newframe.open(sourcep,args,!newfile);
  curframep=&newframe;
}

#define CHECKWORD(A,B)\
  ((STRNCMP(A,B)==0) && (isspace(cptr[sizeof(B)-1]) || !cptr[sizeof(B)-1]))


bool parser_getline(bool executeglobal)
{
  if (!blockname)
    error_abort("parser_getline: no block started");
  if (!curframep)
    return false;
  command_Factory_t& Global_Factory(get_Global_Factory());
  for (;;) {
    if (needfresh && !(curframep->skipcomments())) {
      if (!popinclude())
	return false;
      continue;
    }
    needfresh=true;
    char* cptr=linebuf;
    chewwhite(cptr);
    if (*cptr=='}') {
      if (!insection)
	error_abort("Unexpected } block end");
      cptr++;
      chewwhite(cptr);
      if (*cptr)
	trailingcharacters_warning.raise();
      insection=false;
      if (verbose & VER_PARSE)
	parser_printcontext() << "finishing parsing block: " << blockname << '\n';
      return false; // end block
    }
    if (!*cptr)
      continue;

    if (verbose & VER_PARSE)
      parser_printthread() << "Parsing " << (curframep->state()) << ": " << cptr << std::endl;

    char* tail=cptr+1;
    while (*tail && !isspace(*tail))
      tail++;
    const char storechar=*tail;
    *tail='\0';

//     if (strcmp(cptr,"include")==0) {
//       *tail=storechar;
//       set_curline(tail);
//       parse_include(false);
//       continue;
//     }
//     if (strcmp(cptr,"includeonce")==0) {
//       *tail=storechar;
//       set_curline(tail);
//       parse_include(true);
//       continue;
//     }
//    const bool isglobal=Global_Factory.count(cptr);
    command_Factory_t::iterator curp=Global_Factory.find(cptr);
    const bool isglobal=(curp!=Global_Factory.end());
    if (isglobal) {
      set_curline(storechar ? tail+1 : tail);
      if ((curp->second).parsefunc==parse_setenv) { //!< always execute setenv
	parse_setenv(); //!< execute setenv *AND* return line if in scan mode
	continue;
      }
      else {
	const bool isinclude=((curp->second).parsefunc_qual==parse_include);
	if (isinclude || executeglobal) {
	  (curp->second)();
	  parser_checkfinished();
	  continue;
	}
      }
    }
    *tail=storechar;
    set_curline(cptr);
//     if (CHECKWORD(cptr,"include")) {
//       set_curline(cptr+sizeof("include"));
//       parse_include(false);
//       continue;
//     }
//     if (CHECKWORD(cptr,"includeonce")) {
//       set_curline(cptr+sizeof("includeonce"));
//       parse_include(true);
//       continue;
//     }
    if (!isglobal && !insection) {
      char tmp[128];
      char tmpc;
      if ((sscanf(cptr,"%s %c ",tmp,&tmpc)!=2) || (tmpc!='{'))
	error_abort("Expecting block start");
      if (strcmp(tmp,blockname)==0) {
	have_block=true;
	found_block=true;
	insection=true;
	lastblockwasdefinition=false;
	if (verbose & VER_PARSE)
	  parser_printcontext() << "found expected parsing block: " << blockname << '\n';
	continue;
      }
      //      if (strcmp(tmp,"initialproc")==0)
      //	error_abort("Can't define block named initialproc (functionality withdrawn) - contents should be in proc");

      const funcdefmap_type::const_iterator curp(funcdefmap.find(tmp));
      if (curp!=funcdefmap.end()) {
	if (curp->second)
	  function_redefinition_warning.raise(tmp);
	else {
	  if (isopt) {
	    have_block=true;
	    needfresh=false;
	    return false; //found different block
	  }
	  parser_printcontext() << "expecting " << blockname << " block, but found different built-in block: " << tmp << std::endl;
	  error_abort();
	}
      }
      //define new block
      funcdefmap[tmp]=curframep->sourcep->define(tmp);
      lastblockwasdefinition=true;
      continue;
    }
    return true;
  }
}

char* get_curline() 
{
  char* ptr=placeholder;
  if (ptr)
    chewwhite(ptr);
  return ptr;
}

//bit dodgy
void set_curline(char* ptr) {
  placeholder=ptr;
}

//! return \c true if \a str does not contain lower case characters
bool isenvironment(const char* str)
{
  while (*str) {
    if (islower(*str))
      return false;
    str++;
  }
  return true;
}

inline char* strchrlimit(char* s, char* tail, char c)
{
  char* p=strchr(s,c);
  return (!p || (p>=tail)) ? NMRSIM_NULL : p;
}

namespace {

  bool isvarconstant(const RealVariable* vp) { return vp->isconstant(); }

  bool isvarconstant(const SystemVariableBase* vp)
  {
    const RealVariable* vr=dynamic_cast<const RealVariable*>(vp);
    return vr ? isvarconstant(vr) : true;
  }

  int getuses(const RealVariable* vp) { return vp->uses(); }

  int getuses(const SystemVariableBase* vp)
  {
    const RealVariable* vr=dynamic_cast<const RealVariable*>(vp);
    return vr ? getuses(vr) : 0;
  }
};

ThreadWarning<> nonconstsubs_warning("Formatting non-constant variable outside of main execution loop at $",&NMRsim_once_warning);

char* substitute_string(char* buffer, int nchars, const char* rawcptr, int flags)
{
  int uses=0;
  return substitute_string(buffer, nchars, rawcptr, flags, uses);
}

char* substitute_string(char* buffer, int nchars, const char* rawcptr, int flags, int& accuses)
{
  char* cptr=const_cast<char*>(rawcptr); //!< a bit of a hack, but we only tweak the input string and don't fundamentally violate its constness
  char* curptr=buffer;
  char* endptr=buffer+nchars;
  char tmp[32];
  char tmpbuffer[MAXLINE];

  if ((verbose & VER_PARSE) && (verbose_level>1)) {
    const void* bufferasvoid=buffer;
    const void* sourceasvoid=rawcptr;
    std::cout << "Calling subsitute_string with destination buffer " << bufferasvoid << " (" << nchars << " characters) using " << sourceasvoid << " (" << rawcptr << ")\n";
  }

  const bool ignoreformat=(flags & SUB_FULLPREC);

  *buffer='\0';
  if (!cptr)
    return buffer;

  while (cptr && *cptr) {
    char* dollar=strchr(cptr,'$');
    if (!dollar) {
      pushback(curptr,endptr,cptr);
      return curptr;
    }
    // skip over additional $
    while (dollar[1]=='$')
      dollar++;

    pushback(curptr,endptr,cptr,dollar);

    dollar++; //point to post$

    if ((curptr!=buffer) && (curptr[-1]=='\\')) { //found \$
      if (flags & SUB_ESCAPE)
	curptr--; //move back to write over
      *curptr++='$';
      cptr=dollar;
      continue;
    }
    bool hasparen=false;
    int skip=0;
    char* tail=NMRSIM_NULL;
    if (*dollar=='(') {
      hasparen=true;
//       tail=strchr(++dollar,')');
//       if (!tail)
// 	error_abort("Unbalanced () in $(<variable>)");
      dollar++;
      tail=token_search(dollar,')');
      if (!tail)
	throw InternalError("token_search");
      skip=1;
    }

    const bool isstar=(*dollar=='*');
    if (isstar)
      dollar++;
    char* origcptr=dollar;

    long n=-1;
    if (flags & SUB_NUMERIC) {
      char* ltail;
      n=strtoul(dollar,&ltail,10);
      if (ltail==dollar)
	n=isstar ? 1 : -1; // flag invalid unless $* (equivalent to $1)
      if ((n>=0) && !hasparen)
	tail=ltail;
    }
    // establish end of name
    if (n<0) { //didn't find digit
      if (isstar)
	error_abort("Can only use $* with argument variables e.g. $*1");
      if (!tail) {//!< unless end already known, look for white space or (
	char* endptr=dollar;
	while (isalpha(*endptr) || isdigit(*endptr) || (*endptr=='_'))
	  endptr++;
	if (*endptr)
	  tail=endptr;
      }

      if (tail) {
	//	char* oldptr=curptr;
	char* dollarptr=tmpbuffer; //!< tmpbuffer introduced 25/3/15 otherwise recursion leads to brokeness
	pushback(dollarptr,tmpbuffer+sizeof(tmpbuffer),dollar,tail); //copy into buffer (scratch copy) 
	*dollarptr='\0'; //null terminate
	dollar=tmpbuffer;
	if (!(flags & SUB_NUMERIC))
	  tail+=skip; //!< done later for substitution variables (yuk)
      }      
      // catch bad variable name e.g. $.
      if (*dollar=='\0')
	error_abort("invalid use of $: following variable name can contain only alphanumeric characters (letters and numbers)");
    }
    char* qualifier=NMRSIM_NULL;
    char* haspositive=NMRSIM_NULL;
    bool emptycountsasundefined=false;
    if (hasparen) {            
      qualifier=(n<0) ? strchr(dollar,'?') : strchrlimit(dollar,tail,'?');
      char* altqual=(n<0) ? strchr(dollar,'!') : strchrlimit(dollar,tail,'!');
      if (altqual && (!qualifier || (altqual<qualifier))) {
	qualifier=altqual;
	emptycountsasundefined=true;
      }

      if (qualifier) {
	*qualifier++='\0'; //!< separate name and qualifier
	char tailchar;
	if (tail) {
	  tailchar=*tail;
	  *tail='\0';
	}
	haspositive=token_search(qualifier,':',true);
	if (tail)
	  *tail=tailchar;
	if (haspositive) {
	  *haspositive++='\0';
	  std::swap(qualifier,haspositive);
	}
      }
    }

    if (flags & SUB_NUMERIC) {      
      const char* varval=NMRSIM_NULL;
      if (n>=0) {
	if (curframep->isargdefined(n))
	  varval=curframep->arg(n);	  
      }
      else {
	if (isenvironment(dollar))
	  varval=getstoredenv(dollar);
	else {
	  pushback(curptr,endptr,'$');
	  if (hasparen)
	    *curptr++='(';
	  
	  cptr=origcptr; //not found, continue copy
	  continue;
	}
      }
      if (varval && (*varval=='\0') && emptycountsasundefined)
	varval=NMRSIM_NULL;

      cptr=tail+skip;
      if (!varval) { //not defined?
	if (qualifier) {
	  if (n>=0) {
	    if (tail<qualifier)
	      throw InternalError("substitute_string");
	    *tail='\0';
	  }
	  //	  *parend='\0';
	  int locflags=flags;
	  if (curptr!=buffer)
	    locflags &= ~SUB_ABORTHASH;	  
	  curptr=substitute_string(curptr,endptr-curptr,qualifier,locflags, accuses); //use alternate string (expand recursively)
	  continue;
	}
	//	cptr=tail;
	if ((*buffer=='#') && (flags & SUB_ABORTHASH)) {
	  pushback(curptr,endptr,'$');
	  if (hasparen)
	    pushback(curptr,endptr,'(');
	  if (isstar)
	    pushback(curptr,endptr,'*');
	  pushback(curptr,endptr,dollar);
	  if (hasparen)
	    pushback(curptr,endptr,')');
	  continue;
	}
	parser_printcontext() << "undefined variable: $";
	if (n>=0)
	  std::cerr << n;
	else
	  std::cerr << dollar;
	std::cerr << '\n';
	error_abort();
      }
      if (haspositive) {
	int locflags=flags;
	if (curptr!=buffer)
	  locflags &= ~SUB_ABORTHASH;
	curptr=substitute_string(curptr,endptr-curptr,haspositive,locflags, accuses); //use alternate string (expand recursively)
	continue;
      }
      pushback(curptr,endptr,varval);
      if (isstar) { //also need to push additional arguments
	n++;
	while (n<=curframep->nargs()) {
	  pushback(curptr,endptr,' ');
	  pushback(curptr,endptr,curframep->arg(n++));
	}
      }
      // cptr=tail;
      continue;
    }

    bool isconst=true;
    int uses=0;
    const systemvarmap_type::const_iterator curp(systemvarmap.find(dollar));
    if (curp==systemvarmap.end()) {
      UserVariable* curv=findvariable(dollar);
      isconst=isvarconstant(curv);
      uses=getuses(curv);
      const BaseList<double>& value(curv->get_list());
      if (value.size()==1) {
	snprintf_prec(tmp,sizeof(tmp),value.front(),ignoreformat);
	pushback(curptr,endptr,tmp);
      }
      else {
	pushback(curptr,endptr,"[");
	const char* startlist=curptr;
	for (size_t i=0;i<value.size();i++) {
	  if (curptr-startlist>substitute_maxchars) {
	    pushback(curptr,endptr,",<list truncated>");
	    break;
	  }
	  if (i)
	    pushback(curptr,endptr,",");
	  snprintf_prec(tmp,sizeof(tmp),value(i),ignoreformat);
	  pushback(curptr,endptr,tmp);
	}
	pushback(curptr,endptr,"]");
      }
    }
    else {
      uses=getuses(curp->second);
      isconst=(isvarconstant(curp->second) || isallowed_nonconst_var(curp->second)); //!< make exception for evaluation index which is uniquely non-const overall, but still const outside execution loop
      pushback(curptr,endptr,curp->second->format());
    }
    accuses|=uses;
    if (!isconst && (flags & SUB_NONCONSTWARN))
      nonconstsubs_warning.raise(origcptr);

    cptr=tail;
  }
  *curptr='\0';// ensure terminated
  return curptr;
}

//ContextWarning<> unresolveddollar_warning("$ in text string is not being substituted",&NMRsim_repeat_warning);
ContextWarning<> replacedollarignored_warning("REPLACEDOLLAR flag is being ignored",&NMRsim_repeat_warning);

char* parse_string(int flags)
{
  return parse_string_syntax(NULL,0,flags);
}

char* parse_string_syntax(const char* syn, size_t narg, int flags)
{
  char* cptr=get_token(flags,"string");
  if (!cptr || !(*cptr)) {
    if (flags & F_ALLOWMISSING)
      return NMRSIM_NULL;
    error_abort("Null string argument");
  }
  else {
    if (!(flags & F_IGNORESYNTAX))
      verify_syntax(cptr,syn,narg);
  }

#ifndef NDEBUG
  if (flags & F_REPLACEDOLLAR)
    replacedollarignored_warning.raise();
#endif
//   if (!(flags & F_IGNOREDOLLAR) && ((flags & F_REPLACEDOLLAR) || !nochecks) && strchr(cptr,'$')) {
//     if (flags & F_REPLACEDOLLAR) {
//       char buffer[MAXLINE];
//       substitute_string(buffer,sizeof(buffer),cptr);
//       return strdup(buffer); //memory leak, but not cumulative
//     }	
//     unresolveddollar_warning.raise();
//   }
  if (!(flags & F_IGNOREDOLLAR) && strchr(cptr,'$')) {
    char buffer[MAXLINE];
    substitute_string(buffer,sizeof(buffer),cptr);
    return strdup(buffer); //memory leak, but not cumulative
  }
  return cptr;
}

bool parser_isnormal() 
{
  const char* ptr=get_curline();
  return (ptr && *ptr) ? ((*ptr!='-') || isspace(ptr[1]) || !ptr[1]) : false;
}

bool parser_isalpha()
{
  const char* ptr=get_curline();
  return (ptr && *ptr) ? isalpha(*ptr) : false;
}
  
size_t parse_flag(const flagsmap_type& map, const char* cptr)
{
  if (strcmp(cptr,"?")!=0) {
    const flagsmap_type::const_iterator curp(map.find(cptr));
    if (curp!=map.end()) 
      return curp->second;
  }
  parser_printcontext() << "Allowed flags are ";
  dumpflags(std::cerr,map) << '\n';
  error_abort();
  return 0; //!< dummy to prevent bogus warning
}

size_t parse_flags(const flagsmap_type& map, size_t defflags, bool oneonly)
{
  size_t flags=0;
  for (;;) {
    char* cptr= parser_isflag() ? get_token(F_ALLOWMISSING) : NMRSIM_NULL;
    if (cptr==NMRSIM_NULL)
      return flags ? flags : defflags;
    if (flags & oneonly) {
      parser_printcontext() << "Only one flag allowed out of ";
      dumpflags(std::cerr,map) << '\n';
      error_abort();
    }
    const size_t add=parse_flag(map,++cptr);
    if (flags & add)
      throw InternalError("parse_flags: overlapping flags (or oneonly not set)");
    flags|=add;
  }
}

std::ostream& dumpflags(std::ostream& ostr, const flagsmap_type& map, size_t flags)
{
  flagsmap_type::const_iterator start=map.begin();
  const flagsmap_type::const_iterator end=map.end();
  bool needspace=false;

  while (start!=end) {
    if (flags & (start->second)) {
      //      const std::string& fname=(start->first).name;
      if (needspace)
	ostr << ' ';
      else
	needspace=true;
      ostr << '-' << (start->first);
    }
    ++start;
  }
  return ostr;
}

void parser_checkfinished()
{
  if (are_left())
    error_abort("too many arguments");
}

template<> size_t parse_raw(char* cptr, Variable* vardefp, int flags)
{
  const int val=parse_raw<int>(cptr,vardefp,flags);
  if (val<0) {
    parser_printcontext() << cptr << " evaluated to negative integer\n";
    error_abort();
  }
  return val;
}

int parse_unsigned(Variable* vardefp, int flags)
{
  return parse_unsigned_syntax(NMRSIM_NULL,0,vardefp,flags);
}

int parse_unsigned_syntax(const char* syn, size_t narg, Variable* vardefp, int flags)
{
  char* cptr=get_token(flags,"positive integer");
  if (!cptr) {
    if ((flags & F_ALLOWMISSING)==0)
      throw InternalError("parse_unsigned_syntax");
    return -1;
  }
  verify_syntax(cptr,syn,narg);
  return parse_raw<size_t>(cptr,vardefp,flags);
}
  
size_t parse_unsigned_offset(int offset, int flags)
{
  const int val=parse_unsigned(NULL,flags);
  if (val<offset) {
    parser_printcontext() << "parsed number (" << val << ") is not allowed to be less than " << offset << '\n';
    error_abort();
  }
  return size_t(val-offset);
}

template<> int parse_raw(char* cptr, Variable* vardefp, int flags)
{
  return round_int(parse_raw<double>(cptr,vardefp,flags));
}

int parse_int(Variable* vardefp, int flags)
{
  return parse_raw<int>(get_token(0,"integer"),vardefp,flags);
}

LIST<int> parse_intarray(int flags)
{
  static LIST<double> tmparray;
  parse_array(tmparray,get_token(0,"integer list"),flags);
  LIST<int> destarray(tmparray.size(),mxflag::temporary);
  for (size_t i=tmparray.size();i--;)
    destarray(i)=round_int(tmparray(i));
  return destarray;
}

//! none of this is very efficient - but doesn't need to be
std::ostream& print_syntax_string(std::ostream& ostr, const char* str, size_t n)
{
  const size_t orign=n;
  if (n<1)
    throw InvalidParameter("print_syntax_string");
  const char* curp=str;
  while (n!=1) {
    curp=strchr(curp,'#');
    if (!curp) {
      std::cerr << "Can't find argument " << orign << " in syntax string: " << str << '\n';
      error_abort();
    }
    curp++;
    n--;
  }
  while (*curp) {
    ostr.put( (*curp=='#') ? ' ' : *curp);
    curp++;
  }
  return ostr;
}

void verify_syntax(const char* arg, const char* syn, size_t n)
{
  if (!arg)
    throw InvalidParameter("verify_syntax: NULL argument");
  if (strcmp(arg,"?")!=0)
    return;
  if (syn==NMRSIM_NULL)
    error_abort("Sorry - don't have syntax for this argument!");

  if (n==0)
    parser_printcontext() << "Syntax: " << syn << '\n';
  else {
    parser_printcontext() << "Syntax of remaining argument(s): ";
    print_syntax_string(std::cout,syn,n);
    std::cout << '\n';
  }
  error_abort();
}
  
LIST<size_t> parse_unsignedintarray_syntax(const char* synstr, size_t n, int offset, int flags)
{
  static LIST<double> tmparray;
  char* tok=get_token(0,"positive integer list");
  verify_syntax(tok,synstr,n);
  parse_array(tmparray,tok,flags);
  LIST<size_t> destarray(tmparray.size(),mxflag::temporary);
  for (size_t i=tmparray.size();i--;) {
    const int val=round_int(tmparray(i)-offset);
    if (val<0) {
      parser_printcontext() << "list contains elements less than " << offset << '\n';
      error_abort();
    }
    destarray(i)=val;
  }
  return destarray;
}

LIST<size_t> parse_unsignedintarray(int offset, int flags)
{
  return parse_unsignedintarray_syntax(NULL,0,offset,flags);
}

Variable::Variable(subsid_t subsidv, Setable* ptrv, VariableBase* valuepv)
  : subsid(subsidv), valuep(valuepv)
{
  if (ptrv)
    setptr(ptrv);
  else
    ptr=NMRSIM_NULL;
}

const char* Variable::name() const
{
  if (namestr.empty()) {
    validate();
    std::ostringstream str(std::ostringstream::out);
    ptr->printvariablename(str,subsid);
    namestr=str.str();
    if (namestr.empty())
      throw InternalError("Variable name cannot be empty");
  }
  return namestr.c_str();
}

void Variable::setptr(Setable* ptrv)
{
  if (!ptrv)
    throw InvalidParameter("Variable: pointer not set");
  ptr=ptrv;
  if (valuep) {
    if ((verbose & VER_PARSE) && (verbose_level>1)) {
      std::cout << "Creating variable map entry: " << key() << "\n";
      std::cout << (*this) << '\n';
    }
    varmap[key()]=valuep;
  }
}

 void Variable::validate() const
 {
   if (!ptr || !valuep)
     throw InternalError("Variable: pointers not set");
 }

varkey_t Variable::key() const
{ 
  if (!ptr)
    throw InternalError("Variable::key: Setable pointer not set");
  return varkey_t(ptr,subsid);
}

ThreadWarning<> UserVariable::unlinkedvar_warning("Non-const variable defined but never used - almost certainly an error: ",&NMRsim_repeat_warning);
ThreadWarning<> UserVariable::unlinkedconstvar_warning("Variable defined but not directly used: ",&NMRsim_repeat_warning);

void UserVariable::check_unused()
{
  uservarmap_type::const_iterator curp(vardefmap.begin());
  const uservarmap_type::const_iterator end(vardefmap.end());
  while (curp!=end) {
    const UserVariable& cvar(*((curp->second).get()));
    if (!(cvar.isused()) && !(cvar.ignore_unused())) {
      if (cvar.isconstant())
	unlinkedconstvar_warning.raise(cvar.name());
      else
	unlinkedvar_warning.raise(cvar.name());
    }
    ++curp;
  }
}

double UserVariable::get() const
{ 
  const BaseList<double>& val(value().value());
  if (val.size()!=1)
    throw Failed("Array cannot be used in this context");
  return val.front();
}

void UserVariable::set(double fval, subsid_t subsid)
{
  assert(subsid==S_ARG1); //!< harmless assert
  value().set(fval);
}

void UserVariable::set(const BaseList<double>& fvals, subsid_t subsid)
{
  assert(subsid==S_ARG1); //!< harmless assert
  //   BaseList<double> dest(valuep->value());
  if (fvals.vector()!=value().value().vector()) //!< trap assign to self
    value().set(fvals);
}

// void UserVariable::print(std::ostream& ostr, bool full) const
// {
//   print(ostr,full ? S_NONE : S_ARG1);
// }

void UserVariable::print(std::ostream& ostr) const
{
  ostr << name() << " = ";
  value().print(ostr,true);
  if (isconstant() || (!isused())) {
    ostr << " (";
    if (isconstant())
      ostr << "const" << (isused() ? ")" : ", unused)");
    else
      ostr << "unused)";
  }
  if (RealVariable::uses()) {
    ostr << ' ';
    printattributes(ostr,RealVariable::uses());
  }
}

void UserVariable::printvariablename(std::ostream& ostr, subsid_t) const
{
  ostr << name();
}

template<class T> inline void assign_(BlockedMatrix<T>& to, const BlockedMatrix<T>& from)
{
  //catch assignments that would change matrix size - this must be wrong!
  if ((!!to) && !arematching(to,from))
    throw Failed("Can't assign between matrices of differing size");
  to=from;
}

UserVariable* isvardefined(const char* varchar)
{
  //  const bool islocal=(*varchar=='_');
  //  const uservarmap_type& usemap = islocal ? curframep->localvarmap : vardefmap;
  if (systemvarmap.find(varchar)!=systemvarmap.end())
    error_abort("Can't redefine system variable");
  const uservarmap_type& usemap(vardefmap);
  const uservarmap_type::const_iterator iter=usemap.find(varchar);
  return (iter==usemap.end()) ? NMRSIM_NULL : iter->second.get();
}

UserVariable* makevariable(const char* varchar, VariableBase& var, int flagsv =0)
{
  UserVariable* newvarp = new UserVariable(varchar,var,flagsv);
  vardefmap[varchar].reset(newvarp);
  return newvarp;
}

UserVariable* UserVariable::create(const char* varchar, double val, int userflagsv)
{
  VariableBase* varp=new VariableBase(val);
  return create(varchar,*varp,userflagsv);
}

const char* parse_variable_name()
{
  return clean_variable_name(parse_string());
}

//ContextWarning<> UserVariable::replacingusedwithvar_warning("Replacing already used constant variable with non-const version - this is valid, but not recommended: ",&NMRsim_repeat_warning);
ContextWarning<> UserVariable::redefininguservariable_warning("Replacing already used variable with another definition. Add -ignore_refinition flag to original definition to variable definition if this is intended",&NMRsim_repeat_warning);
ContextWarning<> UserVariable::replacingunused_warning("Replacing unused (redundant?) definition: ",&NMRsim_once_warning);

UserVariable* UserVariable::create(const char* varchar, VariableBase& var, int userflagsv)
{
  //  const bool islocal=(*varchar=='_');
  //  uservarmap_type& usemap = islocal ? curframep->localvarmap : vardefmap;
  //const uservarmap_type::iterator iter=usemap.find(varchar);
  //  const bool found=(iter!=usemap.end());
  UserVariable* newvarp=isvardefined(varchar);
  if (!newvarp)
    return makevariable(varchar,var,userflagsv);
  if (newvarp->userflags() != userflagsv)
    throw InternalError("UserVariable::create - can't change variable flags");
  if (newvarp->isused()) {
    if (newvarp->isconstant()) {
      if (!(newvarp->userflags() & UserVariable::IGNORE_REDEFINITION))
	redefininguservariable_warning.raise(varchar);
      //if (!(var.isconst()))
      //replacingusedwithvar_warning.raise(varchar);
    }
    else {
      parser_printcontext() << "can't redefine variable that is in use: " << varchar << '\n';
      error_abort();
    }
  }
  else
    replacingunused_warning.raise(varchar);
    
  newvarp->reset(var);//,isconst);
  return newvarp;
}

size_t parse_counting(const char* str)
{
  char *end;
  size_t n=strtoul(str,&end,10);
  if ((errno==ERANGE) || (str==end))
    n=0;
  return n;
}

size_t maxslot(const char* str)
{
  size_t maxn=0;
  for (;;) {
    str=strchr(str,'#');
    if (str==NMRSIM_NULL)
      return maxn;
    char* end;
    str++;
    const long n=strtol(str,&end,10);
    if (end!=str) {//!< ignore invalid - if real problem will be caught when expression is parsed
      if (n<1)
	error_abort("Arguments (#<n>) are numbered from 1");
      if (n>maxn)
	maxn=n;
      str=end;
    }
  }
}

Expression* parse_expression(const char* cptr)
{
  expr_grammar& expr_parser(get_expr_parser());
  expr_parser.reset(F_ALLOWSLOT);
  const parse_info<> info(parse(cptr,
			  expr_parser.use_parser<expr_grammar::expression_def>(),
			  space_p));
  if (!info.full) {
    parser_printcontext() << "expression parsing failed at: \"" << info.stop << "\"\n";
    error_abort();
  }
  Expression& expr(expr_parser.expression());
  if (expr.empty())
    throw InternalError("parse_expression");
  Expression* newexprp=new Expression;
  newexprp->swap(expr);
  return newexprp;
}

namespace {
  //! check whether function (or chained) can be redefined
  void checkredefine(const function_spec& spec)   
  {
    Function_Factory_t& factory(get_Function_Factory());    
    std::string curname(spec.first);

    for (;;curname+=NMRSIM_OLD_FUNCTION) {
      const function_spec curspec(curname.c_str(),spec.second);
      Function_Factory_t::iterator iter(factory.find(curspec));
      if (iter==factory.end()) 
	return; //!< not found
    
      const function_def_t& fdef(iter->second);
      if (fdef.used) {
	parser_printcontext() << "can't redefine function previously used (as non const): " << curname << '\n';
	error_abort();
      }    
      function_redefinition_warning.raise(curname.c_str());
    }
  }

  //! do redefinition
  void doredefine(const function_spec& spec, ExpressionUserFunction* funcp)
  {
    function_def_t fdef(funcp);
    Function_Factory_t& factory(get_Function_Factory());    
    std::string curname(spec.first);

    for (;;curname+=NMRSIM_OLD_FUNCTION) {
      const function_spec curspec(curname.c_str(),spec.second);
      Function_Factory_t::iterator iter(factory.find(curspec));
      if (iter==factory.end()) {
	std::string* cnamep=new std::string(curname);	
	factory[function_spec(cnamep->c_str(),spec.second)]=fdef;
	return;
      }
      ::std::swap(fdef,iter->second);
    }
  }
}

ContextWarning<> slotmismatch_warning("mismatch between number of argument types specified and number of arguments used in function: ",&NMRsim_once_warning);

LIST<const char*> argnamestack;

//! helper class to make sure that argnamestack only contains entries when argument names exist
struct argnamestack_guard {
  argnamestack_guard() { argnamestack.clear(); }
  ~argnamestack_guard() { argnamestack.clear(); }
};

void parse_function()
{
  argnamestack_guard argnamestack_guard;
  static LIST<const Expression*> expr_pointer_stack;

  char* funcnamep=parse_string();

  char* startargsp=strchr(funcnamep,'(');
  if (startargsp) {
    *startargsp++='\0';
    char* endargsp=strchr(startargsp,')');
    if ((endargsp==NMRSIM_NULL) || (endargsp[1]!='\0'))
      error_abort("syntax of function arguments is <function>(<args>) - failed to find matching ()");
    *endargsp='\0';
    char* strtokarg=startargsp;
    for (;;) {
      char* arg=strtok(strtokarg,",");
      if (arg==NULL)
	break;
      strtokarg=NMRSIM_NULL;
      check_variable_name(arg);
      if ((strcmp(arg,"L")==0) || (strcmp(arg,"S")==0))
	error_abort("L and S are reserved and can't be used for function argument names");
      for (size_t i=argnamestack.size();i--;) {
	if (strcmp(argnamestack(i),arg)==0)
	  error_abort("argument names need to be distinct");
      }
      argnamestack.push_back(arg);      
    }
    if (argnamestack.empty())
      error_abort("failed to find any function arguments");
    if ( (verbose & VER_PARSE) && (verbose_level>1))
      std::cout << "Found " << argnamestack.size() << " named arguments for function " << funcnamep << '\n';
  }

  const std::string* strnamep=new std::string(funcnamep); //!< need to create new string rather than pointer
  const char* funcname=strnamep->c_str();
  
  const char* exprval=parse_string();

  argdef_t argtype=0;
  size_t argcount=0;
  if (are_left()) {
    const char* argdef=parse_string();
    argcount=strlen(argdef);
    if (!(argnamestack.empty()) && (argnamestack.size()!=argcount))
      error_abort("mismatch between number of function argument names and number of argument type specifiers");

    for (size_t i=argcount;i--;) {
      argtype<<=1;
      switch (argdef[i]) {
      case 'S':
	break;
      case 'L':
	argtype|=1;
	break;
      default:
	parser_printcontext() << "invalid argument type specification (" << argdef << ").  Each argument should be one of S (scalar) or L (list)\n";
	error_abort();
      }
    }
    if (slotmismatch_warning.enabled()) {
      if (maxslot(exprval)!=argcount)
	slotmismatch_warning.raise(funcname,true);
    }
  }
  else {
    if (argnamestack.empty())
      argcount=maxslot(exprval);
    else
      argcount=argnamestack.size();
//     if (argcount==0)
//       nofunctionarguments_warning.raise(funcname,true);
  }
  const function_spec spec(funcname,argcount);
  checkredefine(spec);

  expr_pointer_stack.push_back((const Expression*)NMRSIM_NULL); //!< create space for new indirect pointer
  ExpressionUserFunction* funcp=new ExpressionUserFunction(spec,argtype,expr_pointer_stack.back());
  doredefine(spec,funcp);
  const Expression& expr(*parse_expression(exprval));
//   if (expr.isconstant())
//     constfunction_warning.raise(funcname,true);
  funcp->setexpression(expr);
  if (verbose & VER_PARSE) {
    std::cout << "Creating function " << funcname << " with";
    if (funcp->nargs()) {
      std::cout << " arguments ";
      for (size_t i=0;i<funcp->nargs();i++)
	std::cout << (funcp->islist(i) ? 'L' : 'S');
      std::cout << ": ";
      expr.print();
      std::cout << '\n';
    }
    else
      std::cout << "out arguments\n";
  }
}

void check_variable_name(const char* p)
{
  if (!isalpha(*p++))
    error_abort("Variable/argument names must start with letter");
  for (;*p;p++) {
    if (!isalpha(*p) && !isdigit(*p) && (*p!='_'))
      error_abort("Variable/argument names can only contain letters, digits or _");
  }  
}

void parse_variable()
{
  LIST<const char*> varnames;
  char* varcharall=parse_string();
  Splitter splitter(varcharall,",;",true);
  char* varchar;
  char tupleseparator;
  while ((varchar=splitter.next(tupleseparator))) {
//     const bool islocal=(*varchar=='_');
//     if (islocal) {
//       if ((curframep==topframep) && !nochecks)
// 	parser_printcontext() << "warning: local variable declared at top level and so effectively global (block level scoping intended?)\n-nochecks disables this warning\n";
//     }
//     else {
    check_variable_name(varchar);
    const systemvarmap_type::iterator iter(systemvarmap.find(varchar));
    if (iter!=systemvarmap.end()) {
//       const SystemVariableBase& var(*(iter->second));
//       std::cerr << "Hash: " << iter->first << "  ";
//       var.print(std::cerr); std::cerr << '\n';
//       std::cerr << "Hash of " << varchar << ": " << systemvarmap.hash(varchar) << std::endl;
      error_abort("Cannot redefine existing (system) variable");
    }
    if (isenvironment(varchar))
      error_abort("Variable name must contain lower case letters (to prevent clash with environment variables)");
    
    if ((strcmp(varchar,"x")==0) || (strcmp(varchar,"y")==0))
      xy_variable_warning.raise();

    varnames.push_back(varchar);
  }
  if (varnames.empty())
    error_abort("no variables to define!");
        
  char* cptr=get_token();
  if (checkbracket(cptr,'\'')) {
    if (varnames.size()!=1)
      error_abort("strings can't be used in tuple assignment");    
    const char* varchar(varnames.front());
    if (*varchar=='_')
      error_abort("string variables cannot be local");
    systemvarmap[varchar]=new SystemVariable<std::string>(varchar,cptr);
    return;
  }

  Mark markobj;
  const bool istuple=(varnames.size()>1);
  UserVariable* oldvarp=istuple ? NMRSIM_NULL : isvardefined(varnames.front());
  //Postpone until UserVariable::create ?
//   if (oldvarp && oldvarp->islinked())
//     error_abort("can't redefine variable that is used");
  Variable cvar(S_ARG1,oldvarp);
  parse_double_raw(&cvar,cptr,F_ALLOWLIST | F_ALWAYSNEW); //!< get value
  const bool isconst=(cvar.subsid==S_NONE); //!< will be flagged as non-const even if only one element of list is non-const

  static flagsmap_type vflags;
  if (vflags.empty()) {
    vflags["ignore_unused"]=UserVariable::IGNORE_UNUSED;
    vflags["ignore_refinition"]=UserVariable::IGNORE_REDEFINITION;
  }

  const int userflags=parse_flags(vflags);

  VariableBase* valuep=cvar.valuep;
  if (valuep==NMRSIM_NULL)
    throw InternalError("parse_variable");
  LIST<double> tmpvars;
  const size_t nvalues=valuep->value().size();
  if (istuple) {
    if (isconst) {
      if ((tupleseparator!=';') && (nvalues % varnames.size())) {
	parser_printcontext() << "bad tuple assignment: " << (valuep->value().size()) << " values cannot be split evenly between " << varnames.size() << " variables\n";
	error_abort();
      }
//       if (varnames.size()>nvalues) {
// 	parser_printcontext() << "bad tuple assignment: fewer values (" << nvalues << ") than destination variables (" << varnames.size() << '\n';
// 	error_abort();
//       }
      tmpvars.create((nvalues+varnames.size()-1) / varnames.size());
    }
    else {
      if (!(valuep->isexpr()))
	error_abort("Can only use tuple assignment for constants / expressions");
    }
  }
  static ScratchList<size_t,2> extargs(3U);
  static size_t tuplecount=0;
  UserVariable* tuplevarp=NMRSIM_NULL;
  if (istuple && !isconst) {//!< create base tuple expression
    tuplecount++;
    char tmpname[16];
    sprintf(tmpname,"<tuple%" LCM_PRI_SIZE_T_MODIFIER "u>",tuplecount);
    tuplevarp=UserVariable::create(tmpname,*valuep,UserVariable::IGNORE_UNUSED); //!< internal variable
    tuplevarp->flagused(); //!< always used
    markobj.flush(tuplevarp);
  }
  size_t valuesleft=nvalues;
  const bool commasep=(tupleseparator==',');
  for (size_t j=0;j<varnames.size();j++) {
    VariableBase* usevarp=valuep;
    const char* varchar(varnames(j));
    if (istuple) {
      if (isconst) {
	const BaseList<double> source(valuep->value());
	const size_t n=tmpvars.size();
	if ((n==1) && commasep)
	  usevarp=new VariableBase(source(j));
	else {
	  if (n) {
	    if (commasep)
	      tmpvars=source(libcmatrix::range(n*j,n*j+n-1));
	    else {
	      const size_t items=(valuesleft<n) ? valuesleft : n;
	      valuesleft-=items;
	      if (items) {
		const slice sl(j,items,varnames.size());
		tmpvars=source(sl);
	      }
	      else
		tmpvars.clear();
	    }
	  }
	  usevarp=new VariableBase(tmpvars);
	}
      }
      else { //!< can assume that valuep refers to an expression	
	Expression* newexpr = new Expression;
	extargs(0U)=newexpr->push(new ExpressionVariable(*tuplevarp));
// 	extargs(1U)=newexpr->push(new ExpressionConstant(j+1));
// 	newexpr->push(new FunctionExtract(extargs));
	extargs(1U)=newexpr->push(new ExpressionConstant(varnames.size()));
	extargs(2U)=newexpr->push(new ExpressionConstant(j+1));
	if (tupleseparator==',')
	  newexpr->push(new FunctionSplitRow(extargs));
	else
	  newexpr->push(new FunctionSplitColumn(extargs));
	usevarp=new VariableBase;
	usevarp->initialise(*newexpr);
	Variable var(S_ARG1,isvardefined(varchar),usevarp);
	register_expression(var);
      }
    }
    markobj.flush(UserVariable::create(varchar,*usevarp,userflags));
  }
}

std::ostream& parser_printcontext(std::ostream& ostr, bool printtoken)
{
  parser_printthread(ostr);
  const bool havetoken=((lasttoken!=NMRSIM_NULL) && printtoken);
  if (curframep)
    ostr << "At " << curframep->state();
  if (havetoken)
    ostr << (curframep ? ", parsing " : "While parsing ") << lasttoken;
  if (havetoken || curframep)
    ostr << ": ";
  return ostr;
}

static char formatstr[10]="";
static bool hasprec=false;

std::set< std::string> envlist;

const char* register_and_getenv(const char* varname)
{
  const std::string varnamestr(varname);
  envlist.insert(varnamestr);
  return getenv(varname);
}

bool isregistered(const char* name)
{
  return (envlist.find(std::string(name))!=envlist.end());
}
    
void parse_nmrsim_format(const char* nmrsimformat)
{
  const bool useexplicit=(nmrsimformat!=NMRSIM_NULL) && (*nmrsimformat!='\0');
  if (useexplicit) {
    formatstr[0]='%';
    if (strchr(nmrsimformat,'%'))
      error_abort("NMRSIM_FORMAT should not contain %.  Examples of valid usage: g, f, .*g etc.");
    strncat(formatstr+1,nmrsimformat,sizeof(formatstr)-2);
  }
  else
    strcpy(formatstr,"%.*g");
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Setting numeric formatting to " << formatstr << " (";
    if (useexplicit)
      std::cout << "from NMRSIM_FORMAT)\n";
    else
      std::cout << "default)\n";
  }
  hasprec=(strchr(formatstr,'*')!=NMRSIM_NULL);
}

void snprintf_prec(char* buf, size_t s, double val, bool ignoreformat)
{
  if (ignoreformat) {
    snprintf(buf,s,"%.15g",val);
    return;
  }
  if (*formatstr=='\0') {
    parse_nmrsim_format(register_and_getenv("NMRSIM_FORMAT"));
    if (*formatstr=='\0')
      throw InternalError("snprintf_prec");
  }
  if (hasprec) {
    const ostream_controller& ctrl(cmatrix_ostream_controller());
    const int precision=(ctrl.matrixprecision<0) ? std::cout.precision() : ctrl.matrixprecision;
    snprintf(buf,s,formatstr,precision,val);
  }
  else
    snprintf(buf,s,formatstr,val);
}

void parse_ignored()
{
  set_curline(NMRSIM_NULL);
  ignoreddirective_warning.raise();
}


double Variable::get() const
{
  const BaseList<double>& value(valuep->value());
  if (value.size()!=1)
    throw Failed("Can't return list expression as single value");
  return value.front();
}

const char* Variable::format() const
{
  const BaseList<double>& value(valuep->value());
  if (value.size()>1)
    strcpy(buf,"<list>");
  else
    snprintf_prec(buf,sizeof(buf),value.front());
  return buf;
}

void Variable::printname(std::ostream& ostr) const
{
  if (!ptr)
    throw InternalError("Variable: ptr unset");  
  ptr->printvariablename(ostr,subsid);
}

void Variable::print(std::ostream& ostr, bool full) const
{
  printname(ostr);
  ostr << ": ";
  valuep->print(ostr,full);
}

// namespace {
//   template<class T> T* unwrap(T* v) { return v; }
//   template<class T> T* unwrap(const smartptr<T,false>& v) { return v.get(); }
// }

template<class T> void dumpvariables(const T& varmap, const char* type)
{
  parser_printthread(std::cerr) << "Currently defined " << type << " variables:";
  typename T::const_iterator start(varmap.begin());
  const typename T::const_iterator end(varmap.end());
  if (start==end)
    std::cerr << " <none>";
  else {
    while (start!=end) {
      std::cerr << ' ' << start->second->name();
      ++start;
    }
  }
  std::cerr << std::endl;
}

UserVariable* findvariable(const char* name, bool flagused)
{
  const uservarmap_type::const_iterator curp(vardefmap.find(name));
  const uservarmap_type::const_iterator end(vardefmap.end());
  if (curp==end) {
    if (isenvironment(name)) {
      parser_printcontext() << "unsubstituted environment variable: " << name << "\nUse $(" << name << ") to separate substitution variable from following text\n";
      error_abort();
    }
    parser_printcontext() << "variable " << name << " not defined\n";
    if (noexecute)
      std::cerr << "(NB. user-defined variables are not created with -noexecute)\n";
    else {
      dumpvariables(vardefmap,"user-defined");
      dumpvariables(systemvarmap,"system");
    }
    error_abort();
  }
  UserVariable* retval=(curp->second).get();
  if (flagused)
    retval->flagused();
  return retval;
}

UserVariable* findnewvariable(const char* name)
{
  const uservarmap_type::const_iterator curp(vardefmap.find(name));
  const uservarmap_type::const_iterator end(vardefmap.end());
  if (curp!=end) {
    parser_printcontext() << "variable " << name << " already in use - must use a fresh variable name in this context\n";
    error_abort();
  }
  if (isenvironment(name))
    error_abort("variable name must contain lowercase characters to distinguish from environment variable");
  VariableBase* varp=new VariableBase();
  return makevariable(name,*varp);
}

void process_array_ns()
{
  nacqdims=1+array_dims.get(array_ns,array_skips,array_n0); //extract dimensionality
  ndims=1+array_ns.size();
  //   if ((ndims==1) && (array_n0>1))
  //     ndims++;
  if (verbose & VER_GEN) {
    std::cout << "Data set dimensions: " << ndims;
    std::cout << "\nAcquisition dimensions: " << nacqdims << '\n';
    std::cout << "Total number of rows: " << array_n0 << '\n';
  }
}

dimension_set::dimension_set(size_t n)
  : block_change(true)
{
  dims[0].set(n,true);
}

 void dimension_set::set(size_t dim, size_t n)
{
  if (dim>MAX_DIMENSIONS) {
    parser_printcontext() << "dimension (" << dim << ") out of range.  Maximum is " << MAX_DIMENSIONS << '\n';
    error_abort();
  }
  dims[dim].set(n,block_change);
}

ContextWarning<> rowsskip_warning("total number of rows is not multiple of 'skip factor'",&NMRsim_repeat_warning);

void dimension_set::set(size_t dim, size_t n, size_t skip)
{
  if (dim>MAX_DIMENSIONS) {
    parser_printcontext() << "dimension (" << dim << ") out of range.  Maximum is " << MAX_DIMENSIONS << '\n';
    error_abort();
  }
  if (dim==0)
    throw InvalidParameter("Can't set indirect dimension 0");
  if (block_change)
    error_abort("cannot set indirect dimension at this point (already fixed)");
  if (n<1)
    error_abort("dimension cannot be <1");
  if (skip<1)
    error_abort("row step cannot be <1");
    //    skips(dim-1)=ni_skip;
  if (n % skip)
    rowsskip_warning.raise();
  dims[dim].set(n,skip);
}

void dimension_set::dimension::set(size_t n, size_t skipv)
 {
   if (n==0)
     throw InvalidParameter("Dimension cannot be zero");
   if (hard)
     error_abort("cannot change size of indirect dimension");
   hard=n;
   skip=skipv;
 }

void dimension_set::dimension::set(size_t n, bool block)
 {
   if (n==0)
     throw InvalidParameter("Dimension cannot be zero");
   
   if ((n==dim) || (n==1))
     return; //!< no change

   if (dim && ((n<dim) ? (dim % n) : (n % dim)))
     incommensurate_warning2.raise(state);
   
   if (n>=dim) {
     if (block) {
       parser_printcontext() << "cannot change dimension (from " << dim << " to " << n << ") at this point. Sum, {} variables are best declared in par block\n";
       error_abort();
     }
     dim=n;
     state=curframep->state();
   }
 }

void dimension_set::dimension::error(size_t dim, const char* err)
{
  parser_printcontext() << "dimension " << dim << ": " << err << '\n';
  error_abort();
}

size_t dimension_set::dimension::get(size_t& array_skip, size_t N) const
{
  array_skip=1;
  if (hard) {
    if (dim) { //!< also soft set?
      if ((hard+skip-1)/skip==dim) //!< need skip to make things match?
	array_skip=skip;
      else {
	if (dim!=hard)
	  error(N,"array and nD size are incommensurate");
      }
    }
    else //!< by default follow any "quadrature" setting
      array_skip=skip;
    return hard;
  }
  if (!dim)
    error(N,"not set");
  return dim;
}

size_t dimension_set::get(LIST<size_t>& dest, LIST<size_t>& array_skip, size_t& totdim) const
{  
  size_t maxdim=MAX_DIMENSIONS;
  while (maxdim && !dims[maxdim].isset())
    maxdim--;
  totdim=1;
  dest.create(maxdim);
  array_skip.create(maxdim);
  size_t nhard=0;
  for (size_t dim=maxdim;dim;dim--) { //check for "holes"
    size_t& curskip(array_skip(dim-1));
    size_t n=dims[dim].get(curskip,dim);
    if (dims[dim].ishard())
      nhard++;
    totdim*=n;
    dest(dim-1)=n;
  }
  const size_t dim0=dims[0].dim;
  if (dim0) {
    if (maxdim) {
      bool incom;
      if (dim0>totdim) {
	if (active2D) {
	  parser_printthread(std::cerr) << "Data set size smaller than number of rows set at " << dims[0].state << '\n';
	  error_abort();
	}
	incom=dim0 % totdim;
	totdim=dim0;
      }
      else
	incom=totdim % dim0;
      
      if (incom)
	incommensurate_warning.raise(dims[0].state);
    }
    else //not nD, but arrayed
      totdim=dim0;
  }
  return nhard;
}

void stripleaf(char* fname,const char* strip)
{
  char* poss_start=fname+strlen(fname)-strlen(strip);
  if ( (poss_start>=fname) && (strcmp(poss_start,strip)==0))
    *poss_start='\0';
}

bool read_block(const char* name,command_Factory_t& fac, bool isopt, bool ignoreunk)
{
  parser_newblock(name,isopt);
  bool foundblock=false;
  while (parser_getline()) {
    foundblock=true;
    const char* curcom(parse_string());
    command_Factory_t::iterator curp=fac.find(curcom);
    if (curp==fac.end()) {
      if (ignoreunk) {
	if (verbose & VER_PARSE)
	  parser_printcontext() << "ignoring " << curcom << '\n';
	continue;
      }
      else {	
	if (*curcom=='#')
	  error_abort("directive cannot begin with # - misplaced comment? (# must be first character of line)");
	if (!check_unrecognised(curcom)) {
	  std::cerr << "Available options are: ";
	  dump_map(fac,std::cerr);
	}
	error_abort();
      }
    }
    try {
      (curp->second)();
      parser_checkfinished();
    } catch (MatrixException& exc) {
      parser_printthread(std::cerr) << exc;
      std::cerr.flush();
      error_abort();
    }
    checkdirty();
  }  
  return foundblock;
}

template<class T> struct parser_type_traits {};
template<> struct parser_type_traits<int> {
  static const char* description;
};
template<> struct parser_type_traits<size_t> {
  static const char* description;
};
template<> struct parser_type_traits<double> {
  static const char* description;
};
const char* parser_type_traits<int>::description="integer";
const char* parser_type_traits<size_t>::description="positive integer / index";
const char* parser_type_traits<double>::description="floating point number";

template<class T> void parse_system_variable_syntax(const char* synstr, size_t narg, SystemVariable<T*>& gvar,int flags)
{
  Mark markobj;
  Variable cvar(S_ARG1,&gvar);
  char* tok=get_token(0,parser_type_traits<T>::description);
  verify_syntax(tok,synstr,narg);
  gvar=parse_raw<T>(tok,&cvar,flags);
  const bool isconst=(cvar.subsid==S_NONE);    
  if (isconst)
    gvar.isconstant(true); //! ensure const status
  else {
    gvar.reset(cvar.variable());
    markobj.flush(&gvar);
  }
}

template<class T> void parse_system_variable(SystemVariable<T*>& gvar,int flags)
{
  parse_system_variable_syntax(NULL,0,gvar,flags);
}

template void parse_system_variable(SystemVariable<double*>&, int);
template void parse_system_variable(SystemVariable<size_t*>&, int);
template void parse_system_variable(SystemVariable<int*>&, int);
template void parse_system_variable_syntax(const char*, size_t, SystemVariable<double*>&, int);
template void parse_system_variable_syntax(const char*, size_t, SystemVariable<size_t*>&, int);
template void parse_system_variable_syntax(const char*, size_t, SystemVariable<int*>&, int);

//! round \c double to nearest integer
/** exit with ::ERR_FAILED if not within ::NMRSIM_ROUNDTOL */
int round_int(double f)
{
  const int i=rawround_int(f);
  if (fabs(f-i)>NMRSIM_ROUNDTOL) {
    parser_printcontext() << "quantity (" << f << ") does not round to integer\n";
    error_abort();
  }
  return i;  
}

//NB silently fails if str is not a valid key!
bool havekey(const command_Factory_t& fac, const std::string& str)
{
  const command_Factory_t::const_iterator pos(fac.find(str));
  if (pos==fac.end()) {
    parser_printcontext() << "havekey: key doesn't exist: " << str << '\n';
    error_abort();
  }
  return pos->second.done; 
}

void verify_arraysize(size_t n, size_t dim)
{
  if (n!=1) {
    if (dim)
      have_virtualdimensions=true;
    array_dims.set(dim,n);
  }
}

void ExpressionVariableFunction::create(BaseList<std::string>& args)
{
  if (args.size()!=1)
    throw InternalError("Functions operating on variables only take one argument");
  varp=findvariable(args.front().c_str());
  const VariableBase rawvar(varp->value());
  if (rawvar.isexpr())
    appliedtoexpr_warning.raise();
}
  
//ContextWarning<> ExpressionVariableFunction::outofmainloop_warning("variable contents not well defined in this context: ",&NMRsim_once_warning);
//ContextWarning<> ExpressionVariableFunction::outofmainloop_expression_warning("expression may not be well defined in this context: ",&NMRsim_once_warning);
ContextWarning<> ExpressionVariableFunction::outofmainloop_warning("variable depends on {} array and is not well defined in this context: ",&NMRsim_once_warning);
ContextWarning<> ExpressionVariableFunction::outofmainloop_expression_warning("variable depends on || summation array and is not well defined in this context: ",&NMRsim_once_warning);

void ExpressionVariableFunction::verify_initialised() const
{
  if (varp==NMRSIM_NULL)
    throw InternalError("variable function not initialised");
}

void FunctionValueof::operator()(LIST<double>& dest) const
{
  verify_initialised();
  dest=varp->get_list();
}

bool FunctionValueof::rawisconstant() const
{
  verify_initialised();
  return varp->isconstant();
}

void FunctionValuesof::operator()(LIST<double>& dest) const
{
  verify_initialised();
  const VariableBase val(varp->value());
  if (val.isexpr())
    throw Failed("Valuesof applied to expression");
  if (val.isarray()) {
    //    if (val.isexprarray())
    //  throw Failed("Valuesof applied to non-const array");
    static const int deny_attr = A_ARRAY | A_SUM;
    if (val.arrayuses() & deny_attr)
      error_abort("Valuesof applied to array which can only be evaluated inside calculation loop");
    dest=val.get_array().row();
    if ((verbose & VER_GEN) && (verbose_level>1))
      std::cout << "Evaluated Valuesof to: " << dest << '\n';
  }
  else
    dest=val.value();
}

//ContextWarning<> ExpressionVariableFunction::nonarray_warning("Valuesof applied to non-array variable: ",&NMRsim_once_warning);

bool FunctionValuesof::rawisconstant() const
{
//   assert(varp!=NMRSIM_NULL);
//   if (!(varp->value().isarray()) || varp->value().isexprarray())
//     nonarray_warning.raise(varp->name());
  return true;
}
 
void FunctionErrorof::operator()(LIST<double>& dest) const
{
  if (varp==NMRSIM_NULL)
    throw InternalError("FunctionErrorof");
  double err=varp->get_error();
  if (err<0.0)
    err=0.0;
  dest.create(1U,err);
}

bool FunctionErrorof::rawisconstant() const
{
  if (varp==NMRSIM_NULL)
    throw InternalError("FunctionErrorof::rawisconstant");
  return (varp->value().issimple());
} 
bool FunctionErrorsof::rawisconstant() const { return true; }

void FunctionErrorsof::operator()(LIST<double>& dest) const
{
  if (varp==NMRSIM_NULL)
    throw InternalError("FunctionErrorsof");
  dest=varp->get_errors();
  for (size_t i=dest.size();i--;) {
    if (dest(i)<0.0)
      dest(i)=0.0;
  }
}

void ExpressionVariableFunction::rawprint(std::ostream& ostr) const
{
  ostr << name_;
  if (varp!=NMRSIM_NULL)
    ostr << "(\'" << varp->name() << "\')";
  else
    ostr << "(<UNDEFINED!>)";
}

bool VariableBase::validate_context(const char* name) const
{
  if (isconst() || within_evaluation() || (evaluation_state=CONTEXT_UNCHECKED))
    return true;
  const bool arrayvalid=(evaluation_state==CONTEXT_POSTCALC) || (evaluation_state==CONTEXT_FINALISE);
  if (arrayvalid) { //!< simple arrays are OK within processing blocks
    if (isexpr()) {
      if (sum_n0<2) //!< no summation arrays - always OK
	return true;
    }
    else {
      if (isarray() && (array.size()==1))
	return true;
    }
  }
//   if (isexpr())
//     ExpressionVariableFunction::outofmainloop_expression_warning.raise(name,true); //!< slightly weak - have to give warning in case expression depends on summation variable
//   else {
//     if (isarray())
//       ExpressionVariableFunction::outofmainloop_warning.raise(name,true);    
//   }
  const int use=uses();
  if (use & A_ARRAY)
    ExpressionVariableFunction::outofmainloop_warning.raise(name,true); 
  if (use & A_SUM)
    ExpressionVariableFunction::outofmainloop_expression_warning.raise(name,true);
  return false;
}

ContextWarning<> nonenvironment_warning("setenv variable name will not be recognised in pNMRsim as environment variable (contains lower case characters)",&NMRsim_repeat_warning);
 ContextWarning<> ignoredsetenv_warning("environment variable has already been parsed - this setenv will have no effect",&NMRsim_repeat_warning);

void parse_setenv()
{
#ifdef NOSETENV
  char totalcom[MAXLINE];
#endif
  char* cptr=get_curline();
  chewwhite(cptr);
  char* tail=cptr;
  while (*tail && !isspace(*tail))
    tail++;
  char* storeptr=tail;
  const char storechar=*storeptr;
  if (storechar)
    *tail++='\0';
  if (strchr(cptr,'='))
    error_abort("environment variable name can not contain =");
  chewwhite(tail);
  //if (*tail=='\0')
  //  error_abort("syntax: setenv <variable name> <definition>");
  
  if (!isenvironment(cptr))
    nonenvironment_warning.raise();

  if (isregistered(cptr))
    ignoredsetenv_warning.raise(cptr);

#ifndef NOSETENV
  const int retcode=setenv(cptr,tail,1);
#else
  snprintf(totalcom,sizeof(totalcom),"%s=%s",cptr,tail);  
  const int retcode=putenv(strdup(totalcom)); //!< added strdup 1/8/16. Leak, but not significant
#endif
  if (retcode) { //!< no need to create copy for setenv
    perror("setenv/putenv");
    error_abort();
  }
  
  *storeptr=storechar; //!< restore so original command line preserved
}

void scan_block(const char* name, bool isopt, bool allowmult)
 {
   for (;;) {
     parser_newblock(name,isopt);
     bool done=false;
     while (parser_getline(false)) {
       if (!parser_inblock())
	 std::cout << get_curline() << '\n';
       else {
	 if (!done) {
	   std::cout << name << " {\n";
	   done=true;
	 }
	 std::cout << OUTPUT_SPACER << get_curline() << '\n';
       }
     }     
     if (!done)
       return;
     std::cout << "}\n";
     if (!allowmult)
       return;
   }
 }  
  
