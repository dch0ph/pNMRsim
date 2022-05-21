#ifndef Parser_h_
#define Parser_h_

/*! \file
 \brief  Declarations for input file parsing
*/

#include "NMRsim.h"

#ifndef PATH_SEPARATOR
//! separator between path elements 
#define PATH_SEPARATOR ':'
#endif
#ifndef DIR_SEPARATOR
//! separator charactor within filename
#define DIR_SEPARATOR '/'
#endif

#define MAXLINE 4096 //!< max length of input lines
#define MAXPATH 4096 //!< max length of path name

void parser_init(const char*, int argc, char* argv[], int flags =0); //!< initialise parser
bool parser_inblock(); //!< \c true if parsing lines inside a block (as opposed to global)
const char* parser_getfname(); //!< get input file
void parser_newblock(const char* blockname_, bool, const char* keynames_[] =NMRSIM_NULL, size_t nkeys_ =0); //!< open new {} block
void parser_getinteractive(const char* intro, std::istream& =std::cin); //!< get parse line from input stream with prompt \a intro 
bool parser_foundlast(); //!< return \c true if last block was found
bool parser_isopen(); //!< return \c true if input file open
const char* clean_variable_name(char*); //!< strip any leading $, warn missing literal quotes
bool parser_getline(bool executeglobal =true); //!< get input line (NMRSIM_NULL if no more)
const char* parse_variable_name(); //!< get variable name
size_t parse_counting(const char*); //!< parse counting number (from 1) from string (0 if invalid)
void parser_insertline(const char*);
void parser_flush(); //!< finish with parser
void parser_checkfinished(); //!< generate error if arguments are still left
bool parser_isflag(); //!< return \c true if next token is a flag
bool parser_isnormal(); //!< return \c true if next token is a normal (non-flag) argument
bool parser_isalpha(); //!< return \c true if next token starts with alphabet character (i.e. possible variable name)
bool parser_isvariable_name(); //!< return \c true if next token could be a variable name - must be ' quoted
void parse_ignored(); //!< dummy for ignored instruction
Expression* parse_expression(const char*); //!< create new expression
int parse_index(size_t nspins_cell, size_t nspins_total =0); //!< parse spin number
double parse_double_raw(Variable*, char* cptr, int flags); //!< parse real quantity (raw)
bool read_block(const char* name,command_Factory_t& fac, bool isopt =false, bool ignoreunk =false); //!< parse {} block 
void check_variable_name(const char*); //!< check that variable/argument name is valid

extern ContextWarning<> emptyloop_warning;
extern ThreadWarning<> path_unreadable_warning;
extern ThreadWarning<> path_emptyelement_warning;
extern ContextWarning<> sharedtarget_warning; //!< different variables appear to point to same target
extern ContextWarning<> postcontinuation_whitespace_warning; //!< whitespace after \\ continuation
extern ContextWarning<> extradollar_warning; //!< variable *names* shouldn't have leading $
extern ContextWarning<> namequote_warning; //!< variable names are best quoted
extern ContextWarning<> filetail_ignored_warning; //!< rest of input file is being ignored
extern ContextWarning<> function_redefinition_warning; //!< function is (potentially) being redefined
extern ContextWarning<> misplaced_sum_warning; //!< sum array not meaningful in context
extern ContextWarning<> includeonce_withlist_warning; //!< list arguments used with includeonce
extern ContextWarning<> nofunctionarguments_warning; //!< function with no arguments
extern ContextWarning<> xy_variable_warning; //!< warn about using x or y as variable names
extern ContextWarning<> trailingcharacters_warning; //!< trailing characters after } block end
extern ContextWarning<> ignoreddirective_warning; //!< directive is being ignored
extern ContextWarning<> recursion_warning; //!< excessive recursion warning
extern ContextWarning<> zerolengthscalar_warning; //!< zero length argument when scalar expected
extern ContextWarning<> include_nosubstitution_warning; //!< dynamic include but nothing varies

typedef std::map<std::string, smartptr<UserVariable,false> > uservarmap_type;
extern uservarmap_type vardefmap; //!< registry for user variables

typedef FASTMAPTYPE(size_t) flagsmap_type; //!< type for set of flags - Note that flags must be non-volatile strings for FASTMAP

void declare_builtin_block(const char*); //!< register standard block

std::ostream& dumpflags(std::ostream&, const flagsmap_type&, size_t =0xFFFFFFFF); //!< stream flags set

extern LIST<productoperator_spec*> opspecstack; //!< stack of pointers to all operator specifications
//extern bool nondiagonal_opspec; //!< \c true if created non-diagonal product operator
void stripleaf(char* a, const char* trail); //!< strip \a trail from end of \a a if present
UserVariable* findnewvariable(const char* name); //!< create new (empty) ::UserVariable (error if already defined)
char* parse_string(int flags =0); //!< parse string
char* parse_string_syntax(const char*, size_t, int flags =0); //!< parse string
extern bool iscommonPAS; //!< \c true if all interactions share a common PAS
LIST<int> parse_intarray(int =0); //!< parse list of integers
LIST<size_t> parse_unsignedintarray(int offset =0, int  flags =0); //!< parse list of unsigned integers, with optional offset
LIST<size_t> parse_unsignedintarray_syntax(const char*, size_t, int offset =0, int  flags =0); //!< parse list of unsigned integers, with optional offset
int parse_int(Variable* =NMRSIM_NULL, int =0); //!< parse integer quantity
int parse_unsigned(Variable* =NMRSIM_NULL, int =0); //!< parse unsigned integer
int parse_unsigned_syntax(const char*, size_t, Variable* =NMRSIM_NULL, int =0); //!< parse unsigned integer
size_t parse_unsigned_offset(int offset, int =0); //!< parse unsigned integer after subtracting offset 
double parse_double(Variable* =NMRSIM_NULL, int =0); //!< parse real quantity
double parse_double_syntax(const char*, size_t, Variable* =NMRSIM_NULL, int =0); //!< parse real quantity
double parse_shift(Variable*, double gamrat, int flags =0); //!< parse shift quantity
size_t parse_nucleusname(const char*); //!< parse nucleus name (13C etc.)
VariableBase* parse_double_variable(Variable&, int =0);
VariableBase* parse_double_variable_syntax(Variable&, const char*, size_t, int =0);
template<class T> T parse_raw(char*, Variable*, int); //!< \internal raw parse
size_t parse_flag(const flagsmap_type&, const char*); //!< parse individual flag
size_t parse_flags(const flagsmap_type&, size_t def =0, bool oneonly =false); //!< parse flags
char* get_curline(); //!< get pointer to current parse line
void set_curline(char *); //!< set current parse point
char* get_token(int flags =0, const char* type =NMRSIM_NULL); //!< get next token (and advance)
size_t count_left(bool =true); //!< count remaining tokens
bool are_left(); //!< return \c true if unparsed tokens left
//PhasedSequence* parse_phasedsequence(int flags =0); //!< parse sequence+phase (NB be careful with setable phase issues)
class CycledSequence;
CycledSequence* parse_cycledsequence(int flags =0); //!< parse cycled sequence
CycledSequence* parse_cycledsequence(char*, int flags =0); //!< parse cycled sequence
double get_nmrfreq(size_t qualifier =0); //!< return current (or specified) NMR frequency
void scan_block(const char* name, bool isopt =false, bool allowmult =false);
void dumpenvironment(std::ostream& ostr); 
void process_array_ns(); //!< extract info from dimensionality
void printattributes(std::ostream&, int uses); //!< dump attributes

enum pm_flags { PM_FULL=1, //!< expand matrices
		PM_STRUCTURE=2, //!< output summary of matrix structure
		PM_ONCE=4, //!< only output once
		PM_EIGENBASIS=8, //!< output information on eigenbasis
		PM_STATISTICS=16, //!< output statistics
		PM_STRUCTUREFLAGS=2 | 8 | 16, //!< any structure flag
};

const flagsmap_type& putmatrix_flags();

struct Splitter {
  explicit Splitter(char* source, const char* charsv =",",  bool =false);
  
  char* next(char&);
  char* next() { char junk; return next(junk); } //discard
  bool atend() const { return ((curptr==NMRSIM_NULL) || !(*curptr)); }
  
  char* curptr;
  const char* chars;
  char* s;
  char sepchar;
  bool unique;
};

//! mark point for creation of new variables
class Mark {
public:
  Mark();
  ~Mark();
  void flush(Setable*); //!< flush (set pointer for unassigned variables)
  bool needsflush() const; //!< \c true if unassigned variables
  static ContextWarning<> unflushed_warning; //!< warning if object destroyed before flush
private: 
  size_t stackpos; //!< starting pointer into dirty link stack
};

//! strip leading white space
inline void chewwhite(char*& str) {
  while (isspace(*str)) { str++; } 
}

void verify_arraysize(size_t n, size_t dim =0); //!< register {} array size

template<class T> void verify_arraysizes(const ListList<T>& a, size_t sumtag, const BaseList<size_t>& arraytags, size_t itemlength)
{
  const size_t lists=a.size();
  if (lists>1)
    sum_dims.set(sumtag,lists);
  for (size_t i=lists;i--;) {
    const size_t dim = arraytags.empty() ? 0 : arraytags(i);
    verify_arraysize(a.size(i)/itemlength,dim);
  }
}

//const char* parser_getdirective(char* buf); //!< set parse line to \a buf and extract directive name

template<class Factory> void dump_map(const Factory& factory, std::ostream& ostr =std::cout)
{
  if (factory.empty()) {
    ostr << "<empty>\n";
    return;
  }
  const typename Factory::const_iterator end=factory.cend();
  typename Factory::const_iterator start=factory.cbegin();
  bool needspace=false;
  while (start!=end) {
    if (needspace)
      ostr << ' ';
    else
      needspace=true;
    ostr << (start->first);
    ++start;
  }
  ostr << '\n';
}

//! find key in parse command factory
template<class Factory> typename Factory::mapped_type factory_parse(const char* keyname, const Factory& factory)
{
  const typename Factory::const_iterator curcom=factory.find(keyname);
  if (curcom==factory.end()) {
    // if (verbose & VER_PARSE) {
    //   std::cerr << "Failed to find " << keyname << " in map:\n";
    //   dump_map(factory,std::cerr);
    //}
    if (!check_unrecognised(keyname)) {
      std::cerr << "Available options are: ";
      dump_map(factory,std::cerr);
    }
    error_abort();
  }       
  return curcom->second;
}

//!< call parsing function creating new object
template<class Func,class T> T* factory_call(Func funcp, Type2Type<T>)
{
  T* newac=NMRSIM_NULL;
  Mark markobj;
  try {
    newac=(*funcp)();
    if (newac)
      parser_checkfinished();
  } catch (MatrixException& exc) {
    parser_printcontext() << exc << '\n';
    error_abort();
  }
  if (markobj.needsflush()) {
    if (newac==NMRSIM_NULL)
      throw InternalError("factory_call null pointer");
    newac->isconstant(false); //flag is not const
    markobj.flush(newac);
  }
  return newac;
}

//!< find and call parsing function from key name
template<class Factory, class T> T* factory_parse(const char* keyname, const Factory& factory, Type2Type<T> type)
{
  const typename Factory::mapped_type funcp(factory_parse(keyname,factory));
  return factory_call(funcp,type);
}

template<class StackT,class ControlT> void verify_empty(StackT& pstack, const LIST<size_t>& controlstack, Type2Type<ControlT>)
{
  if (!(controlstack.empty())) {
    ControlT* headp=dynamic_cast<ControlT*>(pstack(controlstack.back()));    
    parser_printcontext() << "block finished before expected end " << (headp->type()) << '\n';  
    error_abort();
  }
}

template<class StackT> void verify_empty(StackT&, const LIST<size_t>&, Type2Type<NullType>) {}

template<class StackT, class ControlT> bool found_control(const char* keyname, StackT& pstack, LIST<size_t>& controlstack, Type2Type<ControlT>)
{
  if (strcmp(keyname,"end")!=0)
    return false;

  const char* rest=parse_string();
  if (controlstack.empty()) {
    parser_printcontext() << "end " << rest << " unmatched\n";	
    error_abort();
  }
  const size_t startindex=controlstack.back();
  ControlT* headp=dynamic_cast<ControlT*>(pstack(controlstack.back()));
  if (headp==NMRSIM_NULL)
    throw InternalError("factor_read: bad control block");
  if (strcmp(headp->type(),rest)!=0) {
    parser_printcontext() << "expecting end " << (headp->type()) << " but found end " << rest << '\n';
    error_abort();
  }
  const size_t items=pstack.size()-startindex-1; //!< subtract one for control block itself
  if (items==0)
    emptyloop_warning.raise();
  headp->items(items);
  controlstack.pop_back();
  if (verbose & VER_PARSE)
    std::cout << "Creating " << headp->type() << " control block with " << items << " items\n";
  return true;
}

template<class StackT> bool found_control(const char* keyname, StackT&, LIST<size_t>& controlstack, Type2Type<NullType>) { return false; }

//!< read complete block given parsing factory
template<class StackT, class Factory, class T,class ControlT> bool factory_read(StackT& pstack, const char* name, const Factory& factory, Type2Type<T> type, Type2Type<ControlT> ctype)
{
  if (!parser_isopen())
    return false;

  LIST<size_t> controlstack;

  parser_newblock(name,true);
  bool found=false;
  while (parser_getline()) {
    const char* keyname=parse_string();
    if (!found_control(keyname,pstack,controlstack,ctype)) {
      T* newproc=factory_parse(keyname,factory,type);
      if (newproc) {
	ControlT* headp=dynamic_cast<ControlT*>(newproc);
	if (headp)
	  controlstack.push_back(pstack.size());
	pstack.push_back(newproc);
	found=true;
      }
    }
  }
  verify_empty(pstack,controlstack,ctype);
  return found;
}

//!< read blocks (which may be repeated) using given parsing factory
template<class StackT, class Factory, class T, class ControlT> size_t multiple_factory_read(LIST<StackT>& mstack, const char* name, const Factory& factory, Type2Type<T> type, Type2Type<ControlT> ctype)
{
  if (!parser_isopen())
    return 0U;
  size_t count=0U;
  for (;;count++) {
    mstack.push_back(); //!< call constructor for set up
    if (!factory_read(mstack.back(),name,factory,type,ctype)) {
      mstack.pop_back();
      return count;
    }
    mstack.back().initialise(); //!< finished parsing, do any tidy up
  }
}

template<class StackT> int overall_attributes(const StackT& stack)
{
  const typename StackT::const_iterator end(stack.end());
  typename StackT::const_iterator start(stack.begin());
  int attr=0;
  while (start!=end) {
    attr |= (*start)->setable_uses();
    ++start;
  }
  return attr;
}

template<class StackT> int multiple_overall_attributes(const BaseList<StackT>& mstack)
{
  int attr=0;
  for (size_t i=mstack.size();i--;)
    attr|=overall_attributes(mstack(i));
  return attr;
}

extern ContextWarning<> inconsistent_const_warning;
extern ContextWarning<> inconsistent_nouses_warning;

extern bool have_virtualdimensions; //!< set to true if virtual dimensions created

#endif
