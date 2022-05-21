/*! \file
 \brief  Parsing for basic objects
*/

#include "parser_common.hpp"
#include "expression_definition.hpp"
#include "ListList.h"
#include "assign_const_action.hpp"

//! limit for truncation of list (avoid by verbose -general and -debug)
#define DEFAULT_TRUNCATIONLIMIT 30

LIST<size_t> suminds;

//! functor to set error step
struct VariableBuilder::set_step
{
  set_step(VariableBuilder& builderv) : builder_(builderv) {}
  //! assign error step
  void operator()(double val) const {
    builder_.var_step=std::fabs(val);
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

//! functor resetting variable builder
struct VariableBuilder::ResetSimple
{
  explicit ResetSimple(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const {
    builder_.resetsimple();
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

//! functor adding array entry
struct VariableBuilder::Add
{
  explicit Add(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const {
    builder_.add();
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

//! functor adding array list to sum array
struct VariableBuilder::AddList
{
  explicit AddList(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const {
    builder_.addlist();
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

//! functor adding simple value
struct VariableBuilder::AddSimple
{
  explicit AddSimple(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const {
    builder_.addsimple();
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

//! functor adding list to sum list
struct VariableBuilder::AddSumList
{
  explicit AddSumList(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const {
    builder_.addsumlist();
  }
  VariableBuilder& builder_; //!< reference to variable builder
};

template<typename IteratorT> void addfiltered(std::string& dest, IteratorT st, IteratorT end, char filter)
{
//     st++; //!< discard ''
//     end--;
  if (end<st)
    throw InternalError("addfiltered");
  dest.reserve(dest.size()+(end-st));
  while (st<end) {
    char c=*st++;
    if ((c=='\\') && (st!=end)) {
      c=*st++;
      if (c!=filter)
	dest.push_back('\\');
    }
    dest.push_back(c);
  }
}

//! functor adding ::ExpressionInclude
/** \c IsDynamic distinguishes between dynamic `` and once-only "" */
template<bool IsDynamic> struct ExpressionInclude::Add
{
  explicit Add(VariableBuilder& exprv)
    : exprbuilder_(exprv) {}
  //! Add new ::ExpressionInclude
  template<typename IteratorT> void operator()(const IteratorT& start, const IteratorT& end) const {
    std::string tmp;
    addfiltered(tmp,start,end, IsDynamic ? '`' : '"');
    exprbuilder_.push(new ExpressionInclude(tmp.c_str(),IsDynamic));
  }
  VariableBuilder& exprbuilder_; //!< reference to expression builder
};

//! functor adding ::ExpressionRange
/** \c HasMiddle is \c true if step is specified */
template<bool HasMiddle> struct ExpressionRange::Add
{
  explicit Add(VariableBuilder& expr) : exprbuilder_(expr) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    const size_t last=exprbuilder_.pop();
    const size_t nextlast=exprbuilder_.pop();
    exprbuilder_.push(HasMiddle 
		      ? new ExpressionRange(exprbuilder_.pop(),nextlast,last)
		      : new ExpressionRange(nextlast,last));
  } 
  VariableBuilder& exprbuilder_;
};

//! functor adding ::ExpressionConstant to expression
struct ExpressionConstant::Add
{
  explicit Add(VariableBuilder& exprv)
    : exprbuilder_(exprv) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    exprbuilder_.push(new ExpressionConstant(exprbuilder_.var_value,exprbuilder_.isint),false); //don't check for const-ness!
  }
  void operator()(double) const {
    exprbuilder_.push(new ExpressionConstant(exprbuilder_.var_value,exprbuilder_.isint),false); //don't check for const-ness!
  }
  VariableBuilder& exprbuilder_;
};

struct VariableBuilder::AddSlot
{
  explicit AddSlot(VariableBuilder& exprv)
    : exprbuilder_(exprv) {}

  template<typename IteratorT> void operator()(const IteratorT& start, const IteratorT& end) const {
    if (!(exprbuilder_.flags & F_ALLOWSLOT))
      error_abort("#<n> argument specifiers only valid in function definitions");
    const std::string name(start,end);
    const size_t slot=ExpressionSlot::findslot(name.c_str());
    exprbuilder_.push(new ExpressionSlot(slot),false); //don't check for const-ness!
  }

  void operator()(unsigned int slot) const {    
    if (!(exprbuilder_.flags & F_ALLOWSLOT))
      error_abort("#<n> argument specifiers only valid in function definitions");
    //    const bool isstar=exprbuilder_.isstar;
    exprbuilder_.push(new ExpressionSlot(slot),false); //don't check for const-ness!
    //    if (slot>exprbuilder_.maxslot)
    //  exprbuilder_.maxslot=slot;
    //exprbuilder_.maxisstar|=isstar;
  }

  VariableBuilder& exprbuilder_;
};
  
//! functor adding ::ExpressionShift
struct ExpressionShift::Add
{
  explicit Add(VariableBuilder& exprv) : exprbuilder_(exprv) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const { doadd(); }
  //  template<typename charT> void operator()(charT) const { doadd(); }
  void doadd() const;
  VariableBuilder& exprbuilder_;
};

void ExpressionShift::Add::doadd() const
{
  //  if (!(exprbuilder_.flags & F_ISSHIFT))
  // error_abort("Trying to create chemical shift in unexpected context - <number>p has confused parser?");
  ExpressionShift* newobjp=new ExpressionShift(exprbuilder_.var_value, exprbuilder_.ppmqual);
  exprbuilder_.push(newobjp);
}

//! functor adding phase quantity
/** \a Quad parameter distinguishes 90 shifts: 0, 1, 2, 3 */  
template<int Quad> struct ExpressionConstant::AddPhase
{
  explicit AddPhase(VariableBuilder& exprv) : exprbuilder_(exprv) {}

  template<typename charT> void operator()(charT) const {
//     if (!(exprbuilder_.flags & F_ISPHASE))
//       throw Failed("Attempt to use phase definitions in non-phase parameter?");
    exprbuilder_.push(new ExpressionConstant(90.0*Quad),false);
  }
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
//     if (!(exprbuilder_.flags & F_ISPHASE))
//       throw Failed("Attempt to use phase definitions in non-phase parameter?");
    exprbuilder_.push(new ExpressionConstant(90.0*Quad),false);
  }
  VariableBuilder& exprbuilder_;
};

//! functor adding unary expression node
template<typename T> struct ExpressionUnary::Add
{
  explicit Add(VariableBuilder& expr) : exprbuilder_(expr) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    exprbuilder_.push(new T(exprbuilder_.pop()));
  }
  VariableBuilder& exprbuilder_;
};

//! functor adding binary expression node
template<typename T> struct ExpressionBinary::Add
{
  explicit Add(VariableBuilder& expr) : exprbuilder_(expr) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    const size_t last=exprbuilder_.pop();
    exprbuilder_.push(new T(exprbuilder_.pop(),last));
  }
  VariableBuilder& exprbuilder_;
};

//! function adding string argument
struct ExpressionStringFunction::AddString {
  explicit AddString(VariableBuilder& expr) : exprbuilder_(expr) {}
  template<typename IteratorT> void operator()(IteratorT st, IteratorT end) const {
    exprbuilder_.stringstack_.push_back("");
    addfiltered(exprbuilder_.stringstack_.back(),st,end,'\'');
    //    exprbuilder_.stringstack_.back().append(st,end-st);
  }
  VariableBuilder& exprbuilder_;
};

//! functor starting list/function
struct ExpressionFunctionBase::Start {
  explicit Start(VariableBuilder& expr) : exprbuilder_(expr) {}
  template<typename IteratorT> void operator()(const IteratorT& st, IteratorT end) const {
    exprbuilder_.mark_push(st,--end); //take off last char (parenthesis)
    exprbuilder_.stringstack_.create(0U);
  }
  template<typename T> void operator()(T) const {
    exprbuilder_.mark_push();
  }
  VariableBuilder& exprbuilder_;
};

//! functor ending list/function
struct ExpressionFunctionBase::End {
  End(VariableBuilder& expr, bool isfunc =true) : exprbuilder_(expr), isfunc_(isfunc) {}
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    exprbuilder_.endfunction(isfunc_);
  }
  template<typename T> void operator()(T) const { exprbuilder_.endfunction(isfunc_); }
  VariableBuilder& exprbuilder_;
  bool isfunc_; //!< \c true if building function (rather than list)
};

//extern bool parser_init;

struct VariableBuilder::Mismatched
{
  explicit Mismatched(char parenv, const char* typev =NMRSIM_NULL)
    : paren_(parenv), type_(typev) {}
  template<typename IteratorT> void operator()(IteratorT st, const IteratorT& en) const {
    parser_printcontext() << "Parsing stopped before finding closing " << paren_;
    if (type_)
      std::cerr << " in " << type_;    
    std::cerr << ": ";
    while (st!=en) {
      std::cerr << *st;
      ++st;
    }
    std::cerr << '\n';
    error_abort();
  }
  const char paren_;
  const char* type_;
};

struct ExpressionVariable::ParseFailed
{
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    error_abort("Failed to parse variable name.  Syntax is $<name> or $(<name>) where <name> contains alphanumeric characters (but does not start with a number");
  }
};

struct VariableBuilder::FailedTag
{
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    error_abort("Failed to parse virtual dimension qualifier.  Syntax is :<dimension> where <dimension> is a dimension number 1,2,...");
  }
};

struct VariableBuilder::FailedSlot
{
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    error_abort("Failed to parse function argument specification.  Syntax is #<argument number> where <argument number> is an integer 1,2,...");
  }
};

//! function adding $ variable node
struct ExpressionVariable::Add 
{
  explicit Add(VariableBuilder& exprv) : exprbuilder_(exprv) {}
  template<typename IteratorT> void operator()(const IteratorT& start, const IteratorT& end) const {
    std::string name(start,end);
    RealVariable& var(parse_variable_name(name.c_str()));
    exprbuilder_.push(new ExpressionVariable(var));
  }
  VariableBuilder& exprbuilder_;
};

template <typename ScannerT> expr_grammar::definition<ScannerT>::definition(const expr_grammar& self)
{
   varname =  lexeme_d[ (alpha_p | '_') >> *(alnum_p | '_') ];

   stepqualdef = lexeme_d[ureal_p[VariableBuilder::set_step(self.builder_)]
      >> !( baseppmrule[assign_const_a(self.builder_.stepqual,'p')] | ch_p('%')[assign_a(self.builder_.stepqual)])];

    Vdef = lexeme_d[ch_p('V')[assign_const_a(self.builder_.var_step,0.0)] || stepqualdef];
    
    sumtag = ch_p(':') >> ( uint_p[assign_a(self.builder_.tagbuilder_.sumtag_)] | eps_p[VariableBuilder::FailedTag()]);

    slot = lexeme_d[ch_p('#') >> (uint_p[VariableBuilder::AddSlot(self.builder_)] | varname[VariableBuilder::AddSlot(self.builder_)] | eps_p[VariableBuilder::FailedSlot()]) ];

    arraytag = (ch_p(':') >> (uint_p[assign_a(self.builder_.tagbuilder_.arraytag_)] | eps_p[VariableBuilder::FailedTag()]))
      | eps_p[assign_const_a(self.builder_.tagbuilder_.arraytag_,0)];
		   
   phasedef
     = ch_p('x')[ExpressionConstant::AddPhase<0>(self.builder_)]
     | ch_p('y')[ExpressionConstant::AddPhase<1>(self.builder_)]
     | str_p("-x")[ExpressionConstant::AddPhase<2>(self.builder_)]
     | str_p("-y")[ExpressionConstant::AddPhase<3>(self.builder_)];

   ppmrule = lexeme_d[ch_p('p')[assign_const_a(self.builder_.ppmqual,0)]
		      >> !(uint_p[assign_a(self.builder_.ppmqual)])
		      >> ~eps_p(alnum_p-ch_p('V'))];

   baseppmrule = lexeme_d[ch_p('p') >> ~eps_p(alnum_p-ch_p('V'))];

   variable
     = lexeme_d[ch_p('$')
		>> ( ( (ch_p('(') >> varname[ExpressionVariable::Add(self.builder_)] >> (ch_p(')') | eps_p[VariableBuilder::Mismatched(')',"$(<variable>) specification")]) ) 
 		      | varname[ExpressionVariable::Add(self.builder_)] ) | eps_p[ExpressionVariable::ParseFailed()])];
  
   //   stringliteral = lexeme_d[ch_p('\'') >> (*(anychar_p - ch_p('\'')))[ExpressionStringFunction::AddString(self.builder_)] >> (ch_p('\'') | eps_p[VariableBuilder::Mismatched('\'')])];
   stringliteral = lexeme_d[confix_p(ch_p('\''), (*anychar_p)[ExpressionStringFunction::AddString(self.builder_)], eps_p(~ch_p('\\')) >> ch_p('\''))];

   function
     = lexeme_d[varname >> ch_p('(')][ExpressionFunctionBase::Start(self.builder_)]
     >> ( list_p(stringliteral, ch_p(','),')') | !list_p( expression, ch_p(','), ')'))
     >> (ch_p(')')[ExpressionFunctionBase::End(self.builder_)] | eps_p[VariableBuilder::Mismatched(')',"function specification")]);
  
   staticinclude = ch_p('"') >> lexeme_d[+((graph_p | blank_p)-ch_p('"'))][ExpressionInclude::Add<false>(self.builder_)] >> (ch_p('"') | eps_p[VariableBuilder::Mismatched('"')]);
   //   staticinclude = confix_p( ch_p('"') , lexeme_d[*(graph_p | blank_p)][ExpressionInclude::Add<false>(self.builder_)] , ~eps_p('\\') >> ch_p('"') );

   //dynamicinclude = confix_p( ch_p('`') , lexeme_d[*(graph_p | blank_p)][ExpressionInclude::Add<true>(self.builder_)] , ~eps_p('\\') >> ch_p('`') );
   //dynamicinclude = confix_p( ch_p('`') , lexeme_d[*( (graph_p | blank_p)-ch_p('`'))][ExpressionInclude::Add<true>(self.builder_)] , ch_p('`') );
  
   dynamicinclude = ch_p('`') >> lexeme_d[+((graph_p | blank_p)-ch_p('`'))][ExpressionInclude::Add<true>(self.builder_)] >> (ch_p('`') | eps_p[VariableBuilder::Mismatched('`')]);
  
   range
     = expression
     || (ch_p(':') >> expression
 	>> ( (ch_p(':') >> expression)[ExpressionRange::Add<true>(self.builder_)] | eps_p[ExpressionRange::Add<false>(self.builder_)]));

   varrange
     = (simplevariable[VariableBuilder::AddSimple(self.builder_)]
        | (expression
 	  || (ch_p(':') >> expression
 	      >> ( (ch_p(':') >> expression)[ExpressionRange::Add<true>(self.builder_)] | eps_p[ExpressionRange::Add<false>(self.builder_)])))[VariableBuilder::Add(self.builder_)]
        );

   //  ppmrule = lexeme_d[(str_p("ppm") | ch_p('p')) >> ~eps_p(alnum_p-ch_p('V'))];

   simplevariable = lexeme_d[eps_p[VariableBuilder::ResetSimple(self.builder_)]
   			    >> real_p[assign_a(self.builder_.var_value)]
 			    >> !(ppmrule[assign_const_a(self.builder_.valisp,true)])
			     >> Vdef
			     //			     >> (ch_p('V')[assign_const_a(self.builder_.var_step,0.0)] || stepqualdef)
 			    ];

   listel
     = staticinclude | dynamicinclude | range;

   varlistel
     = (staticinclude | dynamicinclude)[VariableBuilder::Add(self.builder_)]
     | varrange;
  
   rawlist = eps_p[ExpressionFunctionBase::Start(self.builder_)]
     >> list_p( listel, ch_p(','))[ExpressionFunctionBase::End(self.builder_,false)] >> eps_p[VariableBuilder::Cleanup<VariableBuilder::EXPR>(self.builder_)];
  
   plainlist = ch_p('[')[ExpressionFunctionBase::Start(self.builder_)]
     >> (!list_p( listel, ch_p(',') , ']') >> (ch_p(']') | eps_p[VariableBuilder::Mismatched(']')]) )[ExpressionFunctionBase::End(self.builder_,false)];
  
   arraylist = ch_p('{')
     >> list_p( varlistel , ch_p(','), '}')
     >> (ch_p('}') | eps_p[VariableBuilder::Mismatched('}')])
     >> arraytag[VariableBuilder::AddList(self.builder_)];
  
   sumlistel = arraylist | varlistel[VariableBuilder::AddSumList(self.builder_)];

   sumlist = ch_p('|')
     >> list_p( sumlistel, ch_p(','),'|')
     >> ( (ch_p('|') >> !lexeme_d[Vdef]) | eps_p[VariableBuilder::Mismatched('|',"|| summation list")])
     >> !sumtag;
  
   constarray = sumlist | arraylist;
  
   basefactor_rule
     = lexeme_d[
 	       (strict_real_p[assign_a(self.builder_.var_value)][assign_const_a(self.builder_.isint,false)]
 		| int_p[assign_a(self.builder_.var_value)][assign_const_a(self.builder_.isint,true)])
 	       >> (ppmrule[ExpressionShift::Add(self.builder_)] | eps_p[ExpressionConstant::Add(self.builder_)])
 		    ]
     | lexeme_d[variable]
     | (ch_p('(') >> expression >> (ch_p(')') | eps_p[VariableBuilder::Mismatched(')',"() grouping [or variable quantity placed within ()]")]))
     | plainlist
     | function
     | dynamicinclude
     | lexeme_d[slot]
     | (ch_p('-') >> basefactor_rule)[ExpressionUnary::Add<ExpressionMinus>(self.builder_)]
     ;
    
   factor = basefactor_rule || (ch_p('^') >> basefactor_rule)[ExpressionBinary::Add<ExpressionPower>(self.builder_)];

  //phasedef must be first otherwise x etc. parsed as start of function
    term =
      (phasedef | factor) >>
     *(  (ch_p('*') >> factor)[ExpressionBinary::Add<ExpressionMultiply>(self.builder_)]
 	| (ch_p('/') >> factor)[ExpressionBinary::Add<ExpressionDivide>(self.builder_)]
 	| (ch_p('%') >> factor)[ExpressionBinary::Add<ExpressionModulus>(self.builder_)]
 	);
  
   baseterm =  term >>
     *(  (ch_p('+') >> term)[ExpressionBinary::Add<ExpressionAdd>(self.builder_)]
 	| (ch_p('-') >> term)[ExpressionBinary::Add<ExpressionSubtract>(self.builder_)]
 	);

    expression = baseterm >>
      !( (ch_p('>') >> baseterm)[ExpressionBinary::Add<ExpressionGT>(self.builder_)]
        | (ch_p('<') >> baseterm)[ExpressionBinary::Add<ExpressionLT>(self.builder_)] 
        );

   supervariable_rule
     = (simplevariable >> end_p)[VariableBuilder::Cleanup<VariableBuilder::SIMPLE>(self.builder_)]
     | (constarray >> end_p)[VariableBuilder::Cleanup<VariableBuilder::ARRAY>(self.builder_)]
     | expression[VariableBuilder::Cleanup<VariableBuilder::EXPR>(self.builder_)];
  
  self.isinit=true; //flag rules created
  this->start_parsers(supervariable_rule,rawlist,expression);
}   

parser_assert_t VariableBuilder::stack_empty("VariableBuilder: stack exhausted");

bool VariableBuilder::checkclean()
{
  if (curstack_.size()>1) {
    if (verbose & VER_PARSE)
      parser_printthread(std::cerr) << "VariableBuilder: non-trivial stack: " << curstack_ << '\n';
    return false;
  }
  if (!markstack_.empty()) {
    if (verbose & VER_PARSE)
      parser_printthread(std::cerr) << "VariableBuilder: non-empty mark stack\n";
    return false;
  }
  return true;
}

VariableBuilder::~VariableBuilder()
{
  checkclean();
}

void VariableBuilder::reset()
{
  array_.clear();
  steps_.create(0);
  exprstmp_.create(0);
  offsetstmp_.create(0);
  haveexprs_=false;
  tagbuilder_.reset();
  resetexpr();
  cvar_.initialise(false);
  built=false;
  //  maxslot=0;
  //maxisstar=false;
}

void VariableBuilder::reset(int flagsv)
{
  flags=flagsv;
  reset();
}

void TagBuilder::reset()
{
  hastagspec_=false;
  arraytags_.create(0);
  sumtag_=0;
}

void VariableBuilder::resetsimple()
{
  valisp=false;
  stepqual='\0';
  var_step=-1.0;
}

void VariableBuilder::resetexpr()
{
  expr_.create(0);
  curstack_.create(0);
  markstack_.clear();
  clearevalstacks();  //!< ensure argument stack is clear (this signals expression is being constructed rather than evaluated
}

size_t VariableBuilder::pop()
{
  if (curstack_.empty())
    throw stack_empty;
  const size_t val=curstack_.back();
  curstack_.pop_back();
  return val;
}

const MarkState_t<>& VariableBuilder::mark_back()
{
  if (markstack_.empty())
    throw stack_empty;
  MarkState_t<>& mark=markstack_.back();
  mark.nargs=stringstack_.empty() ? (curstack_.size()-mark.stackpos) : stringstack_.size();
  return mark;
}

bool VariableBuilder::push(ExpressionBase* ptr, bool constcheck) 
{
  const bool isverbose=(parser_verbose_level()>1);
  bool isconst=true;
  if (constcheck) {
    int root=-1;
    isconst=ptr->isconstant(root,expr_);
    if (isconst) { //constant sub-expression
      try {
	ExpressionConstant* newptr=new ExpressionConstant(*ptr,expr_);
	if (isverbose) {
	  std::cout << "Eliminating constant sub-expression back to ";
	  ptr->print(std::cout,expr_);
	  std::cout << ": " << (*newptr)() << '\n';
	}
	if (root>=0)
	  expr_.create(root);
	delete ptr;
	ptr=newptr;
      }
      catch (const MatrixException& exc) {
	if (isverbose)
	  std::cout << "Evaluation of constant sub-expression failed: " << exc;
      }
      catch (const std::exception& exc) {
	if (isverbose)
	  std::cout << "Evaluation of constant sub-expression failed: " << exc.what() << '\n';
      }
    }
  }
  const size_t n=expr_.push_original(ptr);
  if (isverbose)
    std::cout << "Pushed to position: " << n << '\n';    
  if (!curstack_.empty() && (curstack_.back()>=n))
    throw InternalError("Stack crash!");
  curstack_.push_back(n);
  return isconst;
}

//ContextWarning<> VariableBuilder::nonconst_in_array_warning("non-const expression encountered when constructing constant array - behaviour may not be as expected",&NMRsim_repeat_warning);
ContextWarning<> VariableBuilder::zerosizeexpr_in_array_warning("variable expression returning empty result in array - may (incorrectly) depend on variables that are only defined inside evaluation loop?",&NMRsim_repeat_warning);

void VariableBuilder::evaluate(LIST<double>& res)
{
  //  if (!expr_.isconstant())
  //  nonconst_in_array_warning.raise();
  //  throw Failed("Expression must evaluate to constant");    
  expr_.get(res);
  size_t nels=res.size();
  Expression* exprp=NMRSIM_NULL;
  if (expr_.isconstant())
    offsetstmp_.push_back(nels,size_t(-1));
  else {
    if (nels==0)
      zerosizeexpr_in_array_warning.raise();
    else {
      exprp=new Expression();
      exprp->swap(expr_);
      haveexprs_=true;
      for (size_t i=0;i<nels;i++)
	offsetstmp_.push_back(i);
    }
  }
  exprstmp_.push_back(nels,exprp);
  resetexpr(); //clear expression
}

void VariableBuilder::add()
{
  size_t n;
  if (scr_.empty()) {
    evaluate(scr_);
    n=scr_.size();
  }
  else {
    evaluate(tmp_);
    n=tmp_.size();
    scr_.push_back(tmp_);
  }
  steps_.push_back(n,0.0); //not varied
}

void TagBuilder::flushlist()
{
  arraytags_.push_back(arraytag_);
  if (arraytag_)
    hastagspec_=true; //flag virtual dimension present
}

void VariableBuilder::flushlist()
{
  scr_.create(0);
}

void VariableBuilder::addlist()
{
  if (scr_.empty())
    error_abort("List is empty");
  array_.push_back(scr_);
  tagbuilder_.flushlist();
  flushlist();
}

void VariableBuilder::addsumlist()
{
  if (scr_.empty())
    error_abort("Component of sum list is empty");
  for (size_t i=0;i<scr_.size();i++)
    array_.push_back(BaseList<double>(1,&(scr_(i))));
  tagbuilder_.arraytags_.push_back(scr_.size(),0); //individual elements can't have tags
  flushlist();
}

template<> void VariableBuilder::cleanup<VariableBuilder::EXPR>()
{  
  steps_.create(0); //remove any values lying around
  if (expr_.isconstant()) {
    cvar_.initialise(true);
    expr_.get(cvar_.value_);
    if ((cvar_.value_.size()!=1) && !(flags & F_ALLOWLIST))
      error_abort("Can't use list result here");
    return;
  }
  if (flags & F_DENYEXPR)
    error_abort("Can't use non-const expressions here");
  cvar_.initialise(false);
  cvar_.expr_.swap(expr_);
  //  if (!(flags & F_ALLOWSLOT))
  //  cvar_.update(); //!< don't try to evaluate if function
}

double VariableBuilder::checkV(double val)
{
  if (var_step<0) //simple constant
    return 0.0;

  if (flags & F_DENYVAR)
    error_abort("variable quantity not permitted");
  if (flags & F_ISINTEGER)
    error_abort("integer quantity cannot be made variable (V)");
  if ((flags & F_DENYZERO) && (val==0.0))
    error_abort("can't allow this parameter to vary through zero");

  switch (stepqual) {
  case '\0': 
    if (var_step==0.0)
      var_step=std::fabs(val*NMRSIM_DEFAULT_ERROR); //!< apply default error
    break;

  case 'p':
    if (!proton_freq_isconstant())
      error_abort("Can't combine variable proton frequency with variable shifts");
    ::fudge_shift(var_step);
    var_step=std::fabs(var_step);
    break;

  case '%':
    if (var_step==0.0)
      error_abort("can't use 0 as a scaling factor!");
    var_step*=std::fabs(val*0.01);
    break;

  default:
    throw InternalError("Unknown step qualifier");
  }
  if (var_step==0.0)
    error_abort("error cannot be determined from zero parameter - set explicitly");

  return var_step;
}

double VariableBuilder::rawcleanup(bool denyexpr)
{
  if (flags & F_DENYEXPR)
    denyexpr=true;

  if (valisp || (stepqual=='p')) {
//     if (!(flags & F_ISSHIFT))
//       error_abort("'p' specified for non-shift quantity");
    if (!proton_freq_isconstant()) {
      if (var_step>=0.0)
	error_abort("Can't combine variable proton frequency with variable shifts");
      if (denyexpr)	
	error_abort("expression not permitted");
      if (valisp) {
	expr_.push(new ExpressionShift(var_value,ppmqual));
	cleanup<EXPR>();
	return 0.0;
      }	
    }
  }
  if (valisp)
    ::fudge_shift(var_value,ppmqual);
  const double vstep=checkV(var_value);
  if (vstep==0.0) {//simple constant
    cvar_.initialise(var_value,0.0);
    var_step=0.0; //!< reset var_step
    return 0.0;
  }
  return vstep;
}

void VariableBuilder::addsimple()
{
  const double vstep=rawcleanup(true);
  var_step=-1.0;
  scr_.push_back(var_value);
  steps_.push_back(vstep);
  offsetstmp_.push_back(-1);
  exprstmp_.push_back((Expression*)NMRSIM_NULL);
}  

template<> void VariableBuilder::cleanup<VariableBuilder::SIMPLE>()
{
  const double vstep=rawcleanup(false);
  var_step=-1.0;
  if (vstep!=0.0) //!< already initialised/setup if vstep=0
    cvar_.initialise(var_value,vstep);
}

ContextWarning<> ignoredV_warning("Ignored wider scope V specification in favour of local",&NMRsim_repeat_warning);

template<> void VariableBuilder::cleanup<VariableBuilder::ARRAY>()
{
  if ((array_.items()==1) && !haveexprs_) {
    double err=0.0;
    const double val=array_(size_t(0),size_t(0));
    switch (steps_.size()) {
    case 0:
      err=checkV(val);
      break;
    case 1:
      if (var_step>=0.0)
	ignoredV_warning.raise();	
      err=steps_.front(); 
      break;
    default:
      throw InternalError("cleanup<ARRAY>");
    }
    cvar_.initialise(val,err);
    return;
  }
  if (flags & F_DENYARRAY) {
    for (size_t i=array_.size();i--;) {
      if (array_(i).size()>1)
	error_abort("{} arrayed quantity not permitted in this context");
    }
  }
  if ((flags & F_DENYSUM) && (array_.size()>1))
    error_abort("|| arrayed quantity not permitted in this context");
  cvar_.initialise(false);
  cvar_.array.swap(array_);
  cvar_.arraylength_=1;
  cvar_.errors_.swap(steps_);
  if (var_step>=0.0) {
    const size_t n=cvar_.array.items();
    const BaseList<double> asrow(cvar_.array.row());
    if (cvar_.errors_.empty()) {
      cvar_.errors_.create(n);
      for (size_t i=n;i--;)
	cvar_.errors_(i)=checkV(asrow(i));
    }
    else {
      bool foundprob=false;
      for (size_t i=n;i--;) {
	if (cvar_.errors_(i))
	  foundprob=true;
	else
	  cvar_.errors_(i)=checkV(asrow(i));
      }
      if (foundprob)
	ignoredV_warning.raise();
    }
  }
  cvar_.sumtag_=tagbuilder_.sumtag_;
  cvar_.arraytags_.swap(tagbuilder_.arraytags_);
  if (haveexprs_) {
    cvar_.arrayexprs_.swap(exprstmp_);
    cvar_.exproffsets_.swap(offsetstmp_);
  }
  cvar_.update();
  var_step=-1.0; //!< clear any V setting
  steps_.create(0);
}

void ExpressionStringFunction::ExpressionStringFunction::create(BaseList<std::string>& args)
{
  args_.create(args.size());
  for (size_t i=args.size();i--;)
    args_(i).swap(args(i));
}

bool VariableBuilder::endsequence(ExpressionBase* objp)//, bool constcheck)
{
  const MarkState_t<>& fstate=mark_back();
  const size_t nargs=fstate.nargs; 
  if (parser_verbose_level()>1)
    std::cout << "Creating new node (function/list/array) with " << nargs << " children\n";
  if (stringstack_.empty()) {
    objp->children.create(nargs);
    for (size_t i=nargs;i--;)
      objp->children(i)=pop();
  }
  else {
    ExpressionStringFunction* sobjp=dynamic_cast<ExpressionStringFunction*>(objp);
    if (!sobjp)
      error_abort("string arguments passed to function not expecting strings");
    sobjp->create(stringstack_);
    stringstack_.clear();
  }
  return push(objp);
}
 
const ExpressionNamedBase& function_def_t::operator*() const
{
  if (!objp)
    throw InternalError("Missing function pointer");
  return *objp;
}

void VariableBuilder::endfunction(bool isfunc)
{
  const bool isstring=!stringstack_.empty();
  if (!isfunc) {
    if (isstring)
      throw InternalError("list of strings?");
    (void)endsequence(new ExpressionSimpleFunction());
  }
  else {
    const MarkState_t<>& fstate=mark_back();
    const size_t nargs=fstate.nargs;
    
    Function_Factory_t& Function_Factory(get_Function_Factory());
    const Function_Factory_t::const_iterator end(Function_Factory.end());
    const std::string name(fstate.begin,fstate.end);
    Function_Factory_t::iterator iter=Function_Factory.find(function_spec(name.c_str(),-1)); //!< check for multi-variate function
    if (iter==end) {
      iter=Function_Factory.find(function_spec(name.c_str(),nargs));
      if (iter==end) {
	parser_printcontext() << "No function " << name << " with " << nargs << " argument(s)\n";
	Function_Factory_t::const_iterator start(Function_Factory.begin());
	while (start!=end) {
	  std::cout << (start->first).first << " with ";
	  const int n=(start->first).second;
	  if (n<0)
	    std::cout << "variable number of arguments\n";
	  else
	    std::cout << n << " argument(s)\n";
	  start++;
	}
	error_abort();
      }
    }
    function_def_t& funcdef(iter->second);
    ExpressionNamedBase* funcp=dynamic_cast< ExpressionNamedBase* >((*funcdef).clone());
    if (!funcp)
      throw InternalError("function not derived from ExpressionNamedBase");
    if (!endsequence(funcp)) {//!< non-const use of function
      funcdef.used=true;
      if ((verbose & VER_PARSE) && (verbose_level>1))
	std::cout << "Flagging " << funcp->name() << " as used\n";
    }
  }
  mark_pop();
}

void VariableBuilder::fudge_shift(size_t qual)
{
//   if (!(flags & F_ISSHIFT))
//     throw Failed("Can't specify 'p' for non-shift parameter");
  ExpressionConstant* objp=dynamic_cast<ExpressionConstant*>(back());
  if (objp) {
    LIST<double>& vals(objp->value_);
    if (vals.size()==1) {
      ::fudge_shift(vals.front(),qual);
      return;
    }
  }
  throw InternalError("FudgeShift");
}

void VariableBase::set_error(double err, size_t i)
{
  if (isexpr())
    throw InternalError("VariableBase::set_error");
  if (i>=errors_.size())
    throw BadIndex("VariableBase::set_error",i,errors_.size());
  errors_(i)=err;
}

double VariableBase::get_error(size_t i) const
{
  if (isexpr())
    throw InternalError("VariableBase::get_error");
  if (errors_.empty())
    return 0.0;
  if (i>=errors_.size())
    throw BadIndex("VariableBase::get_error",i,errors_.size());
  return errors_(i);
}
  
void RealVariable::reset(VariableBase& valuev)//, bool isconstv)
{
  isconstant(valuev.isconst()); 
  varp_=&valuev; 
}

double VariableBase::get_error() const
{
  if (errors_.empty())
    return 0.0;
  size_t offset=0U;
  if (issimple()) {
    if (errors_.size()!=1)
      throw InternalError("VariableBase::get_error");
  }
  else
    offset=getoffset().first;
  return get_error(offset);
}

 void RealVariable::set_error(double err, size_t i)
 {
   if (varp_==NMRSIM_NULL)
     throw InternalError("RealVariable::set_error");
   varp_->set_error(err,i);
 }

double RealVariable::get_error() const
{
  if (!varp_)
    return 0.0;
  validate();
  return varp_->get_error();
}

const BaseList<double> RealVariable::get_errors() const
{
  if (!varp_)
    return BaseList<double>();
  return varp_->get_errors();
}

bool RealVariable::validate() const
{
  return varp_ ? varp_->validate_context(name_.c_str()) : true;
}
 
const VariableBase& RealVariable::value() const
{
  if (varp_==NMRSIM_NULL)
    throw Undefined("RealVariable");
  return *varp_;
}

VariableBase& RealVariable::value()
{
  if (varp_==NMRSIM_NULL)
    throw Undefined("RealVariable");
  return *varp_;
}

std::ostream& dumparray(std::ostream& ostr, const BaseList<double>& a, char paren)
{
  if (a.size()==1)
    return ostr << a.front();

  ostr << paren;
  const size_t truncationlimit = ((verbose & VER_GEN) && (verbose_level>1)) ? 0 : DEFAULT_TRUNCATIONLIMIT;
  size_t actual=a.size();
  if (truncationlimit && (a.size() > truncationlimit+1))
    actual=truncationlimit;
  for (size_t i=0;i<actual;i++) {
    if (i)
      ostr << ',';
    ostr << a(i);
  }
  if (actual<a.size())
    ostr << ", " << (a.size()-actual) << " items omitted";
  return ostr << matchingbracket(paren);
}

int VariableBase::arrayuses() const
{
  int attr=0;
  for (size_t j=arrayexprs_.size();j--;) {
    const Expression* expr(arrayexprs_(j));
    if (expr)
      attr |= expr->uses();
  }
  return attr;
}

int VariableBase::uses() const 
{
  int attr=0;  
  if (isarray()) {
    if (array.size()>1)
      attr |= A_SUM;
    for (size_t i=array.size();i--;) {
      if (array(i).size()>1)
	attr |= A_ARRAY;
    }
    for (size_t j=errors_.size();j--;) {
      if (errors_(j))
	attr |= A_VAR;
    }
    attr |= arrayuses();
    return attr;
  }
  if (!(arrayexprs_.empty()))
    throw InternalError("VariableBase::uses");

  if (isexpr())
    return expr_.uses();

  if ((errors_.size()==1) && errors_.front())
    return A_VAR;

  return 0;
}

void VariableBase::print(std::ostream& ostr, bool full) const
{
  if (value_.size()==1)
    ostr << value_.front();
  else
    ostr << value_;    

  if (!full)
    return;

  if (isarray()) {
    ostr << ": ";
    const bool havesum=(array.size()>1);
    if (havesum)
      ostr << '|';
    for (size_t i=0;i<array.size();i++) {
      if (i)
	ostr << ',';
      const BaseList<double> a(array(i));
      const bool havearray=(a.size()>1);
      if (havearray)
	ostr << '{';
      const bool havevector=(arraylength_>1);
      bool needcomma=false;
      for (size_t j=0;j<a.size();j+=arraylength_) {
	if (havevector) {
	  ostr << '[';
	  needcomma=false;
	}
	for (size_t k=0;k<arraylength_;k++) {
	  if (needcomma)
	    ostr << ',';
	  else
	    needcomma=true;

	  const size_t offset=array.offset(i,j+k);
	  if (!(arrayexprs_.empty()) && arrayexprs_(offset))
	    ostr << (*arrayexprs_(offset));
	  else {
	    ostr << a(j+k);
	    const double err=errors_(offset);
	    if (err)
	      ostr << " +/- " << err;
	  }
	}
	if (havevector)
	  ostr << ']';
      }
      if (havearray) {
	ostr << '}';
	const size_t tag=arraytags_.empty() ? 0 : arraytags_(i);
	if (tag)
	  ostr << ':' << tag;
      }
    }
    if (havesum) {
      ostr << '|';
      if (sumtag_)
	ostr << ':' << sumtag_;
    }
    return;
  }
  if (isexpr()) {
    ostr << ' ';
    expr_.print(ostr);
    return;
  }
  if (errors_.size()==1) { //!< simple variable
    const double err=errors_.front();
    if (err)
      ostr << " +/- " << err;
  }
}

VariableBase& VariableBase::operator= (const BaseList<double>& vals)
{
  value_=vals;
  array.clear();
  expr_.clear();
  isconst_=true;
  return *this;
}

VariableBase& VariableBase::operator= (const Matrix<double>& vals)
{
  if (!vals)
    error_abort("Can't initialise variable to empty matrix");
  //  value_=vals;
  initialise(true); //!< initialise as constant
  arraylength_=vals.cols();
  array.create(1,vals.size());
  array.push_back(vals.row());
  update();
  return *this;
}

void VariableBase::initialise(size_t rows, size_t cols, bool isconstv)
{
  initialise(isconstv);
  const size_t els=rows*cols;
  if (els==0)
    error_abort("Can't initialise variable to empty matrix");
  arraylength_=cols;
  array.create(BaseList<size_t>(1,const_cast<size_t*>(&els)), 0.0);;
  arraytags_.clear();
  update();
}

void VariableBase::clear()
{
  value_.create(0U); //!< empty without realising memory
  array.clear();
  expr_.clear();
  errors_.clear();
}

VariableBase::VariableBase(double val, double err)
  : value_(1U,val),
    errors_(1U,err),
    isconst_(err==0.0) {
  if (err<0.0)
    throw InvalidParameter("VariableBase: error can't be <0");
}

void VariableBase::initialise(double initvalv, double errv)
{
  if (errv<0.0)
    throw InvalidParameter("VariableBase: error can't be <0");
  clear();
  value_.create(1,initvalv);
  errors_.create(1,errv);
  isconst_=(errv==0.0);
}

void VariableBase::initialise(Expression& expr)
{
  clear();
  expr_.swap(expr);
  isconst_=expr_.isconstant();
  if (isconst_) //!< only get value if const (may trigger error otherwise)
    update();
}

void VariableBase::initialise(bool isconstv) 
{
  clear();
  value_.clear();
  isconst_=isconstv;
}

void VariableBase::swapin(VariableBase& a) //NOT symmetrical (hence swapin)
{
  value_.swap(a.value_);
  isconst_=a.isconst_;
  if (a.isexpr())
    expr_.swap(a.expr_);
  else {
    if (a.isarray()) {
      array.swap(a.array);
      if (a.arraytags_.empty())
	arraytags_.clear();
      else
	arraytags_.swap(a.arraytags_);
      arrayexprs_.swap(a.arrayexprs_);
      exproffsets_.swap(a.exproffsets_);
      sumtag_=a.sumtag_;
      arraylength_=a.arraylength_;
    }
    errors_.swap(a.errors_);
  }
}

// double VariableBase::value() const
// {
//   assert(isexpr());
//   return expr_.get();
// }

void VariableBase::set(double val)
{
  value_.create(1,val);
}

// bool VariableBase::isconst() const
// {
//   if (isexpr())
//     return expr_.isconstant();
//   return !isarray();
// }

class ContextGuard {
public:
  ContextGuard(context_t newv)
    : oldstate_(evaluation_state) { evaluation_state=newv; }

  ~ContextGuard() { evaluation_state=oldstate_; }

private:
  context_t oldstate_;
};

void VariableBase::update(bool allowwarnings)
{
  ContextGuard guard(CONTEXT_UNCHECKED);

  if (isexpr())
    expr_.get(value_);
  else {
    if (isarray())
      value_=BaseList<double>(arraylength_, array.row().vector() ); //!< NB arrays (and expression arrays) are updated by updateindex - this just ensures initial value is valid
  }
}

void Expression::get(LIST<double>& res) const
{
  clearevalstacks();
  try {
    back()->get(res,*this);
  }
  catch (const MatrixException& c) {
    parser_printcontext() << "Error: " << c << "when evaluating: ";
    print(std::cerr);
    std::cerr << '\n';
    error_abort();
  }
  catch (const std::exception& c) {
    parser_printcontext() << "error '" << c.what() << "' evaluating: ";
    print(std::cerr);
    std::cerr << '\n';
    error_abort();
  }
}

shared_rule_t supervariable_rule,basefactor_rule;

//create parser

template expr_grammar::definition< phrase_scanner_t >::definition(const expr_grammar&);

expr_grammar& get_expr_parser()
{
  static expr_grammar expr_parser;
  return expr_parser;
}
