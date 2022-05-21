#ifndef parser_common_hpp_
#define parser_common_hpp_

/*! \file
 \brief  Common functionality for parsing
*/

//! Flag to Spirit that grammars are instantiated only once
#define BOOST_SPIRIT_SINGLE_GRAMMAR_INSTANCE
//! Increase number of parse rules (default is 3)
#define PHOENIX_LIMIT 4
#define BOOST_SPIRIT_GRAMMAR_STARTRULE_TYPE_LIMIT PHOENIX_LIMIT
//#define BOOST_SPIRIT_DEBUG

//! limits.h needed for INT_MAX (problem in Spirit 1.8.5?)
#include <limits.h>
#include <boost/spirit.hpp>
#include <boost/spirit/dynamic/stored_rule.hpp>
#include <boost/spirit/utility/grammar_def.hpp> 

#include "NMRsim_common.h"

using namespace boost::spirit;

//shared rules

typedef rule<phrase_scanner_t> shared_rule_t; //!< type for shared rules

extern shared_rule_t basefactor_rule; //!< rule for multiplicative factor
extern shared_rule_t supervariable_rule; //!< rule for variable/expression

typedef assertion<const char*> parser_assert_t; //!< type for parsing assertions
typedef const char* def_iteratorT; //!< normal input iterator type

//! class for storing parser state at function/list start
template<typename IteratorT =def_iteratorT> struct MarkState_t {
  IteratorT begin; //!< start of name
  IteratorT end; //!< end of name
  size_t stackpos;
  size_t nargs;

  MarkState_t(size_t stackposv)
    : begin(NMRSIM_NULL), end(NMRSIM_NULL), stackpos(stackposv), nargs(-1) {} //!< list (no function name)

  MarkState_t(IteratorT beginv, IteratorT endv, size_t stackposv)
    : begin(beginv), end(endv), stackpos(stackposv), nargs(-1) {} //!< function
};

//! builder object for virtual dimension set
struct TagBuilder {
  bool hastagspec_; //!< \c true if tags have been specified
  size_t sumtag_; //!< sum array tag
  size_t arraytag_; //!< temporary array tag
  LIST<size_t> arraytags_; //!< array tags
  void reset(); //!< clear
  void flushlist(); //!< add new array tag
};

//! builder class for ::VariableBase
struct VariableBuilder {

  VariableBuilder(int flagsv =0) : flags(flagsv), built(false) {}
  ~VariableBuilder();

  VariableBase cvar_; //!< scratch ::VariableBase
  Expression expr_; //!< scratch expression
  int flags; //!< parse flags
  bool built; //!< \c true if variable has been constructed
  LIST<size_t> curstack_; //!< expression mark stack
  LIST< MarkState_t<> > markstack_; //!< function/list start stack
  ListList<double> array_; //!< array variable values
  LIST<double> steps_; //!< error steps
  LIST<Expression*> exprstmp_;
  LIST<size_t> offsetstmp_;
  bool haveexprs_;
  LIST<double> scr_; 
  LIST<double> tmp_;
  TagBuilder tagbuilder_; //!< virtual dimension builder
  LIST<std::string> stringstack_; //!< string arguments - no stack required since not used recursively
  
  //!< variable types
  enum var_t { SIMPLE, //!< scalar constant
	       ARRAY, //!< sum/row array
	       EXPR //!< expression
  };

  double var_value; //!< scalar value temporary
  double var_step; //!< error temporary
  bool valisp; //!< p indicator
  char stepqual; //!< qualifier for error (0 if none)
  bool isint; //!< \c true if integral
  unsigned int ppmqual; //!< shift nucleus qualifier
  
  static parser_assert_t stack_empty; //!< expression mark stack empty assertion
  //  static ContextWarning<> nonconst_in_array_warning;
  static ContextWarning<> zerosizeexpr_in_array_warning;

  void reset(int flags); //!< restart with given input \a flags
  void reset(); //!< reset (without changing flags)
  void resetexpr(); //!< restart expression parsing
  void resetsimple(); //!< restart simple scalar 
  void flushlist();
  double rawcleanup(bool denyexpr); //!< \internal finish variable parsing - return var_step
  bool checkclean(); //!< check stacks are clear
  
  VariableBase& operator()() { return cvar_; } //!< return finished ::VariableBase
  const VariableBase& operator()() const { return cvar_; }

  bool isbuilt() const { return built; } //!< return \c true if build complete

  bool push(ExpressionBase*, bool constcheck =true); //!< add node to expression (return \c true if const)
  size_t pop(); //!< pop expression mark stack
  ExpressionBase* back() { return const_cast<ExpressionBase*>(expr_.back().get()); } //!< last (root) node
  const ExpressionBase* back() const { return expr_.back().get(); } 

  void mark_push() {
    markstack_.push_back(MarkState_t<>(curstack_.size()));
  } //!< add (list) to function/list mark stack
  void mark_push(def_iteratorT begin, def_iteratorT end) {
    markstack_.push_back(MarkState_t<>(begin,end,curstack_.size()));
  } //!< add (function) to function/list mark stack

  const MarkState_t<>& mark_back(); //!< last mark
  void mark_pop() { markstack_.pop_back(); } //!< pop (and discard) mark stack

  bool endsequence(ExpressionBase*); //!< transfer stored end to new object, pop mark state and add to expression stack  
  void endfunction(bool isfunc); //!< finish function/list
  void endstringfunction(); //!< end string function
  void fudge_shift(size_t =0); //!< ppm -> freq
  double checkV(double); //!< process var_step

  void evaluate(LIST<double>&); //!< evaluate expression

  void add(); //!< add entry to {}
  void addlist(); //!< add complete {}
  void addsumlist(); //!< add entry to []
  void addsimple(); //!< add simple scalar

  struct set_step; //!< set error step functor
  struct ResetSimple;
  struct Add;
  struct AddSimple;
  struct AddList;
  struct AddSumList;
  struct AddSlot;
  struct Mismatched;
  struct FailedTag;
  struct FailedSlot;

  template<var_t Type> void cleanup(); //!< cleanup for given type
  template<var_t Type> struct Cleanup; //!< cleanup functor
  template<var_t Type> void fullcleanup() {
    cleanup<Type>();
    checkclean();
    built=true;
  }   
};

template<VariableBuilder::var_t Type> struct VariableBuilder::Cleanup
{
  Cleanup(VariableBuilder& builderv) : builder_(builderv) {}
  template<typename IteratorT> void operator()(IteratorT, IteratorT) const { builder_.fullcleanup<Type>(); }
  VariableBuilder& builder_;
};

//! Spirit grammar for expressions / variables
struct expr_grammar : public grammar<expr_grammar>
{
  expr_grammar() : isinit(false) {}

  //! labels for rules
  enum { variable_def=0,
	 rawlist_def=1,
	 expression_def=2
  };

  template <typename ScannerT> struct definition
    : public grammar_def< rule<ScannerT>, rule<ScannerT>, rule<ScannerT> >
  {
    typedef rule<ScannerT> rule_t;

    rule< typename lexeme_scanner<ScannerT>::type > variable,varname,slot,stepqualdef,Vdef,baseppmrule,ppmrule;
    rule_t expression,factor,function,term,plainlist,staticinclude,dynamicinclude,listel,varlistel,range,varrange,stringliteral,baseterm;
    rule_t arraylist,sumlistel,sumlist,constarray,rawlist,simplevariable,stepqual,sumtag,phasedef,arraytag;

    definition(const expr_grammar&);
  };

  void reset(int flags) { builder_.reset(flags); } //!< reset

  mutable bool isinit; //!< \c true if grammar has been initialised
  mutable VariableBuilder builder_; //!< ::VariableBase builder
  VariableBase& variable() { return builder_.cvar_; } //!< return contructed ::VariableBase
  Expression& expression() { return builder_.expr_; } //!< return temporary expression
  VariableBuilder& builder() { return builder_; }
};

expr_grammar& get_expr_parser(); //!< return base expression grammar object

#endif
