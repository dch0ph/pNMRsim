/*! \file
 \brief  Parsing for higher level objects
*/

#include "parser_common.hpp"
#include "NMRsim_RF.h"
#include "Parser.h"

using namespace boost::spirit;

//bool nondiagonal_opspec=false;
LIST<productoperator_spec*> opspecstack;

#include "assign_const_action.hpp"

//! builder for ::multioperator_spec
struct OperatorSpecBuilder {

  OperatorSpecBuilder(VariableBuilder& varbuilderv)
    : varbuilder_(varbuilderv) {}

  VariableBuilder& varbuilder_; //!< builder for coefficient expressions
  multioperator_spec multispec_; //!< destination object
  productoperator_spec spec_; //!< single product operator term
  LIST<operator_spec> tmp_; //!< temporary list of operator terms
  double scale_; //!< temporary scale factor
  bool isminus_,iscomplex_; //!< flags for i/-i imaginary qualifier

  char op; //!< operator identifier (xyz... or , for ST)
  char optype; //!< type of operator (I, F)
  unsigned int arg; //!< spin/nucleus
  std::string nucname; //!< temporary for nucleus name
  unsigned int r; //!< temporary for ST row
  unsigned int c; //!< temporary for ST col

  void reset(); //!< sub reset
  void reset(int flags); //!< reset builder and parse \a flags
  void add(); //!< add term
  void addlist(); //!< add list of terms
  void addproduct(); //!< add complete ::productoperator_spec

  struct Add;
  struct AddList;
  struct AddProduct;
  struct Reset;
};
  
//! Spirit grammar for complex objects
struct ext_grammar : public grammar<ext_grammar>
{
  ext_grammar(expr_grammar&);
  
  //! labels for rules
  enum { productopspec_def=0, multiopspec_def=1 };

  template <typename ScannerT> struct definition
    : public grammar_def< rule<ScannerT>, rule<ScannerT> >
  {
    typedef rule<ScannerT> rule_t;

    rule< typename lexeme_scanner<ScannerT>::type > opspec;
    rule_t opterm,productopspec,multiopspec;

    definition(const ext_grammar&);
  };

  void reset(int flags =0); //!< reset with parse \a flags

  mutable OperatorSpecBuilder opbuilder_; //!< ::multioperator_spec builder
  multioperator_spec& setableoperator() { return opbuilder_.multispec_; } //!< return constructed ::multioperator_spec
  productoperator_spec& productoperator() { return opbuilder_.spec_; } //!< return last ::productoroperator_spec created
};
  
struct OperatorSpecBuilder::Reset
{
  explicit Reset(OperatorSpecBuilder& builderv)
    : builder_(builderv) {}

  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.reset();
  }
  OperatorSpecBuilder& builder_;
};


template <typename ScannerT> ext_grammar::definition<ScannerT>::definition(const ext_grammar& self)
{  
  /**
     \note Not perfect e.g. will accept 13C::x */
  opspec = lexeme_d[
		    (( ch_p('I') >> ( ch_p('n')[assign_const_a(self.opbuilder_.optype,'F')] | uint_p[assign_const_a(self.opbuilder_.optype,'I')][assign_a(self.opbuilder_.arg)]) )
 		     | (ch_p('C')[assign_a(self.opbuilder_.optype)] >> uint_p[assign_a(self.opbuilder_.arg)])
 		     | ch_p('F')[assign_a(self.opbuilder_.optype)]
 		     | ((+alnum_p)[assign_a(self.opbuilder_.nucname)]
			>> ch_p(':')[assign_a(self.opbuilder_.optype)])
 		     )
		    >> ( 
		    chset<>("xyzpmc+-")[assign_a(self.opbuilder_.op)]    
		    	|  (!ch_p(':') 
			   >> uint_p[assign_a(self.opbuilder_.r)]
			      >> ch_p(',')[assign_a(self.opbuilder_.op)]
			    >> uint_p[assign_a(self.opbuilder_.c)]))
 		    ];
  
  opterm =
    ( !(basefactor_rule >> ch_p('*') >> eps_p[VariableBuilder::Cleanup<VariableBuilder::EXPR>(self.opbuilder_.varbuilder_)])
      >> eps_p[OperatorSpecBuilder::Reset(self.opbuilder_)] 
      >> !( (!(ch_p('-')[assign_const_a(self.opbuilder_.isminus_,true)]) >> ch_p('i') >> ch_p('*'))[assign_const_a(self.opbuilder_.iscomplex_,true)])
      >> lexeme_d[opspec][OperatorSpecBuilder::Add(self.opbuilder_)]
	 >> *(  (ch_p('*') >> lexeme_d[opspec][OperatorSpecBuilder::Add(self.opbuilder_)]) ))[OperatorSpecBuilder::AddList(self.opbuilder_)];
  
  productopspec
    = opterm >>
    // - sign is acknowledged but parsed later
    *(  (ch_p('+') | eps_p('-')) >> opterm);

  multiopspec
    = (ch_p('{') >> list_p( productopspec[OperatorSpecBuilder::AddProduct(self.opbuilder_)], ch_p(','),'}') >> ch_p('}'))
    | productopspec[OperatorSpecBuilder::AddProduct(self.opbuilder_)];

  this->start_parsers(productopspec,multiopspec);
}

void OperatorSpecBuilder::reset()
{
  tmp_.clear();
  isminus_=false;
  iscomplex_=false;
}

struct OperatorSpecBuilder::Add
{ // add element
  explicit Add(OperatorSpecBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.add();
  }
  OperatorSpecBuilder& builder_;
};

struct OperatorSpecBuilder::AddList
{ // add productop term
  explicit AddList(OperatorSpecBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.addlist();
  }
  OperatorSpecBuilder& builder_;
};

struct OperatorSpecBuilder::AddProduct
{ // add productop term
  explicit AddProduct(OperatorSpecBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.addproduct();
  }
  OperatorSpecBuilder& builder_;
};

static bool opspec_const;
static LIST< std::pair<size_t,size_t> > subsidspecs;

double handle_operator_variable(int flags, size_t listoffset, size_t offset)
{
  subsidspecs.push_back(std::pair<size_t,size_t>(listoffset,offset));
  opspec_const=false;
  return handle_variable(flags,subsidspecs.size());
}

size_t parse_nucleusname(const char* label)
{  
  size_t nuc;
  try {
    nuc=labeltonuc(label);
  }
  catch (Failed&) {
    parser_printcontext();
    if (*label)
      std::cerr << "Nucleus name not recognised: " << label << " (should be 13C, 2H etc.)\n";
    else
      std::cerr << "Missing nucleus name\n";
    error_abort();
  }
  return nuc;
}

void OperatorSpecBuilder::add()
{
  if (tmp_.empty()) {
    if (varbuilder_.isbuilt())
      scale_=handle_operator_variable(varbuilder_.flags,multispec_.size(),spec_.size());
    else
      scale_=1.0;
  }
  
  switch (op) {
  case 'm':
    op='-';
    break;
  case 'p':
    op='+';
    break;
  case ',':
    op='\0';
    break;
  }
//   switch (optype) {
// //   case 'n': //convert In -> F
// //     optype='F';
// //     break;
//   case ':':
//     break;
//   }

  const basespin_system* sysp(get_spin_system());

  size_t nuc;
  size_t ind=-1;

  switch (optype) {
  case 'I':
    if (!sysp)
      error_abort("can't use I operators without spin system");
    if ((arg<1) || (arg>sysp->nspins()))
      error_abort("spin index out of range");
    ind=arg-1;
    nuc=(*sysp)(ind).nucleus();
    break;
  case 'F': 
    if (sysp ? !sysp->ishomonuclear() : (nucids.size()>1))
      error_abort("Can't use F<op>/In<op> with heteronuclear spin system");
    nuc=sysp ? (*sysp)(0U).nucleus() : nucids.front();
    break;
  case 'C':
    if ((arg<1) || (arg>nchannels))
      error_abort("channel index out of range");
    nuc=nucids(arg-1);
    break;
  case ':': 
    nuc=parse_nucleusname(nucname.c_str());
    break;
  default:
    throw InternalError("make_operator_spec");
  }

  if (!op) {
    const size_t deg=spin(nuc).deg();
    if ((r<1) || (c<1))
      error_abort("invalid single transition operator specification");
    if ((r>deg) || (c>deg))
      error_abort("single transition operator specification invalid for this nucleus");
    r--;
    c--;
  }

  if (optype=='I') 
    tmp_.push_back(op ? operator_spec(ind,op) : operator_spec(ind,r,c));
  else
    tmp_.push_back(op ? operator_spec(nuclei_spec(nuc),op) : operator_spec(nuclei_spec(nuc),r,c));
}

void OperatorSpecBuilder::addproduct()
{
  if (!spec_)
    throw InternalError("addproduct");
//   if (!(spec_.issimple()  ||  spec_.isdiagonal()))
//     nondiagonal_opspec=true;
  multispec_.push_back(productoperator_spec());
  productoperator_spec& curspec(multispec_.back());  
  curspec.swap(spec_);
  opspecstack.push_back(&curspec);
}
   
void OperatorSpecBuilder::addlist()
{
  const basespin_system* sysp=get_spin_system();
  if (sysp) {
    complex scale(scale_);
    if (iscomplex_)
      scale=complex(0.0,isminus_ ? -scale_ : scale_);
    spec_.push_back(scale,tmp_,*sysp);
  }
  else {
    if (tmp_.size()>1)
      throw Failed("can't use non-trivial product operator here");
    spec_.push_back(complex(scale_),tmp_.front());
  }
  tmp_.clear();
  scale_=1.0;
}

void OperatorSpecBuilder::reset(int flags)
{
  varbuilder_.reset(flags | F_ISPHASE);
  multispec_.clear();
  spec_.clear();
  scale_=1.0;
}

ext_grammar::ext_grammar(expr_grammar& exprgram)
  :
  opbuilder_(exprgram.builder())
{
  if (!exprgram.isinit) //fudge to force initialisation of expr_grammar rules
    (void)::parse("",exprgram.use_parser<expr_grammar::variable_def>(),space_p);
}

template ext_grammar::definition< phrase_scanner_t >::definition(const ext_grammar&);

ext_grammar& get_ext_parser()
{
  static ext_grammar ext_parser(get_expr_parser());
  return ext_parser;
}

void ext_grammar::reset(int flags)
{ 
  //  cyclebuilder_.reset(flags);
  opbuilder_.reset(flags);
}

setableoperator_spec* create_setableoperator(char* curptr, int flags)
{
  ext_grammar& ext_parser(get_ext_parser());
  ext_parser.reset(flags);
  opspec_const=true;
  Mark markobj;
  const parse_info<> info(parse(curptr,
				ext_parser.use_parser<ext_grammar::multiopspec_def>(),
				space_p));
  if (!info.full) {  
    parser_printcontext() << "parsing of product operator list specification (\"" << curptr << "\") failed at: \"" << info.stop << '\"' << std::endl;
    error_abort();
  }
  setableoperator_spec* resultp=new setableoperator_spec(opspec_const);
  resultp->swap(ext_parser.setableoperator());
  const size_t n=resultp->size();
  if (n>1)
    array_dims.set(0,n); //register size
  markobj.flush(resultp);
  return resultp;
}

productoperator_spec* create_productoperator(char* curptr, int flags)
{
  ext_grammar& ext_parser(get_ext_parser());
  ext_parser.reset(flags | F_DENYVAR | F_DENYARRAY);
  const parse_info<> info(parse(curptr,
				ext_parser.use_parser<ext_grammar::productopspec_def>(),
				space_p));
  if (!info.full) {  
    parser_printcontext() << "parsing of product operator specification failed at: \"" << info.stop << "\"\n";
    error_abort();
  }
  productoperator_spec* resultp=new productoperator_spec();
  resultp->swap(ext_parser.productoperator());
  //  if (!(resultp->issimple()  ||  resultp->isdiagonal()))
  //  nondiagonal_opspec=true;
  opspecstack.push_back(resultp);
  return resultp;
}

productoperator_spec* parse_productoperator(int flags)
{
  productoperator_spec* result=create_productoperator(get_curline(),flags);
  set_curline(NMRSIM_NULL);
  return result;
}

setableoperator_spec* parse_setableoperator(int flags)
{
  setableoperator_spec* result=create_setableoperator(get_curline(),flags);
  set_curline(NMRSIM_NULL);
  return result;
}

void setableoperator_spec::set(double val, subsid_t subsid)
{
  size_t which=size_t(subsid)-1;
  const std::pair<size_t,size_t>& offsets(subsidspecs(which));
  productoperator_spec& curspec((*this)(offsets.first));
  curspec.scales(offsets.second)=val;
}

void setableoperator_spec::print(std::ostream& ostr) const
{
  ostr << (*this);
}

//NB unique only within object
void setableoperator_spec::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  size_t which=size_t(subsid)-1;
  const std::pair<size_t,size_t>& offsets(subsidspecs(which));
  if (size()>1)
    ostr << 'i' << (offsets.first+1) << '_';
  ostr << 't' << (offsets.second+1);
//   const productoperator_spec& curspec((*this)(offsets.first));
//   const complex scale=curspec.scales(offsets.second);
//   static const complex one(1.0);
//   if (scale!=one)
//     ostr << scale << ' ';
//   ostr << curspec.specs(offsets.second);
}
