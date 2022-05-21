/*! \file
 \brief  Parsing for objects specific to RF
*/

#include "parser_common.hpp"
#include "NMRsim_RF.h"
#include "Parser.h"

//! builder for ::CycledSequence
struct CycledSeqBuilder {
  CycledSeqBuilder(VariableBuilder& varbuilderv) : varbuilder_(varbuilderv) {}

  VariableBuilder& varbuilder_; //!< builder for phase expressions
  CycledSequence cseq_; //!< destination object

  void reset(int flags); //!< reset builder and parse \a flags
  void add(); //!< add term
  void addlist(); //!< add {} list
  void flush(); //!< finish ::CycledSequence

  struct Add;
  struct AddList;
  struct Flush;

  LIST<PhasedSequence*> tmp_; //!< temporary for {} lists
  size_t arraytag_;
  //  TagBuilder tagbuilder_; //!< virtual dimension builder
};

struct RF_grammar : public grammar<RF_grammar>
{
  RF_grammar(expr_grammar&);
  
  //! labels for rules
  enum { cycledseq_def=0 };

  template <typename ScannerT> struct definition
    : public grammar_def< rule<ScannerT> >
  {
    typedef rule<ScannerT> rule_t;

    rule< typename lexeme_scanner<ScannerT>::type > varname;
    rule_t arraytag,sumtag,phasedseq,seqname,arrayseqlist,cycledseqlist;

    definition(const RF_grammar&);
  };

  mutable std::string laststring_;
  const std::string& laststring() const { return laststring_; }

  void reset(int flags =0); //!< reset with parse \a flags

  mutable CycledSeqBuilder cyclebuilder_; //!< ::CycledSequence builder
  CycledSequence& cycledsequence() { return cyclebuilder_.cseq_; } //!< return constructed ::CycledSequence
};

using namespace boost::spirit;

#include "assign_const_action.hpp"

template <typename ScannerT> RF_grammar::definition<ScannerT>::definition(const RF_grammar& self)
{
  varname
    =  lexeme_d[ +(alnum_p | '_') ];

  seqname = lexeme_d[varname][assign_a(self.laststring_)];

  phasedseq
    = (supervariable_rule >> ch_p('+') >> seqname)
    //NB '-' before term is only checked, parsing is done in term
    | (seqname >> !( (ch_p('+') | eps_p('-')) >> supervariable_rule));
  
  arraytag =
    (ch_p(':') >> uint_p[assign_a(self.cyclebuilder_.arraytag_)]) 
    | eps_p[assign_const_a(self.cyclebuilder_.arraytag_,0)];
    
 //  sumtag = 
//     ch_p(':') >> uint_p[assign_a(self.cyclebuilder_.tagbuilder_.sumtag_)];

//   arrayseqlist
//     = (ch_p('{') >> list_p( phasedseq[CycledSeqBuilder::Add(self.cyclebuilder_,false)], ch_p(','),'}') >> ch_p('}') >> arraytag)
//     | (ch_p('[') >> list_p( phasedseq[CycledSeqBuilder::Add(self.cyclebuilder_,true)], ch_p(','),']') >> ch_p(']'))
//     | phasedseq[CycledSeqBuilder::Add(self.cyclebuilder_,false)]; //default for single entry is {} 
  
//   cycledseqlist
//      = ((ch_p('(') >> list_p( arrayseqlist[CycledSeqBuilder::AddList(self.cyclebuilder_)], ch_p(','),')') >> ch_p(')') >> !sumtag)
//         | arrayseqlist[CycledSeqBuilder::AddList(self.cyclebuilder_)])[CycledSeqBuilder::Flush(self.cyclebuilder_)];

  arrayseqlist
    = (ch_p('[') >> list_p( phasedseq[CycledSeqBuilder::Add(self.cyclebuilder_)], ch_p(','),']') >> ch_p(']'))
    | phasedseq[CycledSeqBuilder::Add(self.cyclebuilder_)];
  
  cycledseqlist
     = ((ch_p('{') >> list_p( arrayseqlist[CycledSeqBuilder::AddList(self.cyclebuilder_)], ch_p(','),'}') >> ch_p('}') >> arraytag)
        | arrayseqlist[CycledSeqBuilder::AddList(self.cyclebuilder_)])[CycledSeqBuilder::Flush(self.cyclebuilder_)];
    
  this->start_parsers(cycledseqlist);
}

void RF_grammar::reset(int flags)
{ 
  cyclebuilder_.reset(flags);
}
 
void CycledSeqBuilder::reset(int flags) 
{ 
  varbuilder_.reset(flags | F_ISPHASE);
  cseq_.clear();
  //  tagbuilder_.reset();
  tmp_.create(0);
}

void CycledSeqBuilder::add()
{
  const subsid_t subsid=CycledSequence::index_to_subsid(cseq_.row().size()+tmp_.size());
  tmp_.push_back(handle_phasedsequence(varbuilder_.flags,subsid));
//   if (cseq_.empty())
//     cseq_.isspecial_=isspec;
//   else {
//     if (cseq_.isspecial_!=isspec)
//       throw Failed("Can't mix [] and {} sequence lists");
//   }
  varbuilder_.reset();
}

void CycledSeqBuilder::addlist()
{
  if (tmp_.empty())
    throw Failed("Empty sequence list");
  cseq_.push_back(tmp_);
//   if (!cseq_.isspecial_)//tag only valid for arrayed ("non special") list
//     tagbuilder_.flushlist();
  tmp_.create(0);
}

struct CycledSeqBuilder::Add
{ // add element
  Add(CycledSeqBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.add();
  }
  CycledSeqBuilder& builder_;
  //  bool isspecial_;
};

struct CycledSeqBuilder::AddList
{ // add list
  AddList(CycledSeqBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.addlist();
  }
  CycledSeqBuilder& builder_;
};

void CycledSeqBuilder::flush()
{
  cseq_.arraytag_=arraytag_;
  //  cseq_.arraytags_.swap(tagbuilder_.arraytags_);
}

struct CycledSeqBuilder::Flush
{
  Flush(CycledSeqBuilder& builderv)
    : builder_(builderv) {}
  
  template<typename IteratorT> void operator()(const IteratorT&, const IteratorT&) const {
    builder_.flush();
  }
  CycledSeqBuilder& builder_;
};

RF_grammar::RF_grammar(expr_grammar& exprgram)
  :
  cyclebuilder_(exprgram.builder())
{
  if (!exprgram.isinit) //fudge to force initialisation of expr_grammar rules
    (void)::parse("",exprgram.use_parser<expr_grammar::variable_def>(),space_p);
}

template RF_grammar::definition< phrase_scanner_t >::definition(const RF_grammar&);

RF_grammar& get_RF_parser()
{
  static RF_grammar RF_parser(get_expr_parser());
  return RF_parser;
}

PhasedSequence* handle_phasedsequence(int flags, subsid_t subsid)
{
  const std::string& namestr(get_RF_parser().laststring());
  CompSequenceBase& seq(*findmap(seqmap,namestr.c_str(),"sequence"));
  return new PhasedSequence(seq,
			    handle_variable(flags,subsid));
}

// PhasedSequence* parse_phasedsequence(char* tptr, int flags)
// {
//   expr_grammar& expr_parser(get_expr_parser());
//   expr_parser.reset(F_ISPHASE | flags);
//   const parse_info<> info(parse(tptr,
// 				get_RF_parser().use_parser<RF_grammar::phasedseq_def>(),
// 			  space_p));

//   if (!info.full) {
//     parser_printcontext() << "parsing of phased sequence failed at: \"" << info.stop << "\"\n";
//     error_abort();
//   }
//   return handle_phasedsequence(flags,S_PHASE);
// }

// PhasedSequence* parse_phasedsequence(int flags)
// {
//   char* cptr=parse_string(flags);
//   return cptr ? parse_phasedsequence(cptr,flags) : NMRSIM_NULL
// }

CycledSequence* parse_cycledsequence(int flags)
{
  return parse_cycledsequence(get_token(0,"sequence specification"),flags);
}

CycledSequence* parse_cycledsequence(char* ptr, int flags)
{
  const cycledseqmap_type::const_iterator curp(cycledseqmap.find(ptr)); // look first for simple CycledSequence with this name
  if (curp!=cycledseqmap.end())
    return curp->second;

  //! if using "smart propagator" optimisation, automatically promote simple sequences to cycled list
  if (optsmartprop()) {
    const seqmap_type::const_iterator curp2(seqmap.find(ptr));
    if (curp2!=seqmap.end()) {
      CompSequenceBase& seq(*(curp2->second));
      if ((verbose & VER_GEN) && (verbose_level>1))
	std::cout << "Creating new sequence list from sequence fragment " << seq.name() << '\n';
      CycledSequence* seqp=new CycledSequence(seq); //create CycledSequence based on simple sequence fragment
      cycledseqmap[seq.name()]=seqp;
      cycledseqlist.push_back(seqp);
      return seqp;
    }
  }

  // new cycled sequence
  RF_grammar& RF_parser(get_RF_parser());
  RF_parser.reset(F_ISPHASE | flags);
  const parse_info<> info(parse(ptr,
			  RF_parser.use_parser<RF_grammar::cycledseq_def>(),
			  space_p));
  if (!info.full) {
    if (flags & F_ALLOWMISSING)
      return NMRSIM_NULL;
    parser_printcontext() << "parsing of cycled sequence failed at: \"" << info.stop << "\"\n";
    error_abort();  
  }
  CycledSequence* seqp=new CycledSequence(RF_parser.cycledsequence());
  seqp->reset(true); //!< ensure that states_ has been created
  //  verify_arraysizes(seqp->aslist(),seqp->sumtag(),seqp->arraytags()); //"register" array sizes
  verify_arraysize(seqp->aslist().size(),seqp->arraytag());
  cycledseqlist.push_back(seqp);
  return seqp;
}
