#include <iostream>
#include "ListList.h"
//#include "parser_definition.hpp"
#include "parser_common.hpp"

int sum_index=0;

using namespace std;

std::ostream& parser_printcontext(std::ostream& ostr) { return ostr; }

char* substitute_string(char* out, int n, const char* in, bool, bool)
{
  strncpy(out,in,n);
  return out;
}

const size_t verbose_level=2;

size_t parser_verbose_level()
{
  return verbose_level;
}

bool proton_freq_isconstant()
{
  return true;
}

double curgrat=1.0;
double proton_freq=500e6;

static const spin_system sys(2,"1H");
size_t nchannels=0;
LIST<size_t> nucids;

const basespin_system& get_spin_system()
{
  return sys;
}

// operator_spec make_operator_spec(char optype, size_t arg, char op, size_t r, size_t c)
// {    
//   if (optype=='I')
//     return op ? operator_spec(size_t(arg),op) : operator_spec(size_t(arg),r,c);

//   nuclei_spec which;
//   switch (optype) {
//   case 'F':
//     which=nuclei_spec("1H");
//     break;
//   case 'C':
//     which=nuclei_spec(arg);
//     break;

//     return operator_spec(nuclei_spec(arg),op);
//   }
//   throw InternalError("make_operator_spec");
// }

  
// struct varval {
//   double value;
//   double step;   
// };

// static varval varval_tmp;

// struct varval_set_value {
//   void operator()(double val) const {
//     varval_tmp.value=val;
//     varval_tmp.step=-1.0;
//   }
// };

// struct varval_isvariable {
//   void operator()(char) const {
//     varval_tmp.step=0.0;
//   }
// };

// struct varval_set_step {
//   void operator()(double val) const {
//     varval_tmp.step=val;
//   }
// };

void fudge_shift(double& value) //500 MHz
{
  value*=proton_freq*curgrat*1e-6;
}

struct Variable : public RealVariable {
  Variable(double valv) : RealVariable(false,true), val_(valv) {}
  const BaseList<double> get_list() const { return BaseList<double>(1,const_cast<double*>(&val_)); }
  void print(std::ostream& ostr) const { ostr << val_; }
  double val_;
};

Variable v_pi(M_PI);  

//only variable recognised is pi!
RealVariable& parse_variable_name(const char* varname)
{
  if (strcmp(varname,"pi")==0)
    return v_pi;
  throw Failed("Variable name not recognised");
}


// namespace {
//   double var_step;

//   struct set_step {
//     void operator()(double val) const {
//       var_step=std::fabs(val);
//     }
//   };
// }

struct PhasedSequence
{
  std::string name_;
  double phase_;
  PhasedSequence(const char* name, double phase) : name_(name), phase_(phase) {}
};

std::ostream& operator<< (std::ostream& ostr, const PhasedSequence& a)
{
  return ostr << a.name_ << '+' << a.phase_;
}

double handle_variable(int,size_t)
{
  return get_expr_parser().variable().value().front();
}

double handle_operator_variable(int flags, size_t listoffset, size_t offset)
{
  return handle_variable(flags,offset);
}

PhasedSequence* handle_phasedsequence(int flags)
{
    return new PhasedSequence(
			      get_ext_parser().laststring().c_str(),
			      handle_variable(flags,1));
}

int main()
{
  const int flags=F_ISPHASE | F_ISSHIFT | F_ALLOWLIST;
  expr_grammar& expr_parser(get_expr_parser());
  ext_grammar& ext_parser(get_ext_parser());

  string str;
  while (getline(cin, str)) {
    if (str.empty() || tolower(str[0]) == 'q')
      break;

    const char* input = str.c_str();   
    expr_parser.reset(flags);
    try {
      if (strncmp(input,"S:",2)==0) {
	ext_parser.reset(flags);
	std::cout << "Parsing as cycled sequence\n";
	const parse_info<> info(parse(input+2,
				      ext_parser.use_parser<ext_grammar::cycledseq_def>(),
				      space_p));
	if (info.full) {
	  cout << "Parsing succeeded\n";
	  cout << ext_parser.cycledsequence() << '\n';
	}
	else
	  cout << "Parsing failed at: \"" << info.stop << "\"\n";
      }
      else {
	if (strncmp(input,"O:",2)==0) {
	  std::cout  << "Parsing as product operator\n";
	  ext_parser.reset(flags);
	  const parse_info<> info(parse(input+2,
					ext_parser.use_parser<ext_grammar::productopspec_def>(),
					space_p));
	  if (info.full) {
	    cout << "Parsing succeeded\n";
	    cout << ext_parser.productoperator() << '\n';
	  }
	  else
	    cout << "Parsing failed at: \"" << info.stop << "\"\n";
	}
	else {
	  std::cout << "Parsing as variable\n";
	  const parse_info<> info(parse(input,
					expr_parser.use_parser<expr_grammar::variable_def>(),
					space_p));
	  if (info.full) {
	    cout << "Parsing succeeded\n";
	    expr_parser.variable().print(cout);
	    cout << '\n';
	  }
	  else
	    cout << "Parsing failed at: \"" << info.stop << "\"\n";
	}
      }
    } catch (MatrixException& exc) {
      std::cerr << exc << '\n';
    }
  }
  return 0;
}


//     cout << "Trying to parse as simple variable\n";
//     if (parse(input, var_rule).full)
//       cout << "Value: " << value << "  Step: " << var_step << '\n'; 
//     else {   
//       cout << "Trying to parse as array\n";

// 	if ((*input=='(') || (*input=='{')) {
// 	  // try to parse as array
// 	  arraybuilder.reset(flags);
// 	  parse_info<> info(parse(input,
// 			    expr_parser.use_parser<expr_grammar::constarray_def>(),
// 				  space_p));
// 	  if (info.full)
// 	    cout << "Parsing succeeded: " << curarray << '\n';
// 	  else
// 	    cout << "Parsing failed at: \"" << info.stop << "\"\n";
// 	}
// 	else {
// 	  cout << "Trying to parse as expression\n";

//   double value=0.0;
//   var_step=-1.0;
//   bool valisp=false;
//   bool stepisp=false;

//   static rule<> var_rule = 
//     lexeme_d[real_p[assign_a(value)] >> !ch_p('p')[assign_a(valisp,true)] ]
//     || ch_p('V')[assign_a(var_step,0.0)]
//     || lexeme_d[ureal_p[set_step()] >> !ch_p('p')[assign_a(stepisp,true)] ];

  //  rule<> var_rule = real_p[varval_set_value()] || ch_p('V')[varval_isvariable()] || ureal_p[varval_set_step()];
