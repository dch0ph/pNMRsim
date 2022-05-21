
#include "NMRsim.h"
#include "Parser.h"
#include "optim.h"

static int min_flags=0;

//framework for external minimisation function

// read any arguments
void parse_minimise_external()
{
}

void parse_minimise()
{
  
}

bool minimise_external(BaseMinFunction* funcp, varpars_t& minvarpars, bool invert)
{
  return minimise(funcp,minvarpars,invert);
}

void prepare_paras(LIST<double>& vals, LIST<double>& errs, varpars_t& lvarpars, int flags)
{
}

bool minimise(LIST<double>& min_vals,BaseMinFunction* funcp, varpars_t& minvarpars, bool invert)
{
  LIST<double> min_vals,min_errs;
  prepare_paras(min_vals,min_errs,minvarpars,min_flags);
  
  //do fitting

  set_variables(min_vals);
  return true;
}
