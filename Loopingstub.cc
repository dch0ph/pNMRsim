// stub for Looping functionality

#include "NMRsim_common.h"

//LIST<size_t> skips; //!< only used for reporting (better merged)
//LIST<int> ns(MAX_DIMENSIONS,0); //!< only used for reporting (better merged)
//LIST<double> sws(MAX_DIMENSIONS,0.0);

std::pair<size_t,bool> VariableBase::getoffset() const
{
  error_abort("support for arrays not enabled");
}

void ensure_array() { throw Failed("ensure_array should not be called in simple simulation"); }

bool getindex(size_t&, size_t, size_t)
{
  error_abort("support for arrays not enabled");
}
