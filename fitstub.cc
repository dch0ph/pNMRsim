// sub for programs not requiring fitting

#include "NMRsim.h"

LIST<Parameter> parameter_list;

ThreadWarning<> notfitting_warning("Creation of fitting variation is being ignored (fitting support not compiled in",&NMRsim_once_warning);
 
void register_fitting_variable(VarVariable*)
{
  notfitting_warning.raise();
}

bool try_find_sw() { return false; }
bool try_find_np() { return false; }

DataStore fit_set;
