/* stub functions if no RF, no pulseq block, no fitting variables or powder arrays */

#include "NMRsim_MasterObj.h"

void reset_pgens() {}
void MasterObj::ensure_channels() {}
//bool hasacquisition_rf() { return false; }
//Variable* create_fitting_variable(const VarVariable&) { throw Failed("create_fitting_variable should not be called in simple simulation"); }
template<class T> void MasterObj::set_Hsystem(const BlockedMatrix<T>&) {}
template void MasterObj::set_Hsystem(const BlockedMatrix<complex>&);
template void MasterObj::set_Hsystem(const BlockedMatrix<double>&);
