// stub for NMRacq.o

#include "timer.h"
#include "NMRsim_common.h"

bool ammaster=true;
bool process2D=true;
dirty_t update_propagators=DIRTY_CLEAN;
void verify_powderaverage() {}
void ensure_sigma0detect() {}

rmatrix vararray,valerrs;
LIST<const char*> varnames;
Matrix<accumulating_timer<> >* timers_arrayp=NMRSIM_NULL;

size_t global_workers=0;
void save_parameters() {}
void update_auxiliary_vars(bool) {}
size_t get_thread_num() { return 0; }
bool transitions_log_active() { return false; }
std::ostream& parser_printthread(std::ostream& ostr) { return ostr; }
