#include "Parser.h"
#include "NMRsim_common.h"

bool update_interactions=true;
double proton_freq=0.0;
context_t evaluation_state=CONTEXT_MAINLOOP; //!< unless over-ridden, in main loop for simple calculations
bool proton_freq_isconstant() { return true; }
bool isallowed_nonconst_var(SystemVariableBase*) { return false; }
size_t parse_nucleusname(const char*) { throw Failed("parsing of nucleus names is not valid here"); }
spin_system* sysp=NMRSIM_NULL;

// flagsmap_type verbose_flags;

// #ifndef NDEBUG
// int def_verbose=VER_GEN | VER_PARSE;
// #else
// int def_verbose=VER_GEN;
// #endif

// void parse_verbose()
// {
//   static bool doneaddmap=false;
//   if (!doneaddmap) {
//     verbose_flags["general"]=VER_GEN;
//     //verbose_flags["powder"]=VER_POWDER;
//     //verbose_flags["optimise"]=VER_OPTIM;
//     verbose_flags["parse"]=VER_PARSE;
//     //verbose_flags["profile"]=VER_PROFILE;
//     doneaddmap=true;
//   }
//   verbose=parse_flags(verbose_flags);
//   if (!verbose)
//     verbose=def_verbose;
// }

double get_nmrfreq(size_t) { return 0.0; }
double getgrat(size_t) { return 0.0; }
