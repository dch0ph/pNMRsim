#ifndef NMRsim_spinsys_h_
#define NMRsim_spinsys_h_

/*! \file
  \brief   Header file defining core entities for spinsys block

*/

#include "NMRsim_common.h"

extern HamiltonianStore<space_T>* interactions_MFp; //!< NMR interactions store (NMRSIM_NULL if unset)

void refresh_interactions(); //!< ensure that tensor information is up-to-date in Hamiltonian store
void dump_interactions(std::ostream& =std::cout); //!< display spin system information

void parse_proton_frequency();
command_Factory_t& get_spinsys_Factory(); //!< return command factory for spinsys
command_Factory_t& initialise_spinsys_Factory(); //!< create base functionality for spinsys
enum { USER_COUPLING=0, USER_SHIFT };
extern LIST<interaction_t> weakints; //!< list of weakly coupled interactions

bool isclassicQ(); //!< \a true if "classic" 2nd order Q treatment is being used

#endif
