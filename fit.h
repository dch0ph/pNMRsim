#ifndef NMRSIM_FIT_H_
#define NMRSIM_FIT_H_

/*! \file
 \brief  Components used in optimisation
*/

#include "NMRsim.h"

//! default maximum fitting iterations
#define NMRSIM_DEFAULT_ITERATIONS 100
//! default stop criterion (fitting)
#define NMRSIM_DEFAULT_FIT_TOLERANCE 1e-4
//! default stop criterion (Minuit)
#define NMRSIM_DEFAULT_OPTIM_TOLERANCE 0.1
#define NMRSIM_DEFAULT_OPTIM_METHOD OPT_GRADIENT
//! Minuit strategy (1 is default)
#define NMRSIM_OPTIM_STRATEGY 1

//! base for optimisation commands
struct FitCommand {
  virtual ~FitCommand() {};
  virtual void exec() const { throw InternalError("FitCommand::exec"); } //!< not pure virtual since some commands cannot be instantiated hence don't define an exec
};

typedef FitCommand* (*FitCommand_function)(); //!< pointer to function creating new ::FitCommand

//! structure 'describing' a FitCommand
struct FitCommand_t {
  FitCommand_function funcp; //!< pointer to create function
  const char* description; //!< pointer to description (help string)
  FitCommand_t(FitCommand_function funcpv =NMRSIM_NULL, const char* descv =NMRSIM_NULL)
    : funcp(funcpv), description(descv) {}
};

void check_extract(cmatrix&); //!< apply any range selection

typedef double (*OptimiseFunction)(const DataStore&); //!< typedef for external optimisation function
typedef OptimiseFunction (*OptimiseParse)();
typedef FASTMAPTYPE(OptimiseParse) Optimise_Factory_t; 
Optimise_Factory_t& get_Optimise_Factory();

typedef FASTMAPTYPE(FitCommand_t) Fit_Factory_t; //!< type of ::FitCommand registry
Fit_Factory_t& get_Fit_Factory(); //!< get registry for (stacked) fit directives
Fit_Factory_t& get_minmax_Factory(); //!< get registry for min/max directives
typedef void (*GlobalFitCommand_function)();
typedef FASTMAPTYPE(GlobalFitCommand_function) GlobalFit_Factory_t;
GlobalFit_Factory_t& get_setup_fit_Factory(); //!< registry for non-interactive directives

extern ThreadWarning<> impossiblesave_warning; //!< save of data that is not available
extern ContextWarning<> callback_override_warning; //!< redefining calculation callback
extern ContextWarning<> fitdirective_ignored_warning; //!< ignored fitting directive
extern ThreadWarning<> fit_allfixed_warning; //!< all parameters fixed
extern ContextWarning<> fit_noparameters_warning; //!< no parameters to fix
extern ThreadWarning<> fit_directiveafterrun_warning; //!< optimisation directive after last run ignored
extern ContextWarning<> fit_expafterprevious_warning; //!< fit exp used after previous data load
extern ContextWarning<> fit_nofilenames_warning; //!< no data loaded in exp
extern ContextWarning<> fit_overwritingdomain_warning; //!< over-writting time/frequency domain specification
extern ContextWarning<> fit_domainunclear_warning; //!< time vs. frequency domain unclear
extern ContextWarning<> fit_ignoringnoiselevel_warning; //!< ignoring previous noiselevel setting
extern ContextWarning<> fit_convergefailed_warning; //!< optimisation failed to converge
extern ThreadWarning<> fit_sw_warning; //!< sw mismatch? warning
extern ThreadWarning<> fit_constraintsactive_warning; //!< warning about error bounds when constraints present
extern ThreadWarning<> fit_nonoiselevel_warning; //!< no noise level specified
extern ThreadWarning<> fit_noscale_warning; //!< no intensity scaling

extern LIST<Parameter> fitting_parameters; //!< fundamental list of fitting parameters

#endif
