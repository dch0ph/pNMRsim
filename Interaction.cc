#include "Parser.h"
#include "Interaction.h"

typedef LIST<Interaction*> InteractionStack;
LIST<interaction_t> weakints;
HamiltonianStore<space_T>* interactions_MFp;

void parse_cells();
//! parse user interaction
/** qualifier distinguishes coupling vs. shift */
void parse_user(int);
void parse_nuclei();
void parse_usernucleus();

typedef std::pair<double,bool> asym_spec;

HamiltonianStore<space_T>& get_Hamiltonian()
{
  if (interactions_MFp)
    return *interactions_MFp;
  error_abort("Spin system hasn't been created");
}

const char tensorerror[]="tensor expression should evaluate to list of six items [<iso> <aniso> <asym> <alpha> <beta> <gamma>]";

namespace {
  InteractionStack inter_stack; //!< stack of interactions
  LIST<interaction_info> interaction_stack;
  inline interaction_info& get_interaction_info(interaction_t id) { return interaction_stack(size_t(id)); } //!< slightly dodgy


  size_t nspins_cell=0;
  int cells=1;

  bool parser_isarray()
  {
    char* ptr=get_curline();
    return (ptr && (*ptr=='{'));
  }
}

void parse_truncate();
void parse_tensorordering();
void parse_interaction(int);
Euler parse_orient(Variable* =NMRSIM_NULL); 
option optgeneralisedQ("generalisedQ","",option::AUTO,option::NOTUSED);
option optclassicQ("classicQ","",option::AUTO,option::NOTUSED);

bool isclassicQ() { 
  static const bool usedboth=(optgeneralisedQ.isenabled() && optclassicQ.isenabled());
  if (usedboth)
    error_abort("Cannot simultaneously enable classicQ and generalisedQ options");  
  return (optgeneralisedQ.isnotstate(option::ON) || optclassicQ.isstate(option::ON));
}

//< can no longer initialise as static in case properties table changed
double get_gamma1H()
{
  static double gamma1H=0.0;
  if (!gamma1H)
    gamma1H=gamma(nuclei_spec(H_NUCLEUS));
  return gamma1H;
}

const double rad_to_deg=180.0/M_PI;

double get_proton_freq()
{
  if (!proton_freq)
    error_abort("proton_frequency unset");
  if (!v_proton_freq.isconstant())
    error_abort("proton_frequency must be fixed");
  return proton_freq;
}

double getgrat(size_t qualifier)
{
  double grat=curgrat;
  if (qualifier) {
    if (sysp==NMRSIM_NULL)
      error_abort("can't qualify shift with nucleus type - no spin system set");
    if ((qualifier<1) || (qualifier>sysp->nspins()))
      error_abort("nucleus index out of range");  
    const size_t ind=qualifier-1;
    grat=(*sysp)(ind).gamma()/get_gamma1H();
    if (grat==0.0)
      throw InternalError("getgrat");
  }
  else {
    if (curgrat==0.0)
      error_abort("no current (or default) nucleus");
  }
  return grat;
}

std::ostream& operator<< (std::ostream& ostr, const Interaction& a)
{
  a.print(ostr);
  return ostr;
}

void dump_interactions(std::ostream& ostr)
{
  ostr << inter_stack << std::endl;
  if (verbose_level>1) {
    ostr << get_Hamiltonian();
    ostr << "Weak/truncated interactions: ";
    if (weakints.size()) {
      for (size_t i=0;i<weakints.size();i++)
	ostr << interaction_name(weakints(i)) << ' ';
      ostr << '\n';
    }
    else
      ostr << "none\n";
  }
}

double get_nmrfreq(size_t qualifier)
{
  return get_proton_freq()*getgrat(qualifier);
}

const interaction_info& interaction_create(interaction_t id, const char* name,bool iscoupling)
{
  command_Factory_t& spinsys_Factory(get_spinsys_Factory());
  interaction_stack.push_back(interaction_info(name,id,iscoupling));  
  const par_t comdesc(&parse_interaction,interaction_stack.size()-1,true);
  if (!(spinsys_Factory.insert(command_Factory_t::value_type(name,comdesc)).second))
    throw InternalError("can't create interaction twice!");    
  return interaction_stack.back();
}

command_Factory_t& get_spinsys_Factory()
{
  static command_Factory_t spinsys_Factory;
  return spinsys_Factory;
}

/* create core spinsys functionality */

command_Factory_t& initialise_spinsys_Factory()
{
  static command_Factory_t& spinsys_Factory(get_spinsys_Factory());
  static bool doneinit=false;
  if (doneinit)
    std::cerr << "Warning: initialise_spinsys_Factory called twice!\n";
  doneinit=true;
  spinsys_Factory["nuclei"]=&parse_nuclei;
  spinsys_Factory["cells"]=&parse_cells;
  spinsys_Factory["usercoupling"]=par_t(&parse_user,USER_COUPLING,true);
  spinsys_Factory["usernucleus"]=par_t(&parse_usernucleus,true);
  spinsys_Factory["usershift"]=par_t(&parse_user,USER_SHIFT,true);
  spinsys_Factory["proton_frequency"]=&parse_proton_frequency;
  spinsys_Factory["precision"]=&parse_precision;
  return spinsys_Factory;
}

ContextWarning<> protonfrequency_warning("proton Larmor frequency is less than 1 MHz.  High field approximation may become questionable, or wrong units (Hz) in proton_frequency?",&NMRsim_once_warning);
ThreadWarning<> exact_warning("'Exact' treatment of quadrupoles is not fully evaluated and will silently fail if used inappropriately (e.g. NQR limit)",&NMRsim_once_warning,BaseWarning::Inherit,std::cout);

void parse_proton_frequency()
{
//   if (!interactions_MFp)
//     error_abort("proton_frequency can't be set before spin system is created");
  
  parse_system_variable(v_proton_freq);
  if (proton_freq<=0)
    error_abort("Proton frequency cannot be <= 0!");
  if (proton_freq<1e6)
    protonfrequency_warning.raise();
}

// class FunctionTensor20 : public ExpressionSimpleFunction {
// public:
//   explicit FunctionTensor20(size_t nargs) : ExpressionSimpleFunction("tensor20",nargs) {}
//   explicit FunctionTensor20(const BaseList<size_t>& childv) : ExpressionSimpleFunction("tensor20",childv) {}
//   ExpressionBase* clone() const { return new FunctionTensor20(*this); }
//   static ExpressionFunctionBase* create(size_t nargs) { return new FunctionTensor20(nargs); }
//   static ExpressionFunctionBase* create(const BaseList<size_t>& childv) { return new FunctionTensor20(childv); }
//   void operator()(LIST<double>&, const BaseList<double>&) const;
// private:
//   mutable space_T A_PAS,A_MF;
//   mutable double lastaniso,lastasym;
//   mutable Euler lasteuler;
//   mutable bool validMF;
// };

// void FunctionTensor20::operator()(LIST<double>& dest, const BaseList<double>& in) const
// {
//   static bool donewarn=nochecks || silent;
//   if (!donewarn && spin_rate) {
//     std::cerr << "Warning: tensor20 only makes sense in context of static simulations\n";
//     donewarn=true;
//   }
//   const double aniso(in.front());
//   const double asym(in(1U));
//   if (!!A_PAS || (aniso!=lastaniso) || (asym!=lastasym)) {
//     A_PAS=spatial_tensor(aniso,asym);
//     lastaniso=aniso;
//     lastasym=asym;
//     validMF=false;
//   }

//   const space_T* usep=&A_PAS;
//   if (in.size()==5) {
//     const Euler PtoM(in(2U),in(3U),in(4U));
//     if (!validMF || (lasteuler!=PtoM)) {
//       A_MF=rotate(A_PAS,PtoM);
//       lasteuler=PtoM;
//       validMF=true;
//     }
//     usep=&A_MF;
//   }
//   dest.push_back(real(rotate(*usep,2,0,global_powder)));
// }

struct Proxy_ {
  Proxy_() {
    interaction_create(I_DIPOLE,"dipole",true);
    interaction_create(I_J,"jcoupling",true);
    interaction_create(I_CS,"shift",false);
    interaction_create(I_QUAD,"quadrupole",false);

    command_Factory_t& spinsys_Factory(get_spinsys_Factory());
    spinsys_Factory["truncate"]=&parse_truncate;
    spinsys_Factory["tensorordering"]=par_t(&parse_tensorordering,true);

    optional_map_t& optional_map(get_optional_map());
    optional_map["generalisedQ"]=&optgeneralisedQ;
    optional_map["classicQ"]=&optclassicQ;
    
//     Function_Factory_t& funcfac(get_Function_Factory());
//     funcfac[function_spec("tensor20",2U)]=function_def_t(new FunctionTensor20(2U));
//     funcfac[function_spec("tensor20",5U)]=function_def_t(new FunctionTensor20(5U));
  }
};

const Proxy_ proxy_;

//bool ishomogeneous=false;
const bool weakcoupling=false;

static const double deg_to_rad=M_PI/180.0;

Interaction::Interaction(const interaction_info& info_,size_t ni_,size_t nj_,double iso_,double aniso_,double asymspecv, bool isxyv, const Euler& PAS_, ordering_convention_t ordering_, subsid_t subsid, double scale_)
  : info(info_),
    ni(ni_), nj(nj_),
    iso(iso_), aniso(aniso_*scale_), asymspec(asymspecv),
    PAS(PAS_),
    ordering(ordering_),
    aniso_scale(scale_),
    isdirty_(true),
    isxy_(isxyv),
    allowtensor_(subsid)
{
  if (isxy_)
    asymspec*=scale_;
}

Interaction::Interaction(const interaction_info& info_,size_t ni_,size_t nj_,double iso_, subsid_t subsid)
  : info(info_),
    ni(ni_), nj(nj_),
    iso(iso_), aniso(0.0), asymspec(0.0),
    ordering(convention_Haeberlen), //!< not used, but complete with sensible value
    aniso_scale(1.0),
   isdirty_(true), isxy_(false),
    allowtensor_(subsid)
{}

bool Interaction::isset() const
{
  return info.iscoupling ? get_Hamiltonian().isset(info.id,ni,nj) : get_Hamiltonian().isset(info.id,ni);
}

void Interaction::ensurevalid()
{
  if (!isdirty_)
    return;
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    std::cout << "Updating tensor information for ";
    print(std::cout);
  }
  update_tensor();
}
  
void Interaction::update_tensor()
{
  //const bool verb = (verbose_level>1) && (verbose & VER_GEN);  
  if (aniso) {
    double eta=asymspec;
    if (isxy_)
      eta/=aniso;
    space_T A_PAS(spatial_tensor(iso,aniso,eta,ordering));
//     if (verb)
//       std::cout << "Tensor in PAS\n" << A_PAS << '\n';
    if (iso==0.0)
      A_PAS.clear(0);
    const space_T A_MF(rotate(A_PAS,PAS));
//     if (verb)
//       std::cout << "Tensor in MF\n" << A_MF << '\n';
    update_tensor(A_MF);
  }
  else {
    space_T A_PAS(0);
    A_PAS(0,0)=iso;
    update_tensor(A_PAS);
  }
  isdirty_=false;
}

void Interaction::update_tensor(const space_T& A)
{
  if (info.iscoupling)
    interactions_MFp->set_coupling(info.id,ni,nj,A);
  else {
    if (info.id==I_QUAD) 
      interactions_MFp->set_quadrupole(ni,A,nj);
    else 
      interactions_MFp->set_shift(ni,A,info.id);
  }
}

void Interaction::set(const BaseList<double>& vals, subsid_t subsid)
{
  switch (vals.size()) {
  case 1:
    set(vals.front(),subsid);
    return;
  case 6:
    if (subsid!=allowtensor_)
      error_abort("can only define full tensor on first item (isotropic value or anisotropy as appropriate for interaction)");
    iso=vals.front();
    aniso=vals(1U)*aniso_scale;
    asymspec=vals(2U);
    if (isxy_)
      asymspec*=aniso_scale;
    PAS.alpha=vals(3U)*deg_to_rad;
    PAS.beta=vals(4U)*deg_to_rad;
    PAS.gamma=vals(5U)*deg_to_rad;
    isdirty_=true;
    update_interactions=true;
    break;
  default:
    error_abort(tensorerror);
  }
}

void Interaction::set(double fval,subsid_t subsid)
{
  switch (subsid) {
  case S_ISO:
    iso=fval;
    break;
  case S_ANISO:
    aniso=fval*aniso_scale;
    break;
//   case S_XY:
//     if (aniso==0)
//       error_abort("detected zero anisotropy while setting asymmetry");
//     eta=fval*aniso_scale/aniso;
//     break;
  case S_ASYM:
    if (isxy_)
      fval*=aniso_scale;
    asymspec=fval;
    break;
  case S_ALPHA:
    PAS.alpha=fval*deg_to_rad;
    break;
  case S_BETA:
    PAS.beta=fval*deg_to_rad;
    break;
  case S_GAMMA:
    PAS.gamma=fval*deg_to_rad;
    break;
  default:
    throw InternalError("unhandled subsid type");
  }
  // update_tensor();
  isdirty_=true;
  update_interactions=true;
}

void Interaction::print(std::ostream& ostr) const
{   
  ostr << interaction_name(info.id) << ' ' << (ni+1);
  if (info.iscoupling)
    ostr << ' ' << (nj+1);
  else {
    if (info.id==I_QUAD)
      ostr << ' ' << nj;  //nj is order
  }
  const bool hasaniso=(info.id==I_QUAD) || (info.id==I_DIPOLE);
  if (!hasaniso)
    ostr << ' ' << iso;

  if (aniso || hasaniso) {
    ostr << ' ' << aniso/aniso_scale << ' ';
    if (asymspec || (info.id!=I_DIPOLE)) {
      if (isxy_) 
	ostr << asymspec/aniso_scale << " [xy] ";
      else
	ostr << asymspec << ' ';
    }
    ostr << PAS;
  }
  switch (ordering) {
  case convention_Haeberlen:
    ostr << " [Haeberlen]";
    break;
  case convention_NQR:
    ostr << " [NQR]";
    break;
  default:
    ostr << " [???]"; //!< no need to fail as reporting only
    break;
  }
  if (verbose_level>1)
    ostr << (isdirty_ ? " [dirty]" : " [clean]");
}

// void Interaction::print(std::ostream& ostr,subsid_t subsid) const
// {   
//   ostr << interaction_name(info.id) << ' ' << (ni+1);
//   if (info.iscoupling)
//     ostr << ' ' << (nj+1);
//   else {
//     if ((info.id==I_QUAD) && (subsid==S_NONE)) 
//       ostr << ' ' << nj;  //nj is order
//   }
//   if (subsid) {
//     switch (subsid) {
//     case S_ISO: ostr << " isotropic=" << iso; break;
//     case S_ANISO: ostr << " anisotropy=" << aniso; break;
//       //    case S_XY: ostr << " xx-yy=" << aniso*eta; break;
//     case S_ASYM:
//       ostr << (isxy_ ? " xx-yy=" : " asymmetry=") << asymspec; break;
//     case S_ALPHA: ostr << " alpha=" << (PAS.alpha*rad_to_deg); break;
//     case S_BETA: ostr << " beta=" << (PAS.beta*rad_to_deg); break;
//     case S_GAMMA: ostr << " gamma=" << (PAS.gamma*rad_to_deg); break;
//     default: throw InternalError("Unknown interaction subsid type");
//     }
//   }
//   else {
//     switch (info.id) {
//     case I_QUAD: case I_DIPOLE:
//       break;
//     default:
//       ostr << ' ' << iso;
//     }
//     if (aniso) {
//       ostr << ' ' << aniso << ' ';
//       if (asymspec || (info.id!=I_DIPOLE))
// 	ostr << asymspec << (isxy_ ? " [xy] " : " ");
//       ostr << PAS;
//     }
//     if (verbose_level>1)
//       ostr << (isdirty_ ? " [dirty]" : " [clean]");
//     ostr << '\n';
//   }
// }

void parse_user(int which)
{
  const char* name(parse_string(F_REPLACEDOLLAR));
  if (*name=='-')
    error_abort("Interaction name cannot start with -");
  const interaction_t id=interaction_register(name);
  interaction_create(id,name,which==USER_COUPLING);
}

static void parse_spinargs(int& ni, int& nj, size_t nargs, size_t nspins_cell, size_t nspins_total)
{
  ni=parse_index(nspins_cell);
  if (ni<0)
    error_abort();
  nj=0;
  if (nargs>1) {
    nj=parse_index(nspins_cell,nspins_total);
    if (nj<0)
      error_abort();
    if (ni==nj) {
      parser_printcontext() << "Can't have interaction between same spin (" << (ni+1) << ")\n";
      error_abort();
    }
  }    
}

ContextWarning<> nonzeroasymmetry_warning("non-zero asymmetry with zero anisotropy",&NMRsim_repeat_warning);

static asym_spec parse_asym(double aniso, double ratX, Variable* cvarp =NMRSIM_NULL)
{
  if (cvarp)
    cvarp->subsid=S_ASYM;
  static flagsmap_type flags;
  if (flags.empty())
    flags["xy"]=1;
  const int isxy=parse_flags(flags);
  double fval=(ratX && isxy)
    ? parse_shift(cvarp,ratX)
    : parse_double(cvarp);
  if (!isxy && (fval<0.0 || fval>1.0)) {
    std::cerr << "Invalid asymmetry: " << fval << std::endl;
    error_abort();
  }
  if (fval && (aniso==0.0)) {
    nonzeroasymmetry_warning.raise();
    fval=0.0;
  }
  return std::pair<double,bool>(fval,isxy);
}

interaction_info& parse_interaction_name_syntax(const char* syntax, size_t narg)
{
  const char* intname(parse_string_syntax(syntax,narg,F_REPLACEDOLLAR));
  command_Factory_t& spinsys_Factory(get_spinsys_Factory());
  const command_Factory_t::iterator curp=spinsys_Factory.find(intname);
  if ((curp==spinsys_Factory.end()) || (curp->second.parsefunc_qual!=&parse_interaction)) {
    parser_printcontext() << "invalid interaction name: " << intname << '\n';
    error_abort();
  }
  return get_interaction_info(curp->second.ident);
}

void parse_truncate()
{
  static const char syntaxstring[]="<list of interaction names>";
  do {
    const interaction_info& info(parse_interaction_name_syntax(syntaxstring,1));
    if (!info.iscoupling)
      parser_printcontext() << "interaction is not a coupling: " << info.name << " (ignored)\n";
    else
      weakints.push_back(info.id);
  } while (are_left());
}

void parse_tensorordering()
{
  static flagsmap_type orderflags;
  enum { ORD_HAEBERLEN=1, ORD_NQR=2 };
  if (orderflags.empty()) {
    orderflags["Haeberlen"]=ORD_HAEBERLEN;
    orderflags["NQR"]=ORD_NQR;
  }
  LIST<interaction_info*> intlist;
  static const char syntaxstring[]="<list of interaction names>#<interaction ordering flag>";
  while (!(parser_isflag())) {
    if (are_left())
      intlist.push_back(&(parse_interaction_name_syntax(syntaxstring,1)));  
    else {
      // output ordering info
      for (size_t i=0;i<intlist.size();i++) {
	const interaction_info& info(*(intlist(i)));
	std::cout << info.name << ' ';
	switch (info.ordering) {
	case convention_Haeberlen:
	  std::cout << "-Haeberlen\n";
	  break;
	case convention_NQR:
	  std::cout << "-NQR\n";
	  break;
	default:
	  throw InternalError("Unknown tensor ordering convention");
	}
      }
      return;
    }
  }      
  if (intlist.empty())
    error_abort("No interaction names given");
  const int ord=parse_flags(orderflags,0,true);
  ordering_convention_t conv=convention_Haeberlen;
  switch (ord) {
  case ORD_HAEBERLEN:
    break;
  case ORD_NQR:
    conv=convention_NQR;
    break;
  default:
    throw InternalError("Unknown tensor ordering convention (2)");
  }
  for (size_t i=intlist.size();i--;) {
    interaction_info& info(*(intlist(i)));
    info.ordering=conv;
  }
}

ContextWarning<> dipolesign_warning("dipolar coupling has incorrect sign (for unaveraged coupling)",&NMRsim_once_warning);

bool parse_first(double& ciso, double &caniso, asym_spec& casym, Euler& PAS, subsid_t& allowtensor, Variable* cvar, subsid_t subsid, double ratX =0.0)
{
  static flagsmap_type pflags;
  if (pflags.empty())
    pflags["xy"]=1;
  cvar->subsid=allowtensor=subsid;
  const int flags=((subsid==S_ANISO) ? F_DENYZERO : 0) | F_ALLOWLIST;
  double v= ratX ? parse_shift(cvar,ratX,flags) : parse_double(cvar,flags);
  VariableBase* valuep=cvar->valuep;
  if (valuep && valuep->isconst()) {
    const BaseList<double>& vals(valuep->value());
    switch (vals.size()) {
    case 1:
      v=vals.front();
      break;
    case 6: {
      const bool isxy=parse_flags(pflags);
      ciso=vals.front();
      caniso=vals(1U);
      casym=asym_spec(vals(2U),isxy);
      PAS.alpha=vals(3U);
      PAS.beta=vals(4U);
      PAS.gamma=vals(5U);	
    }
      return true;
    default:
      error_abort(tensorerror);
    }
  }
  switch (subsid) {
  case S_ISO:
    ciso=v;
    break;
  case S_ANISO:
    caniso=v;
    break;
  default:
    throw InternalError("parse_first");
  }
  if (are_left()) {
    if (parser_isnormal() || (count_left()>1)) {
      allowtensor=S_NONE; //!< disallow tensor specification
      return false;
    }
    const bool isxy=parse_flags(pflags);
    casym.second=isxy;
  }      
  return true;
}

// double Interaction::gamma() const
// {
//   return (*sysp)(ni).gamma();
// }

static void parse(const interaction_info& info)
{
  Mark markobj;
  Variable cvar;
  Interaction* newint=NMRSIM_NULL;
  int nni,nj,order;
  double caniso=0.0;
  double ciso=0.0;
  subsid_t allowtensor=S_NONE;
  asym_spec casym(0.0,false);
  Euler PAS(0,0,0);
  const size_t nspins_cell(get_Hamiltonian().nspins_cell()); //!< use get_Hamiltonian once to ensure system exists
  const size_t nspins(interactions_MFp->nspins());
  const ordering_convention_t conv=info.ordering;

  switch (info.id) {
  case I_DIPOLE: {
    parse_spinargs(nni,nj,2,nspins_cell,nspins);
    if (!parse_first(ciso,caniso,casym,PAS,allowtensor,&cvar,S_ANISO)) {
      size_t argsleft=count_left();
      if (parser_isflag())
	argsleft--; //!< discount -xy flag
      if ((argsleft==1) || (argsleft==4)) {
	casym=parse_asym(caniso,0.0,&cvar); //relevant to dynamically averaged systems
	argsleft--;
      } //don't check sign in dynamically averaged systems (P2 factor may be +ve or -ve)
      else {
	static bool checksign=dipolesign_warning.enabled(); //!< assume not going to change
	if (checksign && caniso*(*sysp)(nni).gamma()*(*sysp)(nj % nspins_cell).gamma()>0.0) {
	  char buf[256];
	  snprintf(buf,sizeof(buf),": dipole %i %i",nni+1,nj+1);
	  dipolesign_warning.raise(buf);
	}    
      }
      switch (argsleft) {
      case 0:
	break;      
      case 3:
	PAS=parse_orient(&cvar);
	break;
      default:
	error_abort("dipole must specify <anisotropy> [<asymmetry>] [<PAS>]");
      }
    }
    newint= new Interaction(info,nni,nj,0.0,caniso,casym.first,casym.second,PAS,conv,allowtensor);
  }
    break;
  case I_QUAD: {
    parse_spinargs(nni,nj,1,nspins_cell,nspins);
    const spin& cspin((*sysp)(nni));
    const size_t deg=cspin.deg();
    if (deg<3)
      error_abort("nucleus is not quadrupolar!");
    
    const double scale=1.0/(2*(deg-2)*(deg-1)); // 4I(2I-1) factor
    order=parse_int();
    switch (order) {
    case 1:
      break;
    case 0: case 2:
      if (optclassicQ.isenabled() && (order==0))
	error_abort("Can't use -enable:classicQ with order '0'");
      if (!proton_freq)
	error_abort("must define field (via proton_frequency) before using 2nd order / exact quadrupoles");
      if (!isclassicQ() || (order==0)) {
	if (cspin.isrestricted())
	  error_abort("cannot combine generalised quadrupole treatment with restriction to central transition");
	//	ishomogeneous=true;
	optgeneralisedQ.set(option::ON); //!< implicitly turn on so that isclassicQ then returns right answer
	optgeneralisedQ.setusage(true);
	exact_warning.raise();
      }
      else {
	if ((nspins>1) && !(optclassicQ.isenabled()))
	  error_abort("classic second-order quadrupole treatment only valid for single (uncoupled) nuclear spin; either enable experimental multiple spin treatment with -enable:generalisedQ or use -enable:classicQ to force classic (SIMPSON) behaviour.");
	optclassicQ.setusage(true);
      }
      break;
    default:
      error_abort("quadrupole order must be 0, 1 or 2");
    }
    if (!parse_first(ciso,caniso,casym,PAS,allowtensor,&cvar,S_ANISO))
      casym=parse_asym(caniso,0.0,&cvar);
    newint = new Interaction(info,nni,order,0.0,caniso,casym.first,casym.second,parse_orient(&cvar),conv,allowtensor,scale);
  }
    break;
  default:
    if (info.iscoupling) { //I_J and other general couplings
      parse_spinargs(nni,nj,2,nspins_cell,nspins);
      if (!parse_first(ciso,caniso,casym,PAS,allowtensor,&cvar,S_ISO)) {
	cvar.subsid=S_ANISO;
	caniso=parse_double(&cvar,F_DENYZERO);
	casym=parse_asym(caniso,0.0,&cvar);
	if (are_left())
	  PAS=parse_orient(&cvar);
      }
      newint = (caniso) ? 
	new Interaction(info,nni,nj,ciso,caniso,casym.first,casym.second,PAS,conv,allowtensor) :
	new Interaction(info,nni,nj,ciso,allowtensor);
    }
    else { //I_CS and other shifts
      parse_spinargs(nni,nj,1,nspins_cell,nspins);
      const double ratX=(*sysp)(nni).gamma()/get_gamma1H();
      if (ratX==0.0)
	throw InternalError("ratX unset");
      if (!parse_first(ciso,caniso,casym,PAS,allowtensor,&cvar,S_ISO,ratX)) {
	cvar.subsid=S_ANISO;
	caniso=parse_shift(&cvar,ratX,F_DENYZERO);
	casym=parse_asym(caniso,ratX,&cvar);
	if (are_left())
	  PAS=parse_orient(&cvar);
      }
      newint = (caniso) ?
	new Interaction(info,nni,0,ciso,caniso,casym.first,casym.second,PAS,conv,allowtensor) :
	new Interaction(info,nni,0,ciso,allowtensor);
    }
  }
  if (newint->isset()) { 
    parser_printcontext() << "Attempt to redefine ";
    newint->print(std::cerr);
    error_abort();
  }
  newint->update_tensor(); //!< need to create so that isset can check
  inter_stack.push_back(newint);
  markobj.flush(newint);
}

void refresh_interactions()
{
  for (size_t i=inter_stack.size();i--;)
    inter_stack(i)->ensurevalid();
}

void parse_interaction(int which)
{
  return parse(interaction_stack(size_t(which)));
}

bool iscommonPAS=true;

Euler parse_orient_raw(Variable* cvarp)
{
  static const subsid_t mapping[3]={ S_ALPHA, S_BETA, S_GAMMA };

  double angles[3];
  for (size_t j=0;j<3;j++) {
    if (cvarp)
      cvarp->subsid=mapping[j];
    angles[j]=parse_double(cvarp)*M_PI/180.0;
  }
  return Euler(angles[0],angles[1],angles[2]);
}

Euler parse_orient(Variable* cvarp)
{
  static bool havePAS=false;
  static Euler Common_PAS;

  if (!are_left())
    return Euler(0,0,0);

  const Euler retval(parse_orient_raw(cvarp));
  if (havePAS) {  // check for unique PAS
    if (retval!=Common_PAS)
      iscommonPAS=false;
  }
  else {
    havePAS=true;
    Common_PAS=retval;
  }
  return retval;
}

// void Interaction::convert_anisotropy(Variable& cvar) const
// {
//   if (cvar.subsid!=S_ETA)
//     throw Failed("Interaction::covert_aniso: Variable is not ETA");
//   if (aniso==0.0)
//     throw InternalError("Can't optimise eta if no anisotropy");
//   cvar.subsid=S_XY;
//   BaseList<double> rawlist(cvar.variable().value());
//   rawlist*=aniso;
// }

void Interaction::printvariablename(std::ostream& dest, subsid_t subsid) const
{
  dest << info.name << '_' << (ni+1);
  if (info.iscoupling)
    dest << '_' << (nj+1);
  dest << '_';
  switch (subsid) {
  case S_ISO: dest << "iso"; break;
  case S_ANISO: dest << "aniso"; break;
  case S_ASYM: dest << (isxy_ ? "xy" : "eta"); break;
  case S_ALPHA: dest << "alpha"; break;
  case S_BETA: dest << "beta"; break;
  case S_GAMMA: dest << "gamma"; break;
  default: 
    throw InternalError("Unhandled interaction subsid type");
  }
}

void parse_usernucleus()
{
  static flagsmap_type flags;
  if (flags.empty())
    flags["gamma"]=1;

  const char* name(parse_string(F_REPLACEDOLLAR));
  char* iname=get_curline();
  size_t deg;
  double gamma;
  try {
    const spin s(iname);
    deg=s.deg();
    gamma=s.gamma();
  }
  catch (MatrixException&) {
    const double I=parse_double();
    if ((I<=0.0) || ((2*I-floor(2*I))>NMRSIM_ROUNDTOL))
      error_abort("Invalid spin quantum number");
    deg=static_cast<size_t>(2*I+1.5);
    gamma=parse_double();
    if (gamma==0.0)
      error_abort("NMR frequency / gamma cannot be zero!");
    if (!parse_flags(flags)) 
      gamma*=get_gamma1H()/get_proton_freq();
  }
  define_nucleus(name,deg,gamma);
}

ContextWarning<> spinhalfrestrict_warning("Central transition restriction (:c) applied to spin-1/2 nucleus is being ignored",&NMRsim_once_warning);

void create_spin_system()
{
  LIST<spin> spinlist;

  char* cptr;
  while ((cptr=parse_string(F_ALLOWMISSING))) {
    bool restrict=false;
    char* tag=cptr+strlen(cptr)-2;
    if ((tag>cptr) && (strcmp(tag,":c")==0)) { //look for :c tag
      *tag='\0';
      restrict=true;
    }
    spinlist.push_back(spin(cptr,restrict));
    if (restrict && (spinlist.back().deg()==2) && !nochecks) //!< :c applied to integer spins already caught
      spinhalfrestrict_warning.raise(); 
  }
  if (spinlist.empty())
    error_abort("empty nuclei specification");

  const size_t lnspins_cell=spinlist.length();
  if (cstructp && cstructp->haspermutation() && (cstructp->permutation().size()!=lnspins_cell))
    error_abort("spins per cell doesn't match previous permutation length");

  sysp=new spin_system(lnspins_cell,spinlist.front());
  for (size_t i=1;i<lnspins_cell;i++)
    (*sysp)(i)=spinlist(i);
  if (verbose & VER_GEN)
    std::cout << "Spin system: " << (*sysp) << '\n';
//   if (cstructp) { //check against CrystalStructure
//     const Permutation& perm(cstructp->permutation());
//     if (!perm.empty()) {
//       if (perm.size()!=lnspins_cell)
// 	error_abort("Permutation in cells doesn't match number of spins in nuclei");
//       for (size_t i=lnspins_cell;i--;) {
// 	if ((*sysp)(i)!=(*sysp)(perm(i)))
// 	  error_abort("Permutation mixes nuclei types");
//       }
//     }
//   }
//   else
  if (!cstructp)
    cstructp=new CrystalStructure(); //no periodic structure
  if (sysp->ishomonuclear())
    curgrat=defgrat=gamma(nuclei_spec(sysp->homonucleus()))/get_gamma1H();
}

void parse_nuclei()
{
  create_spin_system();
  nspins_cell=sysp->nspins();
  interactions_MFp = new HamiltonianStore<space_T>(nspins_cell,cells);
  ScratchList<size_t> usednuc(MAX_NUCLEUS,false);
  for (size_t j=nspins_cell;j--;) {
    const size_t nuc=(*sysp)(j).nucleus();
    if (!usednuc(nuc)) {
      blockingnuclei.push_back(nuc);
      usednuc(nuc)=true;
    }
  }
}

void parse_cells()
{
  if (sysp)
    error_abort("can't set cells after spin system has been created");

  LIST<size_t> lcells;
  Permutation* permp=NMRSIM_NULL;
  while (are_left()) {
    if (parser_isarray()) {
      const LIST<size_t> tmp(parse_unsignedintarray(1)); //!< subtract 1 from indices
      permp=new Permutation(tmp);
      break;
    }
    const size_t n=parse_int();
    if (n<1)
      error_abort("cells cannot be <1");
    lcells.push_back(n);
  }
  cstructp = permp ? new CrystalStructure(lcells,*permp) : new CrystalStructure(lcells);
  cells=cstructp->ncells();
}

int parse_index(size_t nspins_cell, size_t nspins_total)
{
  if (!cstructp)
    throw Failed("parse_index: no periodic structure defined");

  size_t inds[4];
  size_t curind=0;
  Splitter splitter(get_token(),",");
  char* cptr;
  char* tail;
  const size_t periodic_dims(cstructp->dimensions()); //only needs evaluating once;
  const size_t maxinds=nspins_total ? 1+periodic_dims : 1;
  while ((cptr=splitter.next())) {
    if (curind==maxinds) {
      parser_printcontext() << "too many indices\n";
      return -1;
    }
    const long ind=strtol(cptr,&tail,10);
    if (*tail || (ind<1)) {
      parser_printcontext() << "Failed to parse " << cptr << " as index\n";
      return -1;
    }
    inds[curind++]=ind-1;
  }
  if (curind==1) {
    const int max=nspins_total ? nspins_total : nspins_cell;
    if (inds[0]>=max) {
      if (max==1)
	parser_printcontext() << "spin index is >1, when system only contains one spin!\n";
      else
	parser_printcontext() << "spin index (" << (1+inds[0]) << ") is not in range 1 to " << max << '\n';
      return -1;
    }
    return inds[0];
  }
  if (curind!=periodic_dims+1) {
    parser_printcontext() << "number of indices does not match number of perodic dimensions\n";
    return -1;
  }
  const size_t cell=(*cstructp)(BaseList<size_t>(periodic_dims,inds));
  const size_t total=cell*nspins_cell+inds[curind-1];
  if (total>=nspins_total) {
    parser_printcontext() << "spin index out of range\n";
    return -1;
  }
  return total;
}
