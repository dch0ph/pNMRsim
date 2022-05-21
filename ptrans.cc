#include "Action.h"
#include "Parser.h"

using namespace libcmatrix;
using namespace std;

command_Factory_t par_Factory;
command_Factory_t spinsys_Factory;

int verbose_level=1;
bool nochecks=false;
bool need_spinsys=false;
double proton_freq=0.0;

// struct Translator {
//   virtual void init(const char*, const actionstack_t&) =0;
//   virtual void finish() =0;
//   virtual ~Translate() {}
// };

// struct TranslateInfinity {
//   TranslateInfinity(bool isIplusv) : isIplus_(isIplusv), fp_(NMRSIM_NULL {}
//   ~TranslateInfinity() {
//     if (fp_)
//       fclose(fp_);
//   }
//   void init(const char*, const actionstack_t&);
//   void finish();
//   void checkopen() const {
//     if (!fp_)
//       throw InternalError("checkopen");
//   }
//   bool isIplus_;
//   FILE* fp_;
// };

// void TranslateInfinity::init(const char* fname, const actionstack_t&)
// {
//   if (fp_)
//     throw InternalError("init");
// }

// void TranslateInfinity::finish()
// {
//   checkopen();
//   fclose(fp_);
//   fp_=NMRSIM_NULL
// }

void output_line(const char* buffer)
{
  std::cout << buffer << '\n';
}

//void EventID::buildevents() {}
//void EventID::set(double, size_t, subsid_t) {}
//void EventID::duration(double dur) { nomdur_=dur; }

bool proton_freq_isconstant()
{ throw InternalError("proton_freq_isconstant"); }

Variable* create_fitting_variable(const VarVariable&)
{
  error_abort("Can't use fitting variables");
  return NMRSIM_NULL // dummy to avoid warning
}

ActionCommand* ActionFilter::create()
{
  std::cerr << "filter cannot be translated (ignored)\n";
  return NMRSIM_NULL
}

ActionCommand* ActionPutMatrix::create()
{
  return NMRSIM_NULL //!< putmatrix is silently ignored
}

int auto_vars=0;
CrystalStructure* cstructp=NMRSIM_NULL //shouldn't be referred to

double get_nmrfreq()
{
  throw InternalError("get_nmrfreq"); //this should never be required
}

const basespin_system* get_spin_system() { return NMRSIM_NULL }

void post_channels() {}

bool isFop(const setableoperator_spec& mspec, char op)
{
  if (!mspec.isconstant())
    return false;
  const productoperator_spec& pspec(mspec.front());
  if (!pspec.issimple())
    return false;
  const operator_spec spec(pspec.front().front());
  return spec.issumoperator() && (spec.op==op);
}

int main(int argc_, char **argv_)
{
  int count=1;

  if (count==argc_) {
    cerr << "Syntax: ptrans <.in file>\n";
    return 1;
  }
  
  char* fname=argv_[count];
  char fnamebase[256]="";
  if (strcmp(fname,"-")!=0) {    
    strncpy(fnamebase,fname,sizeof(fnamebase));
    stripleaf(fnamebase,".in");
    systemvarmap["name"]=new SystemVariable<std::string>("name",fnamebase);
  }

  size_t argc=argc_-count-1;
  char** argv=argv_+count+1;

  declare_builtin_block("spinsys");
  declare_builtin_block("par");
  declare_builtin_block("pulseq");

  try {

  parser_init(fname,argc,argv);

  //spinsys_Factory["nuclei"]=&parse_nuclei;
  //spinsys_Factory["cells"]=&parse_cells;
  spinsys_Factory["channels"]=&parse_channels;
  spinsys_Factory["variable"]=par_t(&parse_variable,true);
  spinsys_Factory["verbose"]=&parse_verbose;
  spinsys_Factory["time_resolution"]=&parse_time_resolution;
  //spinsys_Factory["transients"]=&parse_transients;
  //spinsys_Factory["usercoupling"]=par_t(&parse_user,USER_COUPLING,true);
  //spinsys_Factory["usershift"]=par_t(&parse_user,USER_SHIFT,true);
  //spinsys_Factory["precision"]=&parse_precision;
  //spinsys_Factory["proton_frequency"]=&parse_proton_frequency;
  
  (void)read_block("spinsys",spinsys_Factory,true,true);
  if (nchannels==0)
    error_abort("No RF channels - nothing to do!");

  make_common_par_variables();

  par_Factory["start_operator"]=&parse_start_operator;
  par_Factory["detect_operator"]=&parse_detect_operator;
  //par_Factory["spin_rate"]=&parse_spin_rate;
  //par_Factory["gamma_angles"]=&parse_gamma_angles;
  //par_Factory["gamma_zero"]=&parse_gamma_zero;
  par_Factory["np"]=&parse_np;
  par_Factory["sw"]=&parse_sw;
  //par_Factory["method"]=par_t(&parse_ignored,true);
  //par_Factory["precision"]=&parse_precision;
  //par_Factory["rotor_angle"]=&parse_rotor_angle;
  //par_Factory["crystal_file"]=&parse_crystal_file;
  //par_Factory["log_file"]=par_t(&parse_log_file,true);
  //par_Factory["histogram"]=&parse_histogram;
  //par_Factory["echo"]=par_t(&parse_par_echo,1,true);
  //par_Factory["puts"]=par_t(&parse_par_echo,0,true);
  //  par_Factory["putmatrix"]=par_t(&parse_par_putmatrix,true);
  par_Factory["delay"]=par_t(&parse_delay,true);
  par_Factory["variable"]=par_t(&parse_variable,true);
  par_Factory["verbose"]=&parse_verbose;
  //par_Factory["matrix"]=par_t(&parse_matrix,true);
  //par_Factory["autoopt"]=&parse_autoopt;
  //par_Factory["fit"]=par_t(&parse_fit,true);
  //par_Factory["minimise"]=par_t(&parse_minimise,true);
  //par_Factory["maximise"]=par_t(&parse_maximise,true);
  //par_Factory["maxdt"]=&parse_maxdt;
  //par_Factory["tolerance"]=&parse_tolerance;

  (void)read_block("par",par_Factory,true,true);

  //  make_pulseq_variables();
  parser_newblock("pulseq",false);
  rebuild_sequences(); //necessary before sequences can be used
  
  char* lbuf;
  while ((lbuf=parser_getline())) { //parse pulseq block
    if (verbose & VER_PARSE)
      std::cout << "Parsing: " << lbuf << '\n';
    set_curline(lbuf);
    if (!actionstack.parse() || !parser_isfinished()) 
      return ERR_INVALID_INPUT; 
  }

  sum_dims.get(sum_ns,array_skips,sum_n0); //array_skips is ignored
  array_dims.get(array_ns,array_skips,array_n0); //extract dimensionality
  if ((array_n0!=1) || (sum_n0!=1))
    error_abort("Can't translate sequence involving arrays");

  if (sigma0_specp) {
    if (!isFop(*sigma0_specp,'z'))
      error_abort("initial density matrix is not z");
  }
  else {
    if (nchannels>1)
      error_abort("intial density matrix must be set explicitly");
  }

  if (detect_specp) {
    if (!isFop(*detect_specp,'+'))
      error_abort("detection operator is not +");
  }
  else {
    if (nchannels==1)
      detect_specp=new setableoperator_spec(operator_spec(nucids.front(),'+'));
    else
      error_abort("detection operator must be set explicitly");
  }

  cout << "Pulse sequence\n" << actionstack << '\n';

//   Translator* transp=new TranslateInfinity(true);
//   transp->init(fnamebase,actionstack);
  
//   const actionstack_t::const_iterator finish(actionstack.end());
//   actionstack_t::const_iterator start(actionstack.begin());
//   while (start!=finish) {
//     dispatch(transp,*start);
//     start++;
//   }
  
  } catch (const std::exception& exc) {
    cerr << exc.what() << '\n';
    return ERR_FAILED;
  }

  return 0;
}

  
