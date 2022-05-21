/* Analyse entanglement build up using pNMRsim framework */

#include "Parser.h"
#include "NMRsim_RF.h"
#include "NMRsim_spinsys.h"
#include "NMRsim_MasterObj.h"

using namespace libcmatrix;
using namespace std;

//! break down into coherence ranks
List<double> analyse(const MasterObj& obj, const BlockedOperator& sigma)
{
  static BlockedMatrix<size_t> coherenceorder;
  static size_t maxrank=0;
  const BlockedMatrix<complex>& sigmar(sigma.row());
  const block_pattern& blkstr(sigma.blockstructure());
  if (!coherenceorder) {
    const SpinOpGenerator& opgen(*(obj.simple_opgenp));
    coherenceorder.duplicate_structure(sigmar);
    const size_t nspins=opgen.nspins();
    List< ListList<bool> > redIzs(nspins);
    for (size_t i=nspins;i--;) {
      ListList<double> Iz;
      opgen.mla_Iz(Iz,1.0,i);
      redIzs(i).duplicate_structure(Iz);
      BaseList<bool> destrow(redIzs(i).row());
      const BaseList<double> sourcerow(Iz.row());
      for (size_t j=sourcerow.size();j--;)
	destrow(j)=(sourcerow(j)<0.0);

      if (verbose & VER_GEN)
	cout << "Iz(" << i << "): " << redIzs(i) << '\n';
    }

    block_pattern::iterator mziter(blkstr);
    size_t bra,ket;
    size_t mzeigSD=0;    
    while (mziter.next(bra,ket)) {
      Matrix<size_t>& curcoher(coherenceorder(mzeigSD));
      for (size_t r=curcoher.rows();r--;) {
	for (size_t s=curcoher.cols();s--;) {
	  size_t coher=0;
	  for (size_t k=redIzs.size();k--;) {
	    const BaseList<bool> brastates(redIzs(k)(bra));
	    const BaseList<bool> ketstates(redIzs(k)(ket));
	    if (brastates(r) ^ ketstates(s))
	      coher++;
	  }
	  if (coher>maxrank)
	    maxrank=coher;
	  curcoher(r,s)=coher;
	}
      }
      mzeigSD++;
    }
    assert(mzeigSD==coherenceorder.size());

    List<size_t> dist(maxrank+1); dist=0U;
    
    for (size_t i=coherenceorder.size();i--;) {
      const BaseList<size_t> curcoher(coherenceorder(i).row());
      for (size_t j=curcoher.size();j--;)
	dist(curcoher(j))++;
    }
    if (verbose & VER_GEN)
      cout << "Coherence order matrix\n" << coherenceorder << '\n';
    const size_t totalstates=sum(dist);
    cout << "Total states available: " << totalstates << '\n';
    cout << "Distribution of states: " << dist << '\n';
    List<double> fdist(dist);
    fdist*=100.0/totalstates;
    cout << "% distribution: " << fdist << '\n' << endl;
  }
  static List<double> cohersum(maxrank+1);
  cohersum=0.0;
  for (size_t i=coherenceorder.size();i--;) {
    const BaseList<size_t> curcoher(coherenceorder(i).row());
    const BaseList<complex> cursigma(sigmar(i).row());
    for (size_t j=curcoher.size();j--;)
      cohersum(curcoher(j))+=norm(cursigma(j));
  }
  if (verbose & VER_GEN) {
    cout << "sigma\n" << sigmar << '\n';
    cout << "Coherence breakdown: " << cohersum << '\n';
  }
  return cohersum;
}

productoperator_spec* simplesigma0_specp=NULL;

void parse_simple_start_operator()
{
  simplesigma0_specp=parse_productoperator();
}

double buildup=-1; //not set by default

void parse_buildup()
{
  buildup=parse_double();
  if (buildup<0.0)
    error_abort("buildup time cannot be <0");
}

int main(int argc, char** argv)
{
  argc-=1;
  argv+=1;

  if (argc && (argv[0][0]=='-'))
    cerr << "Warning: $1 begins with - (" << argv[0] << ").  Misplaced flag?\n";

  declare_builtin_block("spinsys");
  declare_builtin_block("par");

  parser_init("analyseMQ.in",argc,argv,F_SIMPLE);

  command_Factory_t& spinsys_Factory(initialise_spinsys_Factory());
  read_block("spinsys",spinsys_Factory);
  if (!interactions_MFp)
    error_abort("spinsys block failed to defined nuclei / Hamiltonian");
  if (!(interactions_MFp->verify(cout,*cstructp,1e-1)))
    error_abort("Boundary conditions of periodic system incorrect");

  if (verbose)
    dump_interactions();

  if (!(sysp->isspinhalfonly()))
    error_abort("Spin system can only be spin half");

  make_1d_par_variables();
  command_Factory_t& par_Factory(get_par_Factory());  
  par_Factory["start_operator"]=&parse_simple_start_operator;
  par_Factory["buildup"]=&parse_buildup;
  read_par_block();

  const nuclei_spec homonucleus(sysp->homonucleus());
  const List<nuclei_spec> blocknuc(1,homonucleus);
  const bool needrf=(buildup>=0.0);
  MasterObj masterobj(*sysp,*interactions_MFp,needrf ? BaseList<nuclei_spec>() : (BaseList<nuclei_spec>)blocknuc);
  if (needrf) {
    add_channel(homonucleus()); //!< create RF channel (1)
    masterobj.ensure_channels(); //!< create pulse generators
  }
  masterobj.initialise_simple_simulation();

  //!Use Fx if sigma0 not set explicitly
  if (!simplesigma0_specp)
    simplesigma0_specp=new productoperator_spec(operator_spec(homonucleus,'x'));
  else {
    if (buildup>=0.0) //object as we insist on starting with Fx
      error_abort("can't set both buildup time and initial operator");
  }

  //!create initial density matrix
  BlockedOperator sigma(masterobj.make(*simplesigma0_specp));
  const double normfac=norm(sigma.row());
  const double dt=1/sw;
  if (verbose & VER_GEN) {
    cout << "Initial density matrix\n" << sigma << '\n';
    cout << "Normalisation factor: " << normfac << '\n';
  }

  BlockedMatrix<complex> U;
  double t0=0.0;
  if (buildup>0.0) {
    masterobj.propagator(U,0,buildup);
    sigma.unitary_simtrans(U);
    t0+=buildup;    
    Sequence seq;
    const double rf=100e3; //!< nominal RF / Hz
    const double tau90=0.25/rf; //!< 90 length
    seq.push_back(new HardPulse(pulse_generator(),tau90,rf,phase_spec('y')),0.0);
    masterobj.propagator(U,seq,0.0,0.0,0.0);
    if (verbose & VER_GEN)
      cout << "Propagator for 90 pulse:\n" << U << '\n';
    sigma.unitary_simtrans(U);
    if (verbose & VER_GEN)
      cout << "Density matrix after MQ buildup:\n" << sigma << '\n';
  }
  masterobj.propagator(U,t0,t0+dt); //!< dwell time propagator
  if (verbose & VER_GEN)
    cout << "Dwell time propagator\n" << U << '\n';
  
  List<double> resultvec(analyse(masterobj,sigma));
  rmatrix result(np,resultvec.size());
  result.row(0U)=resultvec;
  
  for (size_t i=1;i<np;i++) {
    sigma.unitary_simtrans(U);
    result.row(i)=analyse(masterobj,sigma);
  }
  result*=100.0/normfac;

  cout << "%age entanglement distribution as function of dwell time (";
  prettyprint_time(dt) << ")\n" << result;
  return 0;
}

  
