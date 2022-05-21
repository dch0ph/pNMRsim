// calc statistics on (homonuclear) Hamiltonian

#include "MetaPropagation.h"
#include "CrystalSystem.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

int main(int argc, const char *argv[])
{
  int count=1;

  const size_t M=getint(argc,argv,count,"Spins per cell? ");
  const size_t N=getint(argc,argv,count,"Unit cells (1D)? ",1);
  const size_t totalspins=M*N;
  if ((totalspins<1) || (totalspins>30)) {
    cerr << "Invalid total number of spins: " << totalspins << '\n';
    return 1;
  }
  
  const spinhalf_system Hsys(M,"1H");
  const CrystalStructure cstruct(N);
  const List<nuclei_spec> blocking(1,nuclei_spec(Hsys(0).nucleus()));
  const CrystalOpGenerator opgen(Hsys,cstruct,blocking);
  const size_t eigs=opgen.eigblocks();
  const size_t mzblks=opgen.mzblocks();
  const size_t totsize=((const SpinOpGeneratorBase&)opgen).size();
  cout << "Overall Hilbert space size: " << totsize << '\n';
  cout << "Eigenvalues: " << eigs << '\n';
  cout << "mz values: " << mzblks << '\n';

  size_t elements_diag=0;
  size_t elements_RF=0;
  size_t total=0;
  size_t tottot=0;
  for (size_t eig=0;eig<eigs;eig++) {
    size_t lastr=0;
    size_t tot=0;
    for (size_t mz=0;mz<mzblks;mz++) {
      const size_t r=opgen.size(mz,eig);
      if (r==0)
	continue;
      cout << "Eigenvalue: " << eig << "  Mz: " << mz << ": " << r << " x " << r << '\n';
      tot+=r;
      elements_diag+=r*r;
      elements_RF+=2*r*lastr;
      lastr=r;
    }
    total+=tot*tot;
    tottot+=tot;
    std::cout << '\n';
  }
  if (tottot!=totsize) {
    cerr << "Totals don't match!\n";
    return 1;
  }

  cout << "Total number of elements: " << total << '\n';
  cout << "Elements in free precession H: " << elements_diag << "  " << (100.0*elements_diag)/total << "%\n";
  elements_RF+=elements_diag;
  cout << "Elements in H with RF: " << elements_RF << "  " << (100.0*elements_RF)/total << "%\n";
  
  return 0;
}
    
