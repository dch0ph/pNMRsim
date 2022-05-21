// calculate averaged tensor

#include "space_T.h"
#include "ScratchList.h"
#include "NMR.h"

using namespace libcmatrix;
using namespace std;

double scanfloat(const char* source)
{
  double v;
  if (sscanf(source,"%lf",&v)!=1) {
    cerr << "Invalid floating point number: " << source << '\n';
    exit(1);
  }
  return v;
}
  
double scanlimited(const char* source)
{
  const double v=scanfloat(source);
  if ((v<0.0) || (v>1.0)) {
    cerr << "Invalid order parameter/asymmetry: " << source << '\n';
    exit(1);
  }
  return v;
}

void scanmatrix(Matrix<double>& S, size_t i, size_t j, const char* source)
{  
  S(i,j)=S(j,i)=scanlimited(source);
}
 
// void cartesian_to_PAS_symmetric(double& iso, double& aniso, double& asym, Euler& PAS, const Matrix<double>& A,  bool zeroiso)
// {
//   Matrix<double> V;
//   ScratchList<double,3> eigs(3);
//   hermitian_eigensystem(V,eigs,A);
//   //  std::cout << "Eigenvalues: " << eigs << '\n';
//   //  std::cout << "Eigenvectors\n" << V << '\n';  
//   ScratchList<double,3> tmpeigs(3);
//   const double isotropic=(eigs(0U)+eigs(1U)+eigs(2U))/3.0;
//   for (size_t i=3;i--;)
//     tmpeigs(i)=fabs(eigs(i)-isotropic);
//   //std::cout << "Relative eigenvalues: " << tmpeigs << '\n';
//   ScratchList<size_t,3> order(0U,0U,0U);
//   indexed_sort(order,tmpeigs);
//   std::swap(order(0U),order(1U)); //!< order is y,x,z !
//   //  std::cout << "Order: " << order << '\n';
//   const Matrix<double> Vsorted(V(order,range()));
//   const ScratchList<double,3> sortedeigs(eigs(order));
//   //  std::cout << "Sorted eigenvalues: " << sortedeigs << '\n';
//   //  std::cout << "Sorted eigenvectors\n" << Vsorted << '\n';
//   cartesian_to_anisotropy(iso,aniso,asym,sortedeigs(0U),sortedeigs(1U),sortedeigs(2U));
//   const double sqrt3rd=sqrt(1.0/3.0);
//   const double sqrt2thrd=sqrt(2.0/3.0);
//   iso/=-sqrt3rd;
//   aniso/=sqrt2thrd;
//   //  std::cout << "Isotropic: " << out.iso << "  Anisotropy: " << out.aniso << "  Asymmetry: " << out.asym << '\n';
//   PAS=matrix_to_euler(Vsorted);
//   //std::cout << "Euler angles: " << out.angles << '\n';
// }
 
int main(int argc, const char* argv[])
{
  if (argc!=12) {
    cerr << "Syntax: " << argv[0] << ": <aniso> <asym> <alpha> <beta> <gamma> <Sxx> <Sxy> <Sxz> <Syy> <Syz> <Szz>\n";
    return 1;
  }
  Matrix<double> S(3,3);
  scanmatrix(S,0,0,argv[6]);
  scanmatrix(S,0,1,argv[7]);
  scanmatrix(S,0,2,argv[8]);
  scanmatrix(S,1,1,argv[9]);
  scanmatrix(S,1,2,argv[10]);
  scanmatrix(S,2,2,argv[11]);
  
  const space_T A_PAS(spatial_tensor(scanfloat(argv[1]),scanlimited(argv[2])));
  //cout << "PAS tensor: " << A_PAS << '\n';
  Euler PAS_to_MF(scanfloat(argv[3]),scanfloat(argv[4]),scanfloat(argv[5]));
  PAS_to_MF*=M_PI/180.0; //convert to rad
  const space_T A_MF(rotate(A_PAS,PAS_to_MF));
  Matrix<double> Am_MF;
  tensor_to_matrix(Am_MF,A_MF);
  //cout << "Input tensor: " << A_MF << '\n';
  //cout << "Pre-averaged tensor (molecular frame)\n" << Am_MF << '\n';
  Matrix<double> A_DF;
  multiply(A_DF,S,Am_MF);
  //cout << "Post-averaging tensor\n" << Am_MF << '\n';
  double iso,aniso,asym;
  Euler PAS;
  cartesian_to_PAS_symmetric(iso,aniso,asym,PAS,A_DF);
  PAS*=180/M_PI;
  cout << aniso << ' ' << asym << ' ' << PAS.alpha << ' ' << PAS.beta << ' ' << PAS.gamma << '\n';
  return 0;
}


