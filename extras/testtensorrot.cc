/* Test tensor rotation / Euler angle definitions */

#include "ttyio.h"
#include "space_T.h"
#include "geometry.h"
#include "NMR.h"

using namespace libcmatrix;
using namespace std;

#include "./geometry_utils.h"

ordering_convention_t convention=convention_Haeberlen;
 
void dotest(const Matrix<double>& A_MF_cart, Euler_controller& ctrl)
{
  double iso11,naniso22,nasym33;
  Euler PASorient;

  cartesian_to_PAS_symmetric(iso11,naniso22,nasym33,PASorient,A_MF_cart,convention,ctrl,true); //!< apply weird scaling

  if (ctrl.verbose) {
    if (isorderingIUPAC(convention))
      cout << "11: " << iso11 << "  22: " << naniso22 << "  33: " << nasym33 << '\n';
    else
      cout << "Isotropic value: " << iso11 << "  Anisotropy after rotation: " << naniso22 << "  Asymmetry after rotation: " << nasym33 << '\n';
  }

  cout << "Euler angles describing orientation: " << PASorient << " (degrees)\n";
  if (ctrl.verbose) {
    cout << "Euler angles for reverse transformation: " << -PASorient << " (degrees)\n";

    rmatrix fullrotmat;
    rotation_matrix(fullrotmat,PASorient);
    cout << "Recreated rotation matrix:\n" << fullrotmat << '\n';
  }

  const space_T A_newPF(spatial_tensor(iso11,naniso22,nasym33,convention));
  const space_T A_newMF(rotate(A_newPF,PASorient));
  if (ctrl.verbose) {
    std::cout << "Recreated tensor in PAS (spherical tensor rep)\n" << A_newPF << '\n';
    std::cout << "Recreated tensor in MF (spherical tensor rep)\n" << A_newMF << '\n';
  }
  Matrix<double> nA_PAS_cart;
  tensor_to_matrix(nA_PAS_cart,A_newMF);

  testsame(nA_PAS_cart,A_MF_cart,"Recreated tensor in MF (Cartesian representation)",ctrl.verbose);
}

int main(int argc, const char* argv[])
{
  int count=1;
  const double iso=getfloat(argc,argv,count,"Isotropic component? ",0.0);
  const double aniso=getfloat(argc,argv,count,"Anisotropy? ",1);
  const double asym=getfloat(argc,argv,count,"Asymmetry? ",0);

  Euler_controller ctrl;
  ctrl.verbose=getint(argc,argv,count,"Verbosity? ",1);

  const space_T A_PAS(spatial_tensor(iso,aniso,asym));
  if (ctrl.verbose)
    cout << "Tensor in PAS (spherical tensor representation)\n" << A_PAS << '\n';

  Matrix<double> A_PAS_cart;
  tensor_to_matrix(A_PAS_cart,A_PAS);
  if (ctrl.verbose)
    cout << "Tensor in PAS (Cartesian representation)\n" << A_PAS_cart << '\n';

  const double deg_to_rad=M_PI/180.0;
  const double alpha=getfloat(argc,argv,count,"Euler alpha (degrees)? ",0.0)*deg_to_rad;
  const double beta=getfloat(argc,argv,count,"Euler beta (degrees)? ",0.0)*deg_to_rad;
  const double gamma=getfloat(argc,argv,count,"Euler gamma (degrees)? ",0.0)*deg_to_rad;

  std::cout << "H - Haeberlen\nQ - NQR\nI - IUPAC increasing\nD - IUPAC decreasing\n";
  const size_t mode=getoption(argc,argv,count,"Output type? ","HQID",0);
  switch (mode) {
  case 0:
    convention=convention_Haeberlen;
    break;
  case 1:
    convention=convention_NQR;
    break;
  case 2:
    convention=convention_IUPACIncreasing;
    break;
  case 3:
    convention=convention_IUPACDecreasing;
    break;
  default:
    throw InternalError("");
  }
 
  Euler rotn(alpha,beta,gamma);

  Matrix<double> fullrotmat;
  rotation_matrix(fullrotmat,rotn);

  if (ctrl.verbose)
    std::cout << "Rotation matrix determined from Euler angles:\n" << fullrotmat << '\n';

  const space_T A_MF(rotate(A_PAS,rotn));
  if (ctrl.verbose)
    cout << "Tensor in MF (spherical tensor representation) obtained by rotating tensor in PAS by Euler angles\n" << A_MF << '\n';

  Matrix<double> A_MF_cart;
  tensor_to_matrix(A_MF_cart,A_MF);
  cout << "Test tensor in MF (Cartesian representation) from A_MF:\n" << A_MF_cart << '\n';

  if (ctrl.verbose) {
    Matrix<double> AT;
    active_rotate_DCM_tensor(AT,fullrotmat,A_PAS_cart);
    testsame(AT,A_MF_cart,"Tensor in MF - rotated using rotation matrix (R*A_PAS*R^T)",ctrl.verbose);
  }
    
  ctrl.verify_tolerance=1e-6;
  dotest(A_MF_cart,ctrl);

  return(0);
}

  

  
