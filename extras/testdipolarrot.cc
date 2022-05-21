/* Test tensor rotation / Euler angle definitions */

#include "ttyio.h"
#include "space_T.h"
#include "geometry.h"
#include "NMR.h"

using namespace libcmatrix;
using namespace std;

#include "./geometry_utils.h"

int main(int argc, const char* argv[])
{
  int count=1;
  const double aniso=getfloat(argc,argv,count,"Anisotropy? ",1);
  //const double asym=getfloat(argc,argv,count,"Asymmetry? ",0);

  std::cout.precision(8);
  Euler_controller ctrl;
  ctrl.verbose=getint(argc,argv,count,"Verbosity? ",1);

  const space_T A_PAS(spatial_tensor(aniso));
  if (ctrl.verbose)
    cout << "Dipolar tensor in PAS (spherical tensor representation)\n" << A_PAS << '\n';

  Matrix<double> A_PAS_cart;
  tensor_to_matrix(A_PAS_cart,A_PAS);
  if (ctrl.verbose)
    cout << "Dipolar tensor in PAS (Cartesian representation)\n" << A_PAS_cart << '\n';

  const double deg_to_rad=M_PI/180.0;
  //const double alpha=getfloat(argc,argv,count,"Euler alpha (degrees)? ",0.0)*deg_to_rad;
  const double beta=getfloat(argc,argv,count,"Euler beta (degrees)? ",0.0)*deg_to_rad;
  const double gamma=getfloat(argc,argv,count,"Euler gamma (degrees)? ",0.0)*deg_to_rad;
  Euler rotn(0.0,beta,gamma);

  rmatrix3 R;
  rotation_matrix(R,rotn);
  Matrix<double> fullrotmat;
  rmatrix3_to_rmatrix(fullrotmat,R);
  const vector3 zaxis(0,0,1.0);
  const vector3 newpos(rotate(zaxis,R));

  Matrix<double> A_MF_cart;
  active_rotate_DCM_tensor(A_MF_cart,fullrotmat,A_PAS_cart);

  if (ctrl.verbose) {
    std::cout << "Rotation matrix determined from Euler angles:\n" << fullrotmat << '\n';
    std::cout << "New position of vector initially oriented along z: " << newpos << '\n';
    std::cout << "Dipolar tensor in MF:\n" << A_MF_cart << '\n';
  }

  const spherical sphco(newpos);
  const Euler neweuler(vector_to_Euler(sphco));
  std::cout << "Recreated Euler angles from geometry: " << neweuler << '\n';

  Matrix<double> Rn;
  rotation_matrix(Rn,neweuler);
  testsame(Rn,fullrotmat,"Rotation matrix from Euler angles deduced from geometry",ctrl.verbose);

//   Euler altneweuler(0.0,sphco.theta,sphco.phi);
//   std::cout << "Original euler angles: " << altneweuler << '\n';
//   rotation_matrix(Rn,altneweuler);
//   testsame(Rn,fullrotmat,"Rotation matrix from Euler angles deduced from geometry (original)",ctrl.verbose);

  Euler PASorient;
  double iso,naniso,nasym;
  cartesian_to_PAS_symmetric(iso,naniso,nasym,PASorient,A_MF_cart,convention_Haeberlen,ctrl);
  std::cout << "Recreated Euler angles from rotation tensor: " << PASorient << '\n';  
  
  Matrix<double> nfullrotmat;
  rotation_matrix(nfullrotmat,neweuler);

  Matrix<double> nA_MF_cart;
  active_rotate_DCM_tensor(nA_MF_cart,nfullrotmat,A_PAS_cart);

  const space_T A_MF(rotate(A_PAS,neweuler));
  if (ctrl.verbose)
    cout << "Dipolar tensor in MF obtained by rotating tensor in PAS by Euler angles deduced from geometry\n" << A_MF << '\n';

  testsame(nA_MF_cart,A_MF_cart,"Dipolar tensor in MF obtained by rotating tensor by Euler anglers deduced from geometry",ctrl.verbose);

  return(0);
}
