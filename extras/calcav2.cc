// calculate averaged tensor over 2 fast wobbles

#include "space_T.h"
#include "NMR.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

void wobble_average(double& retaniso, double& retasym, double aniso1, double asym1, double aniso2, double asym2, double wobble, const char* label, double frac1 =0.5)
{
  const space_T A_PAS1(spatial_tensor(aniso1,asym1));
  cout << "Site 1 tensor at start of step " << label << ": " << A_PAS1 << '\n';
  space_T A_PAS2;
  const bool arediff=(aniso1!=aniso2) || (asym1!=asym2);
  if (arediff) {
    A_PAS2=spatial_tensor(aniso2,asym2);
    cout << "Initial site 2 tensor: " << A_PAS2 << '\n';
  }
  Matrix<double> Am_MF,A_tmp;

  for (size_t step=0;step<2;step++) {
    const double beta=(step==0) ? -wobble/2.0 : wobble/2.0;
    const Euler PAS_to_MF(0.0,beta,0.0);
    const space_T& A_PAS_use( (arediff && (step==1)) ? A_PAS2 : A_PAS1);
    const space_T A_MF(rotate(A_PAS_use,PAS_to_MF));
    tensor_to_matrix(A_tmp,A_MF);
    const double scalefac=(step==0) ? frac1 : 1.0-frac1;
    A_tmp*=scalefac;
    Am_MF+=A_tmp;
  }

  cout << "Averaged tensor\n" << Am_MF << '\n';
  double iso;
  Euler PAS;
  cartesian_to_PAS_symmetric(iso,retaniso,retasym,PAS,Am_MF);
  const double fudge=std::sqrt(1.5); 
  retaniso*=fudge;
  
  cout << "After step " << label << " averaging: aniso=" << retaniso << "  asym=" << retasym << "  orientation: " << PAS << '\n';
}
  
int main(int argc, const char* argv[])
{
  int count=1;
  const double deg_to_rad=M_PI/180.0;
  
  const double aniso=getfloat(argc,argv,count,"Anisotropy (kHz)? ",225.0);
  const double asym=getfloat(argc,argv,count,"Asymmetry? ",0.16);
  const double wobble1=getfloat(argc,argv,count,"Wobble (full) angle 1 (degrees)? ",112.0)*deg_to_rad;
  const double wobble2=getfloat(argc,argv,count,"Wobble (full) angle 2 (degrees) [0 if none]? ",0.0)*deg_to_rad;
  double aniso2=aniso;
  double asym2=asym;
  double frac1=0.5;
  if (wobble2==0.0) {
    frac1=getfloat(argc,argv,count,"Site 1 population fraction? ",0.0,1.0,0.5);
    aniso2=getfloat(argc,argv,count,"Anisotropy of site 2 (kHz)? ",aniso2);
    asym2=getfloat(argc,argv,count,"Anisotropy of site 2 (kHz)? ",asym2);
  }
  double newaniso,newasym;
  wobble_average(newaniso,newasym,aniso,asym,aniso2,asym2,wobble1,"1",frac1);
  if (wobble2!=0.0)
    wobble_average(newaniso,newasym,newaniso,newasym,newaniso,newasym,wobble2,"2");
  
  return 0;
}


