/* quick and dirty program to read in a spectrum (SIMPSON format) and measure peak height and linewidth */

#include "simpsonio.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

const double BADNORM=1e-2;

double read_file(List<double>& spec,const char *file)
{
  char tmp[128];
  sprintf(tmp,"%s.spe",file);
  List<complex> cspec;
  simpsonFD fd;
  try {
    read_simpson(cspec,fd,tmp);
  }
  catch (...) {
    std::cerr << "Failed to open " << tmp << '\n';
    exit(1);
  }
  spec=real(cspec);
  return fd.sw;
}

double search(const BaseList<double>& spec,int start,int step, double crit =0.5)
{
  const double searchval=crit*spec(start);
  start+=step;
  for (;(start<spec.size()) && (start>=0);start+=step) {
    if (spec(start)<searchval) { //apply linear interpolation
      const double diff1=spec(start-step)-searchval;
      const double diff2=searchval-spec(start);
      return double(start)-step*diff2/(diff1+diff2);
    }
  }
  return -1.0; //flag invalid
}
      
int main(int argc,const char *argv[])
{
  if (argc!=2) {
    cerr << "Syntax: quantify <filename (without .spe)>\n";
    return 1;
  }
  
  List<double> spec;
  const double sw=read_file(spec,argv[1]);

  const size_t npts=spec.size();
  if (npts<1) {
    std::cerr << "Spectrum is empty!\n";
    return 1;
  }
  double maxval=-1e30;  
  size_t wheremax=-1;
  for (size_t i=spec.size();i--;) {
    if (spec(i)>maxval) {
      wheremax=i;
      maxval=spec(i);
    }
  }
  std::cout << "Maximum value: " << maxval << " at index " << wheremax << '\n';
  
  const double offsetplus=search(spec,wheremax,1);
  const double offsetminus=search(spec,wheremax,-1);
  
  if ((offsetplus<0) || (offsetminus<0))
    std::cerr << "Failed to find FWHM\n";
  else {
    const double lw=(offsetplus-offsetminus)*sw/spec.size();
    std::cout << "Full width at half maximum: " << lw << " Hz\n";
  }
  return 0;
}
