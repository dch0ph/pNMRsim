#include "simpsonio.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

const double BADNORM=1e-2;

void read_file(List<complex> &FID,const char *file)
{
  char tmp[128];
  sprintf(tmp,"%s.fid",file);
  simpsonFD fd;
  try {
    read_simpson(FID,fd,tmp);
  }
  catch (...) {
    std::cerr << "Failed to open " << tmp << '\n';
    exit(1);
  }
}

int main(int argc,const char *argv[])
{
  int count=1;

  const int nfiles=argc-count;
  if (nfiles<1) {
    cerr << "Need more than 1 file!\n";
    return 1;
  }

  List<complex> FIDbase;
  List<complex> FIDtest;
  double normref;
  int failcount=0;

  for (int i=0;i<nfiles;i++) {
    cout << (i==0 ? "Reading " : "Comparing ") << argv[i+count];
    cout.flush();
    if (i==0) {
      read_file(FIDbase,argv[i+count]);
      normref=sqrt(norm(FIDbase));
      cout << endl;
    }
    else {
      read_file(FIDtest,argv[i+count]);
      FIDtest-=FIDbase;
      const double checkval=sqrt(norm(FIDtest))/normref;
      if (checkval>BADNORM)
      	failcount++;
      cout << ": " << checkval << endl;
    }
  }
  cout << "Failures: " << failcount << '\n';
  return (failcount ? 1 : 0);
}
