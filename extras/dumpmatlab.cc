#include "matlabio.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

template<class Ctrl, class T> void dump(const char* name, Ctrl& ctrl, Matrix<T>& a)
{
  ctrl.read(a);
  if (a.rows()>1) {
    if (*name)
      cout << name << '\n';
    cout << a << '\n';
  }
  else {
    if (*name)
      cout << name << ": ";
    if (a.cols()>1)
      cout << a.row() << '\n';
    else
      cout << a(0U,0U) << '\n';
  }
}

template<class T> void dump(T& ctrl)
{
  List< Matrix<double> > blockedr;
  List< Matrix<complex> > blockedc;
  Matrix<double> matrixr;
  Matrix<complex> matrixc;
  List<char> str;
  matlab_controller::header_info info;
  
  while (ctrl.peek(info)) {
    switch (info.type) {
//     case matlab_controller::CELL: //!< assumes cell is being used to store identical objects (use matlab_controller::composite otherwise)
//       if (info.iscomplex) {
// 	ctrl.read(blockedc);
// 	cout << info.name << '\n' << blockedc << '\n';
//       }
//       else {
// 	ctrl.read(blockedr);
// 	cout << info.name << '\n' << blockedr << '\n';
//       }
//       break;
    case matlab_controller::ARRAY:
      if (info.dims.size()>2) {
	cout << info.name << " is multidimensional (skipping)\n";
	ctrl.next();
      }
      else {
	if (info.iscomplex)
	  dump(info.name,ctrl,matrixc);
	else
	  dump(info.name,ctrl,matrixr);
      }
      break;
    case matlab_controller::CHAR:
      ctrl.read(str);
      cout << info.name << ": " << str.vector() << '\n';
      break;
    case matlab_controller::STRUCT: case matlab_controller::CELL: {
      cout << info.name << " (structure/cell):\n";
      matlab_controller::composite comp(ctrl);
      dump(comp);
    }
      break;
    default:
      std::cerr << "<Unhandled object type>\n";
    }
  }
}
  
int main(int argc, const char* argv[])
{
  char fname[256];
  int count=1;
  getstring(argc,argv,count,"File name (no .mat)? ",fname,sizeof(fname));

  try {
    matlab_controller ctrl(fname);
    dump(ctrl);
  } catch (MatrixException& exc) {
    cerr << exc << '\n';
    return 1;
  }
  return 0;
}
