// Calculate 1H 2nd moments etc. from MD trajectories

#include "ttyio.h"
#include "MoleculeStructure.h"
#include "space_T.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

int verbose = 0;

FILE* createoutput(const char* basename, const char* qual)
{
  char scratch[256];
  snprintf(scratch, sizeof(scratch), "%s%s.txt",basename,qual);
  FILE* fp=fopen(scratch, "wa");
  if (fp==NULL) {
    cerr << "Failed to open output file: " << scratch << "\n";
    exit(2);
  }
  return fp;
}

int main(int argc, const char* argv[])
{
  int count = 1;
  
  MoleculeStructure struc;
  const char selatom = 'H';

  Matrix< Matrix<double> > tensors;
  Matrix<bool> included;
  List<double> last_drss;
  List<std::string> labels;
  const double gamma1H(gamma("1H"));
  const double coupling_lim = 50.0; //!< smallest coupling (Hz)
  const double dlim = 1e10*dipolar_coupling_to_r(coupling_lim, gamma1H, gamma1H);
  std::cout << "Limiting distance: " << dlim << " A\n";
  Matrix<double> Dtmp;

  char fname[256];
  getstring(argc,argv,count,"PDB trajectory file? ",fname,sizeof(fname),"Paul-300.pdb");

  FILE* fp=fopen(fname, "r");
  if (fp==NULL) {
    cerr << "Failed to open PDB trajectory file: " << fname << '\n';
    return 1;
  }

  char outfname[256];
  getstring(argc,argv,count,"Output file name? ",outfname,sizeof(outfname),"diamantane_wibble");

  const int startframe = getint(argc,argv,count,"Start frame? ",1);
  const int maxframes = getint(argc,argv,count,"Frames to include (0 for all)? ",0);
  if (maxframes<0) {
    cerr << "Frames can't be <0!\n";
    return 1;
  }
  int endframe = 0;
  if (maxframes)
    endframe=startframe+maxframes;

  size_t allnatoms=0;
  size_t natoms=0;
  size_t cell0_start=0;
  int curframe=startframe;
  
  Euler lF;
  double liso,d,lasym;

  FILE* fpconv = createoutput(outfname,"_converge");
  fprintf(fpconv,"#Frame\trootmean_drss/kHz\n");
  
  int nframes=1;
  bool firstloop=true;
  for (;;curframe++,nframes++,firstloop=false) {
    if ((endframe>0) && (curframe>=endframe))
      break;

    struc.clear();
    try {
      struc.read_PDB_single(fp, curframe, selatom);
      if (struc.empty())
	throw Failed("Empty structure!");
    }
    catch (const MatrixException& exc) {
      if (curframe>1)
	break;

      cerr << "Failed to read any frames from " << fname << ": " << exc;
      return 1;
    }

    cell0_start = struc.cell0_start();
    natoms = struc.cell0_end() - cell0_start;

    //! sanity check on structure size - ideally would check all parameters
    if (firstloop) {
      allnatoms=struc.size();
      cout << selatom << " atoms read: " << allnatoms << "   " << selatom << " atoms in unit cell: " << natoms << '\n';
      cout << "Serial no. of first H is " << struc(cell0_start).serial << "  Last H is " << struc(cell0_start+natoms-1).serial << '\n';
      tensors.create(natoms,allnatoms);
      included.create(natoms,allnatoms);
      labels.create(natoms);
      last_drss.create(natoms,0.0);
      for (size_t i=natoms;i--;)
	labels(i) = struc(cell0_start+i).type;
    }
    else {
      const size_t cursize=struc.size();
      if (allnatoms!=cursize) {
	cerr << "Frame " << curframe << " differs in size from initial frame!\n";
	return 2;
      }
    }    

    size_t inccount = 0;
    for (size_t i=0;i<natoms;i++) {
      const size_t ni = cell0_start+i;
      for (size_t nj = struc.size(); nj--;) {
	if (ni == nj)
	  continue;
	if (!firstloop && !included(i, nj))
	  continue;
	const vector3 diff(struc(ni).coords - struc(nj).coords);
	const spherical sphco(diff);
	if (firstloop) {
	  const bool incl = (sphco.r < dlim);
	  included(i, nj) = incl;
	  if (!incl)
	    continue;
	  inccount++;
	}
	d = dipolar_coupling(gamma1H, gamma1H, sphco.r*1e-10);
	if (verbose>0)
	  cout << "Dipolar coupling between " << ni << " and " << nj << ": " << d << " Hz\n";
	const space_T D_PAS(spatial_tensor(d));
	const Euler leuler(vector_to_Euler(sphco));
	const space_T D_MF(rotate(D_PAS,leuler));
	tensor_to_matrix(Dtmp,D_MF);
	tensors(i, nj) += Dtmp;
      }
    }
    
    if (firstloop)
      std::cout << "Included dipolar couplings: " << inccount << " out of " << natoms*(allnatoms-1) << '\n';

    const double scale=sqrt(3.0/2.0)/nframes; //!< sqrt(3/2) is nastiness from forward and backward conversions of tensors
    double sumsumdss=0.0;
    for (size_t i=0;i<natoms;i++) {
      double sumdss=0.0;
      const size_t ni = cell0_start+i;
      for (size_t nj = allnatoms; nj--;) {
	Matrix<double>& curtensor(tensors(i,nj));
	if (!curtensor)
	  continue;
	cartesian_to_PAS_symmetric(liso,d,lasym,lF,tensors(i,nj));
	d*=scale;
	if (verbose>0)
	  cout << "Averaged dipolar coupling between " << ni << " and " << nj << ": " << d << " Hz\n";
	sumdss+=d*d;
      }
      const double drss = sqrt(sumdss);
      last_drss(i)=drss;
      if (verbose>0)
	cout << "d_rss for " << labels(i) << ": " << (drss*1e-3) << " kHz\n";
      sumsumdss+=sumdss;
    }
    const double rootmean_dss=sqrt(sumsumdss/natoms)*1e-3;
    cout << "root mean d_ss: " << rootmean_dss << " kHz\n";
    fprintf(fpconv,"%i %g\n",curframe,rootmean_dss);

  }
  fclose(fp);
  fclose(fpconv);

  nframes--; //!< last frame was empty - don't count
  cout << "Read " << nframes << " frames\n";

  FILE* fpdrss = createoutput(outfname,"_drss");
  fprintf(fpdrss,"#Label\td_rss/kHz\n");
  for (size_t i=0;i<natoms;i++)
    fprintf(fpdrss,"%s %g\n",labels(i).c_str(),last_drss(i));

  return 0;
}
