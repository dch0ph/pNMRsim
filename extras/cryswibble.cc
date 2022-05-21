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

size_t validatecellstarts(const MoleculeStructure& struc, size_t start, size_t natoms) {
  const int start_index = struc.snumber_to_index(start); //!< may crash if start out of range
  if (start_index<0) {
    cerr << "Serial number " << start << " not present!\n";
    exit(1);
  }
  for (size_t i=natoms;i--;) {
    if (struc(start_index+i).type[0]!='H') {
      cerr << "Atom serial number " << struc(start_index+i).serial << " is not an H!\n";
      exit(1);
    }
  }
  return start_index;
}

int main(int argc, const char* argv[])
{
  int count = 1;
  const double disttolfac=12.0; //!< warn if distance changes by more than this (A)
  
  MoleculeStructure struc(0);
  const char selatom = 'H';

  Matrix< Matrix<double> > tensors;
  Matrix<double> refdist;
  //  Matrix<bool> included;
  List<double> last_drss;
  List<std::string> labels;
  List<size_t> cellstarts;
  const double gamma1H(gamma("1H"));
  const double coupling_lim = 150.0; //!< smallest coupling (Hz)
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
  const int stepframes = getint(argc,argv,count,"Frame step? ",1);
  if (stepframes<1) {
    cerr << "Step frames can't be <1\n";
    return 1;
  }
  const bool multcells = getlogical(argc,argv,count,"Average over multiple? ",false);

  int endframe = 0;
  if (maxframes)
    endframe=startframe+maxframes;

  size_t allnatoms=0;
  size_t natoms=0;
  int curframe=startframe;
  
  Euler lF;
  double liso,d,lasym;
  size_t cellnatoms=0;
  size_t ncells=0;
  List<size_t> nwarnings;

  FILE* fpconv = createoutput(outfname,"_converge");
  fprintf(fpconv,"#Frame\trootmean_drss/kHz\n");
  
  int nframes=1;
  bool firstloop=true;
  for (;;curframe+=stepframes,nframes++,firstloop=false) {
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

    cellnatoms = struc.cell0_end() - struc.cell0_start();

    //! sanity check on structure size - ideally would check all parameters
    if (firstloop) {
      allnatoms=struc.size();
      cout << selatom << " atoms read: " << allnatoms << "   " << selatom << " atoms in unit cell: " << cellnatoms << '\n';
      if (multcells) {
	//!< check that given serial number corresponds to a run of H	
	ncells=4;
	cellstarts.create(4);
	cellstarts(0U)= validatecellstarts(struc, 5693, cellnatoms);
	cellstarts(1U)= validatecellstarts(struc, 5727, cellnatoms);
	cellstarts(2U)= validatecellstarts(struc, 5761, cellnatoms);
	cellstarts(3U)= validatecellstarts(struc, 5795, cellnatoms);
      }
      else {
	const size_t cell0_start = struc.cell0_start();
	cout << "Serial no. of first H is " << struc(cell0_start).serial << "  Last H is " << struc(cell0_start+cellnatoms-1).serial << '\n';
	cellstarts.create(1);
	cellstarts.front()=cell0_start;
	ncells=1;
      }
      nwarnings.create(ncells);
      nwarnings=0U;
      natoms=ncells*cellnatoms;
      tensors.create(natoms,allnatoms);
      refdist.create(natoms,allnatoms);
      labels.create(natoms);
      last_drss.create(natoms,0.0);
      size_t base=0;
      for (size_t cell=0;cell<ncells;cell++) {
	const size_t cell0_start = cellstarts(cell);
	for (size_t i=cellnatoms;i--;)
	  labels(base+i) = struc(cell0_start+i).type;
	base+=cellnatoms;
      }
    }
    else {
      const size_t cursize=struc.size();
      if (allnatoms!=cursize) {
	cerr << "Frame " << curframe << " differs in size from initial frame!\n";
	return 2;
      }
    }

    size_t inccount = 0;
    size_t base=0;
    for (size_t cell=0;cell<ncells;cell++) {      
      for (size_t celli=0;celli<cellnatoms;celli++) {
	const size_t outi=base+celli;
	const size_t ni = cellstarts(cell)+celli;
	for (size_t nj = struc.size(); nj--;) {
	  if (ni == nj)
	    continue;
	  if (!firstloop && (refdist(outi, nj)==0.0))
	    continue;
	  const vector3 diff(struc(ni).coords - struc(nj).coords);
	  const spherical sphco(diff);
	  if (firstloop) {
	    const bool incl = (sphco.r < dlim);
	    refdist(outi, nj) = incl ? sphco.r : 0.0; //!< store distance if included, otherwise 0
	    if (!incl)
	      continue;
	    inccount++;
	  }
	  else {
	    const double r0 = refdist(outi, nj);
	    if (fabs(sphco.r-r0)>disttolfac) {
	      cerr << "Warning: internuclear distance between " << ni << " and " << nj << " has changed significantly (" << r0 << " A vs. " << sphco.r << " A\n";
	      nwarnings(cell)++;
	    }
	  }
	  d = dipolar_coupling(gamma1H, gamma1H, sphco.r*1e-10);
	  if (verbose>0)
	    cout << "Dipolar coupling between " << ni << " and " << nj << ": " << d << " Hz\n";
	  const space_T D_PAS(spatial_tensor(d));
	  const Euler leuler(vector_to_Euler(sphco));
	  const space_T D_MF(rotate(D_PAS,leuler));
	  tensor_to_matrix(Dtmp,D_MF);
	  tensors(outi, nj) += Dtmp;
	}
      }
      base+=cellnatoms;
    }
    if (firstloop)
      std::cout << "Included dipolar couplings: " << inccount << " out of " << natoms*(allnatoms-1) << '\n';

    const double scale=sqrt(3.0/2.0)/nframes; //!< sqrt(3/2) is nastiness from forward and backward conversions of tensors
    base=0;
    fprintf(fpconv,"%i",curframe);
    double sumcells=0.0;
    for (size_t cell=0;cell<ncells;cell++) {      
      double sumsumdss=0.0;
      for (size_t celli=0;celli<cellnatoms;celli++) {
	const size_t outi=base+celli;
	double sumdss=0.0;
	const size_t ni = cellstarts(cell)+celli;
	//cout << "Out atom: " << outi << "\n";
	for (size_t nj = allnatoms; nj--;) {
	  Matrix<double>& curtensor(tensors(outi,nj));
	  if (!curtensor)
	    continue;
	  cartesian_to_PAS_symmetric(liso,d,lasym,lF,tensors(outi,nj));
	  d*=scale;
	  if (verbose>0)
	    cout << "Averaged dipolar coupling between " << ni << " and " << nj << ": " << d << " Hz\n";
	  sumdss+=d*d;
	}
	const double drss = sqrt(sumdss);
	last_drss(outi)=drss;
	if (verbose>0)
	  cout << "d_rss for " << labels(outi) << ": " << (drss*1e-3) << " kHz\n";
	sumsumdss+=sumdss;
      }
      const double rootmean_dss=sqrt(sumsumdss/cellnatoms)*1e-3;
      cout << "root mean d_ss: " << rootmean_dss << " kHz\n";
      fprintf(fpconv," %g",rootmean_dss);
      sumcells+=rootmean_dss;
      base+=cellnatoms;
    }
    if (multcells)
      fprintf(fpconv," %g\n", (sumcells/ncells));
    else
      fprintf(fpconv,"\n");
  }
  fclose(fp);
  fclose(fpconv);

  nframes--; //!< last frame was empty - don't count
  cout << "Read " << nframes << " frames\n";

  FILE* fpdrss = createoutput(outfname,"_drss");
  fprintf(fpdrss,"#Label\td_rss/kHz\n");
  for (size_t i=0;i<natoms;i++)
    fprintf(fpdrss,"%s %g\n",labels(i).c_str(),last_drss(i));

  std::cout << "Distance change warnings per cell " << nwarnings << '\n';
  return 0;
}
