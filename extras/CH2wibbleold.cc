// Calculate 13C,1H dipolar couplings from MD trajectories

#include "ttyio.h"
//#include <cmath>
#include "MoleculeStructure.h"
#include "space_T.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

int verbose = 2;
const double rad2deg=180.0/M_PI;

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

static const double gamma1H(gamma("1H"));
static const double gamma13C(gamma("13C"));

struct Coupling {
  enum coupling_t : int { CH1=0, CH2, HH, ExtraCH };
  
  size_t i,j,cell;
  coupling_t type;
  Matrix<double> tensor;
  double refdist;
  size_t count;
  double gammaX;

  static Matrix<double> Dtmp;

  Coupling(size_t iv, size_t jv, size_t cellv, double refdistv, coupling_t typev) 
    : i(iv), j(jv), cell(cellv), type(typev), refdist(refdistv), count(0U) {
    gammaX = (typev == HH) ? gamma1H : gamma13C;
  }
  
  double add(double& lasym, Euler& lF, const spherical& sphco) {
    static const double sqrt32 = sqrt(3.0/2.0);  //!< sqrt(3/2) is nastiness from forward and backward conversions of tensors

    const double d = dipolar_coupling(gammaX, gamma1H, sphco.r*1e-10);
    if (verbose>0) {
      const char* typelabel = "CH";
      switch (type) {
      case HH:
	typelabel="HH";
	break;
      case ExtraCH:
	typelabel="Long-range CH";
	break;
      default:
	break;
      }
      cout << typelabel << " dipolar coupling between " << (i+1) << " and " << (j+1) << ": " << d << " Hz (r=" << sphco.r << " A)\n";
    }
    const space_T D_PAS(spatial_tensor(d));
    const Euler leuler(vector_to_Euler(sphco));
    const space_T D_MF(rotate(D_PAS,leuler));
    tensor_to_matrix(Dtmp,D_MF);
    tensor += Dtmp;
    if (verbose>1)
      cout << tensor;
    count++;

    const double scale=sqrt32/count;
    double liso,avd;
    cartesian_to_PAS_symmetric(liso,avd,lasym,lF,tensor);
    return avd*scale;
  }
  
};

Matrix<double> Coupling::Dtmp;

int main(int argc, const char* argv[])
{
  int count = 1;
  const double disttolfac=5.0; //!< warn if distance changes by more than this (A)

  const size_t nCatoms=5;
  List<size_t> nCatom(nCatoms);
  nCatom(0U)=0;
  nCatom(1U)=1;
  nCatom(2U)=5;
  nCatom(3U)=8;
  nCatom(4U)=9;

  MoleculeStructure struc(0);

  List< List<Coupling> > couplings(nCatoms);
  List<std::string> labels(nCatoms);

  //  List<size_t> cellstarts;

  //  const double coupling_lim = 150.0; //!< smallest coupling (Hz)
  const double dlim = getfloat(argc, argv, count, "Limiting CH distance (in A)? ",1.6);
  const double dlimCH = 1.5;
  const double dlimHH = 2.0;
  if (dlim < dlimCH) 
    cerr << "Warning: Limiting CH distance cannot be meaningfully less than internal cutoff of " << dlimCH << " A\n"; 
  
  //  std::cout << "Limiting distance: " << dlim << " A\n";
  Matrix<double> Dtmp;

  char fname[256];
  getstring(argc,argv,count,"PDB trajectory file? ",fname,sizeof(fname),"acidch2.pdb");

  FILE* fp=fopen(fname, "r");
  if (fp==NULL) {
    cerr << "Failed to open PDB trajectory file: " << fname << '\n';
    return 1;
  }

  char outfname[256];
  getstring(argc,argv,count,"Output file name? ",outfname,sizeof(outfname),"UICwibble");

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
  enum avmode_t : int { NOAV=0, DOAV, TABALL };
  cout << "N - no average ('central repeat' only)\nA - Average over repeats\nT - Tabulate all averages\n";
  const int avmode = getoption(argc,argv,count,"Averaging mode? ","NAT",NOAV);
  //  const bool multcells = getlogical(argc,argv,count,"Average over repeats? ",false);

  int endframe = 0;
  if (maxframes)
    endframe=startframe+maxframes;

  size_t allnatoms=0;
  //size_t natoms=0;
  int curframe=startframe;
  
  size_t cellnatoms=0;
  size_t ncells=0;

  //List<size_t> nwarnings;
  //List<size_t> Hcounts;
  //List<size_t> H_indices;
  //List<double> sumdssCH2;
  //Matrix<double> sumdCH2;

  FILE* fpconv = NULL;
  FILE* fpconvCH2 = NULL;
  Matrix<double> lastds;
  Matrix<double> lastangles;
  if (avmode!=TABALL) {
    fpconv = createoutput(outfname,"_converge");
    fpconvCH2 = createoutput(outfname,"_CH2converge");
    lastds.create(nCatoms,4);
    lastangles.create(nCatoms,12);
  }

  int nframes=1;
  bool firstloop=true;
  size_t cellstart=0;
  for (;;curframe+=stepframes,nframes++,firstloop=false) {
    if ((endframe>0) && (curframe>=endframe))
      break;

    struc.clear();
    try {
      struc.read_PDB_single(fp, curframe);
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
      cout << "Atoms read: " << allnatoms << "   " << "  Atoms in unit cell: " << cellnatoms << '\n';
      if (avmode!=NOAV) {
	if (allnatoms % cellnatoms) {
	  cerr << "Number of atoms in repeating unit doesn't divide into total number!\n";
	  return 1;
	}
	//!< check that given serial number corresponds to a run of H	
	ncells=allnatoms / cellnatoms;
      }
      else {
	const size_t cell0_start = struc.cell0_start();
	cout << "Serial no. of first atom is " << struc(cell0_start).serial << "  Last atom is " << struc(cell0_start+cellnatoms-1).serial << '\n';
	cellstart = cell0_start;
	ncells=1;
      }
      if (fpconv!=NULL) {
	fprintf(fpconv,"#Frame");
	fprintf(fpconvCH2,"#Frame");
	for (size_t i=0;i<nCatoms;i++) {
	  const char* label = struc(cellstart+nCatom(i)).type;
	  if (label[0] != 'C')
	    throw Failed("Expected to find C atom");
	  fprintf(fpconv,"\t%s_dCH1/kHz\t%s_dCH2/kHz\t%s_dHH/kHz", label,label,label);
	  fprintf(fpconvCH2,"\t%s_dCHrss/kHz", label);
	}
	fprintf(fpconv,"\n");
	fprintf(fpconvCH2,"\n");
      }
    }
    else {
      const size_t cursize=struc.size();
      if (allnatoms!=cursize) {
	cerr << "Frame " << curframe << " differs in size from initial frame!\n";
	return 2;
      }
    }

    if (fpconv!=NULL) {
      fprintf(fpconv,"%i",curframe);
      fprintf(fpconvCH2,"%i",curframe);
    }
    else
      lastds.create(ncells,nCatoms*3);

    for (size_t nC=0;nC<nCatoms;nC++) {
      List<Coupling>& curcouplings(couplings(nC));
	
      if (firstloop) { // first time around search for couplings
	const size_t celli=nCatom(nC);
	labels(nC) = struc(cellstart+celli).type;
	
	for (size_t cell=0;cell<ncells;cell++) {      
	  const size_t ni = cellstart+(cell*cellnatoms)+celli;
	  Coupling::coupling_t coupling_state = Coupling::CH1;
	  List<size_t> CHindices;
	  for (size_t nj = 0; nj<struc.size(); nj++) {
	    if (struc(nj).type[0] != 'H')
	      continue;
	    const vector3 diff(struc(ni).coords - struc(nj).coords);
	    const spherical sphco(diff);
	    if (sphco.r < dlimCH) {
	      if (coupling_state == Coupling::ExtraCH)
		throw Failed("Found >2 matching CH distances");
	      curcouplings.push_back(Coupling(ni, nj, cell, sphco.r, coupling_state));
	      coupling_state =Coupling::coupling_t(int(coupling_state) + 1);
	      CHindices.push_back(nj);
	    }
	    else {
	      if (sphco.r < dlim)
		curcouplings.push_back(Coupling(ni, nj, cell, sphco.r, Coupling::ExtraCH));
	    }
	  }
	  if (CHindices.size()!=2)
	    throw Failed("Didn't find 2 H's per C");	  
	  const size_t nH1=CHindices(0U);
	  const size_t nH2=CHindices(1U);
	  const vector3 diff(struc(nH1).coords - struc(nH2).coords);
	  const spherical sphco(diff);
	  if (sphco.r > dlimHH)
	    throw Failed("Unexpectedly large HH distance");
	  curcouplings.push_back(Coupling(nH1,nH2, cell, sphco.r, Coupling::HH)); 
	}
      }

      List<double> sumdCH2(3U,0.0);
      Matrix<double> sumaux(3U, 4U, 0.0);
      double dss=0.0;
      size_t count=0;
      Euler lF;
      double lasym;
      for (size_t i=curcouplings.size();i--;) {
	Coupling& curcoupling(curcouplings(i));
	const vector3 diff(struc(curcoupling.i).coords - struc(curcoupling.j).coords);
	const spherical sphco(diff);	  
	if (fabs(sphco.r - curcoupling.refdist)>disttolfac) {
	  cerr << "Warnin: internuclear distance between " << (1+curcoupling.i) << " and " << (1+curcoupling.j) << " has changed significantly (" << curcoupling.refdist << " A vs. " << sphco.r << " A\n";
	  //	      nwarnings(cell)++;
	  throw Failed("Bad distance change!");
	}
	const double d = curcoupling.add(lasym, lF, sphco);

	if (curcoupling.type != Coupling::HH)
	  dss += d*d;
	if (curcoupling.type != Coupling::ExtraCH) {
	  sumdCH2(curcoupling.type) += d;
	  sumaux(curcoupling.type, 0U) += lasym;
	  sumaux(curcoupling.type, 1U) += lF.alpha;
	  sumaux(curcoupling.type, 2U) += lF.beta;
	  sumaux(curcoupling.type, 3U) += lF.gamma;
	  if (avmode==TABALL)
	    lastds(curcoupling.cell, 3*nC+size_t(curcoupling.type)) = 1e-3*d;
	  count++;
	}
      }
      if (count != 3*ncells)
	throw Failed("Sanity check failed");
      const double scale = 1e-3/ncells;
      const double drsskHz=1e-3*std::sqrt(dss/ncells);
      if (fpconvCH2) {
	fprintf(fpconvCH2,"\t%g", drsskHz);
	for (size_t i=0;i<3;i++) {
	  const double dkHz= sumdCH2(i)*scale;
	  fprintf(fpconv,"\t%g",dkHz);
	  lastds(nC,i) = dkHz;
	  for (size_t j=0;j<4;j++)
	    lastangles(nC,4*i+j)=sumaux(i,j)/ncells;
	}
	lastds(nC,3U) = drsskHz;
      }
    }
    if (fpconv) {
      fprintf(fpconv,"\n");
      fprintf(fpconvCH2,"\n");
    }
  }
  fclose(fp);
  if (fpconv) {
    fclose(fpconv);
    fclose(fpconvCH2);
  }

  nframes--; //!< last frame was empty - don't count
  cout << "Read " << nframes << " frames\n";

  if (avmode==TABALL) {
    FILE* fpdav = createoutput(outfname,"_dall");
    fprintf(fpdav,"#Cell");
    for (size_t nC=0;nC<nCatoms;nC++) {
      const char* label=labels(nC).c_str();
      fprintf(fpdav,"\t%s_dCH1/kHz\t%s_dCH2/kHz\t%s_dHH/kHz", label,label,label);
    }
    fprintf(fpdav,"\n");
    for (size_t cell=0;cell<ncells;cell++) {
      fprintf(fpdav,"%i",cell+1);
      for (size_t nC=0;nC<nCatoms;nC++)
	for (size_t i=0;i<3;i++)
	  fprintf(fpdav,"\t%g", lastds(cell,3*nC+i));
      fprintf(fpdav,"\n");
    }
  }
  else {
    FILE* fpdav = createoutput(outfname,"_d");
    FILE* fpaux = createoutput(outfname,"_aux");

    fprintf(fpdav,"#Label");
    fprintf(fpaux,"#Label");
    for (size_t i=0;i<3;i++) {
      const char* lab=NULL;
      switch (i) {
      case Coupling::CH1:
	lab="CH1";
	break;
      case Coupling::CH2:
	lab="CH2";
	break;
      case Coupling::HH:
	lab="HH";
      default:
	break;
      }
      fprintf(fpaux,"\t%s_eta\t%s_a/deg\t%s_b/deg\t%s_g/deg", lab, lab, lab, lab);
    }
    fprintf(fpaux,"\n");
    fprintf(fpdav,"#Label\t<dCH1>/kHz\t<dCH2>/kHz\t<dHH>/kHz\t<drssCH>/kHz\n");
    for (size_t nC=0;nC<nCatoms;nC++) {
      fprintf(fpdav,"%s",labels(nC).c_str());
      fprintf(fpaux,"%s",labels(nC).c_str());      
      for (size_t i=0;i<4;i++)
	fprintf(fpdav,"\t%g", lastds(nC,i));
      for (size_t i=0;i<3;i++) {
	fprintf(fpaux,"\t%g", lastangles(nC,4*i));
	for (size_t j=1;j<4;j++) 
	  fprintf(fpaux,"\t%g", rad2deg*lastangles(nC,4*i+j));
      }
      fprintf(fpdav,"\n");
      fprintf(fpaux,"\n");
    }
  }
      
  //  cout << "Distance change warnings per cell " << nwarnings << '\n';
  return 0;
}
