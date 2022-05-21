// Calculate 13C,1H dipolar couplings from MD trajectories

#include "ttyio.h"
//#include <cmath>
#include "MoleculeStructure.h"
#include <utility>
#include <map>
#include "space_T.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

bool recentre = true;
bool includelattice = false;
int verbose = 0;
const double rad2deg=180.0/M_PI;
const double deg2rad=M_PI/180.0;
const double gimbal_tol = 4.0*deg2rad; //! force to z

Euler_controller euctrl;

//void lattice_vectors(Matrix<double>& vecs, double a, double b, double c, double alpha, double beta, double gamma)
void lattice_vectors(Matrix<double>& vecs, const lattice_info_t& lattice)
{
  // adapted from TriclinicFactory in ASE
  double cosa = cos(lattice.alpha*deg2rad);
  if (fabs(cosa)<1e-7)
    cosa = 0.0;
  double cosb = cos(lattice.beta*deg2rad);
  if (fabs(cosb)<1e-7)
    cosb = 0.0;
  const double sinb = sin(lattice.beta*deg2rad);
  const double cosg = cos(lattice.gamma*deg2rad);
  const double sing = sin(lattice.gamma*deg2rad);
  if (sinb==0.0)
    throw InvalidParameter("lattice_vectors: invalid beta angle");
  if (sing==0.0)
    throw InvalidParameter("lattice_vectors: invalid gamma angle");
  const double angfac = (cosa - cosb*cosg)/sing;
  vecs.create(3U,3U,0.0);
  vecs(0U,0U) = lattice.a;
  vecs(1U,0U) = lattice.b*cosg;
  vecs(1U,1U) = lattice.b*sing;
  vecs(2U,0U) = lattice.c*cosb;
  if (angfac) {
    vecs(2U,1U) = lattice.c*angfac;
    vecs(2U,2U) = lattice.c*sqrt(sinb*sinb - angfac*angfac);
  }
  else
    vecs(2U,2U) = lattice.c*fabs(sinb);
}

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

class DipolarWriter 
{
 public:
  typedef std::pair<vector3, std::string>  dipolar_info_t;
  typedef std::map< std::string, List<dipolar_info_t> > map_t;

  DipolarWriter() { clear(); }

  ~DipolarWriter() {
    close();
  }
  void close() {
    if (fp) {
      fclose(fp);
      fp=NULL;
    }
  }
  void clear() {
    data.clear();
    sum_ = vector3(0.0, 0.0, 0.0);
    count_ = 0U;
  }

  void add(const char* label, const vector3& pos, const Matrix<double>& tensor, double scale =1.0)
  {
    //char posstr[128];
    //vector3 npos(pos);
    //npos -= centre_;
    //    snprintf(posstr,sizeof(posstr),"%g %g %g", npos.x, npos.y, npos.z);
    char tenstr[128];
    snprintf(tenstr,sizeof(tenstr),"%g %g %g %g %g %g %g %g %g", scale*tensor(0U,0U), scale*tensor(0U,1U), scale*tensor(0U, 2U),  scale*tensor(1U, 0U), scale*tensor(1U, 1U), scale*tensor(1U, 2U), scale*tensor(0U, 2U), scale*tensor(1U, 2U), scale*tensor(2U, 2U));
    data[label].push_back(dipolar_info_t(pos, tenstr));
    //    std::cout << pos << '\n';
    sum_ += pos;
    count_++;
  }

  void write(const char* fname, const MoleculeStructure& struc) {
    fp=fopen(fname, "wa");
    if (fp==NULL)
      throw Failed("DipolarWriter: file open");
    if (data.empty())
      throw Failed("DipolarWriter: nothing to write!");
    vector3 midp = sum_;
    midp *= 1.0/count_;
    //    std::cout << "Sum: " << sum_ << "  Count: " << count_ << "  Mid point: " << midp << '\n';
    fprintf(fp,"#$magres-abinitio-v1.0\n");
    fprintf(fp,"[atoms]\n");
    if (includelattice) {
      if (struc.has_lattice()) {
	Matrix<double> latticevecs;
	lattice_vectors(latticevecs, struc.lattice_parameters());
	fprintf(fp,"units lattice Angstrom\n");      
	fprintf(fp,"lattice");
	for (size_t i=0;i<3;i++)
	  for (size_t j=0;j<3;j++)
	    fprintf(fp," %g",latticevecs(i,j));
	fprintf(fp,"\n");
      }
      else
	std::cerr << "Warning: no lattice information found\n";
    }
    fprintf(fp,"units atom Angstrom\n");
    for (auto ptr=data.cbegin();ptr!=data.cend();++ptr) {
      const char* label = (ptr->first).c_str();
      const List<dipolar_info_t>& dipinfol(ptr->second);
      for (int i=0;i<dipinfol.length();i++) {
	vector3 vec(dipinfol[i].first);
	if (recentre)
	  vec -= midp;
	fprintf(fp,"atom %c %s %i %g %g %g\n", label[0], label, i+1, vec.x, vec.y, vec.z);
      }
    }
    fprintf(fp,"[/atoms]\n");
    fprintf(fp,"[magres]\n");    
      fprintf(fp,"units ms ppm\n");
    for (auto ptr=data.cbegin();ptr!=data.cend();++ptr) {
      const char* label = (ptr->first).c_str();
      const List<dipolar_info_t>& dipinfol(ptr->second);
      for (int i=0;i<dipinfol.length();i++)
	fprintf(fp,"ms %s %i %s\n", label, i+1, dipinfol[i].second.c_str());
    }
    fprintf(fp,"[/magres]\n");    
    close();
  }

//   const vector3& centre() { return centre_; }
//   void centre(const vector3& centrev) {
//     if (!(data.empty()))
//       throw Failed("Can't move centre part-way through adding data!");
//     centre_ = centrev;
//   }
  
 private:
  FILE* fp;
  map_t data;
  vector3 sum_;
  size_t count_;
};

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
    return diag(lasym,lF,tensor,1.0/count);
  }

  static double diag(double& lasym, Euler& lF, const Matrix<double>& tens, double scale =1.0) {
    static const double sqrt32 = sqrt(3.0/2.0);  //!< sqrt(3/2) is nastiness from forward and backward conversions of tensors
    double liso,avd;
    euctrl.gimbal_tolerance = gimbal_tol;
    cartesian_to_PAS_symmetric(liso,avd,lasym,lF,tens, convention_NQR, euctrl);    
    return avd*scale*sqrt32;
  }
  
  const Matrix<double>& value() const { return tensor; }
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
  //  Matrix<Matrix<double> > sumtensor(nCatoms,3U);

  if (avmode!=TABALL) {
    fpconv = createoutput(outfname,"_converge");
    fpconvCH2 = createoutput(outfname,"_CH2converge");
    lastds.create(nCatoms,4);
    lastangles.create(nCatoms,12);
  }

  DipolarWriter* dipwritep=NULL;
  if (avmode==NOAV)
    dipwritep = new DipolarWriter();

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

    if (dipwritep) {
      dipwritep->clear();
      //  dipwritep->centre(struc.centre());
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
	  fprintf(fpconv,"\t%s_dCH1/kHz\t%s_dCH2/kHz\t%s_dHH/kHz\t%s_dCH1eta\t%s_dCH2eta\t%s_dHHeta", label,label,label, label,label,label);
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
    else {
      lastds.create(ncells,nCatoms*3);
      lastangles.create(ncells,nCatoms*6);
    }

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
      List<double> sumetaCH2(3U,0.0);
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
	  cerr << "Warning: internuclear distance between " << (1+curcoupling.i) << " and " << (1+curcoupling.j) << " has changed significantly (" << curcoupling.refdist << " A vs. " << sphco.r << " A\n";
	  //	      nwarnings(cell)++;
	  throw Failed("Bad distance change!");
	}
	const double d = curcoupling.add(lasym, lF, sphco);	
	if (dipwritep) {
	  char lab[12];

	  switch (curcoupling.type) {
	  case Coupling::CH1: case Coupling::CH2:
	    snprintf(lab,sizeof(lab),"%sH%i", struc(curcoupling.i).type, (curcoupling.type==Coupling::CH1) ? 1 : 2);
	    dipwritep->add(lab, struc(curcoupling.j).coords, curcoupling.value(), 1e-3/curcoupling.count);
	    break;
	  case Coupling::HH: {
	    vector3 pos(struc(curcoupling.i).coords);
	    pos+=struc(curcoupling.j).coords;
	    pos*=0.5;
	    dipwritep->add("HH", pos, curcoupling.value(), 1e-3/curcoupling.count);
	  }
	    break;
	  default:
	    break;
	  }
	}
	if (curcoupling.type != Coupling::HH)
	  dss += d*d;
	if (curcoupling.type != Coupling::ExtraCH) {
	  sumdCH2(curcoupling.type) += d;
	  sumetaCH2(curcoupling.type) += lasym;
	  //	  sumtensor(nC,curcoupling.type) += curcoupling.tensor;
	  //	  cout << curcoupling.type << ": " << lasym << '\n';
	  sumaux(curcoupling.type, 0U) += lasym;
	  sumaux(curcoupling.type, 1U) += lF.alpha;
	  sumaux(curcoupling.type, 2U) += lF.beta;
	  sumaux(curcoupling.type, 3U) += lF.gamma;
	  if (avmode==TABALL) {
	    lastds(curcoupling.cell, 3*nC+size_t(curcoupling.type)) = 1e-3*d;
	    const size_t base = 6*nC+2*size_t(curcoupling.type);
	    lastangles(curcoupling.cell, base) = lasym;
	    lastangles(curcoupling.cell, base+1) = lF.alpha;
	  }
	  count++;
	}
      }
      if (count != 3*ncells)
	throw Failed("Sanity check failed");
      const double scale = 1.0/ncells;
      const double drsskHz=1e-3*std::sqrt(dss/ncells);
      if (fpconvCH2) {
	fprintf(fpconvCH2,"\t%g", drsskHz);
	for (size_t i=0;i<3;i++) {
	  const double dkHz= 1e-3*sumdCH2(i)*scale;
	  fprintf(fpconv,"\t%g",dkHz);
	  lastds(nC,i) = dkHz;
	  for (size_t j=0;j<4;j++)
	    lastangles(nC,4*i+j)=sumaux(i,j)/ncells;
	}
	for (size_t i=0;i<3;i++)
	  fprintf(fpconv,"\t%g",scale*sumetaCH2(i));
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
  if (dipwritep) {
    char scratch[256];
    snprintf(scratch,sizeof(scratch),"%s.magres",outfname);
    dipwritep->write(scratch, struc);
  }

  nframes--; //!< last frame was empty - don't count
  cout << "Read " << nframes << " frames\n";

  if (avmode==TABALL) {
    FILE* fpdav = createoutput(outfname,"_dall");
    FILE* fpaux = createoutput(outfname,"_angall");
    fprintf(fpdav,"#Cell");
    fprintf(fpaux,"#Cell");
    for (size_t nC=0;nC<nCatoms;nC++) {
      const char* label=labels(nC).c_str();
      for (size_t i=0;i<3;i++) {
	const char* tlab = NULL;
	switch (i) {
	case Coupling::CH1:
	  tlab = "CH1";
	  break;
	case Coupling::CH2:
	  tlab = "CH2";
	  break;
	case Coupling::HH:
	  tlab = "HH";
	default:
	  break;
	}
	fprintf(fpdav,"\t%s_%s/kHz", label, tlab);
	fprintf(fpaux,"\t%s_%seta\t%s_%sa/deg", label, tlab, label, tlab);
      }
    }
    fprintf(fpdav,"\n");
    for (size_t cell=0;cell<ncells;cell++) {
      fprintf(fpdav,"%i",cell+1);
      fprintf(fpaux,"%i",cell+1);
      for (size_t nC=0;nC<nCatoms;nC++) {
	for (size_t i=0;i<3;i++) {
	  fprintf(fpdav,"\t%g", lastds(cell,3*nC+i));
	  const size_t base = 6*nC+2*i;
	  fprintf(fpaux,"\t%g\t%g", lastangles(cell,base),rad2deg*lastangles(cell,base+1));
	}
      }
      fprintf(fpdav,"\n");
      fprintf(fpaux,"\n");
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
//       double lasym;
//       Euler lF;

//       for (size_t i=0;i<3;i++) {
// 	(void)Coupling::diag(lasym,lF,sumtensor(nC,i),1.0/(ncells*nframes));
// 	fprintf(fpaux,"\t%g", lasym);
// 	fprintf(fpaux,"\t%g", rad2deg*lF.alpha);
// 	fprintf(fpaux,"\t%g", rad2deg*lF.beta);
// 	fprintf(fpaux,"\t%g", rad2deg*lF.gamma);
//       }
      fprintf(fpdav,"\n");
      fprintf(fpaux,"\n");
    }
  }
      
  //  cout << "Distance change warnings per cell " << nwarnings << '\n';
  return 0;
}
