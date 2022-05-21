
#include "BaseList.h"
#include "geometry.h"
#include "basespin_system.h"
#include "ListList.h"
#include "NMR.h"

using namespace libcmatrix;

#include "geometry_utils.h"

struct default_spin_t {
  const char* elementname;
  size_t massnumber;
};

#define NMRSIM_DEFAULT_SPINS 57
const default_spin_t default_spin_table[NMRSIM_DEFAULT_SPINS] = {
  {"H",1},
             {"C", 13},
             {"F", 19},
             {"He", 3},
             {"Be", 9},
             {"B", 11},
             {"O", 17},
             {"Ne", 21},
             {"Na", 23},
             {"Mg", 25},
             {"Al", 27},
             {"Si", 29},
             {"P", 31},
             {"S", 33},
             {"K", 39},
             {"Ca", 43},
             {"Sc", 45},
             {"V", 51},
             {"Cr", 53},
             {"Mn", 55},
             {"Fe", 57},
             {"Co", 59},
             {"Ni", 61},
             {"Zn", 67},
             {"Ga", 71},
             {"Ge", 73},
             {"As", 75},
             {"Se", 77},
             {"Br", 81},
             {"Kr", 83},
             {"Rb", 87},
             {"Sr", 87},
             {"Y", 89},
             {"Zr", 91},
             {"Nb", 93},
             {"Mo", 95},
             {"Tc", 99},
             {"Rh", 103},
             {"Pd", 105},
             {"Cd", 113},
             {"In", 115},
             {"Sb", 121},
             {"I", 127},
             {"Sn", 119},
             {"Te", 125},
             {"Cs", 133},
             {"Ba", 137},
             {"La", 139},
             {"Ta", 181},
             {"W", 183},
             {"Re", 187},
             {"Ir", 193},
             {"Pt", 195},
             {"Au", 197},
             {"Tl", 205},
             {"Pb", 207},
             {"Bi", 209}
};


void fill(double d[3], const vector3& a)
{
  d[0]=a.x;
  d[1]=a.y;
  d[2]=a.z;
}

void transformation_matrix(rmatrix3& dest, const vector3& x, const vector3& y, const vector3& z)
{
  fill(dest[0],x);
  fill(dest[1],y);
  fill(dest[2],z);
}

void write_PDB_comment(FILE* fp, const char* str)
{
  fprintf(fp,"REMARK  99 %s\n",str);
}

void write_PDB_atom(FILE* fp, size_t serial, const char totlabel[8], const vector3& coords, size_t seq)
{
  fprintf(fp,"HETATM%5" LCM_PRI_SIZE_T_MODIFIER "u %-4s UNK 0 %3" LCM_PRI_SIZE_T_MODIFIER "u    %8.3f%8.3f%8.3f  1.00  0.00\n",serial,totlabel,seq,coords.x,coords.y,coords.z);
}

bool testsame(const Matrix<double>& An, const Matrix<double>& A, const char* head, int verbose)
{
  const double normdiff=norm(An-A);
  const bool fail=(normdiff>1e-6);

  if (verbose || fail)
    std::cout << head << '\n' << An << '\n';
  if (fail)
    std::cout << "FAILED (" << normdiff << ")\n";
  else
    std::cout << "SUCCESS\n";
  return fail;
}

size_t default_nucleus(const char* name)
{
  char buf[10];
  for (size_t i=0;i<NMRSIM_DEFAULT_SPINS;i++) {
    const default_spin_t& curspin(default_spin_table[i]);
    if (strcmp(name,curspin.elementname)==0) {
      sprintf(buf,"%zu%s",curspin.massnumber,name);
      return labeltonuc(buf);
    }
  }
  return NULL_NUCLEUS;
}
  
//check that subsel fits exactly within sel
bool checkoverlap(const BaseList<size_t>& sel, const BaseList<size_t>& subsel)
{
  for (size_t i=subsel.size();i--;) {
    if (!ismember(sel,subsel(i))) {
      std::cerr << "Methyl group " << subsel << " not entirely within selection " << sel << ":  results will be bogus\n";
      return false;
    }
  }
  return true;
}

bool BaseDipolarAverager::check_nooverlap(const BaseList<size_t>& second) const
{
  for (size_t j=second.size();j--;) {
    const int which=which_set(second(j));
    if (which>=0) {
      if (!checkoverlap(second,groups_(which)))
	return false;
    }
  }
  return true;
}

int BaseDipolarAverager::validate_set(const BaseList<size_t>& which) const
{
  if (which.size()==1)
    return -1;
  const int whichset=which_set_(which.front());
  const size_t n=which.size();
  if ((whichset>=0) && (groups_(whichset).size()==n)) {
    bool ok=true;
    for (size_t i=1;i<n;i++) {
      if (which_set_(which(i))!=whichset) {
	ok=false;
	break;
      }
    }
    if (ok)
      return whichset;
  }
  throw InvalidParameter("DipolarAverager: atom group does not match an averaging set");
}
    
BaseDipolarAverager::BaseDipolarAverager(const BaseList<vector3>& positionsv, double gammaval, const ListList<size_t>& groupsv, int verbosev)
  : natoms_(positionsv.size()), 
    positions_(positionsv),
    groups_(groupsv),
    verbose_(verbosev)
{
  gammavals_.create(natoms_,gammaval);
  initialise();
}

BaseDipolarAverager::BaseDipolarAverager(const BaseList<vector3>& positionsv, const BaseList<double>& gammavalsv, const ListList<size_t>& groupsv, int verbosev)
  : natoms_(positionsv.size()), 
    positions_(positionsv),
    groups_(groupsv),
    gammavals_(gammavalsv),
    verbose_(verbosev)
{
  initialise();
}

BaseDipolarAverager::BaseDipolarAverager(const BaseList<vector3>& positionsv, double gammaval, const BaseList< Matrix<size_t> >& groupsv, int verbosev)
  : natoms_(positionsv.size()), 
    positions_(positionsv),
    verbose_(verbosev)
{
  makegroups(groupsv);
  gammavals_.create(natoms_,gammaval);
  initialise();
}

BaseDipolarAverager::BaseDipolarAverager(const BaseList<vector3>& positionsv, const BaseList<double>& gammavalsv, const BaseList< Matrix<size_t> >& groupsv, int verbosev)
  : natoms_(positionsv.size()), 
    positions_(positionsv),
    gammavals_(gammavalsv),
    verbose_(verbosev)
{
  if (gammavals_.size()!=natoms_)
    throw Mismatch("DipolarAverager:: number of gamma values doesn't match number of atoms",gammavals_.size(),natoms_);

  makegroups(groupsv);
  initialise();
}

void BaseDipolarAverager::makegroups(const BaseList< Matrix<size_t> >& groupsv)
{
  size_t items=0;
  size_t blks=0;
  for (size_t i=groupsv.size();i--;) {
    const Matrix<size_t>& curm(groupsv(i));
    blks+=curm.rows();
    items+=curm.size();
  }
  groups_.reserve(blks,items);
  for (size_t i=0;i<groupsv.size();i++) {
    const Matrix<size_t>& curm(groupsv(i));
    for (size_t r=0;r<curm.rows();r++)
      groups_.push_back(curm.row(r));
  }
}

void BaseDipolarAverager::initialise()
{
  which_set_.create(natoms_,-1);
  for (size_t ng=0;ng<groups_.size();ng++) {
    const BaseList<size_t> curgroup(groups_(ng));
    if (curgroup.size()<2)
      throw InvalidParameter("DipolarAverager: group of atoms with <2 members");
    for (size_t na=curgroup.size();na--;) {
      const size_t curi=curgroup(na);
      validate_index(curi);
      if (which_set_(curi)!=-1)
	throw InvalidParameter("DipolarAverager: groups of atoms overlap");
      which_set_(curi)=ng;
    }
  }
}

std::pair<int,int> BaseDipolarAverager::validate_sets(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const
{
  const int seti=validate_set(whichi);
  const int setj=validate_set(whichj);
  if (seti>=0) {
    if (seti==setj)
      throw Failed("DipolarAverager: not appropriate for couplings within the same averaging set");
  }
  else {
    if (setj<0)
      throw Failed("DipolarAverager: not appropriate for couplings between non-averaged nuclei");
  }
  return std::pair<int,int>(seti,setj);
}
  

void TensorDipolarAverager::initialise_tensor()
{
  const size_t nsets=groups_.size();
  //  avposition.create(nsets);
  perpunitvector_.create(nsets);

  for (size_t i=0;i<nsets;i++) {
    const BaseList<size_t> curset(groups_(i));
    const size_t nsites=curset.size();
    if (nsites!=3)
      throw Failed("TensorDipolarAverager: only implemented for 3 site exchange");

    const vector3 planevec1(positions_(curset(size_t(0))) - positions_(curset(size_t(1))) );
    const vector3 planevec2(positions_(curset(size_t(1))) - positions_(curset(size_t(2))) );
    vector3 perp(cross(planevec1,planevec2));
    perp/=perp.length();
      
//     vector3 meanpos(0,0,0);
//     for (size_t c=nsites;c--;)
//       meanpos+=positions_(curset(c));

//     meanpos/=nsites;

    //    avposition(i)=meanpos;
    const spherical sphco(perp);
    perpunitvector_(i)=sphco;
    if (verbose_)
      //      std::cout << "For averaging group " << (i+1) << "  Mean position: " << meanpos << "  Perp unit vector (spherical polar coords): " << sphco << '\n';
      std::cout << "For averaging group " << (i+1) << ": Perp unit vector (spherical polar coords): " << sphco << '\n';
  }
}

void TensorDipolarAverager::d_and_Euler(double& d, double& asym, Euler& euler, const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const
{
  const std::pair<int,int> setij=validate_sets(whichi,whichj);
  const int seti=setij.first;
  const int setj=setij.second;
  const bool inavseti=(seti>=0);
  const bool inavsetj=(setj>=0);
  //vector3 locpositioni,locpositionj;

  size_t source;
  BaseList<size_t> curset;
  if (inavseti) { //!< NB doesn't treat the case of both i and j being averaged (gives dubious results)
    //    locpositioni=avposition(size_t(seti));
    source=whichj.front();
    curset.create(groups_(size_t(seti)));
  }
  else {
    if (!inavsetj)
      throw InternalError("notinavsetj");
    // locpositionj=avposition(size_t(setj));
    source=whichi.front();
    curset.create(groups_(size_t(setj)));
  }
  Matrix<double> Dav_MF,Dtmp;
  const size_t n=curset.size();
  const vector3& sourcepos(positions_(source));
  const double gammasource=unsafe_gamma(source);
  for (size_t k=n;k--;) {
    const size_t curind=curset(k);
    const vector3 ivec(positions_(curind) - sourcepos);	      
    const spherical sphco(ivec);
    const double d=dipolar_coupling(gammasource,unsafe_gamma(curind),sphco.r*1e-10);
    const space_T D_PAS(spatial_tensor(d));
    const Euler leuler(vector_to_Euler(sphco));
    const space_T D_MF(rotate(D_PAS,leuler));
    if (verbose_>1)
      std::cout << "Built tensor with coupling " << d << " Hz, orientation " << leuler << '\n' << D_MF << '\n';

    tensor_to_matrix(Dtmp,D_MF);
    Dav_MF+=Dtmp;
  }
  Dav_MF*=1.0/n;
  if (verbose_>1)
    std::cout << "Averaged tensor\n" << Dav_MF << '\n';
  
  double iso;
  cartesian_to_PAS_symmetric(iso,d,asym,euler,Dav_MF);
  static const double sqrt2thrd=sqrt(2.0/3.0);
  //	    iso/=-sqrt3rd;
  d/=sqrt2thrd; //!< should still be needed  
  if (verbose_>1)
    std::cout << "Final d=" << d << " Hz, orientation: " << euler << '\n';
}

double TensorDipolarAverager::dss(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const
{
  Euler junk;
  double d,asym;
  d_and_Euler(d,asym,junk,whichi,whichj);
  return d*d; 
}

