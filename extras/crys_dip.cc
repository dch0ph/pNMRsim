//Using PDB-file this program calculates root-sum-squared dipolar coupling between specific atoms

#include <fstream>
#include <string.h>
#include <fnmatch.h>
#include "ttyio.h"
#include "simpsonio.h"
#include "geometry.h"
#include "powder.h"
#include "wigner.h"
#include "cmatrix.h"
#include "List.h"
#include "cmatrix_utils.h"
#include <set>
#include <map>
#include "Histogram.h"
#include "matlabio.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

#include "./geometry_utils.h"

class PowderDipolarAverager : public BaseDipolarAverager 
{
public:
  PowderDipolarAverager(const BaseList<vector3>& positionsv, double gammaval, const ListList<size_t>& groupsv, size_t nzcwv, int verbosev)
    : BaseDipolarAverager(positionsv, gammaval, groupsv, verbosev),
      nzcw_(nzcwv)
  {}

  PowderDipolarAverager(const BaseList<vector3>& positionsv, const BaseList<double>& gammavals, const ListList<size_t>& groupsv, size_t nzcwv, int verbosev)
    : BaseDipolarAverager(positionsv, gammavals, groupsv, verbosev),
      nzcw_(nzcwv)
  {}

  double dss(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const;
  void d_and_Euler(double&, double&, Euler&, const BaseList<size_t>&, const BaseList<size_t>&) const { throw Failed("d_and_Euler not possible for PowderDipolarAverager"); }

private:
  size_t nzcw_;
};

double PowderDipolarAverager::dss(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const
{
  validate_sets(whichi,whichj);
  const size_t N=whichi.size()*whichj.size();

  PlanarZCW powdmeth(nzcw_);
  double scale;
  Euler powder;
  rmatrix3 rotmat;
  double ldss=0.0;
  
  while (powdmeth.next(powder,scale)) {
    double meand=0.0;
    rotation_matrix(rotmat,powder);
    for (size_t i=whichi.size();i--;) {
      const size_t seli=whichi(i);
      const double gammai=unsafe_gamma(seli);
      for (size_t j=whichj.size();j--;) {
	const size_t selj=whichj(j);
	const vector3 diff=positions_(seli)-positions_(selj);
	const vector3 rotvec(rotate(diff,rotmat));
	const spherical sphco(rotvec);
	if (sphco.r==0.0)
	  throw Failed("spins are co-incident");
	const double gammaj=unsafe_gamma(selj);
	const double d=dipolar_coupling(gammai,gammaj,sphco.r*1e-10);
	const double dweight=d*legendre(2,sphco.theta);
	if (verbose_>2)
	  std::cout << "P2-scaled d between spins " << (seli+1) << " and " << (selj+1) <<  ": " << dweight << " Hz\n";
	meand+=dweight;
      }
    }
    const double dav=meand/N;
    if (verbose_>1)
      std::cout << "Mean d at orientation " << powder << ": " << dav << " Hz\n";
    ldss+=scale*dav*dav;
  }  
  if (verbose_>1)
    std::cout << "Mean drss (with P2-scaling): " << std::sqrt(ldss) << " Hz\n";
  static const double scalef=20.0/9;
  return ldss*scalef;  
}

double standard_short=1.75; //!< normalise short distances to this (unless zero)

bool knowcastep=false; //!< if true, know that PDB was generated from labelled CASTEP
double global_cutoff=0.0; //!< minimum dipolar coupling to consider

int verbose=1;
bool asymmerge=false;

const size_t MAX_RECORD=100;
const size_t MAX_LINE=100;
const size_t MAX_TYPE_LEN=16;
//bool motionaveraging=false; //!< by default no motional averaging

enum av_t { AV_NONE=0, AV_POWDER, AV_TENSOR };
av_t averaging_method;
enum sort_t { SORT_NONE=0, SORT_ASU, SORT_UNITCELL };
const sort_t def_sorttype=SORT_UNITCELL;
enum output_t { OUT_NORMAL=0, OUT_TABLES, OUT_MAP, OUT_HIST };
output_t outputtype=OUT_NORMAL;

inline bool isorientationmap() { return (outputtype==OUT_MAP) || (outputtype==OUT_HIST); }

size_t nzcw=0;
size_t nbeta=1;
size_t nalpha=1;
char outname[256]="";

// Equal steps over alpha and beta
class Hemisphere : public PowderMethod {
public:
  Hemisphere(size_t asteps, size_t bsteps);
  void index_to_orientation(Euler&, double&, size_t) const;
  PowderMethod* clone() const { return new Hemisphere(*this); }

private:
  size_t asteps;
  double astep,bstep;

  double getbeta(size_t whichb) const { return whichb*bstep; }
};

  Hemisphere::Hemisphere(size_t _asteps,size_t _bsteps) : PowderMethod(_asteps*_bsteps,2), asteps(_asteps)
{
  if (asteps<2 || _bsteps<2)
    throw InvalidParameter("Hemisphere()");

  astep=(2*M_PI)/asteps;
  bstep=(0.5*M_PI)/(_bsteps-1);

  double total=0;
  for (size_t i=_bsteps;i--;)
    total+=sin(getbeta(i));
  normal=1.0/(asteps*total);
}

void Hemisphere::index_to_orientation(Euler& dest,double &w,size_t which) const
{
  dest.alpha=astep*(which % asteps);
  dest.beta=getbeta(which / asteps);
  w=normal*sin(dest.beta);
}

inline char* chew_null(char* start) {
  char* entry=start;
  while (*entry=='\0')
    entry++;
  return entry;
}

inline void push_ifnonnull(List<size_t>& atomlist, char* start)
{
  const size_t atom=atoi(chew_null(start)); //atom.serial number	
  if (atom)
    atomlist.push_back(atom);
}

//! NB sign of output has no meaning (unlike strcmp)
static int strrcmp(const char* a, const char* b)
{
  const int la=strlen(a);
  const int lb=strlen(b);
  return (lb>la) ? 1 : strcmp(a+la-lb,b);
}

//! slow search
size_t searchnuc(const char* name)
{
  size_t nuc=NULL_NUCLEUS;
  for (size_t i=maximumnuc()-1;i!=NULL_NUCLEUS;i--) {
    if (strrcmp(nuctolabel(i),name)==0) {
      if (nuc) {
	std::cerr << "Ambiguous nucleus name: " << name << '\n';
	nuc=NULL_NUCLEUS; //!< ambiguous label
	if (strcmp(name,"H")==0)
	  nuc=H_NUCLEUS;
	else if (strcmp(name,"C")==0)
	  nuc=labeltonuc("13C");
	else if (strcmp(name,"N")==0)
	  nuc=labeltonuc("15N");
	else if (strcmp(name,"F")==0)
	  nuc=labeltonuc("19F");
	if (nuc!=NULL_NUCLEUS)
	  std::cerr << "Defaulting to " << nuctolabel(nuc) << '\n';
	return nuc;
      }
      nuc=i;
    }
  }
  return nuc;	
}

typedef std::map<std::string,size_t> map_t;
static map_t nucmap;

//! cached search
size_t castep2nuc(const char* name)
{
  const std::string strname(name);
  const map_t::iterator iter(nucmap.find(strname));
  return (iter==nucmap.end())
    ? (nucmap[strname]=searchnuc(name))
    : iter->second;
}

//! determine isotope from label using first alpha characters
size_t get_isotope(const char* type)
{
  char buf[8]; 

  size_t i=0; //copy atom type (alphabet chars)
  for (;i<sizeof(buf)-1;i++) {
    if (isalpha(type[i]))
      buf[i]=type[i];
    else
      break;
  }
  buf[i]='\0'; //terminate
  return castep2nuc(buf);
}

template<typename T> bool ismember(const BaseList<T>& vec, const T& v)
{
  for (size_t i=vec.size();i--;) {
    if (vec[i]==v)
      return true;
  }
  return false;
}


ostream& operator<< (ostream& ostr, const set_t& a)
{
  const set_t::const_iterator end(a.end());
  set_t::const_iterator iter(a.begin());
  bool first=true;
  while (iter!=end) {
    if (first)
      first=false;
    else
      ostr << ',';
    ostr << *iter;
    ++iter;
  }
  return ostr;
}

//! contains information about an individual atom
struct mol_data {
  vector3 position; //!< coordinates
  char type[MAX_TYPE_LEN]; //!< atom type, e.g. H or C2
  int snumber; //!< serial number (usually starting from 1)
  size_t aunit; //!< which asymmetric unit
  size_t index; //!< index of element type (if labels not used)
  void rotate(const rmatrix3& R3, const rmatrix&) { position=libcmatrix::rotate(position,R3); }
};

std::ostream& operator<< (std::ostream& ostr, const mol_data& a)
{
  ostr << a.type;
  if (a.aunit)
    ostr << ':' << a.aunit;
  return ostr << " (" << a.snumber << "): " << a.position;
}

template<class T> void rearrange(List<T>& dest, const List<T>& source, const BaseList<size_t>& sizes)
{
  const size_t totsizes=sum(sizes);
  const size_t nels=source.size();
  if (nels % totsizes)
    throw Failed("rearrange: number of elements doesn't divide");
  dest.create(nels);
  const size_t cells=nels / totsizes;
  size_t sptr=0;
  size_t offset=0;
  for (size_t div=0;div<sizes.size();div++) {
    size_t destbase=offset;
    const size_t n=sizes(div);
    for (size_t cell=0;cell<cells;cell++) {
      for (size_t i=0;i<n;i++)
	dest(destbase+i)=source(sptr++);
      destbase+=totsizes;
    }
    offset+=n;
  }
}
  
class ZMolStructure {
public:
  ZMolStructure() : cellsize_(0), ncells_(0) {}

  void readPDB(const char*, sort_t); //!< read PDB-file with given filename analysing HETATM/ATOM and CONECT sections

  size_t findSerialNumber(size_t) const; //!< return internal index of atom with supplied serial number

  List<size_t> findAtomType(const char*, bool incell0 =false) const; //!< return internal indexes of atoms, which type is matched with argument

  void addAtomType(List<size_t>&, const char*, bool incell0 =false) const; //!< add internal indexes of atoms, which type is matched with argument to supplied list

  size_t whichMolecule(size_t m) const { //!< return internal index of molecule, whichs contain atom with given internal index
    if (hasconnections())
      throw Failed("whichMolecule impossible without connection data");
    return whichmolecule_(m);
  }

  bool samemolecule(size_t m, size_t n) const { return (whichMolecule(m)==whichMolecule(n)); }

  void print_atom(size_t, std::ostream&) const;

  bool hasconnections() const { return !(connections.empty()); }
  const mol_data& operator()(size_t i) const { return data(i); }
  mol_data& operator()(size_t i) { return data(i); }
  size_t size() const { return data.size(); }
  void print(std::ostream& =std::cout) const;
  void dumpunitcell(std::ostream& =std::cout) const;
  size_t cell_to_index(size_t ncell, size_t ind) const;

  void write_PDB(const char*) const;
  void set_zprime(size_t, bool uselabels);

private:
  List<mol_data> data; //!< mol_data for all atoms
  List< List<size_t> > connections; //!< serial numbers, grouped by molecules
  List<size_t> whichmolecule_;
  List<size_t> serial_to_internal_;
  vector3 centre_;
  size_t cellsize_,ncells_; //!< note that cells refers to ASU
  sort_t sorttype_;
  List<size_t> cell0contents_; //!< list of atoms in "cell0"

  void readPDB_normal(ifstream&); //!< normal file read (sorted by unit cell)
  void readPDB_resort(ifstream&); //!< resort file read
  bool readline(ifstream& file, char* buffer, size_t nbuf); //!< read line (true if successful)
  bool cleanline(mol_data&, char* buffer); //!< clean input line (false if empty)

  static bool partcmp(const char*, const char*); //!< compare how first atom type from input matches the second
  void addConnections(const BaseList<size_t>&); //!< put entries of CONECT section to the connections superlist (see above)
  void tidy_atoms();
  void tidy_connections();
  void clear_connections() { connections.clear(); }
  void clear() { data.clear(); clear_connections(); }
  sort_t sorttype() const { return sorttype_; }
};

size_t ZMolStructure::findSerialNumber(size_t n) const {
  if ((n<1) || (n>=serial_to_internal_.size())) {
    std::cerr << "Serial number is outside of range: " << n << '\n';
    exit(2);
  }
  const size_t ind=serial_to_internal_(n);
  if (verbose>1)
    std::cout << "Mapped PDB serial number " << n << " to internal index " << ind << '\n';
  return ind;
}

size_t ZMolStructure::cell_to_index(size_t ncell, size_t ind) const
{
  switch (sorttype_) {
  case SORT_ASU:
    return ncell*cellsize_+ind;
  case SORT_UNITCELL:
    return ind*ncells_+ncell;
  default:
    std::cerr << "Can't find atom by unit cell for unsorted type\n";
    exit(2);
  }
}

void ZMolStructure::set_zprime(size_t Z, bool uselabels)
{
  if (Z<1)
    throw InvalidParameter("ZMolStructure::set_zprime");
  if (ncells_ % Z) {
    std::cerr << "Z/Z' (" << Z << ") does not divide into number of cells (" << ncells_ << ")\n";
    exit(2);
  }

  if (!uselabels) { //!< need to number types
    if (sorttype_!=SORT_ASU) {
      std::cerr << "set_zprime only implemented for sorted unit cells (in the absence of :labels)\n";
      exit(2);
    }  

    typedef std::map<char,size_t> map_t;
    map_t counts;
    const size_t tcellsize=cellsize_*Z;

    for (size_t i=0;i<tcellsize;i++) {      
      mol_data& curat(data(i));
      const char el=curat.type[0];
      const map_t::iterator iter(counts.find(el));
      size_t count=(iter==counts.end()) ? 1 : (iter->second)+1;
      curat.index=count;
      curat.aunit=0;
      counts[el]=count;
    }
    for (size_t i=tcellsize;i<size();i++) {
      const size_t redi=i % tcellsize;
      data(i).index=data(redi).index;
      data(i).aunit=0;
    }
  }
  else {
    if (data.size() % cellsize_) 
      throw InternalError("Number of atoms not multiple of ASU size");

    if (sorttype_!=SORT_ASU) {
      const size_t reps = data.size() / cellsize_;
      for (size_t basei=0;basei<data.size();basei+=reps) {
	for (size_t i=0;i<reps;i++)
	  data(basei+i).aunit=1+(i % Z);
      }
    }
    else {
      for (size_t i=data.size();i--;) {    
	const size_t ucell=i / cellsize_;
	data(i).aunit=1+(ucell % Z);
      }
    }
    if (cellsize_>1) {
      for (size_t i=Z;i--;) {
	for (size_t j=i;j--;) {
	  const size_t startind=cell_to_index(i,0);
	  const size_t otherind=cell_to_index(j,0);
	  const vector3 onevec=data(startind).position-data(otherind).position;
	  const size_t startind1=cell_to_index(i,1);
	  const size_t otherind1=cell_to_index(j,1);
	  const vector3 secondvec=data(startind1).position-data(otherind1).position;
	  const vector3 diff(secondvec-onevec);
	  const double difflen=diff.length();
	  if (verbose)
	    std::cout << "ASU " << (i+1) << ", ASU " << (j+1) << " differential distance: " << difflen << " A\n";
	  if (difflen<1e-6) {
	    std::cerr << "Unit cells that should be different asymmetric units seem to be related by translation.\n";
	    std::cerr << "Vector between atom " << (startind+1) << " in ASU " << (i+1) << " and " << (otherind+1) << " in ASU " << (j+1) << ": " << onevec << " A\n";
	    std::cerr << "Vector between atom " << (startind+2) << " in ASU " << (i+1) << " and " << (otherind+2) << " in ASU " << (j+1) << ": " << secondvec << " A\n";
	    exit(2);
	  }
	}
      }
      if (verbose)      
	std::cout << "Verified that different asymmetric units are not related by translation\n";
    }  
  }
}
    
void ZMolStructure::print_atom(size_t i,std::ostream& ostr =std::cout) const
{
  const mol_data& a((*this)(i));
  ostr << a.type;
  if (a.aunit)
    ostr << ':' << a.aunit;
  ostr << '('<< a.snumber;
  if (cellsize_) {
    switch (sorttype_) {
    case SORT_ASU:
      ostr << ':' << ((i % cellsize_)+1) << ", cell " << ((i / cellsize_)+1);
      break;
    case SORT_UNITCELL:
      ostr << ':' << ((i / ncells_)+1) << ", cell " << ((i % ncells_)+1);
      break;
    }
  }
  ostr << ')';
}

void ZMolStructure::print(std::ostream& ostr) const
{
  const bool showmarker=(cellsize_!=0) && (sorttype_==SORT_ASU);
  for (size_t i=0;i<size();i++) {
    if (showmarker && ((i % cellsize_)==0))      
      ostr << "----- cell (ASU) " << ((i/cellsize_)+1) << " -----\n";
    ostr << i << ": " << (*this)(i) << '\n';
  }
}

void ZMolStructure::dumpunitcell(std::ostream& ostr) const
{
  for (size_t i=0;i<cell0contents_.size();i++)
    ostr << data(cell0contents_(i)) << '\n';      
}  

void ZMolStructure::tidy_atoms()
{
  size_t maxs=0;
  for (size_t i=data.size();i--;) {
    if (data(i).snumber>maxs)
      maxs=data(i).snumber;
  }
  List<size_t>& mutable_index=const_cast<ZMolStructure*>(this)->serial_to_internal_;
  mutable_index.create(maxs+1,size_t(-1));
  vector3 total(0,0,0);
  for (size_t i=data.size();i--;)
    mutable_index(data(i).snumber)=i;

  List<vector3> centres(ncells_,vector3(0,0,0));
  const double scalef=1.0/cellsize_;
  for (size_t c=ncells_;c--;) {
    for (size_t i=cellsize_;i--;) {
      const size_t ind=cell_to_index(c,i);      
      const vector3& pos(data(ind).position);
      total+=pos;
      centres(c)+=pos;
    }
    centres(c)*=scalef;
  }
  centre_=total*(1.0/data.size());
  
  if (sorttype_!=SORT_NONE) {
    size_t cell0=0U;
    if (ncells_>1) {

      double mindist=1e30;
      for (size_t i=ncells_;i--;) {
	const double dist=(centres(i)-centre_).length();
	if (verbose>1)
	  std::cout << "Centre of cell " << (i+1) << " at " << centres(i) << " (distance=" << dist << " A)\n";
	if (dist<mindist) {
	  mindist=dist;
	  cell0=i;
	}
      }
      if (verbose) {
	std::cout << "Found " << ncells_ << " 'unit cells' of " << cellsize_ << " atoms\n";
	std::cout << "Cell " << (cell0+1) << " with centre at " << centres(cell0) << " is closest to structure centre at " << centre_ << " A (distance of " << mindist << " A)\n";
      }
    }
    if (sorttype_==SORT_ASU)
      cell0contents_=slice(cell0*cellsize_,cellsize_);
    else
      cell0contents_=slice(cell0,cellsize_,ncells_);

    if (verbose)
      std::cout << "Central cell indices: " << cell0contents_ << '\n';    
  }
  else { // not sorted by cell
    if (verbose)
      std::cout << "Finding atoms closest to origin at " << centre_ << " A\n";
    cell0contents_.create(cellsize_);
    //! find atoms closest to cell origin
    for (size_t typei=0;typei<cellsize_;typei++) { //loop over atom types
      double mindist=1e30;
      size_t mini=-1;
      for (size_t atomi=ncells_;atomi--;) {
	const size_t ind=typei*ncells_+atomi;
	const mol_data& curatom(data(ind));
	const vector3 diff(curatom.position-centre_);
	const double len=diff.norm();
	if (len<mindist) {
	  mindist=len;
	  mini=ind;
	}
      }
      cell0contents_(typei)=mini;
    }
  }
  if (verbose) {
    std::cout << "Central unit cell contents:\n";
    dumpunitcell();
  }
}

void ZMolStructure::tidy_connections()
{
  if (!hasconnections())
    return;
    
  List<size_t>& mutable_index=const_cast<ZMolStructure*>(this)->whichmolecule_;
  for (size_t m=connections.size();m--;) {
    const List<size_t>& curmol(connections(m));
    for (size_t n=curmol.size();n--;)
      mutable_index(serial_to_internal_(curmol(n)))=m;
  }
}

inline bool globcompare(const char* pattern, const char* str)
{
  switch (fnmatch(pattern,str,FNM_NOESCAPE)) {
  case 0:
    return true;
  case FNM_NOMATCH:
    return false;
  }
  perror("globcompare");
  exit(1);
}

bool ZMolStructure::partcmp(const char* ref, const char* second) //compare atom types from user's input and PDB-file (e.g. C and C2 are the same, F2 and F2A are the same etc.)
{
  const int n=strlen(ref)-1;
  
  if (ref[n]!='*')
    return (strcmp(ref,second)==0);

  char p[MAX_TYPE_LEN];
  strcpy(p,ref);
  p[n]='\0'; //get rid of trailing *

  for (size_t i=0; i<n; i++) {
    if (p[i]!=second[i])
      return false;
  }
  if (isdigit(second[n]) && isdigit(second[n-1])) // e.g. F2 and F21 are different
    return false;
  if (!isdigit(second[n]) && !isdigit(second[n-1]) && (strlen(second)!=n)) // e.g. C and Cu are different
    return false;
  return true;
}

List<size_t> ZMolStructure::findAtomType(const char* type, bool incell0) const
{
  List<size_t> result;
  addAtomType(result,type,incell0);
  return result;
}

void ZMolStructure::addAtomType(List<size_t>& result, const char* type, bool incell0) const
{
  if (incell0) {
    for (size_t i=0;i<cell0contents_.size();i++) {
      const size_t ind=cell0contents_(i);
      if (globcompare(type, data(ind).type))
	result.push_back(ind);
    }
  }
  else {
    for (size_t i=0;i<data.size();i++) {
      if (globcompare(type, data(i).type))
	result.push_back(i);
    }
  }
}

void ZMolStructure::addConnections(const BaseList<size_t>& newcon)
{
  for (size_t i=newcon.size();i--;) {  //looking for what molecule to add
    const size_t curatom=newcon(i);
    for (size_t m=connections.size();m--;) {
      List<size_t>& curmol(connections(m));
      if (ismember(curmol,curatom)) {
	for (size_t j=newcon.size();j--;) { //add all new connections
	  if (!ismember(curmol,newcon(j)))
	    curmol.push_back(newcon(j));
	}
	return;
      }
    }
  }
  //molecule is not found i.e. all connections are new
  connections.push_back(newcon);
}

void ZMolStructure::readPDB(const char* fname, sort_t sorttype) //parser for PDB-files (from standard guide)
//see www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
//the parser is not very clever and can fault if PDB-file slightly non-standard
{
  ifstream file(fname);
  if (!file.is_open())
    throw Failed("readPDB: file open failed");

  clear();
  sorttype_=sorttype;
  if (sorttype!=SORT_ASU)
    readPDB_resort(file);
  else
    readPDB_normal(file);
  file.close();

  if (hasconnections() && verbose) {
    std::cout << connections.size() << " molecules found\n";
    if (verbose>1)
      std::cout<< connections<< '\n';
  }
  if (verbose)
    std::cout << "Found " << data.size() << " atoms\n";
}

struct mol_data_type_order {
  bool operator()(const mol_data& a, const mol_data& b) {
    return (strcmp(a.type,b.type)<0);
  }
};

void ZMolStructure::readPDB_resort(ifstream& file)
{
  char buffer[MAX_RECORD];
  mol_data single_entry;
  List<size_t> atomlist;
  while (readline(file,buffer,sizeof(buffer))) {
    if (cleanline(single_entry,buffer))
      data.push_back(single_entry);
  }
  stable_sort_ip(data,mol_data_type_order()); //!< sort by label
  size_t ncells=0; //!< flag number of cells
  size_t lastreset=0;
  std::cout << "Found " << data.size() << " atoms\n";
  if (verbose>1) {
    FILE* fp=fopen("atomdump","w");
    for (size_t i=0;i<data.size();i++)
      fprintf(fp,"%i\t%s\n",data(i).snumber,data(i).type);
    fclose(fp);
  }
  for (size_t i=1;i<data.size();i++) {
    const char* curtype=data(i).type;
    const bool match=(strcmp(data(lastreset).type,curtype)==0);
    if (ncells==0) {
      if (!match) {
	ncells=i;
	if (verbose)
	  std::cout << "Assuming " << ncells <<  " asymmetric unit(s)\n";
	if (data.size() % ncells) {
	  std::cerr << "Number of atoms (" << data.size() << ") is not multiple of number of asymmetric units (" << ncells << ")\n";
	  exit(1);
	}
	lastreset=i;
      }
    }
    else {
      if ((i % ncells)==0) { //at boundary?
	if (match) {
	  std::cerr << "Expecting new atom type at index " << (i+1) << " (serial number " << data(i).snumber << ") but found " << curtype << " again\n";
	  exit(1);
	}
	lastreset=i;
	if (verbose>1)
	  std::cerr << "New atom type: " << data(lastreset).type << '\n';
      }
      else {
	if (!match) {
	  std::cerr << "Expecting repeated atom type " << data(lastreset).type << " at index " << (i+1) << " (serial number " << data(i).snumber << ") but found " << curtype << "\n";
	  exit(1);
	}
      }
    }
  }
  if (ncells==0) {
    ncells_=data.size();
    cellsize_=1;
  }
  else {
    ncells_=ncells;
    cellsize_=data.size() / ncells;
  }
  tidy_atoms();
}

bool ZMolStructure::readline(ifstream& file, char* buffer, size_t nbuf)
{
  if (file.eof())
    return false;
  file.getline(buffer,nbuf);
  return (strcmp(buffer,"ENDM")!=0);
}

bool ZMolStructure::cleanline(mol_data& single_entry, char* buffer)
{
  const size_t len=strlen(buffer);
  if (len<6)
    return false;
  for (size_t i=len;i--;) {
    char& curchr(buffer[i]);
    if (curchr==' ') 
      curchr='\0'; // Substitute all spaces to eol symbols
    else
      curchr=toupper(curchr); //enforce all upper case
  }
  if (strncmp(buffer,"HETATM",6) && strcmp(buffer,"ATOM"))
    return false; //!< no atom data
  if (buffer[13]=='Q')
    return false; //ignore pseudo-atoms

  single_entry.snumber=atoi(chew_null(buffer+6));  //enter (index) (PDB serial number);
  strcpy(single_entry.type,chew_null(buffer+12)); //enter (atom name)
  char* colonptr=strchr(single_entry.type,':');
  if (colonptr) {// strike out :
    knowcastep=true;
    do {
      *colonptr=colonptr[1];
    } while (*(++colonptr));
  }
    
  single_entry.aunit=0; //PDB file has no asymmetric unit info
  single_entry.position.x=atof(chew_null(buffer+30));
  single_entry.position.y=atof(chew_null(buffer+38));
  single_entry.position.z=atof(chew_null(buffer+46));
  return (single_entry.position.x<=9000) || (single_entry.position.y<=9000) || (single_entry.position.z<=9000); //ignore XPLOR pseudo atoms
}

void ZMolStructure::readPDB_normal(ifstream& file)
{
  mol_data single_entry;
  List<size_t> atomlist;
  List<size_t> molsizes;
  char buffer[MAX_RECORD];
  
  cellsize_=0;
  size_t lastreset=0;
  while (readline(file,buffer,sizeof(buffer))) { //analyse each string
    if (cleanline(single_entry,buffer)) {      
      const size_t curoffset=data.size()-lastreset;
      if (cellsize_) {
	if (((curoffset % cellsize_)==0) && (strcmp(data(lastreset).type,single_entry.type)!=0)) { //mismatch from expected
	  if (verbose)
	    std::cout << "Additional molecule detected\n";
	  molsizes.push_back(cellsize_);
	  cellsize_=0; //reset cell size
	  lastreset+=curoffset;
	}
      }
      else {
	if (!data.empty() && (strcmp(data(lastreset).type,single_entry.type)==0)) {//starting to repeat entries?
	  cellsize_=curoffset;
	  if (verbose)
	    std::cout << "Repeated label (" << single_entry.type << ") used to determine unit cell size: " << cellsize_ << " atoms\n";
	}
      }
      
      if (verbose>1) {
	if (cellsize_ && ((curoffset% cellsize_)==0))
	  std::cout << "--- cell " << (1+(curoffset / cellsize_)) << " ---\n";
	std::cout << single_entry.snumber << ": " <<  single_entry.type << " at " << single_entry.position << '\n';
      }

      data.push_back(single_entry);
      if (!strcmp(buffer,"CONECT")) {
	atomlist.create(0); //reset
	push_ifnonnull(atomlist,buffer+6);
	push_ifnonnull(atomlist,buffer+11);
	push_ifnonnull(atomlist,buffer+16);
	push_ifnonnull(atomlist,buffer+21);
	push_ifnonnull(atomlist,buffer+26);
	addConnections(atomlist);
	if (verbose>2)
	  std::cout<< "Connected atoms: " << atomlist << '\n';
      }      
    }
  }
  if (cellsize_==0) {
    cellsize_=data.size()-lastreset;
    ncells_=1;
  }
  if (!molsizes.empty()) { //multiple molecules per cell
    molsizes.push_back(cellsize_);
    if (verbose)
      std::cout << "Molecule sizes: " << molsizes << '\n';
    //re-arrange
    if (hasconnections()) {
      if (verbose)
	std::cout << "Warning: connection information discarded\n";
      clear_connections(); //too much trouble to re-order
    }
    cellsize_=sum(molsizes);
    if (data.size() % cellsize_) {
      std::cerr << data.size() << " atoms does not divide into total cell size of " << cellsize_ << " atoms (atoms at special positions will cause problems!)\n";
      exit(1);
    }
    ncells_ = data.size() / cellsize_;
    List<mol_data> ndata;
    rearrange(ndata,data,molsizes);    
    data.swap(ndata);
    tidy_atoms();
  }    
  else {
    if (cellsize_==0) {
      std::cerr << "Failed to identify unit cell - file corrupt?\n";
      exit(1);
    }
    if (data.size() % cellsize_)
      std::cerr << "Warning: " << data.size() << " atoms does not divide into cell size of " << cellsize_ << " atoms\n";
    else
      ncells_ = data.size() / cellsize_;

    tidy_atoms();
    tidy_connections();
  }
}

size_t resolve_isotope(const ZMolStructure& cryst, size_t, const char* spec, bool safe =false)
{
  char* end;
  const long snum=strtol(spec,&end,10);
  if (*end=='\0') {
    const size_t ind=cryst.findSerialNumber(snum);
    return get_isotope(cryst(ind).type);
  }
  char scratch[20];
//   if (end-spec>sizeof(scratch)-2)
//     throw Failed("buffer overrun");
  strncpy(scratch,spec,end-spec);
  char* source=end;
  char* dest=scratch+(end-spec);
  while (isalpha(*source))
    *dest++=*source++;
  *dest++='\0';
  size_t nuc=default_nucleus(scratch);
  try {
    if (nuc==NULL_NUCLEUS)
      nuc=labeltonuc(scratch);
  } catch (...) {
    if (safe)
      return NULL_NUCLEUS;
    std::cerr << "Failed to identify nucleus from " << scratch << '\n';
    exit(2);
  }
  return nuc;
}

void spec_to_indices_(List<size_t>& dest, const ZMolStructure& cryst, const char* spec, bool incell0 =false)
{
  char* end;
  const long ind=strtol(spec,&end,10);
  if (*end=='\0')
    dest.push_back(cryst.findSerialNumber(ind));
  else {
    const size_t curlen=dest.size();
    cryst.addAtomType(dest,end,incell0);
    if (dest.size()==curlen)
      std::cerr << "Warning: found no atoms matching '" << spec << "'\n";
    //    std::cout << "After adding " << spec << ": " << dest << '\n';
  }
}

inline List<size_t> spec_to_indices(const ZMolStructure& cryst, char* cptr, bool incell0 =false)
{
  List<size_t> dest;
  for (;;) {
    char* tokptr=strtok(cptr,",");
    if (tokptr==NULL)
      break;
    cptr=NULL;
    spec_to_indices_(dest,cryst,tokptr,incell0);
  }
  return dest;
}

size_t spec_to_index(const ZMolStructure& cryst, const char* spec)
{
  char* end;
  const long ind=strtol(spec,&end,10);
  if (*end=='\0')
    return cryst.findSerialNumber(ind);
  List<size_t> tmp(cryst.findAtomType(end,true));
  if (tmp.size()!=1) {
    std::cerr << "Specification (" << spec << ") did not return a single atom: " << tmp << '\n';
    exit(2);
  }
  return tmp.front();
}

struct indexed_type_order {
  indexed_type_order(const ZMolStructure& crystv) : cryst(crystv) {}

  bool operator()(const size_t a, const size_t b) {
    return (strcmp(cryst(a).type,cryst(b).type)<0);
  }
  
  const ZMolStructure& cryst;
};

void buildlists(List< List<size_t> >& out, const ZMolStructure& cryst, char* spec)
{
  List<size_t> matchinginds(spec_to_indices(cryst,spec));
  if (matchinginds.empty()) {
    std::cerr << "Failed to find sites matching: " << spec << '\n';
    exit(1);
  }
  if (isdigit(spec[0])) {
    out.create(1U);    
    out.front()=matchinginds;
    return;
  }
  indexed_type_order indexer(cryst);
  if (verbose>1)
    std::cout << "Indices of matching atoms: " << matchinginds << '\n';
  sort_ip(matchinginds,indexer); //!< sorted by label
  if (verbose>1)
    std::cout << "Sorted by group: " << matchinginds << '\n';
  const char* lasttype;
  out.clear();
  List<size_t> tmp;
  bool donecheck=!asymmerge;
  size_t comparelen=8; //!< maximum label length
  for (size_t curind=0;curind<matchinginds.size();curind++) {
    const size_t atomi=matchinginds(curind);
    if ((lasttype!=NULL) && (strncmp(cryst(atomi).type,lasttype,comparelen)==0))
      tmp.push_back(atomi);
    else {
      if (lasttype) {
	if (verbose)
	  std::cout << "Found group of atoms of type " << lasttype << ": " << tmp << '\n'; 
	out.push_back(tmp); // store last list
      }
      lasttype=cryst(atomi).type;
      if (!donecheck) {	
	comparelen=strlen(lasttype)-1;	
	if (!isalpha(lasttype[comparelen]))
	  std::cerr << "Warning: merging asymmetric units and last character is not alphabetic: " << lasttype << '\n';
	donecheck=true;
      }
      tmp.create(1U,atomi);
    }
  }
  if (verbose)
    std::cout << "Found group of atoms of type " << lasttype << ": " << tmp << '\n'; 
  out.push_back(tmp);    
}

template<class T> void accumulate_intersection(T& dest, const T& source)
{
  T tmp;
  set_intersection(dest.begin(),dest.end(),source.begin(),source.end(),insert_iterator<T>(tmp,tmp.begin()));
  dest.swap(tmp);
}

void print_group(const ZMolStructure& cryst, const BaseList<size_t>& which, std::ostream& ostr =std::cout)
{
  if (which.size()==1)
    cryst.print_atom(which.front(),ostr);
  else
    ostr << "group " << which;
}

double average_dss(const BaseDipolarAverager& averager, const ZMolStructure& cryst, Matrix<bool>& included, BaseList<bool> atomincluded, const BaseList<size_t>& whichi, const BaseList<size_t>& whichj, rmatrix3& rotmat)
{
  if (included(whichi.front(),whichj.front()))
    return 0.0;
  atomincluded(whichi)=true;

  for (size_t i=whichi.size();i--;) {
    const size_t seli=whichi(i);
    for (size_t j=whichj.size();j--;) {
      const size_t selj=whichj(j);
      if (included(seli,selj))
	std::cerr << "Warning coupling " << seli << ',' << selj << " seems to be multiply included!\n";
      included(seli,selj)=included(selj,seli)=true;
    }
  }

  const size_t N=whichi.size()*whichj.size();
  double dss=0.0;

  if (isorientationmap()) {
    double meand=0.0;
    for (size_t i=whichi.size();i--;) {
      const size_t seli=whichi(i);
      const size_t gammai=averager.gamma(seli);
      for (size_t j=whichj.size();j--;) {
	const size_t selj=whichj(j);
	const size_t gammaj=averager.gamma(selj);
	const vector3 diff(cryst(seli).position-cryst(selj).position);
	const vector3 rotvec(rotate(diff,rotmat));
	const spherical sphco(rotvec);
	if (sphco.r==0.0)
	  throw Failed("spins are co-incident");
	double d=dipolar_coupling(gammai,gammaj,sphco.r*1e-10);
	d*=legendre(2,sphco.theta);
	meand+=d;
      }
    }
    dss=meand*meand/N;  
  }
  else {
    dss=averager.dss(whichi,whichj);
  }
  if (verbose) {
    cout << "Averaged d_rss between ";
    print_group(cryst,whichi);
    cout << " and ";
    print_group(cryst,whichj);
    cout << ": " << sqrt(dss) << " Hz\n";
  }
  return dss;
}

void dumptype(FILE* fp, const char* s)
{
  if (asymmerge) //strip asymmetric cell marker
    fwrite(s,1,strlen(s)-1,fp);
  else 
    fputs(s,fp);
}

//! output data table 
void annotated_table(const Matrix<double>& data, const char* sup, const ZMolStructure& cryst, const BaseList<size_t>& first, const BaseList< List<size_t> >& seconds, bool annot =true)
{
  char buffer[256];
  snprintf(buffer, sizeof(buffer), "%s_%s",outname,sup);
  if (annot) {
    FILE* fp=fopen(buffer,"w");
    for (size_t c=0;c<seconds.size();c++) {
      fprintf(fp,"\t");
      dumptype(fp,cryst(seconds(c).front()).type);//!< output name of first in each (should be all same)
    }
    fprintf(fp,"\n");
    for (size_t r=0;r<first.size();r++) {
      dumptype(fp,cryst(first(r)).type);
      for (size_t c=0;c<data.cols();c++)
	fprintf(fp,"\t%g",data(r,c));
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  else
    write_matrix(buffer,data);
}

void doextract(const char* outname, const char* source, const ZMolStructure& cryst, const char* first_type, const char* rsecond_type, int maxn, bool uselabels)
{
  char* second_type=strdup(rsecond_type);
  if (maxn<1) {
    std::cerr << "Bad number of spins\n";
    exit(2);
  }
  List<size_t> includeatoms;
  const size_t origind=spec_to_index(cryst,first_type);  
  const mol_data& sourceat(cryst(origind));
  List<size_t> secondlist;

  if (*second_type!='\0') {
    secondlist=spec_to_indices(cryst,second_type);
    if (secondlist.empty()) {
      std::cerr << "No atoms found of type: " << second_type << '\n';    
      exit(2);
    }
//     if ((minrad==0) && (std::find(secondlist.begin(),secondlist.end(),origind)==secondlist.end()))
//       secondlist.push_back(origind);
    if (verbose)
      std::cout << "Internal indices of selected atoms: " << secondlist << '\n';
  }
  size_t n=secondlist.size();
  if (n==0)
    n=cryst.size();
  List<double> distances(n);
  for (size_t i=0;i<n;i++) {
    const size_t curind= secondlist.empty() ? i : secondlist(i);    
    const mol_data& curat(cryst(curind));
    //    std::cout << curat << '\n';
    const vector3 ivec(curat.position-sourceat.position);
    distances(i)=ivec.length();
    //    if ((d>=minrad) && (d<=maxrad))
    //  includeatoms.push_back(curind);
  }
  const List<size_t> rawsortindex(indexed_sort(distances));
  List<size_t> rindex;
  if (secondlist.size())
    rindex=secondlist(rawsortindex);
  else
    rindex=rawsortindex;

  if (rindex.size()<maxn)
    maxn=rindex.size();
  //  if (verbose) {
  //  const List<size_t> index(rindex(range(maxn)));  
//     const List<double> sdistances(distances(index));
//     std::cout << "Sorted distances of selected atoms: " << sdistances << " A\n";  
//   }
//   if (includeatoms.empty()) {
//     std::cerr << "No selected atoms found!\n";
//     exit(2);
//   }

  FILE* fp=fopen(outname,"w");
  if (fp==NULL) {
    std::cerr << "Failed to open output file: " << outname << '\n';
    exit(2);
  }  
  char comline[MAX_LINE];
  snprintf(comline,sizeof(comline)-1,"Created by crys_dip from %s",source);
  write_PDB_comment(fp,comline);
  snprintf(comline,sizeof(comline)-1,"Starting from %s",first_type);
  write_PDB_comment(fp,comline);
  if (*second_type) {
    snprintf(comline,sizeof(comline)-1,"Selecting %s",rsecond_type);
    write_PDB_comment(fp,comline);
  }

  char label[8];

  for (size_t i=0;i<maxn;i++) {
    const mol_data& curat(cryst(rindex(i)));
    if (uselabels)
      strncpy(label,curat.type,sizeof(label));
    else
      snprintf(label,sizeof(label),"%c%" LCM_PRI_SIZE_T_MODIFIER "u",curat.type[0],curat.index);
    snprintf(comline,MAX_LINE-1,"Atom %" LCM_PRI_SIZE_T_MODIFIER "u (%s) is %.5g A from origin",i+1,label,distances(rawsortindex(i)));
    write_PDB_comment(fp,comline);
  }

  for (size_t i=0;i<maxn;i++) {
    const mol_data& curat(cryst(rindex(i)));
    if (uselabels)
      strncpy(label,curat.type,sizeof(label));
    else
      snprintf(label,sizeof(label),"%c%" LCM_PRI_SIZE_T_MODIFIER "u",curat.type[0],curat.index);
    const vector3 ivec(curat.position-sourceat.position);
    size_t seq=0;
    if (uselabels) {
      seq=curat.aunit;
      if (seq==0)
	seq=1;
    }
    write_PDB_atom(fp,i+1,label,ivec,seq);
  }
  fclose(fp);
}

//! routine that does all the work
void analyse_molecule(const ZMolStructure& cryst, char* first_type, char* second_type, FILE* outfile =NULL)
{  
  const size_t natoms=cryst.size();
  
  const List<size_t> Hlist(cryst.findAtomType("H*")); //!< list of all H atoms
  if (verbose>1)
    std::cout << "H atoms: " << Hlist << '\n';
  const size_t nH=Hlist.size();
  const double cutoff=2.0; //!< cutoff (A) for considering H's as 'neighbours'
  List<set_t> neighbours(nH);
  List<int> reverse(natoms,-1); //!< returns H index for given atom (-1 if not H)
  //!< find all pairs of H's within 2 A
  for (size_t i=nH;i--;) {
    reverse(Hlist(i))=i; //!< add to reverse H index
    //    neighbours(i).insert(Hlist(i)); //!< each H is neighbour to itself?
    const vector3 posi(cryst(Hlist(i)).position); //!< position of current H
    for (size_t j=i+1;j<nH;j++) { //!< loop over remaining H's in list
      const vector3 posj(cryst(Hlist(j)).position);
      const double r=(posi-posj).length(); //!< distance between current H and other
      if (r<cutoff) {
	neighbours(i).insert(Hlist(j)); //!< flag that j is a neighbour of spin i
	neighbours(j).insert(Hlist(i)); //!< flag that i is a neighbour of spin j
      }
    }
  }
  
  List<bool> assigned(natoms,false); //!< flags whether any (H) atom has been asigned to set
  ListList<size_t> methyllist;
  List<int> pairedatom(natoms,-1); //!< non zero if atom is paired with 1 one other
  //  List<int> whichmethyl(natoms,-1);
  for (size_t i=0;i<nH;i++) {
    const set_t& curlist(neighbours(i)); //!< list of neighbours of spin i
    if (assigned(Hlist(i)) || curlist.empty())
      continue; //already assigned or no neighbours
    if (verbose>1)
      cout << "Neighbours of H " << Hlist(i) << ": " << curlist << '\n';

    //! accummulate intersection of all lists - is this necessary?    
    //    set_t intersect(curlist);
    //const set_t::const_iterator end(curlist.end());
    //set_t::const_iterator iter(curlist.begin());
    //while (iter!=end) {
    //  const size_t j=*iter;
    //  ++iter;
    //  if (j!=i) //ignore self
    //	accumulate_intersection(intersect,neighbours(reverse(j)));
    //}

    const size_t curn=curlist.size();
    set_t::const_iterator liter(curlist.begin());
    const size_t i1=Hlist(i); //!< atom number rather than H number
    const size_t i2=*(liter++);
    switch (curn) {
    case 1: {
      pairedatom(i1)=i2; //!< create neighbour pair
      pairedatom(i2)=i1;
      assigned(i1)=true;
      assigned(i2)=true;
      if (verbose)
	std::cout << "Found pair involving " << cryst(i1) << " and " << cryst(i2) << '\n';
    }
      break;
    case 2: {
      const size_t i3=*(liter++);
      const size_t whichm=methyllist.size();
      methyllist.push_back(ScratchList<size_t>(i1,i2,i3));
      if (verbose)
	std::cout << "Found methyl group number " << (whichm+1) << " involving " << cryst(i1) << ", " << cryst(i2) << " and " << cryst(i3) << '\n';
      assigned(i1)=true;
      assigned(i2)=true;
      assigned(i3)=true;
      break;
    }
    default:
      std::cerr << "Unexpected number of neighbours: " << curn << '\n';
    }
  }
  
  List<vector3> positions(natoms);
  for (size_t i=natoms;i--;)
    positions(i)=cryst(i).position;

  if (nH && (verbose>1)) {
    std::cout << "Atom pairings: " << pairedatom << '\n';
    //    std::cout << "In methyl: " << whichmethyl << '\n';
  }

  const List<size_t> first(spec_to_indices(cryst,first_type,true));      
  const size_t nfirst=first.size(); //number of "from" atoms
  if (nfirst==0) {
    std::cerr << "No observe atoms selected!\n";
    return;
  }
  if (verbose)
    std::cout << "Internal indices of observed atoms: " << first << '\n'; 

  const double gam1=libcmatrix::gamma(resolve_isotope(cryst,first.front(),first_type));
  if (verbose>1)
    std::cout << "Gamma value for 'source' atom: " << gam1 << '\n';

  List<double> gammavals(natoms,0.0);
  gammavals(first)=gam1;

  List< List<size_t> > seconds;
  Matrix<double> drsstable,rmintable,dmaxtable;

  if (outputtype==OUT_TABLES) {
    buildlists(seconds,cryst,second_type);

    const size_t ngroups=seconds.size();
    drsstable.create(nfirst,ngroups);
    rmintable.create(nfirst,ngroups);
    dmaxtable.create(nfirst,ngroups);  
  }
  else {
    seconds.create(1U);
    seconds.front()=spec_to_indices(cryst,second_type,false);
  }

  double d2=0; //sum of d^2
  List<rmatrix> drssmap(nfirst);

  for (size_t groupi=0;groupi<seconds.size();groupi++) { //!< loop over master list
    const List<size_t>& second(seconds(groupi));
    const size_t nsecond=second.size(); //number of "to" atoms
    if (nsecond==0) {
      std::cerr << "No remote atoms selected in set " << (groupi+1) << "\n";
      return;
    }
    if (verbose)
      std::cout << "Internal indices of remote atoms in set " << (groupi+1) << ": " << second << '\n'; 

    //!< Don't actually check the sets are homonuclear
    const double gam2=libcmatrix::gamma(resolve_isotope(cryst,second.front(),second_type));

    //    if ((gam1!=gam2) && (averaging_method==AV_POWDER))
    //  throw Failed("Motional averaging using 'powder' method not handled for heteronuclear case");
    gammavals(second)=gam2;
    BaseDipolarAverager* avp=NULL;
    if (averaging_method==AV_TENSOR)
      avp=new TensorDipolarAverager(positions,gammavals,methyllist,verbose);
    else
      avp=new PowderDipolarAverager(positions,gammavals,methyllist,nzcw,verbose);
  
    List<double> d2local(nfirst);
    List<double> rminlocal(nfirst);
    List<double> dmaxlocal(nfirst);
    List<size_t> whichminlocal(nfirst);
    Matrix<bool> included(natoms,natoms);
    List<bool> atomincluded(natoms);
    char mol_info[20];

    if (!(avp->check_nooverlap(second)))
      throw Failed("Remote atom selection partially overlaps with averaging group");
    if (!(avp->check_nooverlap(first)))
      throw Failed("Source atom selection partially overlaps with averaging group");

//     bool ok=true;
//     for (size_t j=nsecond;ok && j--;) {
//       const int which=averager.which_set(second(j));
//       if (which>=0)
// 	ok=checkoverlap(second,methyllist(which)());
//     }
        
    static PowderMethod* powdm=NULL;
    if (!powdm) {
      switch (outputtype) {
      case OUT_MAP:
	powdm=new Hemisphere(nalpha,nbeta);
	break;
      case OUT_HIST:
	powdm=new SphericalZCW(nzcw,hemisphere);
	break;
      default:
	powdm=new PowderSingle; //!< no powder average for simple drss calculation
      }
    }
    const size_t orientations=powdm->orientations();
    if (orientations>1)
      cout << "Orientations: " << orientations << std::endl;
    
    powdm->reset();
    double junkscale;
    Euler powder(0,0,0);
    for (size_t orient=0;orient<orientations;orient++) {
      powdm->index_to_orientation(powder,junkscale,orient);
      
      d2local=0.0;
      rminlocal=0.0;
      included=false;
      atomincluded=false; //!< reset list showing which atoms considered
      
      const bool needrot=(powder.alpha!=0.0) || (powder.beta!=0.0);
      rmatrix3 rotmat;
      rotation_matrix(rotmat,powder);
      
      for (size_t i=0; i<nfirst; i++) { //!< loop over 'source' atoms
	const size_t seli=first(i);
	const mol_data& crysti(cryst(seli));
	const vector3 crysti_coords= needrot ? rotate(crysti.position,rotmat) : crysti.position;
	if (verbose)
	  cout << '\n';
	if (atomincluded(seli)) {
	  if (verbose)
	    cout << "Skipping " << crysti.type << " (equivalent site already considered)\n";
	  continue;
	}
	const bool ismethyli=avp->in_set(seli);
	
	for (size_t j=0; j<nsecond; j++){ //!< loop over 'remote' atoms
	  const size_t selj=second(j);
	  if (seli==selj)
	    continue;
	  const bool ismethylj=avp->in_set(selj);
	  const mol_data& crystj(cryst(second(j)));
	  const vector3 crystj_coords= needrot ? rotate(crystj.position,rotmat) : crystj.position;
	  const spherical sphco(crysti_coords-crystj_coords);
	  if (cryst.hasconnections())
	    strcpy(mol_info,cryst.samemolecule(first(i),first(j)) ? "(intramolecular)" : "(intermolecular)");
	  else
	    mol_info[0]='\0';
	  
	  double r=0.0;
	  double scale=1.0;
	  double locdss=0.0;
	  size_t weight=1;

	  //!< if averaging, then care if either atom is in methyl, otherwise only care about pairing *within* same methyl
	  const bool sameset=(avp->which_set(seli)==avp->which_set(selj));
	  const bool interesting= (averaging_method!=AV_NONE) ? (ismethyli || ismethylj) : (ismethyli && sameset);
	  
	  if (interesting) {
	    if (ismethyli && ismethylj) { //both in methyl
	      if (sameset) { //in *same* methyl
		r = (standard_short!=0.0) ? standard_short : sphco.r;
		strcpy(mol_info,"(intramethyl)");
		if (averaging_method!=AV_NONE)
		  scale=-0.5;
	      }
	      else {
		locdss=average_dss(*avp,cryst,included,atomincluded,methyllist(avp->which_set(seli)),methyllist(avp->which_set(selj)),rotmat);
		weight=3;
	      }
	    }
	    else { //one in, one out
	      if (ismethyli)
		locdss=average_dss(*avp,cryst,included,atomincluded,methyllist(avp->which_set(seli)),ExplicitList<1,size_t>(selj),rotmat);
	      else
		locdss=average_dss(*avp,cryst,included,atomincluded,ExplicitList<1,size_t>(seli),methyllist(avp->which_set(selj)),rotmat);
	      weight=3;
	    }
	  }
	  else { //no methyl business
	    if ((standard_short!=0.0) && (pairedatom(seli)==selj)) {
	      r=standard_short;
	      strcpy(mol_info,"(close H-H pair - using standard distance)");
	    }
	    else
	      r=sphco.r;
	  }
      
	  if (r!=0.0) { // 'real' r
	    if ((rminlocal(i)==0.0) || (r<rminlocal(i))) {
	      rminlocal(i)=r;
	      whichminlocal(i)=selj;
	    }
	    
	    double d=scale*dipolar_coupling(gam1,gam2,r*1e-10); //dipolar coupling

	    if (fabs(d)<global_cutoff)
	      d=0.0;
	    else {
	      if (isorientationmap())
		d*=legendre(2,sphco.theta); //!< include orientation
	      
	      if (verbose) {
		std::cout << "Distance from ";
		cryst.print_atom(seli);
		std::cout << " to ";
		cryst.print_atom(selj);
		std::cout << ": " << r<< " A   d=" << d << " Hz " << mol_info << '\n';
	      }
	    }
	    locdss=d*d;
	  } // otherwise use locdss and weight calculated earlier
	  d2+=weight*locdss;
	  d2local(i)+=weight*locdss;
	}
	if (d2==0) {
	  std::cerr << "No couplings found!\n";
	  exit(1);
	}
	if (isorientationmap()) {
	  rmatrix& curmap(drssmap(i));
	  if (!curmap) {
	    if (outputtype==OUT_HIST)
	      curmap.create(1,orientations);
	    else
	      curmap.create(nbeta,nalpha); //!< assume that alpha is most rapidly varying
	  }
	  BaseList<double> drow(curmap.row());
	  drow(orient)=sqrt(d2local(i));
	}
	else {
	  const double drss=sqrt(d2local(i));
	  dmaxlocal(i)=rminlocal(i) ? dipolar_coupling(gam1,gam2,rminlocal(i)*1e-10) : 0.0;
	  if (outfile)
	    fprintf(outfile,"%s\t%i\t%f\t%f\n",crysti.type,crysti.snumber,drss,rminlocal(i));
	  else  {
	    std::cout << crysti.type << '(' << crysti.snumber << "): d_rss=" << drss << " Hz";
	    if (rminlocal(i)) {
	      cout << "  Shortest (rigid) distance to ";
	      cryst.print_atom(whichminlocal(i));
	      cout << ": " << rminlocal(i) << " A (d=" << dmaxlocal(i) << " Hz)\n";
	    }
	    else
	      cout << '\n';
	  }
	}
      }
    }
    if (outputtype==OUT_TABLES) {
      rmintable(range(),groupi)=rminlocal;
      dmaxtable(range(),groupi)=dmaxlocal;
      for (size_t j=d2local.size();j--;)
	drsstable(j,groupi)=sqrt(d2local(j));
    }
  } //loop over groups
  
  if (isorientationmap()) {
    char tmpname[256];
    matlab_controller ctrl(outname,5);
    for (size_t i=0; i<nfirst; i++) {
      if (!drssmap(i).empty()) {
	snprintf(tmpname,sizeof(tmpname),"drss%" LCM_PRI_SIZE_T_MODIFIER "u",first(i)+1);
	cout << "Writing drss data for " << (first(i)+1) << '\n';
	ctrl.write(drssmap(i),tmpname);
      }
    }
  }
  else {
    if (outputtype==OUT_TABLES) {
      annotated_table(rmintable,"minr",cryst,first,seconds);
      annotated_table(dmaxtable,"maxd",cryst,first,seconds);
      annotated_table(drsstable,"drss",cryst,first,seconds);
    }
    if (!outfile && !d2 && (nfirst>1)) {
      double dav=sqrt(d2/nfirst);
      std::cout <<"\nAverage d_rss="<< dav << " Hz" << std::endl;
    }
  }
}
  
int main(int argc,const char *argv[])
{
  int count=1;
  try {
    
    cout.precision(5);
    char datafile[128]; //filename of PDB file
    const bool interactive=true; //!< disable non-interactive mode
    //    if (argc!=1) {
    //  if (argc<4) {
    //	std::cerr << "Syntax: " << argv[0] << " <observed nuclei> <remote nuclei> <filename> [<filename>...]\nExample nuclei selections: C - all Cs, C1 - carbon 1, C1* - all types of carbon 1 (e.g. C1B)\n";
    //	return 1;
    //  }
    //  interactive=false;
    // }
       
     
    if (interactive) {
      char first_type[MAX_TYPE_LEN],second_type[MAX_TYPE_LEN];

      getstring(argc,argv,count,"PDB file containing expanded (typically 3x3x3) structure? ", datafile,sizeof(datafile));

      std::cout << "N - atoms in no particular order - sort atoms by label in order to find unit cell\nA - atoms ordered by asymmetric unit (O:1 C:1 C:2 O:1 C:1 C:2 etc.)\nU - atoms ordered by unit cell (O:1 O:1 C:1 C:1 C:2 C:2 etc.)\n";
      const sort_t sorttype=(sort_t)getoption(argc,argv,count,"Atom ordering in PDB? ","NAU");

      verbose=getint(argc,argv,count, "Verbose level? ", 0,2,0);  

      ZMolStructure cryst;
 
      cryst.readPDB(datafile,sorttype); //read PDB file
      if (verbose>1) {
	std::cout << "Read in (atoms listed by internal index):\n";
	cryst.print();
      }
   
      cout << "Atom type examples: C* - all carbons; C1 - carbons at position 1; C1* - carbons at position 1 (including C1B, C1C, etc). Add isotope to resolve ambiguities e.g. 1H*, 2H*.\n";
      getstring(argc,argv,count,"From atom type or serial number (e.g. C*, C15, 124)? ", first_type,sizeof(first_type),"1H*");

      getstring(argc,argv,count,"To atom type (e.g. 15N*, 15N17A) [ENTER for unselective]? ",second_type,sizeof(second_type),"");
      
      const bool doradius=getlogical(argc,argv,count,"Select atoms close to site? ",false);

      if (doradius) {
	const int zrat=getint(argc,argv,count,"Ratio between Z and Z'? ",1);
	bool uselabels=true;
	if (!knowcastep)	  
	  uselabels=getlogical(argc,argv,count,"CASTEP output contains labels? ",true);

	cryst.set_zprime(zrat,uselabels);

	const int maxn=getint(argc,argv,count,"Maximum number of spins? ",10);
	
	getstring(argc,argv,count,"Output base filename? ",outname,sizeof(outname)-4);
	strcat(outname,".pdb");
	if (strcmp(outname,datafile)==0) {
	  std::cerr << "Refusing to overwrite input file!\n";
	  exit(2);
	}
	doextract(outname,datafile,cryst,first_type,second_type,maxn,uselabels);
	return 0;
      }
            
      //      if (gam1!=gam2)
      //	std::cerr << "Warning: motional averaging cannot be included for heteronuclear case\n";
      //else {
      std::cout << "N - none\nP - Powder averager\nT - Tensor averager\n";
      averaging_method=(av_t)getoption(argc,argv,count,"Motional averaging method? ","NPT");
      //	motionaveraging=getlogical(argc,argv,count,"Including motional averaging of XH3? ",true);
      standard_short=getfloat(argc,argv,count,"Normalised short internuclear distance (A - 0 if no normalisation)? ",standard_short);
	//}
      
      cout << "N - normal output\nT - output data tables\nM - drss map\nH - drss histogram\n";
      outputtype=(output_t)getoption(argc,argv,count,"Output type? ","NTMH",OUT_NORMAL);
      if (outputtype) {
	switch (outputtype) {
	case OUT_TABLES:
	  asymmerge=getlogical(argc,argv,count,"Merge asymmetric units (last character of atom type)? ",false);
	  break;
	case OUT_HIST:
	  nzcw=getint(argc,argv,count,"ZCW number? ");
	  break;
	default:
	  nalpha=getint(argc,argv,count,"alpha angles? ",30);
	  nbeta=getint(argc,argv,count,"beta angles? ",30);
	}
	getstring(argc,argv,count,"Output base filename? ",outname,sizeof(outname));
      }
      else {
	global_cutoff=getfloat(argc,argv,count,"Smallest coupling to include (Hz)? ",10.0);
	if (averaging_method==AV_POWDER)
	  nzcw=getint(argc,argv,count,"ZCW number? ",8);
      }
      
      analyse_molecule(cryst,first_type,second_type);
    }
    else { //loop through files - currently disabled

      std::cout << "Warning: motional averaging not implemented for command-line output\n";

      for (size_t i=3;i<argc;i++) {
	if (strlen(argv[i])+4>sizeof(outname))
	  throw Failed("Filename too long");
	strcpy(outname,argv[i]);
	char* ptr=strrchr(outname,'.');
	if (ptr && (ptr>strrchr(outname,'/')))
	  strcpy(ptr+1,"txt"); //place extension after .
	else
	  strcat(outname,".txt");
	  
	std::cout << "Analysing: " << argv[i] << std::endl;

	ZMolStructure cryst;
 
	cryst.readPDB(argv[i],def_sorttype); //read PDB file

	FILE* outfile=fopen(outname,"w");
	if (!outfile)
	  std::cerr << "Failed to open: " << outfile << '\n';
	else {
	  fprintf(outfile,"#Output from %s %s %s %s\n",argv[0],argv[1],argv[2],argv[i]);
	  fprintf(outfile,"#<serial no> <drss/Hz> <r_min/A>\n");
	  char* argv1=strdup(argv[1]);
	  char* argv2=strdup(argv[2]);
	  analyse_molecule(cryst,argv1,argv2,outfile);
	  fclose(outfile);
	}
      }
    }
    
  } catch (MatrixException& exc) {
    cerr << exc;
    return 1;
  }
  return 0;
}

// 	char orientstr[MAX_LINE];
// 	getstring(argc,argv,count,"Bond vector and atom specifying orientation (e.g. H6,H5,C1) [ENTER for none]? ",orientstr,sizeof(orientstr));  
// 	if (orientstr[0]) {

// 	  List<size_t> orientindices;
// 	  parse_index_list(orientindices,cryst,orientstr);
// 	  if (verbose)
// 	    std::cout << "Orientation indices: " << orientindices << '\n';
	  
// 	  if (orientindices.size()!=3) {
// 	    std::cerr << "Couldn't parse orientation string as 3 comma-separated atom identifiers\n";
// 	    exit(1);
// 	  }

// 	  cryst.reorient(orientindices,verbose);
// 	}

// void parse_index_list(List<size_t>& dest, const ZMolStructure& cryst, char* cptr)
// {
//   for (;;) {
//     char* tokptr=strtok(cptr,",");
//     if (tokptr==NULL)
//       break;
//     cptr=NULL;
//     dest.push_back(spec_to_index(cryst,tokptr));
//   }
// }

//     PlanarZCW powdmeth(nzcw);
//     double scale;
//     Euler powder;
    
//     while (powdmeth.next(powder,scale)) {
//       double meand=0.0;
//       rotation_matrix(rotmat,powder);
//       for (size_t i=whichi.size();i--;) {
// 	const size_t seli=whichi(i);
// 	for (size_t j=whichj.size();j--;) {
// 	  const size_t selj=whichj(j);
// 	  const vector3 diff=cryst(seli).position-cryst(selj).position;
// 	  const vector3 rotvec(rotate(diff,rotmat));
// 	  const spherical sphco(rotvec);
// 	  if (sphco.r==0.0)
// 	    throw Failed("spins are co-incident");
// 	  double d=dipolar_coupling(gamma1H,gamma1H,sphco.r*1e-10);
// 	  d*=legendre(2,sphco.theta);
// 	  meand+=d;
// 	}
//       }
//       dss+=scale*meand*meand/N;  
//    }  
//    static const double scalef=20.0/9;
//    dss*=scalef;
