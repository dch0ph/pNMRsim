/* Parse magres file to create spin system */

#include "ttyio.h"
//#include "space_T.h"
#include <errno.h>
#include "geometry.h"
#include "smartptr.h"
#include "ScratchList.h"
#include "NMR.h"
#include <map>
#include <string>
#include <set>
#include <sstream>
#include "MoleculeStructure.h"
#define USE_LIBCMATRIX 1
#include "magres.h"

// 2022-02-19  Bug fix at L1250

using namespace libcmatrix;
using namespace std;
using namespace MagRes;

#include "./geometry_utils.h"

static const double asym_zero_tol=1e-4; // asymmetry set to zero if less than this
static const double rad_to_deg=180.0/M_PI;

static const size_t MAX_LINE=1024;

// Ripped from io.f90 of CASTEP
static const double speed_light_si=299792458.0;
static const double planck_si=6.62606876e-34; // J s
static const double electron_mass_si = 9.10938188e-31; // kg
static const double hbar_si=planck_si/(2*M_PI);
static const double elementary_charge_si = 1.602176462e-19; // C
static const double mu_0_si = 4.0e-7*M_PI;
static const double epsilon_0_si = 1.0/(mu_0_si*speed_light_si*speed_light_si);
static const double fine_structure_si = (elementary_charge_si*elementary_charge_si)/(4.0*M_PI*epsilon_0_si*hbar_si*speed_light_si);
static const double joule = 1.0/(fine_structure_si*fine_structure_si*electron_mass_si*speed_light_si*speed_light_si);
static const double hertz = planck_si*joule;
static const double metre = electron_mass_si*speed_light_si*fine_structure_si/hbar_si;
static const double barn=metre*metre*1e-28;

bool asymmetric_unit_only=false;

struct quadrupole_properties {
  const char* label;
  double Q;
};

#define NMRSIM_QUAD_PROPERTIES 111
// Quadrupole moments of all nuclear isotopes, from P. Pyykko, J. Mol. Phys, 2008 106 1965-1974
// Units of mb, Q/10/fm^2

const quadrupole_properties quadrupole_props_Pyykko08[NMRSIM_QUAD_PROPERTIES]={
     {"2H", 2.860},
     {"6Li", -0.808},
     {"7Li", -40.1},
     {"9Be", 52.88},
     {"10B", 84.59},
     {"11B", 40.59},
     {"11C", 33.27},
     {"14N", 20.44},
     {"17O", -25.58},
     {"19F", -94.2},
     {"21Ne", 101.55},
     {"23Na", 104},
     {"25Mg", 199.4},
     {"27Al", 146.6},
     {"33S", -67.8},
     {"35S", 47.1},
     {"35Cl", -81.65},
     {"37Cl", -64.35},
     {"39K", 58.5},
     {"40K", -73},
     {"41K", 71.1},
     {"41Ca", -66.5},
     {"43Ca", -40.8},
     {"45Sc", -220},
     {"47Ti", 302},
     {"49Ti", 247},
     {"50V", 210},
     {"51V", -52},
     {"53Cr", -150},
     {"55Mn", 330},
     {"57Fe", 160},
     {"59Co", 420},
     {"61Ni", 162},
     {"63Cu", -220},
     {"65Cu", -204},
     {"67Zn", 150},
     {"69Ga", 171},
     {"71Ga", 107},
     {"73Ge", -196},
     {"75As", 314},
     {"77Se", 760},
     {"79Br", 313},
     {"81Br", 262},
     {"83Kr", 259},
     {"85Rb", 276},
     {"87Rb", 133.5},
     {"87Sr", 305},
     {"90Y", -125},
     {"91Zr", -176},
     {"93Nb", -320},
     {"95Mo", -22},
     {"97Mo", 255},
     {"99Tc", -129},
     {"99Ru", 79},
     {"101Ru", 457},
     {"105Pd", 660},
     {"113In", 759},
     {"115In", 770},
     {"119Sn", -132},
     {"121Sb", -543},
     {"123Sb", -692},
     {"127I", -696},
     {"131Xe", -114},
     {"133Cs", -3.43},
     {"135Ba", 160},
     {"137Ba", 245},
     {"138La", 450},
     {"139La", 200},
     {"141Pr", -58.9},
     {"143Nd", -630},
     {"145Nd", -330},
     {"147Pm", 740},
     {"147Sm", -259},
     {"149Sm", 75},
     {"151Eu", 903},
     {"153Eu", 2412},
     {"155Gd", 1270},
     {"157Gd", 1350},
     {"159Tb", 1432},
     {"161Dy", 2507},
     {"163Dy", 2648},
     {"165Ho", 3580},
     {"167Er", 3565},
     {"169Tm", -1200}, // Tm missing from Pyykko
     {"173Yb", 2800},
     {"175Lu", 3490},
     {"176Lu", 4970},
     {"177Hf", 3365},
     {"179Hf", 3793},
     {"181Ta", 3170},
     {"185Re", 2180},
     {"187Re", 2070},
     {"189Os", 856},
     {"191Ir", 816},
     {"193Ir", 751},
     {"197Au", 547},
     {"201Hg", 387},
     {"209Pb", -269},
     {"209Bi", -516},
     {"209Rn", 311},
     {"223Fr", 1170},
     {"223Ra", 1210}, // Ra missing from Pyykko
     {"227Ac", 1700},
     {"229Th", 4300},
     {"231Pa", -1720},
     {"233U", 3663},
     {"235U", 4936},
     {"237Np", 3886},
     {"241Pu", 5600},
     {"243Am", 4210},
     {"253Es", 6700},
};

double match_tolerance=0.0;

bool havequads=false;
bool haveshifts=false;
bool makespinsys=false;

typedef std::map<size_t,double> Cqmap_t;
Cqmap_t Cqmap;
std::set<size_t> usednuc;
  
typedef std::map< char, List<string> > labelmap_t;
bool havePDBlabels=false;
bool haveCASTEPlabels=false;
labelmap_t labelmap;

int quadorder=2;
int prec=4;

size_t verbose=1;

inline void dump(std::ostream& ostr, const double d[3])
{
  ostr << d[0] << "  " << d[1] << "  " << d[2] << '\n';
}

const char* makelabel(size_t nuc, size_t index)
{
  char tmp[10];
  const char* name=nuctolabel(nuc);
  while (isdigit(*name))
    name++;
  snprintf(tmp,sizeof(tmp)-1,"%c%" LCM_PRI_SIZE_T_MODIFIER "u",*name,index);
  return strdup(tmp);
}

const char* makelabel(size_t nuc, const char* rest)
{
  char tmp[10];
  const char* name=nuctolabel(nuc);
  while (isdigit(*name))
    name++;
  snprintf(tmp,sizeof(tmp)-1,"%c%s",*name,rest ? rest : ""); //!< fix-up to prevent 'null' being inserted
  return strdup(tmp);
}

std::ostream& operator<< (std::ostream& ostr, const rmatrix3& a)
{
  dump(ostr,a[0]);
  dump(ostr,a[1]);
  dump(ostr,a[2]);
  return ostr;
}

enum { ALL=0, SIMPSON_SINGLE, PNMRSIM_SINGLE, PNMRSIM_SUM, TABLES, GSIM };
int mode=PNMRSIM_SINGLE;

typedef std::map<size_t,std::string> sel_t;
sel_t selnuclei; //!< allowed nuclei (empty if all)

int lico=0;

char* getline(FILE* fp)
{
  static char linebuf[MAX_LINE];
  if (!fgets(linebuf,sizeof(linebuf),fp))
    return NULL;
  lico++;
  return linebuf;
}

void dumpnucleus(FILE* fp, size_t nuc)
{
  fprintf(fp," %s",nuctolabel(nuc));
  const sel_t::const_iterator sel(selnuclei.find(nuc));
  if (sel!=selnuclei.end()) {
    const char* qual=(sel->second).c_str();
    if (*qual)
      fprintf(fp,":%s",qual);
  }
}

void symmetrise(Matrix<double>& A)
{
  if (!issquare(A))
    throw NotSquare("symmetrise");
  for (size_t r=A.rows();r--;) {
    for (size_t c=r;c--;)
      A(r,c)=A(c,r)=0.5*(A(r,c)+A(c,r));
  }
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
      if (selnuclei.find(i)!=selnuclei.end()) //!< if in allowed set, accept automatically
	return i;
      if (nuc) {
	if (verbose)
	  std::cerr << "Ambiguous nucleus name: " << name << '\n';
	nuc=default_nucleus(name); //!< ambiguous label
	if ((nuc!=NULL_NUCLEUS) && verbose)
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
  if (iter!=nucmap.end())
    return iter->second;
  const size_t nuc=searchnuc(name);
  if (nuc)
    nucmap[strname]=nuc;
  return nuc;
}
 
typedef std::pair<size_t,const char*> parsed_atom_t;

parsed_atom_t parse_castep_atom(char* name)
{
  char* colon=strchr(name,':');
  if (colon)
    *colon++='\0';
  const size_t nuc=castep2nuc(name);
  if (nuc)
    return parsed_atom_t(nuc,makelabel(nuc,colon));
  return parsed_atom_t(NULL_NUCLEUS,NULL);
}

// template<class Container> struct IndexCompare : public std::binary_function<size_t,size_t,bool> {
//   const Container& data;
//   IndexCompare(const Container& datav) : data(datav) {}
//   bool operator()(size_t m, size_t n) { return (data(m)<data(n)); }
// };

// template<class Container> List<size_t> indexsort(const Container& data)
// {
//   List<size_t> index(slice(0,data.size()),mxflag::temporary);
//   std::sort(index.begin(),index.end(),IndexCompare<Container>(data));
//   return index;
// }

//typedef std::pair<std::string,size_t> nucident;

const char* nuctoelement(size_t nuc)
{
  const char* p=nuctolabel(nuc);
  while (isdigit(*p))
    p++;
  return p;
}

struct nucident {
  nucident(size_t nucv, const char* labelv, size_t indexv) 
    : nuc(nucv),
      index(indexv) {
    if (!labelv || (*labelv=='\0'))
      throw InvalidParameter("nucident: label can't be empty");
    label=labelv;
  }

   nucident(parsed_atom_t spec, size_t indexv)
     : nuc(spec.first),
       index(indexv) {
     if (!(spec.second) || (*(spec.second)=='\0'))
      throw InvalidParameter("nucident: label can't be empty");
     label=spec.second;
   }

  size_t nuc;
  std::string label;
  size_t index;

  bool istype(const nucident& a) const {
    if (nuc!=a.nuc)
      return false;
    if (*(a.label.c_str()))
      return (label==a.label);
    return (index==a.index);
  }

  bool operator== (const nucident& a) const
  { return (nuc==a.nuc) && (index==a.index) && (label==a.label); }

  bool operator!= (const nucident& a) const
  { return (nuc!=a.nuc) || (index!=a.index) || (label!=a.label); }
};

std::ostream& operator<< (std::ostream& ostr, const nucident& a)
{
  ostr << nuctolabel(a.nuc);
  if (a.label.size())
    ostr << ':' << a.label;
  return ostr << '(' << a.index << ')';
}

class ScaledComparator {
public:
  ScaledComparator(double scalefacv) { scaledtol_=fabs(scalefacv*match_tolerance); }

  bool operator()(double x, double y) const {
    if (scaledtol_)
      return (fabs(x-y)<scaledtol_);
    return (x==y);
  }
private:
  double scaledtol_;
};

struct Comparator {
  bool operator()(double x, double y) const {
    if (match_tolerance) {
      const double compval=std::max(fabs(x),fabs(y))*match_tolerance;
      return compval ? ( fabs(x-y) < compval) : true;
    }
    return (x==y);
  }
};

struct interaction {

  interaction(const Matrix<double>& Av)
    : A(Av), iso(0.0), aniso(0.0), asym(0.0) {}

  interaction()
    : iso(0.0), aniso(0.0), asym(0.0) {}

  void process(bool zeroiso =false);

  enum { ISO =0, ANISO, ASYM, ALPHA, BETA, GAMMA };
  rmatrix A;
  double iso,aniso,asym;
  Euler angles;

  void rotate(const rmatrix&);

  double getvalue(Int2Type<ISO>) const { return iso; }
  double getvalue(Int2Type<ANISO>) const { return aniso; }
  double getvalue(Int2Type<ASYM>) const { return asym; }
  double getvalue(Int2Type<ALPHA>) const { return angles.alpha*rad_to_deg; }
  double getvalue(Int2Type<BETA>) const { return angles.beta*rad_to_deg; }
  double getvalue(Int2Type<GAMMA>) const { return angles.gamma*rad_to_deg; }

  template<int N> void dump(FILE* fp) const { fprintf(fp,"%.*g",prec,getvalue(Int2Type<N>())); }

  void dumpangles(FILE* fp, char sep =' ') const {
    dump<ALPHA>(fp); fputc(sep,fp);
    dump<BETA>(fp); fputc(sep,fp);
    dump<GAMMA>(fp);
  }

  bool operator==(const interaction& a) const { return equalpas(a) && equalangles(a); }
  bool equalpas(const interaction& a) const;
  bool equalangles(const interaction& a) const;
};

class averager {
public:
  averager() {}
  void add(interaction&);
  void perform();
private:
  List<interaction*> storep;
  Matrix<double> summat;
  //  List<bool> negate;
  //vector3 values;
};

//averager::averager()
//  : values(0.0,0.0,0.0) {}

void averager::add(interaction& inter)
{
  storep.push_back(&inter);
  //  double lxx,lyy,lzz;
  //negate.push_back(inter.aniso<0.0);
  // anisotropy_to_cartesian(lxx,lyy,lzz,inter.iso,fabs(inter.aniso),inter.asym);

  //values+=vector3(lxx,lyy,lzz);
  const space_T A_PAS(spatial_tensor(inter.iso,inter.aniso,inter.asym));
  const space_T A_MF(rotate(A_PAS,inter.angles));
  Matrix<double> Am_MF;
  tensor_to_matrix(Am_MF,A_MF);
  summat+=Am_MF;
}

void averager::perform()
{
  if (storep.empty())
    return;
    //    throw Failed("averager: nothing to average over!");
//   values/=storep.size();
  double liso,laniso,lasym;
  Euler lF;
//   cartesian_to_anisotropy(liso,laniso,lasym,values.x,values.y,values.z);
  summat*=1.0/storep.size();
  cartesian_to_PAS_symmetric(liso,laniso,lasym,lF,summat);
  static const double minussqrt3=-sqrt(3.0);
  static const double sqrt3over2=sqrt(3.0/2.0);
  liso*=minussqrt3;
  laniso*=sqrt3over2;

  for (size_t i=storep.size();i--;) {
    interaction& curint(*(storep(i)));
    curint.A=summat; //!< replace cartesian rep
    curint.iso=liso;
    //curint.aniso=(negate(i) ? -laniso : laniso);
    curint.aniso=laniso;
    curint.asym=lasym;
    curint.angles=lF;
  }
}

void interaction::rotate(const rmatrix& R)
{
  A.unitary_simtrans(R);
}

bool interaction::equalpas(const interaction& a) const
{
  static ScaledComparator areequalasym(1.0);
  static Comparator areequal;
  const bool miso=areequal(iso,a.iso);
  const bool maniso=areequal(aniso,a.aniso);
  const bool masym=areequalasym(asym,a.asym);
  if (verbose>1)
    std::cout << miso << maniso << masym;
  return miso && maniso && masym;
}
  
bool interaction::equalangles(const interaction& a) const
{
  static ScaledComparator areequal(2.0*M_PI);
  return areequal(angles.alpha, a.angles.alpha) && areequal(angles.beta, a.angles.beta) && areequal(angles.gamma, a.angles.gamma);
}

  struct nucdesc {
    enum { CSA =0, QUAD };
    nucdesc(nucident identv, const vector3& positionv, const char* extlabelv =NULL, size_t aunitv =0)
      : ident(identv),
	multiplicity(1U),
	position(positionv),
	aunit(aunitv)
	   //	extlabel(extlabelv)
    { isquad=(spin(ident.nuc).deg()>2); 
      if (extlabelv)
	label=extlabelv;
      else {
	std::ostringstream ostr(std::ostringstream::out);
	ostr << nuctolabel(identv.nuc);
	const char* intlabel=identv.label.c_str();
	if (*intlabel)
	  ostr << ':' << intlabel;
	//	ostr << '_' << identv.index;
	label=ostr.str();
      }
    }

    nucident ident;
    bool isquad;
    smartptr<interaction,false> csap,quadp;
    size_t multiplicity;
    //    const char* extlabel;
    std::string label;
    vector3 position;
    size_t aunit;
    //    const char* label() const { return extlabel ? extlabel : ident.label.c_str(); }
    
    //bool operator== (const nucdesc& a) const { return (ident==a.ident); }
    //bool operator!= (const nucdesc& a) const { return (ident!=a.ident); }
    bool operator== (const nucident& a) const { return (ident==a); }
    bool operator!= (const nucident& a) const { return (ident!=a); }

    size_t nucleus() const { return ident.nuc; }
    size_t index() const { return ident.index; }
    void dump(FILE*, int, const char* ="") const;
    void tableintro(FILE*) const;
    void dumpquadtable(FILE* fp) const;
    void dumpshifttable(FILE* fp) const;

    template<int N> void dump(FILE* fp, Int2Type<CSA>, Int2Type<N> which) const { csap->dump<N>(fp); }
    void dump(FILE* fp, Int2Type<CSA>, Int2Type<interaction::ANISO>, bool printp =true) const {
      csap->dump<interaction::ANISO>(fp);
      if (printp)
	fputc('p',fp);
    }
    void dump(FILE* fp, Int2Type<CSA>, Int2Type<interaction::ISO>, bool printp =true) const {
      csap->dump<interaction::ISO>(fp);
      if (printp)
	fputc('p',fp);
    }    
    template<int N> void dump(FILE* fp, Int2Type<QUAD>, Int2Type<N> which) const { quadp->dump<N>(fp); }

    void rotate(const rmatrix3&, const rmatrix&);
  };

typedef List<nucdesc> store_t;
store_t atoms;

bool sortonlabel(const nucdesc& a, const nucdesc& b) { return (a.label.compare(b.label)<0); }

void nucdesc::rotate(const rmatrix3& R3, const rmatrix& R)
{
  position=libcmatrix::rotate(position,R3);
  if (csap.get()) {
    //   std::cout << "Tensor before\n" << A;
    csap->rotate(R);
    //  std::cout << "Tensor after\n" << A;
    csap->process();
  }
  if (quadp.get()) {
    quadp->rotate(R);
    quadp->process(true);
  }
}

struct CompareInteraction : public std::binary_function<interaction,interaction,bool> 
{
  bool operator()(const interaction& a, const interaction& b) const { return (a==b); }
};

struct ComparePAS : public std::binary_function<interaction,interaction,bool> 
{
  bool operator()(const interaction& a, const interaction& b) const { return a.equalpas(b); }
};

template<class Func> struct CompareNuc {
  CompareNuc(const nucdesc& nucv) : nuc(nucv) {}
  bool operator()(const nucdesc& a) const { return isequal(a.csap,nuc.csap,'C') && isequal(a.quadp,nuc.quadp,'Q'); }
  bool isequal(const smartptr<interaction,false>& a, const smartptr<interaction,false>& b, char inter) const
  {
    if ((!a) ^ (!b)) //!< one set and one not?
      return false;
    if (!a)
      return true;
    const bool res=func(*a,*b);
    if (verbose>1) 
      std::cout << inter << ' ' << (res ? "match " : "not matched ");
    return res;
  }  
  Func func;
  const nucdesc& nuc;
};

template<class Func> void compress(List<nucdesc>& atoms)
{
  List<nucdesc> newatoms;
  for (size_t i=0;i<atoms.size();i++) {
    const List<nucdesc>::iterator iter(std::find_if(newatoms.begin(),newatoms.end(),CompareNuc<Func>(atoms(i))));
    if (iter==newatoms.end()) {
      newatoms.push_back(atoms(i));
      if (verbose>1)
	std::cout << "Retaining atom " << (i+1) << "  " << atoms(i).label << '\n';
    }
    else {
      if (verbose>1)
	std::cout << "Merging atom " << (i+1) << "  " << atoms(i).label << '\n';      
      (iter->multiplicity)++; //!< increase multiplicity count
    }
  }	 
  atoms.swap(newatoms);
}

void dump(FILE* fp, const BaseList<size_t>& a)
{
  const bool assum=(mode==PNMRSIM_SUM);
  if (a.size()>1)
    fputc(assum ? '|' : '{',fp);
  for (size_t i=0;i<a.size();i++) {
    if (i)
      fputc(',',fp);
    fprintf(fp,"%" LCM_PRI_SIZE_T_MODIFIER "u",a(i));
  }
  if (a.size()>1)
    fputc(assum ? '|' : '}',fp);
}
      
template<int Int, int N> void dumplist(FILE* fp, const BaseList<nucdesc>& list, Int2Type<Int> whichi, Int2Type<N> which)
{
  const bool assum=(mode==PNMRSIM_SUM);
  try {
    if (list.size()>1)
      fputc(assum ? '|' : '{',fp);
    for (size_t i=0;i<list.size();i++) {
      if (i)
	fputc(',',fp);
      list(i).dump(fp,whichi,which);
    }
    if (list.size()>1)
      fputc(assum ? '|' : '}',fp);
  } catch (Failed&) { //!< catch failure to de-reference smartptr (missing interaction)
    std::cerr << "Interactions must be defined uniformly for all nuclei in this output mode\n";
    exit(1);
  }
}

template<int Int> void dumplistangles(FILE* fp, const BaseList<nucdesc>& list, Int2Type<Int> whichi)
{  
  dumplist(fp,list,whichi,Int2Type<interaction::ALPHA>()); fputc(' ',fp);
  dumplist(fp,list,whichi,Int2Type<interaction::BETA>()); fputc(' ',fp);
  dumplist(fp,list,whichi,Int2Type<interaction::GAMMA>());
}

void nucdesc::dump(FILE* fp, int ind, const char* intro) const
{
  fprintf(fp,"#%s ",label.c_str());
  fprintf(fp,"index %" LCM_PRI_SIZE_T_MODIFIER "u ",ident.index);
  fprintf(fp,"at (%g,%g,%g) A\n",position.x,position.y,position.z);
  if (csap.get()) {
    fprintf(fp,"%s\tshift %i ",intro,ind);
    dump(fp,Int2Type<CSA>(),Int2Type<interaction::ISO>()); fputc(' ',fp);
    dump(fp,Int2Type<CSA>(),Int2Type<interaction::ANISO>()); fputc(' ',fp);
    dump(fp,Int2Type<CSA>(),Int2Type<interaction::ASYM>()); fputc(' ',fp);
    csap->dumpangles(fp);
    fputc('\n',fp);
  }
  if (quadp.get()) {
    fprintf(fp,"%s\tquadrupole %i %i ",intro,ind,quadorder);
    dump(fp,Int2Type<QUAD>(),Int2Type<interaction::ANISO>()); fputc(' ',fp);
    dump(fp,Int2Type<QUAD>(),Int2Type<interaction::ASYM>()); fputc(' ',fp);
    quadp->dumpangles(fp);
    fputc('\n',fp);
  }
}

void nucdesc::tableintro(FILE* fp) const
{
//   const char* llabel=extlabel;
//   if (!llabel)
//     llabel=nuctolabel(nucleus());
  fprintf(fp,"%s\t%" LCM_PRI_SIZE_T_MODIFIER "u\t%" LCM_PRI_SIZE_T_MODIFIER "u\t",label.c_str(),index(),multiplicity);
}

void nucdesc::dumpquadtable(FILE* fpq) const
{
  if (fpq && quadp.get()) {
    tableintro(fpq);
    dump(fpq,Int2Type<QUAD>(),Int2Type<interaction::ANISO>()); fputc('\t',fpq);
    dump(fpq,Int2Type<QUAD>(),Int2Type<interaction::ASYM>()); fputc('\t',fpq);
    quadp->dumpangles(fpq,'\t');
    fputc('\n',fpq);
  }
} 

void nucdesc::dumpshifttable(FILE* fps) const
{
  if (fps && csap.get()) {
    tableintro(fps);
    dump(fps,Int2Type<CSA>(),Int2Type<interaction::ISO>(),false); fputc('\t',fps);
    dump(fps,Int2Type<CSA>(),Int2Type<interaction::ANISO>(),false); fputc('\t',fps);
    dump(fps,Int2Type<CSA>(),Int2Type<interaction::ASYM>()); fputc('\t',fps);
    csap->dumpangles(fps,'\t');
    fputc('\n',fps);
  }
} 

void interaction::process(bool zeroiso)
{
  symmetrise(A);  
  cartesian_to_PAS_symmetric(iso,aniso,asym,angles,A);
  if (zeroiso)
    iso=0.0;
}

void parse_index(List<size_t>& dest, const List<nucdesc>& atoms, char* tokptr, bool isspecific =true)
{
  if (!tokptr)
    throw InvalidParameter("parse_index: NULL input string");
  char* colonptr=strchr(tokptr,':');
  size_t nuc=NULL_NUCLEUS;
  const char* labelptr=NULL;
  unsigned int ind=1;
  bool isvalid=false;
  
  if (colonptr) {
    *colonptr++='\0';
    if (!isdigit(*colonptr))
      labelptr=colonptr;
    else
      isvalid=(sscanf(colonptr,"%u",&ind)==1);
    if (isvalid) {
      nuc=isdigit(*tokptr) ? labeltonuc(tokptr) : castep2nuc(tokptr);
      if (nuc==NULL_NUCLEUS) {
	std::cerr << "Nucleus label is ambiguous: " << tokptr << ".  Use isotope specification e.g. 95Mo\n";
	exit(1);
      }
    }
  }
  else {
    if (sscanf(tokptr,"%u",&ind)==1) {
      if ((ind<1) || (ind>atoms.size())) {
	std::cerr << "Index " << ind << " out range (1-" << atoms.size() << ")\n";
	exit(1);
      }
      dest.push_back(ind-1);
      return;
    }
    else {
      nuc=castep2nuc(tokptr);
      if (nuc==NULL_NUCLEUS) {
	std::cerr << "Nucleus label is ambiguous: " << tokptr << ".  Use isotope specification e.g. 95Mo\n";
	exit(1);
      }
      labelptr=tokptr+1;
      ind=0;
      isvalid=true;
    }
  }

  if (!isvalid) {
    std::cerr << "Failed to parse " << tokptr << " as either <element>:<index>|<label> or <index>\n";
    exit(1);
  }
  if (!labelptr) {
    labelptr=makelabel(nuc,ind);
    ind=0;
  }
  const nucident nucsel(parsed_atom_t(nuc,labelptr),ind);      
  
  if (isspecific) {
    const nucdesc* matchp=std::find(atoms.begin(),atoms.end(),nucsel);      
    if (matchp==atoms.end()) {
      std::cerr << "Failed to find atom " << nucsel << '\n';
      exit(1);
    }
    dest.push_back(matchp-atoms.begin());
  }
  else {
    for (size_t i=0;i<atoms.size();i++) {
      if (atoms(i).ident.istype(nucsel))
	dest.push_back(i);
    }
  }
}

void parse_index_list(List<size_t>& dest, const List<nucdesc>& atoms, char* cptr)
{
  for (;;) {
    char* tokptr=strtok(cptr,",");
    if (tokptr==NULL)
      break;
    cptr=NULL;
    parse_index(dest,atoms,tokptr,true);
  }
}

void parse_generic_index_list(Matrix<size_t>& dest, const List<nucdesc>& atoms, char* cptr)
{
  List<size_t> tmp;
  size_t len=0;
  for (;;) {
    char* tokptr=strtok(cptr,",");
    if (tokptr==NULL)
      break;
    cptr=NULL;
    const size_t cursize=tmp.size();
    parse_index(tmp,atoms,tokptr,false);
    if (len==0)
      len=tmp.size();
    else {
      if ((tmp.size()-cursize)!=len) {
	std::cerr << "Label specifications give lists of different lengths\n";
	exit(1);
      }
    }      
  }
  const size_t ncols=tmp.size() / len;
  dest.create(len,ncols);
  dest.row()=tmp;
}
  
void parse_super_index_list(List< Matrix<size_t> >& superdest, const List<nucdesc>& atoms, char* rawcptr)
{
  char* supersaveptr=NULL;
  for (;;) {
    char* cptr=strtok_r(rawcptr,";",&supersaveptr);
    if (cptr==NULL) {
      //     if (verbose)
      //	std::cout << "Average list: " << superdest << '\n';
      return;
    }
    rawcptr=NULL;
    superdest.push_back();
    parse_generic_index_list(superdest.back(),atoms,cptr);
  }
}

bool isselected(size_t nuc)
{
  if (nuc==NULL_NUCLEUS)
    return false;
  if (selnuclei.empty())
    return true;
  return (selnuclei.find(nuc)!=selnuclei.end());
}

void add_interaction(smartptr<interaction,false>& curint, const Matrix<double>& A, bool isquad)
{
  if (isquad)
    havequads=true;
  else
    haveshifts=true;
  
  interaction* intp=new interaction(A);
  curint.reset(intp);
  intp->process(isquad);
  const space_T A_PAS(spatial_tensor(intp->iso,intp->aniso,intp->asym));
  if (verbose>1)
    std::cout << (isquad ? 'Q' : 'C') << " tensor in PAS:\n" << A_PAS << '\n';
  const space_T A_MF(rotate(A_PAS,intp->angles));
  if (verbose>1)
    std::cout << "Tensor in MF:\n" << A_MF << '\n';
  const space_T A_MFd(A2(A));
  if (verbose>1)
    std::cout << "Tensor in MF (direct):\n" << A_MFd << '\n';      
}

void add_interaction(smartptr<interaction,false>& curint, const MagRes::tensor_t& tens, bool isquad, double scale =1.0)
{
  Matrix<double> A(3,3);

  for (size_t r=3;r--;)
   for (size_t c=3;c--;)
     A(r,c)=scale*tens[r][c];

  add_interaction(curint,A,isquad);
}

double quadscale(const nucdesc& curatom)
{
  if (!(curatom.isquad)) {
    if (verbose>1)
      std::cerr << "Ignoring quadrupole specification for non-quadrupolar nucleus\n";
    return 0;
  }

  const size_t nuc=curatom.nucleus();
  const Cqmap_t::const_iterator iter(Cqmap.find(nuc));
  if (iter==Cqmap.end()) {
    const char* curlabel=nuctolabel(nuc);
    for (size_t i=0;i<NMRSIM_QUAD_PROPERTIES;i++) {
      const quadrupole_properties& curquad(quadrupole_props_Pyykko08[i]);
      if (strcmp(curquad.label,curlabel)==0) {
	const double eQ=curquad.Q*barn/(1000.0*hertz); //!< convert from atomic to SI
	std::cout << "Setting eQ for " << curlabel << ": " << (eQ*1e-6) << " MHz\n";
	Cqmap[nuc]=eQ;
	return eQ;
      }
    }
    std::cerr << "Cq for " << curlabel << " not known!\n";
    return 1;
  }
  return iter->second;
}

//std::map<const char*, size_t> label_usage;

const char* get_PDB_label(const char* nucleus, size_t ind)
{
  const char* PDBlabel=NULL;
  const labelmap_t::const_iterator iter(labelmap.find(nucleus[0]));
  size_t useind=ind;
  if (iter!=labelmap.end()) {
    const List<string>& curlabels(iter->second);
    if (asymmetric_unit_only) {
      useind = ((ind-1) % curlabels.size()) + 1;
      //if (label_usage.empty())
      //	label_usage.create(labelmap.size(),
      //label_usage(ind-1)++;
    }
    if (useind<=curlabels.size())
      PDBlabel=curlabels(useind-1).c_str();
  }
  if (!PDBlabel)
    std::cerr << "Failed to find label from PDB for atom " << (ind+1) << ": " << nucleus << ' ' << ind << '\n';		
  return PDBlabel;
}

void read_old_magres(const char* fname)
{
  FILE* fp=fopen(fname,"r");
  if (!fp) {
    std::cerr << "magres2pNMRsim: failed to open " << fname << '\n';
    exit(1);
  }

  int lico=0;
  enum State { BeforeAtom =1, BeforeData };
  State cstate=BeforeAtom;
	
  nucdesc* curatomp=NULL;
  for (;;) {
    char* cptr;
    if ((cptr=getline(fp))==NULL)
      break;

    char nucleus[MAX_LINE];
    unsigned int ind;
    lico++;
    
    switch (cstate) {
    case BeforeAtom:
//       if (strstr(cptr," Barn")) { //!< found quadrupole spec?
// 	int mass=0;
// 	double Cq;
// 	if (sscanf(cptr," %s Isotope %i Q = %lg Barn",nucleus,&mass,&Cq)!=3) {
// 	  char junk[256];
// 	  if (sscanf(cptr," %s %s %s Q = %lg Barn",nucleus,junk,junk,&Cq)!=4) {
// 	    std::cerr << "Failed to parse quadrupole specification: " << cptr;
// 	    continue;
// 	  }
// 	}
// 	const parsed_atom_t nucspec(parse_castep_atom(nucleus));
// 	const size_t nuc=nucspec.first;
// 	if (nuc) {
// 	  const int mass2=atoi(nuctolabel(nuc)); //!< get mass;
// 	  if ((mass!=0) && (mass!=mass2))
// 	    std::cerr << "Mismatch between CASTEP mass (" << mass2 << ") and assumed mass (" << mass << ") g mol^-1.  Specification ignored\n";
// 	  else {
// 	    Cq*=barn/hertz; //!< convert from atomic to SI
// 	    Cqmap[nuc]=Cq;
// 	    std::cout << "Setting Cq for " << nuctolabel(nuc) << ": " << (Cq*1e-6) << " MHz\n";
// 	  }
// 	}
// 	continue;
//       }
      switch (sscanf(cptr,"Atom: %s %ui",nucleus,&ind)) {
      case 2: {
	parsed_atom_t nucspec(parse_castep_atom(nucleus));
	const size_t nuc=nucspec.first;
	if (nucspec.second)
	  haveCASTEPlabels=true;
	
	const bool allow=isselected(nuc);
	if (allow) {
	  if (!nucspec.second) {
	    nucspec.second=makelabel(nuc,ind);
	    ind=0;
	  }
	  if (verbose) {
	    std::cout << (allow ? "Found " : "Ignoring ") << nucleus << " (" << nuc << ")";
	    if (nucspec.second) 
	      std::cout << ':' << nucspec.second;
	    std::cout << ' ' << ind << '\n';
	  }
	  const nucident cident(nucspec,ind);
	  const store_t::iterator iter(std::find(atoms.begin(),atoms.end(),cident));
	  if (iter==atoms.end()) {
	    if (verbose>1)
	      std::cout << "(New atom)\n";
	    usednuc.insert(nuc);
	    const char* PDBlabel=NULL;
	    if (havePDBlabels && !makespinsys) {
	      static bool donelabelwarn=false;
	      if (haveCASTEPlabels && !donelabelwarn) {
		std::cerr << "Labels supplied both in CASTEP file and PDB - overwriting CASTEP label\n";
		donelabelwarn=true;
	      }
	      PDBlabel=get_PDB_label(nucleus,ind);
// 	      const labelmap_t::const_iterator iter(labelmap.find(nucleus[0]));
// 	      if (iter!=labelmap.end()) {
// 		const List<string>& curlabels(iter->second);
// 		if (ind<=curlabels.size())
// 		  PDBlabel=curlabels(ind-1).c_str();
// 	      }
// 	      if (!PDBlabel)
// 		std::cerr << "Failed to find label from PDB for atom " << (ind+1) << ": " << nucleus << ' ' << ind << '\n';		
	    }
	    vector3 pos;
	    char* linebuf;
	    if ((linebuf=getline(fp)) && (linebuf=getline(fp))) {
	      char label[256];
	      int junki;
	      if (sscanf(linebuf,"%s %i Coordinates %lg %lg %lg",label,&junki,&(pos.x),&(pos.y),&(pos.z))!=5)
		std::cerr << "Error parsing coordinates line: " << linebuf << " (ignored)\n";
	    }
	    else {
	      std::cerr << "File ended before coordinates found\n";
	      exit(1);
	    }
	    atoms.push_back(nucdesc(cident,pos,PDBlabel));
	    curatomp=&(atoms.back());
	  }
	  else {
	    if (verbose>1)
	      std::cout << "(Atom already encountered)\n";
	    curatomp=&(*iter);
	  }
	  cstate=BeforeData;
	}
      }
	break;
      case 0:
	break;
      default:
	std::cerr << "Failed to parse Atom description at line " << lico << '\n';
	exit(1);
      }
      break;
    case BeforeData:
      if (curatomp==NULL)
	throw InternalError("curatomp unset");
      if (strncmp(cptr,"TOTAL ",6))
	break;
      cptr+=6;
      smartptr<interaction,false>* curintp=NULL;
      bool ignore=false;
      double scale=1.0;
      bool isquad=false;
      if (strncasecmp(cptr,"Shielding",9)==0) {
	curintp=&(curatomp->csap);
	scale=-1.0; // "shielding" -> "shift";
      }
      else {
	if (strncasecmp(cptr,"tensor",6)==0) {
	  isquad=true;
	  curintp=&(curatomp->quadp);
	  if (!(curatomp->isquad)) {
	    if (verbose>1)
	      std::cerr << "Ignoring quadrupole specification for non-quadrupolar nucleus\n";
	    ignore=true;
	  }
	  else {
	    const Cqmap_t::const_iterator iter(Cqmap.find(curatomp->nucleus()));
	    if (iter==Cqmap.end()) {
	      std::cerr << "Cq for " << nuctolabel(curatomp->nucleus()) << " not set\n";
	      exit(1);
	    }
	    scale=iter->second;
	  }
	}
      }
      if (!curintp) {
	std::cerr << "Unrecognised interaction: " << cptr << " at line " << lico << '\n';
	exit(1);
      }
      if (!ignore) {
	Matrix<double> A(3,3);
	for (size_t i=0;i<3;) {
	  char* linebuf;
	  if ((linebuf=getline(fp))==NULL) {
	    std::cerr << "End of input file reading tensor\n";
	    exit(1);
	  }
	  BaseList<double> currow(A.row(i));
	  switch (sscanf(linebuf,"%lg %lg %lg",&(currow(0)),&(currow(1)),&(currow(2)))) {
	  case 0: case EOF://!< ignore blank lines
	    break; 
	  case 3:
	    i++;
	    break;
	  default:
	    std::cerr << "Failed to read 3 matrix components at line " << lico << '\n';
	    exit(1);
	  }
	}
	A*=scale;
	add_interaction(*curintp,A,isquad);
      }
      cstate=BeforeAtom;
      curatomp=NULL;
      break;
    }
  }  
  fclose(fp);
}

int main(int argc, const char* argv[])
{
  cmatrix_euler_controller.verify_tolerance=1e-6; //!< actively check Euler angles

  const size_t MAX_FILE=1024;
  char fname[MAX_FILE];

  int count=1;
  getstring(argc,argv,count,"Input filename (including .magres)? ",fname,sizeof(fname));

  char outname[MAX_FILE];
  getstring(argc,argv,count,"Output filename (including .in/.inc extension if relevant)? ",outname,sizeof(outname));
  char baseoutname[MAX_FILE];
  strcpy(baseoutname,outname);
  char* dirptr=strrchr(baseoutname,'\\');
  char* dotptr=strrchr(baseoutname,'.');
  if (dotptr && ((dirptr==NULL) || (dotptr>dirptr)))
    *dotptr='\0';
      
  if (strcmp(outname,fname)==0) {
    std::cerr << "Output and input files must be different\n";
    return 1;
  }
    
  size_t nlabels=0;
  char pdbname[256];

  getstring(argc,argv,count,"PDB file with labels or spin system specification (if any)? ",pdbname,sizeof(pdbname));

  verbose=getint(argc,argv,count,"Verbosity (0 - quiet, 2 - debugging)? ",1);

  MoleculeStructure molstruct(2);

  if (pdbname[0]) {

    enum {  PDB_CELL=0, PDB_MULTIPLE, PDB_SYS, };
    std::cout << "A - PDB contains asymmetric unit\nM - PDB contains multiple repeating asymmetric units\nS - PDB defines spin system\n";
    int pdb_mode=getoption(argc,argv,count,"PDB contains? ","AMS",PDB_MULTIPLE);

    //    if (!table_oriented)
    //makespinsys=getlogical(argc,argv,count,"Use PDB file to specify spin system? ");
    makespinsys=(pdb_mode==PDB_SYS);

    const bool findunitcell=(pdb_mode==PDB_MULTIPLE);
    molstruct.read_PDB(pdbname,'\0',findunitcell);
    if (pdb_mode==PDB_CELL)
      asymmetric_unit_only=true; //!< flag labels need to repeat
    nlabels=molstruct.size();
    if (nlabels==0) {
      std::cerr << "Failed to read PDB file: " << pdbname << '\n';
      exit(1);
    }
    if (makespinsys) {
      bool hasinfo=(molstruct.front().residue!=0);
      for (size_t i=1;i<nlabels;i++) {
	const bool haslinfo=(molstruct(i).residue!=0);
	if (hasinfo!=haslinfo) {
	  std::cerr << "PDB file contains inconsistent information on asymmetric unit (in residue number).\n";
	  exit(2);
	}
      }
      std::cout << "PDB file contains information on asymmetric unit: " << (hasinfo ? "Yes\n" : "No\n");
      havePDBlabels=hasinfo;
    }
    else {
      havePDBlabels=true;
      for (size_t i=0;i<nlabels;i++) {	
	const char* curlabel=molstruct(i).type;
	char curatom=*curlabel;
	if (isdigit(curatom))
	  curatom=curlabel[1];
	if ((curatom<'A') && (curatom>'Z'))
	  std::cerr << "Warning: Bad label in PDB for atom " << i << ": " << curlabel << '\n';
	else {
	  List<string>& curlabels(labelmap[curatom]);
	  if (verbose>1) 
	    std::cout << (i+1) << " " << curatom << " " << molstruct(i).type << '\n';
	  curlabels.push_back(string(curlabel));
	}
      }
    }
  }

  char nucspec[MAX_LINE]="";
  if (!makespinsys) {
    getstring(argc,argv,count,"Comma separated list of nuclei to select (ENTER for all, use 95Mo etc. to specify isotope, O:c for central transition only)? ",nucspec,sizeof(nucspec));

    if (nucspec[0]) {
      char* cptr=nucspec;
      char* tokptr;
      while ((tokptr=strtok(cptr,","))) {
	char* colonptr=strchr(tokptr,':');
	if (colonptr) {
	  *colonptr++='\0';
	  if (*colonptr=='\0') {
	    std::cerr << ": qualifier is empty\n";
	    exit(1);
	  }
	}
	const size_t nuc=isdigit(*tokptr) ? labeltonuc(tokptr) : castep2nuc(tokptr);
	if (nuc==NULL_NUCLEUS) {
	  std::cerr << "Nucleus label is ambiguous: " << tokptr << ".  Use isotope specification e.g. 95Mo\n";
	  exit(1);
	}
	selnuclei[nuc]=colonptr ? colonptr : "";
	cptr=NULL;
      }
    }
  }
  if (verbose>1)
    std::cout << "Number of nucleus types selected: " << selnuclei.size() << '\n';
  
  enum { O_ALL =0, O_UNIQUE, O_PAS };
  int outputsel=O_ALL;
  mode=ALL;

  if (!makespinsys) {
    std::cout << "A - Spin system of all spins\nS - Single spin spectrum (SIMPSON oriented)\nI - Single spin (pNMRsim individual spectra)\nT - Single spin (pNMRsim total spectrum only)\nX - text tables (e.g. for Excel import)\nG - gsim oriented output\n";
    mode=getoption(argc,argv,count,"Output type? ","ASITXG",mode);

    std::cout << "A - All matching spins\nU - Unique parameters\nP - Unique principal components\n";
    outputsel=getoption(argc,argv,count,"Combine entries? ","AUP",O_PAS);
    if (outputsel!=O_ALL)
      match_tolerance=getfloat(argc,argv,count,"Matching tolerance (e.g. 0, 1e-6)? ",0.0); 
  }
  const bool table_oriented=(mode==TABLES) || (mode==GSIM);
  const bool allspins=(mode==ALL) | table_oriented;
         
  //  assingle=getlogical(argc,argv,count,"Single spin per spectrum? ",true);
  try {
    MagresFile magres;
    magres.parse_from_file(fname);
    const bool makePDBlabels=(havePDBlabels && !makespinsys);
    if (makePDBlabels)
      std::cerr << "Warning: labels in new format magres will be ignored in favour of those specified in PDB\n";

    // need to go through magres object and "process" selected atoms
    char labelbuf[20];

    const size_t natoms=magres.num_atoms();
    List<int> revindex(natoms,-1);
    for (size_t i=0;i<natoms;i++) {
      const MagresAtom& curatom(magres.atoms(i));
      const size_t nuc=castep2nuc(curatom.species.c_str());
      if (isselected(nuc)) {
	usednuc.insert(nuc);
	const char* labelp=NULL;
	if (makePDBlabels)
	  labelp=get_PDB_label(curatom.species.c_str(),curatom.index);
	else {
	  labelp=curatom.label.c_str();
	  if (curatom.species==curatom.label) { //!< if no attempt to make label, create one
	    snprintf(labelbuf,sizeof(labelbuf),"%s:%zu",curatom.species.c_str(),curatom.index);
	    labelp=labelbuf;
	  } 
	}
	if (!labelp)
	  throw InternalError("Failed to label atom");
	const nucident cident(nuc,labelp,curatom.index);
	revindex(i)=atoms.size(); 
	atoms.push_back(nucdesc(cident,curatom.position));
      }
    }

    for (size_t i=0;i<magres.num_ms();i++) {
      const MagresMs& curms(magres.ms(i));
      const size_t atomindex=curms.atom - magres.atoms.vector(); //!< ugly hack to get index
      if (atomindex>=natoms)
	throw InternalError("Bad atom index");
      const int ind=revindex(atomindex);
      if (ind>=0) {
	nucdesc& curatom(atoms(size_t(ind)));
	add_interaction(curatom.csap,curms.sigma,false,-1.0);
      }
    }

    for (size_t i=0;i<magres.num_efg();i++) {
      const MagresEfg& curefg(magres.efg(i));
      const size_t atomindex=curefg.atom - magres.atoms.vector(); //!< ugly hack to get index
      if (atomindex>=natoms)
	throw InternalError("Bad atom index");
      const int ind=revindex(atomindex);
      if (ind>=0) {
	nucdesc& curatom(atoms(size_t(ind)));
	const double scale=quadscale(curatom);
	if (scale) {
	  const MagresEfg& curefg(magres.efg(i));
	  add_interaction(curatom.quadp,curefg.V,true,scale);
	}
      }
    }
  }
  catch (exception_t& exc) {
    std::cerr << "Parsing new format failed: " << exc << '\n';
    return 1;
  }
  catch (notmagres_exception_t&) {
    if (verbose)
      std::cerr << "Not a new format magres file - trying old format...\n";
    read_old_magres(fname);
  }
  
  if (atoms.empty()) {
    std::cerr << "Did not find any nuclei to save!\n";
    return 1;
  }

  if (verbose)
    std::cout << "Found " << atoms.size() << " matching nuclei\n";

  List<size_t> multiplicity;
  if (outputsel!=O_ALL) {
    if (outputsel==O_PAS)
      compress<ComparePAS>(atoms);
    else
      compress<CompareInteraction>(atoms);
    multiplicity.create(atoms.size());
    int uniqueval=-1;
    for (size_t i=0;i<atoms.size();i++) {
      const size_t mult=atoms(i).multiplicity;
      multiplicity(i)=mult;
      if (uniqueval<0)
	uniqueval=mult;
      else {
	if (uniqueval!=mult)
	  uniqueval=0;
      }
    }
    if (uniqueval>0)
      multiplicity.clear(); //!< unique multiplicity
    if (verbose)
      std::cout << "Nuclei left after merging: " << atoms.size() << '\n';
  }

  if (makespinsys) {

    if (verbose) {
      for (size_t i=0;i<atoms.size();i++) {
	const nucdesc& curatom(atoms(i));
	std::cout << curatom.ident << " [" << (i+1) << "] at " << curatom.position << '\n';
      }
    }

    List<nucdesc> newatoms;
    // Build spin system 
    char tmplabel[3]={'\0'};
    for (size_t i=0;i<molstruct.size();i++) {
      const PDBatom& curatom(molstruct(i));
      tmplabel[0]=*(curatom.type);
      const size_t nuc=castep2nuc(tmplabel);
      if (nuc==NULL_NUCLEUS) {
	std::cerr << "Unidentified nucleus type: " << tmplabel << " (ignoring)\n";
	continue;
      }
      const nucident cident(nuc,curatom.type,curatom.residue);
      if (haveCASTEPlabels) {
	if (!havePDBlabels) {
	  std::cerr << "Labels found in CASTEP file but no asymmetric units specified in PDB\n";
	  exit(2);
	}
      }
      else {
	if (havePDBlabels) {
	  std::cerr << "No labels in CASTEP file, but asymmetric units specified in PDB - unsupported combination\n";
	  exit(2);
	}
      }

      const store_t::iterator iter(std::find(atoms.begin(),atoms.end(),cident));
      if (iter==atoms.end())
	std::cerr << "Failed to find atom matching " << cident << " - ignoring.\n";
      else {
	newatoms.push_back(*iter);
	nucdesc& cursite(newatoms.back());
	cursite.position=curatom.coords; // overwrite position with coordinates from PDB spin system
      }
    }	
    if (newatoms.empty()) {
      std::cerr << "Failed to find any matching atoms!\n";
      exit(2);
    }
    atoms.swap(newatoms);
  }
  else {
  //  std::sort(atoms.begin(),atoms.end(),sortonlabel);
    sort_ip(atoms,sortonlabel);
  }

  if (verbose) {
    for (size_t i=0;i<atoms.size();i++) {
      const nucdesc& curatom(atoms(i));
      std::cout << curatom.ident << " [" << (i+1) << "] at " << curatom.position << '\n';
      //" in ASU " << curatom.aunit << '\n';
    }
  }

  bool showdipolar=makespinsys;

  char orientstr[MAX_LINE];
  getstring(argc,argv,count,"Bond vector and atom specifying orientation (e.g. H:6,H:5,C:1) [ENTER for none]? ",orientstr,sizeof(orientstr));  
  if (orientstr[0]) {

    List<size_t> orientindices;
    parse_index_list(orientindices,atoms,orientstr);
    if (verbose)
      std::cout << "Orientation indices: " << orientindices << '\n';

    if (orientindices.size()!=3) {
      std::cerr << "Couldn't parse orientation string as 3 comma-separated identifiers (<element><index>)\n";
      exit(1);
    }

    reorient_structure(atoms,orientindices,verbose);
    showdipolar=true;
  }

  //  atoms.front().position.y=1.10;

  char avstr[MAX_LINE];
  getstring(argc,argv,count,"List of spins to average principal values (e.g. H:6,H:5;H1:H2) [ENTER for none, ',' to separate, ';' to group]? ",avstr,sizeof(avstr));  

  List< Matrix<size_t> > savindices;

  if (avstr[0]) {
    parse_super_index_list(savindices,atoms,avstr);

    for (size_t grp=savindices.size();grp--;) {
      const Matrix<size_t>& avindices(savindices(grp));
      for (size_t r=avindices.rows();r--;) {
	averager qav,csav;
	if (verbose>1)
	  std::cout << "Averaging " << avindices.row(r) << '\n';
	size_t asu=0;
	bool mismatch=false;
	for (size_t i=0;i<avindices.cols();i++) {
	  const nucdesc& curatom(atoms(avindices(r,i)));
	  if (i) {
	    if (asu!=curatom.aunit)
	      mismatch=true;
	  }
	  else
	    asu=curatom.aunit;
	  if (curatom.csap.get())
	    csav.add(*(curatom.csap));
	  if (curatom.quadp.get())
	    qav.add(*(curatom.quadp));
	}
	csav.perform();
	qav.perform();
	if (mismatch) {
	  std::cerr << "Warning: averaging atoms that seem to belong in different asymmetric units:\n";
	  for (size_t i=0;i<avindices.cols();i++) {
	    const nucdesc& curatom(atoms(avindices(r,i)));	      
	    std::cerr << curatom.ident << '\n';
	  }
	}
      }
    }
  }

  enum { TENS_NONE =0, TENS_ORTEP, TENS_OTHER, TENS_DELETE };

  std::cout << "N - No tensor output\nO - ORTEP-oriented output\nG - gdis / other\nD - delete anisotropy information (from spin system output)\n";
  int tenstype=getoption(argc,argv,count,"Tensor output? ","NOGD",TENS_NONE);
  if ((tenstype==TENS_DELETE) && (mode==TABLES)) {
    std::cerr << "Delete anisotropy information makes no sense with table output (ignored)\n";
    tenstype=TENS_NONE;
  }
  
  if (tenstype!=TENS_NONE) {
    char buf[256];
    snprintf(buf,sizeof(buf),"%s_tensors.pdb",baseoutname);
    FILE* fp=fopen(buf,"wa");
    fputs("HEADER    DUMMYHEADERINFO\n",fp);
    fprintf(fp,"TITLE     Created by magres2pNMRsim from %s\n",fname);
    List<double> anisopar(6);
    bool usequad=true;
    for (size_t i=0;i<atoms.size();i++) {
      const nucdesc& curatom(atoms(i));
      const vector3& coords(curatom.position);    
      const bool isH=(curatom.ident.nuc==H_NUCLEUS);
      const bool isquad = (usequad && !isH);
      const interaction* intp((isquad ? curatom.quadp : curatom.csap).get());
      const char* source=curatom.label.c_str();
      char totlabel[8];
      const char* colonptr=strchr(source,':');
      if (colonptr) {
	strncpy(totlabel,colonptr+1,sizeof(totlabel)-1);
	if (isH && (tenstype==TENS_ORTEP))
	  *totlabel='F';
      }
      else {
	while (isdigit(*source))
	  ++source;
	size_t j=0;
	if (isH && (tenstype==TENS_ORTEP)) {
	  source++;
	  totlabel[j++]='F'; // fake F atom
	}
	for (;j<sizeof(totlabel);) {
	  if (*source!=':')
	    totlabel[j++]=*source;
	  if (*source++=='\0')
	    break;
	}
      }

      //snprintf(totlabel,sizeof(totlabel)-1,"%s%i",label,nuccount);
      write_PDB_atom(fp,i+1,totlabel,coords,curatom.aunit);
      //    fprintf(fp,"HETATM%5i %-4s UNK 0   1    %8.3f%8.3f%8.3f  1.00  0.00\n",i+1,totlabel,coords.x,coords.y,coords.z);
      if (!intp)
	continue;
      const rmatrix& A(intp->A);
      for (size_t j=3;j--;) 
	anisopar(j)=A(j,j);
      anisopar(size_t(3))=0.5*(A(size_t(0),size_t(1))+A(size_t(1),size_t(0)));
      anisopar(size_t(4))=0.5*(A(size_t(0),size_t(2))+A(size_t(2),size_t(0)));
      anisopar(size_t(5))=0.5*(A(size_t(1),size_t(2))+A(size_t(2),size_t(1)));
      double scale=isquad ? 1e-3 : (isH ? 100 : 10.0);
      if (anisopar(size_t(0))<0)
	anisopar*=-1.0;
      if (!isquad) {
	double iso=(anisopar(0)+anisopar(1)+anisopar(2))/3.0;
	const double aniso=fabs(intp->aniso);
	if (iso>1.5*aniso) {
	  const double subtract=iso-aniso;
	  anisopar(0)-=subtract;
	  anisopar(1)-=subtract;
	  anisopar(2)-=subtract;
	}
      }      
      //		siso=(s11+s22+s33)/3
      //		s11=s11-siso;
      //		s22=s22-siso;
      //		s33=s33-siso;
      fprintf(fp,"ANISOU%5" LCM_PRI_SIZE_T_MODIFIER "u %-4s UNK 0 %3" LCM_PRI_SIZE_T_MODIFIER "u  ",i+1,totlabel,curatom.aunit);
      for (size_t i=0;i<6;i++)
	fprintf(fp,"%7i",int(0.5+scale*anisopar(i)));
      fprintf(fp,"    \n");
    }
    fclose(fp);
  }

  if (table_oriented) {
    char buf[256];
    
    if (haveshifts) {
      if (mode==GSIM)
	snprintf(buf,sizeof(buf),"%s.gsimmagres",baseoutname);
      else
	snprintf(buf,sizeof(buf),"%s_shifts.txt",baseoutname);

      FILE* fpshifts=fopen(buf,"w");
      if (fpshifts==NULL) {
	perror("File open");
	exit(errno);
      }
      if (mode==TABLES) {
	fprintf(fpshifts,"Label\tIndex\tMultiplicity\tIsotropic shift (unref)/ ppm\tAnisotropy / ppm\tAsymmetry\talpha\tbeta\tgamma\n");
	for (size_t i=0;i<atoms.size();i++)
	  atoms(i).dumpshifttable(fpshifts);
      }
      else {
	for (size_t i=0;i<atoms.size();i++) {
	  const nucdesc& curatom(atoms(i));
	  const interaction* csap=curatom.csap.get();	 
	  if (csap) {
	    const char* elname=nuctoelement(curatom.nucleus());
	    // header necessary for gsim to recognise
	    fprintf(fpshifts,"=============\n");
	    fprintf(fpshifts,"Atom: %s  %zu\n",elname,curatom.index());
	    fprintf(fpshifts,"=============\n");
	    fprintf(fpshifts,"%s %zu   Isotropic: %g (ppm)\n",elname,curatom.index(),-(csap->iso)); //!< need to re-reverse sign to convert back shielding!
	  }
	}
      }
      fclose(fpshifts);
    }

    if (havequads && (mode==TABLES)) {
      snprintf(buf,sizeof(buf),"%s_quadrupoles.txt",baseoutname);
      FILE* fpquads=fopen(buf,"w");
      if (fpquads==NULL) {
	perror("File open");
	exit(errno);
      }
      fprintf(fpquads,"Label\tIndex\tMultiplicity\tCq / Hz\tAsymmetry\talpha\tbeta\tgamma\n");
      for (size_t i=0;i<atoms.size();i++)
	atoms(i).dumpquadtable(fpquads);      
      fclose(fpquads);
    }
    return 0;
  }

  FILE* fp=fopen(outname,"w");
  fprintf(fp,"#Created by magres2pNMRsim from %s\n",fname);
  if (pdbname[0])
    fprintf(fp,"#PDB label / selection file: %s\n",pdbname);
  if (havePDBlabels || haveCASTEPlabels) {
    fprintf(fp,"#Site labels:");
    for (size_t i=0;i<atoms.size();i++) {
      const char* label=atoms(i).label.c_str();
      //	fprintf(fp,"%c%s", i ? ' ' : '#', label ? label : "<unknown>");      
      fprintf(fp," %s", *label ? label : "<missing>");
    }
    fputc('\n',fp);
  } 

  fprintf(fp,"\tnuclei");
  const size_t uniquenuc=atoms.front().ident.nuc;
  if (allspins) {
    for (size_t i=0;i<atoms.size();i++)
      dumpnucleus(fp,atoms(i).ident.nuc);
  }
  else {
    if (usednuc.size()!=1) {
      std::cerr << "Must be only one type of nucleus active\n";
      return 1;
    }
    dumpnucleus(fp,uniquenuc);
  }
  fputc('\n',fp);

  const bool notens=(tenstype==TENS_DELETE);

  switch (mode) {
  case PNMRSIM_SINGLE: case PNMRSIM_SUM:
    if (atoms.front().csap.get()) {
      fputs("\tvariable rawshifts ",fp);
      dumplist(fp,atoms,Int2Type<nucdesc::CSA>(),Int2Type<interaction::ISO>()); fputc(' ',fp);
      fputc('\n',fp);
      fputs("\tshift 1 $rawshifts+$(SHIFT",fp);
      const char* llab=nuctolabel(uniquenuc);
      while (*llab)
	fputc(std::toupper(*llab++),fp);
      fputs("?0) ",fp);
      if (!notens) {
	dumplist(fp,atoms,Int2Type<nucdesc::CSA>(),Int2Type<interaction::ANISO>()); fputc(' ',fp);
	dumplist(fp,atoms,Int2Type<nucdesc::CSA>(),Int2Type<interaction::ASYM>()); fputc(' ',fp);
	dumplistangles(fp,atoms,Int2Type<nucdesc::CSA>());  
      }
      fputc('\n',fp);
    }
    if (!notens && atoms.front().quadp.get()) {
      fprintf(fp,"\tquadrupole 1 %i ",quadorder);
      dumplist(fp,atoms,Int2Type<nucdesc::QUAD>(),Int2Type<interaction::ANISO>()); fputc(' ',fp);
      dumplist(fp,atoms,Int2Type<nucdesc::QUAD>(),Int2Type<interaction::ASYM>()); fputc(' ',fp);
      dumplistangles(fp,atoms,Int2Type<nucdesc::QUAD>());      
      fputc('\n',fp);
    }    
    break;

  case ALL:
    for (size_t i=0;i<atoms.size();i++)
      atoms(i).dump(fp,i+1);
    break;
  case SIMPSON_SINGLE:
    for (size_t i=0;i<atoms.size();i++)
      atoms(i).dump(fp,1,i ? "#" : "");
    break;    
  default:
    throw InternalError("Unknown output mode");
  }

  if (tenstype==TENS_DELETE)
    showdipolar=false;
  else {
    if (!showdipolar)
      showdipolar=getlogical(argc,argv,count,"Calculate dipolar couplings? ",true);
  }

  if (showdipolar) {
    const size_t natoms=atoms.size();

    List<vector3> positions(natoms);
    List<double> gammavals(natoms);
    for (size_t i=natoms;i--;) {
      const nucdesc& curatom(atoms(i));
      positions(i)=curatom.position;
      gammavals(i)=spin(curatom.nucleus()).gamma();
    }

    TensorDipolarAverager dipav(positions,gammavals,savindices,verbose);

    for (size_t i=0;i<natoms;i++) {
      const double gammai=gammavals(i);
      const int seti=dipav.which_set(i);
      const bool inavseti=(seti>=0);
      const vector3 positioni=atoms(i).position;

      for (size_t j=i+1;j<natoms;j++) {
	const double gammaj=gammavals(j);
	const int setj=dipav.which_set(j);
	const bool inavsetj=(setj>=0);
	const vector3 positionj=atoms(j).position;

	spherical sphco;
	double scalefac=0.0;
	double d=0.0;
	double asym=0.0;
	
	Euler euler(0.0,0.0,0.0);
	bool haveeuler=false;

	if (inavseti || inavsetj) { //!< one or both in averaging set
	  const bool sameset=(seti==setj);
	  // std::cout << "In different sets? " << indiff << '\n';
	  if (!sameset) { //!< not in same averaging set
	    if (inavseti) {
	      if (inavsetj)
		throw Failed("Not implemented");
	      dipav.d_and_Euler(d,asym,euler,dipav.set(size_t(seti)),ExplicitList<1,size_t>(j));
	    }
	    else
	      dipav.d_and_Euler(d,asym,euler,ExplicitList<1,size_t>(i),dipav.set(size_t(setj)));

	    haveeuler=true;
	  }
	  else { //!< both spins are in the same set
	    //	    assert(whichtype(i)==whichtype(j));
	    if (seti!=setj)
	      throw InternalError("assertion failed: seti != setj");
	    const size_t useti(seti);
	    //	    const Matrix<size_t>& avindices(savindices(whichtype(i)));
	    const BaseList<size_t> curset(dipav.set(useti));
	    const vector3 ivec(positioni-positionj);
	    sphco=dipav.perpunitvector(useti);
	    sphco.r=ivec.length();
	    scalefac=-0.5;
	    if (verbose>1)
	      std::cout << "Setting distance (for spins within same averaging set) to " << sphco.r << " A and scaling factor to " << scalefac << '\n';
	  }
	}
	if (haveeuler)
	  fprintf(fp,"\tdipole %" LCM_PRI_SIZE_T_MODIFIER "u %" LCM_PRI_SIZE_T_MODIFIER "u %g %g %g %g %g\n",i+1,j+1,d,asym,euler.alpha*rad_to_deg,euler.beta*rad_to_deg,euler.gamma*rad_to_deg);
	else {
	  if (scalefac==0.0) {
	    const vector3 ivec(positioni-positionj);
	    sphco=spherical(ivec);
	  }
	  if (sphco.r==0.0)
	    throw Failed("dipolar_parameters: spins are co-incident");	
	  if (d==0.0) {
	    d=dipolar_coupling(gammai,gammaj,sphco.r*1e-10);
	    if (verbose>1) {
	      std::cout << "Calculated unaveraged coupling of " << d << " Hz";
	      if (scalefac)
		std::cout << ", scaled by " << scalefac << '\n';
	      else
		std::cout << '\n';	    
	    }
	    if (scalefac)
	      d*=scalefac;
	  }
	  else {
	    if (verbose>1)
	      std::cout << "Using supplied coupling of " << d << " Hz\n";
	  }
	  euler=vector_to_Euler(sphco);
	  if (verbose>1)
	    std::cout << "Euler angles " << euler << " determined from vector (spherical polar coords): " << sphco << '\n';
	  fprintf(fp,"\tdipole %" LCM_PRI_SIZE_T_MODIFIER "u %" LCM_PRI_SIZE_T_MODIFIER "u %g 0 %g %g\n",i+1,j+1,d,euler.beta*rad_to_deg,euler.gamma*rad_to_deg);
	}
      }
    }
  }

  if (multiplicity.size()) {
    if (mode==ALL) {
      fputs("#NB spectra should be adjusted by the following multiplicities to correctly reproduce the intensity pattern\n#",fp);
      dump(fp,multiplicity);
      fputc('\n',fp);
    }
    else {
      fputs("\tvariable scalefactors ",fp);      
      dump(fp,multiplicity);
      fputs("\n#NB ensure that scale $scalefactors is present in pulseq to correctly reproduce the intensity pattern\n",fp);
    }
    std::cout << "Warning: non-uniform multiplicity factors.  See output file for more information\n";
  }
  else
    fputs("\tvariable scalefactors 1 -ignore_unused\n",fp);

  if ((mode==PNMRSIM_SINGLE) || (mode==PNMRSIM_SUM))
    fprintf(fp,"\tsetenv SAVEFLAGS%s\n",(atoms.size()>1) ? " -sum" : "");
	    
  fclose(fp);
  return 0;
}

//  void cartesian_to_PAS_symmetric(double& iso, double& aniso, double& asym, Euler& PAS, const Matrix<double>& A, bool zeroiso =false)
//  {
//    Matrix<double> V;
//    ScratchList<double,3> eigs(3);
//    hermitian_eigensystem(V,eigs,A);
//    if (verbose>1)
//      std::cout << "Eigenvalues: " << eigs << '\n';
//    //  std::cout << "Eigenvectors\n" << V << '\n';  
//    ScratchList<double,3> tmpeigs(3);
//    const double isotropic=(eigs(0U)+eigs(1U)+eigs(2U))/3.0;
//    for (size_t i=3;i--;)
//      tmpeigs(i)=fabs(eigs(i)-isotropic);
//    //std::cout << "Relative eigenvalues: " << tmpeigs << '\n';
//    ScratchList<size_t,3> order(0U,0U,0U);
//    indexed_sort(order,tmpeigs);
//    std::swap(order(0U),order(1U)); //!< order is y,x,z !
//    if (verbose>1)
//      std::cout << "Order: " << order << '\n';
//    const Matrix<double> Vsorted(V(order,range()));
//    const ScratchList<double,3> sortedeigs(eigs(order));
//    //  std::cout << "Sorted eigenvalues: " << sortedeigs << '\n';
//    //  std::cout << "Sorted eigenvectors\n" << Vsorted << '\n';
//    //  interaction out;
//    cartesian_to_anisotropy(iso,aniso,asym,sortedeigs(0U),sortedeigs(1U),sortedeigs(2U));
//    if (asym<asym_zero_tol)
//      asym=0.0;
//    if (zeroiso)
//      iso=0.0;
//    PAS=matrix_to_euler(Vsorted);

//    if (verbose>1) {
//      if (!zeroiso)
//        std::cout << "Isotropic: " << iso << "  ";
//      std::cout << "Anisotropy: " << aniso << "  Asymmetry: " << asym << '\n';
//      std::cout << "Euler angles: " << PAS << '\n';
//    }
//  }

