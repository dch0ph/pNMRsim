#include "CrystalSystem.h"
#include "ttyio.h"
//#include "Propagation.h"
//#include "MetaPropagation.h"
//#include "simpsonio.h"
//#include "Histogram.h"
#include <fstream>
//#include <set>

using namespace libcmatrix;

#include "./geometry_utils.h"


using namespace std;

const double sw=100e3;
const size_t npts=1000;

const double deg_to_rad=M_PI/180;

const int verbose=1;

vector3 getvector(int argc,const char *argv[], int& count, const vector3& def =vector3(0,0,0))
{
  double x,y,z;
  x=getfloat(argc,argv,count," x: ", def.x);
  y=getfloat(argc,argv,count," y: ", def.y);
  z=getfloat(argc,argv,count," z: ", def.z);
  return vector3(x,y,z);
}

void apply_rotate_ip(BaseList<vector3> a, const rmatrix3& R)
{
  for (size_t i=a.size();i--;)
    a(i)=rotate(a(i),R);
}

void put(std::ostream& ostr, const BaseList<size_t>& ns, const char* sep =",", size_t offset =0)
{
  ostr << (offset+ns.front());
  for (size_t i=1;i<ns.size();i++)
    ostr << sep << (offset+ns(i));
}

//const double normlim=1e-12;

// template<class OpGen> void dosim(OpGen& opgen, const HamiltonianStore<double>& Hstore, const char* fname, const nuclei_spec& nuctype,bool normal =true)
// {
//   BlockedStaticHamiltonian<complex> Ham(opgen,Hstore);
//   if (verbose) {
//     cout << "Hamiltonian\n";
//     //    if (normal)
//     cout << Ham;
//     //else
//     //spy(cout,Ham());
//   }
//   const BlockedOperator Fp(opgen,operator_spec(nuctype,'+'));
//   if (verbose) {
//     cout << "sigma0\n";
//     // if (normal)
//     cout << Fp;
//     //else
//     // spy(cout,Fp);
//   }
//   const BlockedOperator dummy_detect;
//   const StaticSpectrumED obj(verbose);
//   MetaSpectrum<double> speciter(obj,Ham,Fp,dummy_detect,normal ? verbose : 1);
//   double freq;
//   double amp;
//   List<complex> spec(npts,0.0);
//   Histogram<complex> hist(spec,sw,-(sw+sw/npts)/2.0);
//   while (speciter(amp,freq)) {
//     if (norm(amp)>normlim) {
//       hist.add(amp,freq);
//       if (normal)
// 	std::cout << amp << " of " << freq << '\n';
//     }
//   }
//   write_simpson(fname,spec,sw,true);
// }

// void fillstore(HamiltonianStore<double>& Hstore, const spin_system& sys,const CrystalGeometry& geomspec, const Euler& powder)
// {
//   const size_t M=Hstore.nspins_cell();
//   const size_t total=Hstore.nspins();

//   for (size_t j=0; j<M; j++) {
//     for (size_t sk=j+1; sk<total; sk++) {
//       //      if (sys(j)!=sys(geomspec.spin_to_spin(sk))) {
//       const space_T A_PAS(geomspec.dipolar_tensor(sys,j,sk));
//       //	const double d=geomspec.dipolar_coupling(sys, j, sk);
//       const double d=real(rotate(A_PAS,2,0,powder));
//       Hstore.set_coupling(I_DIPOLE,j,sk,0.0,d);
// 	//   }
//     }
//   }
//   cout << "Interactions\n" << Hstore << '\n';
// }

// size_t getnuctype(char nuc)
// {
//   static size_t nuc13C(labeltonuc("13C"));
//   static size_t nuc1H(labeltonuc("1H"));

//   switch (nuc) {
//   case 'C':
//     return nuc13C;
//   case 'H':
//     return nuc1H;
//   }
//   cerr << "Unrecognised nucleus type: " << nuc << '\n';
//   exit(1);
// }

void fillsys(spin_system& sys, const BaseList<size_t>& nuclist, size_t hetnuc)
{
  size_t finish=sys.nspins();
  if (hetnuc!=NULL_NUCLEUS)
    finish--;
  for (size_t base=0;base<finish;base+=nuclist.size()) {
    for (size_t j=nuclist.size();j--;)
      sys.isotope(base+j,spin(nuclist(j)));
  }
  if (hetnuc!=NULL_NUCLEUS)
      sys.isotope(finish,spin(hetnuc));
}

template<typename T> List<T> getnuclist(const char* spinslist)
{
  char tmp[256];
  char* strp=strncpy(tmp,spinslist,sizeof(tmp));
  char *ptr;
  List<T> nuclist;
  while ((ptr=strtok(strp,","))!=NULL) {
    nuclist.push_back(labeltonuc(ptr));
    strp=NULL;
  }
//   const size_t M=strlen(spinslist);

//   for (size_t j=0;j<M;j++)
//     nuclist.push_back(getnuctype(spinslist[j]));
  return nuclist;
}

namespace {

  void write_atom(FILE* fp, size_t nuc, const vector3& coords)
  {
    static size_t serial=0;
    typedef std::map<char,size_t> count_t;
    static count_t countmap;
    serial++;
    const count_t::iterator iter(countmap.find(nuc));
    size_t nuccount=1;
    if (iter==countmap.end())
      countmap[nuc]=1;
    else
      nuccount=(++(iter->second));

    const char* label=nuctolabel(nuc);
    //    bool warn=false;
    //while (isdigit(*label)) {
    //  label++;
    //  if (strlen(label)>2)
    //	warn=true;
    //}
    //  if (warn)
    //   cerr << "Warning: nucleus name > 2 characters long detected.  PDB file name be badly formatted.\n";
    char totlabel[8];
    snprintf(totlabel,sizeof(totlabel-1),"%s%" LCM_PRI_SIZE_T_MODIFIER "u",label,nuccount);
    write_PDB_atom(fp,serial,totlabel,coords);
  }
}

void write_PDB(const char* fname, const CrystalGeometry& crysgeom, const BaseList<size_t>& nuclist, const size_t hetnuc, const vector3& hetpos)
{
  FILE* fp=fopen(fname,"w");
  if (fp==NULL)
    throw Failed("write_PDB: file open");
  
  write_PDB_comment(fp,"created by createdip");

  const size_t M=crysgeom.nspins_cell();
  for (size_t i=0;i<crysgeom.nspins();i++)
    write_atom(fp,nuclist[i % M],crysgeom(i));
  if (hetnuc!=NULL_NUCLEUS)
    write_atom(fp,hetnuc,hetpos);

  fclose(fp);
}

void parselist(List<size_t>& dest, char* mapstr)
{
  char* p=mapstr;
  char* tok;
  while ((tok=strtok(p,","))) {
    p=NULL;
    char* endptr;
    const long val=strtol(tok,&endptr,10);
    if (endptr[0] || (val<1)) {
      fprintf(stderr,"Failed to parse at positive integer: %s\n",tok);
      exit(1);
    }
    dest.push_back(size_t(val-1));
  }
}
  
int main(int argc,const char *argv[])
{
  int count=1;
  const char axis[]="abc";

#ifndef LCM_USE_VECTOR_TO_EULER
  std::cerr << "Warning: not using vector_to_Euler function for determining dipolar PAS - may give Euler angles that are inconsistent with other NMR interactions\n";
#endif
  
  try {
    //  getstring(argc,argv,count,"Isotope name (e.g. 1H) ? ", nuctype,sizeof(nuctype));
    char spinslist[200];
    cout << "Enter number and type of spins per cell (e.g. 13C,1H,1H)\nFirst spin is detected.\n";
    getstring(argc,argv,count,"Number and type of spins? ",spinslist,sizeof(spinslist));
    List<size_t> nuclist(getnuclist<size_t>(spinslist));
    const int M=nuclist.size();

    const char* detectnuc(nuctolabel(nuclist.front()));

    Permutation perm;
//     char permlist[10];
//     if (M>1) {
//       getstring(argc,argv,count,"Permutation vector? ",permlist,sizeof(permlist));
//       if (permlist[0])
// 	perm=Permutation(permlist);
//     }
    const int dims=getint(argc,argv,count,"Number of dimensions? ",0,3,1);
    List<size_t> ncells(dims);
    List<vector3> crys_vecs(dims,vector3(0,0,0));
    int total=M;
    for (size_t dim=0;dim<dims;dim++) {
      cout << "Axis " << axis[dim] << ":\n";
      const int N=getint(argc,argv,count,"Number of cells? ",2);
      ncells(dim)=N;
      total*=N;
      if (N>1) {
	cout << "Translational vector:\n";
	crys_vecs(dim)=getvector(argc,argv,count);
      }
    }
    if (total<2)
      throw Failed("Must have at least 2 spins!");

    List<vector3> baseco(M);
    baseco(0U)=vector3(0,0,0);
    if (M>1) {
      cout << "Coordinates of spin 1 are fixed at " << baseco(0U) << '\n';
      for (size_t i=1; i<M; i++) {
	cout << "Spin no " << (i+1) << '\n';
	baseco(i)=getvector(argc,argv,count);
      }
    }

    char hetname[20];
    bool baseonly=false;
    getstring(argc,argv,count,"Heteronucleus type e.g. 13C, null if none? ",hetname,sizeof(hetname));
    size_t hetnuc=NULL_NUCLEUS;
    if (hetname[0]) {
      baseonly=getlogical(argc,argv,count,"Only include couplings to core unit cell? ",true);
      hetnuc=labeltonuc(hetname);
    }
      
    const bool addhet=(hetnuc!=NULL_NUCLEUS);
	      
    vector3 hetpos(0,0,0);
    if (addhet)
      hetpos=getvector(argc,argv,count);

    if (!addhet && getlogical(argc,argv,count,"Rotate crystal? ",false)) {
      Euler MF_to_CF;
      MF_to_CF.alpha=getfloat(argc,argv,count,"alpha (degrees)? ")*deg_to_rad;
      MF_to_CF.beta=getfloat(argc,argv,count,"beta (degrees)? ")*deg_to_rad;
      MF_to_CF.gamma=getfloat(argc,argv,count,"gamma (degrees)? ",0.0)*deg_to_rad;
      rmatrix3 rotmat;
      rotation_matrix(rotmat,MF_to_CF);
      apply_rotate_ip(baseco,rotmat);
      apply_rotate_ip(crys_vecs,rotmat);
      if (getlogical(argc,argv,count,"Rotate in rotor frame? ")) {
	rotation_matrix(rotmat,Euler(0,MAGIC_ANGLE,0));
	apply_rotate_ip(baseco,rotmat);
	apply_rotate_ip(crys_vecs,rotmat);
      }
    }

    spin_system sys(addhet ? M+1 : M,detectnuc);
    fillsys(sys,nuclist,hetnuc);

    CrystalGeometry geomspec(ncells,crys_vecs,baseco);//,perm);

    char fname[260],fnamebase[256];
    getstring(argc,argv,count,"Output file [without .inc]? ",fnamebase,sizeof(fnamebase));

    //if (fnamebase[0]) {

      double drss=0;
      double drss_theta=0;
      for (size_t sk=1; sk<total; sk++) {
	const double d=geomspec.dipolar_coupling_unscaled(sys,0U,sk);
	drss+=d*d;
	const double d_theta=geomspec.dipolar_coupling(sys,0U,sk);
	drss_theta+=d_theta*d_theta;
      }
      drss=sqrt(drss);
      drss_theta=sqrt(drss_theta);
      cout << "d_rss for spin 1:  without angular dependence = " << drss << " Hz \t with angular dependence = " << drss_theta << " Hz\n";

      sprintf(fname,"%s.pdb",fnamebase);      
      write_PDB(fname,geomspec,nuclist,hetnuc,hetpos);

      const double rss=getfloat(argc,argv,count,"Target rss coupling (kHz)? [0 if none] ",0.0)*1e3;

      const bool exploitp = addhet ? false : getlogical(argc,argv,count,"Exploit periodicity? ",true);
      List<size_t> mapping;
      if (!exploitp) {
	char mapstr[256];
	getstring(argc,argv,count,"Mapping (ENTER for none)? ",mapstr,sizeof(mapstr),"");
	if (mapstr[0])
	  parselist(mapping,mapstr);
      }
      if (mapping.empty()) {
	mapping.create(total);
	for (size_t i=total;i--;)
	  mapping(i)=i;
      }

      const double scale= rss ? rss/drss : 1.0;

      sprintf(fname,"%s.inc",fnamebase);      
      ofstream outfile;
      outfile.open(fname, ios_base::trunc | ios_base::out );    

      outfile << "#created by " << argv[0] << " using " << ncells << " cell(s) of " << M << " spin(s) with coordinates\n";
      for (size_t i=0;i<M;i++)
	outfile << "#Spin " << (i+1) << ": " << baseco(i) << '\n';
      for (size_t i=0;i<dims;i++)
	outfile << "#Translation vector " << axis[i] << ": " << crys_vecs(i) << '\n';
      if (addhet)
	outfile << "#Heteronucleus at " << hetpos << '\n';
      
      if (exploitp && (dims || !perm.empty())) {
	outfile << "cells ";
	if (dims)
	  put(outfile,ncells," ");
	if (!perm.empty()) {
	  for (size_t i=0;i<perm.size();i++)
	    outfile << ((i==0) ? " {" : ",") << (perm(i)+1);
	  outfile << '}';
	}
	outfile << '\n';
      }
      outfile << "nuclei";
      const size_t maxn = exploitp ? M : mapping.size();
      const size_t maxtot = exploitp ? total : maxn;
      for (size_t q=0; q<maxn; q++)
	outfile << ' ' << nuctolabel(nuclist(mapping(q) % M));
      if (addhet)
	outfile << ' ' << hetname;
      outfile << '\n';
      
      List<size_t> ns;
      size_t k;
      Euler PAS;
      double d;

      for (size_t usej=0;usej<maxn;usej++) {
	const size_t j=mapping(usej);
	const size_t base = M*(j/M);
	const size_t bj = j % M;

	for (size_t usesk=usej+1; usesk<maxtot; usesk++) {
	  const size_t sk=mapping(usesk);
	  geomspec.dipolar_parameters(d, PAS, sys, bj, sk-base);
	  //cout << (j+1) << "," << (sk+1) << '\n';
	  geomspec.reverse(ns,k,sk);
	  if ((spinslist[j]!='C') || (spinslist[k]!='C')) { //don't output C-C couplings
	    outfile << "dipole " << usej+1 << ' ';
	    if (exploitp) {
	      if (dims) {
		put(outfile,ns,",",1);
		outfile << ',';
	      }
	      outfile << (k+1);
	    }
	    else
	      outfile << (usesk+1);
	    outfile << ' ' << (d*scale) << ' ' << (PAS.alpha/deg_to_rad) << ' ' << (PAS.beta/deg_to_rad) << ' ' << (PAS.gamma/deg_to_rad) << '\n';
	  }
	}
      }
      if (addhet) {
	const double hetgamma=spin(hetnuc).gamma();
	const size_t maxj = baseonly ? M : total;
	for (size_t usej=0;usej<maxj;usej++) {
	  const size_t bj = usej % M;
	  const vector3 pos(geomspec(usej));
	  const spherical sphco(pos-hetpos);
	  const double d=dipolar_coupling(sys(bj).gamma(),hetgamma,sphco.r*1e-10);
#ifdef LCM_USE_VECTOR_TO_EULER
	  Euler PAS(vector_to_Euler(sphco));
#else
	  Euler PAS(0,sphco.theta,sphco.phi);
#endif
	  outfile << "dipole " << (usej+1) << ' ' << (total+1) << ' ' << d << ' ' << (PAS.alpha/deg_to_rad) << ' ' << (PAS.beta/deg_to_rad) << ' ' << (PAS.gamma/deg_to_rad) << '\n';
	}
      }

      outfile.close();
      //    return 0;
  //   }

//     CrystalStructure cstruct(ncells,perm);
//     cout << "CrystalStructure: " << cstruct << '\n';

//     const Euler powder(M_PI/4,M_PI/3,0.0); //non-special angle
 
//     int flags=0;
//     getstring(argc,argv,count,"Nuclei to block by? ",spinslist,sizeof(spinslist));
//     const List<nuclei_spec> blockingnuc(getnuclist<nuclei_spec>(spinslist));
//     if (!blockingnuc.empty()) {
//       if (getlogical(argc,argv,count,"Exploit mz symmetry? "))
// 	flags|=MetaFlags::UseMzSymmetry;
//     }
//     if (getlogical(argc,argv,count,"Exploit k symmetry? "))
//       flags|=MetaFlags::UseEigSymmetry;    

//     const bool check1d=(dims>1) && perm.empty();
//     enum { CHECK_NONE=0, CHECK_1D, CHECK_EXACT };
//     cout << "N - none\n1 - 1D\nE - exact\n";
//     const size_t checkmethod=getoption(argc,argv,count,"Check method? ","N1E",check1d ? CHECK_1D : CHECK_EXACT);

//     if (checkmethod==CHECK_1D) { //check against 1D symmetrisation

//       const size_t inner= (dims==2) ? ncells(1U) : ncells(1U)*ncells(2U);
//       const size_t lM=M*inner;
//       List<vector3> lbaseco(lM);
//       for (size_t j=lM;j--;)
// 	lbaseco(j)=geomspec(j);
//       CrystalGeometry lgeomspec(ncells.front(),crys_vecs.front(),lbaseco);

//       HamiltonianStore<double> Hstore(lM,ncells.front());
//       spinhalf_system lsys(lM,detectnuc);
//       fillsys(lsys,nuclist);
//       fillstore(Hstore,lsys,lgeomspec,powder);

//       CrystalStructure cstruct1D(ncells(0U));
//       cout << "CrystalStructure 1D: " << cstruct1D << '\n';
      
//       CrystalOpGenerator crystal_opgen(lsys,cstruct1D,blockingnuc,flags,verbose);
//       sprintf(fname,"out%iD.spe",total);
//       dosim(crystal_opgen,Hstore,fname,detectnuc);
//     }

//     HamiltonianStore<double> Hstore(M,cstruct.ncells());
//     fillstore(Hstore,sys,geomspec,powder);
    
//     if (!Hstore.verify(cout,cstruct,1e-4))
//       return 1;
//     cout << "Periodic structure verified\n";

//     if (checkmethod==CHECK_EXACT) {
//       SpinOpGenerator simple_opgen(sys,cstruct,blockingnuc,flags,verbose);
//       sprintf(fname,"out%is.spe",total);
//       dosim(simple_opgen,Hstore,fname,detectnuc,total<6);
//     }

//     CrystalOpGenerator crystal_opgen(sys,cstruct,blockingnuc,flags,verbose);
//     sprintf(fname,"out%ic.spe",total);
//     dosim(crystal_opgen,Hstore,fname,detectnuc);

  } catch (MatrixException& exc) {
    cerr << exc;
    return 1;
  }
  return 0;
}
