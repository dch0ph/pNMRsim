#include <set>

void transformation_matrix(rmatrix3& dest, const vector3& x, const vector3& y, const vector3& z);
  
void write_PDB_comment(FILE* fp, const char* str);
void write_PDB_atom(FILE* fp, size_t serial, const char label[8], const vector3& coords, size_t seq =1);

template<class T> void reorient_structure(BaseList<T> atoms, const BaseList<size_t>& orientindices, int verbose =0)
{
  vector3 ivec(atoms(orientindices.front()).position - atoms(orientindices(size_t(1))).position);
  const double ilength=ivec.length();
  ivec/=ilength;
  const vector3 othervec(atoms(orientindices.front()).position - atoms(orientindices(size_t(2))).position);
  vector3 perp=othervec-dot(ivec,othervec)*ivec;
  perp/=perp.length();
  const double perpcheck=dot(perp,ivec);
  const vector3 third=cross(ivec,perp);
  rmatrix3 rotmat;
  //    transformation_matrix(rotmat,perp,third,ivec); //!< rotate so that ivec is along z, 13C is in xz plane
  transformation_matrix(rotmat,third,ivec,perp); //!< rotate so that ivec is along y, 13C is in yz
  if (verbose) {
    std::cout << "Normalised internuclear vector: " << ivec << " (length: " << ilength << " A)\n";
    std::cout << "Perpendicular vector: " << perp << "  (dot product: " << perpcheck << ")\n";
    std::cout << "Rotation matrix\n" << rotmat;
  }
  Matrix<double> fullrotmat;
  rmatrix3_to_rmatrix(fullrotmat,rotmat);
  
  for (size_t i=atoms.size();i--;)
    atoms(i).rotate(rotmat,fullrotmat);
  
  const vector3 ovec=0.5*(atoms(orientindices.front()).position + atoms(orientindices(size_t(1))).position);
  for (size_t i=atoms.size();i--;)
    atoms(i).position-=ovec;
  
  vector3 nivec(atoms(orientindices.front()).position - atoms(orientindices(size_t(1))).position);
  nivec/=nivec.length();
  const vector3 nothervec(atoms(orientindices.front()).position - atoms(orientindices(size_t(2))).position);
  vector3 nperp=nothervec-dot(nivec,nothervec)*nivec;
  nperp/=nperp.length();
  
  if (verbose) {
    std::cout << "Rotated normalised internuclear vector: " << nivec << '\n';
    std::cout << "Rotated perpendicular vector: " << nperp << '\n';
  }
}

bool testsame(const Matrix<double>& An, const Matrix<double>& A, const char* head, int verbose);
//bool testsame(const Euler&, const Euler&, const char* head, int verbose);

size_t default_nucleus(const char*);

typedef std::set<size_t> set_t;

class BaseDipolarAverager {
public:
  BaseDipolarAverager(const BaseList<vector3>&, double, const ListList<size_t>&, int =0);
  BaseDipolarAverager(const BaseList<vector3>&, const BaseList<double>&, const ListList<size_t>& groupsv, int =0);

  BaseDipolarAverager(const BaseList<vector3>&, double, const BaseList< Matrix<size_t> >&, int =0);
  BaseDipolarAverager(const BaseList<vector3>&, const BaseList<double>&, const BaseList< Matrix<size_t> >& groupsv, int =0);

  bool in_set(size_t i) const {
    validate_index(i);
    return (which_set_(i)>=0);
  }

  int which_set(size_t i) const {
    validate_index(i);
    return which_set_(i);
  }

  const BaseList<size_t> set(size_t i) const {
    if (i>groups_.size())
      throw BadIndex("DipolarAverager: set index out of range",i,groups_.size());
    return groups_(i);
  }

  void validate_index(size_t i) const {
    if (i>=natoms_)
      throw BadIndex("DipolarAverager: atom index out of range",i,natoms_);
  }

  void validate_set_index(size_t i) const {
    if (i>=groups_.size())
      throw BadIndex("DipolarAverager: set index out of range",i,groups_.size());
  }

  bool check_nooverlap(const BaseList<size_t>&) const;

  virtual double dss(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const =0;
  virtual void d_and_Euler(double&, double&, Euler&, const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const =0;

  int validate_set(const BaseList<size_t>&) const;
  std::pair<int,int> validate_sets(const BaseList<size_t>&, const BaseList<size_t>&) const;

  double gamma(size_t i) const {
    validate_index(i);
    return unsafe_gamma(i);
  }

protected:
  size_t natoms_;
  List<vector3> positions_;
  ListList<size_t> groups_;
  List<double> gammavals_;
  List<int> which_set_;
  int verbose_;

 protected:
  double unsafe_gamma(size_t i) const {
    const double v=gammavals_(i);
    if (!v)
      throw Failed("DipolarAverager: gamma of zro detected. Unexpected index reqest");
    return v;
  }

private:
  void makegroups(const BaseList< Matrix<size_t> >&);
  void initialise();
};

class TensorDipolarAverager : public BaseDipolarAverager 
{
public:
  TensorDipolarAverager(const BaseList<vector3>& positionsv, double gammaval, const ListList<size_t>& groupsv, int verbosev =0)
    : BaseDipolarAverager(positionsv, gammaval, groupsv, verbosev)
  { initialise_tensor(); }

  TensorDipolarAverager(const BaseList<vector3>& positionsv, const BaseList<double>& gammavalsv, const ListList<size_t>& groupsv, int verbosev =0)
    : BaseDipolarAverager(positionsv, gammavalsv, groupsv, verbosev)
  { initialise_tensor(); }

  TensorDipolarAverager(const BaseList<vector3>& positionsv, double gammaval, const BaseList< Matrix<size_t> >& groupsv, int verbosev =0)
    : BaseDipolarAverager(positionsv, gammaval, groupsv, verbosev)
  { initialise_tensor(); }

  TensorDipolarAverager(const BaseList<vector3>& positionsv, const BaseList<double>& gammavalsv, const BaseList< Matrix<size_t> >& groupsv, int verbosev =0)
    : BaseDipolarAverager(positionsv, gammavalsv, groupsv, verbosev)
  { initialise_tensor(); }

    double dss(const BaseList<size_t>& whichi, const BaseList<size_t>& whichj) const;

  void d_and_Euler(double&, double&, Euler&, const BaseList<size_t>&, const BaseList<size_t>&) const;
  const spherical& perpunitvector(const size_t setv) const {
    validate_set_index(setv);
    return perpunitvector_(setv);
  }

  private:
  List<spherical> perpunitvector_; //direction vector

  void initialise_tensor();
};
