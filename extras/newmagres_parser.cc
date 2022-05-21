
#include "magres.h"

using namespace libcmatrix;

namespace {
  long parse_integer(const char* p)
  {
    char* tail;
    long val=strtol(p,&tail,10);
    if (*tail!='\0') {
      char errm[MAGRES_ERRM_MAXLEN];
      sprintf(errm,"failed to parse %s as integer", p);
      throw MagRes::exception_t(errm);
    }
    return val;
  }
  
  double parse_double(const char* p)
  {
    char* tail;
    double val=strtod(p,&tail);
    if (*tail!='\0') {
      char errm[MAGRES_ERRM_MAXLEN];
      sprintf(errm,"failed to parse %s as floating point value", p);
      throw MagRes::exception_t(errm);
    }
    return val;
  }

  // Note this duplicates functionality of FileHandle in DynamicList.h
  class FileOpenGuard {
  public:
    FileOpenGuard(const char*);
    ~FileOpenGuard() { fclose(fp); }
    FILE* operator()() { return fp; }
  private:
    FILE* fp;
  };
  
  FileOpenGuard::FileOpenGuard(const char* path)
  {
    fp = fopen(path, "r");          
    char errm[MAGRES_ERRM_MAXLEN];
    if (fp==NULL) {
      sprintf(errm,"Failed to open file for reading: %s",path);
      throw MagRes::exception_t(errm);
    }
  }
 
}

namespace MagRes {

MagresLine::MagresLine(char* p)
{
  char* token;
  while ((token=strtok(p,token_separator))) {
    p=NULL;
    cols.push_back(token);
  }
}

void MagresLine::verify_columns(size_t cols, const char* label) const
{
  if (num_cols() != cols) {
    char errm[MAGRES_ERRM_MAXLEN];
    sprintf(errm,"%s: Wrong number of columns, %zu", label, num_cols());
    throw exception_t(errm);
  }
}

// [atoms] block parser
MagresAtom::MagresAtom(const MagresLine& line)
{
  line.verify_columns(7,"atom");
  species=line.column(1);
  label=line.column(2);
  index = parse_integer(line.column(3));
  position.x = parse_double(line.column(4));
  position.y = parse_double(line.column(5));
  position.z = parse_double(line.column(6));
}

MagresSymmetry::MagresSymmetry(const MagresLine& line)
{
  line.verify_columns(2,"symmetry");
  symmetry_string=line.column(1);
}

void MagresLine::read_tensor(tensor_t& dest, size_t base) const
{
  dest[0][0] = parse_double(column(base));
  dest[0][1] = parse_double(column(base+1));
  dest[0][2] = parse_double(column(base+2));

  dest[1][0] = parse_double(column(base+3));
  dest[1][1] = parse_double(column(base+4));
  dest[1][2] = parse_double(column(base+5));

  dest[2][0] = parse_double(column(base+6));
  dest[2][1] = parse_double(column(base+7));
  dest[2][2] = parse_double(column(base+8));  
}

MagresLattice::MagresLattice(const MagresLine& line)
{
  line.verify_columns(10,"lattice");
  line.read_tensor(lattice,1);
}

void MagresFile::parse_atom(const MagresLine& line)
{
  const char* id=line.column(0);
  if (strcmp(id, "atom") == 0)
    atoms.push_back(MagresAtom(line));
  else if (strcmp(id, "lattice") == 0)
    latticep= new MagresLattice(line);
  else if (strcmp(id, "symmetry") == 0)
    symmetries.push_back(line);
  //ignore everything else, including units
}

const MagresAtom* MagresFile::find_atom(const atomid_t& searchatom) const
{
  LIST_T<MagresAtom>::const_iterator iter=std::find(atoms.begin(),atoms.end(),searchatom);
  if (iter!=atoms.end())
    return &(*iter);

  fprintf(stderr, "Could not find atom %s %zu", searchatom.species.c_str(), searchatom.index);
  return NULL;
}

const MagresAtom* MagresFile::parse_find_atom(const MagresLine& line,size_t base) const
{
  const atomid_t searchatom( line.column(base), parse_integer(line.column(base+1)));
  return find_atom(searchatom);
}

MagresIsc::MagresIsc(const MagresLine& line, MagresFile& magres_file) 
{
  line.verify_columns(14,"isc");

  atom1=magres_file.parse_find_atom(line,1);
  atom2=magres_file.parse_find_atom(line,3);

  line.read_tensor(K,5);
}

MagresEfg::MagresEfg(const MagresLine& line, MagresFile& magres_file)
{
  line.verify_columns(12,"efg");
  atom = magres_file.parse_find_atom(line,1);
  line.read_tensor(V,3);
}

MagresMs::MagresMs(const MagresLine& line, MagresFile& magres_file)
{
  line.verify_columns(12,"ms");
  atom = magres_file.parse_find_atom(line,1);
  line.read_tensor(sigma,3);
}

void MagresFile::parse_magres(const MagresLine& line)
{
  const char* id=line.column(0);
  
  if (strcmp(id, "isc") == 0)
    isc.push_back(MagresIsc(line,*this));
  else if (strcmp(id, "efg") == 0)
    efg.push_back(MagresEfg(line,*this));
  else if (strcmp(id, "ms") == 0)
    ms.push_back(MagresMs(line,*this));
}

char* MagresFile::windforward(char* p)
{
  while( (*p != line_separator) && (*p != 0))
    ++p;
  return p;
}

void MagresFile::parse_lines(char *block, block_t blocktype)
{
  char* lastp=block;

  for (char* p=block;;p++) {
    if ( (*p == line_separator) || (*p == '\0') || (*p == '#') ) {
      if (p>lastp) {
	char storep=*p;	
	*p='\0'; // terminate
        MagresLine line(lastp);
	if (line.num_cols()) {
	  switch (blocktype) {
	  case ATOM:
	    parse_atom(line);
	    break;
	  case MAGRES:
	    parse_magres(line);
	    break;
	  }
	}
	*p=storep; //!< restore
      }

      // If this is a comment, wind forward p to EOL or EOF
      if (*p == '#')
	p=windforward(p);

      lastp = p + 1;
      if (*p == '\0') 
	break;
    }
  }
}

// Parse a magres file from a string
void MagresFile::parse_from_string(char* file)
{
  bool start_block = false; 
  char* block_name = NULL;
  char* block_data_end = NULL;
  char* block_data_start = NULL;

  unsigned int major,minor;
  if (sscanf(file,"#$magres-abinitio-v%u.%u",&major,&minor)!=2)
    throw notmagres_exception_t("ERROR: Not a new-style magres file");
  
  if (major>major_version_limit) {
    char errm[MAGRES_ERRM_MAXLEN];
    sprintf(errm,"ERROR: magres version (%d.%d) exceeds limit (%d)",major,minor,major_version_limit);
    throw exception_t(errm);
  }

  for (char* p= file;*p;p++) {
    switch (*p) {
      case '#': // comment, skip rest of line
        p=windforward(p);
        break;
	
    case '[': // start of block tag
      if(*(p+1) == '/') { // end block tag
	if(!start_block)
	  throw exception_t("ERROR: End block tag with no start block");
	
	start_block = false;
	block_data_end = p - 2;
      } else { // start block tag
	if (start_block)
	  throw exception_t("ERROR: Start block tag inside block");
	
	block_name = p + 1;
	start_block = true;
      }
      break;
      
    case ']': // end of block tag
      if (start_block) {
	*p='\0';
	block_data_start = p + 1;
      } 
      else {
	if ((block_data_end==NULL) || (block_data_start==NULL))
	  throw exception_t("ERROR: block data start/end unset");
	block_data_end[1]='\0'; //!< terminate data block
	if (strcmp(block_name, "atoms") == 0)
	  parse_lines(block_data_start,ATOM);
	else if (strcmp(block_name, "magres") == 0)
	  parse_lines(block_data_start,MAGRES);
	block_name=NULL;
      }
    }
  }
  
  if (block_name != NULL) {
    char errm[MAGRES_ERRM_MAXLEN];
    sprintf(errm,"Unterminated block %s", block_name);
    throw exception_t(errm);
  }
}

 void MagresFile::parse_from_file(const char* path)
   {
     FileOpenGuard FP(path);
     FILE* fp=FP();
     if (fseek(fp, 0L, SEEK_END) == 0) {
       const int buffer_size = ftell(fp);       
       if (buffer_size == -1)
	 throw exception_t("Error reading file size");
       
       std::string contents;     
       contents.resize(buffer_size+1);
       std::rewind(fp);
       const size_t new_length = std::fread(&contents[0], 1, buffer_size, fp);
       if(new_length == 0)
	 throw exception_t("Error reading file");
       
       contents[new_length+1] = '\0'; /* Just to be safe. */
       char* asraw=const_cast<char*>(contents.c_str()); //marginally dodgy - need to assume that parse_from_string doesn't over-run
       parse_from_string(asraw);
     }
     else
       throw exception_t("Error finding file end");
   }

std::ostream& print_tensor(const tensor_t& A, std::ostream& ostr)
   {     
     ostr << "  " << A[0][0] << ' ' << A[0][1] << ' ' << A[0][2] << '\n';
     ostr << "  " << A[1][0] << ' ' << A[1][1] << ' ' << A[1][2] << '\n';
     ostr << "  " << A[2][0] << ' ' << A[2][1] << ' ' << A[2][2] << '\n';
     return ostr;
   }
 
 std::ostream& operator<< (std::ostream& ostr, const MagresAtom& atom)
   {
     return ostr << "  " << atom.index << ' ' << atom.species <<  ' ' << atom.label <<  ' ' << atom.position << '\n';
   }

 std::ostream& operator<< (std::ostream& ostr, const MagresFile& magres_file)
   {
     ostr << magres_file.num_atoms() << " atoms\n";
     ostr << magres_file.num_symmetries() << " symmetries\n";
     
     if (magres_file.latticep) {
       ostr <<  "Lattice:\n";
       print_tensor((magres_file.latticep)->lattice,ostr);
     }
     else
       ostr << "No lattice\n";

     ostr << magres_file.num_isc() << " J-couplings\n";
     ostr << magres_file.num_efg() << " EFG tensors\n";
     ostr << magres_file.num_ms() << " MS tensors\n";

     ostr << "Atoms:\n";
     for (size_t i = 0; i<magres_file.num_atoms(); ++i)
       ostr << magres_file.atoms(i);
     
     printf("Symmetries:\n");
     for (size_t i = 0; i<magres_file.num_symmetries(); ++i)
       ostr << magres_file.symmetries(i).symmetry_string << '\n';

     for (size_t i = 0; i<magres_file.num_isc(); ++i) {
       const MagresIsc& isc(magres_file.isc(i));
       const double K_iso = (isc.K[0][0] + isc.K[1][1] + isc.K[2][2])/3.0;
       ostr << "ISC: ";
       (isc.atom1)->printname(ostr) << " --> ";
       (isc.atom2)->printname(ostr) << " = " << K_iso << '\n';
     }

     for (size_t i = 0; i<magres_file.num_ms(); ++i) {
       const MagresMs& ms(magres_file.ms(i));       
       const double ms_iso = (ms.sigma[0][0] + ms.sigma[1][1] + ms.sigma[2][2])/3.0;
       ostr << "MS: ";
       (ms.atom)->printname(ostr) << " = " << ms_iso << '\n';
     }

     return ostr;
   }

} //namespace MagRes
