/*! \file 
  \brief  Processing commands */

#include "simpsonio.h"
#include "matlabio.h"
#include "cmatrix_utils.h"
#include "NMRsim_Process.h"
#include "NMRsim_logfile.h"
#include "Histogram.h"
#include <sys/stat.h>
#include <sstream>

static const double deg_to_rad=M_PI/180.0;
static const double rad_to_deg=180.0/M_PI;
static const double TWOPI=M_PI*2.0;

BaseList<complex> current_data_row;

void real(ListList<double>& d, const ListList<complex>& s)
{
  d.duplicate_structure(s);
  BaseList<double> drow(d.row());
  real(drow,s.row());
}

// maximum length of Matlab variables
#define NMRSIM_SAVEMAX NMRSIM_MATLAB_VARMAX

static const double sfrq_factor=1e-6; //!< factor for converting pNMRsim sfrq (always in Hz) to external (by default MHz)

double convert_sfrq(double x, int saveflags =0)
{
  if (saveflags & ProcessSave::NOSFRQ)
    return 0.0;
  return sfrq_factor*((saveflags & ProcessSave::POSITIVESFRQ) ? fabs(x) : x);
}

option optsaveprotect("preventoverwrite","prevent file overwrite");

void filestruct::clear()
{
  domain=domain1=D_UNKNOWN;
  sw=sw1=sfrq=sfrq1=ref=ref1=0.0;
}
  
static bool allow_proc_sum=true;

//! flag whether sum allowed in this context
int defsumflags() {
  return allow_proc_sum ? 0 : F_DENYSUM;
}

//! normal flags, allowing sum only if allow_proc_sum set and blocking array / expression variables if 2D active
int defflags() { 
  const static int result=active2D ? (F_DENYEXPR | F_DENYARRAY) : 0;
  return result | defsumflags();
}

ContextWarning<> process_groupsvary_warning("number of processing groups varies between processing stage - this may create very odd results!",&NMRsim_repeat_warning);

bool read_proc_block(procstack_t& pstack, const char* name, const Process_Factory_t& factory, bool allowsum)
{
  allow_proc_sum=allowsum;
  return factory_read(pstack,name,factory,Type2Type<ProcessCommand>(),Type2Type<NullType>());
}

size_t read_proc_blocks(LIST<procstack_t>& pstacks, size_t& accstacks, const char* name, const Process_Factory_t& factory, bool allowsum)
{
  allow_proc_sum=allowsum;
  const size_t nstacks=multiple_factory_read(pstacks,name,factory,Type2Type<ProcessCommand>(),Type2Type<NullType>());
  if (nstacks) {
    if (accstacks<2)
      accstacks=nstacks;
    else {
      if (nstacks!=accstacks)
	process_groupsvary_warning.raise();
    }
  }
  return nstacks;
}

bool try_merge_proc_blocks(LIST<procstack_t>& initialstacks, LIST<procstack_t>& procstacks)
{
  const size_t n=initialstacks.size();
  if ((n==0) || ((n!=procstacks.size()) && !(procstacks.empty())))
    return false;
  const int uses=multiple_overall_attributes(initialstacks);
  const bool mergestacks=((uses & A_SUM)==0);
  if (verbose & VER_PARSE)
    std::cout << "Merging initialproc and proc: " << (mergestacks ? "Yes\n" : "No\n");
  if (!mergestacks)
    return false;
  if (procstacks.empty())
    procstacks.create(n);
  for (size_t i=n;i--;) {
    procstack_t& source(initialstacks(i));
    procstack_t& dest(procstacks(i));
    //    dest.reserve(dest.size()+source.size());
    while (!source.empty()) {
      dest.push_front(source.back());
      source.pop_back();
    }
  }
  initialstacks.clear();
  return true;
}
  
static char scr_fname[2048];
static const int SAVE_FORMATS = ProcessSave::MATLAB | ProcessSave::ASCII | ProcessSave::SIMPSON | ProcessSave::SIMPLOT;
bool ProcessSave::defaultstopoverwrite=false;

Process_Factory_t& get_Process_Factory()
{
  static Process_Factory_t Process_Factory;
  return Process_Factory;
}

Process_Factory_t& get_Finalise_Factory()
{
  static Process_Factory_t Finalise_Factory;
  return Finalise_Factory;
}

typedef std::map<size_t,SaveCommand_function> Save_Factory_t;

Save_Factory_t& get_Save_Factory()
{
  static Save_Factory_t save_Factory;
  return save_Factory;
}

flagsmap_type& get_Save_flags()
{
  static flagsmap_type saveflags;
  return saveflags;
}

void register_save_function(const char* name, SaveCommand_function funcp)
{
  static size_t curflag=ProcessSave::USERFLAGS;
  if (curflag==0)
    throw Failed("register_save_function: flags exhausted");
  get_Save_Factory()[curflag]=funcp;
  get_Save_flags()[name]=curflag;
  curflag<<=1;
}
  
bool havescale=false;

// procflags::procflags()
//   : istimedomain(!isfrequency()),
//     istimedomain1(true),
//     scalevalid(true)
// {}

// void procflags::invalidatescale()
// {
//   scalevalid=false;
// }

processing_state::processing_state(double swv, bool istdv, double sfrqv, double refv)
  : sw(swv),
    istimedomain(istdv),
    sfrq(sfrqv),
    ref(refv)
{}
  
//Warning<> ProcessSave::matlabname_warning("data set name does not start with letter. This may confuse Matlab: ",&NMRsim_once_warning);

namespace {
  smartptr<matlab_controller,false> matlabrawp;
  smartptr<matlab_controller::composite,false> matlabp;
  std::string strbuffer;
  
  void write_simpson_rows(const char* fnamep, const ListList<complex>& FID, const BaseList<processing_state>& pflags, int saveflags)
  {
    simpson_controller ctrl(fnamep,FID.size());
    for (size_t i=FID.size();i--;) {
      simpsonFD spec(FID(i),pflags(i).sw,!(pflags(i).istimedomain));
      spec.sfrq=convert_sfrq(pflags(i).sfrq, saveflags);
      spec.ref=pflags(i).ref;
      ctrl.write(i,spec);
    }
  }

  void scalefirst(cmatrix& FID, double scale, double scale1)
  {
    if (scale1!=1.0) {
      for (size_t n=skips.front();n--;)
	FID.row(n)*=scale1;
    }
    for (size_t n=FID.rows();n--;)
      FID(n,0U)*=scale;
  }

  template<class T> LIST<T> project_row(const Matrix<T>& a)
  {
    const size_t n=a.rows();
    if (!n)
      throw Undefined("project_row");
    LIST<T> d(a.row(0),mxflag::temporary);
    for (size_t j=1;j<n;j++)
      d+=a.row(j);
    return d;
  }
  
  template<class T> LIST<T> project_col(const Matrix<T>& a)
  {
    size_t n=a.rows();
    LIST<T> d(n,T(0),mxflag::temporary);
    for (;n--;) 
      d(n)=sum(a.row(n));
    return d;
  }

  void doshift(BaseList<complex> FID, size_t shift)
  {
    for (size_t c=shift;c--;)
      ::std::swap(FID(c),FID(c+shift));
  }

  void shifthalf(cmatrix& FID, bool is2D)
  {
    size_t shift,r;
    if (FID.cols()>1) {
      if (FID.cols() & 1)
	throw Failed("shifthalf can only be applied to data sets with an even number of points");
      shift=FID.cols()/2;
      for (r=FID.rows();r--;)
	doshift(FID.row(r),shift);
    }
    if (is2D) {
      if ((FID.rows()/skips.front()) & 1)
	throw Failed("shifthalf can only be applied to data sets with an even number of rows");
      shift=FID.rows()/2;
      for (r=shift;r--;) {
	BaseList<complex> FIDr(FID.row(r));
	BaseList<complex> FIDs(FID.row(r+shift));
	for (size_t c=FIDr.size();c--;)
	  ::std::swap(FIDr(c),FIDs(c));
      }
    }
  }
}

void ProcessCommand::print_ftflags(std::ostream& ostr, size_t flags)
{
  if (flags & FT_INV)
    ostr << " -inv";
  if (flags & FT_NOFIRST)
    ostr << " -noscalefirst";
  if (flags & FT_NOSHIFT)
    ostr << " -noshift";
}

bool ProcessSave::openmatlab(const char* rawname) const
{
  const char* basename=getbasename(rawname);
  const size_t rootlength=basename-rawname;
  char tmpname[MAXPATH];
  strncpy(tmpname,rawname,rootlength); //!< copy root
  char* cleanbuf = tmpname+rootlength;
  const char* outp=sanitise_varname(basename,cleanbuf);    
  if (outp!=cleanbuf)
    throw InternalError("ProcessSave::openmatlab");
  //  if (!isalpha(basename[0]))
  //  ProcessSave::matlabname_warning.raise(basename);
  char* endp=strchr(tmpname,'\0');
  strcpy(endp,".mat");
  if (checkoverwrite(tmpname))
    return false;
  *endp='\0'; //suppress temporarily added .mat
  matlabrawp.reset(new matlab_controller(tmpname,5));
  matlabp.reset(new matlab_controller::composite(*matlabrawp,cleanbuf,matlab_controller::STRUCT));
  return true;
}

struct Process_Proxy_ {
  Process_Proxy_() {
    Process_Factory_t& Process_Factory(get_Process_Factory());
    Process_Factory["addlb"]=&ProcessAddLB::create;
    Process_Factory["addnoise"]=&ProcessAddNoise::create;
    Process_Factory["apply"]=&ProcessApply::create;
    Process_Factory["addsignals"]=&ProcessAddSignals::create;
    Process_Factory["conjugate"]=&ProcessConjugate::create;
    Process_Factory["echo"]=&ProcessEcho::create;
    Process_Factory["extract"]=&ProcessExtract::create;
    Process_Factory["fill"]=&ProcessFill::create;
    Process_Factory["ft"]=&ProcessFT::create;
    Process_Factory["ft2d"]=&ProcessFT2d::create;
    Process_Factory["log_file"]=&ProcessLogFile::create;
    Process_Factory["magnitude"]=&ProcessMagnitude::create;
    Process_Factory["newnp"]=&ProcessResample::create;
    Process_Factory["normalise"]=&ProcessNormalise::create;
    Process_Factory["normalize"]=&ProcessNormalise::create;
    Process_Factory["offset"]=&ProcessOffset::create;
    Process_Factory["phase"]=&ProcessPhase::create;
    Process_Factory["resample"]=&ProcessResample::create;
    Process_Factory["rev"]=&ProcessReverse::create;
    Process_Factory["save"]=&ProcessSave::create;
    Process_Factory["scale"]=&ProcessScale::create;    
    Process_Factory["scalefirst"]=&ProcessScaleFirst::create;
    Process_Factory["set"]=&ProcessSet::create;
    Process_Factory["setdomain"]=&ProcessSetDomain::create;
    Process_Factory["shifthalf"]=&ProcessShiftHalf::create;
    Process_Factory["transpose"]=&ProcessTranspose::create;
    //    Process_Factory["sw"]=&ProcessSW::create;
    Process_Factory["zerofill"]=&ProcessZeroFill::create;
    Process_Factory_t& Finalise_Factory(get_Finalise_Factory());
    Finalise_Factory["echo"]=&ProcessEcho::create;
    //    Finalise_Factory["puts"]=&ProcessEcho::create_puts;
    Finalise_Factory["rev"]=&ProcessReverse::create;
    Finalise_Factory["save"]=&ProcessSave::create;
    Finalise_Factory["setdomain"]=&ProcessSetDomain::create;
    Finalise_Factory["system"]=&ProcessSystem::create;
    Finalise_Factory["log_file"]=&ProcessLogFile::create;

    //    optional_map_t& optional_map(get_optional_map());
    add_option(optsaveprotect);
    //    optional_map["preventoverwrite"]=&optsaveprotect;
  }
};

//declare ProcessCommands
static Process_Proxy_ process_proxy_;

void FTObject::negate_alternate(BaseList<complex> FID)
{
  IndirectList<complex,slice> indFID(FID,slice(1,FID.size()/2,2));
  negate_ip(indFID);
}

ThreadWarning<> ProcessFT::doesnothing_warning("FT applied to data set with single point does nothing!",&NMRsim_once_warning);
ThreadWarning<> ProcessFT::nonpowerof2_warning("data set length is not a power of 2 and so not using Fast FT; use zerofill to pad to power of 2?",&NMRsim_once_warning);
ThreadWarning<> ProcessFT::expectingfid_warning("Forward ft applied to spectrum rather than FID may result in reversed data (-nochecks will suppress warning)",&NMRsim_once_warning);
ThreadWarning<> ProcessFT::expectingspectrum_warning("ft -inv applied to FID rather than spectrum may result in reversed data (-nochecks will suppress warning)",&NMRsim_once_warning);

FTObject::FTObject(size_t nv, int flagsv)
  : n(nv), flags(flagsv)
{
  if (n<2)
    ProcessFT::doesnothing_warning.raise();
  doshift=!(flags & ProcessCommand::FT_NOSHIFT);
  doscale=!(flags & ProcessCommand::FT_NOFIRST);
  //  static const int dirmask=(ProcessCommand::FT_INV | ProcessCommand::FT_FID | ProcessCommand::FT_SPECTRUM);
  //switch (flags & dirmask) {
  //case 0: case ProcessCommand::FT_SPECTRUM: case ProcessCommand::FT_INV: case ProcessCommand::FT_FID:
  //  break;
  //default:
  // error_abort("ft: can only specify one of -inv, -fid or -spectrum");
  //}  
  fastft=ispowerof2(n);
  if (!fastft)
    ProcessFT::nonpowerof2_warning.raise();
}

void FTObject::exec(BaseList<complex> FID, const processing_state& pflags) const
{
  if (n<2)
    return;

  int dir=0;
  static const int dirmask=(ProcessCommand::FT_INV | ProcessCommand::FT_FID | ProcessCommand::FT_SPECTRUM);
  switch (flags & dirmask) {
  case 0:
    dir=1;
    if (!nochecks && !(pflags.istimedomain))
      ProcessFT::expectingfid_warning.raise();
    break;

  case ProcessCommand::FT_INV:
    dir=-1;
    if (!nochecks && pflags.istimedomain)
      ProcessFT::expectingspectrum_warning.raise();
    break;
//   case ProcessCommand::FT_SPECTRUM:
//     dir=1;
//     if (pflags.istimedomain)
//       ProcessFT::expectingspectrum_warning.raise();
//     break;
//   case ProcessCommand::FT_FID:
//     dir=-1;
//     if (!(pflags.istimedomain))
//       ProcessFT::expectingfid_warning.raise();
//     break;
  default:
    throw InternalError("FT: unexpected direction flag");
  }  

  if (pflags.istimedomain) {
    if (doshift)
      negate_alternate(FID);
    if (doscale)
      FID.front()*=0.5;
  }
  
  if (fastft)
    fft_ip(FID,dir);
  else
    FID=ft(FID,dir);
  
  if (!(pflags.istimedomain)) { //ensure that TD->FD and FD->TD are inverses
    if (doshift)
      negate_alternate(FID);
    if (doscale)
      FID.front()*=2.0;
  }
}

void ProcessFT::exec(BaseList<complex> FID, processing_state& pflags) const
{
  FTObject obj(FID.size(),ftflags);
  obj.exec(FID,pflags);
  pflags.istimedomain=!pflags.istimedomain;
}

static void split(cmatrix& FIDc, cmatrix& FIDs, const cmatrix& a)
{
  const size_t nrows=a.rows();
  if (!nrows)
    throw Undefined("split");
  if (nrows % 2)
    throw Failed("Hypercomplex signal must have even number of rows");

  const size_t nrows2=nrows/2;
  FIDc.create(nrows2,a.cols());
  FIDs.create(nrows2,a.cols());
  size_t ai=0;
  for (size_t i=0;i<nrows2;i++) {
    FIDc.row(i)=a.row(ai++);
    FIDs.row(i)=a.row(ai++);
  }
}

// void ProcessConjugate::checktd(const processing_state& pflags)
// {
//   static bool donewarn=!allowwarnings();
//   if (!(pflags.istimedomain) && !donewarn) {
//     std::cerr << "Warning (first): conjugate applied to frequency domain data\n";
//     donewarn=true;
//   }
// }

void ProcessMagnitude::exec(BaseList<complex>& FID)
{
  for (size_t i=FID.size();i--;)
    FID(i)=abs(FID(i));
}

void ProcessMagnitude::exec(BaseList<complex> FID, processing_state&) const
{
  exec(FID);
}

void ProcessMagnitude::exec(cmatrix& FID, LIST<processing_state>&) const
{
  // checktd(pflags.back());
  BaseList<complex> asrow(FID.row());
  exec(asrow);
}

void ProcessConjugate::exec(BaseList<complex> FID, processing_state&) const 
{
  //  checktd(pflags);
  exec(FID);
}

void ProcessConjugate::exec(cmatrix& FID, LIST<processing_state>&) const
{
  // checktd(pflags.back());
  FID.conj();
}

void ProcessShiftHalf::exec(BaseList<complex> FID, processing_state& pflag) const
{
  if (FID.size() & 1)
    throw Failed("shifthalf can only be applied to data sets with an even number of points");
  doshift(FID,FID.size()/2);
//   if (pflag.ref) {
//     refreset_warning.raise();
//     pflag.ref=0.0;
//   }
}

void ProcessShiftHalf::exec(cmatrix& FID, LIST<processing_state>& pflags) const
{
  shifthalf(FID,process2D);
//   bool donereset=false;
//   for (size_t i=pflags.size();i--;) {
//     if (pflags(i).ref) {
//       donereset=true;
//       pflags(i).ref=0.0;
//     }
//   }
//   if (donereset)
//     refreset_warning.raise();
}

void process_set_n(int dim, int n, int ni_skip)
{
  raw_set_n(dim,n,ni_skip);
  ensure_array();
}

void ProcessFT2d::exec(cmatrix& FID, LIST<processing_state>& pflags) const
{
  const int dir=(ftflags & FT_INV) ? -1 : 1;

  if (pflags.size()!=2)
    throw Failed("2D FT only valid for 2D data");

  if (!pflags(0U).istimedomain || !pflags(1U).istimedomain)
    throw Failed("2D FT only valid for time-domain data");

  if (!ispowerof2(FID.cols()) || !ispowerof2(FID.rows()))
    error_abort("ft2d only applicable when acquisition dimensions are powers of 2",ERR_FAILED);

  if (!(ftflags & FT_NOFIRST))
    scalefirst(FID,0.5,0.5);

  switch (skips.front()) {
  case 2: {
    cmatrix FIDc,FIDs;
    split(FIDc,FIDs,FID);
    ampfft_ip(FIDc,FIDs,dir);
    if (needphase)
      ampphase_correct_ip(FIDc,FIDs,zero2*deg_to_rad,first2*deg_to_rad,0.0,zero1*deg_to_rad,first1*deg_to_rad,0.0);
    FID.swap(FIDc); //!< 'sin' component is discarded after phasing
    process_set_n(1,FID.rows(),1); //!< change number of rows
  }
    break;
  case 1:
    if (FID.cols()==1) {
      BaseList<complex> FIDr(FID.row());
      fft_ip(FIDr,dir);
    }
    else
      phasefft_ip(FID,dir);
    if (needphase)
      phase_correct_ip(FID,zero2*deg_to_rad,first2*deg_to_rad,0.0,zero1*deg_to_rad,first1*deg_to_rad,0.0);
    break;
  default:
    parser_printthread(std::cerr) << "ft2d only valid with skip of 1 or 2 (" << skips.front() << ")\n";
    error_abort();
  }

  if (!(ftflags & FT_NOSHIFT))
    shifthalf(FID,true);

  pflags(0U).istimedomain=false;
  pflags(1U).istimedomain=false;
}

void ProcessPhase::exec(BaseList<complex> data, processing_state&) const
{
  phase_correct_ip(data,zero*deg_to_rad,first*deg_to_rad,pivot);
}

void
ProcessFT2d::print(std::ostream& ostr) const 
{
  //  assert(subsid==S_NONE);
  ostr << "ft2d ";
  print_ftflags(ostr,ftflags);
  if (needphase)
    ostr << first2 << ' '  << zero2 << ' ' << first1 << ' '  << zero1;
}

void
ProcessPhase::print(std::ostream& ostr) const
{
  ostr << "phase " << first << ' ' << zero << ' ' << "(pivot: " << pivot << ')';
}

void
ProcessPhase::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "phase_";
  switch (subsid) {
  case S_ARG1:
    ostr << "first";
    break;
  case S_ARG2:
    ostr << "zero";
    break;
  default:
    throw InvalidParameter("ProcessPhase::printvariablename");
  }
}

void 
ProcessPhase::set(double v, subsid_t subsid)//, bool)
{
  switch (subsid) {
  case S_ARG1:
    first=v;
    break;
  case S_ARG2:
    zero=v;
    break;
  default:
    throw InternalError("ProcessPhase::set");
  }
}

const char* ProcessSave::makefilename(const char* fnamep, const char* name) const
{
  // const char* bname=name; //!< something odd here - bname not used?
  char* scrp(scr.vector());

  if (name) {
    if (strlen(name)>=NMRSIM_SAVEMAX)
      throw Failed("makefilename: filename too long");
    switch (saveflags_ & SAVE_FORMATS) {
    case SIMPLOT: case SIMPSON:
      strcpy(scrp,insert_simpson(fnamep,name,'_').vector());
      break;
    default:
      if (*fnamep)
		sprintf(scrp,"%s_%s",fnamep,name);
      else
		strcpy(scrp,name);
    }
  }
  else {
    if (fnamep==NMRSIM_NULL)
      throw InternalError("makefilename");
    strncpy(scrp,fnamep,scr.size());
    //    bname=getbasename(fnamep);
  }
  return scrp;
}

bool ProcessSave::checkoverwrite(const char* scrp) const
{
  if ((saveflags_ & STOPOVERWRITE) && isreadable(scrp)) {
    parser_printthread(std::cerr) << "Refusing to overwrite existing file: " << scrp << "\nUse -disable:preventoverwrite to force overwriting\n";
    if (abortonwarning)
      error_abort();
    return true;
  }
  return false;
}

ThreadWarning<> ProcessSave::simpson_notfid_warning("saving freqency-domain data to file with .fid extension",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::simpson_notspec_warning("saving time-domain data to file with .spe extension",&NMRsim_once_warning);

void check_simpson(const char* fname, bool isspec)
{
  if (nochecks)
    return;
  const char* term=fname+strlen(fname)-4;
  if (term<fname)
    return;
  if (strcmp(term,".fid")==0) {
    if (isspec)
      ProcessSave::simpson_notfid_warning.raise();
  }
  else {
    if (!isspec)
      ProcessSave::simpson_notspec_warning.raise();
  }
}

void write_simpson_wrapper(const char* scrp, const BaseList<complex>& data, double sw, bool isspec, double sfrq, double ref)
{
  simpsonFD simpinfo(data,sw,isspec);
  simpinfo.sfrq=sfrq;
  simpinfo.ref=ref;
  check_simpson(scrp,isspec);
  write_simpson(scrp,simpinfo);
}

void write_simpson_wrapper(const char* scrp, const Matrix<complex>& data, double sw, double sw1, bool isspec, double sfrq, double sfrq1, double ref, double ref1)
{
  simpsonFD simpinfo(data,sw,sw1,isspec);
//   if (detect_freq)
//     simpinfo.sfrq=detect_freq;
//   if (sw1)
//     simpinfo.sfrq1=getsfrq1();
  simpinfo.sfrq=sfrq;
  simpinfo.sfrq1=sfrq1;
  simpinfo.ref=ref;
  simpinfo.ref1=ref1;
  check_simpson(scrp,isspec);
  write_simpson(scrp,simpinfo);
}
  
void ProcessSave::save(const BaseList<complex>& data, const processing_state& pflags, const ProcessSave_state& state, const char* name) const
{
  const char* scrp=makefilename(state.fname,name);
  if (!compositesave() && checkoverwrite(scrp))
    return;
  switch (saveflags_ & SAVE_FORMATS) {
  case MATLAB:
    if (matlabp.get())
      matlabp->write(data,name);
    else
      WriteMATLAB(scrp,data,name);
    break;
  case SIMPLOT: case SIMPSON:
    write_simpson_wrapper(scrp,data,pflags.sw,!pflags.istimedomain,convert_sfrq(pflags.sfrq,saveflags_),pflags.ref);
    break;
  default:
    write_vector(scrp,data,mxflag::rmat);
  }
}

void ProcessSave::print(std::ostream& ostr) const
{
  ostr << "save " << fname.vector();
}

namespace {
  inline void write_ascii(const char* name, const rmatrix& a) { write_matrix(name,a); }
  inline void write_ascii(const char* name, const BaseList<double>& a) { write_vector(name,a); }
  inline void write_ascii(const char*, const BaseList<bool>&) { throw InternalError("write_ascii"); }
}

void ProcessSave::makescale(LIST<double>& scale, const processing_state& pflags, size_t np, size_t skip) const
  //, bool istimedomain)
{
  const double locsw=pflags.sw;
  scale.create(np);
  if (pflags.istimedomain) {
    // ignore rev flag - will have generated warning elsewhere
    const double dt=1/(skip*locsw); //fudge dt for skip factor
    for (size_t i=0;i<np;i+=skip) {
      for (size_t j=0;j<skip;j++) 
	scale(i+j)=i*dt*1e6;
    }
  }
  else {
    // const int offset=np/2;
    // const double df=locsw/np;
    const bool rev= (saveflags_ & REVERSEFREQ);
    const double invN=1.0/double(np);
    for (size_t i=0;i<np;i++) {//no skip if transformed
      const size_t desti = rev ? np-1-i : i;
      scale(desti)=locsw*(invN*i-0.5)+pflags.ref;
      //	scale(i)=(i-offset)*df;
    }
  }   
}

ThreadWarning<> ProcessSave::simpsonunsuitable_warning("SIMPLOT/SIMPSON format not suitable for data set (ignored): ",&NMRsim_once_warning);

template<class T> void ProcessSave::save(const T& data, const ProcessSave_state& state, const char* name) const
{
  const char* scrp=makefilename(state.fname,name);
  if (!compositesave() && checkoverwrite(scrp))
    return;
  switch (saveflags_ & SAVE_FORMATS) {
  case MATLAB:
    if (matlabp.get())
      matlabp->write(data,name);
    else
      WriteMATLAB(scrp,data,name);
    break;
  case SIMPLOT: case SIMPSON:
    simpsonunsuitable_warning.raise(name);
    break;
  default:
    write_ascii(scrp,data);
  }
}

namespace {
  size_t save_cols=-1;
  bool save_isspec,save_isrealnD;
}

FILE* ProcessSave::opencomments(const char* outname, const char* mode) const
{
  if (saveflags_ & MATLAB) {
    strbuffer="";
    return NMRSIM_NULL;
  }
  FILE* fp=NMRSIM_NULL;
  if (!checkoverwrite(outname))
    fp=fopen(outname,mode);
  if (fp==NMRSIM_NULL)
    throw Failed("opencomments: failed to open output file");
  return fp;
}

void ProcessSave::writeline(FILE* foutp, const char* str, char comchar) const
{
  if (saveflags_ & MATLAB) {
    strbuffer+=str;
    strbuffer+='\n';
  }
  else {
    if ((saveflags_ & ASCII) && comchar)
      fputc('#',foutp);
    fprintf(foutp,"%s\n",str);
  }
}

void ProcessSave::closecomments(FILE* fp, const char* name) const
{
  if (saveflags_ & MATLAB) {
    matlabp->write(strbuffer.c_str(),name);
    return;
  }
  if (fp)
    fclose(fp);
}

bool ProcessSave::original_isflipped() const
{
  if (original_revguardp)
    return original_revguardp->isflipped();
  return false;
}

void ProcessSave::savesource(const char* outname) const
{  
  const char* origin(parser_getfname());
  FILE* finp=fopen(origin,"r");
  FILE* foutp=NMRSIM_NULL;
  if (!finp)
    throw Failed("ProcessSave::savesource: can't open original input file!");
  try {
    foutp=opencomments(outname);
    std::ostringstream envline(std::ostringstream::out);
    dumpenvironment(envline);
    std::string comline(envline.str());
    comline.insert(0,"Command line: ");
    if (global_argc==0)
      throw InternalError("savesource");
    for (size_t i=0;i<global_argc;i++) {
      comline+=global_argv[i];
      comline+=' ';
    }
    writeline(foutp,comline.c_str());

    comline="Input file: ";
    comline+=origin;
    writeline(foutp,comline.c_str());

    char linebuf[MAXLINE];
    while (fgets(linebuf,MAXLINE,finp)) {
      char* p=linebuf+strlen(linebuf)-1; //!< strip trailing whitespace
      while (isspace(*p) && (p>=linebuf))
	--p;
      p[1]='\0';
      writeline(foutp,linebuf);
    }
    
    closecomments(foutp,"source");

  } catch (...) {
    fclose(finp);
    closecomments(foutp,"source");
    throw;
  }
  fclose(finp);
}
 
ThreadWarning<> ProcessSave::notnd_warning("File format does not support >2D data",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::notirregular_warning("ASCII format does not support irregular data",&NMRsim_once_warning);	    
ThreadWarning<> ProcessSave::nosource_warning("can't save source to this format (ignored)",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::realnotvalid_warning("-realonly flag is ignored for this format",&NMRsim_once_warning);

template<class TC, class TR> void raw_matlab_write(matlab_controller::composite& matlabc, const TC& data, const char* bname, const Type2Type<TR>, bool savereal)
{
  if (savereal) {
    TR rdata;
    real(rdata,data);
    matlabc.write(rdata,bname);
  }
  else
    matlabc.write(data,bname);
}
 
void ProcessSave::save(const DataStore& a, const ProcessSave_state& state, const char* bname) const
{
  if (a.empty()) {
    nothingtosave_warning.raise(bname);
    return;
  }

  //const char* scrp=NMRSIM_NULL;
  // if (bname)
  const char* scrp=makefilename(state.fname,bname);
//   else {
//     scrp=state.fname;
//     bname=getbasename(state.fname);
//   }
  if (!compositesave() && checkoverwrite(scrp))
    return;
  const processing_state statedir(a.states().back());
  const double lsw=statedir.sw;
  const double lsfrq=convert_sfrq(statedir.sfrq, saveflags_);
  const double lref=statedir.ref;
  double save_sw1=0;
  double save_sfrq1=0;
  double save_ref1=0;
  if (a.isnD() && (ndims>1)) {
    const processing_state& stateind(a.state(ndims-2));
    save_sw1=stateind.sw;
    save_sfrq1=convert_sfrq(stateind.sfrq, saveflags_);
    save_ref1=stateind.ref;
  }
  const bool savereal=(saveflags_ & REALONLY);
 
  try {
    switch (saveflags_ & SAVE_FORMATS) {
    case MATLAB: {
      if (a.rows()==1)
	raw_matlab_write(*matlabp,a.row(),bname,Type2Type< List<double> >(),savereal);
      else {
	if (a.isnD()) {
	  switch (ndims) {
	  case 1: case 2:
	    raw_matlab_write(*matlabp,a.matrix(),bname,Type2Type< Matrix<double> >(),savereal);
	    break;
	  case 3:
	    if (save_cols==1)
	      raw_matlab_write(*matlabp,a.matrix(),bname,Type2Type< Matrix<double> >(),savereal);
	    else { //copying out is wasteful, but not worth optimising
	      const MultiMatrix<complex,3> mm(save_cols,array_ns.front(),array_ns(1U),a.row().vector());
	      raw_matlab_write(*matlabp,mm,bname,Type2Type< MultiMatrix<double,3> >(),savereal);
	    }
	    break;
	  case 4:
	    if (save_cols==1) {
	      const MultiMatrix<complex,3> mm(array_ns.front(),array_ns(1U),array_ns(2U),a.row().vector());
	      raw_matlab_write(*matlabp,mm,bname,Type2Type< MultiMatrix<double,3> >(),savereal);
	    }
	    else {
	      const MultiMatrix<complex,4> mm(save_cols,array_ns.front(),array_ns(1U),array_ns(2U),a.row().vector());
	      raw_matlab_write(*matlabp,mm,bname,Type2Type< MultiMatrix<double,4> >(),savereal);
	    }
	    break;
	    // 	  case 4:  {
	    // 	    const MultiMatrix<complex,5> mm(cols,array_ns.front(),array_ns(1U),array_ns(2U),array_ns(3U),a.row());
	    // 	    matlabp->write(mm,bname);
	    // 	  }
	    // 	    break;
	  default:
	    throw Failed("save: unhandled Matlab dimensionality");
	  }	  
	}
	else
	  raw_matlab_write(*matlabp,a.listlist(),bname,Type2Type< ListList<double> >(),savereal);
	//	  matlabp->write(a.listlist(),bname);
      }
    }
      break;
    case SIMPSON:
      if (savereal)
	realnotvalid_warning.raise();

      if (a.isnD()) {
	if (ndims>2)
	  notnd_warning.raise();
	write_simpson_wrapper(scrp,a.matrix(),lsw,save_sw1,save_isspec,lsfrq,save_sfrq1,lref,save_ref1);
      }
      else {
	if (a.rows()==1)
	  write_simpson_wrapper(scrp,a.row(),lsw,save_isspec,lsfrq,lref);
	else
	  write_simpson_rows(scrp,a.listlist(),a.states(),saveflags_);
      }
      break;
    case SIMPLOT:
      if (savereal)
	realnotvalid_warning.raise();

      if (save_isrealnD) {
	write_simpson_rows(scrp,a.matrix(),lsw,save_isspec,lsfrq);
	scrp=NMRSIM_NULL;
      }
      else {
	if (a.isnD() || (a.rows()==1))
	  write_simpson_wrapper(scrp,a.row(),lsw,save_isspec,lsfrq,lref);
	else {
	  write_simpson_rows(scrp,a.listlist(),a.states(),saveflags_);
	  scrp=NMRSIM_NULL; //!< flag can't write
	}	  
      }
      break;
    case ASCII:
      try {
	if (save_isrealnD) {
	  if (ndims>2)
	    notnd_warning.raise();
	  if (savereal) {
	    rmatrix rdata;
	    real(rdata,a.matrix());
	    write_matrix(scrp,rdata);
	  }
	  else
	    write_matrix(scrp,a.matrix());
	}
	else {
	  if (!a.isnD() && (a.rows()>1))
	    notirregular_warning.raise();
	  if (savereal) {
	    List<double> rdata;
	    real(rdata,a.row());
	    write_vector(scrp,rdata);
	  }
	  else
	    write_vector(scrp,a.row());
	}
      } catch (const std::exception& exc) {
	parser_printthread(std::cerr) << "Write failed: " << exc.what() << '\n';
      }
      break;
    default:
      error_abort("ProcessSave: Unknown write format");
    }
    if (saveflags_ & SOURCE) {
      if (scrp)
	savesource(scrp);
      else
	nosource_warning.raise();
    }
  } catch (const std::exception& exc) {
    parser_printthread(std::cerr) << exc.what() << '\n';
  }
}

ThreadWarning<> ProcessSave::noscale_warning("scale save not possible",&NMRsim_repeat_warning);
ThreadWarning<> ProcessSave::sumprojection_warning("-sum / -projection save requires nD data set",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::sumprojection_warning2("-sum / -projection save on single row data set",&NMRsim_once_warning);

void ProcessSave::exec(BaseList<complex> FID, processing_state& pflags) const
{
  substitute_string(scr_fname,sizeof(scr_fname),fname.vector(),SUB_ESCAPE | SUB_NONCONSTWARN);
  ProcessSave_state state(scr_fname);

  //   smartptr<matlab_controller,false> matlabrawp;
  if (matlabp.get()!=NMRSIM_NULL)
    throw InternalError("ProcessSave::exec");
   const bool ismatlab=(saveflags_ & MATLAB);
   if (ismatlab && !openmatlab(state.fname))
     return;

   if (!(saveflags_ & NODATA)) {
     save(FID,pflags,state,ismatlab ? "data" : (const char*)NMRSIM_NULL); //!< use "data" qualifier for Matlab "struct" write
     if (saveflags_ & MATLAB) {
       if (pflags.istimedomain)
	 save(ExplicitList<1,double>(pflags.sw),state,"sw");
       else {
	 const double dt= pflags.sw ? 1.0/pflags.sw : 0.0;
	 save(ExplicitList<1,double>(dt),state,"dt");
       }
     }
   }
   writevars(state,true);

   if (saveflags_ & SCALE) {
     if (pflags.sw) {
       LIST<double> scale;
       //       makescale(scale,pflags.sw,FID.size(),1,pflags.istimedomain);
       makescale(scale,pflags,FID.size(),1);
       save(scale,state,"scale");
     }
     else
       noscale_warning.raise();
   }

   if (saveflags_ & (SUM | PROJ))
     sumprojection_warning.raise();

   writeparametersoptions(state,row_index);
}

ThreadWarning<> ProcessSave::invalidmatlab_warning("invalid characters (replaced by _) or too many characters (truncated) in MATLAB variable name: ",&NMRsim_once_warning);

const char* ProcessSave::sanitise_varname(const char* in, char cleanbuf[NMRSIM_MATLAB_VARMAX+1]) const
{
  if (!(saveflags_ & MATLAB))
    return in;

  static char squash_char='\0';

  if (!squash_char) {
    squash_char=NMRSIM_MATLAB_SQUASH_CHAR;
    if (!isalpha(squash_char))
      throw InternalError("Illegal NMRSIM_MATLAB_SQUASH_CHAR (must be alphabet char)");
  }

  const bool oversize=(strlen(in)>NMRSIM_MATLAB_VARMAX);
  bool problem=oversize;
  if (oversize) {
    strncpy(cleanbuf,in,NMRSIM_MATLAB_VARMAX-1);
    cleanbuf[NMRSIM_MATLAB_VARMAX]='\0'; 
  }
  else
    strcpy(cleanbuf,in);

  for (char* source=cleanbuf;*source;source++) {
    bool ok;
    if (source==cleanbuf)
      ok=isalpha(*source);
    else
      ok=isalnum(*source) || (*source=='_');
    if (!ok) {
      *source=squash_char;
      problem=true;
    }
  }
  if (problem)
    invalidmatlab_warning.raise(in);
  return cleanbuf;
}

ThreadWarning<> ProcessSave::nonconstvariable_warning("Trying to save non-constant variable - explore using log_file mechanism instead? ",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::emptyvariable_warning("Ignoring attempt to save empty variable: ",&NMRsim_once_warning);

void ProcessSave::writevars(const ProcessSave_state& state, bool currentrow) const
{
  char cleanbuf[NMRSIM_MATLAB_VARMAX+1];

  for (size_t nv=0;nv<vars_.size();nv++) {
    const RealVariable& cvar(*(vars_(nv)));
    const char* cleanname=sanitise_varname(cvar.name(),cleanbuf);    
    if (!(cvar.isconstant()) && !currentrow && cvar.hasbase()) {
      const VariableBase& var(cvar.value());
      if ( var.isarray() && (var.get_array().size()==1)) {
	const BaseList<double> vals(var.get_array().row());
	const size_t nvals=vals.size();
	const size_t cols=var.array_item_length();
	if ((cols==0) || (nvals % cols))
	  throw InternalError("ProcessSave::writevars");
	const Matrix<double> mvals(nvals / cols, cols, vals.vector(), mxflag::nondynamic);
	save(mvals,state,cleanname);
	continue;
      }
      nonconstvariable_warning.raise(cvar.name());
    }
    const BaseList<double> vals(cvar.get_list());
    if (vals.empty())
      emptyvariable_warning.raise(cvar.name());
    else
      save(vals,state,cleanname);
  }
}

ThreadWarning<> ProcessSave::noparameters_warning("-parameters save requested but no active simulation parameters",&NMRsim_repeat_warning);

//Warning<> ProcessSave::sumparameters_warning("-parameters save outside summation array will only save last set of parameters (place save in initialproc?)",&NMRsim_once_warning);
ThreadWarning<> ProcessSave::sumparameters_warning("-parameters save outside summation array will only save last set of parameters",&NMRsim_once_warning);

void ProcessSave::writeparametersoptions(const ProcessSave_state& state, int row) const
{
  if (saveflags_ & VARIABLES) {

    if (!vararray && !valerrs)
      noparameters_warning.raise();
    else {
      
      if ( ((saveflags_ & SAVE_FORMATS)==SIMPSON) || ((saveflags_ & SAVE_FORMATS)==SIMPLOT))
	simpsonunsuitable_warning.raise("parameters");
      else {
	const bool usevalerrs = (row<0) && !(valerrs.empty());
	if (!within_evaluation() && (sum_n0>1) && !usevalerrs)
	  sumparameters_warning.raise();
	if (varnames.size()==0)
	  throw InternalError("writeparametersoptions");
	
	if ((saveflags_ & MATLAB) || !usevalerrs || (valerrs.rows()!=varnames.size())) {
	  if (row>=0)
	    save(vararray.row(size_t(row)),state,"parameters");
	  else
	    save(usevalerrs ? valerrs : vararray,state,"parameters");
	  
	  
	  const char* scrp=makefilename(state.fname,"parameternames");
	  FILE* foutp=opencomments(scrp,"w");
	  for (size_t i=0;i<varnames.size();i++)
	    writeline(foutp,varnames(i),'\0');
	  closecomments(foutp,"parameternames");
	}
	else {
	  char buffer[MAXLINE];
	  
	  const char* scrp=makefilename(state.fname,"parameters");
	  FILE* foutp=opencomments(scrp,"w");
	  for (size_t i=0;i<varnames.size();i++) {
	    snprintf(buffer,sizeof(buffer),"%s \t%g \t%g",varnames(i),valerrs(i,size_t(0)),valerrs(i,size_t(1)));	  
	    writeline(foutp,buffer,'\0');
	  }
	  closecomments(foutp,"parameters");	
	}
      }
    }
  }

  for (size_t i=0;i<saveoptions.size();i++) {
    SaveCommand_function funcp(saveoptions(i));
    if (!funcp)
      throw InternalError("ProcessSave");
    (*funcp)(*this,state,row);
  }

  matlabp.clear(); //!< This will flush any composite object (the matlab_controller itself is closed when the function exits)
}

ThreadWarning<> ProcessSave::noscaleforirregular_warning("scale save not possible for irregular data",&NMRsim_once_warning);

void DataStore::reset(bool istd)
{
  if (empty())
    return;
  (*this)=complex(0.0);
  // const bool istd=!isfrequency();
  for (size_t i=procstates_.size();i--;)
    procstates_(i).istimedomain= (!active2D || (i==0)) ? istd : true;
}
      
void DataStore::copy_structure(LIST<size_t>& sizes) const
{
  sizes.create(rows());
  for (size_t i=sizes.size();i--;)
    sizes(i)=row(i).size();
}

//! try to compress listlist into matrix
cmatrix DataStore::trymatrix() const
{
  if (rows()==0)
    throw Undefined("trymatrix");
  if (isnD())
    return cmatrix(matrix(),mxflag::nondynamic);
  size_t commons=0U;
  const ListList<complex>& source(listlist());
  for (size_t n=0;n<source.size();n++) {
    const size_t cursize=source(n).size();
    if (n==0)
      commons=cursize;
    else {
      if (commons!=cursize)
	throw Failed("trymatrix: irregular data set");
    }
  }
  return cmatrix(source.size(),commons,source.row(),mxflag::nondynamic);
}

ReverseGuard::ReverseGuard(DataStore& FIDv, const char* labelv, bool dorev)
  : FID_(FIDv), flipped_(false), label_(labelv) 
  {
    if (dorev)
      flip();
  }

ReverseGuard::~ReverseGuard()
{
  if (flipped_)
    flip();
}

void ReverseGuard::flip()
{
  ProcessReverse revobj;
  revobj.rawexec(FID_);
  flipped_ = !flipped_;
  if ((label_!=NMRSIM_NULL) && (verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Reversed data set: " << label_ << '\n';
}


void ProcessSave::rawexec(DataStore& FID) const
{
  const bool needflip=saveflags_ & REVERSEFREQ;
  ReverseGuard revguard_data(FID, "data", needflip); //!< reverse guard OK here as flip isn't undone until 
  if (needflip && !(fit_set.empty())) {
    ReverseGuard revguard_original(fit_set, "original data set");
    original_revguardp=&revguard_original;
    rawexec_postrev(FID);
  }
  else {
    original_revguardp=NMRSIM_NULL;
    rawexec_postrev(FID);
  }
}

void ProcessSave::rawexec_postrev(DataStore& FID) const
{
  substitute_string(scr_fname,sizeof(scr_fname),fname.vector(),SUB_ESCAPE | SUB_NONCONSTWARN);
  ProcessSave_state state(scr_fname);
  //const char* fnamep=scr_fname;  
  const bool isnD=(FID.isnD() && FID.rows()>1);
  save_cols=isnD ? FID.matrix().cols() : 0;
  save_isrealnD=isnD && (save_cols>1);

  //  smartptr<matlab_controller,false> matlabrawp;
  if (matlabp.get()!=NMRSIM_NULL)
    throw InternalError("ProcessSave::rawexec");
  const bool ismatlab=(saveflags_ & MATLAB);
  if (ismatlab && !openmatlab(state.fname))
    return;
 	       
  const BaseList<processing_state> pflags(FID.states());
  save_isspec=!(pflags.back().istimedomain);

  char scrf[8];  
  if (!(saveflags_ & NODATA)) {
    if ((verbose & VER_GEN) && (verbose_level>1)) {
      const BaseList<complex> rawdata(FID.row());
      std::cout << "Saving data with first point " << rawdata.front() << " and last point " << rawdata(rawdata.size()-1) << '\n';
    }
    save(FID,state,ismatlab ? "data" : (const char*)NMRSIM_NULL); //!< use "data" qualifier for Matlab "struct" write
    if (saveflags_ & MATLAB) {
      const size_t nsws=FID.states().size();
      bool isalltd=true;
      //      bool isallfd=true;
      for (size_t i=nsws;i--;) {
	if (!(pflags(i).istimedomain))
	  isalltd=false;
	//	else
	//  isallfd=false;
      }      
      const bool reverse=FID.isnD(); //!< reverse order of sw etc. for nD data
      ScratchList<double> loc(nsws,0.0);
      ScratchList<double> sfrq(nsws,0.0);
      ScratchList<double> ref(nsws,0.0);
      bool havesfrq=false;
      bool haveref=false;
//       if (detect_freq) {
// 	sfrq.front()=detect_freq;
// 	havesfrq=true;
//       }
//       const double sfrq1=getsfrq1();
//       if ((nsws>1) && sfrq1) {
// 	sfrq(size_t(1))=sfrq1;
// 	havesfrq=true;
//       }
      ScratchList<bool> type(nsws);
      for (size_t i=nsws;i--;) {
	const size_t dest=reverse ? nsws-i-1 : i;
	const double lsw=pflags(i).sw;
	const double lsfrq=convert_sfrq(pflags(i).sfrq);
	const double lref=pflags(i).ref;
	if (lsfrq) {
	  havesfrq=true;
	  sfrq(dest)=lsfrq;
	}
	if (lsw)
	  loc(dest)= isalltd ? 1.0/lsw : lsw;
	type(dest)= pflags(i).istimedomain ? FD_TYPE_FID : FD_TYPE_SPE;
	if (lref) {
	  haveref=true;
	  ref(dest)=lref;
	}
      }
      save(loc,state, isalltd ? "dt" : "sw");
      save(type,state,"type");
      if (havesfrq)
	save(sfrq,state,"sfrq");
      if (haveref)
	save(ref,state,"ref");
    }
  }

  writevars(state,false);

  if (saveflags_ & SCALE) {
    if (FID.isnD()) {
      const cmatrix& FIDm(FID.matrix());
      LIST<double> scale;
      //      bool failed=false;
      const bool simple(ndims==1);
      for (size_t dim=0;dim<ndims;dim++) {
	const double sw=pflags(dim).sw;
	if (sw==0.0)
	  continue;
	//	const bool istd=pflags(dim).istimedomain;
	if (dim==ndims-1) {
	  makescale(scale,pflags(dim),FIDm.cols(),1);
	  //	  makescale(scale,sw,FIDm.cols(),1,istd);
	  save(scale,state,"scale");
	}
	else {
	  sprintf(scrf,"scale%lu",(unsigned long)(dim+1U));
	  //if simple 2D, use actual number of rows (zero-filled?), otherwise use original specification - problematic if this has changed
	  //	  makescale(scale,sw,simple ? FIDm.rows() : array_ns(dim),array_skips(dim),istd);
	  makescale(scale,pflags(dim),simple ? FIDm.rows() : array_ns(dim),array_skips(dim));
	  save(scale,state,scrf);
	}
	// 	else
	// 	  failed=true;
      }
      //else
      //failed=true;
      //if (failed)
      //	std::cerr << "Warning: not all spectral width(s) valid: scale(s) not saved\n"; 
    }
    else
      noscaleforirregular_warning.raise();
  }

  if (saveflags_ & (SUM | PROJ)) {
    try {
      const cmatrix FIDm(FID.trymatrix());
      if (saveflags_ & SUM)
	save(project_row(FIDm),pflags.back(),state,"sum");
    
      if (saveflags_ & PROJ) {
	if (pflags.back().istimedomain) {
	  const LIST<complex> FID0(FIDm(0U,range()));
	  save(FID0,pflags.front(),state,"projection");
	}
	else
	  save(project_col(FIDm),pflags.front(),state,"projection");
      }
      if (FID.rows()==1)
	sumprojection_warning2.raise();	
    } 
    catch (const Failed&) {
      sumprojection_warning.raise();
    }
  }

  writeparametersoptions(state);
//   if (saveflags_ & STATISTICS) {
//     if (fitsaveok()) {
//       dosave(fnamep,"covariance",covar);
//       if (!residuals.empty())
// 	rawexec_(residuals,fnamep,"residuals");
//     }
//   }  

//   if ((saveflags_ & ORIGINAL) && fitsaveok())
//     rawexec_(fit_set,fnamep,"original");

  matlabp.clear(); //!< This will flush any composite object (the matlab_controller itself is closed when the function exits)
}

ProcessEcho::ProcessEcho(const char* str_)
  : ProcessCommand(PROC_HASBOTH),
    BaseEcho(str_)
{
  const int uses=update(); //!< get uses;
  if (uses & A_ARRAY)
    flagvariable();
}

ProcessCommand* ProcessEcho::create()
{
  const char* echoline=get_curline();
  // const bool hasdollar=(strchr(echoline,'$')!=NMRSIM_NULL); //!< if contains $ then indicate that echo string might vary - NASTY HACK
  ProcessCommand* newac = new ProcessEcho(echoline);
  set_curline(NMRSIM_NULL);
  return newac;
}

// ProcessCommand* ProcessEcho::create_puts()
// {
//   return new ProcessEcho(get_token());
// }

ProcessCommand* ProcessLogFile::create()
{
  char* name;
  int logflags;
  logfile_controller::parse_raw(name,logflags);
  return new ProcessLogFile(name,logflags);
}

// NB assumes log file only written per data set rather than row-oriented
void ProcessLogFile::rawexec(DataStore&) const
{
  update_log_file(name.c_str(),logflags);
}

ContextWarning<> ProcessSave::nothingtosave_warning("Data set is empty - nothing to save: ",&NMRsim_repeat_warning);
ContextWarning<> ProcessSave::sfrq_warning("-nosfrq, -positivesfrq flags misused: only valid with SIMPSON/SIMPLOT formats and cannot be sensibly combined",&NMRsim_once_warning);
ContextWarning<> ProcessSave::nosimpsonrevfreq_warning("-reversefrequency incompatible with SIMPSON/SIMPLOT formats since data order is fixed. Use rev if if you really need to flip the data",&NMRsim_once_warning);
ContextWarning<> ProcessSave::revfreqscale_warning("note that -reversefrequency will also reverse order of output frequency scale. Use rev to flip data set without changing normal low to high frequency order of saved scale",&NMRsim_once_warning);

ProcessCommand* ProcessSave::create()
{
  flagsmap_type& saveflags(get_Save_flags());
  static bool doneadd=false;
  if (!doneadd) {
    saveflags["simpson"]=SIMPSON;
    saveflags["simplot"]=SIMPLOT;
    saveflags["matlab"]=MATLAB;
    saveflags["ascii"]=ASCII;
    saveflags["projection"]=PROJ;
    saveflags["nodata"]=NODATA;
    saveflags["sum"]=SUM;
    saveflags["scale"]=SCALE;
    saveflags["parameters"]=VARIABLES;
    saveflags["realonly"]=REALONLY;
    //saveflags["statistics"]=STATISTICS;
    saveflags["source"]=SOURCE;
    //saveflags["original"]=ORIGINAL;
    saveflags["nosfrq"]=NOSFRQ;
    saveflags["positivesfrq"]=POSITIVESFRQ;
    saveflags["reversefrequency"]=REVERSEFREQ;
  }
  const char* name(parse_string(F_IGNOREDOLLAR));

  LIST<RealVariable*> vars;
  while (parser_isnormal()) {
    const char* varname=clean_variable_name(parse_string());
    vars.push_back(&(parse_variable_name(varname)));
  }

  size_t flags=parse_flags(saveflags);
  if (((flags & ~SAVE_FORMATS)==NODATA) && vars.empty()) {
    nothingtosave_warning.raise("data");
    return NMRSIM_NULL;
  }
  if (defaultstopoverwrite)
    flags|=STOPOVERWRITE;

  if (!nochecks && (flags & REVERSEFREQ) && (flags & SCALE))
    revfreqscale_warning.raise();

  return new ProcessSave(name,vars,flags);
}

ProcessSave::ProcessSave(const char* fname_,const BaseList<RealVariable*>& varsv, int flagsv)
  : ProcessCommand(PROC_HAS2D),
    fname(strlen(fname_)+1,const_cast<char*>(fname_)),
    vars_(varsv),
    saveflags_(flagsv),
    scr(strlen(fname_)+NMRSIM_SAVEMAX)
{
  switch (saveflags_ & SAVE_FORMATS) {
  case 0:
    saveflags_|=SIMPLOT;
  case SIMPSON: case SIMPLOT: case ASCII: case MATLAB: break;
//   case MATLAB: {
//     char cleanbuf[NMRSIM_MATLAB_VARMAX+1];
//     const char* cleanname=sanitise_varname(fname_,cleanbuf);    
//     fname=BaseList<char>(strlen(cleanname)+1,cleanbuf);
//     //    if (strchr(fname_,'.')) //catch accidental use of .spe etc.
//     //  error_abort("Matlab filename cannot contain .");
//   }
//    break;
  default:
    error_abort("Can't specify more than one format (Matlab, SIMPSON, ASCII)");
  }
  size_t curmask = USERFLAGS;
  size_t extraflags=saveflags_ & ~(curmask-1);
  Save_Factory_t& save_factory(get_Save_Factory());
  while (extraflags) {
    if (extraflags & curmask) {
      SaveCommand_function funcp(save_factory[curmask]);
      saveoptions.push_back(funcp);
      extraflags-=curmask;
    }
    curmask<<=1;
  }
  const bool isSIMPSON=(saveflags_ & (SIMPLOT | SIMPSON));
  if (isSIMPSON && (saveflags_ & REVERSEFREQ)) {
    nosimpsonrevfreq_warning.raise();
    saveflags_ &= ~REVERSEFREQ;
  }
  if (!isSIMPSON && (saveflags_ & (NOSFRQ | POSITIVESFRQ))) {
    sfrq_warning.raise();
    saveflags_ &= ~(NOSFRQ | POSITIVESFRQ);
  }
  else {
    if ((saveflags_ & NOSFRQ) && (saveflags_ & POSITIVESFRQ))
      sfrq_warning.raise();
  }
}

void ProcessScale::print(std::ostream& ostr) const 
{
  ostr << "scale " << scale;
}

void ProcessScale::printvariablename(std::ostream& ostr, subsid_t) const 
{
  ostr << "scale";
}

void ProcessScale::exec(cmatrix& FID, LIST<processing_state>&) const
{ 
  if (verbose & VER_GEN)
    std::cout << "Scaling 2D FID by " << scale << '\n';
  FID*=scale;
}

void ProcessScale::exec(BaseList<complex> FID, processing_state&) const
{ 
  if (verbose & VER_GEN)
    std::cout << "Scaling 1D FID by " << scale << '\n';
  FID*=scale;
}

void ProcessScaleFirst::set(double v, subsid_t subsid)
{
  switch (subsid) {
  case S_ARG1:
    scale=v;
    break;
  case S_ARG2:
    scale1=v;
    break;
  default:
    throw InternalError("ProcessScaleFirst::set");
  }
}

void ProcessScaleFirst::print(std::ostream& ostr) const
{
  ostr << "scalefirst " << scale << ' ' << scale1;
}

void ProcessScaleFirst::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "scalefirst_";
  switch (subsid) {
  case S_ARG1:
    ostr << "scale";
    break;
  case S_ARG2:
    ostr << "scale1";
    break;
  default:
    throw InvalidParameter("ProcessScaleFirst::printvariablename");
  }
}

ThreadWarning<> ProcessScaleFirst::fd_warning("scalefirst applied to frequency domain data",&NMRsim_once_warning);

void ProcessScaleFirst::exec(BaseList<complex> FID, processing_state& flags) const
{
  if (flags.istimedomain)
    fd_warning.raise();
  FID.front()*=scale;
}

void ProcessScaleFirst::exec(cmatrix& FID, LIST<processing_state>&) const
{ 
  scalefirst(FID,scale,scale1);
}

void rawreverse(BaseList<complex> spec)
{
  const size_t offset=spec.size()-1;
  for (size_t i=spec.size()/2;i--;)
    std::swap(spec(i),spec(offset-i));
}

void ProcessReverse::exec(BaseList<complex> spec)
{
  rawreverse(spec);
}

ThreadWarning<> ProcessReverse::td_warning("rev applied to time domain data",&NMRsim_once_warning);

void ProcessReverse::exec(BaseList<complex> spec, processing_state& pflags) const
{
  if (pflags.istimedomain)
    td_warning.raise();
  //  else {
  //  if (pflags.ref) {
  //    refreset_warning.raise();
  //    pflags.ref=0.0; //!< reset reference to centre
  //  }
  //}
  exec(spec);
}

ProcessCommand* ProcessOffset::create()
{
  Variable cvar(S_ARG1);
//   static flagsmap_type flags;
//   if (flags.empty())
//     flags["imag"]=1;
    //  return new ProcessOffset(parse_double(&cvar,process2D ? F_DENYARRAY : 0));
  VariableBase* varp=parse_double_variable(cvar,defsumflags() | F_ALLOWLIST);
  bool isconst=varp->isconst();
  return new ProcessOffset(*varp,parsecomplexity(ProcessCommand::REAL),isconst);
}

void ProcessOffset::exec(BaseList<complex> FID, processing_state&) const
{
  const BaseList<double> addlist(offsetvar.value());
  const size_t effn=checkcomplexity(addlist.size(),addtype);
  if (effn==1) {
    const double offset(addlist.front());
    switch (addtype) {
    case IMAG:
      FID+=complex(0,offset);
      break;
    case REAL:
      FID+=offset;
      break;
    default:
      FID+=complex(addlist.front(),addlist(1U));
    }
  }
  else {
    if (effn!=FID.size()) {
      parser_printthread(std::cerr) << "offset: mismatch between length of data set (" << FID.size() << ") and expression result (" << effn << ")\n";
      error_abort();
    }
    switch (addtype) {
    case IMAG:
      imags(FID)+=addlist;
      break;
    case REAL:
      FID+=addlist;
      break;
    default:
      FID+=ascomplex(addlist);
    }
  }
}

void ProcessOffset::print(std::ostream& ostr) const
{
  ostr << "offset ";
  offsetvar.print(ostr);
  switch (addtype) {
  case IMAG:
    ostr << " -imag"; break;
  case COMPLEX:
    ostr << " -complex"; break;
  default:
    break;
  }
}


void ProcessOffset::printvariablename(std::ostream& ostr, subsid_t) const
{
  ostr << "offset";
}

ProcessFill::ProcessFill(double valv, complexity_t filltypev, const BaseList<size_t>& col_selv, const BaseList<size_t>& row_selv)
  : ProcessCommand(PROC_HAS2D | (row_selv.empty() ? PROC_HAS1D : 0)),
    val(valv), filltype(filltypev),
    col_sel(col_selv), row_sel(row_selv) {}

void ProcessFill::printvariablename(std::ostream& ostr, subsid_t) const 
{
  ostr << "fill";
}

//ContextWarning<> ProcessFill::fillwillwipe_warning("fill without range or real/imaginary qualifier will wipe entire data set",&NMRsim_repeat_warning);

ContextWarning<> disallowpair_warning("-complexpair is not meaningful in this context (treating as -complex)",&NMRsim_once_warning);

ProcessCommand::complexity_t ProcessCommand::parsecomplexity(complexity_t def, bool allowpair)
{
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["imag"]=IMAG;
    flags["real"]=REAL;
    flags["complex"]=COMPLEX;
    flags["complexpair"]=COMPLEXPAIR;
  }
  complexity_t res=(complexity_t)parse_flags(flags,(size_t)def,true);
  if ((res==COMPLEXPAIR) && !allowpair) {
    disallowpair_warning.raise();
    res=COMPLEX;
  }
  return res;
}

ProcessCommand* ProcessFill::create()
{  
  Variable cvar(S_ARG1);
  const double val=parse_double(&cvar,defflags());
  //  const bool isconst=(cvar.subsid==S_NONE);
  //  VariableBase* varp=parse_double_variable(cvar,F_ALLOWLIST);
  LIST<size_t> col_sel,row_sel;
  if (parser_isnormal()) {
    col_sel=parse_unsignedintarray(1);
    if (parser_isnormal())
      row_sel=parse_unsignedintarray(1);
  }
  const complexity_t flagsval=parsecomplexity(COMPLEX);
//   if ((flagsval==COMPLEX) && col_sel.empty())
//     fillwillwipe_warning.raise();
  return new ProcessFill(val,flagsval,col_sel,row_sel);
}

void ProcessFill::print(std::ostream& ostr) const
{
  ostr << "fill " << val << ' ';
  if (!col_sel.empty()) {
    ostr << ' ' << col_sel;      
    if (!(row_sel.empty()))
      ostr << ' ' << row_sel;
  }
}

void ProcessFill::exec(BaseList<complex> FID, processing_state&) const
{
  exec_(FID);
}

void ProcessFill::exec_(BaseList<complex>& FID) const
{
  const complex cmplxval(val,val);
  if (col_sel.empty()) {
    assert(row_sel.empty()); //!< harmless assert
    switch (filltype) {
    case REAL: reals(FID)=val; break;
    case IMAG: imags(FID)=val; break;
    case COMPLEX: case COMPLEXPAIR: FID=cmplxval; break;
    }
  }
  else {
    for (size_t i=col_sel.size();i--;) {
      complex& cval(FID(col_sel(i)));
      switch (filltype) {
      case REAL: real(cval,val); break;
      case IMAG: imag(cval,val); break;
      case COMPLEX: case COMPLEXPAIR: cval=cmplxval; break;
      }
    }
  }
}

void ProcessFill::exec(cmatrix& FID, LIST<processing_state>&) const
{
  for (size_t i=row_sel.empty() ? FID.rows() : row_sel.size();i--;) {
    BaseList<complex> FIDr(FID.row(row_sel.empty() ? i : row_sel(i)));
    exec_(FIDr);
  }
}

size_t ProcessCommand::checkcomplexity(size_t n, complexity_t type)
{
  if ((type!=COMPLEX) && (type!=COMPLEXPAIR))
    return n;
  if (n % 2) {
    parser_printthread(std::cerr) << "expression returning supposedly complex data returned list with odd number of elements (" << n << ")\n";
    error_abort();
  }
  return n/2;
}

namespace {
  bool isnoncontig(const BaseList<size_t>& sel)
  {
    const size_t n=sel.size();
    if (n>1) {
      size_t prev=sel.front();
      for (size_t i=1;i<n;i++) {
	size_t cur=sel(i);
	if (cur!=(prev+1)) //!< has to be increasing, referencing will still be messed up if reversed
	  return true;
	prev=cur;
      }
    }
    return false;
  }
      
}

ProcessExtract::ProcessExtract(const BaseList<size_t>& col_selv, const BaseList<size_t>& row_selv)
  : ProcessCommand(PROC_HAS2D | (row_selv.empty() ? PROC_HAS1D : 0)),
    col_sel(col_selv), row_sel(row_selv),
    col_noncontiguous(isnoncontig(col_selv)),
    row_noncontiguous(isnoncontig(row_selv))
{
  if (!(col_selv.empty()))
    flags_|=PROC_CHANGES_COLUMNS;
  if (!(row_selv.empty()))
    flags_|=PROC_CHANGES_ROWS;
}

ContextWarning<> ProcessExtract::nullselection_warning("extract: null selection will do nothing - directive ignored",&NMRsim_once_warning);

ThreadWarning<> ProcessExtract::refreset_warning("extract: selection is non-contiguous, so referencing has been reset to give zero at mid spectrum. Use set -ref to explicitly set reference",&NMRsim_repeat_warning);
//ThreadWarning<> ProcessReverse::refreset_warning("rev: referencing has been reset to give zero at mid spectrum. Use set -ref to explicitly set reference",&NMRsim_repeat_warning);
ThreadWarning<> ProcessResample::refreset_warning("resample: referencing has been reset to give zero at mid spectrum. Use set -ref to explicitly set reference",&NMRsim_repeat_warning);
//ThreadWarning<> ProcessShiftHalf::refreset_warning("shifthalf: referencing has been reset to give zero at mid spectrum. Use set -ref to explicitly set reference",&NMRsim_repeat_warning);

ProcessCommand* ProcessExtract::create()
{
  const LIST<size_t> col_sel(parse_unsignedintarray(1));
  if (are_left()) {
    const LIST<size_t> row_sel(parse_unsignedintarray(1));
    if (row_sel.empty() && col_sel.empty()) {
      nullselection_warning.raise();
      return NMRSIM_NULL;
    }
    checkrowchangeok();
    return new ProcessExtract(col_sel,row_sel);
  }
  if (col_sel.empty()) {
    nullselection_warning.raise();
    return NMRSIM_NULL;
  }
  return new ProcessExtract(col_sel);
}

void ProcessExtract::print(std::ostream& ostr) const
{
  ostr << "extract " << col_sel;
  if (!(row_sel.empty()))
    ostr << ' ' << row_sel;
}

void ProcessExtract::scaleswref(processing_state& pflags, size_t startn, size_t newn, size_t oldn, bool isnoncontig)
{
  if (!(pflags.istimedomain)) {
    if (newn==1) { //!< if single point, set sw to 0 and shift reference to frequency of point 
      pflags.ref += pflags.sw*(double(startn)/oldn-0.5);
      pflags.sw =0.0;
    }
    else {
      if (pflags.ref) {
	if (isnoncontig) {
	  refreset_warning.raise();
	  pflags.ref=0;
	}
	else
	  pflags.ref += pflags.sw*( (startn+newn*0.5)/oldn-0.5);
      }
      pflags.sw*=double(newn)/oldn;
    }
  }
}

void ProcessExtract::exec(List<complex>& dest, const BaseList<complex>& source, processing_state& pflags) const
{
  //  pflags.invalidatescale();
  const size_t oldnp=source.size();
  assert(row_sel.empty()); //!< harmless assert
  dest=source(col_sel);
  scaleswref(pflags,col_sel.front(),dest.size(),oldnp,col_noncontiguous);
}

ProcessCommand* ProcessAddNoise::create()
{
  Variable cvar(S_ARG1);
  const double sigmav=parse_double(&cvar,defflags());

  static flagsmap_type flags;
  if (flags.empty())
    flags["reseed"]=1;

  return new ProcessAddNoise(sigmav,parse_flags(flags));
}

ProcessCommand* ProcessScale::create()
{
  Variable cvar(S_ARG1);
  const double scale=parse_double(&cvar,F_DENYZERO | defflags());
  if (cvar.subsid && (cvar.variable().uses() & A_VAR))
    havescale=true;    
  return new ProcessScale(scale);
}

ProcessCommand* ProcessScaleFirst::create()
{
  Variable cvar(S_ARG1);
  const int flags=defflags();
  const double scalef=parse_double(&cvar,flags);
  Variable cvar1(S_ARG2);
  const double scalef1=(process2D && are_left()) ? parse_double(&cvar1,flags) : 1.0;
  return new ProcessScaleFirst(scalef,scalef1);
}

ProcessAddLB::ProcessAddLB(double lb_,double gfrac_, double lb1_, double gfrac1_)
  : ProcessCommand(PROC_HAS1D | PROC_HAS2D),
    gfrac(gfrac_), gfrac1(gfrac1_)
{
  set(lb_,S_ARG1);
  set(lb1_,S_ARG2);
}

void ProcessAddLB::exec(cmatrix& FID, LIST<processing_state>& pflags) const
{
  if (!pflags.back().istimedomain) {
    if (t21)
      error_abort("2D frequency domain version of addlb not (yet?) implemented");
    for (size_t i=FID.rows();i--;) {
      processing_state pstate(pflags.back());
      exec(FID.row(i),pstate);
    }
    return;
  }
//  static bool donewarn=!allowwarnings();
//   if (!pflags.istimedomain && !donewarn) {
//     std::cerr << "Warning (first): addlb applied to frequency domain data\n";
//     donewarn=true;
//   }
  const double lsw=pflags.back().sw;
  if (lsw==0.0)
    error_abort("Can't add linebroadening if spectral width unset");
  const double t2max=FID.cols()/lsw;
  const double sw1=pflags.front().sw;
  const double lbl1= t21 ? FID.rows()/(t21*sw1) : 0.0;
  const size_t ni_skip=skips.front();
  exponential_multiply_ip(FID,lbl1, t2 ? t2max/t2 : 0.0, ni_skip);
  if (lbg || lbg1)
    gaussian_multiply_ip(FID,lbg1,lbg1 ? 1.0/sw1 : 0.0, lbg,1.0/sw, ni_skip);
}

template<class T1,class T2,class T3> void mlaroll(BaseList<T1> dest, T2 scale, const BaseList<T3>& source, size_t roll)
{
  const size_t nbins=dest.size();
  if (nbins!=source.size())
    throw Mismatch("mlaroll");
  if (roll==0) {
    mla(dest,scale,source);
    return;
  }
  if (roll>=nbins)
    throw InvalidParameter("mlaroll: roll");
  BaseList<T1> dest1(dest(range(roll,nbins-1)));  
  mla(dest1,scale,source(range(0U,nbins-roll-1)));
  BaseList<T1> dest2(dest(range(0U,roll-1)));
  mla(dest2,scale,source(range(nbins-roll,nbins-1)));
}

ThreadWarning<> ProcessAddLB::narrow_warning("linewidth in addlb is narrower than histogram bin width - results may be unreliable",&NMRsim_once_warning);

void ProcessAddLB::execf(BaseList<complex>& FID) const
{
  const size_t lnp=FID.size();
  if (lb<(sw/lnp))
    narrow_warning.raise();

  if (!shapegenp)
    shapegenp.reset(new LorentzianGaussian(lb,gfrac));
  else {
    shapegenp->lw(lb);
    shapegenp->fraction(gfrac);
  }
  (void)(*shapegenp)(shapebuf,HistogramSpec(lnp,sw));
  if ((verbose & VER_GEN) && (verbose_level>1))
    std::cout << "Addlb lineshape: " << shapebuf << '\n';
  static LIST<complex> tmp;
  static const complex complexzero(0.0);
  tmp.create(lnp,complexzero);
  for (size_t i=lnp;i--;) {
    if (FID(i)!=complexzero)
      mlaroll(tmp,FID(i),shapebuf,i);
  }
  FID=tmp;
}

void ProcessAddLB::exec(BaseList<complex> FID, processing_state& pflags) const
{
  if (t21)
    error_abort("2D version of addlb inapplicable to 1D data");
  if (sw==0.0)
    throw Failed("Can't add linebroadening if spectral width unset");
  if (!pflags.istimedomain) {
    execf(FID);
    return;
  }
  if (verbose & VER_GEN)
    std::cout << "Add line-broadening with T2: " << t2 << " s  Guass frac: " << lbg << '\n';
  const double t2max=FID.size()/sw;
  exponential_multiply_ip(FID,t2 ? t2max/t2 : 0.0);
  if (lbg)
    gaussian_multiply_ip(FID,lbg,1.0/sw);
}

void ProcessAddNoise::print(std::ostream& ostr) const
{
  ostr << "addnoise " << sigma;
  if (reseed)
    ostr << " -reseed";
}

void ProcessAddNoise::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "addnoise";
}

void ProcessAddNoise::exec(BaseList<complex> FID, processing_state&) const
{
  if (reseed)
    set_seed();
  for (size_t i=FID.size();i--;)
    FID(i)+=complex(gauss(sigma),gauss(sigma));
}

void ProcessAddLB::print_(std::ostream& ostr, double llb, double lfrac)
{
  ostr << llb << " Hz (" << int(100.0*(1.0-lfrac)) << "% Lorentzian)";
}

void ProcessAddLB::print(std::ostream& ostr) const
{
  ostr << "addlb ";
  print_(ostr,lb,gfrac);
  if (lb1) {
    ostr << ' ';
    print_(ostr,lb1,gfrac1);
  }
}

void ProcessAddLB::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "addlb_";
  switch (subsid) {
  case S_ARG2:
    ostr << "lb1";
    break;
  case S_GFRAC1:
    ostr << "gauss_frac1";
    break;
  case S_ARG1:
    ostr << "lb";
    break;    
  case S_GFRAC:
    ostr << "gauss_frac";
    break;
  default:
    throw InvalidParameter("ProcessAddLB::set");
  }
}

void ProcessAddLB::set_(double& lb_, double& gfrac_, double& t2_, double& lbg_, double lbin, double gfracin)
{
  char buf[256];
  if ((gfracin<0.0) || (gfracin>1.0)) {
    snprintf(buf,sizeof(buf),"addlb: Lorentzian/Gaussian fraction must be between 0 and 1 (given %g)",gfracin);
    throw InvalidParameter(buf);
  }
  lb_=lbin;
  gfrac_=gfracin;
  double tmplb=lbin*(1.0-gfracin);
  t2_= tmplb ? (1.0/(M_PI*tmplb)) : 0.0;
  lbg_ = lbin*gfracin; 
}

namespace {
  inline double mod1(double x)
  {
    //    if (x<0)
      //      return 1.0+x+floor(-x);
    //if (x>1)
      // return x-floor(x);
    return x;
  }
}

void ProcessAddLB::set(double val, subsid_t subsid)//, bool is2D)
{
//   if (is2D)
//     std::cerr << warn_proc_in_2D;

  switch (subsid) {
  case S_ARG1:
    set_(lb,gfrac,t2,lbg,fabs(val),gfrac);
    break;
  case S_ARG2:
    set_(lb1,gfrac1,t21,lbg1,fabs(val),gfrac1);
    break;
  case S_GFRAC:
    set_(lb,gfrac,t2,lbg,lb,mod1(val));
    break;
  case S_GFRAC1:
    set_(lb1,gfrac1,t21,lbg1,lb1,mod1(val));
    break;
  default:
    throw InternalError("ProcessAddLB: invalid variable specifier");
  }
}

ProcessCommand* ProcessAddLB::create()
{
  const int parseflags=defflags();  
  Variable cvar(S_ARG1);
  const double lb=parse_double(&cvar,parseflags);
  double gfrac=0.0;
  double gfrac1=0.0;
  double lb1=0.0;
  if (are_left()) {
    cvar.subsid=S_GFRAC;
    gfrac=parse_double(&cvar,parseflags);
    if (process2D && are_left()) {
      cvar.subsid=S_ARG2;
      lb1=parse_double(&cvar,parseflags);
      cvar.subsid=S_GFRAC1;
      gfrac1=parse_double(&cvar,parseflags);
    }
  }
  return new ProcessAddLB(lb,gfrac,lb1,gfrac1);
}

int roundup2(int n)
{
  if (n<=0)
    throw InvalidParameter("roundup2");
  size_t newn=1;
  size_t tmpn=n;
  while ((tmpn>>=1))
    newn<<=1;
  return (newn<n) ? newn<<1 : newn;
}

int zerofillpts(int newpts, int curpts)
{
  if (newpts<curpts) {
    if (newpts>8)
      error_abort("zerofill argument (<new points | fill factor>) is ambiguous");
    if (newpts==0)
      return curpts;
    return roundup2(newpts*curpts);    
  }
  return newpts;
}

ProcessZeroFill::ProcessZeroFill(int npts_, int npts1_)
  : ProcessCommand(PROC_HAS2D | PROC_CHANGES_COLUMNS | (((npts1_==1) && !active2D) ? PROC_HAS1D : 0)),
    npts(npts_), npts1(npts1_) 
{
  if ((npts<0) || (npts1<0))
    error_abort("zerofill: factor/points cannot be negative");
  if (npts1_>1)
    flags_|=PROC_CHANGES_ROWS;
}

ThreadWarning<> ProcessZeroFill::ignoring_warning("zerofill in indirect dimension ignored for irregular data",&NMRsim_once_warning);

void ProcessZeroFill::exec(List<complex>& dest, const BaseList<complex>& source, processing_state&) const
{
  if (npts1)
    ignoring_warning.raise();
  const int curpts=source.size();
  const int newpts=zerofillpts(npts,curpts);
  if (newpts==curpts)
    dest=source;
  else {
    dest.create(newpts,complex(0.0));
    const range colsel(0,curpts-1);
    BaseList<complex> drow(dest(colsel));
    drow=source;
  }
}

void
ProcessZeroFill::exec(cmatrix& FID, LIST<processing_state>&) const
{
  const int curpts=FID.cols();
  const int newpts=zerofillpts(npts,curpts);
  const int curpts1=FID.rows();
  const int newpts1=zerofillpts(npts1,curpts1);
  if ((newpts==curpts) && (newpts1==curpts1))
    return; //do nothing

  cmatrix FIDt(newpts1,newpts,complex(0.0));
  const range colsel(0,curpts-1);
  for (size_t nr=curpts1;nr--;) { //avoid use of IndirectMatrix (inefficient)
    BaseList<complex> drow(FIDt.row(nr));
    drow(colsel)=FID.row(nr);
  }
  FID.swap(FIDt);
}

void
ProcessZeroFill::print(std::ostream& ostr) const
{
  ostr << "zerofill " << npts;
  if (npts1)
    ostr << ' ' << npts1;
}

ProcessCommand*
ProcessZeroFill::create()
{
  int nz=2;
  int nz1=0;
  if (are_left()) {
    nz=parse_int();
    if (process2D && are_left())
      nz1=parse_int();
  }
  return new ProcessZeroFill(nz,nz1);
}

ContextWarning<> ProcessCommand::fidspectrum_warning("-fid/-spectrum flags are undocumented and should be avoided. Niormally use ft for FID->spectrum and ft -inv for spectrum -> FID.",&NMRsim_once_warning);

size_t ProcessCommand::get_ftflags()
{
  static flagsmap_type FTflags;
  if (FTflags.empty()) {
    FTflags["inv"]=FT_INV;
    FTflags["noscalefirst"]=FT_NOFIRST;
    FTflags["noshift"]=FT_NOSHIFT;
    FTflags["fid"]=FT_FID;
    FTflags["spectrum"]=FT_SPECTRUM;
  }
  size_t flags=parse_flags(FTflags);
  
  static const int dirmask=(FT_INV | FT_FID | FT_SPECTRUM);
  switch (flags & dirmask) {
  case 0:
    break;

  case FT_FID:
    fidspectrum_warning.raise();
    flags &= ~FT_FID; //!< clear -fid flag to leave normal "forward" FT
    break;
    
  case FT_SPECTRUM:
    fidspectrum_warning.raise();
    flags &= ~FT_SPECTRUM;
    flags |= FT_INV; //!< replace -spectrum with inverse FT
    break;

  default:
    error_abort("ft: can only specify one of -inv, -fid or -spectrum");
  }
  return flags;
}
  
ProcessCommand*
ProcessFT::create()
{
  const size_t flags=get_ftflags();
//   int dirtype=0;
//   static const int dirmask=(FT_INV | FT_FID | FT_SPECTRUM);
//   switch (flags & dirmask) {
//   case 0: case FT_SPECTRUM:
//     break;
//   case FT_INV: case FT_FID:
//     dirtype=FT_INV;
//     break;
//   default:
//     error_abort("ft: can only specify one of -inv, -fid or -spectrum");
//   }  
//   flags = (flags & ~dirmask) | dirtype;
  return new ProcessFT(flags);
}

void ProcessFT::print(std::ostream& ostr) const
{
  ostr << "ft";
  print_ftflags(ostr,ftflags);
}

ProcessCommand* ProcessPhase::create()
{
  Variable cvar(S_ARG1);
  const double first=parse_double(&cvar,defsumflags() | F_ISPHASE);
  cvar.subsid=S_ARG2;
  const double zero=parse_double(&cvar,defsumflags() | F_ISPHASE);
  return new ProcessPhase(zero,first,0.0); //pivot around zero...
}

ProcessFT2d::ProcessFT2d(size_t flagsv)
  : ProcessCommand(PROC_HAS2D),
    ftflags(flagsv), needphase(false)
{ create_(); }

ProcessFT2d::ProcessFT2d(double zero2_, double first2_, double zero1_, double first1_, size_t ftflagsv)
  : ProcessCommand(PROC_HAS2D), 
    ftflags(ftflagsv), needphase(true), zero2(zero2_), first2(first2_), zero1(zero1_), first1(first1_)
{ create_(); }

ThreadWarning<> ProcessFT2d::invalidskip_warning("ft2d only understands t1 steps of 1 (phase-modulated, default) or 2 (hypercomplex)",&NMRsim_once_warning);

void ProcessFT2d::create_()
{
  if (skips.front()>2)
    invalidskip_warning.raise();
}

ProcessCommand*
ProcessFT2d::create()
{
  if (!process2D)
    error_abort("ft2d only valid for 2D data",ERR_FAILED);

  if (parser_isnormal()) {
    if (!active2D)
      error_abort("ft2d with phasing only valid for true 2D data (ni unset)");
    const double first2=parse_double();
    const double zero2=parse_double();
    const double first1=parse_double();
    const double zero1=parse_double();
    return new ProcessFT2d(zero2,first2,zero1,first1,get_ftflags());
  }
  return new ProcessFT2d(get_ftflags());
}

ContextWarning<> ProcessCommand::rowchange_warning("Processing step changing numbers of rows may cause problems in combination with {} arrays (suppress warning with -nochecks)",&NMRsim_repeat_warning);

void ProcessCommand::checkrowchangeok()
{
  if (array_n0<1) 
    throw InternalError("checkrowchangeok");
  if (!nochecks && !active2D && (array_n0>1))
    rowchange_warning.raise();
}

ProcessCommand* ProcessTranspose::create()
{
  checkrowchangeok();
  return new ProcessTranspose();
}

void ProcessTranspose::rawexec(DataStore& FID) const
{
  if (skips.front()!=1)
      error_abort("transpose not implemented for hypercomplex data (or other step factors <>1)");
  if (!(FID.transpose()))
    error_abort("transpose not possible e.g. data rows are of inconsistent length or parameters (e.g. spectral width");    
}

//Warning<> statemismatch_warning("transpose may be invalid: processing state (e.g. spectral width) inconsistent between rows");

bool DataStore::transpose()
{
  if (isnD()) {
    if (procstates_.size()>1)
      std::swap(procstates_.front(),procstates_(size_t(1))); //!< swap dimension information
    matrix().transpose(); //!< Not exception safe
    return true;
  }
  size_t lnp=row(0).size();
  for (size_t i=1;i<rows();i++) {
    if ((lnp!=row(i).size()) || (procstates_(i)!=procstates_.front()))
      return false;
  }

  List<size_t> sizes(int(lnp),rows());
  ListList<complex>& aslist(listlist());
  ListList<complex> newlist(sizes,aslist.row());
  aslist.swap(newlist);
  return true;
}
 
bool processing_state::operator== (const processing_state& a) const
{
  return (sw==a.sw) && (istimedomain==a.istimedomain) && (sfrq==a.sfrq) && (ref==a.ref);
}

ThreadWarning<> ProcessCommand::varignored_warning("variable quantity used in nD processing (ignore if simple summation variable)?",&NMRsim_once_warning);

ThreadWarning<> ProcessCommand::oneDprocessing_warning("Processing is being applied to 2D data set with 1 data point per row, which may not give expected results. Consider applying transpose first (disable warning with -nochecks).",&NMRsim_once_warning);

void ProcessCommand::rawexec(DataStore& FID) const
{
  bool hasvars=(flags_ & PROC_VARIABLE);
  if (hasvars) {
    if (active2D && (FID.rows()>1))
      varignored_warning.raise();
    
    if (!(flags_ & PROC_HAS1D))
      throw InternalError("Processing method has a variable parameter, but not a 1D method");
    //     if (active2D) {
    //       print(std::cerr);      
    //       error_abort(": nD processing cannot involve variable parameters");
    //     }
  }
  const bool isnD=FID.isnD();

  if (isnD && !hasvars && (flags_ & PROC_HAS2D)) {
    if (!nochecks && (FID.matrix().cols()==1) )
      oneDprocessing_warning.raise();
    //! NB doesn't set current_data_row
    exec(FID.matrix(),FID.states());
    if (flags_ & PROC_CHANGES_ROWS)
      //      ensure_array(); //!< not clear this does anything useful as commands don't change array_n0 etc. ? 16/10/2013
      update_array_loop(FID.rows());
    return;
  }

  array_iter iter(hasvars);
  if (iter.size()!=FID.rows()) {
    std::cout << "Number of rows in data set (" << FID.rows() << ") does not match number of array steps (" << iter.size() << ")\n";
    error_abort();
  }
  processing_state statetmp(0.0,false,0.0); //!< dummy
  processing_state& state0(FID.states().back());
  if (flags_ & (PROC_CHANGES_ROWS | PROC_CHANGES_COLUMNS)) {
    List<complex> rtmp;
    ListList<complex> tmp;
    cmatrix tmpm;

    while (iter.next(&FID)) {
      statetmp=state0; //!< need to copy 
      processing_state& pstate(isnD ? statetmp : FID.state(row_index));
      current_data_row.create(FID.row(row_index));
      exec(rtmp,current_data_row,pstate);
      if (isnD) {
	if (!tmpm)
	  tmpm.create(FID.rows(),rtmp.size());
	else {
	  if (rtmp.size()!=tmpm.cols())
	    error_abort("can't combine changing data size with nD data set");
	}
	tmpm.row(row_index)=rtmp;
      }
      else	    
	tmp.push_back(rtmp);
      if (hasvars)
	save_parameters();
    }
    if (isnD)
      FID.matrix().swap(tmpm);
    else
      FID.listlist().swap(tmp);
  }
  else {
    while (iter.next(&FID)) {
      statetmp=state0; //!< need to copy 
      processing_state& pstate(isnD ? statetmp : FID.state(row_index));
      exec(current_data_row,pstate);
      if (hasvars)
	save_parameters();
    }
  }
  if (isnD)
    state0=statetmp;
}

//! set PROC_VARIABLE is involves non-const var - does not clear it if set explicitly
void ProcessCommand::isconstant(bool isconst)
{
//   if (isconst)
//     flags_&=~PROC_VARIABLE;
//   else
//     flags_|=PROC_VARIABLE;
  if (!isconst)
    flags_|=PROC_VARIABLE;
}

ThreadWarning<> process_nD_warning("nD processing not fully supported...",&NMRsim_once_warning);

void verify_have_FID(const DataStore& FID)
{
  if (FID.empty()) {
    std::cerr << "No FID to process! Need to set at least np and sw." << std::endl;
    error_abort();
  }
}

void procstack_t::exec(DataStore& FID) const
{
  const const_iterator iend(this->end());
  procstack_t::const_iterator start(this->begin());
  if (start==iend)
    return;
  verify_have_FID(FID);
  //  const bool updatevars = (!active2D && !array_ns.empty()) || (optupdate==0);

  const bool isnD=FID.isnD();
  if (isnD && (array_ns.size()>2))
    process_nD_warning.raise();
  const bool verb=(verbose & VER_GEN) && (verbose_level>1);
  if (verb) {
    std::cout << "Initial data structure:\n";
    FID.print_structure();
    std::cout << '\n';
  }
  size_t step=1;
  do {
    const ProcessCommand* procp(*start);
    timer_guard_t guard(procp->stopwatch);
    current_data_row.clear(); //!< reclear for each command - responsibility of rawexec to maintain
    procp->rawexec(FID);
    ++start;
    if (verb) {
      std::cout << "Data structure after processing step " << step << ":\n";
      FID.print_structure();
      std::cout << '\n';
    }
    step++;
  } while (start!=iend);
  current_data_row.clear(); //!< ensure no current data
}

// template<class T> class dataswapper
// {
// public:
//   dataswapper(const BaseList<T>& sourcev)
//     : orig_(sourcev) { cptr_=&orig_; bakptr_=&bak_; }
  
// //   void cleanup() { 
// //     if (cptr_!=&orig_)
// //       orig_.swap(bak_);
// //   }
//   LIST<T>& workspace() { return *bakptr_; }
//   LIST<T>& current() { return *cptr_; }
//   void swap() { std::swap(cptr_,bakptr_); }

// private:
//   LIST<T> orig_,bak_;
//   LIST<T>* cptr_,bakptr_;
// };

void apply_procstacks(const BaseList<procstack_t>& stacks, DataStore& data)
{
  const size_t nstacks=stacks.size();
  switch (nstacks) {
  case 0: return;
  case 1:
    stacks.front().exec(data);
    return;
  }
  verify_have_FID(data);
  if (data.isnD())
    error_abort("can't apply multiple processing blocks to nD data");

  array_iter iter(true);
  ListList<complex>& rdata(data.listlist());
  ListList<complex> tmp;
  List<complex> rtmp1,rtmp2;
  List<complex>* workspacep(&rtmp1);
  List<complex>* currentp(&rtmp2);
  while (iter.next()) {
    const size_t curind=row_index % nstacks;
    const procstack_t& curstack(stacks(curind));
    const procstack_t::const_iterator iend(curstack.end());
    procstack_t::const_iterator start(curstack.begin());
    processing_state& cflags(data.state(row_index));
    *currentp=rdata(row_index);
    save_parameters();
    const bool verb=(verbose & VER_GEN) && (verbose_level>1);
    if (verb)
      std::cout << "Row " << (row_index+1) << ":\nInitial data state: " << cflags << '\n';
    size_t step=1;
    while (start!=iend) {
      const ProcessCommand* procp(*start);
      timer_guard_t guard(procp->stopwatch);
      //      const bool hasvars=(procp->flags_ & PROC_VARIABLE);
      if ((procp->flags()) && (ProcessCommand::PROC_CHANGES_ROWS | ProcessCommand::PROC_CHANGES_COLUMNS)) {
	procp->exec(*workspacep,*currentp,cflags);
	std::swap(workspacep,currentp);
      }
      else	
	procp->exec(*currentp,cflags);
      if (verb)
	std::cout << "Data state after processing step " << step << ": " << cflags << '\n';
      step++;
      ++start;      
    } while (start!=iend);
    tmp.push_back(*currentp);
  }
  rdata.swap(tmp);
}

void ProcessCommand::exec(BaseList<complex>, processing_state&) const
{
  parser_printthread(std::cerr) << "Processing method did not supply a (static) 1D method: ";
  print(std::cerr);
  error_abort();
}

void ProcessCommand::exec(List<complex>&, const BaseList<complex>&, processing_state&) const 
{
  parser_printthread(std::cerr) << "Processing method did not supply a (dynamic) 1D method: ";
  print(std::cerr);
  error_abort();
}

void ProcessCommand::exec(Matrix<complex>&, LIST<processing_state>&) const
{
  parser_printthread(std::cerr) << "Processing method did not supply a 2D method: ";
  print(std::cerr);
  error_abort();
}

void ProcessCommand::set(double,subsid_t)
{ throw InternalError("This function is not setable"); }

void ProcessCommand::printvariablename(std::ostream&, subsid_t) const
{ throw InternalError("This function is not setable and so has no printvariablename"); }

void ProcessSignReverse::exec(BaseList<complex> data, processing_state& pflags) const
{
  if (ifneggamma_ && (detect_freq>=0.0))
    return;
  if (pflags.istimedomain)
    ProcessConjugate::exec(data);
  else
    ProcessReverse::exec(data);
}

void ProcessSignReverse::print(std::ostream& ostr) const 
{ 
  ostr << "sign_reverse";
  if (ifneggamma_)
    ostr << " -if_negative_gamma";
}

void ProcessAddSignals::print(std::ostream& ostr) const
{
  ostr << "addsignals ";
  freqs_->print(ostr);
  ostr << ' ';
  amps_->print(ostr);
  if (!!phases_) {
    ostr << ' ';
    phases_->print(ostr);
  }
}

void ProcessAddSignals::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "addsignals_";
  switch (subsid) {
  case S_OFFSET: ostr << "frequency"; return;
  case S_ARG1: ostr << "intensity"; return;
  case S_PHASE: ostr << "phase"; return;
  default:
    throw InternalError("AddSignals::printvariablename: unknown subsid");
  }
}

namespace {
  template<class T> inline T getlist(const BaseList<T>& a, size_t n) {
    return (a.size()==1) ? a.front() : a(n);
  }
}

//Warning<> ProcessAddSignals::nothistogram_warning("using histogram for frequency domain addsignals, but not in histogram mode!",&NMRsim_once_warning);
ThreadWarning<> ProcessAddSignals::ignoring_phase_warning("non-zero phase parameter has been ignored (frequency domain addsignals)",&NMRsim_repeat_warning);
ThreadWarning<> ProcessAddSignals::timedomainref_warning("non-zero reference frequency is being ignored for time-domain addsignals. Explicitly subtract the reference frequency if required.",&NMRsim_once_warning);

//ThreadWarning<> ProcessAddSignals::histogramloss_warning("addsignals lost:",&NMRsim_once_warning);

void ProcessAddSignals::exec(BaseList<complex> FID, processing_state& pflags) const
{
  if (pflags.sw==0)
    error_abort("Can't apply addsignals when spectral width is undefined");
//   if (!pflags.istimedomain && !errwarn) {
//     std::cerr << "Warning (first): addsignals applied to frequency domain data\n";
//     errwarn=true;
//   }
  const BaseList<double> amps(amps_->value());
  const BaseList<double> freqs(freqs_->value());
  accumulator acclength;
  bool ok= acclength.add(amps.size()) && acclength.add(freqs.size());
  BaseList<double> phases;
  if (!!phases_) {
    phases.create(phases_->value());
    ok&=acclength.add(phases.size());
  }
  //  const double ref=pflags.ref;
  if (!ok)
    error_abort("addsignals: lists have different lengths");
  if (pflags.istimedomain) {
    if (pflags.ref) 
      timedomainref_warning.raise();

    for (size_t i=acclength();i--;) {
      const double fact=getlist(freqs,i);
      const double f=-fact/pflags.sw; //!< reverse sign
      const double a=getlist(amps,i);
      if (phases.empty()) {
	if (verbose & VER_GEN)
	  std::cout << "Adding line of frequency " << (fact*1e-3) << " kHz and amplitude " << a << '\n';
	add_FID(FID,complex(a),f);
      }
      else {
	const double p=getlist(phases,i);
	if (verbose & VER_GEN)
	  std::cout << "Adding line of frequency " << (fact*1e-3) << " kHz, amplitude " << a << " and phase " << p << " degrees\n";
	add_FID(FID,a*expi(p*deg_to_rad),f);
      }
    }
  }
  else {
    //const bool allowwarn=allowwarnings();
//     if (!isfrequency())
//       nothistogram_warning.raise();
    HistogramMaker hist(FID,pflags.sw);
    bool needwarn=false;
    for (size_t i=acclength();i--;) {
      const double f=getlist(freqs,i)-pflags.ref;
      const double a=getlist(amps,i);
      const double p=phases.empty() ? 0.0 : getlist(phases,i);
      if (p)
	needwarn=true;
      hist().add(complex(a),f); //!< ignore phase because lineshape is purely real and cannot be phased
    }
    if (needwarn)
      ignoring_phase_warning.raise();
    const BaseHistogram<complex>& acthist(hist());
    if (!nochecks && acthist.lost()) {
      char errstr[1024];
      const complex lostsum(acthist.lostsum());
      snprintf(errstr,sizeof(errstr),"%" LCM_PRI_SIZE_T_MODIFIER "u frequency(ies) with total amplitude (%g,%g)",acthist.lost(),real(lostsum),imag(lostsum));
      HistogramMaker::loss_warning.raise(errstr);
    }
  }
}

ProcessAddSignals::ProcessAddSignals()
  : ProcessCommand(PROC_HAS1D),
    errwarn(false)
{}

ContextWarning<> ProcessAddSignals::allzero_warning("amplitudes are all zero (or missing).  Note that argument order is <freq> <amp> <phase>.",&NMRsim_once_warning);

bool areallzero(const BaseList<double>& amps)
{
  for (size_t i=amps.size();i--;) {
    if (amps(i))
      return false;
  }
  return true;
}

const char* AddSignals_syntax="<list of frequencies>#<list of amplitudes>#[<list of phases>]";

ProcessCommand* ProcessAddSignals::create()
{
  //Mark markobj; //!< should be redundant since const status checked automatically (11/11/09)
  Variable cvar(S_OFFSET);
  const int flags=defsumflags() | F_ALLOWLIST;
  ProcessAddSignals* newp=new ProcessAddSignals();
  newp->freqs_.reset(parse_double_variable_syntax(cvar,AddSignals_syntax,1,flags));
  cvar.subsid=S_ARG1;
  newp->amps_.reset(parse_double_variable_syntax(cvar,AddSignals_syntax,2,flags));
  if (cvar.subsid && (cvar.variable().uses() & A_VAR))
    havescale=true;    
  if (newp->amps_->isconst() && areallzero(newp->amps_->value()))
    allzero_warning.raise();
  if (are_left()) {
    cvar.subsid=S_PHASE;
    newp->phases_.reset(parse_double_variable_syntax(cvar,AddSignals_syntax,3,flags | F_ISPHASE));
  }
//   if (markobj.needsflush()) {
//     markobj.flush(newp);
//     newp->isconstant(false);
//   }
  return newp;
}

template<class T> void ProcessAddSignals::set_(const T& val, subsid_t subsid)
{
  switch (subsid) {
  case S_OFFSET:
    freqs_->set(val);
    break;
  case S_ARG1:
    amps_->set(val);
    break;
  case S_PHASE:
    phases_->set(val);
    break;
  default:
    throw InternalError("ProcessAddSignals::set");
  }
}

void ProcessAddSignals::set(double val, subsid_t subsid)
{
  set_(val,subsid);
}

void ProcessAddSignals::set(const BaseList<double>& val, subsid_t subsid)
{
  set_(val,subsid);
}


namespace {
  struct Add_ip {
    template<typename T> void operator()(T& a, const T& b) const { a+=b; }
  };

  struct Subtract_ip {
    template<typename T> void operator()(T& a, const T& b) const { a-=b; }
  };

  struct Assign {
    Assign(const complex& v) : v_(v) {}
    template<typename T> void operator()(T& a) const { a=v_; }
    complex v_;
  };
}

const BaseList<complex> DataStore::row() const
{
  switch (type()) {
  case 1:
    return (*this)(Int2Type<1>()).row();
  case 2:
    return (*this)(Int2Type<2>()).row();    
  }
  throw Undefined("DataStore::row");
}

BaseList<complex> DataStore::row()
{
  switch (type()) {
  case 1:
    return (*this)(Int2Type<1>()).row();
  case 2:
    return (*this)(Int2Type<2>()).row();    
  }
  throw Undefined("DataStore::row");
}

const BaseList<complex> DataStore::row(size_t r) const
{
  switch (type()) {
  case 1:
    return (*this)(Int2Type<1>())(r);
  case 2:
    return (*this)(Int2Type<2>()).row(r);    
  }
  throw Undefined("DataStore::row");
}

BaseList<complex> DataStore::row(size_t r)
{
  switch (type()) {
  case 1:
    return (*this)(Int2Type<1>())(r);
  case 2:
    return (*this)(Int2Type<2>()).row(r);    
  }
  throw Undefined("DataStore::row");
}

bool DataStore::empty() const
{
  return (type()==0);//apply(Empty());
}

// size_t DataStore::dimensions() const
// {
//   const size_t n=store_.type();
//   if (n<2)
//     throw Failed("DataStore::dimensions");
//   return n;
// }

size_t DataStore::rows() const 
{
  switch (type()) {
  case 1:
    return (*this)(Int2Type<1>()).size();
  case 2:
    return (*this)(Int2Type<2>()).rows();    
  }
  throw Undefined("DataStore::rows");
  //  return store_.apply(Rows());
}

void DataStore::create(size_t)
{
  set(Int2Type<1>()).clear();  
  //  sws_.create(0); //!< reset spectral 
  procstates_.clear();
}

// double DataStore::sw(size_t r) const
// {
//   if (r>=sws_.size())
//     throw BadIndex("DataStore::sw");
//   return sws_(r);
// }

// double& DataStore::sw(size_t r)
// {
//   if (r>=sws_.size())
//     throw BadIndex("DataStore::sw");
//   return sws_(r);
// }

const processing_state& DataStore::state(size_t r) const
{
  if (r>=procstates_.size())
    throw BadIndex("DataStore::state");
  return procstates_(r);
}

processing_state& DataStore::state(size_t r)
{
  if (r>=procstates_.size())
    throw BadIndex("DataStore::state");
  return procstates_(r);
}

 ListList<complex>& DataStore::listlist()
   {
     return (*this)(Int2Type<1>());
   }

 const ListList<complex>& DataStore::listlist() const
   {
     return (*this)(Int2Type<1>());
   }

cmatrix& DataStore::matrix() 
{
  return (*this)(Int2Type<2>());
}

Warning<> DataStore::mismatchedtype_warning("Combining data sets of different types (nD vs. row-oriented)",&NMRsim_once_warning);

const cmatrix& DataStore::matrix() const
{
  return (*this)(Int2Type<2>());
}

DataStore& DataStore::operator+= (const DataStore& a)
{
  if (type()==0)
    *this=a;
  else
    apply_ip(Add_ip(),a);
  return *this;
}

DataStore& DataStore::operator-= (const DataStore& a)
{
  if (type()==0)
    throw Undefined("DataStore-=");
  apply_ip(Subtract_ip(),a);
  return *this;
}

void DataStore::duplicate_structure_processing(const DataStore& a)
{
  duplicate_structure(a);
  procstates_=a.states();
}

 void DataStore::assign(const cmatrix& v, const BaseList<processing_state>& states)
{
  procstates_=states;
  set(Int2Type<2>())=v;
}

DataStore& DataStore::operator= (const complex& v)
{
  if (row().size()) {
    switch (type()) {
    case 1:
      (*this)(Int2Type<1>())=v;
      break;
    case 2:
      (*this)(Int2Type<2>())=v;
      break;
    default:
      throw Undefined("DataStore=");
    }
  }
  return *this;
}

 void DataStore::createrow(size_t np, const complex& v, const processing_state& state) 
{
  (*this)(Int2Type<1>()).push_back(np,v);
  procstates_.push_back(state);
}

void DataStore::push_back(const BaseList<complex>& v, const processing_state& state)
{
  ListList<complex>& aslistlist(set(Int2Type<1>()));
  aslistlist.push_back(v);
  procstates_.push_back(state);
}

void DataStore::print_processing(std::ostream& ostr) const
{
  for (size_t i=0;i<procstates_.size();i++)
    ostr << "Data row/dimension " << (i+1) << ": " << procstates_(i) << '\n';
}

void DataStore::print_structure(std::ostream& ostr) const
{
  if (isnD()) {
    const cmatrix& mat((*this)(Int2Type<2>()));
    ostr << mat.rows() << " x " << mat.cols();
  }
  else {
    const ListList<complex>& aslistlist((*this)(Int2Type<1>()));
    for (size_t i=0;i<aslistlist.size();i++) {
      if (i)
	ostr << ',';
      ostr << aslistlist.size(i);
    }
  }
  ostr << '\n';
  print_processing(ostr);
}

std::ostream& operator<< (std::ostream& ostr, const processing_state& a)
  {
    ostr << "Spectral width: ";
    if (a.sw)
      ostr << (a.sw*1e-3) << " kHz";
    else
      ostr << "<undefined>";
    if (a.ref) 
      ostr << "  Reference: " << (a.ref*1e-3) << " kHz";
    if (a.sfrq)
      ostr << "  Frequency: " << (a.sfrq*1e-6) << " MHz";
    return ostr << " (" << (a.istimedomain ? "time" : "frequency") << " domain)";    
  }

void DataStore::print(std::ostream& ostr) const
{
  switch (type()) {
  case 0:
    ostr << "<empty>\n";
    return;
  case 1:
    ostr << (*this)(Int2Type<1>()) << '\n';
    return;
  default:
    ostr << (*this)(Int2Type<2>());
  }
  print_processing(ostr);
}

bool arematching(const DataStore& a, const DataStore& b)
{
  if (a.rows()!=b.rows())
    return false;

  for (size_t i=a.rows();i--;) {
    if (!arematching(a.row(i),b.row(i)))
      return false;
  }
  return true;
}

void DataStore::duplicate_structure(const DataStore& a)
{
  //  std::cout << "Store type: " << a.store_.type() << std::endl;
  switch (a.type()) {
  case 0:
    clear();
    break;
  case 1:
    ::libcmatrix::duplicate_structure(set(Int2Type<1>()),a(Int2Type<1>()));
    break;
  case 2:
    ::libcmatrix::duplicate_structure(set(Int2Type<2>()),a(Int2Type<2>()));
    break;
  default:
    throw InternalError("DataStore::duplicate_structure");
  }
}

void ProcessResample::exec(List<complex>& dest, const BaseList<complex>& data, processing_state& pflags) const
{
  if (pflags.istimedomain)
    error_abort("Cannot (or rather won't!) resample time domain data");
  dest.create(newnp_);
  rawresample(dest,data);
  pflags.sw=newsw_; //!< update spectral width
  if (pflags.ref) {
    pflags.ref=0.0;
    refreset_warning.raise();
  }
}

void ProcessResample::exec(cmatrix& data, LIST<processing_state>& pflags) const
{
  if (pflags.back().istimedomain)
    error_abort("cannot (or rather won't!) resample time domain data");
  size_t n=data.rows();
  cmatrix tmp(n,newnp_);
  for (;n--;)
    rawresample(tmp.row(n),data.row(n));
  data.swap(tmp);
  pflags.back().sw=newsw_; //!< update spectral width
  if (pflags.back().ref) {
    pflags.back().ref=0.0;
    refreset_warning.raise();
  }
}
 
void ProcessResample::rawresample(BaseList<complex> dest, const BaseList<complex>& source) const
{
  if (dofold_)
    cubic_interpolate(dest,source);
  else {
    const double olddx=sw/source.size();
    const double oldstart=-(sw+olddx)/2;
    const double newsw=(newsw_<0) ? sw : newsw_;
    const double newdx=newsw/newnp_;
    const double newstart=offset_-(newsw+newdx)/2;
    cubic_interpolate(dest,newstart,newdx,source,oldstart,olddx);
  }
}
 
ProcessResample::ProcessResample(size_t newnpv,bool foldv)
  : ProcessCommand(PROC_HASBOTH | PROC_CHANGES_COLUMNS),
    newnp_(newnpv),
    dofold_(foldv),
    newsw_(-1.0), offset_(0.0)
{}

ProcessResample::ProcessResample(size_t newnpv, double newswv, double offsetv)
  : ProcessCommand(PROC_HASBOTH | PROC_CHANGES_COLUMNS),
    newnp_(newnpv),
    dofold_(false),
    newsw_(newswv), offset_(offsetv)
{}

void ProcessResample::set(double v, subsid_t subsid)
{
  switch (subsid) {
//   case S_ARG1:
//     if (v<0.0) {
//       if (allowwarnings())
// 	std::cerr << "Warning: negative resample spectral encountered (made positive)\n";
//       v=-v;
//     }
//     newsw_=v;
//     break;
  case S_ARG1:
    offset_=v;
    break;
  default:
    throw InternalError("Resample::set");
  }
}   
  
void ProcessResample::print(std::ostream& ostr) const
{
  ostr << "resample " << newnp_;
  if (dofold_)
    ostr << " -fold";
  if (newsw_>=0.0) {
    ostr << ' ' << newsw_;
    if (offset_)
      ostr << ' ' << offset_;
  }
}

void ProcessResample::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  assert(subsid==S_ARG1); //!< harmless assert
  ostr << "resample_offset";
}
  
ProcessCommand* ProcessResample::create()
{
  const size_t newnp=parse_unsigned();
  if (newnp<1)
    error_abort("new number of points cannot be zero");
  double newsw=-1.0;
  double offset=0.0;
  if (are_left() && parser_isnormal()) {
    newsw=parse_double();
    if (newsw<0.0)
      error_abort("spectral width cannot be negative");
    if (are_left()) {
      Variable cvar(S_ARG1);
      offset=parse_double(&cvar,defflags());
    }
  }
  static flagsmap_type flagsmap;
  if (flagsmap.empty())
    flagsmap["fold"]=1;
  const bool dofold=(parse_flags(flagsmap)!=0);
  const bool selrange=(newsw>=0.0);
  if (dofold && selrange)
    error_abort("Cannot combine -fold with frequency range selection");
  return selrange ? new ProcessResample(newnp,newsw,offset) : new ProcessResample(newnp,dofold);
}

bool ProcessNormalise::havealready_=false;
ContextWarning<> ProcessNormalise::multiplenormalise_warning(">1 normalise directives created - this is probably an error",&NMRsim_once_warning);

ProcessNormalise::ProcessNormalise(double valv, metric_t metricv)
    : ProcessCommand(PROC_HASBOTH),
      val_(valv), metric_(metricv) {
  if (havealready_)
    multiplenormalise_warning.raise();
  else
    havealready_=true;
}

const flagsmap_type& ProcessNormalise::get_flags()
{
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["integral"]=INTEGRAL;
    flags["minmax"]=MINMAX;
    flags["abs"]=ABS;
    flags["area"]=AREA;
  }
  return flags;
}

ProcessCommand* ProcessNormalise::create()
{
  double val=1.0;
  if (parser_isnormal()) {
    val=parse_double();
    if (val==0.0)
      error_abort("can't (meaningfully) normalise to zero");
  }
  const metric_t metric=(metric_t)parse_flags(get_flags(),size_t(MINMAX),true);
  return new ProcessNormalise(val,metric);
}

void ProcessNormalise::print(std::ostream& ostr) const
{
  ostr << "normalise " << val_ << ' ';
  switch (metric_) {
  case INTEGRAL: ostr << "-integral"; break;
  case MINMAX: ostr << "-minmax"; break;
  case ABS: ostr << "-abs"; break;
  case AREA: ostr << "-area"; break;
  }
}
  
void ProcessNormalise::exec(BaseList<complex> FID, processing_state& state) const
{
  exec(FID,state.sw);
}

void ProcessNormalise::exec(cmatrix& FID, LIST<processing_state>& state) const
{
  double prodsw=1.0;
  if (metric_==AREA) {
    for (size_t i=state.size();i--;)
      prodsw*=state(i).sw;
  }
  exec(FID.row(),prodsw);
}

void ProcessNormalise::exec(BaseList<complex> FID, double sw) const
{
  double v=0.0;
  switch (metric_) {
  case INTEGRAL: case AREA:
    for (size_t i=FID.size();i--;)
      v+=real(FID(i));
    if (metric_==AREA)
      v*=sw;
    break;

  case MINMAX:
    v=fabs(real(FID.front()));
    for (size_t i=FID.size()-1;i>0;i--) {
      double lv=fabs(real(FID(i)));
      if (lv>v)
	v=lv;
    }
    break;

  case ABS:
    v=norm(FID.front());
    for (size_t i=FID.size()-1;i>0;i--) {
      double lv=norm(FID(i));
      if (lv>v)
	v=lv;
    }
    v=sqrt(v); //!< can delay sqrt to here
    break;
  }
  if (v==0.0)
    error_abort("Normalisation failed - metric was zero");
  FID*=(val_/v);
}

void ProcessApply::set(const BaseList<double>& vs, subsid_t n)
{
  size_t ni=(size_t)n;
  if (ni<1 || (ni>argvals_.size()))
    throw InternalError("ProcessApply::set");
  argvals_(ni-1)=vs;
}

void ProcessApply::print(std::ostream& ostr) const
{
  ostr << "apply " << fname_;
  if (ismultiargument()) {
    ostr << ' ' << argn_;
    for (size_t i=0;i<argvals_.size();i++) {
      if (i!=argn_-1) {
	const BaseList<double>& clist(argvals_(i));
	if (clist.size()!=1)	
	  ostr << ' ' << clist;
	else
	  ostr << ' ' << clist.front();
      }
    }
  }
}

void ProcessApply::printvariablename(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "apply_" << fname_;
  if (ismultiargument())
    ostr << '_' << (size_t)subsid;
}

ProcessApply::ProcessApply(function_spec_t fspecv, complexity_t complexityv, size_t argnv, const BaseList< LIST<double> >& argvalsv, bool isconst)
  : ProcessCommand(PROC_HAS1D | (((complexityv==COMPLEX) || (complexityv==COMPLEXPAIR)) ? PROC_CHANGES_COLUMNS : 0) | (isconst ? 0 : PROC_VARIABLE)),
    fname_(fspecv.name()),
    complexity_(complexityv),
    argn_(argnv),
    argvals_(argvalsv.size()+1),
    expr_(fspecv,argvals_)
{
  for (size_t i=argvalsv.size();i--;) {
    const size_t dest=(i<argnv-1) ? i : i+1;
    argvals_(dest)=argvalsv(i);    
  }
}

ContextWarning<> ProcessApply::singleargument_warning("argument number is redundant when data set is only argument",&NMRsim_once_warning);

ProcessCommand* ProcessApply::create()
{  
  const char* fname=parse_string();
  LIST< LIST<double> > args;
  size_t nargs=0;
  Variable cvar;
  Mark markobj;
  int argn=1;
  bool isconst=true;
  const int flags=defflags();
  if (parser_isnormal()) {
    argn=parse_int();
    if (argn<1)
      error_abort("argument number must be >0");    
    if (!parser_isnormal())
      singleargument_warning.raise();
    while (parser_isnormal()) {
      cvar.subsid=nargs+1;
      VariableBase* valp=parse_double_variable(cvar,F_ALLOWLIST | flags);
      args.push_back(valp->value());
      if (!(valp->isconst()))
	isconst=false;
      nargs++;
    }
  }
  nargs++; //!< count data argument
  if (argn>nargs) {
    parser_printcontext() << "argument number (" << argn << ") is greater than number of supplied arguments (" << nargs << ").\nSyntax for function with more than one argument is apply <function> <arg no. for data set> <remaining arguments>.\n";
    error_abort();
  }
  Function_Factory_t& factory(get_Function_Factory());
  const Function_Factory_t::iterator iter(factory.find(function_spec(fname,nargs)));
  if (iter==factory.end()) {
    parser_printcontext() << "no function named " << fname << " with " << nargs << " argument(s)\n";    
    error_abort();
  }
  const complexity_t complexity=parsecomplexity(COMPLEX,true);  
  return new ProcessApply(*(iter->second), complexity, (size_t)argn, args, isconst);
}

void ProcessApply::zip(LIST<complex>& dest, const BaseList<double>& r, const BaseList<double>& c) const
{
  size_t n=r.size();
  if (n!=c.size()) {
    parser_printthread(std::cerr) << "apply " << fname_ << ": real and complex components have incompatible sizes, " << r.size() << " and " << c.size() << " respectively\n";
    error_abort();
  }
  dest.create(n);
  for (size_t i=n;i--;)
    dest(i)=complex(r(i),c(i));
}

void ProcessApply::apply(LIST<double>& dest, const BaseList<double>& source) const
{
  argvals_(argn_-1)=source; //!< copy data argument 
  expr_.get(dest); //!< evaluate function
}
		
void ProcessApply::exec(LIST<complex>& dest, const BaseList<complex>& source, processing_state&) const
{
  //  (void)checkcomplexity(source.size(),complexity_);
  switch (complexity_) {
  case COMPLEXPAIR:
    apply(bufr_,asdoubles(source));
    (void)checkcomplexity(bufr_.size(),complexity_);
    dest.create(bufr_.size()/2);
    asdoubles(dest)=bufr_;
    break;
  case COMPLEX:
    tmp_=real(source);
    apply(bufr_,tmp_);
    tmp_=imag(source);
    apply(bufc_,tmp_);
    zip(dest,bufr_,bufc_);
    break;
  default:
    throw InternalError("ProcessApply::exec");
  }
}

void ProcessApply::exec(BaseList<complex> FID, processing_state&) const
{
  switch (complexity_) {
  case REAL:
    tmp_=reals(FID);
    apply(bufr_,tmp_);
    reals(FID)=bufr_;
    break;
  case IMAG:
    tmp_=imags(FID);
    apply(bufc_,tmp_);
    imags(FID)=bufc_;
    break;
  default:
    throw InternalError("ProcessApply::exec2");
  }
}

ThreadWarning<> notsimpleelement_warning("MATLAB: dt/sw/type/ref specification is not 1 or 2 element array - using first elements",&NMRsim_once_warning);
ThreadWarning<> repeatdomain_warning("MATLAB: dt/sw/type specified more than once - using last specification",&NMRsim_once_warning);
ThreadWarning<> zerodt_warning("MATLAB: dt/sw is zero - value ignored",&NMRsim_once_warning);
ContextWarning<> containernotstruct_warning("Container for requested data item is not a data structure: ",&NMRsim_repeat_warning);

namespace {
  inline bool matchlast(const char* a, size_t n, const char* b)
  {
    const int which=strlen(a)-n;
    return ((which>=0) && (strcmp(a+which,b)==0));
  }

  enum found_t { NONE, DATA, DT, SW, TYPE, REF, SFRQ, DATAARRAY};
  
  void set(domain_t& dom, double& sw, found_t found, double v)
  { 
    if (found==TYPE) {
      dom= (v==FD_TYPE_FID) ? D_TIME : D_FREQUENCY;
      return;
    }
    if (v==0) 
      zerodt_warning.raise();
    else {
      if (found==DT) {
	dom=D_TIME;
	sw=1.0/v;
      }
      else
	sw=v;
    }
  }

  void readinfo_(double& v, double& v1, const BaseList<double>& datarow, double scalef =1.0)
  {
    v=scalef*datarow.front();
    if (datarow.size()>1)
      v1=scalef*datarow(size_t(1));
  }

  struct varname_search {
    varname_search(const LIST<std::string>* =NMRSIM_NULL);
    LIST<const char*> basenames;
    LIST<const char*> leafnames;
    LIST<const char*> varnames;
    size_t nvars;

    bool checkmatch(matlab_controller::composite&, const char*, const char*);
    void checkleft(std::ostream&) const;
  };

  varname_search::varname_search(const LIST<std::string>* varnamesp)
  {
    if (!varnamesp) {
      nvars=0;
      return;
    }
    nvars=varnamesp->size();
    basenames.create(nvars,(const char*)NULL);
    leafnames.create(nvars,(const char*)NULL);
    varnames.create(nvars,(const char*)NULL);
    for (size_t i=nvars;i--;) {
      char* rawstr=const_cast<char*>((*varnamesp)(i).c_str()); //!< a bit naughty, but we can manipulate in place
      char* colonpos=strchr(rawstr,':');
      if (colonpos) {
	*colonpos++='\0';
	varnames(i)=colonpos;
      }
      char* dotpos=strchr(rawstr,'.');
      if ((*rawstr=='\0') || !dotpos || (dotpos[1]=='\0') || (strchr(dotpos+1,'.')!=NULL))
	error_abort("variables names to be read from data sets need to be of form <base>.<name> e.g. arrays.ct");
      *dotpos++='\0';
      basenames(i)=rawstr;
      leafnames(i)=dotpos;
    }
  }

  bool varname_search::checkmatch(matlab_controller::composite& ctrl, const char* first, const char* second)
  {
    static rmatrix datatmp;
    for (size_t i=basenames.size();i--;) {
      if ( (strcmp(first,basenames(i))==0) && (strcmp(second,leafnames(i))==0)) {
	const char* varname=varnames(i);
	char scratch[NMRSIM_MATLAB_VARMAX];
	if (!varname) {
	  snprintf(scratch,NMRSIM_MATLAB_VARMAX,"%s_%s",basenames(i),leafnames(i));
	  varname=scratch;
	}
	UserVariable& var=*findnewvariable(strdup(varname)); //!< leak but don't care
	ctrl.read(datatmp);
	const BaseList<double> datarow(datatmp.row());
	if ((datarow.size()!=1) && (datarow.size()!=2))
	  notsimpleelement_warning.raise();
	var.set(datarow);
	var.isconstant(true);
	leafnames(i)=NULL; //!< flag read successfully
	return true;
      }
    }
    return false;
  }

  void varname_search::checkleft(std::ostream& ostr) const
  {
    parser_printthread(ostr) << "MATLAB read: failed to read all requested data variables (use -verbose:parse to view available arrays):";	  
    for (size_t i=0;i<nvars;i++) {
      if (leafnames(i)!=NULL)
	ostr << ' ' << basenames(i) << '.' << leafnames(i);
    }
    ostr << std::endl;
  }
	
  template<class T> bool read_matlab_data(T& ctrl, cmatrix& a, filestruct& finfo, varname_search& searchobj)
  {
    try {
      size_t varsleft=searchobj.nvars;

      matlab_controller::header_info info;
      finfo.domain=D_UNKNOWN;
      rmatrix datatmp;
      bool gotdata=false;
      while (ctrl.peek(info)) {
	found_t found=NONE;
	if (strcmp(info.name,"data")==0)
	  found=DATA;
	else if (strcmp(info.name,"dt")==0)
	  found=DT;
	else if (strcmp(info.name,"sw")==0)
	  found=SW;
	else if (strcmp(info.name,"type")==0)
	  found=TYPE;
	else if (strcmp(info.name,"ref")==0)
	  found=REF;
	else if (strcmp(info.name,"sfrq")==0)
	  found=SFRQ;
	else {
	  if ( (info.type==matlab_controller::STRUCT) && (varsleft!=0)) {
	    found=DATAARRAY;
	    matlab_controller::composite array(ctrl);
	    matlab_controller::header_info arrayinfo;
	    const bool dumpinfo=(verbose & VER_PARSE);
	    if (dumpinfo)
	      parser_printthread() << "Available data items:";
	    while (array.peek(arrayinfo)) {
	      if (dumpinfo)
		std::cout << ' ' << info.name << '.' << arrayinfo.name;
	      if (searchobj.checkmatch(array,info.name,arrayinfo.name)) 
		varsleft--;
	      array.next();
	    }
	    if (dumpinfo)
	      std::cout << std::endl;
	  }
	}
	
	if (found!=NONE) {
	  if (found==DATAARRAY)
	    continue; //!< already 'processed'

	  if (info.type!=matlab_controller::ARRAY) {
	    parser_printthread(std::cerr) << "MATLAB read: '" << info.name << "' array is not a simple matrix\n";
	    error_abort();
	  }	  
	  if (found==DATA) {
	    ctrl.read(a);
	    gotdata=true;
	  }
	  else {
	    ctrl.read(datatmp);
	    const BaseList<double> datarow(datatmp.row());
	    if ((datarow.size()!=1) && (datarow.size()!=2))
	      notsimpleelement_warning.raise();
	    switch (found) {
	    case REF:
	      readinfo_(finfo.ref,finfo.ref1,datarow);
	      break;
	    case SFRQ:
	      readinfo_(finfo.sfrq,finfo.sfrq1,datarow,1e6); //!< convert to Hz (23/7/15)
	      break;
	    default:
	      if ((finfo.domain!=D_UNKNOWN) && (found!=TYPE))
		repeatdomain_warning.raise();

	      set(finfo.domain,finfo.sw,found,datarow.front());	      
	      if (datarow.size()>1)
		set(finfo.domain1,finfo.sw1,found,datarow(size_t(1)));
	    }
	  }
	}
	else
	  ctrl.next();
      }
      if (gotdata) {
	if (finfo.domain==D_UNKNOWN)
	  finfo.domain=D_FREQUENCY;

	if (varsleft) {
	  searchobj.checkleft(std::cerr);
	  error_abort();
	}

	return true;
      }
    }
    catch (const std::exception& exc) { //!< problem here indicates file is corrupt
      error_abort(exc.what());
    }
    finfo.clear(); // reset (file pointer will be reset by lock)
    return false;
  }
}

filestruct& filestruct::operator= (const simpsonFD& filedesc)
{    
  sw=filedesc.sw;
  sw1=filedesc.sw1;
  domain=filedesc.istimedomain() ? D_TIME : D_FREQUENCY; //!< NB SIMPSON file format doesn't explicitly specify domain if indirect dimension, so do nothing
  sfrq=filedesc.sfrq;
  sfrq1=filedesc.sfrq1;
  ref=filedesc.ref;
  ref1=filedesc.ref1;
  return *this;
}

bool isdirectory(const char* s)
{
  struct stat buf;
  if (stat(s,&buf))
    return false;
  return S_ISDIR(buf.st_mode);
}

typedef std::map< std::string, std::string> parbuf_t;

void read_spinsight_parfile(parbuf_t& pars, const char* fname)
{
  char linebuf[MAXLINE];

  FILE* fp=file_open(fname,"ra");

  while (fgets(linebuf,MAXLINE,fp)) {
    char* equalptr=strchr(linebuf,'=');
    if (!equalptr)
      continue;
    *equalptr++='\0';
    chewwhite(equalptr);
    pars[linebuf]=equalptr;
  }
  fclose(fp);
}

/*function for reading SPINSIGHT 1D FIDs and spectra. FID stores in "data" file in binary format.
4 bytes per each (real or imag) point in big-endian format (as integer). Spectra - the same but as float. First half of file stores the real part, second half -imaginary
adopted partly from MatNMR of J. vanBeek via gsim. */

inline bool isspaceorend(char c) { return (c=='\0') || isspace(c); }

long tolong(const char* s)
{
  char* end;
  const long n=strtol(s,&end,10);
  if ((end==s) || !isspaceorend(*end)) {
    char buf[MAXLINE];
    snprintf(buf,sizeof(buf)-20,"Failed to parse %s as integer",s);
    throw Failed(buf);
  }
  return n;
}

size_t tochannel(const char* s)
{
  const long n=tolong(s);
  if (n<1 || n>4) {
    char buf[MAXLINE];
    snprintf(buf,sizeof(buf)-20,"%s is not a valid channel (1-4)",s);
    throw Failed(buf);
  }
  return n;
}

double todouble(const char* s)
{
  char* end;
  const double v=strtod(s,&end);
  if (end==s) {
    char buf[MAXLINE];
    snprintf(buf,sizeof(buf)-30,"Failed to parse %s as floating point number",s);
    throw Failed(buf);
  }
  return v;
}

domain_t todomain(const char* s, const char* name)
{
  const long type=tolong(s);
  switch (type) {
  case 1:
    return D_FREQUENCY;
  case 0:
    return D_TIME;
  }
  char buf[MAXLINE];
  snprintf(buf,sizeof(buf)-30,"read_spinsight: %s in 'proc' file is neither 0 nor 1",name);
  throw Failed(buf);
}

double process_ref(double rmp, double rmv, int units, double sfrq)
{
  if (sfrq==0.0)
    throw Failed("process_ref: spectrometer frequency undefined");
  const double base=sfrq-rmp*1e6;

  switch (units) { 
  case 1://ppm
    rmv*=sfrq*1e-6; //convert to Hz
    return rmv+base;

  case 2://kHz
    return rmv*1e3+base;

  case 5://Hz
    return rmv+base;
  }
  throw Failed("read_spinsight: unhandled reference type");
}

  class parbuf_parser {
  public:
    parbuf_parser(const parbuf_t& parbufv)
      : parbuf(parbufv), end(parbufv.end()) {}

    void tolong(size_t& dest, const char* key) const {
      const std::string skey(key);
      const parbuf_t::const_iterator iter(parbuf.find(skey));
      if (iter!=end)
	dest=::tolong( (iter->second).c_str() );
    }

    void todouble(double& dest, const char* key, double scale =1.0) const {
      const parbuf_t::const_iterator iter(parbuf.find(std::string(key)));
      if (iter!=end)
	dest=scale*::todouble( (iter->second).c_str() );
    }

    void tochannel(size_t& dest, const char* key) const {
      const parbuf_t::const_iterator iter(parbuf.find(std::string(key)));
      if (iter!=end)
	dest=::tochannel( (iter->second).c_str() );
    }

    void todomain(domain_t& dest, const char* key) const {
      const parbuf_t::const_iterator iter(parbuf.find(std::string(key)));
      if (iter!=end)
	dest=::todomain( (iter->second).c_str(), key );
    }

  private:
    const parbuf_t& parbuf;
    const parbuf_t::const_iterator end;
  };

void read_spinsight(cmatrix& dest, filestruct& finfo, const char *fname)
{
  char buf[MAXPATH];

  if (strlen(fname)+10>MAXPATH)
    throw Failed("read_spinsight: path name overflow");

  snprintf(buf,sizeof(buf),"%s/data",fname);
  if (!isreadable(buf))
    throw Failed("read_spinsight: failed to open <dir>/data file");

  parbuf_t pars;

  snprintf(buf,sizeof(buf),"%s/acq",fname);
  if (!isreadable(buf))
    throw Failed("read_spinsight: failed to open <dir>/acq file");
  read_spinsight_parfile(pars,buf);

  bool data_istimedomain=true;
  double freq[4]={0.0, 0.0, 0.0, 0.0};
  size_t ch1=0, ch2=0;

  const parbuf_parser acqparser(pars);
  acqparser.todouble(freq[0],"sf1",1e6);
  acqparser.todouble(freq[1],"sf2",1e6);
  acqparser.todouble(freq[2],"sf3",1e6);
  acqparser.todouble(freq[3],"sf4",1e6);
  acqparser.tochannel(ch1,"ch1");
  acqparser.tochannel(ch2,"ch2");
  acqparser.todouble(finfo.sw,"sw",1e3);
  acqparser.todouble(finfo.sw1,"sw2",1e3);
  
  //Determine frequency in both dimensions
  if (ch1)
    finfo.sfrq=freq[ch1-1];
  if (ch2)
    finfo.sfrq1=freq[ch2-1];
  if (!finfo.sfrq1)
    finfo.sfrq1=finfo.sfrq; //probably homonuclear experiment

  size_t ni=1,np;

  //Read proc file
  snprintf(buf,sizeof(buf),"%s/proc",fname);
  pars.clear();
  read_spinsight_parfile(pars,buf);

  double rmp1=0.5, rmp2=0.5;
  double rmv1=0.0, rmv2=0.0;
  size_t rmvunits1=0, rmvunits2=0; //!< flag unset

  parbuf_parser procparser(pars);
  parbuf_t::const_iterator iter(pars.find("datatype"));
  if (iter!=pars.end())
    data_istimedomain=(todomain( (iter->second).c_str(),"datatype")==D_TIME);
  procparser.todomain(finfo.domain,"domain1");
  procparser.todomain(finfo.domain1,"domain2");
  procparser.tolong(np,"current_size1");
  procparser.tolong(ni,"current_size2");
  procparser.todouble(rmp1,"rmp1");
  procparser.todouble(rmp2,"rmp2");
  procparser.todouble(rmv1,"rmv1");
  procparser.todouble(rmv2,"rmv2");
  procparser.tolong(rmvunits1,"rmvunits1");
  procparser.tolong(rmvunits2,"rmvunits2");

  if (rmvunits1>0)
    finfo.ref=process_ref(rmp1,rmv1,rmvunits1,finfo.sfrq);
  if (rmvunits2>0)
    finfo.ref1=process_ref(rmp2,rmv2,rmvunits2,finfo.sfrq1);

  //if (lPar.contains("com"))
  //	par.title=lPar["com"].remove('\n');

  //open data file
  snprintf(buf,sizeof(buf),"%s/data",fname);
  FILE* pFile=file_open(buf, "rb");
  const off_t lSize = file_length(pFile);
  if (lSize & 7) {
    fclose(pFile);
    throw Failed("read_spinsight: file length is not multiple of 8. Corrupt?");
  }

  const size_t len=lSize/8;
  if ((np*ni)!=len) {
    fclose(pFile);
    throw Mismatch("read_spinsight: file length does not match data set size parameters. Corrupt?",len,np*ni);
  }

  List<int32_t> datablock(2*len);    
  int32_t* datat = datablock.vector();
  // copy the file into the buffer.
  fread( datat, 4, 2*len, pFile);
  fclose(pFile);

  // Change byte order if machine is little-endian
  if (!ambigendian()) {
    for (size_t i=2*len;i--;)
      endianswap(datat[i]);
  }
  
  dest.create(ni, np);
  BaseList<complex> asrow(dest.row());
  if (data_istimedomain) {
    for (size_t i=len;i--;)
      asrow(i)=complex(datat[i],datat[i+len]);
  }
  else {
    if (sizeof(float)!=4)
      throw Failed("read_spinsight: 'float' data type on this machine is incompatible");
    const float* asfloat=reinterpret_cast<float*>(datat);
    for (int i=len;i--;)
      asrow(i)=complex(asfloat[i],asfloat[i+len]);
  }
}

ContextWarning<> ignoredvarnames_warning("request to read supporting data arrays ignored (only currently supported for data in MATLAB files)",&NMRsim_repeat_warning);

// Try read as MATLAB (if finishes with .mat) otherwise try reading as SIMPSON file, then as plain ASCII
void read_file_(cmatrix& tmp_set, filestruct& finfo, const char* fname, LIST<std::string>* varnamesp)
  {
    if (isdirectory(fname)) {
      try {
	read_spinsight(tmp_set, finfo, fname);
	if (varnamesp)
	  ignoredvarnames_warning.raise();
	return;
      }
      catch (MatrixException& exc) {
	if (verbose & VER_GEN)
	  parser_printthread(std::cerr) << "Attempt to read " << fname << " as SpinSight data directory failed: " << exc << '\n';
      }
      parser_printcontext() << "Couldn't open data directory: " << fname << " (use verbose -general for more information)\n";
      error_abort();
    }
    if (matchlast(fname,4,".mat")) {
      try {
	varname_search searchobj(varnamesp);

	char fnametmp[MAXPATH];
	strcpy(fnametmp,fname);
	*(strrchr(fnametmp,'.'))='\0';
	matlab_controller ctrl(fnametmp);
	ctrl.padmode=matlab_controller::PADROW;
	matlab_controller::padding_warning.type(BaseWarning::FirstOnly);
	{
	  file_position_lock lock(ctrl.file_pointer());
	  if (read_matlab_data(ctrl,tmp_set,finfo,searchobj)) {
	    lock.unlock();
	    return;
	  }
	}
	matlab_controller::header_info info;
	for (;;) {
	  ctrl.peek(info);
	  switch (info.type) {
	  case matlab_controller::CHAR:
	    break;

	  case matlab_controller::ARRAY:
	    ctrl.read(tmp_set);
	    return;

	  case matlab_controller::CELL: 
	    error_abort("MATLAB file consists of cell array - input must be individual matrix");

	  case matlab_controller::STRUCT: {
	    matlab_controller::composite comp(ctrl);
	    if (read_matlab_data(comp,tmp_set,finfo,searchobj))
	      return;
	    error_abort("MATLAB structure does not contain array called data");	    
	  }
	  default:
	    error_abort("Unknown MATLAB data structure");
	  }
	  if (!(ctrl.next()))
	    error_abort("MATLAB file did not contain data\n");	    
	}
      } catch (...) {}
    }
    else {      
      try {
	simpsonFD filedesc;
	read_simpson(tmp_set,filedesc,fname);
	finfo=filedesc;
	if (varnamesp)
	  ignoredvarnames_warning.raise();
	return;
      }
      catch (MatrixException& exc) {
	if (verbose & VER_GEN)
	  parser_printthread(std::cerr) << "Attempt to read " << fname << " as SIMPSON failed: " << exc << '\n';
      }

      try {
	read_matrix(tmp_set,fname);
	return;
      }
      catch (MatrixException& exc) {
	if (verbose & VER_GEN)
	  parser_printthread(std::cerr) << "Attempt to read " << fname << " as simple matrix failed: " << exc << '\n';
      }
    }
    parser_printcontext() << "Couldn't open data set: " << fname << " (use verbose -general for more information)\n";
    error_abort();
  }

void raw_read_file(cmatrix& tmp_set, filestruct& info, const char* fname, LIST<std::string>* varnamesp)
{
  read_file_(tmp_set,info,fname,varnamesp);
  if (tmp_set.cols()==1)
    tmp_set.transpose();
}

 void raw_build_set(DataStore& fit_set, const BaseList<cmatrix>& fit_sets, const BaseList<processing_state>& fit_procs, int verbose)
 {
  //assemble DataStore to fit
   if (active2D) { //!< restored this 13/2/14
   //   if (!isirregular()) {
    if (fit_sets.size()!=1)
      error_abort("Can't fit nD spectrum to multiple data sets");    
    if (fit_procs.size()==ndims)
      fit_set.assign(fit_sets.front(),fit_procs); //!< if info for all domains provided, use this
    else {
      const bool deftdom=fit_procs.front().istimedomain; //!< default domain matches direct dimension
      ScratchList<processing_state> lstates(ndims,processing_state(0.0,deftdom,0.0)); //!< null processing of other dimensions
      lstates.back()=fit_procs.front();
      fit_set.assign(fit_sets.front(),lstates);
      //    if (fit_procs.front())
      //    fit_set.states().back().sw=fit_sws.front(); //!< copy out sw if supplied
    }
  }
  else {
    if (fit_sets.size()!=array_n0)
      throw Mismatch("Number of supplied data sets and number of rows in calculated data set don't match",fit_sets.size(),array_n0);
    fit_set.create(fit_sets.size());
    //const bool istime=!isfrequency();
    for (size_t i=0;i<fit_sets.size();i++) {
      const cmatrix& source(fit_sets(i));
      if (source.rows()!=1)
	throw Failed("Can't fit irregular data set to matrix");

//       if (allowwarnings() && (fit_procs(i).istimedomain!=istime))
// 	std::cerr << "Warning: time vs. frequency domain status of row " << (i+1) << " may not match simulation\n";
      fit_set.push_back(source.row(),fit_procs(i));
    }
  }
  if (verbose) {
    if (verbose>1)
      std::cout << "Data set\n" << fit_set;
    else {
      std::cout << "Data set stucture: ";
      fit_set.print_structure(std::cout);
    }
  }  			
 }

void ProcessSetDomain::print(std::ostream& ostr) const
{
  ostr << "setdomain ";
  if (dim_)
    ostr << dim_ << ' ';
  dumpflags(ostr,getflags(),flags_);
}

const flagsmap_type& ProcessSetDomain::getflags()
{
  static flagsmap_type flags;
  if (flags.empty()) {
    flags["switch"]=SWITCH;
    flags["time"]=TIME;
    flags["frequency"]=FREQUENCY;
    flags["States"]=STATES;
    flags["noStates"]=NOSTATES;
  }
  return flags;
}

void ProcessSetDomain::exec_(processing_state& state) const
{
  switch (flags_ & dommask) {
  case SWITCH:
    state.istimedomain=!state.istimedomain; 
    break;
  case FREQUENCY:
    state.istimedomain=false;
    break;
  case TIME:
    state.istimedomain=true;
    break;
  case 0:
    break;
  default:
    throw InternalError("ProcessSetDomain::exec_");
  } 
}

void ProcessSetDomain::exec(BaseList<complex>, processing_state& state) const
{
  if (flags_ & (STATES | NOSTATES))
    throw InternalError("ProcessSetDomain::exec (1D)");
  if (flags_==0)
    std::cout << "Domain: " << (state.istimedomain ? "Time" : "Frequency");
  else
    exec_(state);
}

ContextWarning<> ProcessSetDomain::fd_warning("setdomain: changing States flag for frequency domain data",&NMRsim_once_warning);

void ProcessSetDomain::exec(cmatrix&, LIST<processing_state>& states) const
{
  if (flags_==0) {
    const size_t startdim=dim_ ? dim_-1 : 0;
    const size_t enddim=dim_ ? dim_-1 : ndims-1;
    for (size_t k=startdim;k<=enddim;k++) {
      const processing_state& state(states(k));
      std::cout << "Dimension " << (k+1) << ": " << (state.istimedomain ? "Time" : "Frequency");
      if (state.istimedomain && (k!=ndims-1))
	std::cout << "  States: " << ((skips(k)==2) ? "Yes" : "No");
      std::cout << '\n';
    }
    return;
  }
  processing_state& state(states(dim_-1));
  exec_(state);  
  if (flags_ & statesmask) {
    if (!(state.istimedomain))
      fd_warning.raise();
    if (dim_==ndims)
      throw InternalError("ProcessSetDomain::exec (nD)");
    const size_t n=ns(dim_-1);
    size_t skip=1;
    if (flags_ & STATES) {
      if (n & 2)
	error_abort("setdomain -States; number of rows is not even");
      skip=2;
    }
    process_set_n(dim_,n,skip);
  }
}

ProcessCommand* ProcessSetDomain::create()
{
  int dim=0;
  if (parser_isnormal()) {
    dim=parse_int();
    if ((dim<1) || (dim>ndims)) {
      parser_printcontext() << " dimension must be between 1 and " << ndims << '\n';
      error_abort();
    }
  }
  const int flags=parse_flags(getflags());
  if (flags && (dim==0))
    dim=ndims;
  if ((flags & statesmask)==statesmask)
    error_abort("Can't set both States and noStates");
  switch (flags & dommask) {
  case 0: case TIME: case SWITCH:
    break;
  case FREQUENCY:
    if (flags & statesmask)
      fd_warning.raise();
    break;
  default:
    error_abort("can't combine -switch, -time, -frequency");
  }
  return new ProcessSetDomain(dim,flags); 
}    

const int ProcessSet::indmask_=ProcessSet::REF1 | ProcessSet::SW1 | ProcessSet::SFRQ1;
flagsmap_type ProcessSet::flags_;

void ProcessSet::exec(BaseList<complex>, processing_state& cstate) const
{
  if (which_ & indmask_)
    notnd_warning.raise();
  rawset(cstate);
}

void ProcessSet::exec(cmatrix&, LIST<processing_state>& cstates) const
{
  rawset(cstates.back(),0U);
  if (cstates.size()>1)
    rawset(cstates(cstates.size()-2),1U);
  else {
    if (which_ & indmask_)
      notnd_warning.raise();
  }
}
 
ThreadWarning<> ProcessSet::notnd_warning("Cannot set indirect dimension parameters for non-nD data set",&NMRsim_once_warning);

void ProcessSet::rawset(processing_state& state, size_t dim) const
{ 
  switch (dim) {
  case 0:
    if (which_ & SW)
      state.sw=sw_;
    if (which_ & REF)
      state.ref=ref_;
    if (which_ & SFRQ)
      state.sfrq=sfrq_;
    break;
  case 1:
    if (which_ & SW1)
      state.sw=sw1_;
    if (which_ & REF1)
      state.ref=ref1_;
    if (which_ & SFRQ1)
      state.sfrq=sfrq1_;
    break;
  default:
    throw InternalError("ProcessSet does not support >2 dimensions");
  }
}

void ProcessSet::print(std::ostream& ostr) const
{
  ostr << "set";
  flagsmap_type::const_iterator start=flags_.begin();
  const flagsmap_type::const_iterator end=flags_.end();

  while (start!=end) {
    const size_t val=start->second;
    if (which_ & val) {
      ostr << ' ';
      switch (val) {
      case SW:
	ostr << sw_;
	break;
      case SW1:
	ostr << sw1_;
	break;
      case REF:
	ostr << ref_;
	break;
      case REF1:
	ostr << ref1_;
	break;
       case SFRQ:
 	ostr << sfrq_;
 	break;
       case SFRQ1:
 	ostr << sfrq1_;
 	break;
      default:
	throw InternalError("ProcessSet::print");
      }
      ostr << " -" << (start->first);
    }
    ++start;
  }
}

void ProcessSet::set(double v, subsid_t which)
{
  switch (which) {
  case SW:
    sw_=v;
    break;
  case SW1:
    sw1_=v;
    break;
  case REF:
    ref_=v;
    break;
  case REF1:
    ref1_=v;
    break;
  case SFRQ:
    sfrq_=v;
    break;
  case SFRQ1:
    sfrq1_=v;
    break;
  default:
    throw InternalError("ProcessSet::set");
  }
}

ContextWarning<> ProcessSet::repeatset_warning("overriding previous set for ",&NMRsim_once_warning);

ProcessCommand* ProcessSet::create()
{
  static const char syntaxstring[]="invalid number of arguments to set [-sw|-sw1|-ref|-ref1|-sfrq|-sfrq1 <value>] (repeat)";
  static int cumulativewhich=0;

  if (flags_.empty()) {
    flags_["sw"]=SW;
    flags_["sw1"]=SW1;
    flags_["ref"]=REF;
    flags_["ref1"]=REF1;
    flags_["sfrq"]=SFRQ;
    flags_["sfrq1"]=SFRQ1;
    flags_["type"]=0;
  }

  ProcessSet* setcomp=new ProcessSet();
  int which=0;
  while (are_left()) {
    const char* lab=parse_string();
    if (*lab!='-') {
      parser_printcontext() << "Failed to parse as flag: " << lab << '\n';
      which=0;
      break;
    }
    lab++;
    const flagsmap_type::const_iterator flag=flags_.find(lab);
    if (flag==flags_.end()) {
      parser_printcontext() << "unknown set type: " << lab << '\n';      
      error_abort();
    }
    const size_t lwhich=flag->second;
    if (lwhich==0)
      error_abort("use setdomain rather than set -type fid|spe");
    if (cumulativewhich & lwhich)
      repeatset_warning.raise(lab);
    Variable cvar(lwhich);
    const double val=parse_double(&cvar);
    setcomp->set(val,lwhich);    
    which|=lwhich;
  }
  if (which==0)
    error_abort(syntaxstring);

  setcomp->which_=which;
  cumulativewhich|=which;
  
  return setcomp;
}

void write_vector(FILE* fp, const BaseList<bool>& a, int flags)
{   
  char termchar=(flags & mxflag::block) ? ' ' : '\n';

  for (size_t i=0;i<a.size();i++)
    fprintf(fp,"%c%c",(a(i) ? '1' : '0'),termchar);
}

ThreadWarning<> writematrixflag_warning("Flag/format not supported for writing Matrix<bool>",&NMRsim_once_warning);
 
void write_matrix(FILE* fp, const Matrix<bool>& a, const char* comment, int flags)
{
  if (flags & (mxflag::rmat | mxflag::binary | mxflag::doublep))
    writematrixflag_warning.raise();

  write_comment(fp,comment,"%");
  
  if ((flags & mxflag::norowsep) & !(flags & mxflag::block)) {
    write_vector(fp,a.row(),flags);
    return;
  }

  const size_t rows=a.rows();
  for (size_t i=0;i<rows;i++) {
    write_vector(fp,a.row(i),flags);
    if (i<rows-1)
      putc('\n',fp);
  }
  if (flags & mxflag::block)
    putc('\n',fp);
}

ProcessCommand* ProcessSystem::create()
{
  static char comline[MAXLINE];
  const char* restcom=get_curline();
  if (restcom && *restcom) {
    substitute_string(comline,sizeof(comline),restcom,SUB_NUMERIC);
    ProcessCommand* newcom=new ProcessSystem(comline);
    set_curline(NMRSIM_NULL);
    return newcom;
  }
  else
    throw Failed("missing command line");
}

void ProcessSystem::print(std::ostream& ostr) const
{
  ostr << "system " << com_;
}

void ProcessSystem::rawexec(DataStore&) const
{
  char buffer[MAXLINE];
  substitute_string(buffer,sizeof(buffer),com_.c_str(),SUB_ESCAPE);
  system(buffer);
}

void verify_index_list(const BaseList<size_t>& sel, size_t max, const char* name)
{
  const size_t badc=find_bad_index(sel,max);
  if (badc) {
    parser_printthread(std::cerr) << "Bad " << name << " index (" << (badc+1) << ") when data set has " << max << ' ' << name << "s\n";
    error_abort();
  }
}

void get_selection(LIST<complex>& dest, const BaseList<complex>& source, const BaseList<size_t>& inds)
{
  verify_index_list(inds,source.size(),"data point");
  dest=source(inds);
}

void ProcessExtract::exec(cmatrix& data, LIST<processing_state>& pflags) const
{
  const size_t oldcols=data.cols();
  const size_t oldrows=data.rows();
  if (!col_sel.empty())
    verify_index_list(col_sel,oldcols,"column");
  if (oldrows>1) {
    if (row_sel.empty())
      tmp=data(range(),col_sel); //!< don't need to check for empty col_sel
    else {
      verify_index_list(row_sel,oldrows,"row");
      if (col_sel.empty())
	tmp=data(row_sel,range());
      else
	tmp=data(row_sel,col_sel);
    }
    data.swap(tmp);
    if (!(row_sel.empty()))
      scaleswref(pflags.front(),row_sel.front(),data.rows(),oldrows,row_noncontiguous);
  }
  else {
    if (!row_sel.empty())
      error_abort("extract: 1D data set but 2D extract");
    tmpr=data.row()(col_sel);
    data.create(1,tmpr.size(),tmpr.vector());
  }
  if (!(col_sel.empty()))
    scaleswref(pflags.back(),col_sel.front(),data.cols(),oldcols,col_noncontiguous);
}

void check_extract(cmatrix& tmp_set)
{
  cmatrix tmp;
  static const char* synstr=NMRSIM_ROWCOLSTR;
  const LIST<size_t> sel1(parse_unsignedintarray_syntax(synstr,1,1));

  if (parser_isnormal()) {
    //    if (tmp_set.rows()==1)
    //  error_abort("2D submatrix requested but data is only 1D");
    verify_index_list(sel1,tmp_set.rows(),"row");
    const LIST<size_t> col_sel(parse_unsignedintarray_syntax(synstr,2,1));
    if (col_sel.empty()) {
      if (sel1.empty()) //!< nothing to do
	return;
      tmp=tmp_set(sel1,range());
    }
    else {
      verify_index_list(col_sel,tmp_set.cols(),"column");
      if (sel1.empty())
	tmp=tmp_set(range(),col_sel);
      else
	tmp=tmp_set(sel1,col_sel);
    }
  }
  else {
    if (sel1.empty())
      return;
    verify_index_list(sel1,tmp_set.cols(),"column");
    tmp=tmp_set(range(),sel1);
  }
  tmp_set.swap(tmp);
}
