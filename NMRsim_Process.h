#ifndef NMRSIM_PROCESS_H_
#define NMRSIM_PROCESS_H_

/*! \file
 \brief  Processing directives
*/

#include "NMRsim.h"
 //! Need std::list for procstack
#include <list>
#include "Parser.h"
#include "Lineshapes.h"
#include "simpsonio.h"

extern ContextWarning<> process_groupsvary_warning; //!< number of processing groups varies between stages
extern ThreadWarning<> process_nD_warning; //!< nD processing not fully supported

//! (abstract) base class for processing directives
class ProcessCommand : public Setable {
public:
  enum { PROC_HAS1D=1, //!< directive has 1D specialisation
	 PROC_HAS2D=2, //!< directive has 2D specialisation
	 PROC_HASBOTH=3, //!< directive has both 1D and 2D versions
	 PROC_CHANGES_ROWS=4, //!< directive changes no. of columns (split 13/8/15)
	 PROC_VARIABLE=8, //!< directive is not fixed
	 PROC_CHANGES_COLUMNS=16 //!< directive changes no. of rows (indirect dimensions)
  };

  enum {FT_INV=1, //!< inverse Fourier transform
	FT_NOFIRST=2, //!< no first point scaling
	FT_NOSHIFT=4, //!< no half spectrum shifting
	FT_FID=8, //!< from FID
	FT_SPECTRUM=16 //!< from spectrum
  };

  enum complexity_t { 
    REAL =1, //!< set real component
    IMAG =2, //!< set imaginary component
    COMPLEX =3, //!< set both
    COMPLEXPAIR =4 //!< *read* as pair
  };
  static complexity_t parsecomplexity(complexity_t, bool =false);
  static size_t checkcomplexity(size_t n, complexity_t);

  ProcessCommand(int flagsv) : flags_(flagsv) {}
  virtual ~ProcessCommand() {};

  friend class procstack_t;

  virtual void printvariablename(std::ostream&, subsid_t) const;
  virtual void print(std::ostream&) const =0;

  static ContextWarning<> rowchange_warning; //!< row change likely to cause problems
  static ContextWarning<> fidspectrum_warning; //!< -fid/-spectrum deprecated

  static void checkrowchangeok(); //!< warn if changing row size is likely to cause problems
  void isconstant(bool); //!< set const status
  void flagvariable() { flags_|=PROC_VARIABLE; } 
  int flags() const { return flags_; } //!< return flags
  double time() const { return stopwatch(); } //!< return accumulated calculation time
  usage_t usage() const { return usage_t(); } //!< processing commands do not reserve memory

  static ThreadWarning<> varignored_warning; //!< variable quantity ignored for 2D processing
  static ThreadWarning<> oneDprocessing_warning; //!< warn if 2D-capable processing applied to data set with 1 data point per row
  static ContextWarning<> sizechange_warning;

  //! raw processing entry
  /** This can be over-ridden for shape dependent processing */
  virtual void rawexec(DataStore&) const;

  virtual void exec(Matrix<complex>&, LIST<processing_state>&) const; //!< nD processing entry
  virtual void exec(BaseList<complex>, processing_state&) const; //!< row-by-row 1D processing entry
  virtual void exec(LIST<complex>&, const BaseList<complex>&, processing_state&) const; //!< entry for 1D that potentially changes size

  virtual void set(double,subsid_t);

  static ProcessCommand* create(); //!< parse to create new ::ProcessCommand
  static size_t get_ftflags(); //!< get common FT flags
  static void print_ftflags(std::ostream&, size_t); //!< dump FT flags

  mutable accumulating_timer<> stopwatch; //!< accumulates total calculation time
 
protected:
  int flags_; //!< flags
};

inline std::ostream& operator<< (std::ostream& ostr,const ProcessCommand& a)
{
  a.print(ostr);
  return ostr;
}

//! reverse sign on FID or spectrum
class ProcessSignReverse : public ProcessCommand {
public:
  explicit ProcessSignReverse(bool ifneggammav) : 
    ProcessCommand(PROC_HAS1D),
    ifneggamma_(ifneggammav) {}
  void print(std::ostream&) const;
  void exec(BaseList<complex>, processing_state&) const;

private:
  bool ifneggamma_; //!< if \c true, only reverse sign if detection nucleus gamma is negative
};

//! \c echo in process block
struct ProcessEcho : public ProcessCommand, public BaseEcho  {
public:
  ProcessEcho(const char* str_);
  void exec(cmatrix&, LIST<processing_state>&) const { BaseEcho::exec(); }
  void exec(BaseList<complex>, processing_state&) const { BaseEcho::exec(); }
  void print(std::ostream& ostr) const { BaseEcho::print(ostr); }
  static ProcessCommand* create();
};

class ProcessIntegrate : public ProcessCommand {
 public:
  ProcessIntegrate(UserVariable& destv, const BaseList<size_t>& rangesv);

  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream&) const;
  static ProcessCommand* create();
  static ContextWarning<> overlap_warning; //!< integration ranges overlap
  
 private:
  UserVariable& dest_; //!< destination variable
  LIST<size_t> ranges_;
  size_t maxindex_;
  mutable LIST<double> tmp_;
};
    
//! \c log_file directive
struct ProcessLogFile : public ProcessCommand {
public:
  ProcessLogFile(const char* name_, int flagsv)
    : ProcessCommand(PROC_HAS2D), name(name_ ? name_ : ""), logflags(flagsv) {}
  //  void exec(cmatrix&, BaseList<processing_state>) const;
  void rawexec(DataStore&) const;
  void print(std::ostream& ostr) const
  { ostr << "log_file " << name << '\n'; }
  static ProcessCommand* create();

private:
  std::string name; //!< log file name
  int logflags; //!< flags
};

class ProcessSave;

 //! helper structure for ::ProcessSave (summary of input state)
struct ProcessSave_state {
  //  ProcessSave_state(const DataStore& datav, const char* fnamev, int rowv =-1)
//     : data(datav), fname(fnamev), row(rowv) {}
  ProcessSave_state(const char* fnamev) : fname(fnamev) {}
//   ProcessSave_state(const DataStore& datav, const ProcessSave_state& a) //!< change data set
//     : data(datav), fname(a.fname), row(a.row) {}
//   const DataStore& data; //!< reference to data
//   //  processing_state flags; //!< input flags
   const char* fname; //!< base filename
//   int row; //!< data set row (-1 for complete data set)
};

typedef void (*SaveCommand_function)(const ProcessSave&, const ProcessSave_state&, int); //!< save callback
void register_save_function(const char*, SaveCommand_function); //!< register new save callback

// maximum length of Matlab variables
#define NMRSIM_MATLAB_VARMAX 63

extern DataStore fit_set; //!< fit_set needs to be externally visible for ReverseGuard

class ReverseGuard {
public:
  ReverseGuard(DataStore& FIDv, const char* labelv =NMRSIM_NULL, bool dorevv =true);
  ~ReverseGuard(); //!< flip back when guard drops out of scope
  bool isflipped() const { return flipped_; }
  void flip();

private:
  ReverseGuard(const ReverseGuard&);
  ReverseGuard& operator= (const ReverseGuard&); //!< declare these as private to avoid copy
      
  DataStore& FID_;
  bool flipped_;
  const char* label_;
};

//! \c save directive
class ProcessSave : public ProcessCommand {
public:
  enum { SIMPSON=1, //!< SIMPSON format
	 MATLAB=2, //!< Matlab format
	 ASCII=4, //!< ASCII format
	 NODATA=8, //!< don't save data
	 PROJ=16, //!< save projection (sum of columns)
	 SUM=32, //!< sum sums of rows
	 SCALE=64, //!< save scale
	 SIMPLOT=128, //!< SIMPSON 1D only
	 VARIABLES=256, //!< save array variables
	 SOURCE=1024, //!< save input file in comments
	 STOPOVERWRITE=2048, //!< prevent file overwrite
	 NOSFRQ=4096, //!< disable SFRQ
	 POSITIVESFRQ=8192, //!< force SFRQ to be positive
	 REALONLY=16384, //!< save real values only
	 REVERSEFREQ=32768, //!< reverse frequency axis
	 USERFLAGS=65536 //!< start of additional flags
  };
   
  ProcessSave(const char*,const BaseList<RealVariable*>& =BaseList<RealVariable*>(), int =SIMPLOT);

  void rawexec(DataStore& FID) const; //!< override normal exec
  void exec(BaseList<complex>, processing_state&) const; //!< single row save
  void print(std::ostream&) const;
  static ProcessCommand* create();

  bool compositesave() const { return (saveflags_ & MATLAB); } //!< \c true if data saved in single composite file
  template<class T> void save(const T& data, const ProcessSave_state&, const char* qual) const; //!< save raw data object
  void save(const DataStore&, const ProcessSave_state&, const char* qual =NMRSIM_NULL) const; //!< save data set
  void save(const BaseList<complex>& data, const processing_state&, const ProcessSave_state&, const char* qual) const; //!< save row of data set
  bool checkoverwrite(const char*) const; //!< warn and return \c true (or raise error) if overwriting file (if STOPOVERWRITE set)
  static bool defaultstopoverwrite; //!< default stopoverwrite state
  int flags() const { return saveflags_; }
  
  static ThreadWarning<> simpsonunsuitable_warning; //!< SIMPSON/PLOT format unsuitable for data set
  static ThreadWarning<> notnd_warning; //!< File format does not support >2D data
  static ThreadWarning<> notirregular_warning; //!< ASCII format does not support irregalar data
  static ThreadWarning<> nosource_warning; //!< File format does not support source save
  static ThreadWarning<> noscale_warning; //!< scale save not possible
  static ThreadWarning<> sumprojection_warning; //!< sum/projection save requires nD data
  static ThreadWarning<> sumprojection_warning2; //!< sum/projection save for single row
  static ThreadWarning<> noparameters_warning; //!< no active parameters to save
  static ThreadWarning<> noscaleforirregular_warning; //!< can't save scale for irregular data
  static ContextWarning<> nothingtosave_warning; //!< nothing to save
  static ContextWarning<> sfrq_warning; //!< sfrq flags used with inappropriate file format
  static ContextWarning<> nosimpsonrevfreq_warning; //!< can't use -reversefrequency with SIMPSON/SIMPLOT
  static ContextWarning<> revfreqscale_warning; //!< warn about combining -reversefrequency with -scale
  static ThreadWarning<> sumparameters_warning; //!< -parameters save outside summation array
  static ThreadWarning<> invalidmatlab_warning; //!< invalid MATLAB variable name
  static ThreadWarning<> realnotvalid_warning; //!< -realonly not valid for this format
  static ThreadWarning<> nonconstvariable_warning; //!< trying to save non-constant variable
  static ThreadWarning<> emptyvariable_warning; //!< trying to save empty variable
  static ThreadWarning<> simpson_notfid_warning; //!< trying to save spectrum to file with .fid extension
  static ThreadWarning<> simpson_notspec_warning; //!< trying to save spectrum to file with .spe extension

  const char* filename() const { return fname.vector(); }

  bool original_isflipped() const;

private:
  LIST<char> fname; //!< filename
  const LIST<RealVariable*> vars_; //!< set of variable names to save
  int saveflags_; //!< flags
  mutable ScratchList<char> scr; //!< scratch space for filename
  LIST<SaveCommand_function> saveoptions; //!< additional options (from USERFLAGS)
  mutable ReverseGuard* original_revguardp;

  bool openmatlab(const char*) const;
  void writevars(const ProcessSave_state&, bool currentrow =false) const; //!< write variable values
  void writeparametersoptions(const ProcessSave_state&, int row =-1) const; //!< write array parameters and call any optional saves (must be last in sequence)

  void makescale(LIST<double>& scale, const processing_state& pflags, size_t np, size_t skip) const;
  void rawexec_postrev(DataStore&) const; //!< after any flipping of data set

  const char* makefilename(const char* base, const char* qual) const; //!< make qualified filename from \a base and qualifier \a qual (e.g. "residuals")
  // void rawexec_(const DataStore&, const char* base, const char* qual) const; //!< \internal raw save of ::DataStore

  FILE* opencomments(const char* outname, const char* mode ="a") const; //!< create file pointer for comment/source writing (NMRSIM_NULL for Matlab output)
  void writeline(FILE*, const char*, char comchar ='#') const; //!< write comment/source line
  void closecomments(FILE*, const char* name) const; //!< close comment stream
  void savesource(const char* source) const; //!< write file \a source to comment stream
  const char* sanitise_varname(const char*, char buf[NMRSIM_MATLAB_VARMAX+1]) const; //!< cleanup variable name (MATLAB) 
};

//! \c phase directive
class ProcessPhase : public ProcessCommand {
public:
  ProcessPhase(double zero_, double first_, double pivot_)
    : ProcessCommand(PROC_HAS1D),
    zero(zero_), first(first_), pivot(pivot_) {}
  
  void set(double, subsid_t);
  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  static ProcessCommand* create();

private:
  double zero; //!< zero order correction
  double first; //!< first order correction
  double pivot; //!< first order pivot
};

//! \c apply directive
class ProcessApply : public ProcessCommand {
public:
  typedef const ExpressionNamedBase& function_spec_t;
  ProcessApply(function_spec_t, complexity_t, size_t =1U, const BaseList< LIST<double> >& =BaseList< LIST<double> >(), bool isconst =true);

  void set(double v, subsid_t subsid) { set(BaseList<double>(1,&v),subsid); }
  void set(const BaseList<double>&, subsid_t);
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  static ProcessCommand* create();
  bool ismultiargument() const { return (argn_>1); }
  void exec(LIST<complex>&, const BaseList<complex>&, processing_state&) const;
  void exec(BaseList<complex>, processing_state&) const;

  static ContextWarning<> singleargument_warning;
private:
  const char* fname_;
  complexity_t complexity_;
  size_t argn_;
  mutable LIST< LIST<double> > argvals_;
  Expression expr_;
  mutable LIST<double> tmp_,bufr_,bufc_;
  void zip(LIST<complex>&, const BaseList<double>&, const BaseList<double>&) const;
  void apply(LIST<double>&, const BaseList<double>&) const;
};

//! \c offset directive
class ProcessOffset : public ProcessCommand {
public:
  ProcessOffset(VariableBase& offsetv, complexity_t addtypev, bool isconst)
    : ProcessCommand(PROC_HAS1D | (isconst ? 0 : PROC_VARIABLE)),
      offsetvar(offsetv),
      addtype(addtypev)
  {}
  void set(double offsetv, subsid_t) { offsetvar.set(offsetv); }
  void set(const BaseList<double>& offsetl, subsid_t) { offsetvar.set(offsetl); }
  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  static ProcessCommand* create();

private:
  VariableBase& offsetvar; //!< pointer to ::VariableBase holding offset (which may be vector)
  complexity_t addtype; //!< \c component to add to
};

//! \c rev directive
struct ProcessReverse : public ProcessCommand {
  ProcessReverse() : ProcessCommand(PROC_HAS1D) {}
  void print(std::ostream& ostr) const { ostr << "rev"; }
  static void exec(BaseList<complex>);
  void exec(BaseList<complex>, processing_state&) const;
  static ProcessCommand* create() { return new ProcessReverse(); }
  static ThreadWarning<> td_warning; //!< applied to time domain data
  //  static ThreadWarning<> refreset_warning; //!< referencing has been reset
};

//!< \c normalise directive
class ProcessNormalise : public ProcessCommand {
public:
  enum metric_t { MINMAX =1, INTEGRAL, ABS, AREA };

  ProcessNormalise(double, metric_t);

  void print(std::ostream&) const;
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  static ContextWarning<> multiplenormalise_warning; //!< multiple normalise processing directives
  static bool isdefined() { return havealready_; }
private:
  void exec(BaseList<complex>, double) const;
  static const flagsmap_type& get_flags();
  double val_;
  metric_t metric_;
  static bool havealready_;
};

//! \c extract directive
class ProcessExtract : public ProcessCommand {
public:
  ProcessExtract(const BaseList<size_t>& col_selv, const BaseList<size_t>& row_selv =BaseList<size_t>());
  void exec(LIST<complex>&, const BaseList<complex>&, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  static ProcessCommand* create();
  void print(std::ostream&) const;
  static ContextWarning<> nullselection_warning; //!< null selection will do nothing
  static ThreadWarning<> refreset_warning; //!< referencing has been reset

private:
  LIST<size_t> col_sel; //!< column selection
  LIST<size_t> row_sel; //!< row selection
  bool col_noncontiguous; //!< column selection is non-contiguous
  bool row_noncontiguous; //!< row selection is non-contiguous
  mutable cmatrix tmp; //!< scratch space
  mutable LIST<complex> tmpr; //!< scratch space
  static void scaleswref(processing_state&, size_t startn, size_t newn, size_t oldn, bool isnoncontig); //!< adjust sw and ref 
};

//! \c fill directive
class ProcessFill : public ProcessCommand {
public:
  ProcessFill(double val, complexity_t filltypev =COMPLEX, const BaseList<size_t>& col_selv =BaseList<size_t>(), const BaseList<size_t>& row_selv =BaseList<size_t>());
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void set(double vv, subsid_t) { val=vv; }
  static ProcessCommand* create();
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  //  static ContextWarning<> fillwillwipe_warning; //!< command will wipe out data set

private:
  double val; //!< fill value
  complexity_t filltype; //!< real/imag/complex selector
  LIST<size_t> col_sel; //!< column selection
  LIST<size_t> row_sel; //!< row selection  
  void exec_(BaseList<complex>&) const; //!< \internal
};

//! \c set data set directive
// class ProcessSet : public ProcessCommand {
// public:
//   ProcessSet(const Expression&, complexity_t filltypev =COMPLEX);
//   void exec(BaseList<complex>, processing_state&) const;
//   void exec(cmatrix&, BaseList<processing_state>) const;
//   static ProcessCommand* create();
//   void print(std::ostream&) const;
//   static ContextWarning<> constant_warning; //!< set to constant will wipe out data set

// private:
//   const Expression& expr; //!< values expression
//   complexity_t filltype; //!< real/imag/complex selector
//   void exec_(BaseList<complex>&) const; //!< \internal
// };
 
//! \c addnoise directive
class ProcessAddNoise : public ProcessCommand {
public:
  ProcessAddNoise(double sigmav, bool reseedv)
    : ProcessCommand(PROC_HAS1D),
      sigma(sigmav), reseed(reseedv) {}
  void set(double sigmav, subsid_t) { sigma=sigmav; }
  void exec(BaseList<complex> FID, processing_state&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  static ProcessCommand* create();

private:
  double sigma; //!< noise standard deviation
  bool reseed; //!< restart RNG each time if \c true
};

//! \c scale directive
class ProcessScale : public ProcessCommand {
public:
  ProcessScale(double scale_)
    : ProcessCommand(PROC_HASBOTH), scale(scale_) {}

  void set(double scale_, subsid_t =S_NONE) { scale=scale_; }

  void exec(cmatrix&, LIST<processing_state>&) const;
  void exec(BaseList<complex>, processing_state&) const;
  
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;

  static ProcessCommand* create();

private:
  double scale; //!< scale factor
};

//! \c scalefirst directive
class ProcessScaleFirst : public ProcessCommand
{
public:
  ProcessScaleFirst(double scale_, double scale1_ =1.0)
    : ProcessCommand(PROC_HASBOTH), scale(scale_), scale1(scale1_) {}

  void exec(BaseList<complex> FID, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream& ostr, subsid_t) const;
  void set(double, subsid_t);
  static ProcessCommand* create();
  static ThreadWarning<> fd_warning; //!< applied to frequency domain data

private:
  double scale; //!< first column scaling
  double scale1; //!< first row scaling
};

//! \c addsignals directive
class ProcessAddSignals : public ProcessCommand
{
public:
  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  void set(double, subsid_t);
  void set(const BaseList<double>&, subsid_t);
  static ProcessCommand* create();
  //  static Warning<> nothistogram_warning; //!< not in histogram mode
  static ThreadWarning<> ignoring_phase_warning; //!< non-zero phase parameter ignored
  static ContextWarning<> allzero_warning; //!< amplitudes all zero / missing
  static ThreadWarning<> timedomainref_warning; //!< ref ignored for time-domain data

  // static ThreadWarning<> histogramloss_warning; //!< lost frequencies from histogram

private:
  smartptr<VariableBase,false> freqs_; //!< frequencies list
  smartptr<VariableBase,false> amps_; //!< amplitudes list
  smartptr<VariableBase,false> phases_; //!< phases list
  mutable bool errwarn; //!< \internal
  ProcessAddSignals();
  template<class T> void set_(const T&, subsid_t);
};

//! \c resample directive
class ProcessResample : public ProcessCommand
{
public:
  ProcessResample(size_t,bool);
  ProcessResample(size_t,double,double =0.0);
  static ProcessCommand* create();
  void set(double, subsid_t);
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  void exec(LIST<complex>&, const BaseList<complex>&, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  
private:
  size_t newnp_; //!< new number of data points
  bool dofold_; //!< \c true if -fold specified
  double newsw_; //!< new spectral width (<0 if unspecified)
  double offset_; //!< spectrum centre offset
  void rawresample(BaseList<complex>, const BaseList<complex>&) const; //!< \internal
  static ThreadWarning<> refreset_warning; //!< referencing has been reset
};

//! \c addlb directive
class ProcessAddLB : public ProcessCommand
{
public:
  ProcessAddLB(double lb_,double gfrac_, double lb1_, double gfrac1_);
  void set(double, subsid_t);
  void exec(cmatrix&, LIST<processing_state>&) const;
  void exec(BaseList<complex>, processing_state&) const;
  void print(std::ostream&) const;
  void printvariablename(std::ostream&, subsid_t) const;
  static ProcessCommand* create();

private:
  double lb; //!< t2 line-broadening
  double gfrac; //!< guassian fraction (t2)
  double lb1; //!< t1 line-broadening
  double gfrac1; //!< gaussian fraction (t1)

  // enum { S_t1 =1, S_t2, S_gfrac, S_gfrac1 };

  double t2,t21,lbg,lbg1;
  
  static void set_(double&, double&, double&, double&, double, double); //!< \internal
  static void print_(std::ostream& ostr, double lb_, double gfrac_); //!< \internal
  void execf(BaseList<complex>&) const; //!< \internal

  mutable smartptr<LorentzianGaussian,false> shapegenp; //!< lorentzian/gaussian lineshape generator (created as needed)
  mutable LIST<double> shapebuf; //!< cache for lineshape

  static ThreadWarning<> narrow_warning; //!< linewidth is unrealistically narrow
};

//! \c transpose directive
class ProcessTranspose : public ProcessCommand {
public:
  ProcessTranspose() : ProcessCommand(PROC_CHANGES_COLUMNS | PROC_CHANGES_ROWS | PROC_HAS2D) {}
  void rawexec(DataStore&) const;
  void print(std::ostream& ostr) const { ostr << "transpose"; }
  static ProcessCommand* create();
};

//! \c zerofill directive
class ProcessZeroFill : public ProcessCommand {
public:
  ProcessZeroFill(int npts_, int npts1_ =1);
  void exec(cmatrix&, LIST<processing_state>&) const;
  void exec(LIST<complex>&, const BaseList<complex>&, processing_state&) const;
  void print(std::ostream&) const;
  static ProcessCommand* create();
  static ThreadWarning<> ignoring_warning; //!< zerofill in indirect dimension ignored for non-nD data
private:
  const int npts; //!< zerofill factor in t2
  const int npts1; //!< zerofill factor in t1
};

//! Helper object for 1D Fourier Transform
class FTObject {
public:
  FTObject(size_t, int);
  void exec(BaseList<complex>, const processing_state&) const;
private:
  size_t n;
  bool doshift,doscale,fastft;
  //  int dir;
  int flags;
  static void negate_alternate(BaseList<complex>); //!< negate alternate points
};

//! \c ft directive
class ProcessFT : public ProcessCommand {
public:
  ProcessFT(size_t flagsv)
    : ProcessCommand(PROC_HAS1D),
    ftflags(flagsv) {}

  void exec(BaseList<complex> FID, processing_state&) const;
  void print(std::ostream&) const;
  static ProcessCommand* create();
  static ThreadWarning<> doesnothing_warning;
  static ThreadWarning<> nonpowerof2_warning;
  static ThreadWarning<> expectingfid_warning;
  static ThreadWarning<> expectingspectrum_warning;
private:
  size_t ftflags; //!< flags
};

//! \c conjugate directive 
struct ProcessConjugate : public ProcessCommand {
  ProcessConjugate() : ProcessCommand(PROC_HASBOTH) {}
  static void exec(BaseList<complex>& FID) { conj_ip(FID); }
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream& ostr) const { ostr << "conjugate"; }
  static ProcessCommand* create() { return new ProcessConjugate(); }
};

struct ProcessMagnitude : public ProcessCommand {
  ProcessMagnitude() : ProcessCommand(PROC_HASBOTH) {}
  static void exec(BaseList<complex>&);
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream& ostr) const { ostr << "magnitude"; }
  static ProcessCommand* create() { return new ProcessMagnitude(); }
};

//! \c shifthalf directive
struct ProcessShiftHalf : public ProcessCommand {
  ProcessShiftHalf() : ProcessCommand(PROC_HASBOTH) {}
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream& ostr) const { ostr << "shifthalf"; }
  static ProcessCommand* create() { return new ProcessShiftHalf(); }
  //  static ThreadWarning<> refreset_warning; //!< referencing has been reset
};

//!< \v set directive
class ProcessSet : public ProcessCommand {
public:
  ProcessSet(double swv, double sw1v, double refv, double ref1v, double sfrqv, double sfrq1v, int whichv)
    : ProcessCommand(PROC_HASBOTH),
      sw_(swv), sw1_(sw1v),
      ref_(refv), ref1_(ref1v),
      sfrq_(sfrqv), sfrq1_(sfrq1v),
      which_(whichv) {}
  ProcessSet()
    : ProcessCommand(PROC_HASBOTH), sw_(0.0), sw1_(0.0),
      ref_(0.0), ref1_(0.0),
      sfrq_(0.0), sfrq1_(0.0),
      which_(0) {}

  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream&) const;
  void set(double, subsid_t);
  static ProcessCommand* create();
  enum { SW =1, SW1 =2, REF =4, REF1 =8, SFRQ =16, SFRQ1 =32};
  static ThreadWarning<> notnd_warning; //!< attempting to set nD paras in non-nD data set
  static ContextWarning<> repeatset_warning; //!< overriding previous set
private:
  double sw_,sw1_,ref_,ref1_;
  double sfrq_,sfrq1_;
  int which_;
  static const int indmask_;
  static flagsmap_type flags_;
  void rawset(processing_state&, size_t dim =0U) const;
};

//!< \c setdomain directive
class ProcessSetDomain : public ProcessCommand {
public:
  ProcessSetDomain(size_t dimv, size_t flagsv)
    : ProcessCommand(PROC_HASBOTH), dim_(dimv), flags_(flagsv) {}
  void exec(BaseList<complex>, processing_state&) const;
  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream&) const;
  static ProcessCommand* create();
  enum { SWITCH =1, TIME =2, FREQUENCY =4, STATES =8, NOSTATES =16 };
  static ContextWarning<> fd_warning;
private:
  size_t dim_;
  size_t flags_;
  static const flagsmap_type& getflags();
  static const size_t statesmask = STATES | NOSTATES;
  static const size_t dommask =SWITCH | TIME | FREQUENCY;
  void exec_(processing_state&) const;
};

//! \c ft2d directive
class ProcessFT2d : public ProcessCommand {
public:
  ProcessFT2d(size_t flagsv =0);
  ProcessFT2d(double zero2_, double first2_, double zero1_, double first1_, size_t ftflagsv =0);

  void exec(cmatrix&, LIST<processing_state>&) const;
  void print(std::ostream&) const;
  static ProcessCommand* create();
  static ThreadWarning<> invalidskip_warning; //!< skip factor is not 1 or 2
private:
  size_t ftflags; //!< flags
  const bool needphase; //!< apply phase correction if \c true
  double zero2,first2,zero1,first1;

  void create_();
};

class ProcessSystem : public ProcessCommand
{
 public:
  ProcessSystem(const char* comv)
    : ProcessCommand(PROC_HAS2D), com_(comv) {}
    void rawexec(DataStore&) const; //!< override normal exec
    void print(std::ostream&) const;
    static ProcessCommand* create();
 private:
    std::string com_;
};

//! stack for processing actions
struct procstack_t : public std::list<ProcessCommand*> {
  procstack_t() {}
  void initialise() {}
  void exec(DataStore&) const; //!< apply processing
};

void apply_procstacks(const BaseList<procstack_t>&, DataStore&); //!< apply processing stack to data set

bool process_command(const std::string&, char*, Variable* =NMRSIM_NULL); //!< parse processing directive

typedef ProcessCommand* (*ProcessCommand_function)();
typedef FASTMAPTYPE(ProcessCommand_function) Process_Factory_t;
Process_Factory_t& get_Process_Factory(); //!< return registry for processing directives
Process_Factory_t& get_Finalise_Factory(); //!< return registry for objects in \c finalise block

extern LIST<procstack_t> procstacks; //!< per orientation processing actions;
extern LIST<procstack_t> postprocstacks; //!< total spectrum processing actions;
extern LIST<procstack_t> finalisestacks; //!< end-of-calculation actions

size_t read_proc_blocks(LIST<procstack_t>& pstacks, size_t& accstacks, const char* name, const Process_Factory_t& factory, bool allowsum =false); //!< parse (set of) processing block(s)

bool read_proc_block(procstack_t& pstack, const char* name, const Process_Factory_t& factory, bool allowsum =false); //!< parse single processing block

enum domain_t { D_UNKNOWN=0, D_FREQUENCY=1, D_TIME=2 };
struct filestruct {
  filestruct() { clear(); }
  void clear();
  domain_t domain;
  domain_t domain1;
  double sw;
  double sw1;
  double sfrq;
  double sfrq1;
  double ref;
  double ref1;
  filestruct& operator= (const libcmatrix::simpsonFD&);
};
void raw_read_file(cmatrix&, filestruct&, const char*, LIST<std::string>* =NMRSIM_NULL);
void raw_build_set(DataStore&, const BaseList<cmatrix>&, const BaseList<processing_state>&, int verbose =0);
void process_set_n(int,int,int =1);
extern option optsaveprotect; //!< overwrite protection
void get_selection(LIST<complex>&, const BaseList<complex>&, const BaseList<size_t>&);
bool try_merge_proc_blocks(LIST<procstack_t>& source, LIST<procstack_t>& dest); //!< try to combine - return \c true if successful
void rawreverse(BaseList<complex>);
void verify_index_list(const BaseList<size_t>& sel, size_t max, const char* name);

#endif
