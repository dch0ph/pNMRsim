#ifndef NMRsim_logfile_h_
#define NMRsim_logfile_h_

/*! \file 
  \brief  Objects for working with log files */

#include "matlabio.h"

//! class to control logging output
class logfile_controller {
public:
  logfile_controller(const char* fname, int flags); //!< open log file \a fname with flags \a flags
  ~logfile_controller();

  //! flags for log file
  enum { DOUBLE=1, //!< double precision output
	 MATLAB=2, //!< use Matlab format
	 ASCII=4, //!< use ASCII output
	 APPEND=8 //!< append to existing file
  };

  template<class T> void write(const Matrix<T>& a, const char* name) const; //! write matrix named \a name to log 
  template<class T> void write(const BlockedMatrix<T>& a, const char* name) const; //! write blocked matrix named \a name to log
  template<class T> void write(const ListList<T>& a, const char* name) const; //!< write list of lists to log
  template<class T> void write(const BaseList<T>& a, const char* name) const; //!< write list to log

  template<class T> bool try_write_matlab(const T& a, const char* name) const
    {
      if (filep_matlab) {
	if (filep_composite)
	  filep_composite->write(a,name);
	else
	  filep_matlab->write(a,name);
	return true;
      }
      return false;
    }

  void write(const char*, const char* name) const; //!< write string to log
  void flush() const; //!< ensure log file stream is flushed

  static void parse_raw(char*&, int&); //!< \internal
  static void maketitle(char*, int, const char*); //!< \internal
 friend class smartptr<logfile_controller,false>; //!< hack to allow use of smartptr when constructor is private

  static ContextWarning<> ignoring_warning; //!< logging commands will be ignored
  static ContextWarning<> badclose_warning; //!< close request with no open file

  bool supports_composite() const { return (filep_matlab!=NMRSIM_NULL); }
  
  void open_composite(const char*);
  void close_composite();

  class composite_guard {
  public:
    composite_guard(logfile_controller&, const char*);
    ~composite_guard();
  private:
    logfile_controller& logcon;
  };

  static ThreadWarning<> unexpectedcompositeclose_warning;

private:
  logfile_controller(const logfile_controller&); //!< prevent copy and assignment
  logfile_controller& operator= (const logfile_controller&);

  int flags; //!< log file flags
  matlab_controller* filep_matlab; //!< Matlab file controller (NMRSIM_NULL if not Matlab)
  matlab_controller::composite* filep_composite; //!< used for substructures - NMRSIM_NULL if not used;

  FILE* filep_ascii; //!< ASCII file pointer (NMRSIM_NULL if not ASCII)
};

void reset_logcount(); //!< reset counters used to tag log output
void update_log_file(const char* fname, int flags); //!< change log file to \a fname
//extern smartptr<logfile_controller,false> log_filep;
logfile_controller* get_logfile(); //!< pointer to log file controller object (NMRSIM_NULL if non active)

#endif
