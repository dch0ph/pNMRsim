#ifndef NMRsim_MasterObj_h_
#define NMRsim_MasterObj_h_

/*! \file 
  \brief  Header file required by components that are actively involved in simulation */

#include "NMRsim.h"
#include "NMRsim_matrixspec.h"
#include "Propagation.h"
#include "AsynchronousPropagation.h"
#include "MetaPropagation.h"
#include "Sequence.h"

#include "powder.h"
extern range_t rangequal; //!< powder range restriction

//! propagation method types
enum prop_t {
  // PROP_AUTO, //!< automatically choose propagation (not implemented)
       PROP_DIAG, //!< always use diagonalisation
       PROP_CHEBYSHEV //!< always use Chebyshev
};
extern prop_t prop_method; //!< propagation method (a la ::prop_t)

extern LIST<int> prop_flags; //!< propagation flags for each block structure
extern LIST<matrix_partition_set> partitions; //!< ``partitioning'' for each block structure

extern dirty_t update_propagators; //!< dirty status of overall simulation
std::ostream& operator<< (std::ostream& ostr, dirty_t);

enum acc_t {
  ACC_TIME=0, //!< simple time domain accumulation
  ACC_FREQ, //!< frequency domain accumulation
  ACC_PSEUDO //!< frequency domain mode, time domain output
};
extern acc_t acc_mode; //!< accumulation mode

//! Hamiltonian structure types
enum ham_t { H_DIAGONAL, //!< Hamiltonian is diagonal in Zeeman basis
	     H_REAL, //!< Hamiltonian is matrix, but purely real
	     H_COMPLEX, //!< General case of complex Hamiltonian
	     H_NULL//!< Hamiltonian undefined
};

//! Hamiltonian object type
enum hamobj_t { M_NULL, //!< unset
		M_STATIC, //!< static (matrix)
		M_STATICD, //!< static, diagonal
		M_SPIN, //!< spinning (matrix)
		M_SPIND //!< spinning, diagonal
};

//! object defining instantanteous Hamiltonian state - at the moment just current rotor_angle
typedef double hamiltonian_state_t;
//std::ostream& operator<< (std::ostream& ostr, const hamiltonian_state_t&);

class ActionDirectAcq;
class ActionAcqPoint;

/**
   MasterObj contains, in principle, all the "state" about a simulation.
   Originally pNMRsim supported multi-threaded operation (with a
   separate MasterObj per thread).  As pNMRsim only supports MPI, with
   a separate program being run per processor, MasterObj is not strictly
   necessary, although a well-defined "simulation" object would make sense
   in a clean API.
*/
class MasterObj
#ifdef HAVE_PARSYS
  : public BaseThreadFunction //!< In parallel executation, MasterObj must be derived from \c BaseThreadFunction
#endif
 {
 public:
//! initialise from spin system, Hamiltonian specification and list of blocked nuclei
   MasterObj(const basespin_system& sys,const HamiltonianStore<space_T>& Hstorev,const BaseList<nuclei_spec>& blockingnuc,int flags =0, bool =false);
   MasterObj(); //!< constructor in absence of spin system
  
   bool check_powder(const HamiltonianStore<space_T>&, const Euler&) const; //!< check whether Hamiltonian really is invariant if powder range is restricted
   
   BlockedOperator make(const productoperator_spec&) const; //!< create operator from product operator specification
  
  //! accumulate Hamiltonian using appropriate OpGen
  template<class T,class M> void add_Hamiltonian(T& H,const HamiltonianStore<M>& Hstore) const 
  {
#ifndef NOPERIODIC
    if (!simple_opgenp)
      ::libcmatrix::add_Hamiltonian(H,*crystal_opgenp,Hstore);
    else
#endif
      ::libcmatrix::add_Hamiltonian(H,*simple_opgenp,Hstore);
  }

  //! Update operator generators following "change"
  /** \note It is not very obvious \e when an update is required
      (Currently only after changing the proton Larmor frequency for >1 order quadrupoles) */
   void update();
   void ensure_channels(); //!< ensure RF pulse generators have been created

  const SpinOpGeneratorBase& spinop_generator() const {
#ifndef NOPERIODIC
    if (!simple_opgenp)
      return *crystal_opgenp;
#endif
    return *simple_opgenp;
  }
   
   smartptr<HamiltonianStructure> Hstructp; //!< Block structure of problem
   smartptr<SpinOpGenerator> simple_opgenp; //!< (normal) spin operator generator
#ifndef NOPERIODIC
   smartptr<CrystalOpGenerator> crystal_opgenp; //!< periodic spin operator generator initialised if required
#endif
   bool canbuildoperators(const BaseList<productoperator_spec*>&) const; //!< return \c false if can't build list of these spin operators

   BlockedMatrix<complex> U; //!< temporary propagator store
   RotorInfo rinfo; //!< current rotor state
   double rotor_angle; //!< current rotor angle (may change during sequence)
   
   mutable cmatrix dUs; //!< propagator stores
 
   usage_t usage() const; //!< return (total) memory usage
   usage_t Husage() const; //!< return memory usage for Hamiltonian

   bool isspinningH() const { return (hamobjtype==M_SPIN); }

   typedef RealComplexHolder< smartptr< BlockedStaticHamiltonian<double> >, smartptr< BlockedStaticHamiltonian<complex> > > Hstatic_t;
   Hstatic_t Hstaticp;  //!< static Hamiltonian store

   typedef RealComplexHolder< smartptr< BlockedSpinningHamiltonian<double> >, smartptr< BlockedSpinningHamiltonian<complex> > > Hspin_t;
   Hspin_t Hspinp; //!< spinning Hamiltonian store

   bool iscomplex() const { return (Hstaticp.iscomplex() || Hspinp.iscomplex()); } //!< \c true if Hamiltonian is complex
   smartptr<BlockedDiagonalStaticHamiltonian> Hstaticdp; //!< static diagonal Hamilonian store
   smartptr<BlockedDiagonalSpinningHamiltonian> Hspindp; //!< spinning diagonal Hamiltonian store
   
  typedef smartptr<AsynchronousFID,false> AFID_t;
  typedef smartptr<GammaInhomogeneousFID,false> GIFID_t;
  typedef smartptr<InhomogeneousFID,false> IFID_t;
  typedef smartptr<GammaPeriodicFID,false> GPFID_t;
  typedef smartptr<PeriodicFID,false> PFID_t;
  typedef smartptr<StaticFID_H,false> SFIDH_t;
  typedef smartptr<StaticFID_U,false> SFIDU_t;

   UnionHolder< 7, AFID_t, GIFID_t, IFID_t,GPFID_t, PFID_t, SFIDH_t, SFIDU_t> objholder; //!< pointer to one of FID generator types

   typedef smartptr<GammaInhomogeneousSpectrumED,false> GISED_t;
   typedef smartptr<GammaInhomogeneousSpectrum,false> GIS_t;
   typedef smartptr<InhomogeneousSpectrum,false> IS_t;
   typedef smartptr<GammaPeriodicSpectrumED,false> GPSED_t;   
   typedef smartptr<GammaPeriodicSpectrum,false> GPS_t;   
   typedef smartptr<PeriodicSpectrum,false> PS_t;   
   typedef smartptr<StaticSpectrumED,false> SSED_t;   
   typedef smartptr<StaticSpectrum,false> SS_t;   

   UnionHolder<8,GISED_t,GIS_t,IS_t,GPSED_t,GPS_t,PS_t,SSED_t,SS_t> fobjholder; //!< pointer to one of spectrum generator types
  
   void propagator(BlockedMatrix<complex>&, double t1, double t2); //!< calculate propagator for time interval \a t1 to \a t2 (no RF)
   template<class SequenceHolder> void propagator(BlockedMatrix<complex>&, const SequenceHolder& seq, double t1,double t2,double); //!< calculate propagator for time interval \a t1 to \a t2 with given RF sequence

   void add_FID_orientation(DataStore&, size_t i); //!< add orientation with index \a i to data set
   void raw_add(const ActionDirectAcq&, BlockedOperator&, double&, int =-1); //!< do acquisition given current density matrix, time and no. of data points
   void raw_add(const ActionAcqPoint&, BlockedOperator&); //!< acquire single data point without propagation
   void calc(DataStore&, int verbose, bool share =true); //!< calculate complete data set.  If \a share then share result with all processes
   //   void expandoperator(cmatrix&, const BlockedOperator&) const; //!< expand blocked operator into full matrix

   /** 
       \todo These "reset" functions need a real clean up ... */
   //void initrun(); //!< reset at start of new array (inside sum array)
   //   void reset(); //!< reset at start of new orientation
   void initialise_simple_simulation(); //!< initialise simple (1D, const) similation

   //   void exec_action(); //!< apply "pulseq" actions

   void printH(std::ostream& =std::cout, int structureflags =0) const; //!< output current Hamiltonian
   void printfullH(std::ostream& =std::cout) const; //!< stream Hamiltonian (full form)
   void logH(logfile_controller&, const char*, int structureflags =0) const; //!< output Hamiltonian to current "log"

   bool reset_angle(double); //!< set rotor angle
   bool change_angle(double); //!< change angle (during sequence)
   void rotor_phase(double); //!< change rotor phase
   //! return rotor phase at time \a t
   /** Urghh! Clash with ::rotor_phase results in inconsistent naming */
   double phase(double t) const;
   
   hamiltonian_state_t hamiltonian_state() const { return hamiltonian_state_t(rotor_angle); }

   void recreate_Hs(); //!< create/update Hamiltonians
   bool Hempty() const; //!< \c true if Hamiltonian has not be created
   
   //   BlockedFilter* create_filter(const filter_spec&) const; //!< create new coherence filter
   BlockedFilter* create_spinorder(const spinorder_spec&, const block_pattern&) const; //!< create new spinorder / coherence filter

 private:

   void operator()(size_t st,size_t end,size_t nthr) const; //!< accumulate a range of orientations (interface to \c BaseThreadFunction)

   void restart();  //!< reset at start of new row / simple calculation

   MasterObj(const MasterObj&); //!< can't copy objects (due to pointers)
   MasterObj& operator= (const MasterObj&); 
   
   template<class HType> void spin_propagator_(BlockedMatrix<complex>& Udest, const HType& H, const CompSequenceBase& seq, double t1,double t2,double origin);
   template<class HType> void spin_propagator_(BlockedMatrix<complex>& Udest, const HType& H, const Sequence& seq, double t1,double t2,double origin);

   void recreate_Hs(const SpinOpGenerator&, double);
#ifndef NOPERIODIC
   void recreate_Hs(const CrystalOpGenerator&, double);
#endif
   void add_freq(BaseList<complex>, complex, bool, const ActionDirectAcq&);
   template<class PropGen> void add_freq(BaseHistogram<complex>&, complex scale, PropGen&); //!< accumulate single spectrum
   template<class StaticObj,class SpinObj> void add_freq(BaseHistogram<complex>&, complex scale, const StaticObj& Hstaticp, const SpinObj& Hspinp, const ActionDirectAcq&); //!< accumulate single spectrum
   void add_freq(BaseHistogram<complex>&, complex scale, bool iscomplex, const ActionDirectAcq&); //!< accumulate single spectrum using real or complex Hamiltonian
   template<class StaticObj,class SpinObj> void add_td(BaseList<complex> FID, complex scale, const StaticObj& Hstaticp, const SpinObj& Hspinp, const ActionDirectAcq&); //!< accumulate single time-domain FID
   template<class PropGen> void add_td(BaseList<complex> FID, complex scale, PropGen&, double synctime =0);

   template<class T> T& ensurefobj(); //!< ensure spectrum object has been created
   template<class T> T& ensuretobj(); //!< ensure time domain object has been created
 
   template<class OpGen> void create_Hstatic(const OpGen&, Type2Type<double>); //!< create real static Hamiltonian
   template<class OpGen> void create_Hstatic(const OpGen&, Type2Type<complex>); //!< create complex static Hamiltonian
   template<class OpGen> void create_Hspin(const OpGen&, double, Type2Type<double>); //!< create real spinning Hamiltonian
   template<class OpGen> void create_Hspin(const OpGen&, double, Type2Type<complex>); //!< create complex spinning Hamiltonian
   template<typename T> static void set_Hsystem(const BlockedMatrix<T>&);

   template<class HType> void propagator_(BlockedMatrix<complex>&, const HType&, const CompSequenceBase&, double t1,double t2,double origin);
   template<class HType> void propagator_(BlockedMatrix<complex>&, const HType&, const Sequence&, double t1,double t2,double origin);
   template<class StaticObj,class SpinObj> void set_interactions(const Euler&, StaticObj& Hstaticp, SpinObj& Hspinp); //!< update Hamiltonian with NMR interactions for given orientation
   void getgammapars(GammaPeriodicFID&, double& use_period, size_t& steps, double) const;

   ham_t hamtype; //!< Hamiltonian type (full, diagonal etc.)
   hamobj_t hamobjtype; //!< Hamiltonian object type

   DataStore* curstorep; //!< temporary stores (not multi-thread compatible)
   //   double curscale; //!< scaling factor on current orientation/acqusition
};

typedef void (*prefinalise_callback_t)();
typedef void (*precalculation_callback_t)(MasterObj&);
typedef bool (*calculation_callback_t)(MasterObj&, DataStore&);

extern prefinalise_callback_t prefinalise_callback; //!< hook prior to finalise block
extern precalculation_callback_t precalculation_callback; //!< hook for pre-calculation optimisation
extern calculation_callback_t calculation_callback; //!< hook for pre-calculation optimisation
bool default_calculation(MasterObj&, DataStore&); //!< default calculation calculation

void initialise_simulation_environment(); //!< initialise objects / values prior to starting all calculations after reading in set up

#include "NMRsim_matrixspec.h"

//struct filter_def;
typedef std::pair<state_t, const filter_def*>  statistics_t;
std::pair<int,statistics_t> get_putmatrix_flags();

void dump_matrix(const matrixspec&, const char*, MasterObj* =NMRSIM_NULL, int flags =0, statistics_t = statistics_t(0, (const filter_def*)NMRSIM_NULL ) );

extern ThreadWarning<> freq_warning; //!< warn about frequency domain propagation
extern ThreadWarning<> changedsw_warning; //!< spectral width for row has changed
extern ThreadWarning<> nogammaspec_warning; //!< no gamma handling specified
extern ThreadWarning<> gammasetinstatic_warning; //!< gamma_angles set in static simulation
extern ThreadWarning<> arrayambiguity_warning; //!< ambiguity in interpretation of array lengths
extern ThreadWarning<> arrayincommensurate_warning;
extern ThreadWarning<> gammacompute_warning; //!< couldn't use gamma COMPUTE
extern ThreadWarning<> ignoringrev_warning; //!< ignoring rev flag for powder averaging
extern ThreadWarning<InvalidParameter> invalidweighting_warning; //!< found invalid powder weighting factor
extern ContextWarning<> ignoringrange_warning; //!< ignoring range qualifier on powder averaging
extern ThreadWarning<> orientationsfixed_warning; //!< number of powder orientations can't be changed
extern ThreadWarning<> gamma3angle_warning; //!< gamma angle integration combined with 3 angle set
extern ContextWarning<> rowsskip_warning; //!< number of rows in data not multiple of "skip"
#endif
