#ifndef Interaction_h_
#define Interaction_h_

/*! \file
 \brief Defines NMR interactions for \c spinsys 
*/

#include "NMRsim_spinsys.h"

//! description of interaction
struct interaction_info {
  interaction_info(const char* namev, interaction_t id_, bool iscoupling_, ordering_convention_t ordering_ =convention_Haeberlen) 
    : name(namev), id(id_), iscoupling(iscoupling_), ordering(ordering_) {}

  std::string name; //!< name (pNMRsim)
  interaction_t id; //!< interaction identifier (libcmatrix)
  bool iscoupling; //!< \c true if coupling (rather than shift)
  ordering_convention_t ordering;
};

//! NMR interaction object
class Interaction : public Setable {
public:
  //initialise anisotropic interaction
  Interaction(const interaction_info& info_,size_t ni_,size_t nj_,double iso_,double aniso_,double xy_,bool isxy, const Euler&, ordering_convention_t, subsid_t =S_NONE, double scale =1.0);
  //initialise isotropic interaction
  Interaction(const interaction_info& info_,size_t ni_,size_t nj_,double iso_, subsid_t =S_NONE);

  void print(std::ostream&) const; //!< stream contents
  void set(double fval,subsid_t subsid); //!< set parameter \a subsid to \a fval
  void set(const BaseList<double>&, subsid_t); //!< whole tensor set
  void printvariablename(std::ostream&, subsid_t) const; //!< create variable name
  void ensurevalid();
  void update_tensor(); //!< recreate tensor 
  bool isset() const; //!< \c true if tensor has been set already
  bool isxy() const { return isxy_; } //!< \c true if asymmetry is being handled as xx-yy

 private:
  const interaction_info& info; //!< reference to interaction type
  const size_t ni; //!< spin number
  const size_t nj; //!< 'other' spin (if coupling)
  double iso; //!< isotropic component
  double aniso; //!< anisotropy
  double asymspec; //!< asymmetry/xx-yy
  Euler PAS; //!< Euler angles defining PAS
  const ordering_convention_t ordering; //!< tensor ordering convention
  const double aniso_scale; //!< scale factor between external representation of anisotropy and frequency
  bool isdirty_; //!< \c true if parameters changed but tensor not updated
  bool isxy_; //!< \c true if interaction is defined in terms of xx-yy rather than eta
  subsid_t allowtensor_; //!< allow full tensor definition if this subsid is being changed (S_NONE for none)
  // double gamma() const; //!< return gamma (for first spin)

  bool set_tensor(); //!< after initial set of tensor
  void update_tensor(const space_T&); //!< change tensor
};  

std::ostream& operator<< (std::ostream&, const Interaction&);

extern ContextWarning<> protonfrequency_warning; //!< proton Larmor frequency seems small
extern ContextWarning<> nonzeroasymmetry_warning; //!< non-zero assymetry with zero anisotropy
extern ContextWarning<> dipolesign_warning; //!< dipole coupling appears to have incorrect sign

#endif
