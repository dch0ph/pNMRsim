
// Shift spectrum / FID

class ProcessShift : public ProcessCommand {
  ProcessShift(double,double =0.0);
  static ProcessCommand* create();
  void exec(BaseList<complex>, procflags&) const;
  void exec(cmatrix&, procflags&) const;
  void set(double, subsid_t);
  void print(std::ostream&, subsid_t =S_NONE) const;
private:
  double shift_,shift1_;
  static void doshift1(cmatrix&, int); //!< integer roll in indirect dimension
  static void doshift(BaseList<complex>, int); //!< integer roll in direct dimension
  static size_t confine(size_t,size_t,int);
  static void doshift(BaseList<complex>&, const complex&); //!< apply arb roll to direct dimension 
  static int validateshift(double shf, double lsw, size_t lnp, size_t skip =1);
  static complex getshiftfactor(double shf, double lsw); //!< return complex multiplier factor
};

ProcessCommand* ProcessShift::create()
{
  const int flags=process2D ? F_DENYARRAY : 0;
  Variable cvar(S_ARG1);
  const double shiftv=parse_double(&cvar,flags);
  Variable cvar1(S_ARG2);
  const double shift1v=(process2D && are_left()) ? parse_double(&cvar1,flags) : 1.0;
  return new ProcessShift(shiftv,shift1v);
}

void ProcessShift::set(double v, subsid_t subsid)
{
  switch (subsid) {
  case S_ARG1:
    shift_=v;
    break;
  case S_ARG2:
    shift1_=v;
    break;
  default:
    throw InternalError("ProcessShift::set");
  }
}

void ProcessShift::print(std::ostream& ostr, subsid_t subsid) const
{
  ostr << "shift ";
  switch (subsid) {
  case S_NONE:
    ostr << shift_;
    if (shift1_)
      ostr << ' ' << shift1_;
    ostr << '\n';
    break;
  case S_ARG1:
    ostr << "shift=" << shift_;
    break;
  case S_ARG2:
    ostr << "shift1=" << shift1_;
    break;
  default:
    throw InvalidParameter("Processhift::print");
  }
}

ProcessShift::ProcessShift(double shiftv, double shift1v)
  : ProcessCommand(PROC_HASBOTH),
    shift_(shiftv), shift1_(shift1v) {}

size_t ProcessShift::confine(size_t i, size_t max, int roll)
{
  const int newi=int(i)+roll;
  if (newi<0)
    return size_t(newi+max);
  return (newi>=max) ? size_t(i-int(max)) : newi;
}

void ProcessShift::doshift(BaseList<complex> a, int roll)
{
  static LIST<complex> tmp; //!< N.B. Not re-entrant
  const size_t cs(a.size());
  tmp.create(cs);
  for (size_t i=cs;i--;)
    tmp(i)=a(confine(i,cs,roll));
  a=tmp;
}

void ProcessShift::doshift(BaseList<complex>& a, const complex& shiftfac)
{
  complex curfac(1.0);
  const size_t len=a.size();
  for (size_t i=1;i<len;i++) {
    curfac*=shiftfac;
    a(i)*=curfac;
  }
}
 
complex ProcessShift::getshiftfactor(double shf, double lsw) { return expi(TWOPI*shf/lsw); }

int ProcessShift::validateshift(double shf, double lsw, size_t lnp, size_t skip)
{
  if (lsw==0.0)
    throw InternalError("ProcessShift: spectral width not defined");
  if ((skip!=1) && (lnp % skip))
    error_abort("shift cannot be applied to indirect dimension - number of rows is not a multiple of the skip factor\n");
  const double df=skip*lsw/lnp;
  const double froll=shf/df;
  const int roll=rawround_int(froll);
  if (fabs(froll-roll)>NMRSIM_ROUNDTOL) {
    std::cerr << "shift in frequency domain (" << shf << " Hz) must be multiple of frequency resolution (" << df << " Hz)\n";
    error_abort();
  }
  return skip*roll;
}
    
void ProcessShift::exec(BaseList<complex> a, procflags& pflags) const
{
  if (shift_) {
    if (pflags.istimedomain)
      doshift(a,getshiftfactor(shift_,sw));
    else
      doshift(a,validateshift(shift_,sw,np));
  }
}

void ProcessShift::doshift1(cmatrix& a, int roll)
{
  const size_t rs=a.rows();
  cmatrix d(rs,a.cols());

  for (size_t i=rs;i--;)
    d.row(i)=a.row(confine(i,rs,roll));

  a.swap(d);
}

void ProcessShift::exec(cmatrix& a, procflags& pflags) const
{
  
  if (shift1_) {
    const size_t skip=skips.front();
    const double sw1=sws.front();
    if (pflags.istimedomain1) {
      const complex shiftfac1(getshiftfactor(shift1_,sw1));
      complex curshift1(1.0);
      for (size_t i=0;i<a.rows();) {
	BaseList<complex> ai(a.row(i));
	if (i)
	  ai*=curshift1;
	exec(ai,pflags);
	i++;
	if ((i % skip)==0)
	  curshift1*=shiftfac1;
      }
      return;
    }
    else
      doshift1(a,validateshift(shift_,sw1,a.rows(),skip));    
  }
  if (shift_) {
    for (size_t i=a.rows();i--;)
      exec(a.row(i),pflags);
  }
}

