// Components for looping over arrays

#include "NMRsim.h"
#include "Parser.h"

static size_t array_loop=0;
bool update_vars=true;

size_t ndims=-1;
size_t nacqdims=-1;
LIST<size_t> sum_ns,array_ns;
size_t sum_n0=-1;
size_t array_n0=-1;
//size_t nD_n0=-1;
LIST<size_t> array_skips,arrayinds;
bool active2D=false; //!< true if nD sequence

LIST<size_t> skips(MAX_DIMENSIONS,1U); //!< only used for reporting (better merged)
LIST<int> ns(MAX_DIMENSIONS,0); //!< only used for reporting (better merged)
LIST<double> sws(MAX_DIMENSIONS,0.0);
//LIST<double> sfrqs(MAX_DIMENSIONS,0.0);

SystemVariable<int*> v_ni("ni",&(ns.front()),V_ISFIXED | V_ISCONST);
SystemVariable<double*> v_sw1("sw1",&(sws.front()),1.0,V_ISFIXED | V_ISCONST);

int row_index=-1;
int var_index=-1;
int sum_index=-1;

void ensure_array()
{
  if (array_ns.empty() && (array_n0>1))
    array_ns.create(1,array_n0);
  array_loop=array_n0 ? array_n0 : 1;
}

//!< almighty fudge - changing data set size doesn't sit naturally with looping etc.
void update_array_loop(size_t nr)
{
  if ( (array_ns.size()>1) || have_virtualdimensions)
    error_abort("Processing which changes number of rows is currently incompatible with multi-dimensional data and/or virtual dimensions");
  array_loop=array_n0=nr; 
  if (array_ns.size()==1)
    array_ns.front()=nr;
}

size_t array_iter::size() const
{
  if (!array_loop)
    throw InternalError("array_iter::size");
  return array_loop;
}

array_iter::array_iter(bool updatevarsv)
  : updatevars_(updatevarsv), iter(array_ns)
{ 
  var_index=-1;
  row_index=-1;
  finished=false;
  if (array_loop==0)
    throw InternalError("array_iter: array_loop not set");
  //  if (updatevars_)
  //  varinds.create(ndims-1);
}

ThreadWarning<> arrayambiguity_warning("list length is ambiguous; values will be cycled with increments",&NMRsim_once_warning);
//Warning<> arrayincommensurate_warning("variable array length incommensurate with data array length",&NMRsim_once_warning);

static void gettranslation(size_t& skipfac,size_t& modfac,size_t dim,size_t arraylen)
{
  if (dim) { //virtual dimension
    skipfac=array_skips(dim-1);
    modfac=array_ns(dim-1);
  }
  else {
    skipfac=active2D ? array_skips.front() : 1; // if no explicit tag use indirect dimension skip if nD
    modfac=array_n0;
  }
  //const bool suppress=!allowwarnings();
  char buf[256];
  if (skipfac>1) {
    const size_t actrows=(modfac+skipfac-1)/skipfac; //!< rather wasteful re-evaluating this each time...
    if ((arraylen==skipfac) && (arraylen==actrows)) {
      std::cerr << "Ambiguity in dimension " << dim << ": " << "array of length " << arraylen << " could cycle with quadrature factor or number of increments - resolve by specifying values for each row\n";
      error_abort();
    }
    if (arraylen==skipfac) {
      modfac=skipfac;
      skipfac=1;
    }
    else {
      if (arraylen==modfac)
	skipfac=1;
      else {
	if (arraylen==actrows)
	  modfac=actrows;
	else {
	  snprintf(buf,sizeof(buf)," (array length=%lu)",(unsigned long)arraylen);	  
	  arrayambiguity_warning.raise(buf);
	  skipfac=1;
	  return; //!< not sure about this...
	}
      }
    }	
  }
  // if (suppress)
  //  return;
  bool incom=false;
  //  const size_t totfac=modfac*skipfac;
  if (arraylen<modfac)
    incom=(modfac % arraylen);
  else
    incom=(arraylen % modfac);
  if (incom) {
    std::cerr << "Variable array length (" << arraylen << ") doesn't match data array size (" << modfac << ").  Lists of incompatible length (check for earlier warnings)?\n";
    error_abort();
    //    arrayincommensurate_warning.raise(buf,true);
  }
  if (arraylen<modfac)
    modfac=arraylen;
    //    modfac=skipfac*arraylen; //!< no wrap-around
}

bool getindex(size_t& transindex, size_t tag, size_t arraylen)
{
  if (arraylen==1) {
    transindex=0;
    return false;
  }
  size_t usei;
  if (tag) { 
	if (arrayinds.empty())
		throw InternalError("getindex: called before row iterator set up");
    usei=arrayinds(tag-1);
  }
  else {
    if (row_index<0) //!< change from var_index 17/10/13 - var index not always updated
      throw InternalError("getindex: row_index not valid");
    usei=row_index;
  }
  size_t skipfac,modfac;
  gettranslation(skipfac,modfac,tag,arraylen);
  bool update=true;
  if (skipfac>1) {
    if (usei % skipfac)
      update=false;
    usei/=skipfac;
  }
  transindex=usei % modfac;
  if ((verbose & VER_GEN) && (verbose_level>1)) {
    if (update)
      std::cout << "Translated global row index " << usei << " into parameter index " << transindex << " (actual index: " << usei << "  modulation factor: " << modfac << "  skip factor: " << skipfac << ")\n";
    else
      std::cout << "No update for current row index " << usei << '\n';
  }
  return update;
}

size_t VariableBase::summationindex() const
{
  if (!isarray())
    throw InternalError("VariableBase::summationindex");
  size_t which=0;
  if (array.size()>1) {
    if (sum_index<0)
      //      error_abort("attempt to use sum array || outside valid region (processing command belongs in initialproc?)");
      error_abort("attempt to use sum array || outside valid region");
    which=(sumtag_ ? suminds(sumtag_-1) : sum_index) % array.size();
  }
  return which;
}

// const BaseList<double> VariableBase::getsummation() const
// {
//   assert(isarray());
//   size_t which=0;
//   if (array.size()>1) {
//     if (sum_index<0)
//       throw Failed("sum_index used outside valid region");
//     which=(sumtag_ ? suminds(sumtag_-1) : sum_index) % array.size();
//   }
//   return array(which);
// }

std::pair<size_t,bool> VariableBase::getoffset() const
{
  const size_t sblock=summationindex();
  const size_t bsize=array.size(sblock);
  if (bsize==arraylength_)
    return std::pair<size_t,bool>(array.offset(sblock),true);
  size_t transindex;
  const size_t tag = arraytags_.empty() ? 0 : arraytags_(sblock);
  const bool needup=getindex(transindex,tag,bsize);
  return std::pair<size_t,bool>(array.offset(sblock,transindex),needup);
}

const BaseList<double> VariableBase::get_row(size_t i) const
{
  return BaseList<double>(arraylength_, const_cast<double*> (&(array.row()(i*arraylength_))));
}

BaseList<double> VariableBase::get_row(size_t i)
{
  if (!nochecks) {
    const size_t offset=(i+1)*arraylength_;
    if (offset>array.items())
      throw BadIndex("VariableBase::get_row",offset,array.items());
  }
  return BaseList<double>(arraylength_, &(array.row()(i*arraylength_)));
}

void VariableBase::set_current_row(const BaseList<double>& vs)
{
  if (!isarray())
    throw InternalError("VariableBase::set_current_row");
  std::pair<size_t,bool> res(getoffset());
  BaseList<double> currow(get_row(res.first));
  currow=vs;
  value_=vs; //!< make sure current value is also updated
}

bool VariableBase::updateindex()
{
  static LIST<double> tmp;
  if (!isarray())
    throw InternalError("VariableBase::updateindex");
  if ((array.size()>1) && (sum_index<0))
    return false; //!< don't update summation arrays outside region where valid
  std::pair<size_t,bool> res(getoffset());
//   if ((verbose & VER_GEN) & (verbose_level>1))
//     std::cout << "Index: " << res.first << "  Update: " << (res.second ? "Yes\n" : "No\n");
  if (res.second) {
    if (arraylength_==1) {
      double& val(array.row()(res.first));
      if (!(arrayexprs_.empty())) {
	const Expression* exprp=arrayexprs_(res.first);
	if (exprp) {
	  exprp->get(tmp);
	  const size_t ind=exproffsets_(res.first);
	  if (ind>=tmp.size()) {
	    std::cerr << "Expression in array has returned fewer elements (" << tmp.size() << ") than required (" << (ind+1) << ")\n";
	    error_abort();
	  }
	  val=tmp(ind);
	}
      }
      value_.create(1,val);
    }
    else {
      if (!(arrayexprs_.empty()))
	throw InternalError("can't have expression array with array_length!=1");
      value_=get_row(res.first);
    }
    return true;
  }
  return false;
}
    
size_t getarrayindex(size_t i)
{
  size_t res=arrayinds(i);
  if (active2D)
    res/=array_skips(i);
  return res;
}
      
bool array_iter::next(const DataStore* datap)
{
  if (finished)
    return false;

  row_index++;
  if (datap)
    current_data_row.create(datap->row(row_index));
  else
    current_data_row.clear(); //!< can't be valid

  //  const bool forceupdate=(optupdate==0) || !active2D;
  update_vars=(row_index==0) || !active2D || !optupdate;
  if (updatevars_) {
      
    NMRSIM_EXPECT(iter.next(arrayinds));
    
//     if (active2D) {
//       for (size_t i=ndims-1;i--;) {
// 	// 	if ((arrayinds(i) % array_skips(i))==0) {
// 	varinds(i)=arrayinds(i)/array_skips(i);
//  	  //update_vars=true;
// 	  //}
//       }
//     }
//     else
//       varinds=arrayinds;

//    bool updateexpr=forceupdate;
    var_index++;
    //if ((verbose & VER_GEN) && (verbose_level>1))
    //  std::cout << "Parameter index: " << var_index << "  indices: " << varinds << '\n';
    //    update_variables_expressions(varpars,var_index,varinds);

    const varpars_t::iterator end=varpars.end();
    varpars_t::iterator start=varpars.begin();
    while (start!=end) {
      if (start->variable().updateindex()) {
	start->rawupdate();
	update_vars=true; //!< flag that at least one variable has changed
      }
      ++start;
    }
    //    update_variables(varpars,var_index,varinds);
    // "Updating expressions" redundant if followed by updates
    if ((verbose & VER_GEN) && (!update_vars || (verbose_level>1)))
      std::cout << "Updating expressions: " << (update_vars ? "Yes\n" : "No\n");
    if (update_vars)
      update_expressions();
  }
  if (update_vars)
    update_auxiliary_vars();

  finished=(row_index==array_loop-1);
  return true;
}

bool accumulator::add(size_t n)
{
  if ((length_==1) || (length_<0))
    length_=n;
  else {
    if ((n!=1) && (n!=length_)) {
      failed_=false;
      return false;
    }
  }
  return true;
}

size_t accumulator::operator()() const
{ 
  if (length_<0)
    throw Undefined("accumulator");
  if (failed_)
    throw Failed("accumulator: incompatible list sizes");
  return length_;
}

void ensure_dimension(size_t i)
{
  ScratchList<bool,MAX_DIMENSIONS> done(MAX_DIMENSIONS,false);
  if (done(i))
    return;
  char vname[6];
  sprintf(vname,"n%lu",(unsigned long)i);
  add_systemvarmap(*(new SystemVariable<int*>(vname,&(ns(i-1)),V_ISFIXED | V_ISCONST)));
  sprintf(vname,"sw%lu",(unsigned long)i);
  add_systemvarmap(*(new SystemVariable<double*>(vname,&(sws(i-1)),1.0,V_ISFIXED | V_ISCONST)));
}

void raw_set_n(int dim, int n, int ni_skip)
{
  ns(dim-1)=n; //!< only used for reporting
  skips(dim-1)=ni_skip;
}

void set_n(int dim, int n, int ni_skip)
{
  ensure_dimension(dim);
  array_dims.set(dim,n,ni_skip);
  raw_set_n(dim,n,ni_skip);
  active2D=true;
}

void parse_swn(int n)
{
  ensure_dimension(n);
  double& cursw(sws(n-1));
  if (cursw)
    error_abort("cannot change indirect dimension spectral width");  
  cursw=parse_double();
  if (cursw<=0.0)
    error_abort("spectral width cannot be <=0!");
//   if (are_left())
//     sfrqs(n-1)=parse_double();
  if (n==1)
    v_sw1.isconstant(true); //!< flag that sw1 is now fixed
}
