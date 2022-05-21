#ifndef ref_const_value_actor_hpp_
#define ref_const_value_actor_hpp_

/*! \file 
 \brief  declares assign-to-const actor

 Adapted from boost/spirit/actor/ref_const_ref_actor.hpp
*/

template<
  typename T,
  typename ValueT,
  typename ActionT
  >
class ref_const_value_actor : public ActionT
{
private:
  T& ref; //!< reference to destination
  ValueT value; //!< value
public:
  ref_const_value_actor(
		  T& ref_,
		  ValueT const& value_
		  )
    :
    ref(ref_),
    value(value_)
  {}
    
  //!< do assign
  template<typename T2>
  void operator()(T2 const& ) const
  {
    this->act(ref,value);
  }

  //!< do assign  
  template<typename IteratorT>
  void operator()(
		  IteratorT const&,
		  IteratorT const&
		  ) const
  {
    this->act(ref,value);
  }
};

template<
  typename T,
  typename ValueT
  >
inline ref_const_value_actor<T,ValueT,assign_action> assign_const_a(
								  T& ref_,
								  ValueT const& value_
								  )
{
  return ref_const_value_actor<T,ValueT,assign_action>(ref_,value_);
}

#endif
