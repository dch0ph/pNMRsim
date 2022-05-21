#ifndef LCM_rational_h_
#define LCM_rational_h_

/*! \file
  \brief  Defines simple class for working with rational numbers */

#include <iostream>
#include <cmath>

namespace libcmatrix {

  //! class to store rational or non-rational number 
  class rational {
  public:
    rational(double rv) //!< initialise with non-rational value
      : M(0), N(0), ratio(rv) {}
    rational(int Mv,int Nv) //!< initialise rational number \a Mv / \a Nv
      : M(Mv), N(Nv), ratio(double(M)/N) {}
    operator double() const { return ratio; } //!< return as floating point
    bool operator! () const { return (N==0); } //!< \c true if not storing a rational
    rational reciprocal() const { //!< return reciprocal
      return N ? rational(N,M) : rational(1/ratio);
    }
    int numerator() const { return M; } //!< getter for numerator (0 if non-rational)
    int denominator() const { return N; } //!< getter for denominator (0 if non-rational)

    //! \c true if equal to \a a (within tolerance)
    /*! \note must compare floating point representation otherwise 3/9 will not equal 1/3 etc.
     */
    bool operator== (const rational& a) const {
      return (fabs(a.ratio-ratio)<=tolerance);
    }
    bool operator!= (const rational& a) const { //!< \c true if distinct from \a a (within tolerance)
      return (fabs(a.ratio-ratio)>tolerance);
    }
 
    friend std::ostream& operator<< (std::ostream&, const rational&);

  private:
    int M; //!< numerator (0 if non-rational) 
    int N; //!< denominator (0 if non-rational)
    double ratio; //!< as floating point
    static const double tolerance;  //!< tolerance within which two ratios are regarded as equal [Not defined in class - candidate for constexpr]
  };

  inline std::ostream& operator<< (std::ostream& ostr, const rational& a)
  {
    return ostr << a.M << '/' << a.N << " (" << a.ratio << ')';
  }

}

#endif

