#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
using boost::math::cyl_bessel_j;
using boost::math::ellint_1;
using boost::math::ellint_2;

namespace dolfin {
  class phi : public Expression
  {
    
  public :  
    phi () : Expression () {}
    
    void eval (Array <double >& values , const Array <double >& x) const
    {
      values [0] = ellint_1 ( x[0] ) ;
    }
  };

}
