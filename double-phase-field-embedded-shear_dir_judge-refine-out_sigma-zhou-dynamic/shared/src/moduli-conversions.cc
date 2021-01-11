#include <geocentric/utilities/moduli-conversions.h>


namespace convert_moduli
{

  // If bulk modulus and Poisson's ratio are known
  namespace bulk_poisson
  {
    double mu(const double bulk, const double poisson)
    {
      return 3*bulk*(1-2*poisson)/(2+2*poisson);
    }

    double lambda(const double bulk, const double poisson)
    {
      return 3*bulk*poisson/(1+poisson);
    }

    double youngs(const double bulk, const double poisson)
    {
      return 3*bulk*(1-2*poisson);
    }
  }

  // If Young's modulus and Poisson's ratio are known
  namespace youngs_poisson
  {
    double mu(const double youngs, const double poisson)
    {
      return youngs/(2+2*poisson);
    }

    double lambda(const double youngs, const double poisson)
    {
      return youngs*poisson/(1+poisson)/(1-2*poisson);
    }
    double bulk(const double youngs, const double poisson)
    {
      return youngs/(3*(1-2*poisson)); 
    }
  }

  // If shear modulus and Poisson's ratio are known
  namespace shear_poisson
  {
    double lambda(const double mu, const double poisson)
    {
      return (2*mu*poisson)/(1-2*poisson);
    }
  }

}
