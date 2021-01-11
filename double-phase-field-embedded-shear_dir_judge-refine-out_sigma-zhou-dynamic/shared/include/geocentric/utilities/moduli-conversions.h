#ifndef MODULI_CONVERSIONS_H
#define MODULI_CONVERSIONS_H


namespace convert_moduli
{

  namespace bulk_poisson
  {
    double mu(const double bulk, const double poisson);
    double lambda(const double bulk, const double poisson);
    double youngs(const double bulk, const double poisson);
  }

  namespace youngs_poisson
  {
    double mu(const double youngs, const double poisson);
    double lambda(const double youngs, const double poisson);
    double bulk(const double youngs, const double poisson); 
  }

  // If shear modulus and Poisson's ratio are known
  namespace shear_poisson
  {
    double lambda(const double mu, const double poisson);
  }

}


#endif
