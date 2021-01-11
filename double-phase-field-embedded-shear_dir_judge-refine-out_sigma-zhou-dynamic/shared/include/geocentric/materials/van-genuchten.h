#ifndef VAN_GENUCHTEN_H
#define VAN_GENUCHTEN_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


class VanGenuchten
{
  public:
    VanGenuchten();
    ~VanGenuchten();

    void declare_parameters(ParameterHandler &prm);
    void initialize(ParameterHandler &prm);

    void update(const double p, const bool sat_only = false);

    double sat();
    double sat(const double p);
    double sat_deriv_p();
    double sat_deriv_p(const double p);
    double rel_perm();
    double rel_perm(const double p);
    double rel_perm_deriv_p();
    double rel_perm_deriv_p(const double p);

  private:
    double Sr;
    double Ss;
    double sa;
    double exp_n;
    double exp_m;

    double S;
    double S_deriv_p;
    double kr;
    double kr_deriv_p;
};


#endif
