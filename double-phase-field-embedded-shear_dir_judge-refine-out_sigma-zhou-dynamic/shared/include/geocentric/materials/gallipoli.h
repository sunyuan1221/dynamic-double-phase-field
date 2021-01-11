#ifndef GALLIPOLI_H
#define GALLIPOLI_H

#include <deal.II/base/parameter_handler.h>

using namespace dealii;


class Gallipoli
{
  public:
    Gallipoli();
    ~Gallipoli();

    void declare_parameters(ParameterHandler &prm);
    void initialize(ParameterHandler &prm);

    void update(const double p, const double J, const bool sat_only=false);

    void set_init_porosity(const double new_phi_f_zero);

    double sat();
    double sat(const double p, const double J);
    double sat_deriv_J();
    double sat_deriv_p();
    double sat_deriv_p(const double p, const double J);
    double rel_perm();
    double rel_perm(const double p, const double J);
    double rel_perm_deriv_J();
    double rel_perm_deriv_p();
    double rel_perm_deriv_p(const double p, const double J);

  private:
    double a1;
    double a2;
    double exp_m;
    double exp_n;
    double v_zero;

    double S;
    double S_deriv_J;
    double S_deriv_p;
    double kr;
    double kr_deriv_J;
    double kr_deriv_p;
    double kr_deriv_S;
};


#endif
