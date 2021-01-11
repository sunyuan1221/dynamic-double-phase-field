#ifndef POROUS_MEDIA_H
#define POROUS_MEDIA_H

#include <deal.II/base/parameter_handler.h>

using namespace dealii;


class PorousMedia
{
  public:
    PorousMedia();
    ~PorousMedia();

    virtual void declare_parameters(ParameterHandler &prm);
    virtual void initialize(ParameterHandler &prm);

    double solid_density();
    double fluid_density();
    double biot_coefficient();
    double porosity();
    double porosity(const double J);
    double porosity_deriv_J(const double J);
    double mobility();
    double mobility(const double J);
    double mobility_deriv_J(const double J);
    double mixture_density();
    double mixture_density(const double J);
    double mixture_density_deriv_J(const double J);
    double abs_perm_output(); 
    double viscosity_output(); 
    double fluid_modulus(); 

    // thermal
    double thermal_conductivity();
    double specific_heat_capacity();
    double thermal_expansion_coefficient();

  private:
    double rho_s;
    double rho_f;
    double biot;
    double phi;
    double phi_s_zero;
    double kappa;
    double rho_bar;
    double abs_perm; 
    double viscosity; 
    double fluid_bulk; 

    bool   const_perm;
    double kappa_zero_KC;

    // thermal
    double kappa_t;
    double c_t;
    double alpha_t;
};


#endif
