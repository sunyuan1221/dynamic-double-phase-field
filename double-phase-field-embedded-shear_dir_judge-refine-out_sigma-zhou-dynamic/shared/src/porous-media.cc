#include "geocentric/materials/porous-media.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
PorousMedia::PorousMedia()
{}

PorousMedia::~PorousMedia()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
void PorousMedia::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("BasicProperties");
    prm.declare_entry("SolidDensity","0.0",Patterns::Double());
    prm.declare_entry("FluidDensity","0.0",Patterns::Double());
    prm.declare_entry("BiotCoefficient","0.0",Patterns::Double());
    prm.declare_entry("Porosity","0.0",Patterns::Double());
    prm.declare_entry("AbsolutePerm","0.0",Patterns::Double());
    prm.declare_entry("DynamicViscosity","0.0",Patterns::Double());
    prm.declare_entry("ConstantPermeability","false",Patterns::Bool());
    prm.declare_entry("FluidBulkModulus","0.0",Patterns::Double()); 
    // thermal
    prm.declare_entry("ThermalConductivity","0.0",Patterns::Double());
    prm.declare_entry("SpecificHeatCapacity","0.0",Patterns::Double());
    prm.declare_entry("ThermalExpansionCoeff","0.0",Patterns::Double());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
void PorousMedia::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("BasicProperties");
    rho_s = prm.get_double("SolidDensity");
    rho_f = prm.get_double("FluidDensity");
    biot  = prm.get_double("BiotCoefficient");
    phi   = prm.get_double("Porosity");
    abs_perm = prm.get_double("AbsolutePerm"); 
    viscosity = prm.get_double("DynamicViscosity");
    const_perm = prm.get_bool("ConstantPermeability");
    fluid_bulk = prm.get_double("FluidBulkModulus"); 
    // thermal
    kappa_t = prm.get_double("ThermalConductivity");
    c_t     = prm.get_double("SpecificHeatCapacity");
    alpha_t = prm.get_double("ThermalExpansionCoeff");
  prm.leave_subsection();

  kappa = abs_perm/viscosity; 

  // Mixture density (constant under infinitesimal deformation/full saturation)
  rho_bar = phi*rho_f + (1 - phi)*rho_s;

  // For finite strain
  phi_s_zero = 1 - phi;
  kappa_zero_KC = kappa * pow(1 - phi, 2) / pow(phi, 3); // phi_f_zero = 1 - phi_s_zero
}


// -----------------------------------------------------------------------------
// MOBILITY AND DERIVATIVE (DIRECTLY FROM JACOBIAN)
// -----------------------------------------------------------------------------
double PorousMedia::porosity(const double J)
{
  return 1 - phi_s_zero/J;
}

double PorousMedia::porosity_deriv_J(const double J)
{
  return phi_s_zero/(J*J);
}

double PorousMedia::mobility(const double J)
{
  if (const_perm)
  {
    return kappa;
  }
  else
  {
    double phi_tmp = 1 - phi_s_zero/J;
    return kappa_zero_KC * pow(phi_tmp, 3) / pow(1 - phi_tmp, 2);
  }
}

double PorousMedia::mobility_deriv_J(const double J)
{
  if (const_perm)
  {
    return 0;
  }
  else
  {
    double phi_tmp = 1 - phi_s_zero/J;

    double kappa_deriv_phi = kappa_zero_KC * (3 - phi_tmp) * pow(phi_tmp, 2) / pow(1 - phi_tmp, 3);
    double phi_deriv_J     = phi_s_zero/(J*J);

    return kappa_deriv_phi * phi_deriv_J;
  }
}

double PorousMedia::mixture_density(const double J)
{
  return porosity(J)*rho_f + (1 - porosity(J))*rho_s;
}

double PorousMedia::mixture_density_deriv_J(const double J)
{
  return porosity_deriv_J(J)*(rho_f - rho_s);
}


// -----------------------------------------------------------------------------
// GET PHYSICAL PROPERTIES
// -----------------------------------------------------------------------------
double PorousMedia::solid_density() { return rho_s; }
double PorousMedia::fluid_density() { return rho_f; }
double PorousMedia::biot_coefficient() { return biot; }
double PorousMedia::porosity() { return phi; }
double PorousMedia::mobility() { return kappa; }
double PorousMedia::mixture_density() { return rho_bar; }
double PorousMedia::abs_perm_output(){ return abs_perm; }
double PorousMedia::viscosity_output(){ return viscosity; }
double PorousMedia::fluid_modulus(){ return fluid_bulk; }
// thermal
double PorousMedia::thermal_conductivity() { return kappa_t; }
double PorousMedia::specific_heat_capacity() { return c_t; }
double PorousMedia::thermal_expansion_coefficient() { return alpha_t; }
