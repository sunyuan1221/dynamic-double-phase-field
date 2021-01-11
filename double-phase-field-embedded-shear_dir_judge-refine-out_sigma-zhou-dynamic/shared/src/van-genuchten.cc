#include "geocentric/materials/van-genuchten.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
VanGenuchten::VanGenuchten()
{}

VanGenuchten::~VanGenuchten()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
void VanGenuchten::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VanGenuchten");
    prm.declare_entry("MaximumSat","1.0",Patterns::Double());
    prm.declare_entry("MinimumSat","0.0",Patterns::Double());
    prm.declare_entry("ScalingPressure","0.0",Patterns::Double());
    prm.declare_entry("ExponentN","0.0",Patterns::Double());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
void VanGenuchten::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("VanGenuchten");
    Ss = prm.get_double("MaximumSat");
    Sr = prm.get_double("MinimumSat");
    sa = prm.get_double("ScalingPressure");
    exp_n = prm.get_double("ExponentN");
    exp_m = (exp_n-1)/exp_n;
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// SATURATION AND ITS DERIVATIVES
// -----------------------------------------------------------------------------
void VanGenuchten::update(const double p, const bool sat_only)
{
  if (p >= 0)
  {
    S = Ss;
    S_deriv_p = 0.0;
    kr = 1.0;
    kr_deriv_p = 0.0;
  }
  else
  {
    S = Sr + (Ss - Sr)*pow(1 + pow(-p/sa,exp_n),-exp_m);

    if (sat_only)
      return;

    double ratio = -p/sa;
    double rn = pow(ratio,exp_n);
    S_deriv_p = (Ss - Sr)/sa*exp_m*exp_n*pow(1 + pow(ratio,exp_n),-exp_m-1)*pow(ratio,exp_n-1);
    kr = pow(1 - pow(ratio,exp_n-1)*pow(1 + pow(ratio,exp_n),-exp_m),2)
         * pow(1 + pow(ratio,exp_n),-exp_m/2);
    kr_deriv_p = 0.5 / sa * pow(ratio,exp_n-3) * pow(1 + rn,-5*exp_m/2-1)
                 * (rn - ratio*pow(1+rn,exp_m))
                 * (4*(1 + rn) + exp_n*(-4 + (5*exp_m - 4)*rn - exp_m*ratio*pow(1 + rn,exp_m)));
  }
}


// -----------------------------------------------------------------------------
// SATURATION AND RELATIVE PERMEABILITY (DIRECTLY FROM PRESSURE)
// -----------------------------------------------------------------------------
double VanGenuchten::sat(double p)
{
  if (p >= 0)
    return Ss;
  else
    return Sr + (Ss-Sr)*pow(1+pow(-p/sa,exp_n),-exp_m);
}


double VanGenuchten::sat_deriv_p(double p)
{
  if (p >= 0)
    return 0.0;
  else
  {
    double ratio = -p/sa;
    return (Ss-Sr)/sa*exp_m*exp_n*
           pow(1+pow(ratio,exp_n),-exp_m-1)*
           pow(ratio,exp_n-1);
  }
}


double VanGenuchten::rel_perm(double p)
{
  if (p >= 0)
    return 1.0;
  else
  {
    double ratio = -p/sa;
    return pow( 1-pow(ratio,exp_n-1)*pow(1+pow(ratio,exp_n),-exp_m) ,2)*
             pow( 1+pow(ratio,exp_n) ,-exp_m/2);
  }
}


double VanGenuchten::rel_perm_deriv_p(double p)
{
  if (p >= 0)
    return 0.0;
  else
  {
    double ratio = -p/sa;
    double rn = pow(ratio,exp_n);

    return 0.5 / sa * pow(ratio,exp_n-3) *
             pow(1+rn,-5*exp_m/2-1)*
             ( rn-ratio*pow(1+rn,exp_m))*
             ( 4*(1+rn)+exp_n *(-4+(5*exp_m-4)*rn-exp_m*ratio*pow(1+rn,exp_m)));
  }
}


// -----------------------------------------------------------------------------
// "GET" METHOD
// -----------------------------------------------------------------------------
double VanGenuchten::sat(){ return S; }
double VanGenuchten::sat_deriv_p(){ return S_deriv_p; }
double VanGenuchten::rel_perm(){ return kr; }
double VanGenuchten::rel_perm_deriv_p(){ return kr_deriv_p; }

