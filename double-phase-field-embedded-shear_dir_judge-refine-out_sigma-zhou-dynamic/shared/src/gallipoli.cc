#include "geocentric/materials/gallipoli.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
Gallipoli::Gallipoli()
{}

Gallipoli::~Gallipoli()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
void Gallipoli::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Gallipoli");
    prm.declare_entry("Parameter_a1","0.0",Patterns::Double());
    prm.declare_entry("Parameter_a2","0.0",Patterns::Double());
    prm.declare_entry("ExponentM","0.0",Patterns::Double());
    prm.declare_entry("ExponentN","0.0",Patterns::Double());
    prm.declare_entry("InitialPorosity","0.0",Patterns::Double());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
void Gallipoli::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("Gallipoli");
    a1     = prm.get_double("Parameter_a1");
    a2     = prm.get_double("Parameter_a2");
    exp_m  = prm.get_double("ExponentM");
    exp_n  = prm.get_double("ExponentN");
    v_zero = 1.0/(1.0 - prm.get_double("InitialPorosity"));
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// UPDATE
// -----------------------------------------------------------------------------
void Gallipoli::update(const double p, const double J, const bool sat_only)
{
  if (p >= 0)
  {
    S = 1.0;
    S_deriv_J = 0.0;
    S_deriv_p = 0.0;

    kr = 1.0;
    kr_deriv_J = 0.0;
    kr_deriv_p = 0.0;
  }
  else
  {
    double S_tmp = pow(-a1*pow(J*v_zero - 1, a2)*p, exp_n);
    S = pow(1 + S_tmp, -exp_m);

    if (sat_only)
      return;

    double S_deriv_S_tmp = -exp_m*pow(1 + S_tmp, -exp_m-1);

    S_deriv_J = S_deriv_S_tmp * ( a2*exp_n*v_zero*S_tmp/(J*v_zero - 1) );
    S_deriv_p = S_deriv_S_tmp * ( exp_n*S_tmp/p );

    double kr_tmp_1 = 1 - pow(S, 1/exp_m);
    double kr_tmp_1_deriv_S = -(1/exp_m)*pow(S, 1/exp_m - 1);
    double kr_tmp_2 = 1 - pow(kr_tmp_1,exp_m);
    double kr_tmp_2_deriv_kr_tmp_1 = -exp_m*pow(kr_tmp_1, exp_m - 1);

    kr = sqrt(S)*pow(kr_tmp_2, 2);

    kr_deriv_S = 0.5*pow(S,-0.5)*pow(kr_tmp_2, 2)
                 + 2*sqrt(S)*kr_tmp_2*kr_tmp_2_deriv_kr_tmp_1*kr_tmp_1_deriv_S;

    kr_deriv_J = kr_deriv_S*S_deriv_J;
    kr_deriv_p = kr_deriv_S*S_deriv_p;
  }
}


// -----------------------------------------------------------------------------
// SATURATION/PERMEABILITY FOR GIVEN PRESSURE/JACOABIAN
// -----------------------------------------------------------------------------
double Gallipoli::sat(const double p, const double J)
{
  if (p >= 0)
  {
    return 1.0;
  }
  else
  {
    return pow(1 + pow(-a1*pow(J*v_zero - 1, a2)*p, exp_n), -exp_m);
  }
}
double Gallipoli::sat_deriv_p(const double p, const double J)
{
  if (p >= 0)
  {
    return 0.0;
  }
  else
  {
    double S_tmp = pow(-a1*pow(J*v_zero - 1, a2)*p, exp_n);
    double S_deriv_S_tmp = -exp_m*pow(1 + S_tmp, -exp_m-1);

    return S_deriv_S_tmp * ( exp_n*S_tmp/p );
  }
}
double Gallipoli::rel_perm(const double p, const double J)
{
  if (p >= 0)
  {
    return 1.0;
  }
  {
    double this_S = sat(p,J);

    double kr_tmp_1 = 1 - pow(this_S, 1/exp_m);
    double kr_tmp_2 = 1 - pow(kr_tmp_1,exp_m);

    return sqrt(this_S)*pow(kr_tmp_2, 2);
  }
}
double Gallipoli::rel_perm_deriv_p(const double p, const double J)
{
  if (p >= 0)
  {
    return 0.0;
  }
  {
    double this_S = sat(p,J);

    double kr_tmp_1 = 1 - pow(this_S, 1/exp_m);
    double kr_tmp_1_deriv_S = -(1/exp_m)*pow(this_S, 1/exp_m - 1);
    double kr_tmp_2 = 1 - pow(kr_tmp_1,exp_m);
    double kr_tmp_2_deriv_kr_tmp_1 = -exp_m*pow(kr_tmp_1, exp_m - 1);

    double this_kr_deriv_S = 0.5*pow(this_S,-0.5)*pow(kr_tmp_2, 2)
                             + 2*sqrt(this_S)*kr_tmp_2*kr_tmp_2_deriv_kr_tmp_1*kr_tmp_1_deriv_S;

    return this_kr_deriv_S*sat_deriv_p(p,J);
  }
}


// -----------------------------------------------------------------------------
// SET INITIAL SOLID VOLUME
// -----------------------------------------------------------------------------
void Gallipoli::set_init_porosity(const double new_phi_f_zero)
{
  v_zero = 1.0/(1.0 - new_phi_f_zero);
}


// -----------------------------------------------------------------------------
// "GET" METHODS
// -----------------------------------------------------------------------------
double Gallipoli::sat(){ return S; }
double Gallipoli::sat_deriv_J(){ return S_deriv_J; }
double Gallipoli::sat_deriv_p(){ return S_deriv_p; }
double Gallipoli::rel_perm(){ return kr; }
double Gallipoli::rel_perm_deriv_J(){ return kr_deriv_J; }
double Gallipoli::rel_perm_deriv_p(){ return kr_deriv_p; }
