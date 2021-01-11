#include "geocentric/materials/phasefield-shear-monolithic.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
PhaseFieldShearMono<dim>::PhaseFieldShearMono()
{ this->name_ = "Phase-field Shear Monolithic"; }

template <int dim>
PhaseFieldShearMono<dim>::~PhaseFieldShearMono()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldShearMono<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldShearMono");
    prm.declare_entry("BulkMod","0.0",Patterns::Double());
    // prm.declare_entry("YoungsMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergy","0.0",Patterns::Double());
    prm.declare_entry("LengthScale","0.0",Patterns::Double());
    prm.declare_entry("CriticalStress","0.0",Patterns::Double());
    prm.declare_entry("ResidualStress","0.0",Patterns::Double()); 
    prm.declare_entry("Xi","0.0",Patterns::Double()); 
    prm.declare_entry("ShapePara","0.0",Patterns::Double()); 
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldShearMono<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldShearMono");
    K      = prm.get_double("BulkMod");
    // E      = prm.get_double("YoungsMod");
    nu     = prm.get_double("Poisson");
    G_c    = prm.get_double("FractureEnergy");
    l_zero = prm.get_double("LengthScale");
    sig_c  = prm.get_double("CriticalStress");
    sig_r  = prm.get_double("ResidualStress"); 
    xi     = prm.get_double("Xi"); 
    a2     = prm.get_double("ShapePara"); 
  prm.leave_subsection();

  mu     = convert_moduli::bulk_poisson::mu     (K,nu);
  lambda = convert_moduli::bulk_poisson::lambda (K,nu);
  // E      = convert_moduli::bulk_poisson::youngs (K,nu);
  // mu     = convert_moduli::youngs_poisson::mu     (E,nu);
  // lambda = convert_moduli::youngs_poisson::lambda (E,nu);
  // K      = convert_moduli::youngs_poisson::bulk   (E,nu); 
  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  // psi_c = 0.5*sig_c*sig_c*(1/E); // for tensile cohesive fracture 
  psi_c = 0.5*sig_c*sig_c*1/mu; 

  // For Quasi-quadrature
  b  = l_zero/2; 
  c0 = Pi; 
  a1 = xi/(c0*b) * G_c/psi_c; 
 
  phi_1 = pow(sig_c/sig_r, 2) * a1/xi; 
  g1    = 1/(1 + phi_1); 
  // g1 = 0; 

  new_cto    = elastic_cto;
  new_H_plus = psi_c;

  degrade = 1;
  degrade_deriv = 0;
  degrade_deriv_2 = 0; 

  lagrange_mult = 0; 

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldShearMono<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_strain = invert(new_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <>
void PhaseFieldShearMono<3>::update(const SymmetricTensor<2,3> &incr_strain)
{
  run_3D_update(incr_strain);
}
template <>
void PhaseFieldShearMono<2>::update(const SymmetricTensor<2,2> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  run_3D_update(incr_strain_2d);
}


template <int dim>
void PhaseFieldShearMono<dim>::update_phasefield(const double pf, const double pf_extra)
{
  // quasi-linear
  // degrade         = (1 - pf)/(1 - pf + M*pf);
  // degrade_deriv   = -M/pow(1 + pf*(M - 1), 2);
  // degrade_deriv_2 = (2*M*(M - 1))/pow(1 + pf*(M - 1), 3);

  // quasi-quadratic
  // degrade = pow(1 - pf, 2) / ( pow(1 - pf,2) + M*pf*(1 + pf) );
  // degrade_deriv = -(1 - pf)*(M + 3*M*pf) / pow(1 + pf*(-2 + pf + M + M*pf), 2);
  // degrade_deriv_2 = 2*M*(-3 + M - 3*pf*(-1 + (-1 + pf)*pf)*(1 + M)) /
  //                   pow(1 + pf*(-2 + pf + M + pf*M), 3);

  // quasi-linear for frictional shear fracture
  // degrade = ((1 - g1)*(1 - pf) + M*g1*pf)/((1 - g1)*(1 - pf) + M*pf);
  // degrade_deriv = -M*pow(1 - g1,2)/pow((1 - g1)*(1 - pf) + M*pf, 2); 
  // degrade_deriv_2 = 2*M*pow(1-g1,2)*(M + g1 -1)/pow((1 - g1)*(1 - pf) + M*pf, 3);

  // Wu's model for cohesive fracture 
  // degrade         = pow(1 - pf, 2)/( pow(1 - pf, 2) + a1*pf*(1 - 0.5*pf) ); 
  // degrade_deriv   = -a1*(4 - 4*pf)/ pow( (a1 - 2)*pf*pf + (4 - 2*a1)*pf - 2, 2); 
  // degrade_deriv_2 = -4*a1*( a1*(3*pf*pf - 6*pf + 4) - 6*(pf*pf - 2*pf + 1) )/ 
  //                   pow( a1*(pf - 2)*pf - 2*(pf*pf - 2*pf + 1), 3 );

  // Slip Weakening Model (general form)
  degrade         = ((1 - g1)*pow(1 - pf, 2) + a1*g1*pf*(1 + a2*pf)) / ( (1 - g1)*pow(1 - pf,2) + a1*pf*(1 + a2*pf)); 
  degrade_deriv   = (a1*pow(g1 - 1, 2)*(pf - 1)*(2*a2*pf + pf +1)) / pow(pf*pf*(a1*a2 + 1) + (a1 - 2)*pf - g1*(pf - 1)*(pf - 1) + 1, 2); 
  degrade_deriv_2 = (2*a1*pow(g1-1,2)*(a1*(a2*pf*(a2*(2*pf-3)*pf+pf*pf-3)-1)-(g1-1)*pow(pf-1,2)*(2*a2*pf+a2+pf+2))) / 
                    pow(-pf*(a1*a2*pf+a1+pf-2)+g1*pow(pf-1,2)-1, 3); 

  if (pf_extra > 0) 
    degrade_extra = ((1 - g1)*pow(1 - pf_extra, 2) + a1*g1*pf_extra*(1 + a2*pf_extra)) / ( (1 - g1)*pow(1 - pf_extra,2) + a1*pf_extra*(1 + a2*pf_extra)); 
  else 
    degrade_extra = degrade; 
  // degrade         = std::max((1-pf)*(1-pf), 1e-6);
  // degrade_deriv   = -2*(1-pf);
  // degrade_deriv_2 = 2;
}

template <int dim> // for deviatoric-volumetric decomposition 
void PhaseFieldShearMono<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain)
{
  // Update total strain 
  new_strain = old_strain + incr_strain; 
  // new_strain = incr_strain; // for adaptive remeshing 

  // Vol-Dev decomposition of the strain 
  double tr_eps = trace(new_strain);

  SymmetricTensor<2,3> new_strain_dev; 

  new_strain_dev = new_strain - 1./3.*tr_eps*eye; 

  // Construct stress from volumetric strain and deviatoric strain 
  SymmetricTensor<2,3> new_stress_dev = 2*mu*new_strain_dev;

  // if (tr_eps >= 0) 
  //  new_stress = K*degrade_extra*tr_eps*eye + new_stress_dev*degrade_extra; 
  // else 
   new_stress = K*tr_eps*eye + new_stress_dev*degrade_extra; 

  // Update strain histroy variable if current shear energy exceeds it  
  double strain_energy_degraded = 0; 

  // if (tr_eps >=0)
  //   strain_energy_degraded = 0.5*new_strain_dev*new_stress_dev + 0.5*K*tr_eps*tr_eps;
  // else 
    strain_energy_degraded = 0.5*new_strain_dev*new_stress_dev; 
  
  new_cto_pf = 0;

  // if (tr_eps > 0)
  //   new_cto_pf = degrade_deriv*(K*tr_eps*eye + 2*mu*new_strain_dev);
  // else
    new_cto_pf = degrade_deriv*2*mu*new_strain_dev;  

  new_H_plus = std::max(std::max(strain_energy_degraded, psi_c), old_H_plus); // for cohesive and shear fracture 

  new_H_plus_deriv_strain = 0; 

  if (new_H_plus > old_H_plus)
  {
    // if (tr_eps > 0)
    //   new_H_plus_deriv_strain = K*tr_eps*eye + 2*mu*new_strain_dev;
    // else 
      new_H_plus_deriv_strain = 2*mu*new_strain_dev; 
  }

  // Construct CTO 
  // if (tr_eps > 0) 
  //   new_cto = degrade_extra*K*eye_dyad_eye + 2*mu*degrade_extra*(big_eye - 1./3.*eye_dyad_eye); 
  // else 
    new_cto = K*eye_dyad_eye + 2*mu*degrade_extra*(big_eye - 1./3.*eye_dyad_eye); 
}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void PhaseFieldShearMono<dim>::save_state()
{
  old_stress = new_stress;
  old_strain = new_strain;
  old_H_plus = new_H_plus;
}

// Stress
template <>
SymmetricTensor<2,3> PhaseFieldShearMono<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> PhaseFieldShearMono<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldShearMono<dim>::full_stress()
{
  return new_stress;
}

// Strain
template <>
SymmetricTensor<2,3> PhaseFieldShearMono<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> PhaseFieldShearMono<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldShearMono<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double PhaseFieldShearMono<dim>::equiv_plastic_strain(){ return 0.; }

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> PhaseFieldShearMono<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> PhaseFieldShearMono<2>::cto()
{
  SymmetricTensor<4,2> new_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    new_cto_2d[i][j][k][l] = new_cto[i][j][k][l];

  return new_cto_2d;
}

template <> 
SymmetricTensor<2,3> PhaseFieldShearMono<3>::cto_pf()
{
  return new_cto_pf; 
}

template <> 
SymmetricTensor<2,2> PhaseFieldShearMono<2>::cto_pf()
{
  SymmetricTensor<2,2> new_cto_pf_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_cto_pf_2d[i][j] = new_cto_pf[i][j];

  return new_cto_pf_2d;
}

// Bulk modulus
template <int dim>
double PhaseFieldShearMono<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double PhaseFieldShearMono<dim>::shear_mod(){ return mu; }

// Critical fracture energy (critical energy release rate)
template <int dim>
double PhaseFieldShearMono<dim>::fracture_energy(){ return G_c; }

// Length scale parameter
template <int dim>
double PhaseFieldShearMono<dim>::length_scale(){ return l_zero/2.; }

// c0 coefficient for degradation function 
template <int dim> 
double PhaseFieldShearMono<dim>::c_0(){ return c0; }

// xi for crack density function 
template <int dim> 
double PhaseFieldShearMono<dim>::xi_output(){ return xi; }

// Crack driving forces
template <int dim>
double PhaseFieldShearMono<dim>::crack_driving_force()
{
  return new_H_plus;
}

// for monolithic 
template <>
SymmetricTensor<2,3> PhaseFieldShearMono<3>::crack_driving_force_deriv_strain()
{
  return new_H_plus_deriv_strain;
}

template <>
SymmetricTensor<2,2> PhaseFieldShearMono<2>::crack_driving_force_deriv_strain()
{
  SymmetricTensor<2,2> new_H_plus_deriv_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_H_plus_deriv_strain_2d[i][j] = new_H_plus_deriv_strain[i][j];

  return new_H_plus_deriv_strain_2d;
}

// Lagrange multiplier for irreversibility by augmented Lagrangian
template <int dim>
double PhaseFieldShearMono<dim>::lagrange_multiplier()
{
  return lagrange_mult;
}

template <int dim>
void PhaseFieldShearMono<dim>::set_lagrange_multiplier(double new_lagrange_mult)
{
  lagrange_mult = new_lagrange_mult;
}

template <int dim>
double PhaseFieldShearMono<dim>::degradation(){ return degrade; }

template <int dim> 
double PhaseFieldShearMono<dim>::degradation_deriv()
{
  return degrade_deriv; 
}

template <int dim>
double PhaseFieldShearMono<dim>::degradation_deriv_2()
{
  return degrade_deriv_2; 
}

// EXPLICIT INSTANTIATION
template class PhaseFieldShearMono<2>;
template class PhaseFieldShearMono<3>;
