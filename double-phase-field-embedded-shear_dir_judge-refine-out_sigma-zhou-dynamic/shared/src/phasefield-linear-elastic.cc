#include "geocentric/materials/phasefield-linear-elastic.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
PhaseFieldLinearElastic<dim>::PhaseFieldLinearElastic()
{ this->name_ = "Phase-field Linear Elastic"; }

template <int dim>
PhaseFieldLinearElastic<dim>::~PhaseFieldLinearElastic()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldLinearElastic<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldLinearElastic");
    prm.declare_entry("BulkMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergy","0.0",Patterns::Double());
    prm.declare_entry("LengthScale","0.0",Patterns::Double());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldLinearElastic<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldLinearElastic");
    K      = prm.get_double("BulkMod");
    nu     = prm.get_double("Poisson");
    G_c    = prm.get_double("FractureEnergy");
    l_zero = prm.get_double("LengthScale");
  prm.leave_subsection();

  mu     = convert_moduli::bulk_poisson::mu     (K,nu);
  lambda = convert_moduli::bulk_poisson::lambda (K,nu);
  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  new_cto    = elastic_cto;
  new_H_plus = 0;

  degrade = 1;
  degrade_deriv = 0;
  degrade_deriv_2 = 0;

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldLinearElastic<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_strain = invert(new_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <>
void PhaseFieldLinearElastic<3>::update(const SymmetricTensor<2,3> &incr_strain)
{
  run_3D_update(incr_strain);
}
template <>
void PhaseFieldLinearElastic<2>::update(const SymmetricTensor<2,2> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  run_3D_update(incr_strain_2d);
}


template <int dim>
void PhaseFieldLinearElastic<dim>::update_phasefield(const double pf)
{
  degrade         = std::max((1-pf)*(1-pf), 1e-6);
  degrade_deriv   = -2*(1-pf);
  degrade_deriv_2 = 2;
}


// -----------------------------------------------------------------------------
// ACTUAL 3D UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldLinearElastic<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain)
{
  // Update total strain
  new_strain = old_strain + incr_strain;
  // new_strain = incr_strain; // for adaptive remeshing

  // // Spectral decomposition of the strain
  // Tensor<1,3>               principal_strains;
  // std::vector<Tensor<1,3> > principal_directions(3);
  // eigen_decompose<3>(new_strain, principal_strains, principal_directions);

  // // Additively decompose into tensile and compressive strains
  // Tensor<1,3> principal_strains_plus, principal_strains_minus;

  // for (unsigned a=0; a<3; ++a)
  // {
  //   double eps_a = principal_strains[a];
  //   if (eps_a > 0) principal_strains_plus[a]  = eps_a;
  //   else           principal_strains_minus[a] = eps_a;
  // }

  // SymmetricTensor<2,3> new_strain_plus, new_strain_minus;
  // Tensor<2,3> m_aa;
  // for (unsigned a=0; a<3; ++a)
  // {
  //   m_aa = outer_product(principal_directions[a],principal_directions[a]);
  //   new_strain_plus  += principal_strains_plus[a]*symmetrize(m_aa);
  //   new_strain_minus += principal_strains_minus[a]*symmetrize(m_aa);
  // }

  // // Construct stress from tensile stress + compressive stress
  // double tr_eps = trace(new_strain);
  // double tr_eps_plus  = std::max(0.0, tr_eps);
  // double tr_eps_minus = tr_eps - tr_eps_plus;

  // Tensor<1,3> principal_stresses_plus, principal_stresses_minus;
  // for (unsigned a=0; a<3; ++a)
  // {
  //   principal_stresses_plus[a]  = lambda*tr_eps_plus  + 2*mu*principal_strains_plus[a];
  //   principal_stresses_minus[a] = lambda*tr_eps_minus + 2*mu*principal_strains_minus[a];
  // }
  // Tensor<1,3> principal_stresses = degrade*principal_stresses_plus + principal_stresses_minus;

  // new_stress = 0;
  // for (unsigned a=0; a<3; ++a)
  // {
  //   m_aa = outer_product(principal_directions[a],principal_directions[a]);
  //   new_stress += principal_stresses[a]*symmetrize(m_aa);
  // }

  // // Update strain history variable if current tensile strain energy exceeds it
  // double strain_energy_plus = 0.5 * principal_stresses_plus * principal_strains;
  // new_H_plus = (strain_energy_plus > old_H_plus) ? strain_energy_plus : old_H_plus;
  // // new_H_plus = strain_energy_plus; // for adaptive remeshing

  // // Construct consistent tangent operator (CTO)
  // SymmetricTensor<2,3> hessPsi_plus, hessPsi_minus;

  // for (unsigned i=0; i<3; ++i)
  // for (unsigned j=i; j<3; ++j)
  // {
  //   hessPsi_plus[i][j]  = lambda*(tr_eps >  0.0) + 2*mu*(i==j and principal_strains_plus[i]  > 0);
  //   hessPsi_minus[i][j] = lambda*(tr_eps <= 0.0) + 2*mu*(i==j and principal_strains_minus[i] < 0);
  // }

  // SymmetricTensor<2,3> elastic_a = degrade*hessPsi_plus + hessPsi_minus;

  // Tensor<2,3> m_bb, m_ab, m_ba;
  // Tensor<4,3> m_aabb, m_abab, m_abba;

  // new_cto = 0;

  // for (unsigned a=0; a<3; ++a)
  // for (unsigned b=0; b<3; ++b)
  // {
  //   m_aa = outer_product(principal_directions[a],principal_directions[a]);
  //   m_bb = outer_product(principal_directions[b],principal_directions[b]);
  //   m_ab = outer_product(principal_directions[a],principal_directions[b]);
  //   m_ba = outer_product(principal_directions[b],principal_directions[a]);

  //   m_aabb = outer_product(m_aa,m_bb);
  //   m_abab = outer_product(m_ab,m_ab);
  //   m_abba = outer_product(m_ab,m_ba);

  //   for (unsigned i=0; i<3; ++i)
  //   for (unsigned j=i; j<3; ++j)
  //   for (unsigned k=0; k<3; ++k)
  //   for (unsigned l=k; l<3; ++l)
  //     new_cto[i][j][k][l] += elastic_a[a][b] * m_aabb[i][j][k][l];

  //   if (a != b)
  //   {
  //     double strain_diff = principal_strains[b]  - principal_strains[a];
  //     double stress_diff = principal_stresses[b] - principal_stresses[a];
  //     double coeff;

  //     if (fabs(strain_diff) < tol)
  //       coeff = elastic_a[b][b] - elastic_a[a][b];
  //     else
  //       coeff = stress_diff / strain_diff;

  //     for (unsigned i=0; i<3; ++i)
  //     for (unsigned j=i; j<3; ++j)
  //     for (unsigned k=0; k<3; ++k)
  //     for (unsigned l=k; l<3; ++l)
  //       new_cto[i][j][k][l] += 0.5 * coeff * (m_abab[i][j][k][l] + m_abba[i][j][k][l]);
  //   }
  // } 

  // // dev-vol decomposition 
  // double tr_eps = trace(new_strain); 
  
  // double tr_eps_plus = std::max(0.0, tr_eps);

  // SymmetricTensor<2,3> new_strain_dev; 

  // new_strain_dev = new_strain - 1./3.*tr_eps*eye; 

  // // Construct stress from volumetric strain and deviatoric strain 
  // SymmetricTensor new_stress_dev = 2*degrade*mu*new_strain_dev;

  // if (tr_eps >= 0) 
  //  new_stress = K*degrade*tr_eps*eye + new_stress_dev; 
  // else 
  //  new_stress = K*tr_eps*eye + new_stress_dev; 

  // // Update strain histroy variable if current shear energy exceeds it  
  // // double strain_energy_degraded = 0.5*new_strain_dev*new_stress_dev; // only for shear degradation 
  // double strain_energy_degraded = 0.5*new_strain_dev*new_stress_dev + 0.5*K*tr_eps_plus*tr_eps_plus*degrade; // for tensile and shear degradation 
  // // new_H_plus = std::max(strain_energy_degraded, psi_c); // for cohesive fracture 
  // new_H_plus = strain_energy_degraded; 

  // if (new_H_plus < old_H_plus)
  //  new_H_plus = old_H_plus; 
  //   // new_H_plus = strain_energy_degraded; // for adaptive remeshing

  // // Construct CTO 
  // if (tr_eps_plus > 0) 
  //  new_cto = degrade*K*eye_dyad_eye + 2*mu*degrade*(big_eye - 1./3.*eye_dyad_eye); 
  // else 
  //  new_cto = K*eye_dyad_eye + 2*mu*degrade*(big_eye - 1./3.*eye_dyad_eye); 

  // Isotropic Degradation 
  new_cto = degrade * elastic_cto; 

  new_stress = new_cto * new_strain; 

  double strain_energy_degraded = 0.5*new_strain*new_stress; // for tensile and shear degradation 
  // new_H_plus = std::max(strain_energy_degraded, psi_c); // for cohesive fracture 
  new_H_plus = strain_energy_degraded; 

  if (new_H_plus < old_H_plus)
   new_H_plus = old_H_plus; 

}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void PhaseFieldLinearElastic<dim>::save_state()
{
  old_stress = new_stress;
  old_strain = new_strain;
  old_H_plus = new_H_plus;
}

// Stress
template <>
SymmetricTensor<2,3> PhaseFieldLinearElastic<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> PhaseFieldLinearElastic<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldLinearElastic<dim>::full_stress()
{
  return new_stress;
}

// Strain
template <>
SymmetricTensor<2,3> PhaseFieldLinearElastic<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> PhaseFieldLinearElastic<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldLinearElastic<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double PhaseFieldLinearElastic<dim>::equiv_plastic_strain(){ return 0.; }

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> PhaseFieldLinearElastic<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> PhaseFieldLinearElastic<2>::cto()
{
  SymmetricTensor<4,2> new_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    new_cto_2d[i][j][k][l] = new_cto[i][j][k][l];

  return new_cto_2d;
}


// Bulk modulus
template <int dim>
double PhaseFieldLinearElastic<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double PhaseFieldLinearElastic<dim>::shear_mod(){ return mu; }

// Critical fracture energy (critical energy release rate)
template <int dim>
double PhaseFieldLinearElastic<dim>::fracture_energy(){ return G_c; }

// Length scale parameter
template <int dim>
double PhaseFieldLinearElastic<dim>::length_scale(){ return l_zero; }

// Degradation value 
template <int dim> 
double PhaseFieldLinearElastic<dim>::g_d(){ return degrade; }

template <int dim> 
double PhaseFieldLinearElastic<dim>::gd_deriv(){ return degrade_deriv; }

// Crack driving forces
template <int dim>
double PhaseFieldLinearElastic<dim>::crack_driving_force()
{
  return degrade_deriv*new_H_plus;
}

template <int dim>
double PhaseFieldLinearElastic<dim>::crack_driving_force_deriv()
{
  return degrade_deriv_2*new_H_plus;
}

template <int dim>
void PhaseFieldLinearElastic<dim>::set_strain_energy_plus(double H_plus)
{
  new_H_plus = H_plus;
  old_H_plus = H_plus;
}

// EXPLICIT INSTANTIATION
template class PhaseFieldLinearElastic<2>;
template class PhaseFieldLinearElastic<3>;
