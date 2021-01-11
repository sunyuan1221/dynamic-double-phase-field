#include "geocentric/materials/phasefield-cohesive.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
PhaseFieldCohesive<dim>::PhaseFieldCohesive()
{ this->name_ = "Phase-field Cohesive Zone"; }

template <int dim>
PhaseFieldCohesive<dim>::~PhaseFieldCohesive()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldCohesive<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldCohesive");
    // prm.declare_entry("BulkMod","0.0",Patterns::Double());
    prm.declare_entry("YoungsMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergy","0.0",Patterns::Double());
    prm.declare_entry("LengthScale","0.0",Patterns::Double());
    prm.declare_entry("CriticalStress","0.0",Patterns::Double());
    prm.declare_entry("ShapeParameter","1.0",Patterns::Double());
    prm.declare_entry("CohesiveFracture","true",Patterns::Bool());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldCohesive<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldCohesive");
    // K      = prm.get_double("BulkMod");
    E     = prm.get_double("YoungsMod");
    nu    = prm.get_double("Poisson");
    G_c   = prm.get_double("FractureEnergy");
    L     = prm.get_double("LengthScale");
    sig_c = prm.get_double("CriticalStress");
    exp_p = prm.get_double("ShapeParameter");
    cohesive_fracture = prm.get_bool("CohesiveFracture");
  prm.leave_subsection();

  // When K is given
  // {
  //   mu     = convert_moduli::bulk_poisson::mu     (K,nu);
  //   lambda = convert_moduli::bulk_poisson::lambda (K,nu);
  // }

  // When E is given
  {
    K      = convert_moduli::youngs_poisson::bulk   (E,nu);
    mu     = convert_moduli::youngs_poisson::mu     (E,nu);
    lambda = convert_moduli::youngs_poisson::lambda (E,nu);
  }
  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  psi_c = (0.5/E)*sig_c*sig_c;

  M = 3*G_c/(4*L*psi_c);

  if (cohesive_fracture)
  {
    double L_fpz = E*G_c*pow(sig_c,-2);
    if (L > (1.5/(exp_p+2))*L_fpz)
      std::cout << "ERROR: The regularization length exceeds the upper bound" << std::endl;
  }

  new_cto = elastic_cto;

  // new_H_plus = (cohesive_fracture) ? psi_c : 0.0;
  new_H_plus = (cohesive_fracture) ? psi_c : 0.0;

  g_d = 1;
  g_d_deriv = 0;
  g_d_deriv_2 = 0;

  K_damaged = K;
  mu_damaged = mu;

  pressure_contribution = 0;

  lagrange_mult = 0;

  initial_crack_point = false;

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldCohesive<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_strain = invert(new_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <>
void PhaseFieldCohesive<3>::update(const SymmetricTensor<2,3> &incr_strain)
{
  run_3D_update(incr_strain);
}
template <>
void PhaseFieldCohesive<2>::update(const SymmetricTensor<2,2> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  run_3D_update(incr_strain_2d);
}


template <int dim>
void PhaseFieldCohesive<dim>::update_phasefield(const double pf, const double pf_extra)
{
  if (cohesive_fracture)
  {
    // Quasi-quadratic (Lorentz)
    g_d = pow(1 - pf, 2) / ( pow(1 - pf, 2) + M*pf*(1 + exp_p*pf) );
    g_d_deriv = (-1 + pf)*M*(1 + pf + 2*pf*exp_p)
                * pow(1 + pf*(-2 + pf + M + pf*M*exp_p), -2);
    g_d_deriv_2 = -2*M*(2 - M + exp_p + pf*(-3 + pf*pf + pf*(-3 + 2*pf)*exp_p)*(1 + M*exp_p))
                  * pow(1 + pf*(-2 + pf + M + pf*M*exp_p), -3);

    if (pf_extra > 0)
      g_d_extra = pow(1 - pf_extra, 2) / ( pow(1 - pf_extra,2) + M*pf_extra*(1 + exp_p*pf_extra) ); // for quasi-monolithic
    else
      g_d_extra = g_d;
  }
  else
  {
    // Standard phase-field
    g_d = pow(1 - pf, 2);
    g_d_deriv = -2*(1 - pf);
    g_d_deriv_2 = 2;

    if (pf_extra > 0)
      g_d_extra = pow(1 - pf_extra, 2); // for quasi-monolithic
    else
      g_d_extra = g_d;
  }
}

template <int dim>
void PhaseFieldCohesive<dim>::update_pressure_contribution(const double pressure_term)
{
  pressure_contribution = pressure_term;
}


// -----------------------------------------------------------------------------
// ACTUAL 3D UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldCohesive<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain)
{
  // Update total strain
  new_strain = old_strain + incr_strain;
  // new_strain = incr_strain; // for adaptive remeshing

  // 1) vol--dev split (Amor et al., 2009)
  // TODO: clean below

        // double vol_strain = trace(new_strain);
        // double strain_energy_plus = 0;

        // SymmetricTensor<2,3> dev_strain = new_strain - one_third*vol_strain*eye;

        // if (vol_strain > 0)
        // {
        //   strain_energy_plus = 0.5*K*vol_strain*vol_strain + mu*scalar_product(dev_strain,dev_strain);

        //   K_damaged  = g_d_extra*K;
        //   mu_damaged = g_d_extra*mu;
        // }
        // else
        // {
        //   strain_energy_plus = mu*scalar_product(dev_strain,dev_strain);

        //   K_damaged  = K;
        //   mu_damaged = g_d_extra*mu;
        // }

        // new_stress = K_damaged*vol_strain*eye + 2*mu_damaged*dev_strain;

        // if (cohesive_fracture)
        // {
        //   new_H_plus = std::max( std::max(strain_energy_plus + pressure_contribution, psi_c), old_H_plus);
        // }
        // else
        // {
        //   new_H_plus = std::max( strain_energy_plus + pressure_contribution, old_H_plus);
        // }


        // // Construct consistent tangent operator (CTO)
        // if (fabs(g_d - 1.0) < tol)
        // {
        //   new_cto = elastic_cto;
        // }
        // else
        // {
        //   new_cto = (K_damaged - two_thirds*mu_damaged)*eye_dyad_eye + 2*mu_damaged*big_eye;
        // }

        // new_cto_pf = 0;

        // if (vol_strain > 0)
        // {
        //   new_cto_pf = g_d_deriv*(K*vol_strain*eye + 2*mu*dev_strain);
        // }
        // else
        // {
        //   new_cto_pf = g_d_deriv*2*mu*dev_strain;
        // }

        // new_H_plus_deriv_strain = 0;

        // if (new_H_plus > old_H_plus)
        // {
        //   if (vol_strain > 0)
        //   {
        //     new_H_plus_deriv_strain = K*vol_strain*eye + 2*mu*dev_strain;
        //   }
        //   else
        //   {
        //     new_H_plus_deriv_strain = 2*mu*dev_strain;
        //   }
        // }



  // 2) Principal strain split (Miehe)

  // Spectral decomposition of the strain
  Tensor<1,3>               principal_strains;
  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(new_strain, principal_strains, principal_directions);

  // Additively decompose into tensile and compressive strains
  Tensor<1,3> principal_strains_plus, principal_strains_minus;

  for (unsigned a=0; a<3; ++a)
  {
    double eps_a = principal_strains[a];
    if (eps_a > 0) principal_strains_plus[a]  = eps_a;
    else           principal_strains_minus[a] = eps_a;
  }

  SymmetricTensor<2,3> new_strain_plus, new_strain_minus;
  Tensor<2,3> m_aa;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_strain_plus  += principal_strains_plus[a]*symmetrize(m_aa);
    new_strain_minus += principal_strains_minus[a]*symmetrize(m_aa);
  }

  // Construct stress from tensile stress + compressive stress
  double tr_eps = trace(new_strain);
  double tr_eps_plus  = std::max(0.0, tr_eps);
  double tr_eps_minus = tr_eps - tr_eps_plus;

  Tensor<1,3> principal_stresses_plus, principal_stresses_minus;
  for (unsigned a=0; a<3; ++a)
  {
    principal_stresses_plus[a]  = lambda*tr_eps_plus  + 2*mu*principal_strains_plus[a];
    principal_stresses_minus[a] = lambda*tr_eps_minus + 2*mu*principal_strains_minus[a];
  }
  Tensor<1,3> principal_stresses = g_d_extra*principal_stresses_plus + principal_stresses_minus;

  new_stress = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_stress += principal_stresses[a]*symmetrize(m_aa);
  }

  // Update strain history variable if current tensile strain energy exceeds it
  double strain_energy_plus = 0.5 * principal_stresses_plus * principal_strains;

  if (cohesive_fracture)
  {
    new_H_plus = std::max( std::max(strain_energy_plus + pressure_contribution, psi_c), old_H_plus);
    // new_H_plus = std::max(strain_energy_plus + pressure_contribution, psi_c);
  }
  else
  {
    new_H_plus = std::max( strain_energy_plus + pressure_contribution, old_H_plus);
  }

  // Derivatives for monolithic phase-field solvers
  // new_H_plus_deriv_strain = 0;

  // if (new_H_plus > old_H_plus)
  // {
  //   for (unsigned a=0; a<3; ++a)
  //   {
  //     m_aa = outer_product(principal_directions[a],principal_directions[a]);
  //     new_H_plus_deriv_strain += principal_stresses_plus[a]*symmetrize(m_aa);
  //   }
  // }

  // new_cto_pf = 0;

  // if (fabs(g_d - 1.0) < tol)
  // {
  //   for (unsigned a=0; a<3; ++a)
  //   {
  //     m_aa = outer_product(principal_directions[a],principal_directions[a]);
  //     new_cto_pf += g_d_deriv*principal_stresses_plus[a]*symmetrize(m_aa);
  //   }
  // }

  // Construct consistent tangent operator (CTO)

  // return undamaged part early
  if (fabs(g_d - 1.0) < tol)
  {
    new_cto = elastic_cto;
    return;
  }

  SymmetricTensor<2,3> hessPsi_plus, hessPsi_minus;

  for (unsigned i=0; i<3; ++i)
  for (unsigned j=i; j<3; ++j)
  {
    hessPsi_plus[i][j]  = lambda*(tr_eps >  0.0) + 2*mu*(i==j and principal_strains_plus[i]  > 0);
    hessPsi_minus[i][j] = lambda*(tr_eps <= 0.0) + 2*mu*(i==j and principal_strains_minus[i] < 0);
  }

  SymmetricTensor<2,3> elastic_a = g_d_extra*hessPsi_plus + hessPsi_minus;

  Tensor<2,3> m_bb, m_ab, m_ba;
  Tensor<4,3> m_aabb, m_abab, m_abba;

  new_cto = 0;

  for (unsigned a=0; a<3; ++a)
  for (unsigned b=0; b<3; ++b)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    m_bb = outer_product(principal_directions[b],principal_directions[b]);
    m_ab = outer_product(principal_directions[a],principal_directions[b]);
    m_ba = outer_product(principal_directions[b],principal_directions[a]);

    m_aabb = outer_product(m_aa,m_bb);
    m_abab = outer_product(m_ab,m_ab);
    m_abba = outer_product(m_ab,m_ba);

    for (unsigned i=0; i<3; ++i)
    for (unsigned j=i; j<3; ++j)
    for (unsigned k=0; k<3; ++k)
    for (unsigned l=k; l<3; ++l)
      new_cto[i][j][k][l] += elastic_a[a][b] * m_aabb[i][j][k][l];

    if (a != b)
    {
      double strain_diff = principal_strains[b]  - principal_strains[a];
      double stress_diff = principal_stresses[b] - principal_stresses[a];
      double coeff;

      if (fabs(strain_diff) < tol)
        coeff = elastic_a[b][b] - elastic_a[a][b];
      else
        coeff = stress_diff / strain_diff;

      for (unsigned i=0; i<3; ++i)
      for (unsigned j=i; j<3; ++j)
      for (unsigned k=0; k<3; ++k)
      for (unsigned l=k; l<3; ++l)
        new_cto[i][j][k][l] += 0.5 * coeff * (m_abab[i][j][k][l] + m_abba[i][j][k][l]);
    }
  }
}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void PhaseFieldCohesive<dim>::save_state()
{
  old_stress = new_stress;
  old_strain = new_strain;
  old_H_plus = new_H_plus;
}

// Stress
template <>
SymmetricTensor<2,3> PhaseFieldCohesive<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> PhaseFieldCohesive<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldCohesive<dim>::full_stress()
{
  return new_stress;
}

// Strain
template <>
SymmetricTensor<2,3> PhaseFieldCohesive<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> PhaseFieldCohesive<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldCohesive<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double PhaseFieldCohesive<dim>::equiv_plastic_strain(){ return 0.; }

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> PhaseFieldCohesive<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> PhaseFieldCohesive<2>::cto()
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
SymmetricTensor<2,3> PhaseFieldCohesive<3>::cto_pf()
{
  return new_cto_pf;
}
template <>
SymmetricTensor<2,2> PhaseFieldCohesive<2>::cto_pf()
{
  SymmetricTensor<2,2> new_cto_pf_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_cto_pf_2d[i][j] = new_cto_pf[i][j];

  return new_cto_pf_2d;
}


// Bulk modulus
template <int dim>
double PhaseFieldCohesive<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double PhaseFieldCohesive<dim>::shear_mod(){ return mu; }

// Critical fracture energy (critical energy release rate)
template <int dim>
double PhaseFieldCohesive<dim>::fracture_energy(){ return G_c; }

// Length scale parameter
template <int dim>
double PhaseFieldCohesive<dim>::length_scale(){ return L; }

// Crack driving forces
template <int dim>
double PhaseFieldCohesive<dim>::crack_driving_force()
{
  return new_H_plus;
}

// Crack driving forces over strain (for monolithic)
template <>
SymmetricTensor<2,3> PhaseFieldCohesive<3>::crack_driving_force_deriv_strain()
{
  return new_H_plus_deriv_strain;
}

template <>
SymmetricTensor<2,2> PhaseFieldCohesive<2>::crack_driving_force_deriv_strain()
{
  SymmetricTensor<2,2> new_H_plus_deriv_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_H_plus_deriv_strain_2d[i][j] = new_H_plus_deriv_strain[i][j];

  return new_H_plus_deriv_strain_2d;
}

template <int dim>
double PhaseFieldCohesive<dim>::old_crack_driving_force()
{
  return old_H_plus;
}

template <int dim> 
double PhaseFieldCohesive<dim>::H_plus_output(){ return new_H_plus; }


// Degradation function values (used for hydrofrac)
template <int dim>
double PhaseFieldCohesive<dim>::degrade()
{
  return g_d;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade_extra()
{
  return g_d_extra;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade_deriv()
{
  return g_d_deriv;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade_deriv_2()
{
  return g_d_deriv_2;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade(const double pf)
{
  double tmp = pow(1 - pf, 2) / ( pow(1 - pf, 2) + M*pf*(1 + exp_p*pf) );

  return tmp;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade_deriv(const double pf)
{
  double tmp = (-1 + pf)*M*(1 + pf + 2*pf*exp_p)
               * pow(1 + pf*(-2 + pf + M + pf*M*exp_p), -2);

  return tmp;
}

template <int dim>
double PhaseFieldCohesive<dim>::degrade_deriv_2(const double pf)
{
  double tmp = -2*M*(2 - M + exp_p + pf*(-3 + pf*pf + pf*(-3 + 2*pf)*exp_p)*(1 + M*exp_p))
                * pow(1 + pf*(-2 + pf + M + pf*M*exp_p), -3);

  return tmp;
}


// Set
template <int dim>
void PhaseFieldCohesive<dim>::set_fracture_energy(double new_G_c)
{
  G_c = new_G_c;
}

template <int dim>
void PhaseFieldCohesive<dim>::set_youngs_mod(double new_E)
{
  E = new_E;

  K      = convert_moduli::youngs_poisson::bulk   (E,nu);
  mu     = convert_moduli::youngs_poisson::mu     (E,nu);
  lambda = convert_moduli::youngs_poisson::lambda (E,nu);

  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;
}

template <int dim>
void PhaseFieldCohesive<dim>::set_crack_driving_force(double new_crack_driving_force)
{
  new_H_plus = new_crack_driving_force;
  old_H_plus = new_H_plus;
}

// Lagrange multiplier for irreversibility by augmented Lagrangian
template <int dim>
double PhaseFieldCohesive<dim>::lagrange_multiplier()
{
  return lagrange_mult;
}

template <int dim>
void PhaseFieldCohesive<dim>::set_lagrange_multiplier(double new_lagrange_mult)
{
  lagrange_mult = new_lagrange_mult;
}


// EXPLICIT INSTANTIATION
template class PhaseFieldCohesive<2>;
template class PhaseFieldCohesive<3>;
