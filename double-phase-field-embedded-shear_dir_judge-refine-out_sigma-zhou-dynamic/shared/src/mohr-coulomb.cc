#include "geocentric/materials/mohr-coulomb.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
MohrCoulomb<dim>::MohrCoulomb()
{ this->name_ = "Mohr-Coulomb"; }

template <int dim>
MohrCoulomb<dim>::~MohrCoulomb()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void MohrCoulomb<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("MohrCoulomb");
    prm.declare_entry("BulkMod","1000",Patterns::Double());
    prm.declare_entry("Poisson","0.25",Patterns::Double());
    prm.declare_entry("Cohesion","0.0",Patterns::Double());
    prm.declare_entry("FrictionAngle","30.0",Patterns::Double());
    prm.declare_entry("DilationAngle","10.0",Patterns::Double());
    prm.declare_entry("ShapeFactor","0.0",Patterns::Double());
    prm.declare_entry("Verbose","false",Patterns::Bool());
    prm.declare_entry("FiniteStrain","false",Patterns::Bool());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// INITIALIZE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void MohrCoulomb<dim>::initialize(ParameterHandler &prm)
{
  // Parse Parameters
  prm.enter_subsection("MohrCoulomb");
    K        = prm.get_double("BulkMod");
    nu       = prm.get_double("Poisson");
    cohesion = prm.get_double("Cohesion");
    phi      = prm.get_double("FrictionAngle")*M_PI/180.;
    psi      = prm.get_double("DilationAngle")*M_PI/180.;
    shape    = prm.get_double("ShapeFactor");
    verbose  = prm.get_bool("Verbose");
    finite_strain = prm.get_bool("FiniteStrain");
  prm.leave_subsection();

  // Elastic Tangent Operator
  mu     = convert_moduli::bulk_poisson::mu(K,nu);
  lambda = convert_moduli::bulk_poisson::lambda(K,nu);
  youngs = convert_moduli::bulk_poisson::youngs(K,nu);

  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  // Elastic Tangent for Principal Values
  for (unsigned i=0; i<3; ++i)
  for (unsigned j=0; j<3; ++j)
  {
    elastic_a[i][j]     = ( (i==j) ? lambda+2*mu : lambda);
    elastic_a_inv[i][j] = ( (i==j) ? 1.0 : -1*nu) / youngs;
  }

  // Plasticity parameters
  cosphi    = std::cos(phi);
  cospsi    = std::cos(psi);
  sinphi    = std::sin(phi);
  sinpsi    = std::sin(psi);
  shape_phi = std::pow(shape*cohesion*cosphi,2);
  shape_psi = std::pow(shape*cohesion*cospsi,2);

  // Plastic Internal Variable / Strain
  new_multiplier = 0;
  new_plastic_strain = 0;
  new_equiv_plastic_strain = 0;

  // Kinematics
  if (finite_strain)
  {
    new_elastic_b = eye;
    new_J = 1.;
  }

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void MohrCoulomb<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_elastic_strain = invert(elastic_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE
// -----------------------------------------------------------------------------
// infinitesimal strain
template <int dim>
void MohrCoulomb<dim>::update(const SymmetricTensor<2,dim> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_3d;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=i; j<dim; ++j)
    incr_strain_3d[i][j] = incr_strain[i][j];

  new_elastic_strain = old_elastic_strain + incr_strain_3d;

  run_3D_update();
}

// finite strain
template <int dim>
void MohrCoulomb<dim>::update(const Tensor<2,dim> &rel_def_grad)
{
  Tensor<2,3> rel_def_grad_3d;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
    rel_def_grad_3d[i][j] = rel_def_grad[i][j];

  if (dim == 2)
    rel_def_grad_3d[2][2] = 1;

  new_elastic_b = rel_def_grad_3d * old_elastic_b * transpose(rel_def_grad_3d);

  Tensor<1,3> principal_stretches_square;
  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(new_elastic_b, principal_stretches_square, principal_directions);

  Tensor<2,3> m_aa;
  new_elastic_strain = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_elastic_strain += 0.5*log(principal_stretches_square[a])*symmetrize(m_aa);
  }

  run_3D_update();

  new_stress = (1/new_J) * new_stress; // since logarithmic strains gives kirchhoff stress

  eigen_decompose<3>(new_elastic_strain, principal_e_strains, principal_directions);

  new_elastic_b = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_elastic_b += exp(2*principal_e_strains[a])*m_aa;
  }
}


// -----------------------------------------------------------------------------
// ACTUAL 3D STRESS UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void MohrCoulomb<dim>::run_3D_update()
{
  // Compute a trial state in which the increment is assumed to be fully elastic
  new_plastic_strain       = old_plastic_strain;
  new_equiv_plastic_strain = old_equiv_plastic_strain;
  new_strain               = new_elastic_strain + new_plastic_strain;
  new_multiplier           = old_multiplier;

  if (finite_strain)
    new_J = exp(trace(new_strain));

  new_stress = new_cto*new_elastic_strain;

  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(new_stress, principal_stresses, principal_directions);

  double yield = F(principal_stresses);

  // If Elastic (yield F <= 0, elastic response)
  // -------------------------------------------
  if (yield < tol)
  {
    new_cto = elastic_cto;

    return;
  }

  // Else Plastic
  // Need to perform plastic iterations using Newton algorithm
  // ---------------------------------------------------------

  // Return Mapping in Principal Axes
  //
  // Eigen (spectral) decompose stress
  // and assign initial state vectors
  // slot 0 = Sigma 1
  // slot 1 = Sigma 2
  // slot 2 = Sigma 3
  // slot 3 = Delta lambda (plastic multiplier)
  // See Sect. 4.11 in Plasticity Modeling and Computation

  // Save trial strains (trial: not updated during Newton iteration)
  principal_e_strains = elastic_a_inv * principal_stresses;
  Tensor<1,3> trial_principal_e_strains = principal_e_strains;

  Tensor<1,3> gradF, gradG;
  SymmetricTensor<2,3> hessG;

  // Newton iteration
  FullMatrix<double> jacobian(4,4);
  FullMatrix<double> jacobian_inv(4,4);

  Vector<double>     xdelta(4);
  Vector<double>     residual(4);
                     residual(3) = yield;

  double dLambda = 0;

  double resid     = residual.l2_norm();
  double resid_ini = 1 + resid;

  ConvergenceTable convergence_table;

  for (unsigned iter=0; iter<15; ++iter)
  {
    // Derivatives of F and G w.r.t. principal stresses
    gradF = grad_F(principal_stresses);
    gradG = grad_F(principal_stresses,/* potential = */true);

    // Residual
    for (unsigned i=0; i<3; ++i)
      residual(i) = principal_e_strains[i]
                    - trial_principal_e_strains[i]
                    + dLambda*gradG[i];

    residual(3) = F(principal_stresses);

    resid = residual.l2_norm();

    if (verbose)
    {
      convergence_table.add_value("Iter",iter);
      convergence_table.add_value("RHS Norm",resid);
      convergence_table.add_value("Reduction",resid/resid_ini);
    }

    if (resid/resid_ini < tol or resid < tol)
      break;

    // Jacobian matrix for local return mapping
    for (unsigned i=0; i<3; ++i)
    {
      jacobian(i,3) = gradG[i];
      jacobian(3,i) = gradF[i];

      for (unsigned j=0; j<3; ++j)
        jacobian(i,j) = elastic_a_inv[i][j]; // hessian is zero for M-C
    }
    jacobian(3,3) = 0;

    // Newton update: xnew = xcurrent - xdelta, xdelta = jacobian_inv*residual
    jacobian_inv.invert(jacobian);
    jacobian_inv.vmult(xdelta,residual);

    // Update principal stresses and strains
    principal_stresses[0] -= xdelta(0);
    principal_stresses[1] -= xdelta(1);
    principal_stresses[2] -= xdelta(2);
    dLambda               -= xdelta(3);

    principal_e_strains = elastic_a_inv * principal_stresses;
  } // end of return mapping

  if (verbose)
  {
    convergence_table.set_scientific("RHS Norm",true);
    convergence_table.set_scientific("Reduction",true);
    std::cout << "Convergence during local return mapping of " << this->name() << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    convergence_table.write_text(std::cout,TableHandler::TextOutputFormat::org_mode_table);
  }

  // Upon convergence, reconstruct stresses, strains, cto, etc.
  // ----------------------------------------------------------

  // Spectral directions
  Tensor<2,3> m_aa, m_bb, m_ab, m_ba;
  Tensor<4,3> m_aabb, m_abab, m_abba;

  // Reconstruct stress
  new_stress = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_stress += principal_stresses[a]*symmetrize(m_aa);
  }

  // Derivative of potential w.r.t. stress
  SymmetricTensor<2,3> dGdsigma;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    dGdsigma += gradG[a]*symmetrize(m_aa);
  }

  // Reconstruct strain
  new_elastic_strain       -= dLambda*dGdsigma;
  new_plastic_strain       += dLambda*dGdsigma;
  new_equiv_plastic_strain += root23*dLambda*dGdsigma.norm();
  new_multiplier           += dLambda;

  // Reconstruct consistent tangent opterator (see Sect. 4.12)
  new_cto = 0;

  FullMatrix<double> elastic_a_trial(3,3);
  for (unsigned a=0; a<3; ++a)
  for (unsigned b=0; b<3; ++b)
    for (unsigned i=0; i<4; ++i)
    {
      if (i == b)
        elastic_a_trial[a][b] += jacobian_inv[a][i];
    }

  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);

    for (unsigned b=0; b<3; ++b)
    {
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
        new_cto[i][j][k][l] += elastic_a_trial[a][b] * m_aabb[i][j][k][l];

      if (a != b)
      {
        if (not finite_strain)
        {
          double strain_diff = trial_principal_e_strains[b] - trial_principal_e_strains[a];
          double stress_diff = principal_stresses[b] - principal_stresses[a];
          double coeff;

          if (fabs(strain_diff) < tol)
            coeff = 0.5 * (elastic_a_trial[b][a] - elastic_a_trial[a][a]);
          else
            coeff = 0.5 * stress_diff / strain_diff;

          for (unsigned i=0; i<3; ++i)
          for (unsigned j=i; j<3; ++j)
          for (unsigned k=0; k<3; ++k)
          for (unsigned l=k; l<3; ++l)
            new_cto[i][j][k][l] += coeff * (m_abab[i][j][k][l] + m_abba[i][j][k][l]);
        }
        else
        {
          double stretch_sq_a = exp(2*trial_principal_e_strains[a]);
          double stretch_sq_b = exp(2*trial_principal_e_strains[b]);
          double stretch_diff = stretch_sq_b - stretch_sq_a;
          double stress_diff  = principal_stresses[b] - principal_stresses[a];
          double coeff;

          if (fabs(stretch_diff) < tol)
            coeff = (0.5/stretch_sq_a)*(elastic_a_trial[b][a] - elastic_a_trial[a][a]);
          else
            coeff = stress_diff / stretch_diff;

          for (unsigned i=0; i<3; ++i)
          for (unsigned j=i; j<3; ++j)
          for (unsigned k=0; k<3; ++k)
          for (unsigned l=k; l<3; ++l)
            new_cto[i][j][k][l] += coeff * (stretch_sq_b*m_abab[i][j][k][l] + stretch_sq_a*m_abba[i][j][k][l]);
        }
      }
    }
  }

}


// -----------------------------------------------------------------------------
// YIELD AND POTENTIAL SURFACES AND THEIR DERIVATIES
// -----------------------------------------------------------------------------
// Yield/potential function
// (angle is phi for yield, psi for potential)
template <int dim>
double MohrCoulomb<dim>::F(const Tensor<1,3> &sigma, bool potential)
{
  // Note: already sorted as sigma[0] < sigma[1] < sigma[2]
  // (sig[0] is the maximum principal stress in compression)
  const double sigma_m = 0.5*(sigma[2] + sigma[0]);
  const double tau_m   = 0.5*(sigma[2] - sigma[0]);

  if (not potential)
    return (sqrt(tau_m*tau_m + shape_phi) + sigma_m*sinphi - cohesion*cosphi);
  else
    return (sqrt(tau_m*tau_m + shape_psi) + sigma_m*sinpsi - cohesion*cospsi);
}

// Derivative of yield/potential function w.r.t. principal stresses
// (angle is phi for yield, psi for potential)
template <int dim>
Tensor<1,3> MohrCoulomb<dim>::grad_F(const Tensor<1,3> &sigma, bool potential)
{
  const double tau_m = 0.5*(sigma[2] - sigma[0]);

  Tensor<1,3> gradF;
  if (not potential)
  {
    const double tmp = tau_m/sqrt(tau_m*tau_m + shape_phi);
    gradF[0] = 0.5*(-tmp + sinphi);
    gradF[2] = 0.5*( tmp + sinphi);
  }
  else
  {
    const double tmp = tau_m/sqrt(tau_m*tau_m + shape_psi);
    gradF[0] = 0.5*(-tmp + sinpsi);
    gradF[2] = 0.5*( tmp + sinpsi);
  }

  return gradF;
}


// -----------------------------------------------------------------------------
// COMMON METHODS: SAVE, GET, ETC.
// -----------------------------------------------------------------------------

// Save state
// -----------------------------------------------------------------------------
template <int dim>
void MohrCoulomb<dim>::save_state()
{
  old_stress               = new_stress;
  old_elastic_strain       = new_elastic_strain;
  old_plastic_strain       = new_plastic_strain;
  old_strain               = new_strain;
  old_equiv_plastic_strain = new_equiv_plastic_strain;
  old_multiplier           = new_multiplier;

  if (finite_strain)
  {
    old_elastic_b = new_elastic_b;
    old_J         = new_J;
  }
}


// Stresses
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> MohrCoulomb<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> MohrCoulomb<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];
  return new_stress_2d;
}
template <int dim>
SymmetricTensor<2,3> MohrCoulomb<dim>::full_stress()
{
  return new_stress;
}


// Strains
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> MohrCoulomb<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> MohrCoulomb<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];
  return new_strain_2d;
}
template <int dim>
SymmetricTensor<2,3> MohrCoulomb<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double MohrCoulomb<dim>::equiv_plastic_strain() { return new_equiv_plastic_strain; }

// Volumetric strain
template <int dim>
double MohrCoulomb<dim>::vol_strain(){ return trace(new_plastic_strain); }

// Deviatoric strain
template <int dim>
double MohrCoulomb<dim>::dev_strain()
{
  auto dev_strain = new_plastic_strain - one_third*trace(new_plastic_strain)*eye;
  return (root23*dev_strain.norm());
}

template <int dim>
double MohrCoulomb<dim>::old_jacobian(){ return old_J; }

// Moduli
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<4,3> MohrCoulomb<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> MohrCoulomb<2>::cto()
{
  SymmetricTensor<4,2> new_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    new_cto_2d[i][j][k][l] = new_cto[i][j][k][l];
  return new_cto_2d;
}

template <int dim>
Tensor<4,dim> MohrCoulomb<dim>::cto_kirchhoff()
{
  Tensor<4,dim> cto_kirchhoff;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
  for (unsigned k=0; k<dim; ++k)
  for (unsigned l=0; l<dim; ++l)
    cto_kirchhoff[i][j][k][l] = new_cto[i][j][k][l] - (j==k)*new_J*new_stress[i][l];

  return cto_kirchhoff;
}

// Bulk Modulus
template <int dim>
double MohrCoulomb<dim>::bulk_mod(){ return K; }

// Shear Modulus (TODO: Plastic shear modulus)
template <int dim>
double MohrCoulomb<dim>::shear_mod(){ return mu; }


// Explicit Instantiation
template class MohrCoulomb<2>;
template class MohrCoulomb<3>;
