#include <geocentric/materials/drucker-prager.h>

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
DruckerPrager<dim>::DruckerPrager()
  :
  xtrial(3),
  xcurrent(3),
  residual(3),
  jacobian(3,3),
  jacobian_inv(3,3),
  deriv(2,2)
{ this->name_ = "Drucker-Prager"; }

template <int dim>
DruckerPrager<dim>::~DruckerPrager()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("DruckerPrager");
    prm.declare_entry("BulkMod","1000",Patterns::Double());
    prm.declare_entry("Poisson","0.25",Patterns::Double());
    prm.declare_entry("Cohesion","0.0",Patterns::Double());
    prm.declare_entry("FrictionAngle","30.0",Patterns::Double());
    prm.declare_entry("DilationAngle","10.0",Patterns::Double());
    prm.declare_entry("CohesionSoftening","0.0",Patterns::Double());
    prm.declare_entry("ShapeFactor","0",Patterns::Double());
    prm.declare_entry("Corners","compression",Patterns::Anything());
    prm.declare_entry("Verbose","false",Patterns::Bool());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// INITIALIZE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::initialize(ParameterHandler &prm)
{
  // Parse Parameters
  prm.enter_subsection("DruckerPrager");
    K        = prm.get_double("BulkMod");
    nu       = prm.get_double("Poisson");
    c_zero   = prm.get_double("Cohesion");
    phi      = prm.get_double("FrictionAngle")*M_PI/180.;
    psi      = prm.get_double("DilationAngle")*M_PI/180.;
    k1       = prm.get_double("CohesionSoftening");
    shape    = prm.get_double("ShapeFactor");
    corners  = prm.get("Corners");
    verbose  = prm.get_bool("Verbose");
  prm.leave_subsection();

  // Elastic Tangent Operator
  mu     = convert_moduli::bulk_poisson::mu(K,nu);
  lambda = convert_moduli::bulk_poisson::lambda(K,nu);

  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  // Yield & Potential Function Parameters
  //   Yield function     F = root23*Q-(AF-BF*P)
  //   Dilatancy function G = root23*Q-(AG-BG*P)
  cosphi = std::cos(phi);
  cospsi = std::cos(psi);
  sinphi = std::sin(phi);
  sinpsi = std::sin(psi);

  if (corners == "compression")
  {
    AF = 2*root6*c_zero*cosphi/(3-sinphi);
    BF = 2*root6*sinphi/(3-sinphi);
    AG = 2*root6*c_zero*cospsi/(3-sinpsi);
    BG = 2*root6*sinpsi/(3-sinpsi);
  }
  else if (corners == "tension")
  {
    AF = 2*root6*c_zero*cosphi/(3+sinphi);
    BF = 2*root6*sinphi/(3+sinphi);
    AG = 2*root6*c_zero*cospsi/(3+sinpsi);
    BG = 2*root6*sinpsi/(3+sinpsi);
  }
  AF0 = AF;

  // Plastic Internal Variable / Strain
  new_multiplier = 0;
  old_multiplier = 0;

  new_plastic_strain = 0;
  old_plastic_strain = 0;

  new_equiv_plastic_strain = 0;
  old_equiv_plastic_strain = 0;
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress         = in_situ_stress;
  old_elastic_strain = invert(elastic_cto)*in_situ_stress;
  new_elastic_strain = old_elastic_strain;
  new_strain         = new_elastic_strain;
  new_plastic_strain = 0;
  new_multiplier     = 0;
}


// -----------------------------------------------------------------------------
// UPDATE
// -----------------------------------------------------------------------------
template <>
void DruckerPrager<3>::update(const SymmetricTensor<2,3> &incr_strain)
{
  run_3D_update(incr_strain);
}
template <>
void DruckerPrager<2>::update(const SymmetricTensor<2,2> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  run_3D_update(incr_strain_2d);
}


// -----------------------------------------------------------------------------
// ACTUAL 3D STRESS UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain)
{
  // Compute a trial state in which the increment is assumed to be fully elastic
  new_elastic_strain       = old_elastic_strain + incr_strain;
  new_cto                  = elastic_cto;
  new_stress               = new_cto*new_elastic_strain;
  new_plastic_strain       = old_plastic_strain;
  new_equiv_plastic_strain = old_equiv_plastic_strain;
  new_multiplier           = old_multiplier;

  // Return Mapping in Invariant Space
  //
  // Decompose into P and Q invariants
  // and assign initial state vectors
  // slot 0 = P (mean effective stress)
  // slot 1 = Q (deviatoric stress)
  // slot 2 = delta lambda (plastic multiplier)
  // See Sect. 4.13 in Plasticisty Modeling and Computation
  double& P       = xcurrent(0);
  double& Q       = xcurrent(1);
  double& dLambda = xcurrent(2);

  // Calculate stress invariants
  P = one_third*trace(new_stress);
  SymmetricTensor<2,3> n_hat = new_stress - P*eye; // s (deviatoric)
  double s_norm = n_hat.norm();
  Q = root32*s_norm;
  if (s_norm > 1e-15) n_hat /= s_norm;
  else                n_hat = 0;

  // Compute Yield Stress
  //
  // Call assembly routine in "quick mode," just to get yield value.
  // Residual is organized as:
  // slot 0 = P consistency eqn
  // slot 1 = Q consistency eqn
  // slot 2 = yield F value
  assemble_residual_and_jacobian(true);

  // If Elastic (yield F <= 0, elastic response)
  // -------------------------------------------
  if (residual(2) < tol)
  {
    new_strain = new_elastic_strain;
    return;
  }

  // Else Plastic
  // Need to perform plastic iterations using Newton algorithm
  // ---------------------------------------------------------

  // Reset plastic multiplier
  dLambda = 0;

  // Save trial stresses (xtrial: not update during Newton iteration)
  xtrial = xcurrent;

  double resid     = residual.l2_norm();
  double resid_ini = 1 + resid;

  ConvergenceTable convergence_table;

  for (unsigned iter=0; iter<15; ++iter)
  {
    assemble_residual_and_jacobian(false);

    resid = residual.l2_norm();

    if (verbose)
    {
      convergence_table.add_value("Iter",iter);
      convergence_table.add_value("RHS Norm",resid);
      convergence_table.add_value("Reduction",resid/resid_ini);
    }

    if (resid/resid_ini < tol || resid < tol)
      break;

    // Newton's method: xcurrent = xcurrent - jacobian_inv*residual
    residual *= -1;
    jacobian_inv.invert(jacobian);
    jacobian_inv.vmult_add(xcurrent,residual);
  } // end of return mapping

  if (verbose)
  {
    convergence_table.set_scientific("RHS Norm",true);
    convergence_table.set_scientific("Reduction",true);
    std::cout << "Convergence during return mapping of " << this->name() << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    convergence_table.write_text(std::cout,TableHandler::TextOutputFormat::org_mode_table);
  }

  // Upon convergence, reconstruct stresses, strains, cto, etc.
  // ----------------------------------------------------------

  // Reconstruct new stress using the converged state
  new_stress = P*eye + root23*Q*n_hat;

  // Reconstruct new elastic strain
  SymmetricTensor<2,3> dFdsigma = one_third*deriv(0,0)*eye + root32*deriv(0,1)*n_hat;
  SymmetricTensor<2,3> dGdsigma = one_third*deriv(1,0)*eye + root32*deriv(1,1)*n_hat;

  new_elastic_strain       -= dLambda*dGdsigma;
  new_plastic_strain       += dLambda*dGdsigma;
  new_strain                = new_elastic_strain + new_plastic_strain;
  new_equiv_plastic_strain += root23*dLambda*dGdsigma.norm();
  new_multiplier           += dLambda;

  // Reconstruct consistent tangent operator
  SymmetricTensor<2,3> A1 = elastic_cto*dGdsigma;
  SymmetricTensor<2,3> A2 = dFdsigma*elastic_cto;

  new_cto -= outer_product(A1,A2) / (dFdsigma*A1);
}


// -----------------------------------------------------------------------------
// ASSEMBLE RESIDUAL AND JACOBIAN
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::assemble_residual_and_jacobian(const bool quick_check)
{
  const double& P = xcurrent(0);
  const double& Q = xcurrent(1);

  // Yield surface parameters (second term: hyperbolic smoothening)
  double omega_F = two_thirds*Q*Q + shape*shape*AF*AF;

  residual(0) = 0;
  residual(1) = 0;
  residual(2) = sqrt(omega_F) - (AF - BF*P); // yield function

  // If we just need yield value, return now
  if (quick_check)
    return;

  // Else, assemble full Newton system

  // Trial stresses
  const double& dLambda  = xcurrent(2);
  const double& Ptr      = xtrial(0);
  const double& Qtr      = xtrial(1);

  // Hardening/softening and update yield surface parameters
  double cumulative_multiplier = old_multiplier + dLambda;
  cohesion = c_zero * exp( -std::pow(cumulative_multiplier*k1,2) );

  if (corners == "compression")
  {
    AF = 2*root6*cohesion*cosphi/(3-sinphi);
    AG = 2*root6*cohesion*cospsi/(3-sinpsi);
  }
  else if (corners == "tension")
  {
    AF = 2*root6*cohesion*cosphi/(3+sinphi);
    AG = 2*root6*cohesion*cospsi/(3+sinpsi);
  }

               omega_F = two_thirds*Q*Q + shape*shape*AF*AF;
  const double omega_G = two_thirds*Q*Q + shape*shape*AG*AG;

  // Derivatives of F and G w.r.t. invariants (P and Q)
  deriv(0,0) = BF;  // dF/dP
  deriv(1,0) = BG;  // dG/dP
  deriv(0,1) = two_thirds*Q/sqrt(omega_F); // dF/dQ
  deriv(1,1) = two_thirds*Q/sqrt(omega_G); // dG/dQ

  // Residual
  residual(0) = P - Ptr + K*dLambda*deriv(1,0);
  residual(1) = Q - Qtr + 3*mu*dLambda*deriv(1,1);
  residual(2) = sqrt(two_thirds*Q*Q + shape*shape*AF*AF) - (AF-BF*P);

  // Jacobian matrix for local return mapping
  // (0: P; 1: Q; 2: dLambda)
  const double dGdQQ = two_thirds/sqrt(omega_G) - four_ninths*Q*Q*pow(omega_G,-1.5);

  jacobian(0,0) = 1;
  jacobian(0,1) = 0;
  jacobian(0,2) = K*deriv(1,0);

  jacobian(1,0) = 0;
  jacobian(1,1) = 1 + 3*mu*dLambda*dGdQQ;
  jacobian(1,2) = 3*mu*deriv(1,1);

  jacobian(2,0) = deriv(0,0);
  jacobian(2,1) = deriv(0,1);
  jacobian(2,2) = (1-shape)*2*AF0*k1*k1*cumulative_multiplier
                   * exp( -std::pow( cumulative_multiplier*k1,2) );
}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------

// Save state
// -----------------------------------------------------------------------------
template <int dim>
void DruckerPrager<dim>::save_state()
{
  old_stress               = new_stress;
  old_elastic_strain       = new_elastic_strain;
  old_plastic_strain       = new_plastic_strain;
  old_strain               = new_strain;
  old_equiv_plastic_strain = new_equiv_plastic_strain;
  old_multiplier           = new_multiplier;
}


// Stresses
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> DruckerPrager<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> DruckerPrager<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];
  return new_stress_2d;
}
template <int dim>
SymmetricTensor<2,3> DruckerPrager<dim>::full_stress()
{
  return new_stress;
}


// Strains
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> DruckerPrager<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> DruckerPrager<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];
  return new_strain_2d;
}
template <int dim>
SymmetricTensor<2,3> DruckerPrager<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double DruckerPrager<dim>::equiv_plastic_strain(){ return new_equiv_plastic_strain; }

template <int dim>
double DruckerPrager<dim>::vol_strain(){ return trace(new_plastic_strain); }

template <int dim>
double DruckerPrager<dim>::dev_strain()
{
  auto dev_strain = new_plastic_strain - one_third*trace(new_plastic_strain)*eye;
  return dev_strain.norm();
}


// Moduli
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<4,3> DruckerPrager<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> DruckerPrager<2>::cto()
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
double DruckerPrager<dim>::bulk_mod(){ return K; }

template <int dim>
double DruckerPrager<dim>::shear_mod(){ return mu; }

template <int dim>
double DruckerPrager<dim>::poisson(){ return nu; }


// Explicit Instantiation
template class DruckerPrager<2>;
template class DruckerPrager<3>;
