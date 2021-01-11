#include <geocentric/materials/threeinvariant-cam-clay.h>

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
ThreeInvariantCamClay<dim>::ThreeInvariantCamClay()
{ this->name_ = "three-invariant Cam-Clay"; }

template <int dim>
ThreeInvariantCamClay<dim>::~ThreeInvariantCamClay()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void ThreeInvariantCamClay<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("ThreeInvariantCamClay");
    // Hyperelasticity
    prm.declare_entry("ReferencePressure","1",Patterns::Double());
    prm.declare_entry("ReferencePoisson","0.3",Patterns::Double());
    prm.declare_entry("ShearModulusEvolution","0.0",Patterns::Double());
    // Cam-Clay plasticity
    prm.declare_entry("PreconsolidationPressure","1",Patterns::Double());
    prm.declare_entry("VirginCompressionIndex","0.1",Patterns::Double());
    prm.declare_entry("RecompressionIndex","0.01",Patterns::Double());
    prm.declare_entry("CSLSlope","1.0",Patterns::Double());
    prm.declare_entry("Associativity","1.0",Patterns::Double());
    prm.declare_entry("Verbose","false",Patterns::Bool());
    // Lode's angle enhancement
    prm.declare_entry("Ellipticity","1.0",Patterns::Double());
    // Flag for finite strain
    prm.declare_entry("FiniteStrain","false",Patterns::Bool());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// INITIALIZE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void ThreeInvariantCamClay<dim>::initialize(ParameterHandler &prm)
{
  // Parse Parameters
  prm.enter_subsection("ThreeInvariantCamClay");
    // hyperelasticity
    P_zero   = prm.get_double("ReferencePressure")*(-1);
    nu       = prm.get_double("ReferencePoisson");
    alpha    = prm.get_double("ShearModulusEvolution");
    // cam-clay plasticity
    new_Pc   = prm.get_double("PreconsolidationPressure")*(-1);
    Cc       = prm.get_double("VirginCompressionIndex");
    Cr       = prm.get_double("RecompressionIndex");
    M        = prm.get_double("CSLSlope");
    beta     = prm.get_double("Associativity");
    verbose  = prm.get_bool("Verbose");
    // Lode's angle enhancement
    scaling_function.set_ellipticity(prm.get_double("Ellipticity"));
    // Finite strain
    finite_strain = prm.get_bool("FiniteStrain");
  prm.leave_subsection();

  // Hyperelasticity parameters
  P          = P_zero;
  Q          = 0;
  K          = -P_zero/Cr;
  mu         = convert_moduli::bulk_poisson::mu(K,nu);
  mu_zero    = mu + alpha*P_zero;
  new_stress = P_zero*eye;
  eps_v_zero = 0;

  // Plastic internal variable / strain
  new_plastic_strain = 0;
  new_equiv_plastic_strain = 0;
  new_multiplier = 0;

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
void ThreeInvariantCamClay<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  if (old_stress == in_situ_stress)
    return;

  new_stress = in_situ_stress;

  P = one_third*trace(new_stress);
  Q = root32*(new_stress - P*eye).norm();

  K  = -P/Cr;
  mu = mu_zero - alpha*P;
  nu = (3*K-2*mu) / (2*(3*K+mu));

  // Reset strain (note: bulk and shear moduli are assumed to be uncoupled, i.e. alpha = 0)
  eps_v = -Cr*log(P/P_zero);
  eps_s = Q/(3*mu);

  SymmetricTensor<2,3> n_hat = new_stress - P*eye;
  double e_norm = n_hat.norm();
  if (e_norm > 1e-15) n_hat /= e_norm;
  else                n_hat = 0;

  new_elastic_strain       = one_third*eps_v*eye + root32*eps_s*n_hat;
  new_plastic_strain       = 0;
  new_strain               = new_elastic_strain;
  new_equiv_plastic_strain = 0;
  new_multiplier           = 0;

  if (finite_strain)
  {
    new_J = exp(eps_v);

    std::vector<Tensor<1,3> > principal_directions(3);
    eigen_decompose<3>(new_elastic_strain, principal_e_strains, principal_directions);

    Tensor<2,3> m_aa;
    new_elastic_b = 0;
    for (unsigned a=0; a<3; ++a)
    {
      m_aa = outer_product(principal_directions[a],principal_directions[a]);
      new_elastic_b += pow(exp(principal_e_strains[a]),2)*m_aa;
    }
  }

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE
// -----------------------------------------------------------------------------
// infinitesimal strain
template <int dim>
void ThreeInvariantCamClay<dim>::update(const SymmetricTensor<2,dim> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_3d;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=i; j<dim; ++j)
    incr_strain_3d[i][j] = incr_strain[i][j];

  new_elastic_strain = old_elastic_strain + incr_strain_3d;

  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(new_elastic_strain, principal_e_strains, principal_directions);

  run_3D_update(principal_e_strains, principal_directions);
}

// finite strain
template <int dim>
void ThreeInvariantCamClay<dim>::update(const Tensor<2,dim> &rel_def_grad)
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

  principal_e_strains[0] = 0.5*log(principal_stretches_square[0]);
  principal_e_strains[1] = 0.5*log(principal_stretches_square[1]);
  principal_e_strains[2] = 0.5*log(principal_stretches_square[2]);

  run_3D_update(principal_e_strains, principal_directions);

  new_stress = (1/new_J) * new_stress; // since logarithmic strains gives kirchhoff stress

  new_elastic_b = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_elastic_b += exp(2*principal_e_strains[a])*m_aa;
  }
}


// -----------------------------------------------------------------------------
// ACTUAL 3D UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void ThreeInvariantCamClay<dim>::run_3D_update(Tensor<1,3> &principal_e_strains,
                                               const std::vector<Tensor<1,3> > &principal_directions)
{
  new_plastic_strain       = old_plastic_strain;
  new_strain               = new_elastic_strain + new_plastic_strain;
  new_equiv_plastic_strain = old_equiv_plastic_strain;
  new_multiplier           = old_multiplier;
  new_Pc                   = old_Pc;

  if (finite_strain)
    new_J = exp(trace(new_strain));

  // Get principal stresses from derivative of the free energy function
  const Tensor<1,3> trial_principal_e_strains = principal_e_strains;
  principal_stresses = grad_Psi(trial_principal_e_strains);
  scaling_function.update(principal_stresses);
  Zeta = scaling_function.value();

  // Compute Yield Stress
  double yield = F(old_Pc);

  // If Elastic (yield F <= 0, elastic response)
  // -------------------------------------------
  if (yield < tol)
  {
    // Stress
    new_stress = 0;
    Tensor<2,3> m_aa;
    for (unsigned a=0; a<3; ++a)
    {
      m_aa = outer_product(principal_directions[a],principal_directions[a]);
      new_stress += principal_stresses[a]*symmetrize(m_aa);
    }

    // Elastic CTO (Eq. 6.15 in Plasticity Modeling and Computation)
    SymmetricTensor<2,3> n_hat = new_elastic_strain - one_third*eps_v*eye;
    double e_norm = n_hat.norm();
    if (e_norm > 1e-15) n_hat /= e_norm;
    else                n_hat = 0;

    SymmetricTensor<4,3> eye_dyad_nhat = outer_product(eye,n_hat);
    SymmetricTensor<4,3> nhat_dyad_eye = outer_product(n_hat,eye);

    double c1 = 2*mu;
    double c2 = K - two_thirds*mu;
    double c3 = root23*coupling_mod;

    new_cto = c1 * big_eye +
              c2 * eye_dyad_eye +
              c3 * (eye_dyad_nhat + nhat_dyad_eye);

    return;
  }

  // Else Plastic
  // Need to perform plastic iterations using Newton algorithm
  // ---------------------------------------------------------

  // Return Mapping in Principal Axes
  //
  // Eigen (spectral) decompose elastic strains
  // and assign initial state vectors
  // slot 0 = Epsilon 1 consistency (equation) / Epsilon 1 (unknown)
  // slot 1 = Epsilon 2 consistency (equation) / Epsilon 2 (unknown)
  // slot 2 = Epsilon 3 consistency (equation) / Epsilon 3 (unknown)
  // slot 3 = Yield function (equation) / Plastic multiplier increment (unknown)
  // See Sect. 6.9 in Plasticity Modeling and Computation

  double eps_v_tr = 0;
  for (unsigned i=0; i<3; ++i)
    eps_v_tr += trial_principal_e_strains[i];

  Tensor<1,3> gradF, gradG;
  SymmetricTensor<2,3> hessG, hessPsi;

  // Newton iteration
  FullMatrix<double> jacobian(4,4);
  FullMatrix<double> jacobian_inv(4,4);

  Vector<double> xdelta(4);
  Vector<double> residual(4);
                 residual(3) = yield;

  double dLambda = 0;

  double resid = residual.l2_norm();
  double resid_ini = 1 + resid;

  ConvergenceTable convergence_table;

  for (unsigned iter=0; iter<20; ++iter)
  {
    // Hardening
    new_Pc = old_Pc*exp( (eps_v-eps_v_tr)/(Cc-Cr) );

    // Derivatives of F and G w.r.t. principal stresses
    gradF   = grad_F(new_Pc);
    gradG   = grad_F(new_Pc,beta);
    hessG   = hess_F(beta);
    hessPsi = hess_Psi(principal_e_strains);

    // Residual
    for (unsigned i=0; i<3; ++i)
      residual(i) = principal_e_strains[i]
                    - trial_principal_e_strains[i]
                    + dLambda*gradG[i];
    residual(3) = F(new_Pc);

    resid = residual.l2_norm();

    if (verbose)
    {
      convergence_table.add_value("Iter",iter);
      convergence_table.add_value("RHS Norm",resid);
      convergence_table.add_value("Reduction",resid/resid_ini);
    }

    Tensor<1,3>          gradF_eps;
    SymmetricTensor<2,3> hessG_eps;

    // Jacobian matrix for local return mapping
    for (unsigned a=0; a<3; ++a)
      for (unsigned i=0; i<3; ++i)
      {
        gradF_eps[i] += gradF[a]*hessPsi[a][i];

        for (unsigned j=i; j<3; ++j)
          hessG_eps[i][j] += hessG[i][a]*hessPsi[a][j];
      }

    Pc_deriv_eps = new_Pc/(Cc-Cr);

    for (unsigned i=0; i<3; ++i)
    {
      jacobian(i,3) = gradG[i];
      jacobian(3,i) = gradF_eps[i] - P*Pc_deriv_eps;

      for (unsigned j=0; j<3; ++j)
        jacobian(i,j) = (i==j) + dLambda*(hessG_eps[i][j] - one_third*Pc_deriv_eps);
    }
    jacobian(3,3) = 0;

    if (resid/resid_ini < tol or resid < tol)
      break;

    // Newton update: xnew = xcurrent - xdelta, xdelta = jacobian_inv*residual
    jacobian_inv.invert(jacobian);
    jacobian_inv.vmult(xdelta,residual);

    // Update principal stresses and strains
    principal_e_strains[0] -= xdelta(0);
    principal_e_strains[1] -= xdelta(1);
    principal_e_strains[2] -= xdelta(2);
    dLambda                -= xdelta(3);

    // Stress and scaling function
    principal_stresses = grad_Psi(principal_e_strains);
    scaling_function.update(principal_stresses);
    Zeta = scaling_function.value();
  } // end of return mapping

  if (verbose)
  {
    convergence_table.set_auto_fill_mode(true);
    convergence_table.set_scientific("RHS Norm",true);
    convergence_table.set_scientific("Reduction",true);
    std::cout << "Convergence during local return mapping of " << this->name() << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    convergence_table.write_text(std::cout,TableHandler::TextOutputFormat::org_mode_table);
  }

  // Upon convergence, reconstruct strains, stresses, cto, etc.
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

  // Reconstruct strain
  SymmetricTensor<2,3> dGdsigma;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    dGdsigma += gradG[a]*symmetrize(m_aa);
  }

  new_elastic_strain       -= dLambda*dGdsigma;
  new_plastic_strain       += dLambda*dGdsigma;
  new_strain                = new_elastic_strain + new_plastic_strain;
  new_equiv_plastic_strain += root23*dLambda*dGdsigma.norm();
  new_multiplier           += dLambda;

  // Reconstruct consistent tangent opterator (see Sect. 6.9)
  // Note: derivative w.r.t. trial strains
  new_cto = 0;

  Pc_deriv_epstr = -Pc_deriv_eps;

  FullMatrix<double> residual_deriv(4,3);
  for (unsigned i=0; i<4; ++i)
  for (unsigned b=0; b<3; ++b)
  {
    if (i < 3)
      residual_deriv[i][b] = -(i==b) + dLambda*(-one_third)*Pc_deriv_epstr;
    else if (i == 3)
      residual_deriv[i][b] = -P*Pc_deriv_epstr;
  }

  jacobian_inv.invert(jacobian);
  FullMatrix<double> x_deriv(3,3);
  for (unsigned c=0; c<3; ++c)
  for (unsigned b=0; b<3; ++b)
    for (unsigned i=0; i<4; ++i)
      x_deriv[c][b] -= jacobian_inv[c][i]*residual_deriv[i][b];

  FullMatrix<double> elastic_a_trial(3,3);
  for (unsigned a=0; a<3; ++a)
  for (unsigned b=0; b<3; ++b)
    for (unsigned c=0; c<3; ++c)
      elastic_a_trial[a][b] += hessPsi[a][c]*x_deriv[c][b];

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
template <int dim>
double ThreeInvariantCamClay<dim>::F(const double Pc,
                                     const double beta)
{
  // Commented out because we have updated the three invariants before calling them.
  // However, below should be uncommented when used independently.
/*
  Tensor<1,3> &sig = principal_stresses;
  P = one_third * (sig[0] + sig[1] + sig[2]);
  Q = sqrt( 0.5 * (   (sig[0]-sig[1]) * (sig[0]-sig[1])
                    + (sig[1]-sig[2]) * (sig[1]-sig[2])
                    + (sig[0]-sig[2]) * (sig[0]-sig[2]) ) );
  Zeta = scaling_function.value();
*/
  return ( Zeta*Zeta*beta*Q*Q/(M*M) + P*(P-Pc) );
}

// Derivative of yield/potential function w.r.t. principal stresses
template <int dim>
Tensor<1,3> ThreeInvariantCamClay<dim>::grad_F(const double Pc,
                                               const double beta)
{
  double dFdP    = 2*P - Pc;
  double dFdQ    = Zeta*Zeta*beta/(M*M); // 2Q is omitted because gradQ contains 1/2Q
  double dFdZeta = 2*Zeta*beta*Q*Q/(M*M);

  Tensor<1,3> &sig = principal_stresses;
  Tensor<1,3> gradQ;
  gradQ[0] = 2*sig[0] - sig[1] - sig[2];
  gradQ[1] = 2*sig[1] - sig[0] - sig[2];
  gradQ[2] = 2*sig[2] - sig[0] - sig[1];
  Tensor<1,3> gradZeta = scaling_function.grad();

  Tensor<1,3> gradF;
  gradF[0] = dFdP*one_third + dFdQ*gradQ[0] + dFdZeta*gradZeta[0];
  gradF[1] = dFdP*one_third + dFdQ*gradQ[1] + dFdZeta*gradZeta[1];
  gradF[2] = dFdP*one_third + dFdQ*gradQ[2] + dFdZeta*gradZeta[2];

  return gradF;
}

// Second derivative of yield/potential function w.r.t. principal stresses
template <int dim>
SymmetricTensor<2,3> ThreeInvariantCamClay<dim>::hess_F(const double beta)
{
  Tensor<1,3> S;
  S[0] = principal_stresses[0] - P;
  S[1] = principal_stresses[1] - P;
  S[2] = principal_stresses[2] - P;
  double trS2 = S[0]*S[0] + S[1]*S[1] + S[2]*S[2];
  double chi = sqrt(trS2);
  if (fabs(chi) < tol)
    chi = tol;
  double const1 = root32/chi;

  Tensor<1,3> gradQ = const1*S;
  SymmetricTensor<2,3> hessQ;
  for (unsigned a=0; a<3; ++a)
  for (unsigned b=a; b<3; ++b)
    hessQ[a][b] = const1*( (a==b) - one_third - S[a]*S[b]/(chi*chi) );

  Tensor<1,3>          gradZeta = scaling_function.grad();
  SymmetricTensor<2,3> hessZeta = scaling_function.hess();

  Tensor<2,3> gradQ_dyad_gradQ       = outer_product(gradQ,   gradQ);
  Tensor<2,3> gradZeta_dyad_gradZeta = outer_product(gradZeta,gradZeta);
  Tensor<2,3> gradQ_dyad_gradZeta    = outer_product(gradQ,   gradZeta);
  Tensor<2,3> gradZeta_dyad_gradQ    = outer_product(gradZeta,gradQ);

  double dFdQ      = 2*Zeta*Zeta*beta*Q/(M*M);
  double d2FdQ2    = 2*Zeta*Zeta*beta/(M*M);
  double dFdZeta   = 2*Zeta*beta*Q*Q/(M*M);
  double d2FdZeta2 = 2*beta*Q*Q/(M*M);
  double d2FdZetaQ = 4*Zeta*beta*Q/(M*M);

  // Eq. (6.42) in Plasticity Modeling and Computation
  SymmetricTensor<2,3> hessF = dFdQ*hessQ
                               + dFdZeta*hessZeta
                               + d2FdQ2*symmetrize(gradQ_dyad_gradQ)
                               + d2FdZeta2*symmetrize(gradZeta_dyad_gradZeta)
                               + d2FdZetaQ*(symmetrize(gradQ_dyad_gradZeta) +
                                            symmetrize(gradZeta_dyad_gradQ) );

  for (unsigned i=0; i<3; ++i)
  for (unsigned j=i; j<3; ++j)
    hessF[i][j] += two_ninths;

  return hessF;
}


// -----------------------------------------------------------------------------
// DERIVATIES OF FREE ENERGY FUNCTION
// -----------------------------------------------------------------------------
// Derivative of free energy function w.r.t. principal strains
template <int dim>
Tensor<1,3> ThreeInvariantCamClay<dim>::grad_Psi(const Tensor<1,3> &eps)
{
  eps_v = eps[0] + eps[1] + eps[2];
  eps_s = sqrt( two_ninths * (   (eps[0]-eps[1]) * (eps[0]-eps[1])
                               + (eps[1]-eps[2]) * (eps[1]-eps[2])
                               + (eps[0]-eps[2]) * (eps[0]-eps[2]) ) );

  double exp_omega = exp( (eps_v_zero - eps_v)/Cr );
  P                = P_zero * (1 + 3*alpha*eps_s*eps_s/(2*Cr)) * exp_omega;
  K                = -P/Cr;
  mu               = mu_zero - alpha*P_zero*exp_omega;
  nu               = (3*K-2*mu) / (2*(3*K+mu));
  coupling_mod     = 3*P_zero*alpha*eps_s*exp_omega/Cr;
  Q                = 3*mu*eps_s;

  Assert(nu >= 0, ExcMessage("Negative Poisson's ratio produced!"));

  double const1 = P;
  double const2 = two_thirds*mu;
  Tensor<1,3> gradPsi;

  gradPsi[0] = const1 + const2*(2*eps[0] - eps[1] - eps[2]);
  gradPsi[1] = const1 + const2*(2*eps[1] - eps[0] - eps[2]);
  gradPsi[2] = const1 + const2*(2*eps[2] - eps[0] - eps[1]);

  return gradPsi;
}

// Second derivative of free energy function w.r.t. principal strains
template <int dim>
SymmetricTensor<2,3> ThreeInvariantCamClay<dim>::hess_Psi(const Tensor<1,3> &eps)
{
  double exp_omega = exp( (eps_v_zero - eps_v)/Cr );
  double const1 = (-1/Cr) * P_zero * (1 + 3*alpha*eps_s*eps_s/(2*Cr)) * exp_omega;
  double const2 = two_thirds*(mu_zero - alpha*P_zero*exp_omega);
  double const3 = two_thirds*alpha*P_zero*exp_omega/Cr;

  SymmetricTensor<2,3> hessPsi;

  hessPsi[0][0] = const1 + 2*const2 + 2*const3*( (2*eps[0] - eps[1] - eps[2]) );
  hessPsi[0][1] = const1 -   const2 +   const3*( (2*eps[0] - eps[1] - eps[2]) + (2*eps[1] - eps[0] - eps[2]) );
  hessPsi[0][2] = const1 -   const2 +   const3*( (2*eps[0] - eps[1] - eps[2]) + (2*eps[2] - eps[0] - eps[1]) );
  hessPsi[1][1] = const1 + 2*const2 + 2*const3*( (2*eps[1] - eps[0] - eps[2]) );
  hessPsi[1][2] = const1 -   const2 +   const3*( (2*eps[1] - eps[0] - eps[2]) + (2*eps[2] - eps[0] - eps[1]) );
  hessPsi[2][2] = const1 + 2*const2 + 2*const3*( (2*eps[2] - eps[0] - eps[1]) );

  return hessPsi;
}


// -----------------------------------------------------------------------------
// SET TO CERTAIN OCR
// -----------------------------------------------------------------------------
template <int dim>
void ThreeInvariantCamClay<dim>::set_OCR(double target_OCR)
{
  if (fabs(P) < 1e-15) new_Pc = target_OCR * (Zeta*Zeta*Q*Q/(M*M*(-1e-15) - 1e-15));
  else                 new_Pc = target_OCR * (Zeta*Zeta*Q*Q/(M*M*P) + P);
  old_Pc = new_Pc;
}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------

// Save state
// -----------------------------------------------------------------------------
template <int dim>
void ThreeInvariantCamClay<dim>::save_state()
{
  old_stress               = new_stress;
  old_elastic_strain       = new_elastic_strain;
  old_plastic_strain       = new_plastic_strain;
  old_strain               = new_strain;
  old_equiv_plastic_strain = new_equiv_plastic_strain;
  old_multiplier           = new_multiplier;
  old_Pc                   = new_Pc;

  if (finite_strain)
  {
    old_elastic_b = new_elastic_b;
    old_J         = new_J;
  }
}


// Stresses
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> ThreeInvariantCamClay<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> ThreeInvariantCamClay<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> ThreeInvariantCamClay<dim>::full_stress()
{
  return new_stress;
}

template <int dim>
double ThreeInvariantCamClay<dim>::vol_stress(){ return P; }

template <int dim>
double ThreeInvariantCamClay<dim>::dev_stress(){ return Q; }

template <int dim>
double ThreeInvariantCamClay<dim>::stress_ratio(){ return Q/P; }


// Strains
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<2,3> ThreeInvariantCamClay<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> ThreeInvariantCamClay<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];
  return new_strain_2d;
}
template <int dim>
SymmetricTensor<2,3> ThreeInvariantCamClay<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double ThreeInvariantCamClay<dim>::equiv_plastic_strain(){ return new_equiv_plastic_strain; }

template <int dim>
double ThreeInvariantCamClay<dim>::vol_strain(){ return trace(new_strain); }

template <int dim>
double ThreeInvariantCamClay<dim>::dev_strain()
{
  SymmetricTensor<2,3> dev_strain = new_strain - one_third*trace(new_strain)*eye;
  return root23*dev_strain.norm();
}

template <int dim>
double ThreeInvariantCamClay<dim>::old_jacobian(){ return old_J; }


// Moduli
// -----------------------------------------------------------------------------
template <>
SymmetricTensor<4,3> ThreeInvariantCamClay<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> ThreeInvariantCamClay<2>::cto()
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
Tensor<4,dim> ThreeInvariantCamClay<dim>::cto_kirchhoff()
{
  Tensor<4,dim> cto_kirchhoff;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
  for (unsigned k=0; k<dim; ++k)
  for (unsigned l=0; l<dim; ++l)
    cto_kirchhoff[i][j][k][l] = new_cto[i][j][k][l] - (j==k)*new_J*new_stress[i][l];

  return cto_kirchhoff;
}

template <int dim>
double ThreeInvariantCamClay<dim>::bulk_mod(){ return K; }

template <int dim>
double ThreeInvariantCamClay<dim>::shear_mod(){ return mu; }

template <int dim>
double ThreeInvariantCamClay<dim>::poisson(){ return nu; }


// History variables
// -----------------------------------------------------------------------------

// Preconsolidation stress
template <int dim>
double ThreeInvariantCamClay<dim>::preconsolidation_pressure(){ return -new_Pc; }

// OCR
template <int dim>
double ThreeInvariantCamClay<dim>::OCR()
{
  if (abs(P) < 1e-15) return new_Pc / (Zeta*Zeta*Q*Q/(M*M*(-1e-15) - 1e-15));
  else                return new_Pc / (Zeta*Zeta*Q*Q/(M*M*P) + P);
}


// Explicit Instantiation
template class ThreeInvariantCamClay<2>;
template class ThreeInvariantCamClay<3>;


