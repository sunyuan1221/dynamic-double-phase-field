#include "geocentric/materials/linear-elastic.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
LinearElastic<dim>::LinearElastic()
{ this->name_ = "Linear Elastic"; }

template <int dim>
LinearElastic<dim>::~LinearElastic()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void LinearElastic<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("LinearElastic");
    prm.declare_entry("BulkMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FiniteStrain","false",Patterns::Bool());
    prm.declare_entry("NeoHookean","false",Patterns::Bool());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void LinearElastic<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("LinearElastic");
    K  = prm.get_double("BulkMod");
    nu = prm.get_double("Poisson");
    finite_strain = prm.get_bool("FiniteStrain");
    neo_hookean   = prm.get_bool("NeoHookean");
  prm.leave_subsection();

  mu     = convert_moduli::bulk_poisson::mu     (K,nu);
  lambda = convert_moduli::bulk_poisson::lambda (K,nu);
  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  if (finite_strain)
  {
    new_b = eye;
    new_J = 1.;
  }

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void LinearElastic<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_strain = invert(elastic_cto)*in_situ_stress;
  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <int dim>
void LinearElastic<dim>::update(const SymmetricTensor<2,dim> &incr_strain)
{
  SymmetricTensor<2,3> incr_strain_3d;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=i; j<dim; ++j)
    incr_strain_3d[i][j] = incr_strain[i][j];

  new_strain = old_strain + incr_strain_3d;
  new_stress = elastic_cto*new_strain;
}

template <int dim>
void LinearElastic<dim>::update(const Tensor<2,dim> &rel_def_grad)
{
  Tensor<2,3> rel_def_grad_3d;

  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
    rel_def_grad_3d[i][j] = rel_def_grad[i][j];

  if (dim == 2)
    rel_def_grad_3d[2][2] = 1;

  new_b = rel_def_grad_3d * old_b * transpose(rel_def_grad_3d);

  Tensor<1,3> principal_stretches_square;
  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(new_b, principal_stretches_square, principal_directions);

  Tensor<2,3> m_aa;
  new_strain = 0;
  for (unsigned a=0; a<3; ++a)
  {
    m_aa = outer_product(principal_directions[a],principal_directions[a]);
    new_strain += 0.5*log(principal_stretches_square[a])*symmetrize(m_aa);
  }

  double eps_v = trace(new_strain);
  new_J = exp(eps_v);

  // Update stresses
  SymmetricTensor<4,3> new_cto;
  if (neo_hookean)
  {
    double new_lambda = lambda*eps_v;
    new_stress = (1/new_J) * ( mu*symmetrize(new_b) - (mu - new_lambda)*eye );

    new_cto = elastic_cto - 2*new_lambda*big_eye;
  }
  else // hencky elasticity
  {
    new_stress = (1/new_J) * elastic_cto*new_strain;
    new_cto = elastic_cto;
  }
}


// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void LinearElastic<dim>::save_state()
{
  old_stress = new_stress;
  old_strain = new_strain;

  if (finite_strain)
  {
    old_b   = new_b;
    old_J   = new_J;
  }
}

// Stress
template <>
SymmetricTensor<2,3> LinearElastic<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> LinearElastic<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

// Strain
template <>
SymmetricTensor<2,3> LinearElastic<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> LinearElastic<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> LinearElastic<3>::cto()
{
  return elastic_cto;
}
template <>
SymmetricTensor<4,2> LinearElastic<2>::cto()
{
  SymmetricTensor<4,2> elastic_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    elastic_cto_2d[i][j][k][l] = elastic_cto[i][j][k][l];

  return elastic_cto_2d;
}

// Bulk modulus
template <int dim>
double LinearElastic<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double LinearElastic<dim>::shear_mod(){ return mu; }

// Equivalent plastic strain
template <int dim>
double LinearElastic<dim>::equiv_plastic_strain(){ return 0.; }

template <int dim>
double LinearElastic<dim>::old_jacobian(){ return old_J; }


// EXPLICIT INSTANTIATION
template class LinearElastic<2>;
template class LinearElastic<3>;
