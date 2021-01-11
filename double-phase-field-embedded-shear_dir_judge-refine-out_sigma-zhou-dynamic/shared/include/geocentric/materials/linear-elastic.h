#ifndef LINEAR_ELASTIC_H
#define LINEAR_ELASTIC_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class LinearElastic : public ConstitutiveModel<dim>
{
  public:
    LinearElastic();
    ~LinearElastic();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update(const Tensor<2,dim> &rel_def_grad);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,dim> strain();
    SymmetricTensor<4,dim> cto();
    double                 bulk_mod();
    double                 shear_mod();
    double                 equiv_plastic_strain();
    double                 old_jacobian();

  private:
    // Input Parameters
    double K;
    double nu;
    double lambda;
    double mu;

    bool finite_strain;
    bool neo_hookean;

    // History variable
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;

    SymmetricTensor<4,3> elastic_cto;

    // Additional variables for finite strain
    Tensor<2,3> new_b; // left Cauchy-Green tensor
    Tensor<2,3> old_b;
    double new_J;
    double old_J;

    Tensor<4,3> new_cto_kirchhoff;

    void run_3D_update (const SymmetricTensor<2,3> &incr_strain);
    void run_3D_update (const Tensor<2,3> &ref_def_grad);
};


#endif

