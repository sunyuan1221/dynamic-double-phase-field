#ifndef MOHR_COULOMB_H
#define MOHR_COULOMB_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class MohrCoulomb : public ConstitutiveModel<dim>
{
  public:
    MohrCoulomb();
    ~MohrCoulomb();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress (const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update(const Tensor<2,dim> &rel_def_grad);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();

    SymmetricTensor<2,dim> strain();
    SymmetricTensor<2,3>   full_strain();
    double                 vol_strain();
    double                 dev_strain();
    double                 equiv_plastic_strain();
    double                 old_jacobian();

    SymmetricTensor<4,dim> cto();
    Tensor<4,dim>          cto_kirchhoff();
    double                 bulk_mod();
    double                 shear_mod();

  private:
    // Input parameters
    double K;
    double nu;
    double cohesion;
    double phi;
    double psi;
    double shape;

    bool verbose;

    bool finite_strain;

    // Computed Parameters
    //   Elastic moduli
    double lambda;
    double mu;
    double youngs;
    //   Yield and potential surfaces parameters
    double cosphi;
    double cospsi;
    double sinphi;
    double sinpsi;
    double shape_phi;
    double shape_psi;

    Tensor<1,3> principal_stresses;
    Tensor<1,3> principal_e_strains;

    double F           (const Tensor<1,3> &sigma, bool potential = false);
    Tensor<1,3> grad_F (const Tensor<1,3> &sigma, bool potential = false);

    // Stresses, Strains, Moduli
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;
    SymmetricTensor<2,3> new_elastic_strain;
    SymmetricTensor<2,3> old_elastic_strain;
    SymmetricTensor<2,3> new_plastic_strain;
    SymmetricTensor<2,3> old_plastic_strain;
    double new_equiv_plastic_strain;
    double old_equiv_plastic_strain;

    SymmetricTensor<4,3> new_cto;
    SymmetricTensor<4,3> elastic_cto;

    Tensor<2,3> elastic_a; // tangent for principal values
    Tensor<2,3> elastic_a_inv;

    // Plastic multiplier
    double new_multiplier;
    double old_multiplier;

    // Additional variables for finite strain
    Tensor<2,3> new_elastic_b; // elastic part of left Cauchy-Green tensor
    Tensor<2,3> old_elastic_b;
    double new_J;
    double old_J;

    // Update
    void run_3D_update();
};

#endif

