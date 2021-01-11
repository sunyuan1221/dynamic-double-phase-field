#ifndef THREEINVARIANT_CAM_CLAY_H
#define THREEINVARIANT_CAM_CLAY_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class ThreeInvariantCamClay : public ConstitutiveModel<dim>
{
  public:
    ThreeInvariantCamClay();
    ~ThreeInvariantCamClay();

    void declare_parameters(ParameterHandler &prm);
    void initialize(ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update(const Tensor<2,dim> &rel_def_grad);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();
    double                 vol_stress();
    double                 dev_stress();
    double                 stress_ratio();

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
    double                 poisson();

    double                 preconsolidation_pressure();
    double                 OCR();

    void                   set_OCR(double target_OCR);

  private:
    // Input Parameters
    //   Hyperelasticity
    double P_zero;        // hyperelastic reference pressure (<0)
    double mu_zero;       // hyperelastic reference shear modulus
    double alpha;         // hyperelastic shear parameter
    //   Cam-clay plasticity
    double new_Pc,old_Pc; // preconsolidation stress
    double Cr;            // hyperelastic compressibility kappa-tilde
    double Cc;            // compressibility parameter lambda-tilde
    double M;             // slope of CSL
    double beta;          // non-associateive flow rule parameter

    bool verbose;

    bool finite_strain;

    // Lode's angle enhancement
    GudehusArgyris scaling_function;

    // Computed Parameters
    double K;             // D_11
    double mu;            // D_22/3
    double coupling_mod;  // D_12 = D_21
    double nu;            // Poisson's ratio (calculated)
    double eps_v_zero;    // hyperelastic reference volumetric strain (<0);
    double P,Q,Zeta;      // stress invariants
    double eps_v,eps_s;   // strain invariants

    Tensor<1,3> principal_stresses;
    Tensor<1,3> principal_e_strains;

    double                      F(const double Pc, const double beta=1.0);
    Tensor<1,3>            grad_F(const double Pc, const double beta=1.0);
    SymmetricTensor<2,3>   hess_F(const double beta=1.0);
    Tensor<1,3>          grad_Psi(const Tensor<1,3> &eps);
    SymmetricTensor<2,3> hess_Psi(const Tensor<1,3> &eps);
    double               Pc_deriv_eps, Pc_deriv_epstr;

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

    // Plastic Multiplier
    double new_multiplier;
    double old_multiplier;

    // Additional variables for finite strain
    Tensor<2,3> new_elastic_b; // elastic part of left Cauchy-Green tensor
    Tensor<2,3> old_elastic_b;
    double new_J;
    double old_J;
    // Tensor<4,3> new_cto_kirchhoff;

    // Update
    void run_3D_update(Tensor<1,3> &pricipal_e_strains,
                       const std::vector<Tensor<1,3> > &principal_directions);
};


#endif

