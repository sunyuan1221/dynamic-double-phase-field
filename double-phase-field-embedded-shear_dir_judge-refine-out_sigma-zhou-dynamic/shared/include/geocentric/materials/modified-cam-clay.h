#ifndef MODIFIED_CAM_CLAY_H
#define MODIFIED_CAM_CLAY_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class ModifiedCamClay : public ConstitutiveModel<dim>
{
  public:
    ModifiedCamClay();
    ~ModifiedCamClay();

    void declare_parameters(ParameterHandler &prm);
    void initialize(ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_str);

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

    SymmetricTensor<4,dim> cto();
    double                 bulk_mod();
    double                 shear_mod();
    double                 poisson();

    double                 preconsolidation_pressure();
    double                 OCR();

    void                   set_OCR(double target_OCR);

  private:
    // Input Parameters
    // hyperelasticity
    double P_zero;        // hyperelastic reference pressure (<0)
    double mu_zero;       // hyperelastic reference shear modulus
    double alpha;         // hyperelastic shear parameter
    // cam-clay plasticity
    double new_Pc,old_Pc; // preconsolidation stress
    double Cr;            // hyperelastic compressibility kappa-tilde
    double Cc;            // compressibility parameter lambda-tilde
    double M;             // slope of CSL
    double beta;          // non-associateive flow rule parameter

    bool verbose;

    // Computed Parameters
    double K;             // D_11
    double mu;            // D_22/3
    double coupling_mod;  // D_12 = D_21
    double nu;            // Poisson's ratio (calculated)
    double eps_v_zero;    // hyperelastic reference volumetric strain (<0);
    double P,Q;           // stress invariants
    FullMatrix<double> D,B,HD,DB;
    double Pc_deriv_epsv, Pc_deriv_epsvtr;

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

    // Plastic Multiplier
    double new_multiplier;
    double old_multiplier;

    // Local Newton System
    Vector<double>     xtrial,xcurrent,residual;
    FullMatrix<double> jacobian,jacobian_inv,deriv;

    // Update
    void run_3D_update(const SymmetricTensor<2,3> &incr_strain);
    void assemble_residual_and_jacobian(const bool quick_check);
};


#endif

