#ifndef PHASEFIELD_COHESIVE_H
#define PHASEFIELD_COHESIVE_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class PhaseFieldCohesive : public ConstitutiveModel<dim>
{
  public:
    PhaseFieldCohesive();
    ~PhaseFieldCohesive();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update_phasefield(const double pf, const double pf_extra = 0);
    void update_pressure_contribution(const double pressure_term);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();

    SymmetricTensor<2,dim> strain();
    SymmetricTensor<2,3>   full_strain();
    double                 equiv_plastic_strain();

    SymmetricTensor<4,dim> cto();
    SymmetricTensor<2,dim> cto_pf();

    double                 bulk_mod();
    double                 shear_mod();

    double                 fracture_energy();
    double                 length_scale();
    double                 critical_stress();

    double                 crack_driving_force();
    double                 crack_driving_force_deriv();
    SymmetricTensor<2,dim> crack_driving_force_deriv_strain();
    double                 old_crack_driving_force();
    double                 H_plus_output();

    double                 strain_energy_plus();
    double                 strain_energy_minus();

    double                 degrade();
    double                 degrade_extra();
    double                 degrade_deriv();
    double                 degrade_deriv_2();
    double                 degrade(const double pf);
    double                 degrade_deriv(const double pf);
    double                 degrade_deriv_2(const double pf);

    double                 lagrange_multiplier();
    void                   set_lagrange_multiplier(double new_lagrange_mult);

    void                   set_fracture_energy(double new_G_c);
    void                   set_youngs_mod(double new_E);
    void                   set_crack_driving_force(double new_crack_driving_force);

  private:
    // Input Parameters
    double E;
    double K;
    double nu;
    double G_c;
    double L;
    double sig_c;
    double psi_c;
    double exp_p;
    bool   cohesive_fracture;

    // Computed Parameters
    //   Elastic moduli
    double lambda;
    double mu;
    double K_damaged;
    double mu_damaged;

    // History variable
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;

    SymmetricTensor<4,3> new_cto;
    SymmetricTensor<4,3> elastic_cto;

    SymmetricTensor<2,3> new_cto_pf;

    double g_d;
    double g_d_deriv;
    double g_d_deriv_2;
    double g_d_extra;
    double M;
    double new_H_plus, old_H_plus; // irreversibility parameter
    double pressure_contribution;
    double lagrange_mult;

    SymmetricTensor<2,3> new_H_plus_deriv_strain;

    bool initial_crack_point;

    // Update
    void run_3D_update (const SymmetricTensor<2,3> &incr_strain);
};


#endif
