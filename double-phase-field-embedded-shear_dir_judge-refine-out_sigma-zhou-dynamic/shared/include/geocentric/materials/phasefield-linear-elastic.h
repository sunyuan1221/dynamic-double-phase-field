#ifndef PHASEFIELD_LINEAR_ELASTIC_H
#define PHASEFIELD_LINEAR_ELASTIC_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class PhaseFieldLinearElastic : public ConstitutiveModel<dim>
{
  public:
    PhaseFieldLinearElastic();
    ~PhaseFieldLinearElastic();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update_phasefield(const double pf);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();

    SymmetricTensor<2,dim> strain();
    SymmetricTensor<2,3>   full_strain();
    double                 equiv_plastic_strain();

    SymmetricTensor<4,dim> cto();

    double                 bulk_mod();
    double                 shear_mod();

    double                 fracture_energy();
    double                 length_scale();

    double                 crack_driving_force();
    double                 crack_driving_force_deriv();
    double                 g_d();
    double                 gd_deriv();

    double                 strain_energy_plus();
    double                 strain_energy_minus();

    void                   set_strain_energy_plus(double H_plus);

  private:
    // Input Parameters
    double K;
    double nu;
    double G_c;
    double l_zero;

    // Computed Parameters
    //   Elastic moduli
    double lambda;
    double mu;

    // History variable
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;

    SymmetricTensor<4,3> new_cto;
    SymmetricTensor<4,3> elastic_cto;

    double degrade;
    double degrade_deriv;
    double degrade_deriv_2;
    double new_H_plus, old_H_plus; // irreversibility parameter

    // Update
    void run_3D_update (const SymmetricTensor<2,3> &incr_strain);
};


#endif
