#ifndef PHASEFIELD_SHEAR_MONO_H
#define PHASEFIELD_SHEAR_MONO_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class PhaseFieldShearMono : public ConstitutiveModel<dim>
{
  public:
    PhaseFieldShearMono();
    ~PhaseFieldShearMono();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain);
    void update_phasefield(const double pf, const double pf_extra = 0);

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
    double                 c_0(); 
    double                 xi_output(); 
    double                 critical_stress();

    double                 crack_driving_force();
    double                 crack_driving_force_deriv();
    SymmetricTensor<2,dim> crack_driving_force_deriv_strain(); 

    double                 strain_energy_plus();
    double                 strain_energy_minus();

    double                 lagrange_multiplier();
    void                   set_lagrange_multiplier(double new_lagrange_mult);
   
    double                 degradation();
    double                 degradation_deriv(); 
    double                 degradation_deriv_2();  

  private:
    // Input Parameters
    double K;
    // double E; 
    double nu;
    double G_c;
    double l_zero;
    double sig_c;
    double psi_c;
    double sig_r; 
    double g1; 
    double phi_1; 

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

    SymmetricTensor<2,3> new_cto_pf; 

    double degrade;
    double degrade_deriv;
    double degrade_deriv_2;
    double degrade_extra; 
    double a1; 
    double a2; 
    double b;
    double c0; 
    double xi;
    double new_H_plus, old_H_plus; // irreversibility parameter
    double lagrange_mult;

    SymmetricTensor<2,3> new_H_plus_deriv_strain; 

    // Update
    void run_3D_update (const SymmetricTensor<2,3> &incr_strain);
};


#endif
