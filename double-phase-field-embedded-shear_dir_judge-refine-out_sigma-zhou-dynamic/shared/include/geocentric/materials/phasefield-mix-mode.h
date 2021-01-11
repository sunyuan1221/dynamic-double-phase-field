#ifndef PHASEFIELD_MIX_MODE_H
#define PHASEFIELD_MIX_MODE_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;

template <int dim>
class PhaseFieldMixMode : public ConstitutiveModel<dim>
{
  public:
    PhaseFieldMixMode();
    ~PhaseFieldMixMode();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

	void update(const SymmetricTensor<2, dim> &incr_strain, const Point<dim> quadrature);
    void update_phasefield_1(const double pf);
    void update_phasefield_2(const double pf);

    void update_mode(const SymmetricTensor<2,3> strain); 
    void update_direction(const SymmetricTensor<2,3> strain);

    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();

    SymmetricTensor<2,dim> strain();
    SymmetricTensor<2,3>   full_strain();

    SymmetricTensor<4,dim> cto();

    double                 bulk_mod();
    double                 shear_mod();

    double                 fracture_energy_1();
    double                 fracture_energy_2();

    double                 length_scale();
    double                 c0_output();
    double                 xi_output();

    double                 crack_driving_force_1();
    double                 crack_driving_force_deriv_1();

    double                 crack_driving_force_2();
    double                 crack_driving_force_deriv_2();

    double                 degrade_1_output(); 
    double                 degrade_1_output(const double pf);
    double                 degrade_deriv_1_output(); 
    double                 degrade_deriv_1_output(const double pf); 

    double                 degrade_2_output(); 
    double                 degrade_2_output(const double pf);
    double                 degrade_deriv_2_output(); 
    double                 degrade_deriv_2_output(const double pf); 

    double                 H_plus_1_output(); 
    double                 H_plus_2_output();
    double                 sigma_n_output();//modify-201227

  private:
    // Input parameters
    // double K;
    double E; 
    double nu;
    double G_I;
    double G_II;
    double sig_p; 
    double phi_p;
    double phi_r;  
    double cohesion; 

    // Input modeling parameters
    double l_zero;
    std::string modeltype;

    // Modeling parameters
    double c0; 
    double xi;
    double M1;
    double M2; 
    double p; 
    double n_p;

    // Computed parameters
    // Elastic moduli
    double K; 
    double lambda;
    double mu;

    double M; // P-wave modulus 

    double tau_p; 
    double tau_r; 

    // Friction coefficient 
    double tanphi_p;
    double tanphi_r;  

    double pf_1, pf_2; 

    // Threshold crack driving force
    double H_p_1; 
    double H_p_2; 

    // Crack direction angle
    double theta;

    double sintheta; 
    double costheta; 

    bool open_flag; // = true is Mode I damage; = false is Mode II damage 

    double tr_eps; 

    double p_N;

    SymmetricTensor<2,3> stress_bulk; 

    // Directional vectors & tensors 
    Tensor<1,3>   normal, m1, m2; 
    SymmetricTensor<2,3> n_dyad_n, m1_dyad_m1, m2_dyad_m2;
    SymmetricTensor<2,3> sym_n_dyad_m1, sym_n_dyad_m2, sym_m1_dyad_m2; 

    // History variable
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;

    SymmetricTensor<4,3> new_cto;
    SymmetricTensor<4,3> elastic_cto;

    double degrade_1;
    double degrade_deriv_1;
    double degrade_deriv_deriv_1;

    double degrade_2;
    double degrade_deriv_2;
    double degrade_deriv_deriv_2;

    double new_H_plus_1, old_H_plus_1; 
    double new_H_plus_2, old_H_plus_2; 

    double new_H_surplus_2, old_H_surplus_2; 

    double H_surplus_2_max; 

    // Update
    void run_3D_update (const SymmetricTensor<2,3> &incr_strain);
};


#endif
