#ifndef PHASEFIELD_CONTACT_H
#define PHASEFIELD_CONTACT_H

#include <geocentric/materials/constitutive-model.h>
#include <geocentric/utilities/line-segment.h>

using namespace dealii;


template <int dim>
class PhaseFieldContact : public ConstitutiveModel<dim>
{
  public:
    PhaseFieldContact();
    ~PhaseFieldContact();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update(const SymmetricTensor<2,dim> &incr_strain,const bool check,const Point<dim> q_point);
    void update_phasefield(const double pf);
    void update_pf_grad_normal(const Tensor<1,dim> grad_pf_normal);

    void update_frac_normal(const std::vector<Point<dim>> crack_point_set, 
                            const std::vector<LineSegment<dim>> line_set, 
                            const Point<dim> quadrature);

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
    double                 c_0(); 
    double                 xi_output(); 

    double                 crack_driving_force();
    double                 crack_driving_force_deriv();
    void                   set_crack_driving_force(double H_plus); 

    double                 degrade_output(); 
    double                 degrade_output(const double pf);
    double                 degrade_deriv_output(); 
    double                 degrade_deriv_output(const double pf);  

    double                 degrade_deriv_2_output(); 

    double                 H_plus_output(); 
    double                 yield_output();
    double                 contact_pressure_output();
    double                 shear_bulk_output();

  private:
    // Input Parameters
    // Material Parameters
    // double K; 
    double E; // if E is provided 
    double nu;
    double G_c;
    double phi;
    double sig_c;

    std::string modeltype;

    // Modeling Parameters
    double l_zero;
    double c0; 
    double xi;
    double a1; 
    double a2; 
    double n_p;

    double strain_energy_threshold;

    // Computed Parameters
    // Elastic moduli
    double lambda;
    double mu;
    double K; // If E is provided 

    double fric_coeff;
    double cohesion; 

    // Contact Parameters
    double contact_pressure;
    double shear_stress_bulk; 
    double new_friction;
    double old_friction;  
    double yield;

    bool stick; 

    Tensor<1,3>   frac_normal; 
    Tensor<1,dim> pf_gradient_norm; 
    double        pf_quadrature; // to store pf value at this quadratrue 

    // History variable
    SymmetricTensor<2,3> new_stress;
    SymmetricTensor<2,3> old_stress;

    SymmetricTensor<2,3> new_stress_bulk;
    SymmetricTensor<2,3> old_stress_bulk;

    SymmetricTensor<2,3> new_stress_interface;  
    SymmetricTensor<2,3> old_stress_interface;
    
    SymmetricTensor<2,3> stress_no_penetration; 
    SymmetricTensor<2,3> stress_friction; 

    SymmetricTensor<2,3> new_strain;
    SymmetricTensor<2,3> old_strain;

    SymmetricTensor<4,3> elastic_cto;
    SymmetricTensor<4,3> elastic_cto_inv; 

    SymmetricTensor<4,3> new_cto;
    SymmetricTensor<4,3> cto_interface; 
    SymmetricTensor<4,3> cto_no_penetration;
    SymmetricTensor<4,3> cto_friction;

    double tr_eps;

    // For Phase-field 
    double degrade;
    double degrade_deriv;
    double degrade_deriv_2;

    double new_H_plus, old_H_plus; // irreversibility parameter

    // Update
    void run_3D_update(const SymmetricTensor<2,3> &incr_strain,const bool check,const Point<dim> q_point);   
};


#endif
