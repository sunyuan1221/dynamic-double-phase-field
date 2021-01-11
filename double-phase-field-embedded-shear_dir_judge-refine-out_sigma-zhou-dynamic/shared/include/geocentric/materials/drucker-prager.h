#ifndef DRUCKER_PRAGER_H
#define DRUCKER_PRAGER_H

#include <geocentric/materials/constitutive-model.h>

using namespace dealii;


template <int dim>
class DruckerPrager : public ConstitutiveModel<dim>
{
  public:
    DruckerPrager();
    ~DruckerPrager();

    void declare_parameters(ParameterHandler &prm);
    void initialize        (ParameterHandler &prm);

    void set_stress(const SymmetricTensor<2,3> &in_situ_stress);

    void update    (const SymmetricTensor<2,dim> &incr_strain);
    void save_state();

    SymmetricTensor<2,dim> stress();
    SymmetricTensor<2,3>   full_stress();

    SymmetricTensor<2,dim> strain();
    SymmetricTensor<2,3>   full_strain();
    double                 vol_strain();
    double                 dev_strain();
    double                 equiv_plastic_strain();

    SymmetricTensor<4,dim> cto();
    double                 bulk_mod();
    double                 shear_mod();
    double                 poisson();

  private:
    // Input Parameters
    double K;
    double nu;
    double c_zero;
    double phi;
    double psi;
    double shape;
    double k1;
    std::string corners;

    bool verbose;

    // Computed Parameters
    //   Elastic moduli
    double lambda;
    double mu;
    //   Yield and potential surfaces parameters
    double cohesion;
    double cosphi;
    double cospsi;
    double sinphi;
    double sinpsi;
    double AF;
    double BF;
    double AG;
    double BG;
    double AF0;

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

    // Plastic multiplier
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

