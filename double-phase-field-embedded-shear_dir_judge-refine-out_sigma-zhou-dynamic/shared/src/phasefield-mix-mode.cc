#include "geocentric/materials/phasefield-mix-mode.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
PhaseFieldMixMode<dim>::PhaseFieldMixMode()
{ this->name_ = "Phase-field Mix Mode"; }

template <int dim>
PhaseFieldMixMode<dim>::~PhaseFieldMixMode()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldMixMode<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldMixMode");
    //prm.declare_entry("BulkMod","0.0",Patterns::Double());
    prm.declare_entry("YoungsMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergyModeI","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergyModeII","0.0",Patterns::Double());
    prm.declare_entry("LengthScale","0.0",Patterns::Double());
    prm.declare_entry("TensilePeakStress","0.0",Patterns::Double());
    prm.declare_entry("PeakFrictionAngle","0.0",Patterns::Double());
    prm.declare_entry("ResidualFrictionAngle","0.0",Patterns::Double());
    prm.declare_entry("Cohesion","0.0",Patterns::Double()); 
    prm.declare_entry("ModelType","Brittle",Patterns::Selection("Brittle|Lorentz|Wu"));
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldMixMode<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldMixMode");
    //K          = prm.get_double("BulkMod");
    E          = prm.get_double("YoungsMod");
    nu         = prm.get_double("Poisson");
    G_I        = prm.get_double("FractureEnergyModeI");
    G_II       = prm.get_double("FractureEnergyModeII");
    l_zero     = prm.get_double("LengthScale");
    sig_p      = prm.get_double("TensilePeakStress");
    phi_p      = prm.get_double("PeakFrictionAngle")/180*M_PI; 
    phi_r      = prm.get_double("ResidualFrictionAngle")/180*M_PI; 
    cohesion   = prm.get_double("Cohesion"); 
    modeltype  = prm.get("ModelType");
  prm.leave_subsection();

  // // If K is provided 
 // mu     = convert_moduli::bulk_poisson::mu     (K,nu);
 // lambda = convert_moduli::bulk_poisson::lambda (K,nu);

  // if E is provided 
   mu     = convert_moduli::youngs_poisson::mu       (E,nu); 
   lambda = convert_moduli::youngs_poisson::lambda   (E,nu); 
   K      = convert_moduli::youngs_poisson::bulk     (E,nu); 

  // P-wave modulus 
  M = lambda + 2*mu; 

  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;

  // Friction coefficient 
  tanphi_p = std::tan(phi_p);
  tanphi_r = std::tan(phi_r);  

  theta = 0.5*(0.5*M_PI - phi_r); 

  sintheta = std::sin(theta); 
  costheta = std::cos(theta); 

  // Initialization of contact pressure, peak and residual shear stresses
  p_N = 0; 
  tau_p = cohesion + p_N*tanphi_p; 
  tau_r = p_N*tanphi_r; 

  // Initialize modeling parameters 
  if (modeltype == "Wu")
  {
    xi  = 2; 
    p   = -0.5; 
    n_p = 2;
    c0  = M_PI;   
  }
  else if (modeltype == "Lorentz")
  {
    xi  = 1; 
    p   = 1; 
    n_p = 2; 
    c0  = 8./3.; 
  }
  else if (modeltype == "Brittle")
  {
    xi  = 0; 
    n_p = 2;
    c0  = 2; 
  }

  if (modeltype == "Brittle")
  {
    H_p_1 = 0;
    H_p_2 = 0; 
  }
  else
  {
    H_p_1 = (0.5/M)*sig_p*sig_p;
    double L_fpz_1 = M*G_I*pow(sig_p,-2);

    M1 = 2*xi/(c0*l_zero)*L_fpz_1;

    H_p_2 = (0.5/mu)*(tau_p - tau_r)*(tau_p - tau_r);
    double L_fpz_2 = mu*G_II*pow(tau_p - tau_r,-2);

    M2 = 2*xi/(c0*l_zero)*L_fpz_2;

    if (modeltype == "Lorentz")
    if (l_zero > 0.5*L_fpz_1 or l_zero > 0.5*L_fpz_2)
      std::cout << "ERROR: The regularization length exceeds the upper bound" << std::endl;
    
    if (modeltype == "Wu")
    if (l_zero > 0.8488*L_fpz_1 or l_zero > 0.8488*L_fpz_2)
    {
      std::cout << "ERROR: The regularization length exceeds the upper bound" << std::endl;
    }
     // std::cout << "l_zero=" <<l_zero << std::endl;
     // std::cout << "0.8488*L_fpz_1=" <<0.8488*L_fpz_1 << std::endl;
     // std::cout << "0.8488*L_fpz_2=" <<0.8488*L_fpz_2 << std::endl;
  }

  update_phasefield_1(0.0);
  update_phasefield_2(0.0);

  new_H_plus_1 = H_p_1; 
  old_H_plus_1 = new_H_plus_1;

  new_H_plus_2 = H_p_2; 
  old_H_plus_2 = new_H_plus_2;

  new_H_surplus_2 = 0; 
  old_H_surplus_2 = new_H_surplus_2; 

  new_cto = elastic_cto;

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldMixMode<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress = in_situ_stress;
  new_strain = invert(new_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <>
void PhaseFieldMixMode<3>::update(const SymmetricTensor<2,3> &incr_strain, const Point<3> quadrature)
{
  run_3D_update(incr_strain);
}
template <>
void PhaseFieldMixMode<2>::update(const SymmetricTensor<2,2> &incr_strain, const Point<2> quadrature)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  //if (quadrature[1] > 60) // for rigid footing
  //{
  //	  E = 3e9;
  //	  nu = 0.0;

  //	  mu = convert_moduli::youngs_poisson::mu(E, nu);
  // lambda = convert_moduli::youngs_poisson::lambda(E, nu);
  //	  K = convert_moduli::youngs_poisson::bulk(E, nu);

  // elastic_cto = lambda*eye_dyad_eye + 2 * mu*big_eye;

  //}

  run_3D_update(incr_strain_2d);
}

template <int dim>
void PhaseFieldMixMode<dim>::update_phasefield_1(const double pf)
{
  if (modeltype == "Brittle")
  {
    // Degradation function for brittle fracture 
    degrade_1             = pow(1 - pf, n_p); 
    degrade_deriv_1       = - n_p*pow(1 - pf, n_p - 1); 
    degrade_deriv_deriv_1 = n_p*(n_p - 1)*pow(1 - pf, n_p - 2); 
  }
  else 
  {
    // Degradation function for cohesive fracture 
    double tmp = M1*pf*(p*pf + 1) + pow(1 - pf, n_p);

    degrade_1             = pow(1 - pf, n_p) / tmp;
    degrade_deriv_1       = - M1*pow(1 - pf, n_p - 1)*(1 + pf*(2*p + n_p - 1) + p*(n_p - 2)*pf*pf) / pow(tmp, 2);
    degrade_deriv_deriv_1 = (n_p - 1)*n_p*pow(1 - pf, n_p - 2) / tmp + 2*n_p*pow(1 - pf, n_p - 1)*(M1*p*pf + M1*(p*pf + 1) - n_p*pow(1 - pf, n_p - 1))/pow(tmp, 2)
                            + pow(1 - pf, n_p)*(2*pow(M1*p*pf + M1*(p*pf + 1) - n_p*pow(1 - pf, n_p - 1), 2) / pow(tmp, 3) - (2*M1*p + (n_p - 1)*n_p*pow(1 - pf, n_p - 2)) / pow(tmp, 2)); 
  }

  pf_1 = pf; 
}

template <int dim>
void PhaseFieldMixMode<dim>::update_phasefield_2(const double pf)
{
  if (modeltype == "Brittle")
  {
    // Degradation function for brittle fracture 
    degrade_2             = pow(1 - pf, n_p); 
    degrade_deriv_2       = - n_p*pow(1 - pf, n_p - 1); 
    degrade_deriv_deriv_2 = n_p*(n_p - 1)*pow(1 - pf, n_p - 2); 
  }
  else 
  {
    // Degradation function for cohesive fracture 
    double tmp = M2*pf*(p*pf + 1) + pow(1 - pf, n_p);

    degrade_2             = pow(1 - pf, n_p) / tmp;
    degrade_deriv_2       = - M2*pow(1 - pf, n_p - 1)*(1 + pf*(2*p + n_p - 1) + p*(n_p - 2)*pf*pf) / pow(tmp, 2);
    degrade_deriv_deriv_2 = (n_p - 1)*n_p*pow(1 - pf, n_p - 2) / tmp + 2*n_p*pow(1 - pf, n_p - 1)*(M2*p*pf + M2*(p*pf + 1) - n_p*pow(1 - pf, n_p - 1))/pow(tmp, 2)
                            + pow(1 - pf, n_p)*(2*pow(M2*p*pf + M2*(p*pf + 1) - n_p*pow(1 - pf, n_p - 1), 2) / pow(tmp, 3) - (2*M2*p + (n_p - 1)*n_p*pow(1 - pf, n_p - 2)) / pow(tmp, 2)); 
  }

  pf_2 = pf; 
}

// -----------------------------------------------------------------------------
// UPDATE NORMAL VECTOR
// -----------------------------------------------------------------------------
template <int dim> 
void PhaseFieldMixMode<dim>::update_mode(const SymmetricTensor<2,3> strain)
{
  Tensor<1,3>               principal_strains;
  std::vector<Tensor<1,3> > principal_directions(3);
  eigen_decompose<3>(strain,principal_strains,principal_directions);// maximum compressive strain at first 

  Tensor<1,3> max_principal_direction; 

  double max_principal_strain; 
  
  if (dim == 2) // For 2D plane strain condition 
  {
    if (principal_directions.at(2)[2] == 0)
    {
      max_principal_strain = principal_strains[2];

      max_principal_direction = principal_directions.at(2); 
    }
    else 
    {
      max_principal_strain = principal_strains[1];

      max_principal_direction = principal_directions.at(1);  
    }
  }
  else // For 3D condition 
  {
    max_principal_strain = principal_strains[2];
    max_principal_direction = principal_directions.at(2); //modify-201218
  }
  
  double max_principal_stress = lambda*tr_eps + 2*mu*max_principal_strain; 

  if (max_principal_stress > 0.0)
    open_flag = true;
  else  
  {
    open_flag = false; 

    if (dim == 2)
    {
      m2[0] = 0; 
      m2[1] = 0; 
      m2[2] = 1; 
    }
    else 
    {
      m2 = principal_directions.at(1);
    }

    // By Rodrigue's formula
    // when inclination of flaws to axial >=90 degree(otherwise, switch m2 and max_principal_direction)
    //modify-201222
    //normal = max_principal_direction*costheta + cross_product_3d(max_principal_direction,m2)*sintheta;
    //normal = max_principal_direction*costheta + cross_product_3d(m2,max_principal_direction)*sintheta;
      // By Rodrigue's formula
      // when inclination of flaws to axial >=90 degree(otherwise, switch m2 and max_principal_direction)
      //modify-201229
        Tensor<1,3> normal1, normal2, crack_dir_normal;
        double cos1,cos2;
        crack_dir_normal[0] = 0.5;
        crack_dir_normal[1] = 0.0;
        crack_dir_normal[2] = 0.8660254;
        normal1 = max_principal_direction*costheta + cross_product_3d(max_principal_direction,m2)*sintheta;
        normal2 = max_principal_direction*costheta + cross_product_3d(m2,max_principal_direction)*sintheta;
        cos1=abs(scalar_product(normal1,crack_dir_normal)/(normal1.norm()*crack_dir_normal.norm()));
        cos2=abs(scalar_product(normal2,crack_dir_normal)/(normal2.norm()*crack_dir_normal.norm()));
        if(cos1<cos2)//modify-201224
        normal = normal1;
        else
        normal = normal2;
       //modify-201223
  }
}

// -----------------------------------------------------------------------------
// DIRECTION VECTORS UPDATE
// -----------------------------------------------------------------------------
template <int dim> 
void PhaseFieldMixMode<dim>::update_direction(const SymmetricTensor<2,3> strain)
{
  if (pf_1 < tol_pf and pf_2 < tol_pf) // only update mode & direction when intact
    update_mode(strain); 

  Tensor<1,3> max_principal_direction; 

  if (open_flag) // Update directional vectors to principal directions for Mode I damage 
  {
    Tensor<1,3>               principal_strains;
    std::vector<Tensor<1,3> > principal_directions(3);
    eigen_decompose<3>(old_strain,principal_strains,principal_directions);// Semi-implicit method for getting the directional vectors

    if (dim == 2) // 2D plane strain condition 
    {
      if (principal_directions.at(2)[2] == 0)
        max_principal_direction = principal_directions.at(2);
      else 
        max_principal_direction = principal_directions.at(1);

      m2[0] = 0; 
      m2[1] = 0; 
      m2[2] = 1; 
    }
    else // 3D condition 
    {
      max_principal_direction = principal_directions.at(2);
      m2 = principal_directions.at(1);
//      //modify-201222
//      if(principal_directions.at(1)[1]<0)
//      {
//          m2[0] = -principal_directions.at(1)[0];
//          m2[1] = -principal_directions.at(1)[1];
//          m2[2] = -principal_directions.at(1)[2];
//      }
//      else
//        m2 = principal_directions.at(1);
//      //modify-201222
    }
    normal = max_principal_direction; 
  }

  m1 = cross_product_3d(normal, m2); 

  n_dyad_n   = symmetrize(outer_product(normal,normal)); 
  m1_dyad_m1 = symmetrize(outer_product(m1,m1)); 
  m2_dyad_m2 = symmetrize(outer_product(m2,m2));

  sym_n_dyad_m1  = symmetrize(outer_product(normal,m1) + outer_product(m1,normal)); 
  sym_n_dyad_m2  = symmetrize(outer_product(normal,m2) + outer_product(m2,normal)); 
  sym_m1_dyad_m2 = symmetrize(outer_product(m1,m2) + outer_product(m2,m1)); 
}

// -----------------------------------------------------------------------------
// ACTUAL 3D UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldMixMode<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain)
{
  // Update total strain
  new_strain = old_strain + incr_strain; 

  stress_bulk = elastic_cto*new_strain; 

  SymmetricTensor<2,3> stress_interface_tr = elastic_cto*new_strain; 

  tr_eps = trace(new_strain);

  update_direction(new_strain); 

  double eps_nn = new_strain*n_dyad_n; 
  double sig_nn = stress_bulk*n_dyad_n; 

  double tau_bulk = std::fabs(mu*new_strain*sym_n_dyad_m1); 

  p_N = - sig_nn;
  
  if (p_N < 0) p_N = 0; // no contact pressure when open

  // Update peak & residual shear stress
  tau_r = p_N*tanphi_r; 
  tau_p = p_N*tanphi_p + cohesion; 

  if (open_flag) // Mode I
  {
    if (pf_1 < tol_pf)
    {
      new_stress = stress_bulk; 

      new_cto = elastic_cto; 
    }
    else 
    {
      double lateral_coeff = lambda / (lambda + 2*mu); 

      double sig_nn_plus; 

      sig_nn_plus = std::max(0.0, sig_nn);

      SymmetricTensor<2,3> stress_plus, stress_minus; 

      SymmetricTensor<2,3> stress_n_m1_sym, stress_n_m2_sym; 

      stress_n_m1_sym = 0.5*stress_bulk*outer_product(sym_n_dyad_m1,sym_n_dyad_m1);
      stress_n_m2_sym = 0.5*stress_bulk*outer_product(sym_n_dyad_m2,sym_n_dyad_m2); 

      stress_plus = sig_nn_plus*n_dyad_n + lateral_coeff*sig_nn_plus*(m1_dyad_m1 + m2_dyad_m2); 

      new_stress = stress_bulk + (degrade_1 - 1)*stress_plus; 

      SymmetricTensor<4,3> cto_plus;

      SymmetricTensor<2,3> tmp_1 = n_dyad_n + lateral_coeff*(m1_dyad_m1 + m2_dyad_m2); 
      SymmetricTensor<2,3> tmp_2 = lambda*eye + 2*mu*n_dyad_n; 

      cto_plus = (sig_nn > 0)*outer_product(tmp_1,tmp_2);

      new_cto = elastic_cto + (degrade_1 - 1)*cto_plus;
    }

    // Update crack driving force 
    double strain_energy_degraded_1 = (sig_nn > 0)*0.5*sig_nn*sig_nn/M; 

    new_H_plus_1 = std::max(old_H_plus_1,strain_energy_degraded_1); 

    // No update of Mode II crack driving force as no shear stress on principal plane
    new_H_plus_2 = old_H_plus_2;
  }
  else // Mode II (same as phase-field shear)
  {
    // Update trial shear stress 
    SymmetricTensor<2,3> in_plane_shear_tr = 0.5*stress_interface_tr*outer_product(sym_n_dyad_m1,sym_n_dyad_m1); 
    SymmetricTensor<2,3> out_plane_shear_tr = 0.5*stress_interface_tr*outer_product(sym_n_dyad_m2,sym_n_dyad_m2);

    double tau_in_plane = sqrt(0.5)*in_plane_shear_tr.norm();
    double tau_out_plane= sqrt(0.5)*out_plane_shear_tr.norm();

    SymmetricTensor<2,3> n_in_plane; 
    
    if (tau_in_plane < tol) n_in_plane = 0; 
    else                    n_in_plane = in_plane_shear_tr/tau_in_plane;

    double gamma_increment = std::fabs(incr_strain*sym_n_dyad_m1);

    double tau_yield; 

    if (pf_2 > tol_pf) // Threshold to ensure reach peak stress
      tau_yield = tau_r; 
    else 
      tau_yield = tau_p; 

    double yield = tau_bulk - tau_yield; 

    if (yield < 0) // Stick case 
    {
      new_stress = old_stress + elastic_cto*incr_strain; 

      new_cto = elastic_cto; 

      if (pf_2 < tol_pf)
      {
        H_p_2 = (0.5/mu)*(tau_p - tau_r)*(tau_p - tau_r);
        double L_fpz_2 = mu*G_II*pow(tau_p - tau_r,-2);

        M2 = 2*xi/(c0*l_zero)*L_fpz_2;

        update_phasefield_2(pf_2); 
      }

      new_H_plus_2 = std::max(old_H_plus_2,H_p_2);
    }
    else // Slip case
    {
      // Calculate the increment of crack driving force 
      double strain_energy_increment_2 = (tau_bulk - tau_r)*gamma_increment; 

      new_H_surplus_2 = old_H_surplus_2 + strain_energy_increment_2; 

      if (new_H_surplus_2 > H_surplus_2_max) // To consider the case with initial phase-field crack
        H_surplus_2_max = new_H_surplus_2; 

      new_H_plus_2 = H_p_2 + H_surplus_2_max; 

      SymmetricTensor<4,3> cto_interface = elastic_cto 
                                           - mu*outer_product(sym_n_dyad_m1,sym_n_dyad_m1) 
                                           - mu*outer_product(sym_n_dyad_m2,sym_n_dyad_m2);                       

      new_cto = degrade_2*elastic_cto + (1 - degrade_2)*cto_interface; 

      double sig_nn_tmp = (elastic_cto*old_strain)*n_dyad_n; 
      double tau_r_tmp = (sig_nn_tmp < 0)? -sig_nn_tmp*tanphi_r:0.0;

      SymmetricTensor<2,3> stress_friction = tau_r_tmp*n_in_plane; 
      SymmetricTensor<2,3> stress_interface = stress_interface_tr - in_plane_shear_tr - out_plane_shear_tr + stress_friction;

      new_stress = degrade_2*stress_bulk + (1 - degrade_2)*stress_interface; 
    }

    // No update of Mode I crack driving force as in closed condition 
    new_H_plus_1 = old_H_plus_1; 
  }
}

// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void PhaseFieldMixMode<dim>::save_state()
{
  old_stress   = new_stress;
  old_strain   = new_strain;

  old_H_plus_1 = new_H_plus_1;
  old_H_plus_2 = new_H_plus_2;

  old_H_surplus_2 = new_H_surplus_2; 
}

// Stress
template <>
SymmetricTensor<2,3> PhaseFieldMixMode<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> PhaseFieldMixMode<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldMixMode<dim>::full_stress()
{
  return new_stress;
}

// Strain
template <>
SymmetricTensor<2,3> PhaseFieldMixMode<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> PhaseFieldMixMode<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldMixMode<dim>::full_strain()
{
  return new_strain;
}

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> PhaseFieldMixMode<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> PhaseFieldMixMode<2>::cto()
{
  SymmetricTensor<4,2> new_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    new_cto_2d[i][j][k][l] = new_cto[i][j][k][l];

  return new_cto_2d;
}


// Bulk modulus
template <int dim>
double PhaseFieldMixMode<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double PhaseFieldMixMode<dim>::shear_mod(){ return mu; }

// Critical fracture energy (critical energy release rate)
template <int dim>
double PhaseFieldMixMode<dim>::fracture_energy_1(){ return G_I; }

template <int dim>
double PhaseFieldMixMode<dim>::fracture_energy_2(){ return G_II; }

// Length scale parameter
template <int dim>
double PhaseFieldMixMode<dim>::length_scale(){ return l_zero; }

// c0 coefficient for degradation function 
template <int dim> 
double PhaseFieldMixMode<dim>::c0_output(){ return c0; }

// xi for crack density function 
template <int dim>
double PhaseFieldMixMode<dim>::xi_output(){ return xi; }

// Degradation value and its derivatives 
template <int dim> 
double PhaseFieldMixMode<dim>::degrade_1_output(){ return degrade_1; }

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_1_output(const double pf)
{
  double gd;

  if (modeltype == "Brittle")
  {
    gd = pow(1 - pf, n_p);
  }
  else 
  {
    double tmp = M1*pf*(p*pf + 1) + pow(1 - pf, n_p);

    gd = pow(1 - pf, n_p) / tmp;
  }

  return gd; 
}

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_2_output(){ return degrade_2; }

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_2_output(const double pf)
{
  double gd;

  if (modeltype == "Brittle")
  {
    gd = pow(1 - pf, n_p);
  }
  else 
  {
    double tmp = M2*pf*(p*pf + 1) + pow(1 - pf, n_p);

    gd = pow(1 - pf, n_p) / tmp;
  }

  return gd; 
}

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_deriv_1_output() { return degrade_deriv_1; }

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_deriv_1_output(const double pf)
{
  double gd_deriv; 

  if (modeltype == "Brittle")
  {
    gd_deriv = - n_p*pow(1 - pf, n_p - 1); 
  }
  else 
  {
    double tmp = M1*pf*(p*pf + 1) + pow(1 - pf, n_p);

    gd_deriv = -M1*pow(1 - pf, n_p - 1)*(1 + pf*(2*p + n_p -1) + p*(n_p - 2)*pf*pf) / pow(tmp, 2);
  }

  return gd_deriv; 
}

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_deriv_2_output() { return degrade_deriv_2; }

template <int dim> 
double PhaseFieldMixMode<dim>::degrade_deriv_2_output(const double pf)
{
  double gd_deriv; 

  if (modeltype == "Brittle")
  {
    gd_deriv = - n_p*pow(1 - pf, n_p - 1); 
  }
  else 
  {
    double tmp = M2*pf*(p*pf + 1) + pow(1 - pf, n_p);

    gd_deriv = -M2*pow(1 - pf, n_p - 1)*(1 + pf*(2*p + n_p -1) + p*(n_p - 2)*pf*pf) / pow(tmp, 2);
  }

  return gd_deriv; 
}

// Crack driving force 
template <int dim> 
double PhaseFieldMixMode<dim>::H_plus_1_output(){ return new_H_plus_1; }

template <int dim> 
double PhaseFieldMixMode<dim>::H_plus_2_output(){ return new_H_plus_2; }

//modify-201227
template <int dim>
double PhaseFieldMixMode<dim>::sigma_n_output(){ return stress_bulk*n_dyad_n; }
//modify-201227

// Crack driving forces times degradation derivatives 
template <int dim>
double PhaseFieldMixMode<dim>::crack_driving_force_1()
{
  return degrade_deriv_1*new_H_plus_1;
}

template <int dim>
double PhaseFieldMixMode<dim>::crack_driving_force_deriv_1()
{
  return degrade_deriv_deriv_1*new_H_plus_1;
}

template <int dim>
double PhaseFieldMixMode<dim>::crack_driving_force_2()
{
  return degrade_deriv_2*new_H_plus_2;
}

template <int dim>
double PhaseFieldMixMode<dim>::crack_driving_force_deriv_2()
{
  return degrade_deriv_deriv_2*new_H_plus_2;
}

// EXPLICIT INSTANTIATION
template class PhaseFieldMixMode<2>;
template class PhaseFieldMixMode<3>;
