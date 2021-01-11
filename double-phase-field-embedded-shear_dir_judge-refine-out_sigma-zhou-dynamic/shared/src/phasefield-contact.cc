#include "geocentric/materials/phasefield-contact.h"

using namespace dealii;


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
PhaseFieldContact<dim>::PhaseFieldContact()
{ this->name_ = "Phase-field Contact"; }

template <int dim>
PhaseFieldContact<dim>::~PhaseFieldContact()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldContact<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldContact");
    // prm.declare_entry("BulkMod","0.0",Patterns::Double());
    prm.declare_entry("YoungsMod","0.0",Patterns::Double());
    prm.declare_entry("Poisson","0.0",Patterns::Double());
    prm.declare_entry("FractureEnergy","0.0",Patterns::Double());
    prm.declare_entry("LengthScale","0.0",Patterns::Double()); 
    prm.declare_entry("FrictionAngle","0.0",Patterns::Double());
    prm.declare_entry("Cohesion","0.0",Patterns::Double());
    prm.declare_entry("PeakStress","0.0",Patterns::Double());
    prm.declare_entry("ModelType","Lorentz",Patterns::Anything());
  prm.leave_subsection();
}


// -----------------------------------------------------------------------------
// PARSE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldContact<dim>::initialize(ParameterHandler &prm)
{
  prm.enter_subsection("PhaseFieldContact");
    // K      = prm.get_double("BulkMod");
    E          = prm.get_double("YoungsMod");
    nu         = prm.get_double("Poisson");
    G_c        = prm.get_double("FractureEnergy");
    l_zero     = prm.get_double("LengthScale");
    phi        = prm.get_double("FrictionAngle")/180*M_PI;  
    cohesion   = prm.get_double("Cohesion"); 
    sig_c      = prm.get_double("PeakStress");
    modeltype  = prm.get("ModelType"); 
  prm.leave_subsection();

  // // If K is provided 
  // mu     = convert_moduli::bulk_poisson::mu     (K,nu);
  // lambda = convert_moduli::bulk_poisson::lambda (K,nu);

  // if E is provided 
  mu     = convert_moduli::youngs_poisson::mu       (E,nu); 
  lambda = convert_moduli::youngs_poisson::lambda   (E,nu); 
  K      = convert_moduli::youngs_poisson::bulk     (E,nu); 

  elastic_cto = lambda*eye_dyad_eye + 2*mu*big_eye;
  elastic_cto_inv = invert(elastic_cto);

  new_cto    = elastic_cto;

  new_stress_bulk = 0; 
  old_stress_bulk = new_stress_bulk; 

  new_stress_interface = 0; 
  old_stress_interface = new_stress_interface;

  fric_coeff = std::tan(phi);

  contact_pressure = 0; 

  new_friction = contact_pressure*fric_coeff + cohesion; 
  old_friction = new_friction; 

  // Initialize modeling parameters 
  if (modeltype == "Wu")
  {
    xi = 2; 
    n_p = 2; 
    a2 = -0.5; 

    c0 = M_PI;   
  }
  else if (modeltype == "Lorentz")
  {
    xi = 1; 
    n_p = 2; 
    a2 = 1; 

    c0 = 8./3.; 
  }
  else if (modeltype == "Brittle")
  {
    xi = 0; 

    c0 = 2; 
  }
  else 
  {
    std::cerr << "ERROR: Crack density function type does not exist"<< std::endl;
    exit(1);
  } 

  if (modeltype == "Brittle")
  {
    strain_energy_threshold = 0;
  }
  else
  {
    strain_energy_threshold = 0.5*sig_c*sig_c/E; 
    double L_fpz = E*G_c*pow(sig_c,-2);
    a1 = 2*xi/(c0*l_zero)*L_fpz;

    if (modeltype == "Lorentz" and l_zero > 0.5*L_fpz)
      std::cout << "ERROR: The regularization length exceeds the upper bound" << std::endl;
    else if (modeltype == "Wu" and l_zero > 0.8488*L_fpz)
      std::cout << "ERROR: The regularization length exceeds the upper bound" << std::endl;
  } 

  new_H_plus = strain_energy_threshold;

  update_phasefield(0);

  // frac_normal[0] = - 0.6585046079; 
  // frac_normal[1] = 0.7525766947;
  // frac_normal[2] = 0.0;

  // frac_normal[0] = - 0.5734623444; 
  // frac_normal[1] = 0.8192319205;
  // frac_normal[2] = 0.0;

  // frac_normal[0] = -0.7071067812;
  // frac_normal[1] = 0.7071067812; 
  // frac_normal[2] = 0.0; 

  // frac_normal[0] = - 0.1961161351; 
  // frac_normal[1] = 0.9805806757;
  // frac_normal[2] = 0.0;

  // frac_normal[0] = 0.0; 
  // frac_normal[1] = 1.0;
  // frac_normal[2] = 0.0;

  frac_normal[0] = 0.0; 
  frac_normal[1] = 0.0;
  frac_normal[2] = 0.0;

  save_state();
}


// -----------------------------------------------------------------------------
// SET STRESS
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldContact<dim>::set_stress(const SymmetricTensor<2,3> &in_situ_stress)
{
  new_stress         = in_situ_stress;
  new_strain         = invert(new_cto)*in_situ_stress;

  save_state();
}


// -----------------------------------------------------------------------------
// UPDATE STATE
// -----------------------------------------------------------------------------
template <>
void PhaseFieldContact<3>::update(const SymmetricTensor<2,3> &incr_strain,const bool check,const Point<3> q_point)
{
  run_3D_update(incr_strain,check,q_point); 
}
template <>
void PhaseFieldContact<2>::update(const SymmetricTensor<2,2> &incr_strain,const bool check,const Point<2> q_point)
{
  SymmetricTensor<2,3> incr_strain_2d;

  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    incr_strain_2d[i][j] = incr_strain[i][j];

  run_3D_update(incr_strain_2d,check,q_point); 
}

template <int dim>
void PhaseFieldContact<dim>::update_phasefield(const double pf)
{
  if (modeltype == "Brittle")
  {
    // Degradation function for brittle fracture 
    degrade         = (1 - pf)*(1 - pf); 
    degrade_deriv   = 2*pf - 2; 
    degrade_deriv_2 = 2; 
  }
  else 
  {
    // Degradation function for cohesive fracture 
    double tmp = a1*pf*(a2*pf + 1) + pow(1 - pf, n_p);

    degrade         = pow(1 - pf, n_p) / tmp;
    degrade_deriv   = -a1*pow(1 - pf, n_p - 1)*(1 + pf*(2*a2 + n_p -1) + a2*(n_p - 2)*pf*pf) / pow(tmp, 2);
    degrade_deriv_2 = (n_p - 1)*n_p*pow(1 - pf, n_p - 2) / tmp + 2*n_p*pow(1 - pf, n_p - 1)*(a1*a2*pf + a1*(a2*pf + 1) - n_p*pow(1 - pf, n_p - 1))/pow(tmp, 2)
                      + pow(1 - pf, n_p)*(2*pow(a1*a2*pf + a1*(a2*pf + 1) - n_p*pow(1 - pf, n_p - 1), 2) / pow(tmp, 3) - (2*a1*a2 + (n_p - 1)*n_p*pow(1 - pf, n_p - 2)) / pow(tmp, 2)); 
  }

  pf_quadrature = pf;
}

template <int dim> 
void PhaseFieldContact<dim>::update_pf_grad_normal(const Tensor<1,dim> grad_pf_normal)
{
  pf_gradient_norm = grad_pf_normal;
}

// -----------------------------------------------------------------------------
// FRACTURE NORMAL VECTOR 
// -----------------------------------------------------------------------------
template <int dim> 
void PhaseFieldContact<dim>::update_frac_normal(const std::vector<Point<dim>> crack_point_set, 
                                                const std::vector<LineSegment<dim>> line_set, 
                                                const Point<dim> quadrature)
{
  if (pf_quadrature < tol)
  {
    frac_normal = 0; 

    return; 
  }

  Tensor<1,dim> crack_orientation; 

  double min_distance = 1e8; 

  unsigned line_index; 
  unsigned point_index; 

  Point<dim> point_coord; 

  bool master_side = false; 

  double size_point = crack_point_set.size(); 
  double size_line  = line_set.size(); 

  for (unsigned n=0; n<size_point; ++n)
  {
    double tmp_distance = quadrature.distance(crack_point_set.at(n)); 

    if (tmp_distance < min_distance)
    {
      min_distance = tmp_distance; 
      point_index = n; 
      point_coord = crack_point_set.at(n); 
    }
  }

  if (point_index == 0) 
  {
    line_index = 0; 
    LineSegment<dim> tmp_line = line_set.at(line_index);

    crack_orientation = tmp_line.direction_vector(); 
    master_side = tmp_line.crack_side(quadrature); 
  }
  else if (point_index == size_point - 1)
  {
    line_index = size_line - 1; 
    LineSegment<dim> tmp_line = line_set.at(line_index);

    crack_orientation = tmp_line.direction_vector();  
    master_side = tmp_line.crack_side(quadrature);
  }
  else 
  {
    LineSegment<dim> line_q; 
    line_q.input_points(quadrature,point_coord); 

    LineSegment<dim> tmp_line_front = line_set.at(point_index - 1);
    LineSegment<dim> tmp_line_back = line_set.at(point_index); 

    double proj_front = tmp_line_front.projection(line_q); 
    double proj_back  = tmp_line_back.projection(line_q); 

    if (proj_front >= proj_back)
    {
      line_index = point_index - 1; 
      LineSegment<dim> tmp_line = line_set.at(line_index);

      crack_orientation = tmp_line.direction_vector(); 
      master_side = tmp_line.crack_side(quadrature); 
    }
    else 
    {
      line_index = point_index; 
      LineSegment<dim> tmp_line = line_set.at(line_index); 

      crack_orientation = tmp_line.direction_vector(); 
      master_side = tmp_line.crack_side(quadrature); 
    }
  }

  if (master_side)
  {
    frac_normal[0] = -crack_orientation[1];
    frac_normal[1] = crack_orientation[0]; 
    frac_normal[2] = 0; 
  }
  else 
  {
    frac_normal[0] = crack_orientation[1];
    frac_normal[1] = -crack_orientation[0]; 
    frac_normal[2] = 0; 
  }
}

// -----------------------------------------------------------------------------
// ACTUAL 3D UPDATE
// -----------------------------------------------------------------------------
template <int dim>
void PhaseFieldContact<dim>::run_3D_update(const SymmetricTensor<2,3> &incr_strain,const bool check,const Point<dim> q_point)
{
  new_strain = old_strain + incr_strain;
  double tr_eps = trace(new_strain);

  unsigned fault_index; 

  if (pf_quadrature > 1e-4)
  {
    Tensor<1,3> frac_normal_1, frac_normal_2; 

    frac_normal_1[0] = 0; 
    frac_normal_1[1] = 1.0; 
    frac_normal_1[2] = 0; 

    frac_normal_2[0] = - 0.8660254038; 
    frac_normal_2[1] = 0.5; 
    frac_normal_2[2] = 0; 

    double a = 0; 
    double b = 1.0; 
    double c = - 5.0;  

    double d = - 1.732050808; 
    double e = 1.0; 
    double f = - 1.241449748;

    double x_coord = q_point[0]; 
    double y_coord = q_point[1];

    double distance_1 = std::fabs(a*x_coord + b*y_coord + c) * pow(a*a + b*b, -0.5);
    double distance_2 = std::fabs(d*x_coord + e*y_coord + f) * pow(d*d + e*e, -0.5); 

    // if (distance_1 > distance_2/* or x_coord < 2.17*/)
    //   frac_normal = frac_normal_2;
    // else 
    //   frac_normal = frac_normal_1;

    if (5.0 - l_zero <= y_coord and y_coord <= 5.0 + l_zero)
    {
      frac_normal = frac_normal_1;

      fault_index = 1; 
    }
    else if (distance_1 < distance_2)
    {
      frac_normal = frac_normal_1; 

      fault_index = 1;
    }
    else 
    {
      frac_normal = frac_normal_2;

      fault_index = 2;
    }
  }
  
  // Calculate trial stress tensor
  new_stress_bulk = elastic_cto*new_strain;

  // SymmetricTensor<2,3> stress_interface_tr = old_stress_interface + elastic_cto*incr_strain; 
  SymmetricTensor<2,3> stress_interface_tr = elastic_cto*new_strain; 

  // Directional stress decomposition 
  Tensor<1,3> tangential_s, tangential_t; 

  // plane strain condition 
  tangential_t[0] = 0; 
  tangential_t[1] = 0; 
  tangential_t[2] = 1; 

  tangential_s[0] = frac_normal[1];
  tangential_s[1] = - frac_normal[0]; 
  tangential_s[2] = 0; 

  SymmetricTensor<2,3> m_rr = symmetrize(outer_product(frac_normal,frac_normal)); 
  SymmetricTensor<2,3> m_ss = symmetrize(outer_product(tangential_s,tangential_s)); 
  SymmetricTensor<2,3> m_tt = symmetrize(outer_product(tangential_t,tangential_t));

  SymmetricTensor<2,3> sym_m_rs = symmetrize(outer_product(frac_normal,tangential_s) + outer_product(tangential_s,frac_normal)); 
  SymmetricTensor<2,3> sym_m_rt = symmetrize(outer_product(frac_normal,tangential_t) + outer_product(tangential_t,frac_normal)); 
  SymmetricTensor<2,3> sym_m_st = symmetrize(outer_product(tangential_t,tangential_s) + outer_product(tangential_s,tangential_t)); 

  // Update contact pressure 
  double normal_eps = new_strain*m_rr;
  double normal_sig = new_stress_bulk*m_rr;
  // contact_pressure = std::abs(normal_sig); 
  contact_pressure = -normal_sig;

  if (contact_pressure < 0)
    contact_pressure = 0;

  // Update bulk shear stress
  shear_stress_bulk = mu*new_strain*sym_m_rs;
  double out_plane_shear_tr = 0.5*stress_interface_tr*sym_m_rt;
  double in_plane_shear_tr = 0.5*stress_interface_tr*sym_m_rs;

  // Case #1: undamaged 
  if (pf_quadrature < tol) 
  {
    new_cto = elastic_cto; 
    new_stress = new_stress_bulk;

    return;
  }

  // Opening  
  if (/*contact_pressure < 0*/false) // Check for physically consistent crack orientation dependent degradation later 
  {
    new_cto = degrade*elastic_cto; 
    new_stress = new_cto*new_strain; 

    return; 
  }

  // // Check yield function in crack 
  if (check and pf_quadrature > 1e-4)
  {   
    new_friction = contact_pressure*fric_coeff + cohesion;  
  }

  yield = (shear_stress_bulk < 0)? (- shear_stress_bulk - new_friction):(shear_stress_bulk - new_friction); 

  // Case #3: stick condition 
  if (yield < tol or fault_index == 1)
  { 
    new_cto = elastic_cto; 

    new_stress_interface = stress_interface_tr;

    new_stress = degrade*new_stress_bulk + (1 - degrade)*new_stress_interface; 

    new_H_plus = old_H_plus; 

    return; 
  }

  // Case #4: slip condition 
  else 
  {
    // Phase-field contact formulation for stress 
    stress_no_penetration = stress_interface_tr - in_plane_shear_tr*sym_m_rs - out_plane_shear_tr*sym_m_rt; 
    stress_friction       = (in_plane_shear_tr < 0)? - new_friction*sym_m_rs:new_friction*sym_m_rs;

    new_stress_interface = stress_no_penetration + stress_friction; 

    new_stress = degrade*new_stress_bulk + (1 - degrade)*new_stress_interface;

    // Update cto 
    cto_no_penetration = elastic_cto - mu*outer_product(sym_m_rs,sym_m_rs) - mu*outer_product(sym_m_rt,sym_m_rt);
    // cto_friction       = (in_plane_shear_tr < 0)? fric_coeff*(lambda*outer_product(sym_m_rs,eye) + 2*mu*outer_product(sym_m_rs,m_rr))
    //                                               : -fric_coeff*(lambda*outer_product(sym_m_rs,eye) + 2*mu*outer_product(sym_m_rs,m_rr));

    cto_interface = cto_no_penetration/* + cto_friction*/;

    new_cto = degrade*elastic_cto + (1 - degrade)*cto_interface;

    new_H_plus = old_H_plus; // static crack

  }
}

// -----------------------------------------------------------------------------
// COMMON METHODS: GET, SAVE, ETC.
// -----------------------------------------------------------------------------
// Save state
template <int dim>
void PhaseFieldContact<dim>::save_state()
{
  old_stress = new_stress;
  old_strain = new_strain;
  old_H_plus = new_H_plus;

  old_stress_bulk = new_stress_bulk; 
  old_stress_interface = new_stress_interface;

  old_friction = new_friction; 
}

// Stress
template <>
SymmetricTensor<2,3> PhaseFieldContact<3>::stress()
{
  return new_stress;
}
template <>
SymmetricTensor<2,2> PhaseFieldContact<2>::stress()
{
  SymmetricTensor<2,2> new_stress_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_stress_2d[i][j] = new_stress[i][j];

  return new_stress_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldContact<dim>::full_stress()
{
  return new_stress;
}

// Strain
template <>
SymmetricTensor<2,3> PhaseFieldContact<3>::strain()
{
  return new_strain;
}
template <>
SymmetricTensor<2,2> PhaseFieldContact<2>::strain()
{
  SymmetricTensor<2,2> new_strain_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
    new_strain_2d[i][j] = new_strain[i][j];

  return new_strain_2d;
}

template <int dim>
SymmetricTensor<2,3> PhaseFieldContact<dim>::full_strain()
{
  return new_strain;
}

template <int dim>
double PhaseFieldContact<dim>::equiv_plastic_strain(){ return 0; }

// Consistent tangent operator (CTO)
template <>
SymmetricTensor<4,3> PhaseFieldContact<3>::cto()
{
  return new_cto;
}
template <>
SymmetricTensor<4,2> PhaseFieldContact<2>::cto()
{
  SymmetricTensor<4,2> new_cto_2d;
  for (unsigned i=0; i<2; ++i)
  for (unsigned j=i; j<2; ++j)
  for (unsigned k=0; k<2; ++k)
  for (unsigned l=k; l<2; ++l)
    new_cto_2d[i][j][k][l] = new_cto[i][j][k][l];

  return new_cto_2d;
}

// MATERIAL AND MODEL PARAMETERS
// Bulk modulus
template <int dim>
double PhaseFieldContact<dim>::bulk_mod(){ return K; }

// Shear modulus
template <int dim>
double PhaseFieldContact<dim>::shear_mod(){ return mu; }

// Critical fracture energy (critical energy release rate)
template <int dim>
double PhaseFieldContact<dim>::fracture_energy(){ return G_c; }

// Length scale parameter
template <int dim>
double PhaseFieldContact<dim>::length_scale(){ return l_zero; }

// c0 coefficient for degradation function 
template <int dim> 
double PhaseFieldContact<dim>::c_0(){ return c0; }

// xi for crack density function 
template <int dim>
double PhaseFieldContact<dim>::xi_output(){ return xi; }

// OUTPUT FUNCTION
// ---------------
// Degradation 
template <int dim> 
double PhaseFieldContact<dim>::degrade_output(){ return degrade; }

template <int dim> 
double PhaseFieldContact<dim>::degrade_output(const double pf)
{
  double gd;

  if (modeltype == "Brittle")
  {
    gd = std::max((1 - pf)*(1 - pf), 1e-6);
  }
  else
  {
    double tmp = a1*pf*(a2*pf + 1) + pow(1 - pf, n_p);

    gd = pow(1 - pf, n_p) / tmp;
  }

  return gd; 
}

// Degradation derivatives
template <int dim> 
double PhaseFieldContact<dim>::degrade_deriv_output(){ return degrade_deriv; }

template <int dim> 
double PhaseFieldContact<dim>::degrade_deriv_output(const double pf)
{
  double gd_deriv; 

  if (modeltype == "Brittle")
  {
    gd_deriv = -2*(1 - pf);
  }
  else
  {
    double tmp = a1*pf*(a2*pf + 1) + pow(1 - pf, n_p);

    gd_deriv = -a1*pow(1 - pf, n_p - 1)*(1 + pf*(2*a2 + n_p -1) + a2*(n_p - 2)*pf*pf) / pow(tmp, 2);
  }

  return gd_deriv; 
}

// Degradation second derivatives 
template <int dim> 
double PhaseFieldContact<dim>::degrade_deriv_2_output(){ return degrade_deriv_2; }

template <int dim> 
double PhaseFieldContact<dim>::H_plus_output(){ return new_H_plus; }

template <int dim>
double PhaseFieldContact<dim>::yield_output(){ return yield; }

template <int dim> 
double PhaseFieldContact<dim>::contact_pressure_output(){ return contact_pressure; }

template <int dim> 
double PhaseFieldContact<dim>::shear_bulk_output(){ return shear_stress_bulk; }

// Crack driving forces
template <int dim>
double PhaseFieldContact<dim>::crack_driving_force()
{
  return degrade_deriv*new_H_plus;
}

template <int dim>
double PhaseFieldContact<dim>::crack_driving_force_deriv()
{
  return degrade_deriv_2*new_H_plus;
}

// SET
// ---
// Set crack driving force to initialize crack 
template <int dim> 
void PhaseFieldContact<dim>::set_crack_driving_force(double H_plus)
{
  new_H_plus = H_plus; 
  old_H_plus = new_H_plus; 
}


// EXPLICIT INSTANTIATION
template class PhaseFieldContact<2>;
template class PhaseFieldContact<3>;
