#include <geocentric/materials/lode-angle.h>

using namespace dealii;


// =============================================================================
// Gudehus--Argyris Scaling Function (Gudehus 1973, Argyris et al. 1974)
// =============================================================================
GudehusArgyris::GudehusArgyris()
  :
  rho(1.0)
{}

GudehusArgyris::GudehusArgyris(double ellipticity)
  :
  rho(ellipticity)
{}

GudehusArgyris::~GudehusArgyris()
{}


// -----------------------------------------------------------------------------
// SET/GET METHODS
// -----------------------------------------------------------------------------
void GudehusArgyris::set_ellipticity(double ellipticity) { rho = ellipticity; }

double GudehusArgyris::value(){ return Zeta; }
Tensor<1,3> GudehusArgyris::grad(){ return gradZeta; }
SymmetricTensor<2,3> GudehusArgyris::hess(){ return hessZeta; }


// -----------------------------------------------------------------------------
// UPDATE SCALING FUNCTION AND ITS DERIVATIES
// -----------------------------------------------------------------------------
void GudehusArgyris::update(const Tensor<1,3> &sig)
{
  // Lode's angle
  LodeAngle lode_angle;
            lode_angle.update(sig);

  double theta = lode_angle.value();
  double cos3theta = std::cos(3*theta);
  double sin3theta = sqrt(1 - cos3theta*cos3theta);
  Tensor<1,3> gradTheta = lode_angle.grad();
  SymmetricTensor<2,3> hessTheta = lode_angle.hess();
  Tensor<2,3> gradTheta_dyad_gradTheta = outer_product(gradTheta,gradTheta);

  // Scaling function and its derivatives
  const double one_over_tworho = 0.5/rho;
  Zeta = one_over_tworho*((1 + rho) + (1 - rho)*cos3theta);

  double dZeta_dTheta = one_over_tworho*(-3*(1 - rho)*sin3theta);
  gradZeta = dZeta_dTheta*gradTheta;

  double d2Zeta_dTheta2 = one_over_tworho*(-9*(1 - rho)*cos3theta);
  hessZeta = dZeta_dTheta*hessTheta + d2Zeta_dTheta2*symmetrize(gradTheta_dyad_gradTheta);
}



// =============================================================================
// Hyperbolic Mohr-Coulomb Approximation Function (Abbo and Sloan, 1995)
// =============================================================================
HyperbolicMohrCoulomb::HyperbolicMohrCoulomb()
  :
  Theta_t(25*M_PI/180)
{}

HyperbolicMohrCoulomb::HyperbolicMohrCoulomb(double limit_lode_angle)
  :
  Theta_t(limit_lode_angle*M_PI/180)
{}

HyperbolicMohrCoulomb::~HyperbolicMohrCoulomb()
{}


// -----------------------------------------------------------------------------
// SET METHODS
// -----------------------------------------------------------------------------
void HyperbolicMohrCoulomb::set_limit_lode_angle(double limit_lode_angle)
{
  Theta_t = limit_lode_angle;
}


// -----------------------------------------------------------------------------
// UPDATE SCALING FUNCTION AND ITS DERIVATIES
// -----------------------------------------------------------------------------
void HyperbolicMohrCoulomb::update(const Tensor<1,3> &sig)
{
  // Lode's angle
  LodeAngle lode_angle;
            lode_angle.update(sig);

  Theta = lode_angle.value() - M_PI/6; // translate to the range between -6/pi to 6/pi
  gradTheta = lode_angle.grad();
  hessTheta = lode_angle.hess();
  gradTheta_dyad_gradTheta = symmetrize(outer_product(gradTheta,gradTheta));

  // // Friction angle value
  // const double sinphi = std::sin(phi);
  // const double cosphi = std::cos(phi);

  // // Scaling function and its derivatives
  // double dZeta_dTheta = 0;
  // double d2Zeta_dTheta2 = 0;
  // double dZeta_dTheta_deriv_phi = 0;

  // if (std::abs(Theta) < theta_t)
  // {
  //   double costheta = std::cos(theta);
  //   double sintheta = std::sin(theta);

  //   Zeta = costheta - root13*sinphi*sintheta;
  //   dZeta_dTheta = -sintheta - root13*sinphi*costheta;
  //   d2Zeta_dTheta2 = -costheta + root13*sinphi*sintheta;

  //   Zeta_deriv_phi = -root13*cosphi*sintheta;
  //   dZeta_dTheta_deriv_phi = -root13*cosphi*costheta;
  // }
  // else
  // {
  //   double sin3theta = std::sin(3*theta);
  //   double cos3theta = std::cos(3*theta);

  //   double costheta_t  = std::cos(theta_t);
  //   double sintheta_t  = std::sin(theta_t);
  //   double tantheta_t  = std::tan(theta_t);
  //   double cos3theta_t = std::cos(3*theta_t);
  //   double tan3theta_t = std::tan(3*theta_t);

  //   double sgn_theta = (theta < 0) ? -1 : (theta >= 0);

  //   double A = one_third*costheta_t*(3 + tantheta_t*tan3theta_t + root13*sgn_theta*(tan3theta_t - 3*tantheta_t)*sinphi);
  //   double B = (1./(3*cos3theta_t))*(sgn_theta*sintheta_t + root13*sinphi*costheta_t);

  //   Zeta = A - B*sin3theta;
  //   dZeta_dTheta = -3*B*cos3theta;
  //   d2Zeta_dTheta2 = 3*B*sin3theta;

  //   double A_deriv_phi = one_third*costheta_t*(root13*sgn_theta*(tan3theta_t - 3*tantheta_t)*cosphi);
  //   double B_deriv_phi = (1./(3*cos3theta_t))*(root13*cosphi*costheta_t);
  //   Zeta_deriv_phi = A_deriv_phi - B_deriv_phi*sin3theta;
  //   dZeta_dTheta_deriv_phi = -3*B_deriv_phi*cos3theta;
  // }

  // gradZeta = dZeta_dTheta*gradTheta;
  // hessZeta = dZeta_dTheta*hessTheta + d2Zeta_dTheta2*symmetrize(gradTheta_dyad_gradTheta);

  // gradZeta_deriv_phi = dZeta_dTheta_deriv_phi*gradTheta;
}


double HyperbolicMohrCoulomb::value(const double Angle)
{
  // Friction angle value
  const double sinAngle = std::sin(Angle);

  if (std::abs(Theta) < Theta_t)
  {
    double cosTheta = std::cos(Theta);
    double sinTheta = std::sin(Theta);

    return cosTheta - root13*sinAngle*sinTheta;
  }
  else
  {
    double sin3Theta = std::sin(3*Theta);

    double cosTheta_t  = std::cos(Theta_t);
    double sinTheta_t  = std::sin(Theta_t);
    double tanTheta_t  = std::tan(Theta_t);
    double cos3Theta_t = std::cos(3*Theta_t);
    double tan3Theta_t = std::tan(3*Theta_t);

    double sgn_Theta = (Theta < 0) ? -1 : (Theta >= 0);

    double A = one_third*cosTheta_t*(3 + tanTheta_t*tan3Theta_t + root13*sgn_Theta*(tan3Theta_t - 3*tanTheta_t)*sinAngle);
    double B = (1./(3*cos3Theta_t))*(sgn_Theta*sinTheta_t + root13*sinAngle*cosTheta_t);

    return A - B*sin3Theta;
  }
}

Tensor<1,3> HyperbolicMohrCoulomb::grad(const double Angle)
{
  const double sinAngle = std::sin(Angle);

  if (std::abs(Theta) < Theta_t)
  {
    double cosTheta = std::cos(Theta);
    double sinTheta = std::sin(Theta);

    return (-sinTheta - root13*sinAngle*cosTheta)*gradTheta;
  }
  else
  {
    double cos3Theta = std::cos(3*Theta);

    double cosTheta_t  = std::cos(Theta_t);
    double sinTheta_t  = std::sin(Theta_t);
    double cos3Theta_t = std::cos(3*Theta_t);

    double sgn_Theta = (Theta < 0) ? -1 : (Theta >= 0);

    double B = (1./(3*cos3Theta_t))*(sgn_Theta*sinTheta_t + root13*sinAngle*cosTheta_t);

    return (-3*B*cos3Theta)*gradTheta;
  }
}

SymmetricTensor<2,3> HyperbolicMohrCoulomb::hess(const double Angle)
{
  const double sinAngle = std::sin(Angle);

  if (std::abs(Theta) < Theta_t)
  {
    double cosTheta = std::cos(Theta);
    double sinTheta = std::sin(Theta);

    return (-cosTheta + root13*sinAngle*sinTheta)*gradTheta_dyad_gradTheta
           + (-sinTheta - root13*sinAngle*cosTheta)*hessTheta;
  }
  else
  {
    double cos3Theta = std::cos(3*Theta);
    double sin3Theta = std::sin(3*Theta);

    double cosTheta_t  = std::cos(Theta_t);
    double sinTheta_t  = std::sin(Theta_t);
    double cos3Theta_t = std::cos(3*Theta_t);

    double sgn_Theta = (Theta < 0) ? -1 : (Theta >= 0);

    double B = (1./(3*cos3Theta_t))*(sgn_Theta*sinTheta_t + root13*sinAngle*cosTheta_t);

    return (3*B*sin3Theta)*gradTheta_dyad_gradTheta
           + (-3*B*cos3Theta)*hessTheta;
  }
}

double HyperbolicMohrCoulomb::value_deriv_angle(const double Angle)
{
  const double cosAngle = std::cos(Angle);

  if (std::abs(Theta) < Theta_t)
  {
    double sinTheta = std::sin(Theta);

    return -root13*cosAngle*sinTheta;
  }
  else
  {
    double sin3Theta = std::sin(3*Theta);

    double cosTheta_t  = std::cos(Theta_t);
    double tanTheta_t  = std::tan(Theta_t);
    double cos3Theta_t = std::cos(3*Theta_t);
    double tan3Theta_t = std::tan(3*Theta_t);

    double sgn_Theta = (Theta < 0) ? -1 : (Theta >= 0);

    double A_deriv_Angle = one_third*cosTheta_t*(root13*sgn_Theta*(tan3Theta_t - 3*tanTheta_t)*cosAngle);
    double B_deriv_Angle = (1./(3*cos3Theta_t))*(root13*cosAngle*cosTheta_t);

    return A_deriv_Angle - B_deriv_Angle*sin3Theta;
  }
}

Tensor<1,3> HyperbolicMohrCoulomb::grad_deriv_angle(const double Angle)
{
  const double cosAngle = std::cos(Angle);

  if (std::abs(Theta) < Theta_t)
  {
    double cosTheta = std::cos(Theta);

    return (-root13*cosAngle*cosTheta)*gradTheta;
  }
  else
  {
    double cos3Theta = std::cos(3*Theta);

    double cosTheta_t  = std::cos(Theta_t);
    double cos3Theta_t = std::cos(3*Theta_t);

    double B_deriv_Angle = (1./(3*cos3Theta_t))*(root13*cosAngle*cosTheta_t);

    return (-3*B_deriv_Angle*cos3Theta)*gradTheta;
  }
}


// =============================================================================
// LODE ANGLE (common for many scaling functions)
// =============================================================================
LodeAngle::LodeAngle()
{}

LodeAngle::LodeAngle(const Tensor<1,3> &sig)
{ update(sig); }

LodeAngle::~LodeAngle()
{}


// -----------------------------------------------------------------------------
// SET/GET METHODS
// -----------------------------------------------------------------------------
double LodeAngle::value(){ return Theta; }
Tensor<1,3> LodeAngle::grad(){ return gradTheta; }
SymmetricTensor<2,3> LodeAngle::hess(){ return hessTheta; }


// -----------------------------------------------------------------------------
// UPDATE LODE ANGLE AND ITS DERIVATIVES
// -----------------------------------------------------------------------------
void LodeAngle::update(const Tensor<1,3> &sig)
{
  // Lode's angle
  // ------------
  Tensor<1,3> S;
  double P = one_third * (sig[0] + sig[1] + sig[2]);
  S[0] = sig[0] - P;
  S[1] = sig[1] - P;
  S[2] = sig[2] - P;

  double trS2 = S[0]*S[0]      + S[1]*S[1]      + S[2]*S[2];
  double trS3 = S[0]*S[0]*S[0] + S[1]*S[1]*S[1] + S[2]*S[2]*S[2];

  const double chi = (fabs(trS2) > tol) ? sqrt(trS2) : sqrt(tol);
  const double one_over_chi  = 1./chi;
  const double one_over_chi2 = pow(chi,-2);
  const double one_over_chi3 = pow(chi,-3);
  // const double one_over_chi5 = pow(chi,-5);

  // For Lode's angle defined from 0 to pi/3
  double cos3theta = root6*trS3*one_over_chi3;

  // bound Lode's angle to the limits (0 and pi/3)
  if (fabs(cos3theta + 1) < tol)
    cos3theta = -1;
  else if (fabs(cos3theta - 1) < tol)
    cos3theta =  1;

  Theta = one_third*std::acos(cos3theta);

  // Gradient of Lode's angle w.r.t. principal stresses
  // --------------------------------------------------
  gradTheta = 0;

  // Borja (2013)
  /*
  double const01 = 3*one_over_chi3;
  double const02 = -3*trS3*one_over_chi5;
  double const03 = -1*one_over_chi;

  // Derivative of Lode's angle w.r.t. Y (Eq. (6.47))
  double dTheta_dY = -2/(root6*sin3theta);

  Tensor<1,3> gradTheta;
  gradTheta[0] = dTheta_dY * (const01*S[0]*S[0] + const02*S[0] + const03);
  gradTheta[1] = dTheta_dY * (const01*S[1]*S[1] + const02*S[1] + const03);
  gradTheta[2] = dTheta_dY * (const01*S[2]*S[2] + const02*S[2] + const03);
  */

  // Nodargi and Bisegna (2016)
  double Theta_plus_23pi  = Theta + two_thirds*M_PI;
  double Theta_minus_23pi = Theta - two_thirds*M_PI;

  Tensor<1,3> gradTheta;
  gradTheta[0] = -root23*one_over_chi*std::sin(Theta_plus_23pi);
  gradTheta[1] = -root23*one_over_chi*std::sin(Theta_minus_23pi);
  gradTheta[2] = -root23*one_over_chi*std::sin(Theta);


  // Hessian of Lode's angle w.r.t. principal stresses
  // -------------------------------------------------
  hessTheta = 0;

  // Borja (2013)
  /*
  double const11 = 2*const01;
  double const12 = const02;
  double const13 = -5*one_over_chi2;
  double const14 = 1*one_over_chi3;
  double const15 = -9*one_over_chi5;

  SymmetricTensor<2,3> hessY;
  for (unsigned a=0; a<3; ++a)
  for (unsigned b=a; b<3; ++b)
    hessY[a][b] = const11*S[a]*(a==b)
                  + const12*((a==b) - one_third + const13*S[a]*S[b])
                  + const14*(S[b] + S[a])
                  + const15*(S[a]*S[b]*S[b] + S[a]*S[a]*S[b]);

  double three_cot3theta = 3*(cos3theta/sin3theta);
  SymmetricTensor<2,3> hessTheta = dTheta_dY*hessY - three_cot3theta*gradTheta_dyad_gradTheta;
  */

  // Nodargi and Bisegna (2016)
  // (Note: it is singular when theta is 0 or pi/3)
  const double perturb = 1e-8;
  if (fabs(Theta - one_third*M_PI) < perturb)
    Theta = one_third*M_PI - perturb;
  else if (fabs(Theta) < perturb)
    Theta = perturb;

  SymmetricTensor<2,3> hessTheta;
  hessTheta[0][0] = two_thirds*one_over_chi2*std::sin(2*Theta_plus_23pi);
  hessTheta[1][1] = two_thirds*one_over_chi2*std::sin(2*Theta_minus_23pi);
  hessTheta[2][2] = two_thirds*one_over_chi2*std::sin(2*Theta);
  hessTheta[0][1] = hessTheta[2][2];
  hessTheta[0][2] = hessTheta[1][1];
  hessTheta[1][2] = hessTheta[0][0];
}
