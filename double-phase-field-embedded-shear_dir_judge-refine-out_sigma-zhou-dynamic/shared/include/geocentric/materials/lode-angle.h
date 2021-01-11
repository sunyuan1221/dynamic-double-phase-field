#ifndef LODE_ANGLE_H
#define LODE_ANGLE_H

#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <geocentric/utilities/constants.h>

using namespace dealii;


// =============================================================================
// Gudehus--Argyris Scaling Function (Gudehus 1973, Argyris et al. 1974)
// =============================================================================
class GudehusArgyris
{
  public:
    GudehusArgyris();
    GudehusArgyris(double ellipticity);
    ~GudehusArgyris();

    void set_ellipticity(double ellipticity);

    double               value();
    Tensor<1,3>          grad();
    SymmetricTensor<2,3> hess();

    // Update scaling function according to current stress (Lode's angle)
    void update(const Tensor<1,3> &sig);

  private:
    // Ellipticity: should be in between 7/9(0.778) and 1
    double rho;

    // Scaling function and its derivatives
    double               Zeta;
    Tensor<1,3>          gradZeta;
    SymmetricTensor<2,3> hessZeta;
};


// =============================================================================
// Hyperbolic Mohr-Coulomb Approximation Function (Abbo and Sloan, 1995)
// =============================================================================
class HyperbolicMohrCoulomb
{
  public:
    HyperbolicMohrCoulomb();
    HyperbolicMohrCoulomb(double limit_lode_angle);
    ~HyperbolicMohrCoulomb();

    void set_limit_lode_angle(double limit_lode_angle);

    double               value(const double Angle);
    Tensor<1,3>          grad(const double Angle);
    SymmetricTensor<2,3> hess(const double Angle);

    double      value_deriv_angle(const double Angle);
    Tensor<1,3> grad_deriv_angle(const double Angle);

    // Update scaling function according to current stress (Lode's angle)
    void update(const Tensor<1,3> &sig);

  private:
    double Theta_t; // Limit load angle (should be close to 30 degrees)

    // Lode angle and its derivatives
    double Theta;
    Tensor<1,3> gradTheta;
    SymmetricTensor<2,3> hessTheta;
    SymmetricTensor<2,3> gradTheta_dyad_gradTheta;
};


// =============================================================================
// LODE ANGLE (common for many scaling functions)
// =============================================================================
class LodeAngle
{
  public:
    LodeAngle();
    LodeAngle(const Tensor<1,3> &sig);
    ~LodeAngle();

    double               value();
    Tensor<1,3>          grad();
    SymmetricTensor<2,3> hess();

    // Update Lode's angle according to current stress (Lode's angle)
    void update(const Tensor<1,3> &sig);

  private:
    // Lode angle and its derivatives
    double               Theta;
    Tensor<1,3>          gradTheta;
    SymmetricTensor<2,3> hessTheta;
};

#endif

