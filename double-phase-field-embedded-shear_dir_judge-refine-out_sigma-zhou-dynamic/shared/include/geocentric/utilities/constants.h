#ifndef GEOCENTRIC_CONSTANTS_H
#define GEOCENTRIC_CONSTANTS_H

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

static const double tol         = 1e-9;
static const double tol_pf      = 1e-4; 

static const double one_third   = 1./3.;
static const double two_thirds  = 2./3.;
static const double one_ninth   = 1./9.;
static const double two_ninths  = 2./9.;
static const double four_ninths = 4./9.;
static const double root13      = sqrt(1./3.);
static const double root23      = sqrt(2./3.);
static const double root32      = sqrt(3./2.);
static const double root6       = sqrt(6.);
static const double Pi          = 3.14159265; 

static const double p_atm       = 101.3;

static const SymmetricTensor<2,3> zero         = SymmetricTensor<2,3>();
static const SymmetricTensor<2,3> eye          = unit_symmetric_tensor<3>();
static const SymmetricTensor<2,2> eye_2d       = unit_symmetric_tensor<2>();
static const SymmetricTensor<4,3> big_eye      = identity_tensor<3>();
static const SymmetricTensor<4,3> eye_dyad_eye = outer_product(eye,eye);


#endif
