#ifndef CONSTITUTIVE_MODEL_H
#define CONSTITUTIVE_MODEL_H

#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <geocentric/materials/lode-angle.h>
#include <geocentric/utilities/constants.h>
#include <geocentric/utilities/moduli-conversions.h>
#include <geocentric/utilities/eigen-decomposition.h>


using namespace dealii;

template <int dim>
class ConstitutiveModel
{
  public:
    ConstitutiveModel();

    std::string name();

  protected:
    std::string name_;
};


// Constructor
template <int dim>
ConstitutiveModel<dim>::ConstitutiveModel()
{}

// Methods
template <int dim>
std::string ConstitutiveModel<dim>::name()
{ return name_; }


#endif
