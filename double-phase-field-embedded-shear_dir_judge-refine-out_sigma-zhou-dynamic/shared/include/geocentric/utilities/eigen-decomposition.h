#ifndef EIGEN_DECOMPOSITION_H
#define EIGEN_DECOMPOSITION_H

#include <stdio.h>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <deal.II/base/symmetric_tensor.h>

using namespace dealii;

template <int dim>
void eigen_decompose(const SymmetricTensor<2,dim> &tensor,
                     Tensor<1,dim>                &eigenvalues,
                     std::vector<Tensor<1,dim> >  &eigenvectors);

template <int dim>
void eigen_decompose(const Tensor<2,dim>          &tensor,
                     Tensor<1,dim>                &eigenvalues,
                     std::vector<Tensor<1,dim> >  &eigenvectors);


#endif

