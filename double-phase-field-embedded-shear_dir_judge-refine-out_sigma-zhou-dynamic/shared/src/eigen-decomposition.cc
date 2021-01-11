#include <geocentric/utilities/eigen-decomposition.h>


// -----------------------------------------------------------------------------
// FUNCTION FOR EIGEN-DECOMPOSITION (SPECTRAL DECOMPOSITION)
// -----------------------------------------------------------------------------
template <int dim>
void eigen_decompose(const SymmetricTensor<2,dim> &input_tensor,
                     Tensor<1,dim>                &eigenvalues,
                     std::vector<Tensor<1,dim> >  &eigenvectors)
{
  // Copy input tensor to matrix of Eigen library
  Eigen::MatrixXd tensor(dim,dim);
  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
    tensor(i,j) = input_tensor[i][j];

  // Compute eigenvalues and eigenvectors
  // (Note: automatically sorted in an ascendeing order)
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(tensor);

  // Copy back to deal.II tensors
  for (unsigned i=0; i<dim; ++i)
  {
    eigenvalues[i] = es.eigenvalues()[i];
    for (unsigned j=0; j<dim; ++j)
      eigenvectors[i][j] = es.eigenvectors().col(i)[j];
  }
}


template <int dim>
void eigen_decompose(const Tensor<2,dim>          &input_tensor,
                     Tensor<1,dim>                &eigenvalues,
                     std::vector<Tensor<1,dim> >  &eigenvectors)
{
  // Copy input tensor to matrix of Eigen library
  Eigen::MatrixXd tensor(dim,dim);
  for (unsigned i=0; i<dim; ++i)
  for (unsigned j=0; j<dim; ++j)
    tensor(i,j) = input_tensor[i][j];

  // Compute eigenvalues and eigenvectors
  // (Note: automatically sorted in an ascendeing order)
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(tensor);

  // Copy back to deal.II tensors
  for (unsigned i=0; i<dim; ++i)
  {
    eigenvalues[i] = es.eigenvalues()[i];
    for (unsigned j=0; j<dim; ++j)
      eigenvectors[i][j] = es.eigenvectors().col(i)[j];
  }
}


// EXPLICIT INSTANTIATIONS
template
void eigen_decompose(const SymmetricTensor<2,2> &tensor,
                     Tensor<1,2>                &eigenvalues,
                     std::vector<Tensor<1,2> >  &eigenvectors);
template
void eigen_decompose(const SymmetricTensor<2,3> &tensor,
                     Tensor<1,3>                &eigenvalues,
                     std::vector<Tensor<1,3> >  &eigenvectors);

template
void eigen_decompose(const Tensor<2,2>          &tensor,
                     Tensor<1,2>                &eigenvalues,
                     std::vector<Tensor<1,2> >  &eigenvectors);
template
void eigen_decompose(const Tensor<2,3>          &tensor,
                     Tensor<1,3>                &eigenvalues,
                     std::vector<Tensor<1,3> >  &eigenvectors);
