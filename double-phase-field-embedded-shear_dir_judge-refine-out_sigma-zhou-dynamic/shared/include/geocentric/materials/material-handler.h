#ifndef MATERIAL_HANDLER_H
#define MATERIAL_HANDLER_H

// Deal.II
#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/distributed/tria.h>

// Basic materials
#include <geocentric/materials/porous-media.h>

// Solid models
#include <geocentric/materials/linear-elastic.h>
#include <geocentric/materials/phasefield-linear-elastic.h>
#include <geocentric/materials/phasefield-cohesive.h>
#include <geocentric/materials/phasefield-contact.h>
#include <geocentric/materials/drucker-prager.h>
#include <geocentric/materials/mohr-coulomb.h>
#include <geocentric/materials/phasefield-mix-mode.h>


using namespace dealii;
namespace PD = parallel::distributed;


template <int dim, class Material>
class MaterialHandler
{
  public:
    MaterialHandler(const MPI_Comm &mpi_comm,
                    const PD::Triangulation<dim> &tria,
                    const Quadrature<dim> &quadrature);
    ~MaterialHandler();

    void declare_parameters(ParameterHandler &prm);
    void initialize(ParameterHandler &prm);
    void get_viz_data(DataOut<dim> &data_out);
    std::vector<Material> materials_at_cell(unsigned cell_index);

    typedef typename std::vector<std::vector<Material> >::iterator iterator;
    iterator begin();
    iterator end();

  private:
    const MPI_Comm         &mpi_comm;
    const Quadrature<dim>  &quadrature;
    const unsigned          max_material_id;

    const parallel::distributed::Triangulation<dim> &tria;

    std::vector<std::vector<Material> > data;
};


#endif
