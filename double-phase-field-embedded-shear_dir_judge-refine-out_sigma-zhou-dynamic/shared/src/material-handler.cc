#include <geocentric/materials/material-handler.h>


// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESTRUCTOR
// -----------------------------------------------------------------------------
template <int dim, class Material>
MaterialHandler<dim,Material>::MaterialHandler
                              (const MPI_Comm &mpi_comm,
                               const PD::Triangulation<dim> &tria,
                               const Quadrature<dim> &quadrature)
                               :
                               mpi_comm(mpi_comm),
                               quadrature(quadrature),
                               max_material_id(20),
                               tria(tria)
{}

template <int dim, class Material>
MaterialHandler<dim,Material>::~MaterialHandler()
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim, class Material>
void MaterialHandler<dim,Material>::declare_parameters(ParameterHandler &prm)
{
  // Instantiate a material model
  Material material;

  // Declare parameters for all material ids
  for (unsigned m=0; m<=max_material_id; ++m)
  {
    prm.enter_subsection("Material_"+Utilities::int_to_string(m));
      material.declare_parameters(prm);
    prm.leave_subsection();
  }
}


// -----------------------------------------------------------------------------
// INITIALIZE
// -----------------------------------------------------------------------------
template <int dim, class Material>
void MaterialHandler<dim,Material>::initialize(ParameterHandler &prm)
{
  // Clear materials data (needed when re-initialize in adaptive remeshing)
  data.clear();

  // Instantiate a material model
  Material material;

  // Initialize an array for containing materials per cell
  unsigned n_q_points = quadrature.size();
  std::vector<Material> point_objects(n_q_points);

  // Set cell iterators
  typename PD::Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();

  // Loop over cells in the local triangulation
  for (; cell != endc; cell++)
  if (cell->is_locally_owned())
  {
    // Get material id of current cell
    unsigned m = cell->material_id();

    // Initialize material properties of current cell from parameter input
    prm.enter_subsection("Material_"+Utilities::int_to_string(m));
      material.initialize(prm);
    prm.leave_subsection();

    // Assign material objects per quadrature point
    for (unsigned q=0; q<n_q_points; ++q)
      point_objects[q] = material;

    // Save point objects (size q) per cell
    data.push_back(point_objects);
  }
}


// -----------------------------------------------------------------------------
// GET MATERIAL OBJECTS AT A SPECIFIC CELL
// -----------------------------------------------------------------------------
template <int dim, class Material>
std::vector<Material> MaterialHandler<dim,Material>::materials_at_cell(unsigned cell_index)
{
  return data.at(cell_index);
}


// -----------------------------------------------------------------------------
// ITERATORS
// -----------------------------------------------------------------------------
template <int dim, class Material>
typename MaterialHandler<dim,Material>::iterator MaterialHandler<dim,Material>::begin()
{
  return data.begin();
}

template <int dim, class Material>
typename MaterialHandler<dim,Material>::iterator MaterialHandler<dim,Material>::end()
{
  return data.end();
}


// EXPLICIT INSTANTIATIONS

// Basic materials
template class MaterialHandler<2,PorousMedia>;
template class MaterialHandler<3,PorousMedia>;

// Solid constitutive models
template class MaterialHandler<2,LinearElastic<2> >;
template class MaterialHandler<3,LinearElastic<3> >;

template class MaterialHandler<2,PhaseFieldLinearElastic<2> >;
template class MaterialHandler<3,PhaseFieldLinearElastic<3> >;

template class MaterialHandler<2,PhaseFieldCohesive<2> >;
template class MaterialHandler<3,PhaseFieldCohesive<3> >;

template class MaterialHandler<2,PhaseFieldContact<2> >;
template class MaterialHandler<3,PhaseFieldContact<3> >;

template class MaterialHandler<2,PhaseFieldMixMode<2> >; 
template class MaterialHandler<3,PhaseFieldMixMode<3> >; 

template class MaterialHandler<2,DruckerPrager<2> >; 
template class MaterialHandler<3,DruckerPrager<3> >; 

template class MaterialHandler<2,MohrCoulomb<2> >; 
template class MaterialHandler<3,MohrCoulomb<3> >; 