#include <geocentric/utilities/function-handler.h>


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
FunctionHandler<dim>::FunctionHandler()
  :
  pcout(std::cout,Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{}


// -----------------------------------------------------------------------------
// DECLARE PARAMETERS
// -----------------------------------------------------------------------------
template <int dim>
void FunctionHandler<dim>::declare_parameters(ParameterHandler               &prm,
                                              const std::vector<std::string> &user_tags,
                                              const std::vector<unsigned>    &n_comp)
{
  // Modify the size of records vector to the number of tags
  records.resize(user_tags.size());

  // Declare functions and IDs for each tag
  for (unsigned c=0; c<records.size(); ++c)
  {
    records[c].tag = user_tags[c];

    if (n_comp.size() == 0)
      records[c].n_components = 1;
    else
      records[c].n_components = n_comp[c];

    prm.enter_subsection(records[c].tag);
      prm.declare_entry("Function","0.0",Patterns::Anything());
      prm.declare_entry("IDs","-1",Patterns::Anything());
    prm.leave_subsection();
  }
}


// -----------------------------------------------------------------------------
// INITIALIZE
// -----------------------------------------------------------------------------
template <int dim>
void FunctionHandler<dim>::initialize(ParameterHandler &prm)
{
  for (unsigned c=0; c<records.size(); ++c)
  {
    prm.enter_subsection(records[c].tag);

      records[c].function_string = prm.get("Function");
      records[c].ids_string      = prm.get("IDs");

      if (records[c].ids_string != "-1")
      {
        std::istringstream ids_stream(records[c].ids_string);

        while (true)
        {
          int n;
          ids_stream >> n;
          if (!ids_stream) break;
          records[c].ids.push_back(n);
        }

      }

    prm.leave_subsection();
  }
}


// -----------------------------------------------------------------------------
// PRINT
// -----------------------------------------------------------------------------
template <int dim>
void FunctionHandler<dim>::print()
{
  for (unsigned c=0; c<records.size(); ++c)
  if (records[c].ids.size() > 0)
  {
    pcout << "  " << records[c].tag << std::endl
          << "   Function   :: " << records[c].function_string << std::endl
          << "   IDs        :: " << records[c].ids_string << std::endl;

    if (records[c].n_components > 1)
    pcout << "   Components :: " << records[c].n_components << " (vector-valued)" << std::endl;
  }
}


// -----------------------------------------------------------------------------
// GET FUNCTIONS
// -----------------------------------------------------------------------------
template <int dim>
void FunctionHandler<dim>::get(const std::string     &tag,
                               std::vector<unsigned> &ids,
                               FunctionParser<dim>   &function)
{
  unsigned r = 100;
  for (unsigned c=0; c<records.size(); ++c)
    if (records[c].tag == tag)
      r = c;

  AssertThrow(r < records.size(),
              StandardExceptions::ExcMessage("Function tag not found"));

  ids = records[r].ids;

  std::string                  variable_names;
  std::map<std::string,double> constants;
  bool                         time_dependent = true;

  switch (dim)
  {
    case 1:
      variable_names = "x,t";
      break;

    case 2:
      variable_names = "x,y,t";
      break;

    case 3:
      variable_names = "x,y,z,t";
      break;
  }

  function.initialize(variable_names,
                      records[r].function_string,
                      constants,
                      time_dependent);
}


template <int dim>
void FunctionHandler<dim>::get_id(const std::string     &tag,
                                  std::vector<unsigned> &ids)
{
  unsigned r = 100;
  for (unsigned c=0; c<records.size(); ++c)
    if (records[c].tag == tag)
      r = c;

  AssertThrow(r < records.size(),
              StandardExceptions::ExcMessage("Function tag not found"));

  ids = records[r].ids;
}


template <int dim>
double FunctionHandler<dim>::get_value(const std::string  &tag,
                                       const unsigned     &id,
                                       const double       &now,
                                       const Point<dim>   &pt,
                                       const unsigned     &comp,
                                       bool               &found_flag)
{
  unsigned r = 100;
  for (unsigned c=0; c<records.size(); ++c)
    if (records[c].tag == tag)
      r = c;

  AssertThrow(r < records.size(),
              StandardExceptions::ExcMessage("Function tag not found"));

  found_flag = false;

  for (unsigned c=0; c<records[r].ids.size(); ++c)
  {
    if (records[r].ids[c] == id)
    {
      found_flag = true;
    }
  }

  if (found_flag == false)
    return 0.0;

  std::string                  variable_names;
  std::map<std::string,double> constants;
  bool                         time_dependent = true;

  switch (dim)
  {
    case 1:
      variable_names = "x,t";
      break;

    case 2:
      variable_names = "x,y,t";
      break;

    case 3:
      variable_names = "x,y,z,t";
      break;
  }

  FunctionParser<dim> function(records[r].n_components);

  function.initialize(variable_names,
                      records[r].function_string,
                      constants,
                      time_dependent);

  function.set_time(now);

  return function.value(pt,comp);
}


// EXPLICIT INSTANTIATIONS
template class FunctionHandler<1>;
template class FunctionHandler<2>;
template class FunctionHandler<3>;

