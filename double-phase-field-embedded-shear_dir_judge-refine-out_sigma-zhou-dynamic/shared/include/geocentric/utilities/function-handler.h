#ifndef FUNCTION_HANDLER_H
#define FUNCTION_HANDLER_H

#include <stdio.h>
#include <vector>
#include <string>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/mpi.h>
// #include <deal.II/base/mpi.templates.h>

using namespace dealii;


template <int dim>
class FunctionHandler
{
  public:
    FunctionHandler();

    void declare_parameters(ParameterHandler               &prm,
                            const std::vector<std::string> &user_tags,
                            const std::vector<unsigned>    &n_comp = std::vector<unsigned>());

    void initialize(ParameterHandler &prm);

    void print();

    void get(const std::string     &tag,
             std::vector<unsigned> &ids,
             FunctionParser<dim>   &function);

    void get_id(const std::string     &tag,
                std::vector<unsigned> &ids);

    double get_value(const std::string &tag,
                     const unsigned    &id,
                     const double      &now,
                     const Point<dim>  &pt,
                     const unsigned    &comp,
                     bool              &found_flag);

  private:
    struct Record
    {
      std::string           tag;
      std::string           ids_string;
      std::string           function_string;
      std::vector<unsigned> ids;
      unsigned              n_components;
    };

    std::vector<Record> records;

    ConditionalOStream pcout;
};


#endif

