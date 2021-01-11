#ifndef LINE_SEGMENT_H
#define LINE_SEGMENT_H 

// C++
#include <cmath>

// deal.II 
#include <deal.II/base/point.h>
#include <deal.II/base/tensor_function.h>

using namespace dealii; 

template <int dim> 
class LineSegment
{
  public:
  	LineSegment(); 
  	~LineSegment(); 

  	void input_points(const Point<dim> &pt1, const Point<dim> &pt2); 

    double projection(LineSegment<dim> &l); 
    bool   crack_side(const Point<dim> &pt); 

  	Point<dim>    start_point_output(); 
  	Point<dim>    end_point_output(); 
  	double        segment_length(); 
    Tensor<1,dim> direction_vector(); 

  private:
  	Point<dim>    start_point;
    Point<dim>    end_point; 
    double        length; 
    Tensor<1,dim> unit_vector; 
};

// -----------------------------------------------------------------------------
// CONSTRUCTOR/DESCRUCTOR
// -----------------------------------------------------------------------------
template <int dim>
LineSegment<dim>::LineSegment()
{}

template <int dim> 
LineSegment<dim>::~LineSegment()
{}

// -----------------------------------------------------------------------------
// INPUT FUNCTION 
// -----------------------------------------------------------------------------
template <int dim>
void LineSegment<dim>::input_points(const Point<dim> &pt1, const Point<dim> &pt2)
{
  start_point = pt1; 
  end_point   = pt2; 

  length = pt1.distance(pt2); 

  Tensor<1,dim> tmp_line; 
  for (unsigned d=0; d<dim; ++d)
    tmp_line[d] = pt2[d] - pt1[d]; 

  unit_vector = tmp_line/length; 
}

// -----------------------------------------------------------------------------
// PROJECTION ON LINE SEGMENT
// -----------------------------------------------------------------------------
template <int dim> 
double LineSegment<dim>::projection(LineSegment<dim> &l)
{
  double length_proj; 
  double L_length = l.segment_length(); 
  Tensor<1,dim> L_direction_vector = l.direction_vector(); 

  length_proj = L_length*fabs(L_direction_vector*unit_vector); 

  return length_proj; 
}

// -----------------------------------------------------------------------------
// JUDGING THE POINT ON MASTER SIDE OF SLAVE SIDE
// -----------------------------------------------------------------------------
template <int dim> 
bool LineSegment<dim>::crack_side(const Point<dim> &pt)
{
  bool side = false; // true for master side, false for slave side

  if (fabs(end_point[0] - start_point[0]) < tol)
  {
    if (pt[0] < end_point[0])
    {
      side = true; 
    }
  }
  else 
  {
    double slope = (end_point[1] - start_point[1])/(end_point[0] - start_point[0]); 

    double crack_y = slope*(pt[0] - start_point[0]) + start_point[1]; 

    if (pt[1] < crack_y)
    {
      side = true; 
    }
  }

  return side; 
}

// -----------------------------------------------------------------------------
// OUTPUT FUNCTION 
// -----------------------------------------------------------------------------
template <int dim>
Point<dim> LineSegment<dim>::start_point_output(){ return start_point; }

template <int dim> 
Point<dim> LineSegment<dim>::end_point_output(){ return end_point; }

template <int dim> 
double LineSegment<dim>::segment_length(){ return length; }

template <int dim> 
Tensor<1,dim> LineSegment<dim>::direction_vector(){ return unit_vector; }

#endif