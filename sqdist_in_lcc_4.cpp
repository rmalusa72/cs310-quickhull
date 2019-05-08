#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <random>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Linear_algebraHd.h>
#include <CGAL/constructions_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>

const int dim = 4; 
bool select_furthest = true;

// Define kernel and its geometric objects
typedef CGAL::Gmpzf Gmpzf;
typedef CGAL::Homogeneous_d<double> K;
typedef CGAL::Vector_d<K> Vector;
typedef CGAL::Segment_d<K> Segment;
typedef CGAL::Hyperplane_d<K> Plane;
typedef CGAL::Point_d<K> Point;

// Define vector of points for use in attributes
typedef std::vector<Point> p_vector;

// Define attributes of linear cell complex: 
// 0. a point, by default, plus an int used for indexing when outputting to .off 
// dim-1. (here, 2) a tuple:
//     * a plane (the supporting plane of the facet)
//     * a vector of points (the outside set of the facet)
//     * a boolean (whether the facet has been processed as part of a visible set and should be deleted)
struct Lcc_attributes_normal
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point<Refs, int> Point_attribute;
    typedef CGAL::Cell_attribute<Refs, std::tuple<Plane, p_vector*, bool>> Facet_attribute;
    typedef CGAL::cpp11::tuple<Point_attribute, void, void, Facet_attribute, void> Attributes;
  };
};

// Define linear cell complex with desired kernel and attributes
typedef CGAL::Linear_cell_complex_for_combinatorial_map<dim, dim, CGAL::Linear_cell_complex_traits<dim, K>, Lcc_attributes_normal> LCC;
typedef LCC::Dart_handle                                 Dart_handle;
typedef std::list<Dart_handle> dart_list;

double squared_distance(Point y1, Point q, Plane* plane_ptr);

// String class used in writing file
typedef std::string string;

int main(){
	std::vector<double> origin = {0,0,0,0};
  std::vector<double> coords1 = {1,1,1,1};
  std::vector<double> coords2 = {10,10,1,10};
  std::vector<double> coords3 = {20, 5, 1, 4};
  std::vector<double> coords4 = {7, -8, 1, -10};
  Point a1 = Point(4, coords1.begin(), coords1.end());
  Point b = Point(4, coords2.begin(), coords2.end());
  Point c = Point(4, coords3.begin(), coords3.end());
  Point d1 = Point(4, coords4.begin(), coords4.end());
  Point o = Point(4, origin.begin(), origin.end());
  Point points[4] = {a1,b,c,d1};
  K::Hyperplane_d plane = K::Hyperplane_d(points, &points[4], o, CGAL::ON_NEGATIVE_SIDE);

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      std::cout << points[i].cartesian(j) << " ";
    }
    std::cout << std::endl;
  }

  for(int i=0; i<5; i++){
    std::cout << plane.coefficient(i) << std::endl;
  }

  Vector orthogonal = plane.orthogonal_vector();
  for(int i=0; i<4; i++){
    std::cout << orthogonal.cartesian(i) << " " << std::endl;
  }

  // Finding closest point to an arbitrary point y
  double ycoords[4] = {20, 30, 40, 50};
  Vector y = Vector(4, ycoords, &ycoords[4]);
  Point y1 = Point(4, ycoords, &ycoords[4]);

  // p is a point on the hyperplane, but we want it as a vector
  Vector p = a1 - o;

  // a is the normal vector
  Vector a = orthogonal;

  // d = p dot a
  CGAL::Quotient<double> d = p*a;
  std::cout << to_double(d) << std::endl;

  // distance = y dot a - d / |a|
  //double squared_distance = to_double(pow(to_double(y*a-d), 2.0)/a.squared_length());
  //std::cout << squared_distance << std::endl;

  std::cout << "dist:" << std::to_string(squared_distance(y1, a1, &plane)) << std::endl; 
}

// Given a point y, a point q, and a hyperplane that q is on, 
// find the square of the distance from y to the hyperplane 
double squared_distance(Point y1, Point q, Plane* plane_ptr){
  std::vector<double> origin = {0,0,0,0};
  Point o = Point(4, origin.begin(), origin.end());

  Vector a = (*plane_ptr).orthogonal_vector();
  Vector p = q - o;
  Vector y = y1 - o;
  CGAL::Quotient<double> d = p*a;
  return(to_double(pow(to_double(y*a-d), 2.0)/a.squared_length()));
}