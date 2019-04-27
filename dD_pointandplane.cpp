#include <iostream>
#include <cstdlib>
#include <cmath>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Vector_d<K> Vector;
typedef CGAL::Segment_d<K> Segment;
//typedef CGAL::Hyperplane_d<K> Hyperplane;
typedef K::Squared_distance_d Squared_distance;
typedef K::Orientation_d Orientation;

int main()
{

  double origin[3] = {0,0,0};
  double coords1[3] = {1,1,1};
  double coords2[3] = {10,10,1};
  double coords3[3] = {20, 5, 1};
  Point a1 = Point(3, coords1, &coords1[3]);
  Point b = Point(3, coords2, &coords2[3]);
  Point c = Point(3, coords3, &coords3[3]);
  Point o = Point(3, origin, &origin[3]);
  Point points[3] = {a1,b,c};
  K::Hyperplane_d plane = K::Hyperplane_d(points, &points[3], o, CGAL::ON_NEGATIVE_SIDE);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      std::cout << points[i].cartesian(j);
    }
    std::cout << std::endl;
  }

  for(int i=0; i<4; i++){
    std::cout << plane.coefficient(i) << std::endl;
  }

  Vector orthogonal = plane.orthogonal_vector();
  for(int i=0; i<3; i++){
    std::cout << orthogonal.cartesian(i) << " " << std::endl;
  }

  // Finding closest point to an arbitrary point y
  double ycoords[3] = {20, 30, 40};
  Vector y = Vector(3, ycoords, &ycoords[3]);

  // p is a point on the hyperplane, but we want it as a vector
  Vector p = a1 - o;

  // a is the normal vector
  Vector a = orthogonal;

  // d = p dot a
  double d = p*a;
  std::cout << d << std::endl;

  // distance = y dot a - d / |a|
  double squared_distance = pow(y*a - d,2.0) / a.squared_length();
  std::cout << squared_distance << std::endl;

}