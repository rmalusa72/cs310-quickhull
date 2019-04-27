#include <iostream>
#include <cstdlib>
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
  Point a = Point(1,1,1, 1);
  Point b = Point(10,10,1, 1);
  Point c = Point(20,5,1, 1);
  Point d = Point(20,20,10, 1);
  Point e = Point(10, 10, 5, 1);
  Point f = Point(4, 4, 20, 1);
  Point pts[4] = {a,b,c,d};
  std::cout << Squared_distance()(a,b) << std::endl;
  std::cout << Orientation()(pts, &pts[4]);
  //std::cout << K::Hyperplane_d(2,pts,&pts[2]); This is trying to get coefficients u dumbass
  K::Hyperplane_d plane = K::Hyperplane_d(pts, &pts[2], pts[3]); // hyperplane thru first few points 
  //std::cout << Squared_distance()(plane, pts[3]) << std::endl; Can't automatically get signed distance to hypeprlane. SIGH
  Vector orthogonal = plane.orthogonal_vector(); // Get normal vector to plane
  std::cout << orthogonal.cartesian(0) << std::endl;
  std::cout << orthogonal.cartesian(1) << std::endl;
  std::cout << orthogonal.cartesian(2) << std::endl;
}