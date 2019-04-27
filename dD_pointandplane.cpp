#include <iostream>
#include <cstdlib>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Segment_d<K> Segment;
typedef CGAL::Hyperplane_d<K> Hyperplane;

int main()
{
  Point a = Point(1,1,1, 1);
  Point b = Point(10,10,1, 1);
  Point c = Point(20,5,1, 1);
  Point d = Point(20,20,10, 1);
  Point e = Point(10, 10, 5, 1);
  Point f = Point(4, 4, 20, 1);
  Point pts[4] = {a,b,c,d};
  std::cout << K::Squared_distance_d()(a,b) << std::endl;
  std::cout << K::Orientation_d()(pts, &pts[4]);
  //std::cout << Hyperplane(3, pts, pts+3) << std::endl;

}