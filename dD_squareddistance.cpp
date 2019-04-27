#include <iostream>
#include <cstdlib>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Segment_d<K> Segment;

int main()
{
  Point a(0,0), b(10,10);
  std::cout << K::Squared_distance_d()(a,b);
}