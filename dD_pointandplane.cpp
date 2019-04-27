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

  double origin[3] = {0,0,0};
  double coords1[3] = {1,1,1};
  double coords2[3] = {10,10,1};
  double coords3[3] = {20, 5, 1};
  Point a = Point(3, coords1, &coords1[3]);
  Point b = Point(3, coords2, &coords2[3]);
  Point c = Point(3, coords3, &coords3[3]);
  Point o = Point(3, origin, &origin[3]);
  Point points[3] = {a,b,c};
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


  // Point a = Point(1,1,1, 1);
  // Point b = Point(10,10,1, 1);
  // Point c = Point(20,5,1, 1);
  // Point d = Point(20,20,10, 1);
  // Point e = Point(10, 10, 5, 1);
  // Point f = Point(4, 4, 20, 1);
  // Point pts[4] = {a,b,c,d};
  // std::cout << Squared_distance()(a,b) << std::endl;
  // std::cout << Orientation()(pts, &pts[4]);
  // //std::cout << K::Hyperplane_d(2,pts,&pts[2]); This is trying to get coefficients u dumbass
  // K::Hyperplane_d plane = K::Hyperplane_d(pts, &pts[2], pts[3]); // hyperplane thru first few points 
  // //std::cout << Squared_distance()(plane, pts[3]) << std::endl; Can't automatically get signed distance to hypeprlane. SIGH
  // Vector orthogonal = plane.orthogonal_vector(); // Get normal vector to plane
  // std::cout << orthogonal.cartesian(0) << std::endl;
  // std::cout << orthogonal.cartesian(1) << std::endl;
  // std::cout << orthogonal.cartesian(2) << std::endl;
}