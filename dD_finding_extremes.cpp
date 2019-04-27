#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Vector_d<K> Vector;
typedef CGAL::Segment_d<K> Segment;
//typedef CGAL::Hyperplane_d<K> Hyperplane;
typedef K::Squared_distance_d Squared_distance;
typedef K::Orientation_d Orientation;
typedef std::list<Point> list; 
typedef K::Affinely_independent_d independent;
typedef std::vector<Point> vector;

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

  list points = list();
  points.push_back(a1);
  points.push_back(b);
  points.push_back(c);
  points.push_back(o);

  int dim = 3;

  vector extremes;
  double extreme_coordinates[2*dim];

  Point curr; 
  list::iterator it = points.begin();
  for(int i=0; i<dim; i++){
    extremes.push_back(*it);
    extremes.push_back(*it);
    extreme_coordinates[2*i]=(*it)[i];
    extreme_coordinates[2*i+1]=(*it)[i]; 
  }

  // get extremes in various dimensions (we have 2d and need d+1) and then make sure they're linearly independent
  for(++it; it != points.end(); ++it){
    curr = *it;
    for(int i=0; i<dim; i++){
      if(curr[i] <= extreme_coordinates[2*i]){
        extremes[2*i] = curr;
        extreme_coordinates[2*i] = curr[i];
      }
      if(curr[i] >= extreme_coordinates[2*i+1]){
        extremes[2*i+1]=curr;
        extreme_coordinates[2*i+1] = curr[i];
      }
    }
  }

  // Deal with possibility that they are not independent. Later
  if(!(independent()(extremes.begin(), extremes.end()))){
    std::cout << "error";
    exit(0);
  }


}