#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/Cartesian.h>
#include <list>
#include <vector>

const int dim = 3; 

typedef CGAL::Cartesian<double> K;
typedef CGAL::Vector_3<K> Vector;
typedef CGAL::Segment_3<K> Segment;
typedef K::Compute_squared_distance_3 Squared_distance;
typedef CGAL::Orientation Orientation;
typedef K::Coplanar_3 not_independent;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<dim, dim, CGAL::Linear_cell_complex_traits<dim, K>> LCC;
typedef LCC::Dart_handle                                 Dart_handle;
typedef LCC::Point                                       Point;
typedef std::vector<Point> vector;
typedef std::list<Point> list; 

// typedef CGAL::Cartesian_d<double> K;
// typedef CGAL::Vector_d<K> Vector;
// typedef CGAL::Segment_d<K> Segment;
// typedef K::Squared_distance_d Squared_distance;
// typedef K::Orientation_d Orientation;
// typedef K::Affinely_independent_d independent;
// typedef CGAL::Linear_cell_complex_for_combinatorial_map<dim, dim, CGAL::Linear_cell_complex_traits<dim, K>> LCC;
// typedef LCC::Dart_handle                                 Dart_handle;
// typedef LCC::Point                                       Point;
// typedef std::vector<Point> vector;
// typedef std::list<Point> list; 

const LCC lcc; 

void quickhull(vector points);
vector find_initial_points(vector points);

int main(){
  // Load points into a list
  double origin[3] = {0,0,0};
  double coords[6][3] = {{1,1,1}, {2,3,5}, {0,10,5}, {6,7,7}, {15,20,10}, {80, 90, -10}};
  Point o = Point(0,0,0);
  Point p1 = Point(1,1,1);
  Point p2 = Point(2,3,5);
  Point p3 = Point(0,10,5);
  Point p4 = Point(6,7,7);
  Point p5 = Point(15,-20,10);
  Point p6 = Point(80, 90, -10);
  vector points;
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  points.push_back(p4);
  points.push_back(p5);
  points.push_back(p6);

  // Call quickhull
  quickhull(points);
}

// Given a set of d-dimensional points, constructs their convex hull in lcc
void quickhull(vector points){
  // Find initial (independent, ideally extreme) points
  vector extremes = find_initial_points(points);
  // Construct simplex
  for(int i=0; i<extremes.size(); i++){
    std::cout << extremes[i] << std::endl;
  }

  // Store darts from each face in list of faces

  // Sort other points into outside sets

  // Iterate thru face list until it is empty
}

// Given a list of d-dimensional points,
// return d+1 linearly independent points, ideally extreme 
vector find_initial_points(vector points){
  vector extremes;
  double extreme_coordinates[2*dim];
  Point curr; 
  vector::iterator it = points.begin();
  for(int i=0; i<dim; i++){
    extremes.push_back(*it);
    extremes.push_back(*it);
    extreme_coordinates[2*i]=(*it)[i];
    extreme_coordinates[2*i+1]=(*it)[i]; 
  }

  // get extremes in various dimensions
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

  vector initials;
  initials.push_back(extremes[0]);
  for(int i=1; i<extremes.size(); i++){
    if(std::find(initials.begin(), initials.end(), extremes[i])==initials.end()){
        initials.push_back(extremes[i]);
        if(initials.size()==3 && collinear(initials[0], initials[1], initials[2])){
            initials.pop_back();
        }
        if(initials.size()==4 && coplanar(initials[0], initials[1], initials[2], initials[3])){
            initials.pop_back();
        }
        if(initials.size()==4){
            return initials;
        }
    }

  }

  for(int i=0; i<points.size(); i++){
    if(std::find(initials.begin(), initials.end(), points[i])==initials.end()){
        initials.push_back(points[i]);
        if(initials.size()==3 && collinear(initials[0], initials[1], initials[2])){
            initials.pop_back();
        }
        if(initials.size()==4 && coplanar(initials[0], initials[1], initials[2], initials[3])){
            initials.pop_back();
        }
        if(initials.size()==4){
            return initials;
        }
    }
  }

  return initials;
}