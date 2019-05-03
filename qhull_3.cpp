#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/enum.h>
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
typedef CGAL::Plane_3<K> Plane; 
typedef K::Compute_squared_distance_3 Squared_distance;
typedef K::Orientation_3 Orientation;
typedef K::Coplanar_3 not_independent;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<dim, dim, CGAL::Linear_cell_complex_traits<dim, K>> LCC;
typedef LCC::Dart_handle                                 Dart_handle;
typedef LCC::Point                                       Point;
typedef std::vector<Point> p_vector;

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

LCC lcc; 

class facet{
    public:
    Dart_handle handle;
    p_vector vertices;
    p_vector outside_set;
    Plane plane; 
};

typedef std::list<facet> facet_list; 

void quickhull(p_vector points);
p_vector find_initial_points(p_vector points);
p_vector get_cell_vertices(Dart_handle d);

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
  p_vector points;
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
void quickhull(p_vector points){
  // Find initial (independent, ideally extreme) points
  p_vector extremes = find_initial_points(points);
  // Construct simplex
  for(int i=0; i<extremes.size(); i++){
    std::cout << extremes[i] << std::endl;
  }
  lcc.make_tetrahedron(extremes[0], extremes[1], extremes[2], extremes[3]);

  // Store darts from each face in list of facets
  facet_list facets; 
  for(LCC::One_dart_per_cell_range<2>::iterator it = lcc.one_dart_per_cell<2>().begin(), itend = lcc.one_dart_per_cell<2>().end(); it != itend; ++it){
    facet f;
    f.handle = it; 
    f.vertices = get_cell_vertices(it);
    // Need to make sure this is the outward-facing plane - that the other hull point is outside it 
    // In dD can provide a point on other side, but in 3D have to do it ourselves by grabbing the point at the vertex of the simplex
    // that isn't in this face
    f.plane = Plane(f.vertices[0], f.vertices[1], f.vertices[2]);
    Point other = lcc.point(lcc.beta(lcc.beta(it, 2), 0));
    if(f.plane.oriented_side(other) != CGAL::ON_NEGATIVE_SIDE){
        std::cout << "switching plane"; 
        f.plane = f.plane.opposite(); // Remember to free this later
    }

    facets.push_back(f);
  }

  // Sort other points into outside sets
  // Make sure orientation here is correct
  for(int i=0; i<points.size(); i++){
    Point curr_point = points[i];
    if(std::find(extremes.begin(), extremes.end(), curr_point) == extremes.end()){
        std::cout << "curr point: " << curr_point << std::endl;
        for(facet_list::iterator it = facets.begin(), itend = facets.end(); it!=itend; ++it){
            if((*it).plane.oriented_side(curr_point) == CGAL::ON_POSITIVE_SIDE){
                std::cout << "point on pos side according to plane" << std::endl;
                std::cout << i << std::endl; 
                (*it).outside_set.push_back(curr_point);
                std::cout << "outside size: " << (*it).outside_set.size();
                break;
            }          
        }   
    }
  }

  for(facet_list::iterator it = facets.begin(), itend = facets.end(); it!=itend; ++it){
    std::cout << "vertices: ";
    for(int i=0; i < (*it).vertices.size(); i++){
        std::cout << (*it).vertices[i] << "/"; 
    }
    std::cout << std::endl << "outside set: ";
    for(int i = 0; i < (*it).outside_set.size(); i++){
        std::cout << (*it).outside_set[i] << "/";
    }
    std::cout << std::endl; 
  }

  // Iterate thru face list until it is empty
  while(facet_list.size() != 0){
    curr_facet = facet_list.pop_front();
    if (curr_facet.outside_set.size() != 0){
      // Find furthest point
      Point max_p = curr_facet.outside_set[0]; 
      double max_distance = Squared_distance()(max_p, curr_facet.plane);
      for(int i=1; i<curr_facet.outside_set.size(); i++){
        // For 4d will have to use my squared_distance and not the built-in one
        double curr_distance = Squared_distance()(curr_facet.outside_set[i], curr_facet.plane);
        if(curr_distance > max_distance){
          max_distance = curr_distance;
          max_p = curr_facet.outside_set[i];
        }
      }

      
      // Find visible set
      // Find boundary
      // Join new point to boundary with facets
      // Glue together facets 
      // Delete visible set (when you delete a dart, make sure not to try to process that face later)
    }
  }
}

// Given a list of d-dimensional points,
// return d+1 linearly independent points, ideally extreme 
p_vector find_initial_points(p_vector points){
  p_vector extremes;
  double extreme_coordinates[2*dim];
  Point curr; 
  p_vector::iterator it = points.begin();
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

  p_vector initials;
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
  // Points are coplanar
  std::cout << "Points are coplanar" << std::endl;
  exit(-1);
}

p_vector get_cell_vertices(Dart_handle handle){
    p_vector p; 
    for(LCC::Dart_of_cell_range<dim-1>::iterator it = lcc.darts_of_cell<dim-1>(handle).begin(), itend = lcc.darts_of_cell<dim-1>(handle).end(); 
        it != itend; ++it){
        p.push_back(lcc.point(it));
    }   
    return p; 
}


