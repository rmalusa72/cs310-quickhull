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
typedef CGAL::Point_3<K> Point;
typedef K::Compute_squared_distance_3 Squared_distance;
typedef K::Orientation_3 Orientation;
typedef K::Coplanar_3 not_independent;
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

struct Lcc_attributes_normal
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point<Refs, void> Point_attribute;
    typedef CGAL::Cell_attribute<Refs, CGAL::cpp11::tuple<Plane, p_vector*>> Facet_attribute;
    typedef CGAL::cpp11::tuple<Point_attribute, void, Facet_attribute, void> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_for_combinatorial_map<dim, dim, CGAL::Linear_cell_complex_traits<dim, K>, Lcc_attributes_normal> LCC;
typedef LCC::Dart_handle                                 Dart_handle;
typedef std::list<Dart_handle> dart_list;
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
p_vector get_ridge_vertices(Dart_handle d);
Dart_handle get_matching_dart(Dart_handle d1, Dart_handle d2);

int main(){
  // Load points into a list
  Point o = Point(0,0,0);
  Point p1 = Point(0,-2,0);
  Point p2 = Point(0,0,5);
  Point p3 = Point(0,5,0);
  Point p4 = Point(5,0,0);
  Point p5 = Point(-5,1,1);
  //Point p6 = Point(80, 90, -10);
  p_vector points;
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  points.push_back(p4);
  points.push_back(p5);
  //points.push_back(p6);

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

  for(int i=0; i<extremes.size(); i++){
    std::cout << extremes[i] << std::endl;
  }

  lcc.make_tetrahedron(extremes[0], extremes[1], extremes[2], extremes[3]);

  // Create & associate 2-attributes to all darts
  for (LCC::Dart_range::iterator
       it=lcc.darts().begin(), itend=lcc.darts().end();
       it!=itend; ++it)
  {
    if ( lcc.attribute<dim-1>(it)==NULL )
      lcc.set_attribute<dim-1>(it, lcc.create_attribute<dim-1>());
  }

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
        std::cout << "switching plane" << std::endl; 
        f.plane = f.plane.opposite(); // Remember to free this later
    }
    facets.push_back(f);
    // Associate facet plane to this facet
    lcc.info<dim-1>(it) = f.plane;
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

  facet curr_facet; 
  // Iterate thru face list until it is empty
  while(facets.size() != 0){
    curr_facet = facets.front();
    facets.pop_front();
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
      dart_list visible; 
      dart_list to_visit;
      dart_list boundary; 
      dart_list new_facets;
      to_visit.push_back(curr_facet.handle);
      visible.push_back(curr_facet.handle);

      LCC::size_type m = lcc.get_new_mark(); 
      while(to_visit.size() != 0){
        // Pop next dart 
        Dart_handle curr = to_visit.front();
        to_visit.pop_front();
        std::cout << "CURR: " << std::endl;
        std::cout << lcc.point(curr) << std::endl;
        std::cout << lcc.point(lcc.beta(curr, 1)) << std::endl;

        // Mark darts of this cell visited
        // & Grab adjacent darts
        for(LCC::Dart_of_cell_range<dim-1>::iterator it = lcc.darts_of_cell<dim-1>(curr).begin(), itend = lcc.darts_of_cell<dim-1>(curr).end(); it != itend; ++it){
          lcc.mark(it, m);
          std::cout << "IT: " << std::endl;
          std::cout << lcc.point(it) << std::endl;
          std::cout << lcc.point(lcc.beta(it, 1)) << std::endl;
          
          // Dart in adjacent cell to current dart
          Dart_handle cross = lcc.beta(it, dim-1);
          std::cout << "CROSS: " << std::endl;
          std::cout << lcc.point(cross) << std::endl;
          std::cout << lcc.point(lcc.beta(cross, 1)) << std::endl;
          // If dart is marked - ignore
          // If dart plane is visible - add to to_visit
          // If dart plane is not visible - add this ridge to boundary
          if (!(lcc.is_marked(cross, m))){
            if(lcc.info<dim-1>(cross).oriented_side(max_p) == CGAL::ON_POSITIVE_SIDE){
              to_visit.push_back(cross);
              visible.push_back(cross);
            } else {
              boundary.push_back(cross);
              std::cout << "Pushing to boundary: " << lcc.point(cross) << std::endl;
            } // TODO: Account for coplanar case
          }
        }        
      }
      // Now boundary should contain one dart from each ridge on boundary, and visible should contain one dart from each visible facet
      //Free mark
      lcc.free_mark(m);

      // Join new point to boundary with new facets
      while (boundary.size() != 0){
        Dart_handle curr = boundary.front();
        boundary.pop_front(); 
        // build new simplex with vertices = points of curr and max_p, making sure it is oriented so it can be linked
        // Note: for shapes of higher dimension than triangles this may require more complex fiddling
        p_vector ridge_points = get_ridge_vertices(curr);
        for(int i=0; i<ridge_points.size(); i++){
          std::cout << "ridge points " << i << ":" << ridge_points[i] << std::endl;
        }

        Dart_handle new_facet = lcc.make_triangle(ridge_points[1], ridge_points[0], max_p);

        std::cout << "new_facet: " << lcc.point(new_facet) << std::endl;
        std::cout << "new_facet next: " << lcc.point(lcc.beta(new_facet, 1)) << std::endl;
        std::cout << "new_facet next next : " << lcc.point(lcc.beta(lcc.beta(new_facet, 1), 1)) << std::endl;
               
        // Unlink from current (visible) facet and link to new facet along the boundary ridge
        lcc.unsew<dim-1>(curr);
        lcc.sew<dim-1>(curr, new_facet);

        // Add to list of new facets for later linking
        new_facets.push_back(new_facet);

        // Find facet normal and add it to list, at some point
        // use any other point on hull to determine which side is positive and which side is negative (? check that this point works)
        facet f; 
        f.handle = new_facet; 
        f.vertices = get_cell_vertices(new_facet);
        f.plane = Plane(f.vertices[0], f.vertices[1], f.vertices[2]);
        Point other = lcc.point(lcc.beta(curr, 0));
        if(f.plane.oriented_side(other) != CGAL::ON_NEGATIVE_SIDE){
          std::cout << "switching plane" << std::endl;
          f.plane = f.plane.opposite(); // Remember to free this later
        }
        facets.push_back(f);

        // Create & associate 2-attributes to new darts
        for (LCC::Dart_of_cell_range::iterator
          n=lcc.darts_of_cell<dim-1>(new_facet).begin(), n_end=lcc.darts_of_cell<dim-1>(new_facet).end();
          n!=n_end; ++n)
        {
          if ( lcc.attribute<dim-1>(n)==NULL )
          lcc.set_attribute<dim-1>(n, lcc.create_attribute<dim-1>());
        }

        // Set those attributes to the plane of the new facet
        lcc.info<dim-1>(new_facet) = f.plane; 

      }

      // Glue together new facets along matching edges
      // This loop iterates through pairs of facets 
      for(dart_list::iterator it = new_facets.begin(), itend = new_facets.end(); it!=itend; ++it){
        dart_list::iterator it2 = it; 
        it2++; 
        for(dart_list::iterator itend2 = new_facets.end(); it2!=itend2; ++it2){
          bool match_found = false;
          Dart_handle match1;
          Dart_handle match2;
          std::cout << "PAIR: " << std::endl; 
          std::cout << "it1 1: " << lcc.point(*it) << std::endl;
          std::cout << "it1 2: " << lcc.point(lcc.beta(*it, 1)) << std::endl; 
          std::cout << "it2 1: " << lcc.point(*it2) << std::endl;
          std::cout << "it2 2: " << lcc.point(lcc.beta(*it2, 1)) << std::endl;

          // Now, for each pair of facets, iterate through each pair of ridges!!!!
          for(LCC::One_dart_per_incident_cell_range<dim-2,dim-1>::iterator r = lcc.one_dart_per_incident_cell<dim-2, dim-1>(*it).begin(), r_end = lcc.one_dart_per_incident_cell<dim-2,dim-1>(*it).end(); r!=r_end; ++r){
            for(LCC::One_dart_per_incident_cell_range<dim-2,dim-1>::iterator r1 = lcc.one_dart_per_incident_cell<dim-2, dim-1>(*it2).begin(), r1_end = lcc.one_dart_per_incident_cell<dim-2,dim-1>(*it2).end(); r1!=r1_end; ++r1){
              std::cout << "ridge 1:" << lcc.point(r) << "/" << lcc.point(lcc.beta(r, 1)) << std::endl;
              std::cout << "ridge 2:" << lcc.point(r1) << "/" << lcc.point(lcc.beta(r1, 1)) << std::endl;
              Dart_handle match = get_matching_dart(r, r1);
              if(match != r){
                match_found = true; 
                match1 = r; 
                match2 = match;
                break;
              }
            }
            if(match_found){
              break;
            }           
          }

          if(match_found){
            std::cout << "sewing along match : " << std::endl;
            std::cout << "match1: " << lcc.point(match1) << "/" << lcc.point(lcc.beta(match1, 1)) << std::endl;
            std::cout << "match2: " << lcc.point(match2) << "/" << lcc.point(lcc.beta(match2, 1)) << std::endl;
            lcc.sew<dim-1>(match1, match2);
          }

        }

      }

      // Resort outside sets of visible set & delete them
      // Emptying their outside sets will prevent us from processing them later, since we only
      // process facets with nonempty outside sets! 
      // ........... oh my god i need to associate outside sets with faces TOO???? 

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


    }
  }

  // Export finished hull as .OFF file 

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

p_vector get_ridge_vertices(Dart_handle handle){
  p_vector p; 
  if(dim == 3){
    // Ridge is a single edge; manually grab both points
    p.push_back(lcc.point(handle));
    p.push_back(lcc.point(lcc.beta(handle, 1)));
  } else if (dim > 3){
    for(LCC::Dart_of_cell_range<dim-2>::iterator it = lcc.darts_of_cell<dim-2>(handle).begin(), itend = lcc.darts_of_cell<dim-2>(handle).end(); 
        it != itend; ++it){
        p.push_back(lcc.point(it));
    }       
  }
  return p;
}

// Check if d1 and d2 are on ridges that match (have same length, same coords, opposite direction); if so, return the dart that 
// corresponds to d1 in d2's ridge
// This is arguably simpler in 3 and 4d because the ridges are lines and triangles, which have a pretty natural ordering
// Right now returning d1 if there is no match TODO: make this nicer
Dart_handle get_matching_dart(Dart_handle d1, Dart_handle d2){
  if(dim==3){
    Point p1 = lcc.point(d1);
    Point p2 = lcc.point(lcc.beta(d1, 1));
    if(lcc.point(d2) == p2 && lcc.point(lcc.beta(d2, 1)) == p1){
      return d2;
    }
  }
  return d1; 
}


