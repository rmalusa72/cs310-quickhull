#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <random>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Gmpzf.h>

const int dim = 4; 
bool select_furthest = true;

// Define kernel and its geometric objects
typedef CGAL::Gmpzf Gmpzf;
typedef CGAL::Homogeneous_d<Gmpzf> K;
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

// String class used in writing file
typedef std::string string;

void glue_matching_facets(dart_list* facets);
Dart_handle get_matching_dart(Dart_handle d1, Dart_handle d2);

// The LCC that is used for all of our operations
LCC lcc; 

int main(){
  exit(0);
}

// Given a pointer to a list containing handles on facets, glue them together
// where any two share a ridge with the same coordinates
void glue_matching_facets(dart_list* facets){
  // This loop iterates through pairs of facets 
  for(dart_list::iterator it = (*facets).begin(), itend = (*facets).end(); it!=itend; ++it){
    dart_list::iterator it2 = it; 
    it2++; 
    for(dart_list::iterator itend2 = (*facets).end(); it2!=itend2; ++it2){
      bool match_found = false;
      Dart_handle match1;
      Dart_handle match2;

      // Now, for each pair of facets, iterate through each pair of ridges
      for(LCC::One_dart_per_incident_cell_range<dim-2,dim-1>::iterator r = lcc.one_dart_per_incident_cell<dim-2, dim-1>(*it).begin(), r_end = lcc.one_dart_per_incident_cell<dim-2,dim-1>(*it).end(); r!=r_end; ++r){
        for(LCC::One_dart_per_incident_cell_range<dim-2,dim-1>::iterator r1 = lcc.one_dart_per_incident_cell<dim-2, dim-1>(*it2).begin(), r1_end = lcc.one_dart_per_incident_cell<dim-2,dim-1>(*it2).end(); r1!=r1_end; ++r1){
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
        std::cout << "Sewing new facets together" << std::endl;
        lcc.sew<dim-1>(match1, match2);
      }
    }
  }
}

// Given two dart_handles d1 and d2, check if they are on opposite ridges with matching coordinates, 
// and if so, return the dart in d2's ridge that corresponds to d1 and can be sewn to it. 
// If not, return d1. 
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