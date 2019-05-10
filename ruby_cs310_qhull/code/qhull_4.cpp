#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <random>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>

const int dim = 4; 
bool select_furthest = true;

// Define kernel and its geometric objects
typedef CGAL::Homogeneous_d<double> K;
typedef CGAL::Vector_d<K> Vector;
typedef CGAL::Segment_d<K> Segment;
typedef CGAL::Hyperplane_d<K> Plane;
typedef CGAL::Point_d<K> Point;

// Define functions from kernel
typedef K::Affinely_independent_d aff_independent;

// Define vector of points for use in attributes
typedef std::vector<Point> p_vector;

// Define attributes of linear cell complex: 
// 0. a point, by default, plus an int used for indexing when outputting to .off 
// dim-1. (here, 3) a tuple:
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

// The LCC that is used for all of our operations
LCC lcc; 

std::vector<double> origin_coords = {0,0,0,0};
Point origin = Point(4, origin_coords.begin(), origin_coords.end());

// Function headers
void write_points(p_vector* plist_ptr);
void quickhull(p_vector points);
p_vector find_initial_points(p_vector points);
void make_simplex(int dim, p_vector* points_ptr);
void make_all_facets(dart_list* flist_ptr);
void sort_into_outside_sets(p_vector* plist_ptr, dart_list* flist_ptr, bool ignore_point=false, Point point_to_ignore=origin);
double squared_distance(Point y1, Point q, Plane* plane_ptr);
Point get_furthest_point(Dart_handle d);
Point get_first_point(Dart_handle d);
p_vector get_cell_vertices(Dart_handle d);
p_vector get_ridge_vertices(Dart_handle d);
void make_facet_cone(dart_list* boundary_ptr, Point* furthest_p_ptr, dart_list* new_facets);
void glue_matching_facets(dart_list* facets);
Dart_handle get_matching_dart(Dart_handle d1, Dart_handle d2);
Plane* face_plane(Dart_handle dh);
p_vector* outside_set(Dart_handle dh);
void set_deleted(Dart_handle dh, bool val); 
bool get_deleted(Dart_handle dh);
void write_off();
void write_cube_markers_off(p_vector* plist_ptr);

// Load points into a vector, and call quickhull on it
int main(int argc, char *argv[]){

  p_vector points;

  if(argc == 1){
    std::cout << "Please provide parameters:" << std::endl; 
    std::cout << "c to demonstrate a cube" << std::endl; 
    std::cout << "r followed by an integer n to demonstrate on n points uniformly distributed in a cube (or no number for default 100)" << std::endl; 
    std::cout << "f followed by absolute path to a file to read points from that file" << std::endl;
    std::cout << "h for information on reading files" << std::endl; 
    exit(0);
  } else if(argv[1][0] == 'c'){
    // Load points with vertices of a hypercube
    std::vector<double> coords1 = {0,0,0,0};
    std::vector<double> coords2 = {0,0,0,5};
    std::vector<double> coords3 = {0,0,5,0};
    std::vector<double> coords4 = {0,5,0,0};
    std::vector<double> coords5 = {0,5,0,5};
    std::vector<double> coords6 = {0,5,5,0};
    std::vector<double> coords7 = {0,0,5,5};
    std::vector<double> coords8 = {0,5,5,5};
    std::vector<double> coords9 = {5,0,0,0};
    std::vector<double> coords10 = {5,0,0,5};
    std::vector<double> coords11 = {5,0,5,0};
    std::vector<double> coords12 = {5,5,0,0};
    std::vector<double> coords13 = {5,5,0,5};
    std::vector<double> coords14 = {5,5,5,0};
    std::vector<double> coords15 = {5,0,5,5};
    std::vector<double> coords16 = {5,5,5,5};

    points.push_back(Point(4, coords1.begin(), coords1.end()));
    points.push_back(Point(4, coords2.begin(), coords2.end()));
    points.push_back(Point(4, coords3.begin(), coords3.end()));
    points.push_back(Point(4, coords4.begin(), coords4.end()));
    points.push_back(Point(4, coords5.begin(), coords5.end()));
    points.push_back(Point(4, coords6.begin(), coords6.end()));
    points.push_back(Point(4, coords7.begin(), coords7.end()));
    points.push_back(Point(4, coords8.begin(), coords8.end()));      
    points.push_back(Point(4, coords9.begin(), coords9.end()));
    points.push_back(Point(4, coords10.begin(), coords10.end()));
    points.push_back(Point(4, coords11.begin(), coords11.end()));
    points.push_back(Point(4, coords12.begin(), coords12.end()));
    points.push_back(Point(4, coords13.begin(), coords13.end()));
    points.push_back(Point(4, coords14.begin(), coords14.end()));
    points.push_back(Point(4, coords15.begin(), coords15.end()));
    points.push_back(Point(4, coords16.begin(), coords16.end()));      

  } else if(argv[1][0] == 'r'){
    int n;  
    if(argc >= 2){
      std::istringstream ss(argv[2]);
      if(!(ss >> n)){
        std::cerr << "Argument not accepted as integer\n";
        exit(1);
      }
    } else {
      n = 100;
    }

    // Load points with uniformly distributed points
    double lower_bound = -100;
    double upper_bound = 100;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re(seed);
    for(int i=0; i<n; i++){
      std::vector<double> coords = {unif(re), unif(re), unif(re), unif(re)};
      points.push_back(Point(dim, coords.begin(), coords.end()));
    }
  } else if(argv[1][0] == 'f'){
    // Read points from a file
    if(argc == 2){
      std::cout << "Please provide a file path as an argument" << std::endl;
      exit(0);
    } 
    std::ifstream in;
    in.open(argv[2]);
    if (!in) {
      std::cerr << "Unable to open file provided\n";
      exit(1); 
    }   

    int dim_1;
    in >> dim_1;
    int num_points;
    in >> num_points;
    double curr_coordinate;
    std::vector<double> coords;
    for(int i=0; i<num_points; i++){
      coords.clear();
      for(int j=0; j<dim_1; j++){
        in >> curr_coordinate;
        coords.push_back(curr_coordinate);
      }
      points.push_back(Point(dim, coords.begin(), coords.end()));
    }
  } else if(argv[1][0] == 'h'){
    std::cout << "Format for point list file:" << std::endl;
    std::cout << "First line: dimension value\n";
    std::cout << "Second line: number of points\n";
    std::cout << "Subsequent lines: space separated coordinates (decimals accepted)\n";
    std::cout << "Example file for a 3d cube provided in cube.txt\n";
  } else {
    std::cout << "Please provide parameters:" << std::endl; 
    std::cout << "c to demonstrate a cube" << std::endl; 
    std::cout << "r to demonstrate on points uniformly distributed in a cube" << std::endl; 
    std::cout << "f followed by a filename to read points" << std::endl;
    std::cout << "h for information on reading files" << std::endl;     
  }

  // Write points to a file for double checking with qhull
  write_points(&points);

  // The quickhull function constructs the hull in lcc
  // and writes the polytope to the hull_output.off file
  write_cube_markers_off(&points);
  quickhull(points);
}

// Given a set of d-dimensional points, this function constructs their convex hull
// in the linear cell complex lcc, and then writes it to a .off file. 
void quickhull(p_vector points){

  // Find initial points - independent, ideally extreme
  p_vector initial_points = find_initial_points(points);
  // Construct simplex 
  make_simplex(4, &initial_points);

  // Associate a new dim-1-attribute to each dart
  for (LCC::Dart_range::iterator
       it=lcc.darts().begin(), itend=lcc.darts().end();
       it!=itend; ++it)
  {
    if (lcc.attribute<dim-1>(it)==NULL){
      lcc.set_attribute<dim-1>(it, lcc.create_attribute<dim-1>());
    }
  }  

  write_off();

  // Initialize the list of facets to be processed with the newly created facets
  dart_list facets; 
  make_all_facets(&facets); 

  // Create a copy of the points vector without the extreme points 
  p_vector remaining_points;
  for(int i=0; i<points.size(); i++){
    if(std::find(initial_points.begin(), initial_points.end(), points[i]) == initial_points.end()){
      remaining_points.push_back(points[i]);
    }
  }

  // Sort the remaining points into outside sets 
  sort_into_outside_sets(&remaining_points, &facets);

  // Iterate through list of facets to be processed until it is empty 

  while(facets.size() != 0){

    // Pop first handle from facets list 
    Dart_handle curr_facet = facets.front();
    facets.pop_front();

    // Facet may be a remnant that has already been processed as part of a visible set
    // If so, its boolean get_deleted attribute will be true 
    if(get_deleted(curr_facet)){
      delete outside_set(curr_facet);
      lcc.remove_cell<dim-1>(curr_facet);
      continue; 
    }

    // If facet has points in outside set, execute quickhull step 
    if ((*(outside_set(curr_facet))).size() != 0){
      Point furthest_p;
      if(select_furthest){
        furthest_p = get_furthest_point(curr_facet);
      } else {
        furthest_p = get_first_point(curr_facet);
      }
      
      dart_list visible; // Will contain a handle on each visible facet
      dart_list to_visit; // Used for breadth-first search
      dart_list boundary; // Will contain a handle on each ridge bordering the visible set

      to_visit.push_back(curr_facet);
      LCC::size_type m = lcc.get_new_mark(); 

      while(to_visit.size() != 0){

        // Pop next dart 
        Dart_handle curr = to_visit.front();
        to_visit.pop_front();
        if(!lcc.is_marked(curr, m)){

          // Check if furthest_p can see curr
          if((*face_plane(curr)).oriented_side(furthest_p) == CGAL::ON_POSITIVE_SIDE || (*face_plane(curr)).oriented_side(furthest_p) == CGAL::ON_ORIENTED_BOUNDARY){
            // If yes, push curr to visible set and add neighbors to to_visit
            // This will add duplicates, but each will only be processed once due to marks
            // Mark every dart of facet visited
            visible.push_back(curr);
            for(LCC::Dart_of_cell_range<dim-1>::iterator it = lcc.darts_of_cell<dim-1>(curr).begin(), itend = lcc.darts_of_cell<dim-1>(curr).end(); it != itend; ++it){
              lcc.mark(it, m);
              Dart_handle cross = lcc.beta(it, dim-1); 
              to_visit.push_back(cross);
            }
          } else {
            boundary.push_back(curr);
            // Mark every dart OF RIDGE visited
            // (Multiple ridges of one facet might be on boundary)
            for(LCC::Dart_of_cell_range<dim-2>::iterator it = lcc.darts_of_cell<dim-2>(curr).begin(), itend = lcc.darts_of_cell<dim-2>(curr).end(); it != itend; ++it){
              lcc.mark(it, m);
            }            
          }
        }
      }
      lcc.free_mark(m);

      // Join new point to boundary with new facets
      dart_list new_facets;
      make_facet_cone(&boundary, &furthest_p, &new_facets);

      // Add new facets to facet list to be processed
      for(dart_list::iterator nf = new_facets.begin(), nf_end = new_facets.end(); nf!=nf_end; nf++){
        facets.push_back(*nf);
      }

      // Glue together new facets along matching edges
      glue_matching_facets(&new_facets);

      // Resort outside sets of visible set
      // For each facet in visible set 
      for(dart_list::iterator v = visible.begin(), v_end = visible.end(); v != v_end; ++v){
        // Pass furthest point as a point to ignore, as it has been processed and should not get re-sorted
        sort_into_outside_sets(outside_set(*v), &new_facets, true, furthest_p);
        // Set boolean in all visible facets to false, so they will be deleted when they are processed later
        set_deleted(*v, true);   
      }
      // Remove the facet we are processing directly, since it will not get re-added and processed
      delete outside_set(curr_facet);
      lcc.remove_cell<dim-1>(curr_facet);
    }
  }
  // Export finished hull as .OFF file 
  write_off();
}

void write_points(p_vector* plist_ptr){
  std::ofstream pt_output; 
  pt_output.open("input_pts.txt");
  pt_output << dim << std::endl; 
  pt_output << (*plist_ptr).size() << std::endl; 
  for(int i=0; i<(*plist_ptr).size(); i++){
    for(int j=0; j<dim; j++){
      pt_output << to_double((*plist_ptr)[i][j]) << " ";
    }
    pt_output << std::endl; 
  }
  pt_output.close();
}

// Given a list of d-dimensional points,
// return d+1 linearly independent points, ideally extreme 
p_vector find_initial_points(p_vector points){
  if(points.size()==0){
    std::cout << "No points provided!" << std::endl;
    std::exit(1);
  }
  p_vector extremes;
  CGAL::Quotient<double> extreme_coordinates[2*dim];
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
      if(!aff_independent()(initials.begin(), initials.end())){
        initials.pop_back();
      }
      if(initials.size()==dim+1){
        return initials; 
      }
    }
  }

  for(int i=0; i<points.size(); i++){
    if(std::find(initials.begin(), initials.end(), points[i])==initials.end()){
      initials.push_back(points[i]);
      if(!aff_independent()(initials.begin(), initials.end())){
        initials.pop_back();
      }
      if(initials.size()==dim+1){
        return initials; 
      }
    }
  }
  // Points are coplanar
  std::cout << "Points are coplanar or set is too small" << std::endl;
  exit(-1);
}

void make_simplex(int dim, p_vector* points_ptr){
  p_vector points = *points_ptr; 
  if(dim==3){
    lcc.make_tetrahedron(points[0], points[1], points[2], points[3]);
  } else if (dim==4){
    // Make a tetrahedron with each subset of four points, in matching orientations,
    // then glue together
    dart_list tetrahedra; 
    tetrahedra.push_back(lcc.make_tetrahedron(points[0], points[1], points[2], points[4]));
    tetrahedra.push_back(lcc.make_tetrahedron(points[1], points[0], points[3], points[4]));
    tetrahedra.push_back(lcc.make_tetrahedron(points[3], points[0], points[2], points[4]));
    tetrahedra.push_back(lcc.make_tetrahedron(points[1], points[0], points[2], points[3]));
    tetrahedra.push_back(lcc.make_tetrahedron(points[1], points[2], points[4], points[3]));
    glue_matching_facets(&tetrahedra);
  }
}

// For each facet in the lcc, associate its appropriate attributes;
// push a handle of each facet into the provided list
void make_all_facets(dart_list* flist_ptr){
  // Iterate over one dart of each of the dim-1-cells in the LCC
  for(LCC::One_dart_per_cell_range<dim-1>::iterator it = lcc.one_dart_per_cell<dim-1>().begin(), itend = lcc.one_dart_per_cell<dim-1>().end(); it != itend; ++it){
    p_vector vertices = get_cell_vertices(it);
    Point other = lcc.point(lcc.beta(lcc.beta(lcc.beta(it, dim-1),dim-2),0));
    Plane plane = Plane(vertices.begin(), vertices.end(), other, CGAL::ON_NEGATIVE_SIDE);

    // Attributes: the supporting plane, an empty outside set, and 
    // false to represent that the facet has not been processed and should not be deleted
    lcc.info<dim-1>(it) = std::make_tuple(plane, new p_vector(), false);
    (*flist_ptr).push_back(it);
  }

}

// Given a pointer to a list of points and a pointer to a list of facets, 
// sort the points into the outside sets of the facets
void sort_into_outside_sets(p_vector* plist_ptr, dart_list* flist_ptr, bool ignore_point/*=false*/, Point point_to_ignore/*=Point(0,0,0)*/){
  p_vector points = *plist_ptr;
  dart_list facets = *flist_ptr; 
  for(int i=0; i<points.size(); i++){
    Point curr_point = points[i];
    if(!ignore_point || (curr_point != point_to_ignore)){
      for(dart_list::iterator it = facets.begin(), itend = facets.end(); it!=itend; ++it){
        if((*face_plane(*it)).oriented_side(curr_point) == CGAL::ON_POSITIVE_SIDE || (*face_plane(*it)).oriented_side(curr_point) == CGAL::ON_ORIENTED_BOUNDARY){
            (*(outside_set(*it))).push_back(curr_point);
            break;
        }             
      }
    }
  }
}

// Given a point y, a point q, and a hyperplane that q is on, 
// find the square of the distance from y to the hyperplane 
double squared_distance(Point y1, Point q, Plane* plane_ptr){
  Vector a = (*plane_ptr).orthogonal_vector();
  Vector p = q - origin;
  Vector y = y1 - origin;
  CGAL::Quotient<double> d = p*a;
  return(to_double(pow(to_double(y*a-d), 2.0)/a.squared_length()));
}

// Iterate through d's outside set in search of the point furthest from d's supporting plane
Point get_furthest_point(Dart_handle d){
  Point furthest_p = (*(outside_set(d)))[0]; 
  CGAL::Quotient<double> max_distance = squared_distance(furthest_p, lcc.point(d), (face_plane(d)));
  for(int i=1; i<(*(outside_set(d))).size(); i++){
    // For 4d will have to use my squared_distance and not the built-in one
    CGAL::Quotient<double> curr_distance = squared_distance((*(outside_set(d)))[i], lcc.point(d), (face_plane(d)));
    if(curr_distance > max_distance){
      max_distance = curr_distance;
      furthest_p = (*(outside_set(d)))[i];
    }
  }
  return furthest_p; 
}

// Return the first point in d's outside set
// (Alternate selection step)
Point get_first_point(Dart_handle d){
  return (*(outside_set(d)))[0];
}

p_vector get_cell_vertices(Dart_handle handle){
    p_vector p; 
    for(LCC::Dart_of_cell_range<dim-1>::iterator it = lcc.darts_of_cell<dim-1>(handle).begin(), itend = lcc.darts_of_cell<dim-1>(handle).end(); 
        it != itend; ++it){
        if(std::find(p.begin(), p.end(), lcc.point(it)) == p.end()){
          p.push_back(lcc.point(it));         
        }
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
      if(std::find(p.begin(), p.end(), lcc.point(it)) == p.end()){
        p.push_back(lcc.point(it));         
      }
    }       
  }
  return p;
}

// Join point furthest_p to the ridges of the boundary list with a set of new facets; 
// put a handle to each new facet in the new_facets list.
// Note: leaves visible facets floating/unconnected
void make_facet_cone(dart_list* boundary_ptr, Point* furthest_p_ptr, dart_list* new_facets){
  dart_list boundary = *boundary_ptr; 
  Point furthest_p = *furthest_p_ptr; 
  while (boundary.size() != 0){
    Dart_handle curr = boundary.front();
    boundary.pop_front(); 
    // build new simplex with vertices = points of curr and furthest_p, making sure it is oriented so it can be linked
    // Note: for shapes of higher dimension than triangles this may require more complex fiddling
    // This version is written specifically for 4D
    p_vector ridge_points = get_ridge_vertices(curr);

    Dart_handle new_facet = lcc.make_tetrahedron(ridge_points[1], ridge_points[0], ridge_points[2], furthest_p);
           
    // Unlink from current (visible) facet and link to new facet along the boundary ridge
    lcc.unsew<dim-1>(curr);
    lcc.sew<dim-1>(curr, new_facet);

    // Add to list of new facets for later linking
    (*new_facets).push_back(new_facet);

    // Find supporting plane 
    // Any other point on the hull can be used to determine which way the plane should point (as it is convex)
    p_vector vertices = get_cell_vertices(new_facet);
    Point other = lcc.point(lcc.beta(lcc.beta(lcc.beta(new_facet, dim-1),dim-2),0));
    Plane new_facet_plane = Plane(vertices.begin(), vertices.end(), other, CGAL::ON_NEGATIVE_SIDE);

    // Create & associate dim-1-attributes to new darts
    for (LCC::Dart_of_cell_range<dim-1>::iterator
      n=lcc.darts_of_cell<dim-1>(new_facet).begin(), n_end=lcc.darts_of_cell<dim-1>(new_facet).end();
      n!=n_end; ++n)
    {
      if ( lcc.attribute<dim-1>(n)==NULL ){
        lcc.set_attribute<dim-1>(n, lcc.create_attribute<dim-1>());
      }
    }

    // Assign appropriate attributes (plane, empty outside set, not deleted yet)
    lcc.info<dim-1>(new_facet) = std::make_tuple(new_facet_plane, new p_vector(), false); 
  }
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
    // Looking for matching line segments
    Point p1 = lcc.point(d1);
    Point p2 = lcc.point(lcc.beta(d1, 1));
    if(lcc.point(d2) == p2 && lcc.point(lcc.beta(d2, 1)) == p1){
      return d2;
    }
  } else if (dim==4){
    // Looking for matching triangles
    for(int i=0; i<3; i++){
      if(lcc.point(d1) == lcc.point(lcc.beta(d2, 1)) && lcc.point(lcc.beta(d1, 1)) == lcc.point(d2) && lcc.point(lcc.beta(d1, 0)) == lcc.point(lcc.beta(d2, 0))){
        return d2; 
      }  
      d2 = lcc.beta(d2, 1);
    }
  }
  return d1; 
}

// Convenience function to grab the plane attribute of a facet
Plane* face_plane(Dart_handle dh){
  return &(std::get<0>(lcc.info<dim-1>(dh)));
}

// Convenience function to grab the outside set attribute of a facet
p_vector* outside_set(Dart_handle dh){
  return std::get<1>(lcc.info<dim-1>(dh));
}

void set_deleted(Dart_handle dh, bool val){
  std::get<2>(lcc.info<dim-1>(dh)) = val; 
}

bool get_deleted(Dart_handle dh){
  return std::get<2>(lcc.info<dim-1>(dh));
}

// Writes the contents of the linear cell complex to a .off file
void write_off(){
  std::ofstream hull_output; 
  hull_output.open("hull_output.off"); /*+ std::to_string(facets_processed)*/
  hull_output << "nOFF\n4\n";
  
  std::list<string> vertices;
  std::list<string> faces;

  int num_vertices = 0;
  Point p; 
  for(LCC::One_dart_per_cell_range<0>::iterator it = lcc.one_dart_per_cell<0>().begin(), itend = lcc.one_dart_per_cell<0>().end(); it!=itend; ++it){
    string vertex;
    if(!get_deleted(it)){
      lcc.info<0>(it) = num_vertices;
      p = lcc.point(it);
      for(int i=0; i<dim; i++){
        vertex = vertex + std::to_string(to_double(p[i])) + " ";
      }
      vertices.push_back(vertex);
      num_vertices++; 
    }
  }

  int num_faces = 0;
  int vertices_of_face;
  string face = "";
  for(LCC::One_dart_per_cell_range<2>::iterator it = lcc.one_dart_per_cell<2>().begin(), itend = lcc.one_dart_per_cell<2>().end(); it!=itend; ++it){
    if(!get_deleted(it)){
      num_faces++;
      vertices_of_face = 0; 
      for(LCC::Dart_of_orbit_range<1>::iterator it2 = lcc.darts_of_orbit<1>(it).begin(), itend2 = lcc.darts_of_orbit<1>(it).end(); 
        it2 != itend2; ++it2){
        vertices_of_face++; 
        face = face + std::to_string(lcc.info<0>(it2)) + " "; 
      }
      faces.push_back(std::to_string(vertices_of_face) + " " + face);
      face = "";
    }
  }

  hull_output << num_vertices << " " << num_faces << " 0\n";

  for(std::list<string>::iterator it = vertices.begin(), itend = vertices.end(); it!=itend; ++it){
    hull_output << *it << std::endl; 
  }

  for(std::list<string>::iterator it = faces.begin(), itend = faces.end(); it!=itend; ++it){
    hull_output << *it << std::endl; 
  }

  hull_output.close();
}

// Writes a hypercube (hyperrectangular prism?) around each point that extends along the x axis
void write_cube_markers_off(p_vector* plist_ptr){
  std::ofstream hull_output; 
  hull_output.open("hypercubes.off"); /*+ std::to_string(facets_processed)*/
  hull_output << "nOFF\n4\n";
  p_vector points = *plist_ptr;
  hull_output << (points.size()*16) << " " << (points.size()*24) << " 0" << std::endl;
  for(int i=0; i<points.size(); i++){
    
    Point p1 = points[i];
    std::vector<double> p;
    for(int j=0; j<dim; j++){
      p.push_back(to_double(p1[j]));
    }

    hull_output << p[0]+10 << " " << p[1] << " " << p[2] << " " << p[3] << std::endl; 
    hull_output << p[0]+10 << " " << p[1] + 5 << " " << p[2] << " " << p[3] << std::endl; 
    hull_output << p[0]+10 << " " << p[1] << " " << p[2]+5 << " " << p[3] << std::endl; 
    hull_output << p[0]+10 << " " << p[1] << " " << p[2] << " " << p[3]+5 << std::endl; 

    hull_output << p[0]+10 << " " << p[1] << " " << p[2]+5 << " " << p[3]+5 << std::endl; 
    hull_output << p[0]+10 << " " << p[1] + 5 << " " << p[2] << " " << p[3]+5 << std::endl; 
    hull_output << p[0]+10 << " " << p[1]+5 << " " << p[2]+5 << " " << p[3] << std::endl; 
    hull_output << p[0]+10 << " " << p[1] + 5 << " " << p[2]+5 << " " << p[3]+5 << std::endl; 
    
    hull_output << p[0]-10 << " " << p[1] << " " << p[2] << " " << p[3] << std::endl; 
    hull_output << p[0]-10 << " " << p[1] + 5 << " " << p[2] << " " << p[3] << std::endl; 
    hull_output << p[0]-10 << " " << p[1] << " " << p[2]+5 << " " << p[3] << std::endl; 
    hull_output << p[0]-10 << " " << p[1] << " " << p[2] << " " << p[3]+5 << std::endl; 

    hull_output << p[0]-10 << " " << p[1] << " " << p[2]+5 << " " << p[3]+5 << std::endl; 
    hull_output << p[0]-10 << " " << p[1] + 5 << " " << p[2] << " " << p[3]+5 << std::endl; 
    hull_output << p[0]-10 << " " << p[1]+5 << " " << p[2]+5 << " " << p[3] << std::endl; 
    hull_output << p[0]-10 << " " << p[1] + 5 << " " << p[2]+5 << " " << p[3]+5 << std::endl; 
  }

  for(int i=0; i<(16*points.size()); i+=16){
    hull_output << "4" << " " << i << " " << i+1 << " " << i+2 << " " << i+3 << std::endl; 
    hull_output << "4" << " " << i << " " << i+4 << " " << i+5 << " " << i+1 << std::endl; 
    hull_output << "4" << " " << i << " " << i+8 << " " << i+9 << " " << i+1 << std::endl; 
    hull_output << "4" << " " << i+1 << " " << i+2 << " " << i+10 << " " << i+9 << std::endl; 
    hull_output << "4" << " " << i+1 << " " << i+2 << " " << i+6 << " " << i+5 << std::endl; 
    hull_output << "4" << " " << i << " " << i+3 << " " << i+7 << " " << i+4 << std::endl; 
    hull_output << "4" << " " << i << " " << i+3 << " " << i+11 << " " << i+8 << std::endl; 
    hull_output << "4" << " " << i+3 << " " << i+11 << " " << i+10 << " " << i+2 << std::endl; 
    hull_output << "4" << " " << i+2 << " " << i+3 << " " << i+7 << " " << i+6 << std::endl; 
    hull_output << "4" << " " << i+2 << " " << i+6 << " " << i+14 << " " << i+10 << std::endl;
    hull_output << "4" << " " << i+7 << " " << i+15 << " " << i+14 << " " << i+6 << std::endl; 
    hull_output << "4" << " " << i+7 << " " << i+6 << " " << i+5 << " " << i+4 << std::endl; 
    hull_output << "4" << " " << i+7 << " " << i+3 << " " << i+11 << " " << i+15 << std::endl;  
    hull_output << "4" << " " << i+7 << " " << i+4 << " " << i+12 << " " << i+15 << std::endl; 
    hull_output << "4" << " " << i << " " << i+4 << " " << i+8 << " " << i+12 << std::endl;
    hull_output << "4" << " " << i+4 << " " << i+12 << " " << i+5 << " " << i+13 << std::endl;
    hull_output << "4" << " " << i+1 << " " << i+9 << " " << i+3 << " " << i+5 << std::endl;
    hull_output << "4" << " " << i+5 << " " << i+6 << " " << i+14 << " " << i+13 << std::endl;
    hull_output << "4" << " " << i+8 << " " << i+9 << " " << i+10 << " " << i+11 << std::endl;
    hull_output << "4" << " " << i+12 << " " << i+13 << " " << i+14 << " " << i+15 << std::endl;
    hull_output << "4" << " " << i+8 << " " << i+11 << " " << i+15 << " " << i+12 << std::endl;
    hull_output << "4" << " " << i+9 << " " << i+13 << " " << i+14 << " " << i+10 << std::endl;
    hull_output << "4" << " " << i+8 << " " << i+9 << " " << i+13 << " " << i+12 << std::endl;
    hull_output << "4" << " " << i+10 << " " << i+11 << " " << i+15 << " " << i+14 << std::endl;
  }
  hull_output << std::endl; 
  hull_output.close();
}

