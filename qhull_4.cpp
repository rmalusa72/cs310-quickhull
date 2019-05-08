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

// The LCC that is used for all of our operations
LCC lcc; 

// Global int used for writing repeated .off files because I am lazy
int facets_processed; 
int points_processed;
std::vector<double> origin_coords = {0,0,0,0};
Point origin = Point(4, origin_coords.begin(), origin_coords.end());

// Function headers TODO: move these to a .h file
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
    // Load points with vertices of a cube
    // TODO replace with hypercube
    /*points.push_back(Point(0,0,0,0));
    points.push_back(Point(0,0,5,0));
    points.push_back(Point(0,5,0,0));
    points.push_back(Point(5,0,0,0));
    points.push_back(Point(0,5,5,0));
    points.push_back(Point(5,0,5,0));
    points.push_back(Point(5,5,0,0));
    points.push_back(Point(5,5,5,0));
    */
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

    // In 3d version, the dimension value is ignored (as points are constructed from individual values, not arrays)
    int dim_1;
    in >> dim_1;
    int num_points;
    in >> num_points;
    // TODO: Check that this works 
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
  quickhull(points);
  std::cout << "Points processed: " << points_processed << std::endl;
}

// Given a set of d-dimensional points, this function constructs their convex hull
// in the linear cell complex lcc, and then writes it to a .off file. 
void quickhull(p_vector points){

  // Find initial points - independent, ideally extreme
  p_vector initial_points = find_initial_points(points);
  // Construct simplex 
  make_simplex(4, &initial_points);

  // PRINT each point of the initial set 
  for(int i=0; i<initial_points.size(); i++){
    std::cout << initial_points[i] << std::endl;
  }

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
  for(LCC::One_dart_per_cell_range<3>::iterator it=lcc.one_dart_per_cell<3>().begin(), itend = lcc.one_dart_per_cell<3>().end(); it!=itend; ++it){
    std::cout << "Handle: " << lcc.info<0>(it) << std::endl;
    std::cout << "points to " << lcc.info<0>(lcc.beta(it, 1)) << std::endl;
    std::cout << "points to " << lcc.info<0>(lcc.beta(lcc.beta(it, 1),1)) << std::endl;
    std::cout << "handle cross: " << lcc.info<0>(lcc.beta(it, 2)) << "/" << lcc.info<0>(lcc.beta(lcc.beta(it,2),1)) << "/" << lcc.info<0>(lcc.beta(lcc.beta(lcc.beta(it,2),1),1)) << std::endl;
  }

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

  // PRINT info on facets
  for(dart_list::iterator it = facets.begin(), itend = facets.end(); it!=itend; ++it){
    std::cout << (*(outside_set(*it))).size() << std::endl; 
    std::cout << (*(face_plane(*it))) << std::endl; 
  }

  // Iterate through list of facets to be processed until it is empty 
  facets_processed = 0; 
  points_processed = 5;
  write_off(); 
  while(facets.size() != 0){

    // PRINT info about current state of facets list 
    std::cout << "Num facets in list: " << facets.size() << std::endl;
    for(dart_list::iterator f = facets.begin(), fend = facets.end(); f!=fend; ++f){
      std::cout << "Facet handle: " << lcc.point(*f) << "/" << lcc.point(lcc.beta(*f, 1)) << "/" << lcc.point(lcc.beta(lcc.beta(*f, 1), 1))  << std::endl;
      std::cout << "To be deleted: " << get_deleted(*f) << std::endl; 
      std::cout << "Outside set size: " << (*(outside_set(*f))).size() << std::endl;
    }
    lcc.display_characteristics(std::cout);
    std::cout << std::endl; 

    // Pop first handle from facets list 
    Dart_handle curr_facet = facets.front();
    facets.pop_front();

    // Facet may be a remnant that has already been processed as part of a visible set
    // If so, its boolean get_deleted attribute will be true 
    if(get_deleted(curr_facet)){
      lcc.remove_cell<dim-1>(curr_facet);
      std::cout << "Removing previously processed facet!" << std::endl;
      facets_processed++;
      write_off();
      continue; 
    }

    // PRINT information about this facet's outside set 
    std::cout << "outside list size: " << (*(outside_set(curr_facet))).size() << std::endl;
    for(int i=0; i<(*(outside_set(curr_facet))).size(); i++){
      std::cout << i << ": " <<  (*(outside_set(curr_facet)))[i] << std::endl; 
    }

    // If facet has points in outside set, execute quickhull step 
    if ((*(outside_set(curr_facet))).size() != 0){
      Point furthest_p;
      if(select_furthest){
        furthest_p = get_furthest_point(curr_facet);
      } else {
        furthest_p = get_first_point(curr_facet);
      }
      
      points_processed++;

      // TODO: only mark darts of ridge, not facet, when dart is considered and is not visible - 
      // then boundary can be found on multiple ridges of that facet ! 

      std::cout << "BREADTH FIRST SEARCH:" << std::endl;

      dart_list visible; // Will contain a handle on each visible facet
      dart_list to_visit; // Used for breadth-first search
      dart_list boundary; // Will contain a handle on each ridge bordering the visible set

      to_visit.push_back(curr_facet);
      LCC::size_type m = lcc.get_new_mark(); 

      while(to_visit.size() != 0){

        std::cout << "VISIBLE:";
        for(dart_list::iterator vis = visible.begin(), vis_end = visible.end(); vis!=vis_end; ++vis){
          std::cout << lcc.info<0>(*vis) << "->" << lcc.info<0>(lcc.beta(*vis, 1)) << "->" << lcc.info<0>(lcc.beta(lcc.beta(*vis, 1),1)) << std::endl;
        }

        std::cout << "TO_VISIT:";
        for(dart_list::iterator tvis = to_visit.begin(), tvis_end = to_visit.end(); tvis!=tvis_end; ++tvis){
          std::cout << lcc.info<0>(*tvis) << "->" << lcc.info<0>(lcc.beta(*tvis, 1)) << "->" << lcc.info<0>(lcc.beta(lcc.beta(*tvis, 1),1)) << std::endl;
        }

        // Pop next dart 
        Dart_handle curr = to_visit.front();
        to_visit.pop_front();
        if(!lcc.is_marked(curr, m)){

          // Check if furthest_p can see curr
          if((*face_plane(curr)).oriented_side(furthest_p) == CGAL::ON_POSITIVE_SIDE || (*face_plane(curr)).oriented_side(furthest_p) == CGAL::ON_ORIENTED_BOUNDARY){
            std::cout << "Adding to visible set & pushing neighbors to visit" << std::endl;
            // If yes, push curr to visible set and add neighbors to to_visit
            // This will add duplicates, but each will only be processed once due to marks (hopefully?)
            // Mark every dart of facet visited
            visible.push_back(curr);
            for(LCC::Dart_of_cell_range<dim-1>::iterator it = lcc.darts_of_cell<dim-1>(curr).begin(), itend = lcc.darts_of_cell<dim-1>(curr).end(); it != itend; ++it){
              lcc.mark(it, m);
              Dart_handle cross = lcc.beta(it, dim-1); 
              to_visit.push_back(cross);
            }
          } else {
            std::cout << "Pushing to boundary" << std::endl;
            boundary.push_back(curr);
            // Mark every dart OF RIDGE visited
            // (Multiple ridges of one facet might be on boundary)
            for(LCC::Dart_of_cell_range<dim-2>::iterator it = lcc.darts_of_cell<dim-2>(curr).begin(), itend = lcc.darts_of_cell<dim-2>(curr).end(); it != itend; ++it){
              lcc.mark(it, m);
            }            
          }
        } else { 
          std::cout << "dart is already marked" << std::endl; 
        }    
      }
      lcc.free_mark(m);

      std::cout << "boundary:\n";
      for(dart_list::iterator b = boundary.begin(), b_end = boundary.end(); b!=b_end; ++b){
        std::cout << lcc.info<0>(*b) << "->" << lcc.info<0>(lcc.beta(*b, 1)) << "(->" << lcc.info<0>(lcc.beta(lcc.beta(*b, 1),1)) << ")" << std::endl; 
      }

      // Join new point to boundary with new facets
      dart_list new_facets;
      make_facet_cone(&boundary, &furthest_p, &new_facets);

      // Add new facets to facet list to be processed
      for(dart_list::iterator nf = new_facets.begin(), nf_end = new_facets.end(); nf!=nf_end; nf++){
        facets.push_back(*nf);
      }

      std::cout << "Number of visible facets: " << visible.size() << std::endl; 
      for(dart_list::iterator v = visible.begin(), v_end = visible.end(); v!=v_end; ++v){
        std::cout << lcc.info<0>(*v) << "->" << lcc.info<0>(lcc.beta(*v, 1)) << "->" << lcc.info<0>(lcc.beta(lcc.beta(*v, 1),1)) << std::endl; 
      }

      // Glue together new facets along matching edges
      glue_matching_facets(&new_facets);

      // Resort outside sets of visible set
      // For each facet in visible set 
      for(dart_list::iterator v = visible.begin(), v_end = visible.end(); v != v_end; ++v){
        std::cout << "Processing visible facet:" << lcc.point(*v) << std::endl; 
        // Pass furthest point as a point to ignore, as it has been processed and should not get re-sorted
        sort_into_outside_sets(outside_set(*v), &new_facets, true, furthest_p);
        // Set boolean in all visible facets to false, so they will be deleted when they are processed later
        set_deleted(*v, true);   
      }
      // Remove the facet we are processing directly, since it will not get re-added and processed
      lcc.remove_cell<dim-1>(curr_facet);
      lcc.display_characteristics(std::cout);
    }

    facets_processed++;
    write_off();

  }
  // Export finished hull as .OFF file 
  //write_off();
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
// TODO: Check that dD version works
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
    Point other = lcc.point(lcc.beta(lcc.beta(lcc.beta(it, dim-1),dim-2),0)); // TODO check that this actually gives dart outside facet
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
    // PRINT info about darts we are trying to sew
    
    lcc.unsew<dim-1>(curr);
    std::cout << "Sewing new facet in" << std::endl;
    lcc.sew<dim-1>(curr, new_facet);

    // Add to list of new facets for later linking
    (*new_facets).push_back(new_facet);

    // Find supporting plane 
    // Any other point on the hull can be used to determine which way the plane should point (as it is convex)
    p_vector vertices = get_cell_vertices(new_facet);
    Point other = lcc.point(lcc.beta(lcc.beta(lcc.beta(new_facet, dim-1),dim-2),0));
    std::cout << "vertices:";
    for(p_vector::iterator v = vertices.begin(), v_end = vertices.end(); v!=v_end; ++v){
      std::cout << *v << "/";
    } 
    std::cout << "\n" << "other:" << other;
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
// TODO: adapt for dD
void write_off(){
  std::ofstream hull_output; 
  hull_output.open("hull_output_" + std::to_string(facets_processed) + ".off"); /*+ std::to_string(facets_processed)*/
  hull_output << "OFF\n";
  
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
    std::cout << "face 2-free status:" << std::endl;
    if(!get_deleted(it)){
      num_faces++;
      vertices_of_face = 0; 
      for(LCC::Dart_of_orbit_range<1>::iterator it2 = lcc.darts_of_orbit<1>(it).begin(), itend2 = lcc.darts_of_orbit<1>(it).end(); 
        it2 != itend2; ++it2){
        vertices_of_face++; 
        face = face + std::to_string(lcc.info<0>(it2)) + " "; 
        std::cout << lcc.highest_nonfree_dimension(it2) << "/";
      }
      std::cout << std::endl; 
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



