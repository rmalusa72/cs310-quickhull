#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Segment_d<K> Segment;
typedef std::vector<Point> vector;
typedef K::Equal_d equal;

const int dim = 3;

struct Vertex_point_item{
    template<class CMap>
    struct Dart_wrapper{
        typedef CGAL::Cell_attribute<CMap, Point> Vertex_attribute;
        typedef CGAL::cpp11::tuple<Vertex_attribute, void, void, void> Attributes;
    };
};

typedef CGAL::Combinatorial_map<dim, Vertex_point_item> CMap; 
typedef CMap::Dart_handle Dart_handle; 
CMap cm; 

void sew_on_match(Dart_handle dh1, Dart_handle dh2){
  for (CMap::One_dart_per_incident_cell_range<dim-2, dim-1>::iterator
       it=cm.one_dart_per_incident_cell<dim-2,dim-1>(dh1).begin(),
       itend=cm.one_dart_per_incident_cell<dim-2,dim-1>(dh1).end(); it!=itend; ++it)
  { 
    for (CMap::One_dart_per_incident_cell_range<dim-2, dim-1>::iterator
         it2=cm.one_dart_per_incident_cell<dim-2,dim-1>(dh2).begin(),
         itend2=cm.one_dart_per_incident_cell<dim-2,dim-1>(dh2).end(); it2!=itend2; ++it2)
    { 

      bool ridge_match = true;
      Dart_handle dart1 = it;
      Dart_handle dart2; 
      // ITERATE THROUGH THE VERTICES OF EACH I-1 CELL
      vector vertices;
      // Add all vertices from first i-1-cell to vertices
      // TODO: RIGHT NOW THIS IS ADDING ONLY ONE VERTEX
      std::cout<<"Vertices"<<std::endl;

      if(dim>3){
        for(CMap::One_dart_per_incident_cell_range<0, dim-2>::iterator it3 = cm.one_dart_per_incident_cell<0, dim-2>(it).begin(), itend3 = cm.one_dart_per_incident_cell<0, dim-2>(it).end(); it3!=itend3; ++it3){
          std::cout << cm.info<0>(it3) << std::endl;
          vertices.push_back(cm.info<0>(it3));
        }        
      }

      // Check that vertices of second i-1-cell are all in there
      // TODO: Is there a possibility they will have different lengths? 
      if(dim==3){
        if(equal()(cm.info<0>(dart1), cm.info<0>(cm.beta(it2,1))) && equal()(cm.info<0>(it2), cm.info<0>(cm.beta(dart1,1)))){
          dart2 = it2; 
        } else {
          ridge_match = false;
        }
      } else {
        for(CMap::One_dart_per_incident_cell_range<0, dim-2>::iterator it4 = cm.one_dart_per_incident_cell<0, dim-2>(it2).begin(), itend4 = cm.one_dart_per_incident_cell<0, dim-2>(it2).end(); it4!=itend4; ++it4){
          if (std::find(vertices.begin(), vertices.end(), cm.info<0>(it4)) == vertices.end()){
            ridge_match = false;
            // Exit this check
          }
          if(equal()(cm.info<0>(cm.beta(it4,1)), cm.info<0>(dart1)) && equal()(cm.info<0>(cm.beta(dart1, 1)), cm.info<0>(it4))){
            dart2 = it4; 
          }
        }
      }

      std::cout << (ridge_match?"yes":"no") << std::endl;

      // If they do match, find aligned darts, sew, and exit this loop
      if(ridge_match){
        // Edge1 and edge2 originate at the same vertex; need to find the compatriot of edge2 that is the opposite of edge1
        // Iterate through darts that share that vertex 

        std::cout << "dart1:" <<cm.info<0>(dart1) << " " << cm.info<0>(cm.beta(dart1, 1)) << std::endl;
        std::cout << "dart2:" <<cm.info<0>(dart2) << " " << cm.info<0>(cm.beta(dart2, 1)) << std::endl;

        if(dim==3){
          std::cout << "What the fuck" << std::endl;
          cm.link_beta<2>(dart1, dart2);
          std::cout << "What the fuck" << std::endl;
        } else {
          std::cout << "What the fuck2" << std::endl;
          cm.sew<dim-1>(dart1, dart2);
        }
        std::cout << "what even" << std::endl;
        return;
      }    
    }
  }
}

int main(){
  Dart_handle dh1 = cm.make_combinatorial_polygon(3);
  Dart_handle dh2 = cm.make_combinatorial_polygon(3);

  double origin[3] = {0,0,0};
  double coords[6][3] = {{1,1,1}, {2,3,5}, {0,10,5}, {6,7,7}, {15,20,10}, {80, 90, -10}};

  Point o = Point(3, origin, &origin[3]);
  Point p1 = Point(3, coords[0], &coords[0][3]);
  Point p2 = Point(3, coords[1], &coords[1][3]);
  Point p3 = Point(3, coords[2], &coords[2][3]);
  Point p4 = Point(3, coords[3], &coords[3][3]);
  Point p5 = Point(3, coords[4], &coords[4][3]);
  Point p6 = Point(3, coords[5], &coords[5][3]);
  Point points[6] = {p1, p2, p3, p4, p3, p2};
  
 // 1) Give every dart empty attributes
  for (CMap::Dart_range::iterator
       it=cm.darts().begin(), itend=cm.darts().end();
       it!=itend; ++it)
  {
    if ( cm.attribute<0>(it)==NULL )
      cm.set_attribute<0>(it, cm.create_attribute<0>());
  }

  // 2) Assign each vertex a point
  cm.info<0>(dh1) = points[0];
  cm.info<0>(cm.beta(dh1,1))=points[1];
  cm.info<0>(cm.beta(dh1,0))=points[2];
  cm.info<0>(dh2) = points[3];
  cm.info<0>(cm.beta(dh2,1))=points[4];
  cm.info<0>(cm.beta(dh2,0))=points[5];

  sew_on_match(dh1, dh2);
  std::cout << "is happening" << std::endl;
}



