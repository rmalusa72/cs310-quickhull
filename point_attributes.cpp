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

const int dim = 4;

struct Vertex_point_item{
    template<class CMap>
    struct Dart_wrapper{
        typedef CGAL::Cell_attribute<CMap, Point> Vertex_attribute;
        typedef CGAL::cpp11::tuple<Vertex_attribute, void, void, void> Attributes;
    };
};

typedef CGAL::Combinatorial_map<dim, Vertex_point_item> CMap; 
typedef CMap::Dart_handle Dart_handle; 

int main(){

  double origin[3] = {0,0,0};
  double coords1[3] = {1,1,1};
  double coords2[3] = {10,10,1};
  double coords3[3] = {20, 5, 1};
  Point a1 = Point(3, coords1, &coords1[3]);
  Point b = Point(3, coords2, &coords2[3]);
  Point c = Point(3, coords3, &coords3[3]);
  Point o = Point(3, origin, &origin[3]);
  Point points[4] = {a1,b,c,o};

  CMap cm; 
  cm.make_combinatorial_tetrahedron();

  // 1) Give every dart empty attributes
  for (CMap::Dart_range::iterator
       it=cm.darts().begin(), itend=cm.darts().end();
       it!=itend; ++it)
  {
    if ( cm.attribute<0>(it)==NULL )
      cm.set_attribute<0>(it, cm.create_attribute<0>());
  }

  // 2) Assign each vertex a point
  int i=0; 
  for(CMap::One_dart_per_cell_range<0>::iterator it = cm.one_dart_per_cell<0>().begin(),
                                                 itend = cm.one_dart_per_cell<0>().end(); 
      it!=itend; ++it){

    cm.info<0>(it) = points[i];
    i++;

  } 
}
