#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/draw_linear_cell_complex.h>


const int dim = 4;

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<4> LCC;
typedef LCC::Dart_handle                                 Dart_handle;
typedef LCC::Point                                       Point;

LCC lcc;

int main(){

  double origin[3] = {0,0,0};
  double coords[6][3] = {{1,1,1}, {2,3,5}, {0,10,5}, {6,7,7}, {15,20,10}, {80, 90, -10}};

  Point o = Point(3, origin, &origin[3]);
  Point p1 = Point(3, coords[0], &coords[0][3]);
  Point p2 = Point(3, coords[1], &coords[1][3]);
  Point p3 = Point(3, coords[2], &coords[2][3]);
  Point p4 = Point(3, coords[3], &coords[3][3]);
  Point p5 = Point(3, coords[4], &coords[4][3]);
  Point p6 = Point(3, coords[5], &coords[5][3]);

  Dart_handle dh1 = lcc.make_triangle(p1, p2, p3);
  Dart_handle dh2 = lcc.make_triangle(p2, p1, p4);
  lcc.sew<2>(dh1, dh2);
  Dart_handle dh3 = lcc.make_triangle(p3, p2, p4);
  lcc.sew<2>(lcc.beta(dh1, 1), dh3);
  lcc.sew<2>(lcc.beta(dh2, 0), lcc.beta(dh3, 1));
  Dart_handle dh4 = lcc.make_triangle(p4, p1, p3);
  lcc.sew<2>(lcc.beta(dh1, 0), lcc.beta(dh4,1));
  lcc.sew<2>(lcc.beta(dh2, 1), dh4);
  lcc.sew<2>(lcc.beta(dh3, 0), lcc.beta(dh4,0));
  for(LCC::Dart_of_cell_range<3>::iterator it = lcc.darts_of_cell<3>(dh1).begin(), itend = lcc.darts_of_cell<3>(dh1).end(); it!=itend; ++it){
    std::cout << lcc.point(it) << std::endl;
  }  

  std::cout << std::endl << std::endl;;

  Dart_handle dh5 = lcc.make_triangle(p3, p2, p1);
  Dart_handle dh6 = lcc.make_triangle(p2, p3, p5);
  lcc.sew<2>(dh5, dh6);
  Dart_handle dh7 = lcc.make_triangle(p1, p2, p5);
  lcc.sew<2>(lcc.beta(dh5, 1), dh7);
  lcc.sew<2>(lcc.beta(dh6, 0), lcc.beta(dh7, 1));
  Dart_handle dh8 = lcc.make_triangle(p5, p3, p1);
  lcc.sew<2>(lcc.beta(dh5, 0), lcc.beta(dh8,1));
  lcc.sew<2>(lcc.beta(dh6, 1), dh8);
  lcc.sew<2>(lcc.beta(dh7, 0), lcc.beta(dh8,0));
  for(LCC::Dart_of_cell_range<3>::iterator it = lcc.darts_of_cell<3>(dh5).begin(), itend = lcc.darts_of_cell<3>(dh5).end(); it!=itend; ++it){
    std::cout << lcc.point(it) << std::endl;
  }  

  lcc.sew3_same_facets();
  std::cout << std::endl; 

  for(LCC::Dart_of_cell_range<4>::iterator it = lcc.darts_of_cell<4>(dh1).begin(), itend = lcc.darts_of_cell<4>(dh1).end(); it!=itend; ++it){
     std::cout << lcc.point(it) << std::endl;
  }
}