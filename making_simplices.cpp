#include <CGAL/Combinatorial_map.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <list>
#include <vector>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>
//#include "making_simplices.h"

// Dimension is set as a constant right now because I need to parametrize CMap. 
const int dim = 4; 

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

typedef CGAL::Combinatorial_map<dim> CMap;
typedef CMap::Dart_handle Dart_handle; 

// Get rid of this forward declaration later and use header files, pls
Dart_handle make_simplex(int simplex_dimension, CMap cm);

int main(){
    CMap cm;
    make_simplex(4,cm);
}

Dart_handle make_simplex(int simplex_dimension, CMap cm){
    if (simplex_dimension == 2){
        return cm.make_combinatorial_polygon(3);
    }
    if (simplex_dimension == 3){
        return cm.make_combinatorial_tetrahedron();
    }
    if (simplex_dimension == 4){
        std::vector<Dart_handle> components;
        components.push_back(cm.make_combinatorial_tetrahedron());
        components.push_back(cm.make_combinatorial_tetrahedron());
        components.push_back(cm.make_combinatorial_tetrahedron());
        components.push_back(cm.make_combinatorial_tetrahedron());
        components.push_back(cm.make_combinatorial_tetrahedron());

        // Ideally should not write in all simplices by hand! Will not work for higher dimensions. But for now this is what's working. 
        // for(int i=0; i<simplex_dimension+1; i++){
        //     components.push_back(make_simplex(simplex_dimension-1, cm));
        // }
        std::cout << "test1" << std::endl;
        cm.sew<3>(components[0], components[1]); //0, 1
        std::cout << "test2" << std::endl;
        cm.sew<3>(cm.beta(components[0], 2), components[2]); //0,2
        std::cout << "test3" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[0], 1), 2), components[3]); //0,3
        std::cout << "test4" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[0], 0), 2), components[4]); //0,4
        std::cout << "test5" << std::endl;
        cm.sew<3>(cm.beta(components[2], 2), cm.beta(components[1], 2)); //2,1
        std::cout << "test6" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[2], 0), 2), cm.beta(cm.beta(components[4], 1), 2)); //2,4
        std::cout << "test7" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[2], 1), 2), cm.beta(cm.beta(components[3], 0), 2)); //2,3
        std::cout << "test8" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[1], 0), 2), cm.beta(components[3], 2)); //1,3
        std::cout << "test9" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[1], 1), 2), cm.beta(components[4], 2)); //1,4
        std::cout << "test10" << std::endl;
        cm.sew<3>(cm.beta(cm.beta(components[3], 1), 2), cm.beta(cm.beta(components[4], 0), 2)); //3, 4
        std::cout << "test10" << std::endl;
        CGAL_assertion(cm.is_valid());
        
    }
    return NULL; 

    // Building simplices by hand for now
    // else if (simplex_dimension >= 4 && simplex_dimension <= dim){
        
    //     // Make d+1 smaller dimensional simplices

    //     std::vector<Dart_handle> components;

    //     for(int i=0; i<=simplex_dimension; i++){
    //         components.push_back(make_simplex(simplex_dimension-1, cm));
    //     }
    //     // Sew them together
    //     // get handle on each, and iterate through, sewing first one to other d, second one to other d-1, etc? 
    //     // stepping iterator forward once each time so that facets are not reused 
    //     // what about face orientation..


    //     for (CMap::Dart_of_orbit_range<>::iterator it = cm.; i < count; ++i)
    //     {
    //         /* code */
    //     }

    // }
    return NULL; 
}