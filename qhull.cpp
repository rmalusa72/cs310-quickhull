#include <iostream>
#include <cstdlib>
#include <cmath>
#include <CGAL/Cartesian_d.h>
#include <CGAL/constructions_d.h>

typedef CGAL::Cartesian_d<double> K;
typedef CGAL::Point_d<K> Point;
typedef CGAL::Vector_d<K> Vector;
typedef CGAL::Segment_d<K> Segment;
typedef CGAL::Hyperplane_d<K> Hyperplane;
typedef K::Squared_distance_d Squared_distance;
typedef K::Orientation_d Orientation;

double point_plane_sqdistance(Hyperplane plane, Point on_plane, Point y){
    int dim = plane.dimension();

    double origin[dim];
    for(int i=0; i<dim; i++){
        origin[i] = 0;
    }
    Point o = Point(dim, origin, &origin[dim]);

    Vector a = plane.orthogonal_vector();
    Vector yv = y-o;
    Vector pv = on_plane-o;

    // d = p dot a
    double d = pv*a;
    // distance = y dot a - d / |a|
    double squared_distance = pow(yv*a - d,2.0) / a.squared_length();
    return squared_distance;
}

int main()
{

}

