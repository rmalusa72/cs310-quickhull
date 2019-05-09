Compiling qhull_3.cpp or qhull_4.cpp requires CGAL and thus CGAL's dependencies, as well as gmp and mpfr. If the executables do not work, they should be compilable with 

c++ -std=c++11 -lcgal qhull_3.cpp -o qhull_3 -lboost_serialization -lgmp
c++ -std=c++11 -lcgal qhull_4.cpp -o qhull_4 -lboost_serialization -lgmp

after which they can be run with 

./qhull_3 (flags)
./qhull_4 (flags)

For an explanation of each flag, run 

./qhull_3 h