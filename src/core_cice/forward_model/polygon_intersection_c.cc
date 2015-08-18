/* -*- C++ -*- */

#include "polygon_intersection_c.hh"
// #include <core/core> WHL - replaced by the line below
#include "./core/core" 
#include <vector>

#ifdef UNDERSCORE
#define polygon_intersection_c polygon_intersection_c_
#else
#ifdef DOUBLEUNDERSCORE
#define polygon_intersection_c polygon_intersection_c__
#endif
#endif

extern "C" {

void polygon_intersection_c(
    double *x1,
    double *y1,
    int    *n1p,
    double *x2,
    double *y2,
    int    *n2p,
    double *xInt,
    double *yInt,
    int    *nIntMaxp
    )
{
    typedef core::point<double> Point;
    typedef std::vector<Point> Vec_Point;

    int &n1      = *n1p;
    int &n2      = *n2p;
    int &nIntMax = *nIntMaxp;

    //std::cout << "n1 "      << n1 << std::endl;
    //std::cout << "n2 "      << n2 << std::endl;
    //std::cout << "nIntMax " << nIntMax << std::endl;

    Vec_Point p1(n1), p2(n2);

    for ( int i = 0; i < n1; ++i ) {
        p1[i].x = x1[i];
        p1[i].y = y1[i];
    }

    for ( int i = 0; i < n2; ++i ) {
        p2[i].x = x2[i];
        p2[i].y = y2[i];
    }

    Vec_Point intersect;

    // for now, call core's version with the default parameter values
    core::polygon_intersection(p1, p2, intersect);

    if ( intersect.size() > nIntMax ) {
      std::core_cerr <<
          "[polygon_intersection_c] ERROR: not enough room for "
          "intersection\npolygon; have " << intersect.size() << "points, but "
          "you gave space for only " << nIntMax << " points"
                     << std::endl;
      exit(1);
    }

    nIntMax = intersect.size();

    for ( int i = 0; i < nIntMax; ++i ) {
        xInt[i] = intersect[i].x;
        yInt[i] = intersect[i].y;
    }
}

} // extern "C"
