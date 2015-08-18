
#ifndef INCLUDED_POLYGON_INTERSECTION_C
#define INCLUDED_POLYGON_INTERSECTION_C

#ifdef UNDERSCORE
#define polygon_intersection_c polygon_intersection_c_
#else
#ifdef DOUBLEUNDERSCORE
#define polygon_intersection_c polygon_intersection_c__
#endif
#endif

extern "C" {

/*
  Interface to core::polygon_intersection(), using bare-bones C-style
  arrays.

  x1       x-coordinates of polygon 1
  y1       y-coordinates of polygon 1
  n1p      number of points in polygon 1
  x2       x-coordinates of polygon 2
  y2       y-coordinates of polygon 2
  n2p      number of points in polygon 2
  xInt     returned x-coordinates of polygon intersection
  yInt     returned y-coordinates of polygon intersection
  nIntMaxp On input, max allowable points in the polygon intersection.
           On output, the number of points in the polygon intersection.
*/
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
    );

}

#endif // INCLUDED_POLYGON_INTERSECTION_C

