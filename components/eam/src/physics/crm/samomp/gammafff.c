/*
 gamma-function for Fortran
 (C) Marat Khairoutdinov */

#include <math.h>
#include <stdio.h>

#ifdef __cplusplus 
extern "C" {
#endif

//float gammafff(float *x) {return (float)exp(lgamma(*x));}
double gammafff(double *x) {return (double)exp(lgamma(*x));}

//float gammafff_(float *x) {return (float)exp(lgamma(*x));}
double gammafff_(double *x) {return (double)exp(lgamma(*x));}

#ifdef __cplusplus
}
#endif
