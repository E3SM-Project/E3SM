/*******************************************************************
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2010; all rights reserved
 *-----------------------------------------------------------------
 * This is the implementation file for an implementation
 * of the NVECTOR package, specifically suited to interface with
 * a external [Fortran] data structure and vector implementation.
 *
 * It contains the N_Vector kernels listed in nvector_external.h,
 * as well as implementations of the interfaces to the
 * create and destroy routines N_VNew_EXT and N_VDestroy_EXT.
 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "nvector_external.h"
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)


/********************* Fortran Interface Functions ***************/

/* define global N_Vector variables for Fortran interface */
N_Vector F2C_CVODE_vec;
N_Vector F2C_IDA_vec;
N_Vector F2C_KINSOL_vec;
N_Vector F2C_ARKODE_vec;


/* FNVEXT_INIT is the Fortran interface routine to initialize the
   template NVector */
void FNVEXT_INIT(int *code, int *ier)
{
  /* initialize return value to success */
  *ier = 0;

  /* Create solver-specific global NVector */
  /*   code = FCMIX_CVODE  (1) -> use with CVODE */
  /*   code = FCMIX_IDA    (2) -> use with IDA */
  /*   code = FCMIX_KINSOL (3) -> use with KINSOL */
  /*   code = FCMIX_ARKODE (4) -> use with ARKODE */
  /*   [others not implemented] */
  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_EXT();
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = NULL;
    F2C_IDA_vec = N_VNewEmpty_EXT();
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_EXT();
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_vec = NULL;
    F2C_ARKODE_vec = N_VNewEmpty_EXT();
    if (F2C_ARKODE_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


/********************* Exported Functions ************************/


/* Returns vector type ID.interface */
N_Vector_ID N_VGetVectorID_EXT(N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}


/* Function to create a new, empty EXT vector */
N_Vector N_VNewEmpty_EXT()
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_EXT content;

  /* Initialize vector pointers to NULL */
  v = NULL;
  ops = NULL;
  content = NULL;

  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) {free(v); return(NULL); }

  /* Attach custom vector routines to N_Vector_Ops structure */
  ops->nvgetvectorid     = N_VGetVectorID_EXT;
  ops->nvclone           = N_VClone_EXT;
  ops->nvcloneempty      = N_VCloneEmpty_EXT;
  ops->nvdestroy         = N_VDestroy_EXT;
  ops->nvspace           = N_VSpace_EXT;
  ops->nvgetarraypointer = N_VGetArrayPointer_EXT;
  ops->nvsetarraypointer = N_VSetArrayPointer_EXT;
  ops->nvlinearsum       = N_VLinearSum_EXT;
  ops->nvconst           = N_VConst_EXT;
  ops->nvprod            = N_VProd_EXT;
  ops->nvdiv             = N_VDiv_EXT;
  ops->nvscale           = N_VScale_EXT;
  ops->nvabs             = N_VAbs_EXT;
  ops->nvinv             = N_VInv_EXT;
  ops->nvaddconst        = N_VAddConst_EXT;
  ops->nvdotprod         = N_VDotProd_EXT;
  ops->nvmaxnorm         = N_VMaxNorm_EXT;
  ops->nvwrmsnorm        = N_VWrmsNorm_EXT;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_EXT;
  ops->nvmin             = N_VMin_EXT;
  ops->nvwl2norm         = N_VWL2Norm_EXT;
  ops->nvl1norm          = N_VL1Norm_EXT;
  ops->nvcompare         = N_VCompare_EXT;
  ops->nvinvtest         = N_VInvTest_EXT;
  ops->nvconstrmask      = N_VConstrMask_EXT;
  ops->nvminquotient     = N_VMinQuotient_EXT;

  /* Create content */
  content =
    (N_VectorContent_EXT) malloc(sizeof(struct _N_VectorContent_EXT));
  if (content == NULL) {free(ops); free(v); return(NULL);}

  /* Initialize content structure members */
  content->data     = NULL;
  content->own_data = SUNFALSE;

  /* Attach content and ops to generic N_Vector */
  v->content = content;
  v->ops     = ops;

  return(v);
}



/* N_VMake_EXT (or nvmake) creates a EXT N_Vector with 'content'
   pointing to user-provided structure */
N_Vector N_VMake_EXT(void* v_data)
{
  /* initialize v to NULL */
  N_Vector v = NULL;

  /* Create vector */
  v = N_VNewEmpty_EXT();
  if (v == NULL) return(NULL);

  /* Attach data if it is non-NULL */
  if ( v_data != NULL ) {
    NV_OWN_DATA_EXT(v) = SUNFALSE;
    NV_DATA_EXT(v)     = v_data;
  }

  return(v);
}



/* N_VPrint_EXT prints the N_Vector v to stdout.  This routine is
   provided to aid in debugging code using this vector package. */
void N_VPrint_EXT(N_Vector v)
{
  /* initialize vd to NULL */
  void *vd = NULL;

  /* access content field from vector; return if NULL */
  vd = NV_DATA_EXT(v);
  if (vd == NULL) return;

  /* call Fortran subroutine */
  FNVEXT_PRINT(vd);
}



/***************************************************************************/
/* BEGIN implementation of vector operations */

/* N_VCloneEmpty_EXT returns a new N_Vector of the same form as the
   input N_Vector, but with empty data container */
N_Vector N_VCloneEmpty_EXT(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_EXT content;

  /* initialize pointers to NULL */
  v = NULL;
  ops = NULL;
  content = NULL;

  /* Check that w has been created */
  if (w == NULL) return(NULL);

  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  /* Attach operations */
  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;
  ops->nvprod            = w->ops->nvprod;
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale;
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */
  content =
    (N_VectorContent_EXT) malloc(sizeof(struct _N_VectorContent_EXT));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Initialize content structure members */
  content->own_data = SUNFALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}



/* N_VClone_EXT returns a new N_Vector of the same form as the
   input N_Vector. */
N_Vector N_VClone_EXT(N_Vector w)
{
  N_Vector v;
  void *wdata, *vdata;
  int ierr;

  /* initialize pointers to NULL */
  v = NULL;
  wdata = NULL;
  vdata = (void **) malloc(sizeof(void *));

  /* Create vector */
  v = N_VCloneEmpty_EXT(w);
  if (v == NULL) return(NULL);

  /* Access data pointer for w */
  wdata = NV_DATA_EXT(w);

  /* Call Fortran routine to allocate data */
  FNVEXT_CLONE(wdata, vdata, &ierr);

  /* Fail gracefully if data is still NULL */
  if(ierr != 0) { N_VDestroy_EXT(v); return(NULL); }

  /* Attach data */
  NV_OWN_DATA_EXT(v) = SUNTRUE;
  NV_DATA_EXT(v)     = vdata;

  return(v);
}



/* N_VDestroy_EXT frees the data storage for a N_Vector */
void N_VDestroy_EXT(N_Vector v)
{
  void *vdata = NULL;
  if ( (NV_OWN_DATA_EXT(v) == SUNTRUE) && (NV_DATA_EXT(v) != NULL) ) {
    vdata = NV_DATA_EXT(v);
    FNVEXT_DESTROY(vdata);
    NV_DATA_EXT(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}



/* N_VSpace_EXT should return the space requirements for one N_Vector;
   however this routine is not utilized by SUNDIALS unless the user asks a
   particular integrator about its storage requirements.  As such, this is
   essentially a dummy routine, and so we return zero for both arguments. */
void N_VSpace_EXT(N_Vector v, long int* lrw, long int* liw)
{
  *lrw = ZERO;
  *liw = 0;
  return;
}



/* N_VGetArrayPointer_EXT (or nvgetarraypointer) extracts the
   data component array from the N_Vector v.  While the resulting pointer
   is cast to have 'realtype' type, this will just be re-cast to the
   appropriate pointer type by the user routines.

   **Warning: do not attempt to use SUNDIALS' direct linear solvers
     from this interface, since they will attempt to access specific
     entries of the data 'array', which will in actuality point to
     random memory locations in the Fortran vector data structure. */
realtype *N_VGetArrayPointer_EXT(N_Vector v)
{
  return((realtype *) NV_DATA_EXT(v));
}



/* N_VSetArrayPointer_EXT or (nvsetarraypointer) attaches the
   data component array v_data to the N_Vector v.  While the input
   pointer is assumed to have 'realtype' type, this will just be
   re-cast from the actual Fortran pointer type. */
void N_VSetArrayPointer_EXT(realtype* v_data, N_Vector v)
{
  /* if (v_data != NULL)  NV_DATA_EXT(v) = v_data; */
  if (v_data != NULL)  NV_DATA_EXT(v) = (void *) v_data;
}



/* N_VLinearSum_EXT (or nvlinearsum) calculates z = a*x + b*y */
void N_VLinearSum_EXT(realtype a, N_Vector x, realtype b,
		      N_Vector y, N_Vector z)
{
  /* extract data array handles from vectors */
  void *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_EXT(x);
  yd = NV_DATA_EXT(y);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_LINSUM(&a, xd, &b, yd, zd);

  return;
}



/* N_VConst_EXT (or nvconst) calculates z[i] = c for all i */
void N_VConst_EXT(realtype c, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *zd = NULL;
  zd = NV_DATA_EXT(z);

  /* call fortran routine to set z to c */
  FNVEXT_CONST(&c, zd);

  return;
}



/* N_VProd_EXT (or nvprod) calculates z[i] = x[i]*y[i] */
void N_VProd_EXT(N_Vector x, N_Vector y, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_EXT(x);
  yd = NV_DATA_EXT(y);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do product */
  FNVEXT_PROD(xd, yd, zd);

  return;
}



/* N_VDiv_EXT (or nvdiv) calculates z[i] = x[i]/y[i] */
void N_VDiv_EXT(N_Vector x, N_Vector y, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_EXT(x);
  yd = NV_DATA_EXT(y);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do division */
  FNVEXT_DIV(xd, yd, zd);

  return;
}



/* N_VScale_EXT (or nvscale) calculates z = c*x */
void N_VScale_EXT(realtype c, N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_SCALE(&c, xd, zd);

  return;
}



/* N_VAbs_EXT or (nvabs) calculates z[i] = |x[i]| */
void N_VAbs_EXT(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_ABS(xd, zd);

  return;
}



/* N_VInv_EXT (or nvinv) calculates z[i] = 1/x[i].
   Note: it does not check for division by 0.  It should be called only
   with an N_Vector x which is guaranteed to have all non-zero components. */
void N_VInv_EXT(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_INV(xd, zd);

  return;
}



/* N_VAddConst_EXT (or nvaddconst) calculates z[i] = x[i] + b */
void N_VAddConst_EXT(N_Vector x, realtype b, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_ADDCONST(&b, xd, zd);

  return;
}



/* N_VDotProd_EXT (or nvdotprod) returns the value of the
   ordinary dot product of x and y, i.e. sum (i=0 to N-1) {x[i] * y[i]} */
realtype N_VDotProd_EXT(N_Vector x, N_Vector y)
{
  /* extract data array handles from N_Vectors */
  void *xd, *yd;
  realtype sum = ZERO;
  xd = yd = NULL;
  xd = NV_DATA_EXT(x);
  yd = NV_DATA_EXT(y);

  /* call fortran routine to do operation */
  FNVEXT_DOTPROD(xd, yd, &sum);
  return(sum);
}



/* N_VMaxNorm_EXT (or nvmaxnorm) returns the maximum norm of x,
   i.e. max(i=1 to N-1) |x[i]|  */
realtype N_VMaxNorm_EXT(N_Vector x)
{
  /* extract data array handles from N_Vectors */
  realtype maxval = ZERO;
  void *xd  = NULL;
  xd = NV_DATA_EXT(x);

  /* call fortran routine to do operation */
  FNVEXT_MAXNORM(xd, &maxval);
  return(maxval);
}



/* N_VWrmsNorm_EXT (or nvwrmsnorm) returns the weighted root
   mean square norm of x with weight factor w,
   i.e. sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N] */
realtype N_VWrmsNorm_EXT(N_Vector x, N_Vector w)
{
  /* extract data array handles from N_Vectors */
  void *xd, *wd;
  realtype wrmsval = ZERO;
  xd = wd = NULL;
  xd = NV_DATA_EXT(x);
  wd = NV_DATA_EXT(w);

  /* call fortran routine to do operation */
  FNVEXT_WRMSNORM(xd, wd, &wrmsval);
  return(wrmsval);
}



/* N_VWrmsNormMask_EXT or (nvwrmsnormmask) returns wrms norm over
   indices indicated by id */
realtype N_VWrmsNormMask_EXT(N_Vector x, N_Vector w, N_Vector id)
{
  /* extract data array handles from N_Vectors */
  void *xd, *wd, *idd;
  realtype wrmsval = ZERO;
  xd = wd = idd = NULL;
  xd  = NV_DATA_EXT(x);
  wd  = NV_DATA_EXT(w);
  idd = NV_DATA_EXT(id);

  /* call fortran routine to do operation */
  FNVEXT_WRMSNORMMASK(xd, wd, idd, &wrmsval);
  return(wrmsval);
}



/* N_VMin_EXT (or nvmin) returns the smallest element of x */
realtype N_VMin_EXT(N_Vector x)
{
  /* extract data array handles from N_Vectors */
  realtype minval = ZERO;
  void *xd = NULL;
  xd = NV_DATA_EXT(x);

  /* call fortran routine to do operation */
  FNVEXT_MIN(xd, &minval);
  return(minval);
}



/* N_VWL2Norm_EXT (or nvwl2norm) returns the weighted
   Euclidean L2 norm of x with weight factor w,
   i.e. sqrt [(sum (i=0 to N-1) {(x[i]*w[i])^2}) ] */
realtype N_VWL2Norm_EXT(N_Vector x, N_Vector w)
{
  /* extract data array handles from N_Vectors */
  void *xd, *wd;
  realtype wl2val = ZERO;
  xd = wd = NULL;
  xd = NV_DATA_EXT(x);
  wd = NV_DATA_EXT(w);

  /* call fortran routine to do operation */
  FNVEXT_WL2NORM(xd, wd, &wl2val);
  return(wl2val);
}



/* N_VL1Norm_EXT (or nvl1norm) returns the L1 norm of x,
   i.e. sum (i=0 to N-1) {|x[i]|} */
realtype N_VL1Norm_EXT(N_Vector x)
{
  /* extract data array handles from N_Vectors */
  realtype l1norm = ZERO;
  void *xd = NULL;
  xd = NV_DATA_EXT(x);

  /* call fortran routine to do operation */
  FNVEXT_L1NORM(xd, &l1norm);
  return(l1norm);
}



/* N_VCompare_EXT (or nvcompare) calculates
   z[i] = 1 if |x[i]| > c, z[i] = 0 otherwise */
void N_VCompare_EXT(realtype c, N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_COMPARE(&c, xd, zd);

  return;
}



/* N_VInvTest_EXT (or nvinvtest) computes z[i] = 1/x[i]
   with a test for x[i] == 0 before inverting x[i].  This routine
   returns SUNTRUE if all components of x are nonzero (successful
   inversion) and returns SUNFALSE otherwise. */
booleantype N_VInvTest_EXT(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  void *xd, *zd;
  int testval = 0;
  xd = zd = NULL;
  xd = NV_DATA_EXT(x);
  zd = NV_DATA_EXT(z);

  /* call fortran routine to do operation */
  FNVEXT_INVTEST(xd, zd, &testval);

  if (testval == ZERO)  return(SUNTRUE);
  else  return(SUNFALSE);
}



/* N_VConstrMask_EXT (or nvconstrmask) returns a boolean SUNFALSE
   if any element fails the constraint test, and SUNTRUE if all passed.  The
   constraint test is as follows:
         if c[i] =  2.0, then x[i] must be >  0.0
         if c[i] =  1.0, then x[i] must be >= 0.0
         if c[i] = -1.0, then x[i] must be <= 0.0
         if c[i] = -2.0, then x[i] must be <  0.0
   It also sets a mask vector m, with elements equal to 1.0 where the
   corresponding constraint test failed, and equal to 0.0 where the
   constraint test passed.  This routine is specialized in that it is
   used only for constraint checking. */
booleantype N_VConstrMask_EXT(N_Vector c, N_Vector x, N_Vector m)
{
  /* extract data array handles from N_Vectors */
  int testval = 0;
  void *cd, *xd, *md;
  cd = xd = md = NULL;
  cd = NV_DATA_EXT(c);
  xd = NV_DATA_EXT(x);
  md = NV_DATA_EXT(m);

  /* call fortran routine to do operation */
  FNVEXT_CONSTRMASK(cd, xd, md, &testval);

  if (testval == ZERO)  return(SUNTRUE);
  else  return(SUNFALSE);
}



/* N_VMinQuotient_EXT (or nvminquotient) returns
   min(num[i]/denom[i]) over all i such that denom[i] != 0. */
realtype N_VMinQuotient_EXT(N_Vector num, N_Vector denom)
{
  /* extract data array handles from N_Vectors */
  void *nd, *dd;
  realtype minquot = ZERO;
  nd = dd = NULL;
  nd = NV_DATA_EXT(num);
  dd = NV_DATA_EXT(denom);

  /* call fortran routine to do operation */
  FNVEXT_MINQUOTIENT(nd, dd, &minquot);
  return(minquot);
}


/*********************** END OF FILE ***********************/
