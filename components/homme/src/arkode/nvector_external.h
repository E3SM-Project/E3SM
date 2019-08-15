/*******************************************************************
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2010; all rights reserved
 *------------------------------------------------------------------
 * This is the header file for an implementation of an NVECTOR
 * package, specifically suited to interface with an external
 * [Fortran] implementation of the vector operations. In these 
 * vectors, local data consists of an arbitrary array of a given 
 * length.
 *
 * Part I of this file contains definitions and prototypes of the
 * Fortran interface functions for this implementation.
 *
 * Part II of this file contains declarations which are specific to 
 * this particular vector specification.  This includes the typedef 
 * for the 'content' field of the vector.
 *
 * Part III of this file defines accessor macros that allow the
 * user to efficiently use N_Vector type without making explicit 
 * references to its underlying representation.
 *
 * Part IV of this file contains the prototypes for
 * initialization and printing routines specific to this
 * implementation.
 *
 * Part V of this file contains prototypes for the vector kernels
 * which operate on the N_Vector.  These prototypes are unique to
 * this particular implementation of the vector package, and are
 * attached to the generic N_Vector structure for use within the
 * various SUNDIALS solvers.
 *
 * NOTES:
 *
 * The definitions of the generic N_Vector structure is in the
 * SUNDIALS header file nvector.h.
 *
 * The definition of the type realtype is in the SUNDIALS header
 * file sundialstypes.h and may be changed according to the user's
 * needs. The SUNDIALS file sundialstypes.h also contains the
 * definition for the type booleantype.
 *
 *******************************************************************/


#ifndef _NVECTOR_EXTERNAL_H
#define _NVECTOR_EXTERNAL_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif



/****************************************************************
 * PART I:
 * Fortran callable wrappers for the interface routines, as
 * well as prototypes of the supplied Fortran vector operations
 * to be called by C routines.
 ****************************************************************/

/* Define Fortran/C interface wrappers */
#if defined(SUNDIALS_F77_FUNC)

#define FNVEXT_INIT         SUNDIALS_F77_FUNC(fnvextinit,         FNVEXTINIT)
#define FNVEXT_PRINT        SUNDIALS_F77_FUNC(fnvextprint,        FNVEXTPRINT)
#define FNVEXT_CLONE        SUNDIALS_F77_FUNC(fnvextclone,        FNVEXTCLONE)
#define FNVEXT_DESTROY      SUNDIALS_F77_FUNC(fnvextdestroy,      FNVEXTDESTROY)
#define FNVEXT_LINSUM       SUNDIALS_F77_FUNC(fnvextlinearsum,    FNVEXTLINEARSUM)
#define FNVEXT_CONST        SUNDIALS_F77_FUNC(fnvextconst,        FNVEXTCONST)
#define FNVEXT_PROD         SUNDIALS_F77_FUNC(fnvextprod,         FNVEXTPROD)
#define FNVEXT_DIV          SUNDIALS_F77_FUNC(fnvextdiv,          FNVEXTDIV)
#define FNVEXT_SCALE        SUNDIALS_F77_FUNC(fnvextscale,        FNVEXTSCALE)
#define FNVEXT_ABS          SUNDIALS_F77_FUNC(fnvextabs,          FNVEXTABS)
#define FNVEXT_INV          SUNDIALS_F77_FUNC(fnvextinv,          FNVEXTINV)
#define FNVEXT_ADDCONST     SUNDIALS_F77_FUNC(fnvextaddconst,     FNVEXTADDCONST)
#define FNVEXT_DOTPROD      SUNDIALS_F77_FUNC(fnvextdotprod,      FNVEXTDOTPROD)
#define FNVEXT_MAXNORM      SUNDIALS_F77_FUNC(fnvextmaxnorm,      FNVEXTMAXNORM)
#define FNVEXT_WRMSNORM     SUNDIALS_F77_FUNC(fnvextwrmsnorm,     FNVEXTWRMSNORM)
#define FNVEXT_WRMSNORMMASK SUNDIALS_F77_FUNC(fnvextwrmsnormmask, FNVEXTWRMSNORMMASK)
#define FNVEXT_MIN          SUNDIALS_F77_FUNC(fnvextmin,          FNVEXTMIN)
#define FNVEXT_WL2NORM      SUNDIALS_F77_FUNC(fnvextwl2norm,      FNVEXTWL2NORM)
#define FNVEXT_L1NORM       SUNDIALS_F77_FUNC(fnvextl1norm,       FNVEXTL1NORM)
#define FNVEXT_COMPARE      SUNDIALS_F77_FUNC(fnvextcompare,      FNVEXTCOMPARE)
#define FNVEXT_INVTEST      SUNDIALS_F77_FUNC(fnvextinvtest,      FNVEXTINVTEST)
#define FNVEXT_CONSTRMASK   SUNDIALS_F77_FUNC(fnvextconstrmask,   FNVEXTCONSTRMASK)
#define FNVEXT_MINQUOTIENT  SUNDIALS_F77_FUNC(fnvextminquotient,  FNVEXTMINQUOTIENT)

#else

#define FNVEXT_INIT         fnvextinit_
#define FNVEXT_PRINT        fnvextprint_
#define FNVEXT_CLONE        fnvextclone_
#define FNVEXT_DESTROY      fnvextdestroy_
#define FNVEXT_LINSUM       fnvextlinearsum_
#define FNVEXT_CONST        fnvextconst_
#define FNVEXT_PROD         fnvextprod_
#define FNVEXT_DIV          fnvextdiv_
#define FNVEXT_SCALE        fnvextscale_
#define FNVEXT_ABS          fnvextabs_
#define FNVEXT_INV          fnvextinv_
#define FNVEXT_ADDCONST     fnvextaddconst_
#define FNVEXT_DOTPROD      fnvextdotprod_
#define FNVEXT_MAXNORM      fnvextmaxnorm_
#define FNVEXT_WRMSNORM     fnvextwrmsnorm_
#define FNVEXT_WRMSNORMMASK fnvextwrmsnormmask_
#define FNVEXT_MIN          fnvextmin_
#define FNVEXT_WL2NORM      fnvextwl2norm_
#define FNVEXT_L1NORM       fnvextl1norm_
#define FNVEXT_COMPARE      fnvextcompare_
#define FNVEXT_INVTEST      fnvextinvtest_
#define FNVEXT_CONSTRMASK   fnvextconstrmask_
#define FNVEXT_MINQUOTIENT  fnvextminquotient_

#endif

  /* Declarations of global variables */
  extern N_Vector F2C_CVODE_vec;
  extern N_Vector F2C_IDA_vec;
  extern N_Vector F2C_KINSOL_vec;
  extern N_Vector F2C_ARKODE_vec;

  /* Prototype of initialization function */
  void FNVEXT_INIT(int *code, int *ier);
  
  /* Prototypes of the Fortran routines -- these are provided by 
     the user, and define the variables in the supplied vector 
     depending on the physical problem under consideration */
  void FNVEXT_PRINT(void*);
  void FNVEXT_CLONE(void*, void*, int*);
  void FNVEXT_DESTROY(void*);
  void FNVEXT_LINSUM(realtype*, void*, realtype*, void*, void*);
  void FNVEXT_CONST(realtype*, void*);
  void FNVEXT_PROD(void*, void*, void*);
  void FNVEXT_DIV(void*, void*, void*);
  void FNVEXT_SCALE(realtype*, void*, void*);
  void FNVEXT_ABS(void*, void*);
  void FNVEXT_INV(void*, void*);
  void FNVEXT_ADDCONST(realtype*, void*, void*);
  void FNVEXT_DOTPROD(void*, void*, realtype*);
  void FNVEXT_MAXNORM(void*, realtype*);
  void FNVEXT_WRMSNORM(void*, void*, realtype*);
  void FNVEXT_WRMSNORMMASK(void*, void*, void*, realtype*);
  void FNVEXT_MIN(void*, realtype*);
  void FNVEXT_WL2NORM(void*, void*, realtype*);
  void FNVEXT_L1NORM(void*, realtype*);
  void FNVEXT_COMPARE(realtype*, void*, void*);
  void FNVEXT_INVTEST(void*, void*, int*);
  void FNVEXT_CONSTRMASK(void*, void*, void*, int*);
  void FNVEXT_MINQUOTIENT(void*, void*, realtype*);



/****************************************************************
 * PART II: implementation
 ****************************************************************/

/* The FNVEXT implementation of the N_Vector 'content' structure 
   contains the a void** pointer to the actual Fortran vector 
   structure memory location, and a flag indicating whether this 
   vector was allocated from SUNDIALS calls (versus by the user program) */
struct _N_VectorContent_EXT {
  void *data;              /* ptr to Fortran vector structure location */
  booleantype own_data;    /* flag for ownership of data */
};
typedef struct _N_VectorContent_EXT *N_VectorContent_EXT;


/****************************************************************
 * PART III: Macros
 *    NV_CONTENT_EXT, NV_DATA_EXT, NV_OWN_DATA_EXT
 *--------------------------------------------------------------
 * (1) NV_CONTENT_EXT
 *
 *     This macro gives access to the 'content' structure of the
 *     N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_EXT(v) sets
 *     v_cont to be a pointer to the EXT N_Vector
 *     content structure.
 *
 * (2) NV_DATA_EXT
 *
 *     This macro gives access to the actual void** Fortran 
 *     vector object portion held in the 'content' structure.
 *
 *     The assignment v_data = NV_DATA_EXT(v) sets v_data to
 *     be a pointer to the void** data component from v.
 *
 *     The assignment NV_DATA_EXT(v) = v_data sets the 
 *     void** data component in v by storing the pointer v_data.
 *
 * (3) NV_OWN_DATA_EXT
 *
 *     This macro gives access to the flag within the N_Vector 
 *     content structure indicating ownership of the data.
 *
 ****************************************************************/ 

#define NV_CONTENT_EXT(v) ( (N_VectorContent_EXT)(v->content) )
#define NV_DATA_EXT(v) ( NV_CONTENT_EXT(v)->data )
#define NV_OWN_DATA_EXT(v) ( NV_CONTENT_EXT(v)->own_data )



/****************************************************************
 * PART III: Functions exported by nvector_ext
 *
 * CONSTRUCTORS:
 *    N_VNewEmpty_EXT
 *    N_VMake_EXT
 * EXTRA OPERATIONS:
 *    N_VPrint_EXT
 *--------------------------------------------------------------*/

/* N_VNewEmpty_EXT */ 
/* This function creates a new EXT vector with an empty (NULL) content pointer */
N_Vector N_VNewEmpty_EXT();

/* N_VMake_EXT */
/* This function creates an EXT vector with user-provided content pointer */
N_Vector N_VMake_EXT(void* FVec);

/* N_VPrint_EXT */
/* This function prints the content of a EXT vector to stdout */
void N_VPrint_EXT(N_Vector v);



/****************************************************************
 * PART VI: vector operations in nvector_ext
 *--------------------------------------------------------------*/

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_EXT(N_Vector);
SUNDIALS_EXPORT N_Vector    N_VCloneEmpty_EXT(N_Vector);
SUNDIALS_EXPORT N_Vector    N_VClone_EXT(N_Vector);
SUNDIALS_EXPORT void        N_VDestroy_EXT(N_Vector);
SUNDIALS_EXPORT void        N_VSpace_EXT(N_Vector, long int *, long int *);
SUNDIALS_EXPORT realtype*   N_VGetArrayPointer_EXT(N_Vector);
SUNDIALS_EXPORT void        N_VSetArrayPointer_EXT(realtype *, N_Vector);
SUNDIALS_EXPORT void        N_VLinearSum_EXT(realtype, N_Vector, realtype, 
                                             N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VConst_EXT(realtype, N_Vector);
SUNDIALS_EXPORT void        N_VProd_EXT(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VDiv_EXT(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VScale_EXT(realtype, N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VAbs_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VInv_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT void        N_VAddConst_EXT(N_Vector, realtype, N_Vector);
SUNDIALS_EXPORT realtype    N_VDotProd_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VMaxNorm_EXT(N_Vector);
SUNDIALS_EXPORT realtype    N_VWrmsNorm_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VWrmsNormMask_EXT(N_Vector, N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VMin_EXT(N_Vector);
SUNDIALS_EXPORT realtype    N_VWL2Norm_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT realtype    N_VL1Norm_EXT(N_Vector);
SUNDIALS_EXPORT void        N_VCompare_EXT(realtype, N_Vector, N_Vector);
SUNDIALS_EXPORT booleantype N_VInvTest_EXT(N_Vector, N_Vector);
SUNDIALS_EXPORT booleantype N_VConstrMask_EXT(N_Vector, N_Vector, N_Vector);   
SUNDIALS_EXPORT realtype    N_VMinQuotient_EXT(N_Vector, N_Vector);


#ifdef __cplusplus
}
#endif

#endif



/*********************** END OF FILE ***********************/
