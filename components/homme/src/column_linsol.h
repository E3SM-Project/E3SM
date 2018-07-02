/*******************************************************************
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2017; all rights reserved
 * -----------------------------------------------------------------
 * This is the header file for a Fortran columnwise linear solver 
 * package for HOMME+ARKode.  This is specifically created
 * to pair with the NVECTOR_EXTERNAL vector implementation.  It is 
 * assumed that the linear solver 'setup' is empty, and that the 
 * 'solve' routine internally sets up the linear system matrix 
 * on-the-fly and uses it in the solve process.  
 *******************************************************************/

#ifndef _COLUMNSOL_H
#define _COLUMNSOL_H

#include <sundials/sundials_nvector.h>
#include <arkode/arkode_impl.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*-----------------------------------------------------------------
  Part I: Fortran callable wrappers for the interface routines, as
  well as prototypes of the supplied Fortran vector operations
  to be called by C routines.
  -----------------------------------------------------------------*/

/* Define Fortran/C interface wrappers */
#if defined(SUNDIALS_F77_FUNC)

#define FCOLUMNSOL_INIT   SUNDIALS_F77_FUNC(fcolumnsolinit,  FCOLUMNSOLINIT)
#define FCOLUMNSOL_SOLVE  SUNDIALS_F77_FUNC(fcolumnsolsolve, FCOLUMNSOLSOLVE)

#else

#define FCOLUMNSOL_INIT   fcolumnsolinit_
#define FCOLUMNSOL_SOLVE  fcolumnsolsolve_

#endif

  /* Prototype of initialization function -- called from Fortran to initialize 
     this linear solver interface, and attach it to the ARKode memory structure */
  void FCOLUMNSOL_INIT(int *ier);

  /* Prototype of the Fortran-supplied 'solve' routine: given a RHS vector B, the 
     current time T, the current solution, Y, and the current time step scale factor, 
     GAMMA, this solves the linear systems (defined internally), storing the 
     solution back in the vector B.  

     The return value should be 0 on success, a positive value on a recoverable 
     failure (e.g. an illegal value that might be fixed with a smaller time step 
     size), and a negative value on a non-recoverable failure.
  */
  void FCOLUMNSOL_SOLVE(void* B, double* T, void* Y, double *GAMMA, int* IER);

  /* Prototype of the C interface to FCOLUMNSOL_SOLVE -- this is the 'lsolve' 
     routine attached to ARKode's memory structure */
  int ColumnSolSolve(struct ARKodeMemRec *ark_mem, N_Vector b,
                     N_Vector ycur, N_Vector fcur);


#ifdef __cplusplus
}
#endif

#endif

