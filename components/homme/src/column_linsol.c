/*******************************************************************
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2017; all rights reserved
 * -----------------------------------------------------------------
 * This is the implementation file for a Fortran columnwise linear 
 * solver package for HOMME+ARKode.  This is specifically created
 * to pair with the NVECTOR_EXTERNAL vector implementation.  It is 
 * assumed that the linear solver 'setup' is empty, and that the 
 * 'solve' routine internally sets up the linear system matrix 
 * on-the-fly and uses it in the solve process.
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "column_linsol.h"
#include "nvector_external.h"


/********************* Fortran Interface Functions ***************/

/* define global matrix and linear solver variables for Fortran interface */
extern void *ARK_arkodemem;


/* FCOLUMNSOL_INIT is the Fortran interface routine to initialize the 
   Fortran linear solver object. */
void FCOLUMNSOL_INIT(int *ier)
{
  ARKodeMem ark_mem;

  /* initialize return value to success */
  *ier = 0;

  /* verify that ARKode FCMIX interface has been initialized */
  if (ARK_arkodemem == NULL) {
    *ier = -1;
    return;
  }
  ark_mem = (ARKodeMem) ARK_arkodemem;

  /* free any existing system solver attached to ARKode */
  if (ark_mem->ark_lfree)  ark_mem->ark_lfree(ark_mem);

  /* Set system linear solver components of the ARKode memory structure -- 
     only lsolve is required (and that is all that is used here).
     Also specify that this is a 'custom' linear solver */
  ark_mem->ark_lsolve = ColumnSolSolve;
  ark_mem->ark_linit  = NULL;
  ark_mem->ark_lsetup = NULL;
  ark_mem->ark_lfree  = NULL;
  ark_mem->ark_lsolve_type = 2;
  ark_mem->ark_lmem = NULL;
}


/* ColumnSolSolve is the C interface routine to call the Fortran-supplied 
   FCOLUMNSOL_SOLVE routine.
 */
int ColumnSolSolve(struct ARKodeMemRec *ark_mem, N_Vector b,
                   N_Vector ycur, N_Vector fcur) {

  /* declare local data */
  void *bd, *yd;
  double tcur, gamma;
  int ier;

  /* initialize return flag to success */
  ier = 0;

  /* access current time */
  tcur = ark_mem->ark_tn;

  /* extract N_Vector data pointers */
  bd = yd = NULL;
  bd  = NV_DATA_EXT(b);
  yd  = NV_DATA_EXT(ycur);
  gamma = ark_mem->ark_gamma;

  /* call Fortran routine to do operation */
  FCOLUMNSOL_SOLVE(bd, &tcur, yd, &gamma, &ier);
  return(ier);
}


/*********************** END OF FILE ***********************/
