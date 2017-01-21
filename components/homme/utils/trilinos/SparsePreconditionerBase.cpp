#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "SparsePreconditionerBase.hpp"

#include "trilinosDefines.h" // Preconditioner Options
#include "trilinosUtils.hpp" // Helpful functions

#include "Teuchos_RCP.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

//#define DEBUG_PRINT_ON

using std::cout;
using std::endl;
using std::flush;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// ==============================================================================
// Prototypes for Fortran functions
// ==============================================================================
extern "C" {

  void Weights_Vel(double*, int, void*);
  void Weights_Ps(double*, int, void*);
  void Weights_T(double*, int, void*);

  // global IDs for building vector and matrix maps
  void GlobalIDs_Vel(int*, int, void*);
  void GlobalIDs_Ps(int*, int, void*);
  void GlobalIDs_T(int*, int, void*);

  // zero out matrix functions
  void Zero_Ablock(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Zero_Bblock(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Zero_Cblock(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Zero_Dblock(int, int, int*, int*, int*, double*, void*);
  void Zero_Eblock(int, int, int*, int*, int*, double*, void*);

  void Zero_Fblock(int, int, int*, int*, int*, double*, void*);
  void Zero_Gblock(int, int, int*, int*, int*, double*, void*);
  void Zero_Hblock(int, int, int*, int*, int*, double*, void*);
  
  // matrix fill functions
  void Jac_Ablock(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Jac_Ablock_VertAdv(int, int, int*, int*, int*, double*, int*, int*, double*, void*);

  void Jac_Bblock_GP(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Jac_Bblock_gradP(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Jac_Bblock_VertAdv(int, int, int*, int*, int*, double*, int*, int*, double*, void*);

  void Jac_Cblock_GP(int, int, int*, int*, int*, double*, int*, int*, double*, void*);
  void Jac_Cblock_gradP(int, int, int*, int*, int*, double*, int*, int*, double*, void*);

  void Jac_Dblock(int, int, int*, int*, int*, double*, void*);

  void Jac_Eblock(int, int, int*, int*, int*, double*, void*);

  void Jac_Fblock_HorizAdv(int, int, int*, int*, int*, double*, void*);
  void Jac_Fblock_VertAdv(int, int, int*, int*, int*, double*, void*);
  void Jac_Fblock_PVertVel(int, int, int*, int*, int*, double*, void*);

  void Jac_Gblock_VertAdv(int, int, int*, int*, int*, double*, void*);
  void Jac_Gblock_PVertVel(int, int, int*, int*, int*, double*, void*);

  void Jac_Hblock_Time(int, int, int*, int*, int*, double*, void*);
  void Jac_Hblock_HorizAdv(int, int, int*, int*, int*, double*, void*);
  void Jac_Hblock_VertAdv(int, int, int*, int*, int*, double*, void*);
  void Jac_Hblock_PVertVel(int, int, int*, int*, int*, double*, void*);
}


SparsePreconditionerBase::SparsePreconditionerBase
( const int np_, const int nlev_, const int nelemd_, void* StateData_,
  RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_) 
  : 
  hommePreconditionerBase(xMap_),
  np(np_), nlev(nlev_), nelemd(nelemd_), StateData(StateData_)
{
  // get process ID, root process does output
  const Epetra_Comm& Comm = xVec_->Comm();
  myid = Comm.MyPID();
  printproc = (myid == 0);

  nState = (3*np*np*nlev + np*np)*(nelemd);
  nVel   = 2*np*np*nlev*nelemd;
  nPs    = np*np*nelemd;
  nT     = np*np*nlev*nelemd;

  // zero timers
  ZeroATime = TimeMonitor::getNewCounter("PC: Zero A");
  ZeroBTime = TimeMonitor::getNewCounter("PC: Zero B");
  ZeroCTime = TimeMonitor::getNewCounter("PC: Zero C");
  ZeroDTime = TimeMonitor::getNewCounter("PC: Zero D");
  ZeroETime = TimeMonitor::getNewCounter("PC: Zero E");
  ZeroFTime = TimeMonitor::getNewCounter("PC: Zero F");
  ZeroGTime = TimeMonitor::getNewCounter("PC: Zero G");
  ZeroHTime = TimeMonitor::getNewCounter("PC: Zero H");

  // fill timers 
  FillATime = TimeMonitor::getNewCounter("PC: Fill A");
  FillBTime = TimeMonitor::getNewCounter("PC: Fill B");
  FillCTime = TimeMonitor::getNewCounter("PC: Fill C");
  FillDTime = TimeMonitor::getNewCounter("PC: Fill D");
  FillETime = TimeMonitor::getNewCounter("PC: Fill E");
  FillFTime = TimeMonitor::getNewCounter("PC: Fill F");
  FillGTime = TimeMonitor::getNewCounter("PC: Fill G");
  FillHTime = TimeMonitor::getNewCounter("PC: Fill H");

  // ============================================================================
  // compute local and global node IDs
  // ============================================================================

  // get weights (1/repeates) for scaling duplicated vector entries 
  VelWgt = new double[nVel];
  PsWgt  = new double[nPs];
  TWgt   = new double[nT];

  Weights_Vel(VelWgt, nVel, StateData);
  Weights_Ps(PsWgt,   nPs,  StateData);
  Weights_T(TWgt,     nT,   StateData);  
  
  // node IDs including duplicates
  VelIdx = new int[nVel];
  PsIdx  = new int[nPs];
  TIdx   = new int[nT];

  GlobalIDs_Vel(VelIdx, nVel, StateData);
  GlobalIDs_Ps(PsIdx,   nPs,  StateData);
  GlobalIDs_T(TIdx,     nT,   StateData);

  // sorted node IDs
  VelIdx_unsorted = new int[nVel];
  PsIdx_unsorted  = new int[nPs];
  TIdx_unsorted   = new int[nT];

  for (int i=0; i<nVel; i++) {
    VelIdx_unsorted[i] = VelIdx[i];
  }

  for (int i=0; i<nPs; i++) {
    PsIdx_unsorted[i] = PsIdx[i];
  }

  for (int i=0; i<nT; i++) {
    TIdx_unsorted[i] = TIdx[i];
  }
  
  std::sort(VelIdx, VelIdx + nVel);
  std::sort(PsIdx,  PsIdx + nPs);
  std::sort(TIdx,   TIdx + nT); 
  
  // number of unique node IDs on this process (not globally unique)
  int* VelEnd = std::unique(VelIdx, VelIdx + nVel);
  int* PsEnd  = std::unique(PsIdx,  PsIdx + nPs);
  int* TEnd   = std::unique(TIdx,   TIdx + nT);   

  nVel_unique = VelEnd - (&VelIdx[0]);
  nPs_unique  = PsEnd  - (&PsIdx[0]);
  nT_unique   = TEnd   - (&TIdx[0]);

  nState_unique = nVel_unique + nPs_unique + nT_unique;

  // count duplicates, allocate map array, and store location of duplicates
  int temp[10]; // <<< need to replace with max number of repeats for regionally refined cases

  // ----------------------------------------------------------------------------
  // velocity index map
  // ----------------------------------------------------------------------------
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Velocity index map" << endl;
#endif

  // allocate arrays to hold the number of times an ID occurs
  nVelIdxMap = new int[nVel_unique];
  
  // initialize duplicate count
  for (int i=0; i < nVel_unique; i++) {
    nVelIdxMap[i] = 0;
  }
  
  VelIdxMap = new int*[nVel_unique];
  
  for (int i=0; i < nVel_unique; i++) {
    for (int j=0; j < nVel; j++) {
      
      // count number of time index occurs and store locations
      if (VelIdx[i] == VelIdx_unsorted[j]) {
	temp[nVelIdxMap[i]] = j; 
	nVelIdxMap[i]++;	
      }
    }
    
    // allocate map for duplicate entries 
    VelIdxMap[i] = new int[nVelIdxMap[i]];
    
    // store locations of duplicate entries
    for (int k=0; k < nVelIdxMap[i]; k++) {
      VelIdxMap[i][k] = temp[k];
    }
  }
  
  // ----------------------------------------------------------------------------
  // surface pressure index map
  // ----------------------------------------------------------------------------
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Surface pressure index map" << endl;
#endif
  
  // allocate and initialize number of times an ID occurs
  nPsIdxMap = new int[nPs_unique];
  
  for (int i=0; i < nPs_unique; i++) {
    nPsIdxMap[i] = 0;
  }
  
  // allocate 
  PsIdxMap = new int*[nPs_unique];
  
  for (int i=0; i < nPs_unique; i++) {
    for (int j=0; j < nPs; j++) {
      
      // count number of time index occurs and store locations
      if (PsIdx[i] == PsIdx_unsorted[j]) {
	temp[nPsIdxMap[i]] = j + nVel; 
	nPsIdxMap[i]++;	
      }
    }    
    
    // allocate map for duplicate entries 
    PsIdxMap[i] = new int[nPsIdxMap[i]];
    
    // store locations of duplicate entries
    for (int k=0; k < nPsIdxMap[i]; k++) {
      PsIdxMap[i][k] = temp[k];
    }
  }
  
  // ----------------------------------------------------------------------------
  // temperature index map
  // ----------------------------------------------------------------------------
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Temperature index map" << endl;
#endif

  // allocate arrays to hold the number of times an ID occurs
  nTIdxMap = new int[nT_unique];

  // initialize duplicate count
  for (int i=0; i < nT_unique; i++) {
    nTIdxMap[i] = 0;
  }
  
  TIdxMap = new int*[nT_unique];
  
  for (int i=0; i < nT_unique; i++) {
    for (int j=0; j < nT; j++) {
      
      // count number of time index occurs and store locations
      if (TIdx[i] == TIdx_unsorted[j]) {
	temp[nTIdxMap[i]] = j + nVel + nPs; 
	nTIdxMap[i]++;	
      }
    }
    
    // allocate map for duplicate entries 
    TIdxMap[i] = new int[nTIdxMap[i]];
    
    // store locations of duplicate entries
    for (int k=0; k < nTIdxMap[i]; k++) {
      TIdxMap[i][k] = temp[k];
    }
  }
  
  // number of non-unique nodes on this process
  nRows_3D = np*np*nlev;
  nRows_2D = np*np;
  
  // approximate number of non-zeros per row (too small?) 
  nnz_2D = np+np-1;       // 2d var
  nnz_3D = nnz_2D * nlev; // 3d var
  
  // maximum number of nonzeros on this process
  nnz_max = 2*nnz_3D*nRows_3D;
  
  // actual number of nonzeros per row
  row_nnz_varies = new int[nRows_3D];
  
  // arrays for filling matrices
  Urows = new int[nRows_3D];
  Ucols = new int[nnz_max];
  Uvals = new double[nnz_max]; 
  
  Vrows = new int[nRows_3D];
  Vcols = new int[nnz_max];
  Vvals = new double[nnz_max]; 
}


// ========================================================================
// A block 
// ========================================================================
int SparsePreconditionerBase::ZeroAblock(RCP<Epetra_FECrsMatrix> Ablock)
{    
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing A Block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*ZeroATime);

  // initialize error flag
  int ierr = 0;
  
  if (Ablock->Filled()) {
    
    Ablock->PutScalar(0.0);
    
  } else {
          
    for (int ie=1; ie <= nelemd; ie++) {

      Zero_Ablock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, 
                  Vrows, Vcols, Vvals, StateData);
      
      // U equations
      for (int i=0; i < nRows_3D; i++) {
        ierr = Ablock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);      
        check_flag(&ierr, "Ablock U, InsertGlobalValues", 1, myid);
      }
      
      // V equations
      for (int i=0; i < nRows_3D; i++) {
        ierr = Ablock->InsertGlobalValues(Vrows[i], row_nnz_fixed, 
                                          Vvals + i*row_nnz_fixed, 
                                          Vcols + i*row_nnz_fixed);	
        check_flag(&ierr, "Ablock V, InsertGlobalValues", 1, myid);
      }
    }
  }  
  return ierr;
}


int SparsePreconditionerBase::FillAblock(RCP<Epetra_FECrsMatrix> Ablock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling A Block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*FillATime);

  // initialize error flag
  int ierr = 0;
  
  for (int ie=1; ie <= nelemd; ie++) {
    
    // time derivative, vorticity, coriolis, kinetic energy   
    Jac_Ablock(ie, nnz_max, &row_nnz_fixed, 
               Urows, Ucols, Uvals, 
               Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Ablock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Ablock U, SumIntoGlobalValues", 1, myid);
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Ablock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);	
      check_flag(&ierr, "Ablock V, SumIntoGlobalValues", 1, myid);
    }
    
    // vertical advection   
    Jac_Ablock_VertAdv(ie, nnz_max, &row_nnz_fixed, 
                       Urows, Ucols, Uvals, 
                       Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {
      
      ierr = Ablock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Ablock U Vert Adv, SumIntoGlobalValues", 1, myid);	
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Ablock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);
      check_flag(&ierr, "Ablock V Vert Adv, SumIntoGlobalValues", 1, myid);
    }
  } // element loop
  return ierr;
}


// ==========================================================================
// B Block
// ==========================================================================
int SparsePreconditionerBase::ZeroBblock(RCP<Epetra_FECrsMatrix> Bblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing Out B Matrix" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*ZeroBTime);

  // initialize error flag
  int ierr = 0;
  
  // zero out B matrix
  if (Bblock->Filled()) {
    
    Bblock->PutScalar(0.0);
    
  } else {

    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Bblock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, 
                  Vrows, Vcols, Vvals, StateData);
      
      // u equations 
      for (int i=0; i < nRows_3D; i++) {
        ierr = Bblock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);
        check_flag(&ierr, "Bblock U, InsertGlobalValues", 1, myid);
      }
      
      // v equations
      for (int i=0; i < nRows_3D; i++) {
        ierr = Bblock->InsertGlobalValues(Vrows[i], row_nnz_fixed, 
                                          Vvals + i*row_nnz_fixed, 
                                          Vcols + i*row_nnz_fixed); 	
        check_flag(&ierr, "Bblock V, InsertGlobalValues", 1, myid);
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillBblock(RCP<Epetra_FECrsMatrix> Bblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling B Matrix" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*FillBTime);
  
  // initialize error flag
  int ierr = 0;
  
  // fill B block
  for (int ie=1; ie <= nelemd; ie++) {
    
    // geopotential
    Jac_Bblock_GP(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, Vrows, Vcols, Vvals, StateData);
	
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock U geopotential, SumIntoGlobalValues", 1, myid);
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock V geopotential, SumIntoGlobalValues", 1, myid);
    }
    
    // vertical advection
    Jac_Bblock_VertAdv(ie, nnz_max, &row_nnz_fixed, 
                       Urows, Ucols, Uvals, Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock U Vert Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock V Vert Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // pressure gradient
    Jac_Bblock_gradP(ie, nnz_max, &row_nnz_fixed, 
                     Urows, Ucols, Uvals, Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock U gradp, SumIntoGlobalValues", 1, myid);
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Bblock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);
      check_flag(&ierr, "Bblock V gradp, SumIntoGlobalValues", 1, myid);
    } 
  } // element loop
  return ierr;
}


// ==========================================================================
// C Block
// ==========================================================================
int SparsePreconditionerBase::ZeroCblock(RCP<Epetra_FECrsMatrix> Cblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing out C block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*ZeroCTime);

  // initialize error flag
  int ierr = 0;
  
  // zero out C block  
  if (Cblock->Filled()) {
    
    Cblock->PutScalar(0.0);
    
  } else {

    for (int ie=1; ie <= nelemd; ie++) {
    
      Zero_Cblock(ie, nnz_max, row_nnz_varies, 
                  Urows, Ucols, Uvals, 
                  Vrows, Vcols, Vvals, StateData);
      
      // U equations
      for (int i=0, offset=0; i < nRows_3D; i++) {
        ierr = Cblock->InsertGlobalValues(Urows[i], row_nnz_varies[i], 
                                          Uvals + offset, Ucols + offset);
        check_flag(&ierr, "Cblock U, InsertGlobalValues", 1, myid);
        
        offset += row_nnz_varies[i];
      }
      
      // V equations
      for (int i=0, offset=0; i < nRows_3D; i++) {
        ierr = Cblock->InsertGlobalValues(Vrows[i], row_nnz_varies[i], 
                                          Vvals + offset, Vcols + offset); 	
        check_flag(&ierr, "Cblock V, InsertGlobalValues", 1, myid);
        
        offset += row_nnz_varies[i];
      }
    }
  }
  return ierr; 
}


int SparsePreconditionerBase::FillCblock(RCP<Epetra_FECrsMatrix> Cblock)
{ 
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Fill C Matrix" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*FillCTime);

  // initialize error flag
  int ierr = 0;
  
  // fill C block
  for (int ie=1; ie <= nelemd; ie++) {
    
    // geopotential
    Jac_Cblock_GP(ie, nnz_max, row_nnz_varies, 
                  Urows, Ucols, Uvals, Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0,offset=0; i < nRows_3D; i++) {

	ierr = Cblock->SumIntoGlobalValues(Urows[i], row_nnz_varies[i], 
                                           Uvals + offset, 
                                           Ucols + offset);
	check_flag(&ierr, "Cblock U geopotential, SumIntoGlobalValues", 1, myid);
	
	offset += row_nnz_varies[i];
    }
    
    // sum in entries for V equations
    for (int i=0,offset=0; i < nRows_3D; i++) {

      ierr = Cblock->SumIntoGlobalValues(Vrows[i], row_nnz_varies[i], 
                                         Vvals + offset, 
                                         Vcols + offset);
      check_flag(&ierr, "Cblock V geopotential, SumIntoGlobalValues", 1, myid);
      
      offset += row_nnz_varies[i];
    }
    
    // pressure gradient
    Jac_Cblock_gradP(ie, nnz_max, &row_nnz_fixed, 
                     Urows, Ucols, Uvals, Vrows, Vcols, Vvals, StateData);
    
    // sum in entries for U equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Cblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Cblock U gradp, SumIntoGlobalValues", 1, myid);
    }
    
    // sum in entries for V equations
    for (int i=0; i < nRows_3D; i++) {

      ierr = Cblock->SumIntoGlobalValues(Vrows[i], row_nnz_fixed, 
                                         Vvals + i*row_nnz_fixed, 
                                         Vcols + i*row_nnz_fixed);
      check_flag(&ierr, "Cblock V gradp, SumIntoGlobalValues", 1, myid);
    }
  } 
  return ierr;
}


// ==========================================================================
// D Block
// ==========================================================================
int SparsePreconditionerBase::ZeroDblock(RCP<Epetra_FECrsMatrix> Dblock)
{ 
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing Out D Matrix" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*ZeroDTime);

  // initialize error flag
  int ierr = 0;
  
  // zero out D block  
  if (Dblock->Filled()) {
    
    Dblock->PutScalar(0.0);
    
  } else {
    
    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Dblock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, StateData);
      
      for (int i=0; i < nRows_2D; i++) {
        ierr = Dblock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);
        check_flag(&ierr, "Dblock, InsertGlobalValues", 1, myid);
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillDblock(RCP<Epetra_FECrsMatrix> Dblock)
{ 
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling D Block" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*FillDTime);
  
  // initialize error flag
  int ierr = 0;

  // fill D block 
  for (int ie=1; ie <= nelemd; ie++) {
    
    Jac_Dblock(ie, nnz_max, &row_nnz_fixed, 
               Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_2D; i++) {

      ierr = Dblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Dblock, SumIntoGlobalValues", 1, myid);
    }
  }
  return ierr; 
}


// ==========================================================================
// E Block
// ==========================================================================
int SparsePreconditionerBase::ZeroEblock(RCP<Epetra_FECrsMatrix> Eblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing Out E Block" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*ZeroETime);

  // initialize error flag
  int ierr = 0;
  
  // zero out E block
  if (Eblock->Filled()) {
    
    Eblock->PutScalar(0.0);

  } else {

    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Eblock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, StateData);
      
      for (int i=0; i < nRows_2D; i++) {
        ierr = Eblock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);      
        check_flag(&ierr, "Eblock, InsertGlobalValues", 1, myid);
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillEblock(RCP<Epetra_FECrsMatrix> Eblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling E Block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*FillETime);
  
  // initialize error flag
  int ierr = 0;
  
  // fill E block 
  for (int ie=1; ie <= nelemd; ie++) {
    
    Jac_Eblock(ie, nnz_max, &row_nnz_fixed, 
               Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_2D; i++) {

      ierr = Eblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Eblock, SumIntoGlobalValues", 1, myid);
    }
  }
  return ierr;
}


// ==========================================================================
// F Block
// ==========================================================================
int SparsePreconditionerBase::ZeroFblock(RCP<Epetra_FECrsMatrix> Fblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing Out F Block" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*ZeroFTime);

  // initialize error flag
  int ierr = 0;
  
  if (Fblock->Filled()) {

    Fblock->PutScalar(0.0);

  } else {
         
    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Fblock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, StateData);
      
      for (int i=0; i < nRows_3D; i++) {
        ierr = Fblock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);
        check_flag(&ierr, "Fblock, InsertGlobalValues", 1, myid);
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillFblock(RCP<Epetra_FECrsMatrix> Fblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling F Block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*FillFTime);
  
  // initialize error flag
  int ierr = 0;
             
  // fill F block     
  for (int ie=1; ie <= nelemd; ie++) {
    
    // horizontal advection
    Jac_Fblock_HorizAdv(ie, nnz_max, &row_nnz_fixed, 
                        Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Fblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Fblock Horiz Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // vertical advection
    Jac_Fblock_VertAdv(ie, nnz_max, &row_nnz_fixed, 
                       Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Fblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Fblock Vert Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // pressure vertical velocity
    Jac_Fblock_PVertVel(ie, nnz_max, row_nnz_varies, 
                        Urows, Ucols, Uvals, StateData);
    
    for (int i=0,offset=0; i < nRows_3D; i++) {

      ierr = Fblock->SumIntoGlobalValues(Urows[i], row_nnz_varies[i], 
                                         Uvals + offset, 
                                         Ucols + offset);
      check_flag(&ierr, "Fblock P Vert Vel, SumIntoGlobalValues", 1, myid);
      
      offset += row_nnz_varies[i];
    }
  }
  return ierr; 
}


// ==========================================================================
// G Block
// ==========================================================================
int SparsePreconditionerBase::ZeroGblock(RCP<Epetra_FECrsMatrix> Gblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing G block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*ZeroGTime);

  // initialize error flag
  int ierr = 0;
  
  // zero out G block
  if (Gblock->Filled()) {

    Gblock->PutScalar(0.0);

  } else {
         
    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Gblock(ie, nnz_max, &row_nnz_fixed, 
                  Urows, Ucols, Uvals, StateData);
      
      for (int i=0; i < nRows_3D; i++) {
        ierr = Gblock->InsertGlobalValues(Urows[i], row_nnz_fixed, 
                                          Uvals + i*row_nnz_fixed, 
                                          Ucols + i*row_nnz_fixed);
        check_flag(&ierr, "Gblock, InsertGlobalValues", 1, myid);
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillGblock(RCP<Epetra_FECrsMatrix> Gblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling G Block" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*FillGTime);
  
  // initialize error flag
  int ierr = 0;
  
  // fill G block 
  for (int ie=1; ie <= nelemd; ie++) {
    
    // vertical advection
    Jac_Gblock_VertAdv(ie, nnz_max, &row_nnz_fixed, 
                       Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Gblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Gblock Vert Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // pressure vertical velocity
    Jac_Gblock_PVertVel(ie, nnz_max, &row_nnz_fixed, 
                        Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Gblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Gblock P Vert Vel, SumIntoGlobalValues", 1, myid);
    }
  }
  return ierr;
}


// ==========================================================================
// H Block
// ==========================================================================
int SparsePreconditionerBase::ZeroHblock(RCP<Epetra_FECrsMatrix> Hblock)
{  
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Zeroing H block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*ZeroHTime);

  // initialize error flag
  int ierr = 0;
  
  if (Hblock->Filled()) {
    
    Hblock->PutScalar(0.0);

  } else {
    
    for (int ie=1; ie <= nelemd; ie++) {
      
      Zero_Hblock(ie, nnz_max, row_nnz_varies, 
                  Urows, Ucols, Uvals, StateData);
      
      for (int i=0,offset=0; i < nRows_3D; i++) {
        ierr = Hblock->InsertGlobalValues(Urows[i], row_nnz_varies[i], 
                                          Uvals + offset, 
                                          Ucols + offset);
        check_flag(&ierr, "Hblock, InsertGlobalValues", 1, myid);
        
        offset += row_nnz_varies[i];
      }
    }
  }
  return ierr;
}


int SparsePreconditionerBase::FillHblock(RCP<Epetra_FECrsMatrix> Hblock)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Filling H block" << endl;
#endif

  // start timer
  TimeMonitor LocalTimer(*FillHTime);

  // initialize error flag
  int ierr = 0;
    
  // fill H block
  for (int ie=1; ie <= nelemd; ie++) {
    
    // time derivative
    Jac_Hblock_Time(ie, nnz_max, &row_nnz_fixed, 
                    Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Hblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Hblock time, SumIntoGlobalValues", 1, myid);
    }
    
    // horizontal advection
    Jac_Hblock_HorizAdv(ie, nnz_max, &row_nnz_fixed, 
                        Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Hblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Hblock Horiz Adv, SumIntoGlobalValues", 1, myid);
    }
    
    // vertical advection
    Jac_Hblock_VertAdv(ie, nnz_max, row_nnz_varies, 
                       Urows, Ucols, Uvals, StateData);
    
    for (int i=0,offset=0; i < nRows_3D; i++) {

      ierr = Hblock->SumIntoGlobalValues(Urows[i], row_nnz_varies[i], 
                                         Uvals + offset, 
                                         Ucols + offset);
      check_flag(&ierr, "Hblock Vert Adv, SumIntoGlobalValues", 1, myid);
      
      offset += row_nnz_varies[i];
    }
    
    // pressure vertical velocity
    Jac_Hblock_PVertVel(ie, nnz_max, &row_nnz_fixed, 
                        Urows, Ucols, Uvals, StateData);
    
    for (int i=0; i < nRows_3D; i++) {

      ierr = Hblock->SumIntoGlobalValues(Urows[i], row_nnz_fixed, 
                                         Uvals + i*row_nnz_fixed, 
                                         Ucols + i*row_nnz_fixed);
      check_flag(&ierr, "Hblock P Vert Vel, SumIntoGlobalValues", 1, myid);
    }
  } // element loop
  return ierr; 
}



// ========================================================================
// Build RHS
// ========================================================================
// int SparsePreconditionerBase::BuildRHS(RCP<Epetra_FECrsMatrix> Ablock)
// {
//   TimeMonitor LocalTimer(*BuildRHSTime);
  
//   StateVec_local->PutScalar(0.0);
  
//   // Velocities
//   double VelVals[nVel];
  
//   for (int i=0; i < nVel; i++)
//     VelVals[i] = VelWgt[i] * X[0][i]; 
  
//   ierr = StateVec_local->SumIntoGlobalValues(nVel, VelVals, VelIdx_unsorted);
//   check_flag(&ierr,"SumIntoGlobalValues: VelVec_local Velocities",1);
  
//   // Surface pressure
//   double PsVals[nPs];
  
//   for (int i = 0; i < nPs; i++)
//     PsVals[i] = PsWgt[i] * X[0][i+nVel];
  
//   ierr = StateVec_local->SumIntoGlobalValues(nPs, PsVals, PsIdx_unsorted);
//   check_flag(&ierr,"SumIntoGlobalValues: StateVec_local Surface Pressure",1);
  
//   // Temperatures
//   double TVals[nT];
  
//   for (int i = 0; i < nT; i++)
//     TVals[i] = TWgt[i] * X[0][i+nVel+nPs]; 
  
//   ierr = StateVec_local->SumIntoGlobalValues(nT, TVals, TIdx_unsorted);
//   check_flag(&ierr,"SumIntoGlobalValues: StateVec_local Temperature",1);
  
//   // Export to local to global
//   StateVec->PutScalar(0.0);
  
//   ierr = StateVec->Export(*StateVec_local, *State_Exporter, Add);
//   check_flag(&ierr,"Export: StateVec",1);

//   return ierr;
// }



// ==============================================================================
// Ouput Matrix to file 
// ==============================================================================

void SparsePreconditionerBase::print_matrix(RCP<Epetra_FECrsMatrix> Matrix, 
                                            const std::string &MatrixName)
{    
  std::string MyPID;
  std::ostringstream temp;
  std::string fname;
  std::ofstream outFile;
  
  temp << std::setw(2) << std::setfill('0') << myid + 1;
  MyPID = temp.str(); 
  
  // write global matrix to file
  fname = MatrixName+"block_proc"+MyPID+".txt";
  outFile.open(fname.c_str());
  outFile.precision(16);
  
  cout << "Writing Global matrix to disk:" << MatrixName << endl;
  Matrix->Print(outFile);
  outFile.close();
}
