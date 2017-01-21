#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosDefines.h" // Preconditioner Options
#include "trilinosUtils.hpp" // Helpful functions

#include "BlockPreconditioner_Sparse.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

// extra output for script to read residual values
//#define PRINT_SOLVE_INFO 

//#define PREC_IDENT        // Identity preconditioner for testing
#define PREC_SIMPLE_2D      // 2D version of the SIMPLE preconditioner
//#define PREC_SIMPLEC_2D   // 2D version of the SIMPLE C preconditioner

//#define DEBUG_PRINT_ON

using std::cout;
using std::endl;
using std::flush;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// ==============================================================================
// Constructor for Block Preconditioner
// ==============================================================================
BlockPreconditioner_Sparse::BlockPreconditioner_Sparse
( const int np_, const int nlev_, const int nelemd_,
  RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, 
  void* StateData_, 
  void (*precUpdateFunction_)(double*, int, void*),
  const RCP<ParameterList>& ASolvePL_,
  const RCP<ParameterList>& SSolvePL_,
  const RCP<ParameterList>& PSolvePL_,
  int* ATotalIters_,
  int* STotalIters_,
  int* PTotalIters_)
  : 
  SparsePreconditionerBase(np_, nlev_, nelemd_, StateData_, xVec_, xMap_), // Required Base Class construction 
  precUpdateFunction(precUpdateFunction_),
  ASolvePL(ASolvePL_),
  SSolvePL(SSolvePL_),
  PSolvePL(PSolvePL_),
  ATotalIters(ATotalIters_),
  STotalIters(STotalIters_),
  PTotalIters(PTotalIters_),
  firstfill(true)
{

  if (printproc) 
    cout << " >>> Constructing BlockPreconditioner_Sparse <<< " << endl;

  // MPI communicator
  const Epetra_Comm& Comm = xVec_->Comm();

  // ============================================================================
  // create timers
  // ============================================================================
  ComputePCTime = TimeMonitor::getNewCounter("PC: Compute PC");
  UpdateStateTime = TimeMonitor::getNewCounter("PC: Update State");

  ABlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build A");
  BBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build B");
  CBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build C");
  DBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build D");
  EBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build E");
  FBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build F");
  GBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build G");
  HBlockAssemblyTime = TimeMonitor::getNewCounter("PC: Build H");

  BuildAhatInv   = TimeMonitor::getNewCounter("PC: Build AhatInv"); 
  BuildAhatInvB	 = TimeMonitor::getNewCounter("PC: Build AhatInv B");		
  BuildAhatInvC	 = TimeMonitor::getNewCounter("PC: Build AhatInv C");		
  BuildDAhatInvB = TimeMonitor::getNewCounter("PC: Build D AhatInv B"); 
  BuildDAhatInvC = TimeMonitor::getNewCounter("PC: Build D AhatInv C"); 
  BuildFAhatInvB = TimeMonitor::getNewCounter("PC: Build F AhatInv B"); 
  BuildFAhatInvC = TimeMonitor::getNewCounter("PC: Build F AhatInv C"); 
  BuildSchurS	 = TimeMonitor::getNewCounter("PC: Build Schur S"); 
  BuildShatInv	 = TimeMonitor::getNewCounter("PC: Build ShatInv"); 
  BuildX	 = TimeMonitor::getNewCounter("PC: Build X"); 
  BuildR	 = TimeMonitor::getNewCounter("PC: Build R"); 
  BuildRX	 = TimeMonitor::getNewCounter("PC: Build RX"); 
  BuildSchurP	 = TimeMonitor::getNewCounter("PC: Build Schur P"); 
  
  ApplyTime = TimeMonitor::getNewCounter("PC: Apply PC");

  BuildRHSTime = TimeMonitor::getNewCounter("PC: Build RHS Vec");

  ASolveTime = TimeMonitor::getNewCounter("PC: Solve A");
  SSolveTime = TimeMonitor::getNewCounter("PC: Solve S");
  PSolveTime = TimeMonitor::getNewCounter("PC: Solve P");

  BuildYvecTime = TimeMonitor::getNewCounter("PC: Build Y Vec");

  // ============================================================================
  // create Epetra objects
  // ============================================================================

  // ----------------------------------------------------------------------------
  // local and global Epetra Maps
  // ----------------------------------------------------------------------------
  
  // create local maps 
  VelMap_local = rcp(new Epetra_Map(-1, nVel_unique, VelIdx, 0, Comm)); 
  PsMap_local  = rcp(new Epetra_Map(-1, nPs_unique,  PsIdx,  0, Comm));
  TMap_local   = rcp(new Epetra_Map(-1, nT_unique,   TIdx,   0, Comm));

  // create global maps
  VelMap = rcp(new Epetra_Map(Epetra_Util::Create_OneToOne_Map(*VelMap_local)));
  PsMap  = rcp(new Epetra_Map(Epetra_Util::Create_OneToOne_Map(*PsMap_local)));
  TMap   = rcp(new Epetra_Map(Epetra_Util::Create_OneToOne_Map(*TMap_local)));

  // ----------------------------------------------------------------------------
  // Epetra export and import objects
  // ----------------------------------------------------------------------------

  // exporters (source map = local, target map = global)
  Vel_Exporter = rcp(new Epetra_Export(*VelMap_local, *VelMap));
  Ps_Exporter  = rcp(new Epetra_Export(*PsMap_local,  *PsMap));
  T_Exporter   = rcp(new Epetra_Export(*TMap_local,   *TMap));

  // importers (target map = local, source map = global)
  Vel_Importer = rcp(new Epetra_Import(*VelMap_local, *VelMap));
  Ps_Importer  = rcp(new Epetra_Import(*PsMap_local,  *PsMap));
  T_Importer   = rcp(new Epetra_Import(*TMap_local,   *TMap));

  // ----------------------------------------------------------------------------
  // local and global Epetra solution vectors
  // ----------------------------------------------------------------------------
  bool zeroout=true;

  VelVec_local = rcp(new Epetra_Vector(*VelMap_local, zeroout));
  PsVec_local  = rcp(new Epetra_Vector(*PsMap_local,  zeroout));
  TVec_local   = rcp(new Epetra_Vector(*TMap_local,   zeroout));

  VelVec = rcp(new Epetra_Vector(*VelMap, zeroout));
  PsVec  = rcp(new Epetra_Vector(*PsMap,  zeroout));
  TVec   = rcp(new Epetra_Vector(*TMap,   zeroout));

  // ----------------------------------------------------------------------------
  // preconditioner matrices and vectors
  // ----------------------------------------------------------------------------

  // use vector to store hat{A}^{inv} since it is a diagonal matrix
  Ahatinv = rcp(new Epetra_Vector(*VelMap, zeroout)); 

  // solution vectors in preconditioner linear solves
  Ax = rcp(new Epetra_Vector(*VelMap, zeroout));
  Sx = rcp(new Epetra_Vector(*PsMap,  zeroout));
  Px = rcp(new Epetra_Vector(*TMap,   zeroout)); 

  // temporary vectors used in preconditioner
  DAx = rcp(new Epetra_Vector(*PsMap, zeroout));
  w1  = rcp(new Epetra_Vector(*VelMap, zeroout));

  // ----------------------------------------------------------------------------
  // Jacobian blocks
  // ----------------------------------------------------------------------------
  Ablock = rcp(new Epetra_FECrsMatrix(Copy, *VelMap, 8*nnz_3D));
  Bblock = rcp(new Epetra_FECrsMatrix(Copy, *VelMap, 4*nnz_2D));
  Cblock = rcp(new Epetra_FECrsMatrix(Copy, *VelMap, 4*nnz_3D));

  Dblock = rcp(new Epetra_FECrsMatrix(Copy, *PsMap,  8*nnz_3D));
  Eblock = rcp(new Epetra_FECrsMatrix(Copy, *PsMap,  4*nnz_2D));

  Fblock = rcp(new Epetra_FECrsMatrix(Copy, *TMap,   8*nnz_3D));
  Gblock = rcp(new Epetra_FECrsMatrix(Copy, *TMap,   4*nnz_2D));
  Hblock = rcp(new Epetra_FECrsMatrix(Copy, *TMap,   nnz_3D));

  // AhatinvB = rcp(new Epetra_FECrsMatrix(*Bblock));
  Wmatrix = rcp(new Epetra_FECrsMatrix(Copy, *PsMap, 2*nnz_3D));
  Smatrix = rcp(new Epetra_CrsMatrix(Copy, *PsMap, 2*nnz_3D));
  
  // ============================================================================
  // create Belos linear problems (Op, LHS, RHS) and solver managers
  // ============================================================================  
  AProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Ablock, Ax, VelVec));
  ASolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(AProblem, ASolvePL));

  SProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Smatrix, Sx, PsVec));
  SSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(SProblem, SSolvePL));
  
  PProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Hblock, Px, TVec));
  PSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(PProblem, PSolvePL));
}



// ==============================================================================
// Update preconditioner
// ==============================================================================
int BlockPreconditioner_Sparse::computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "Computing Preconditioner" << endl;
#endif
  
  // start timer
  TimeMonitor LocalTimer(*ComputePCTime);

  // error flag
  int ierr = 0;
  
  try{

    {
      TimeMonitor LocalTimer(*UpdateStateTime);
      
      // update pointer to HOMME data
      StateData = StateData_;
      
      // update HOMME data with current iterate values
      precUpdateFunction(xVec->Values(), nState, StateData); 
    }
   
    // build Jacobian blocks
    ZeroAblock(Ablock);
    FillAblock(Ablock);
    {
      TimeMonitor LocalTimer(*ABlockAssemblyTime);
      
      ierr = Ablock->GlobalAssemble(firstfill, Add, true);  
      check_flag(&ierr, "GlobalAssemble: A Block", 1);
    }

    ZeroBblock(Bblock);
    FillBblock(Bblock);
    {
      TimeMonitor LocalTimer(*BBlockAssemblyTime);

      ierr = Bblock->GlobalAssemble(*PsMap, *VelMap, firstfill, Add, true);  
      check_flag(&ierr, "Bblock, Global FillComplete", 1);     
    }

    // ZeroCblock(Cblock);
    // FillCblock(Cblock);
    // {
    //   TimeMonitor LocalTimer(*CBlockAssemblyTime);

    //   ierr = Cblock->GlobalAssemble(*TMap, *VelMap, firstfill, Add, true);
    //   check_flag(&ierr, "Cblock, Global FillComplete", 1);
    // }
 
    ZeroDblock(Dblock);
    FillDblock(Dblock);
    {
      TimeMonitor LocalTimer(*DBlockAssemblyTime);

      ierr = Dblock->GlobalAssemble(*VelMap, *PsMap, firstfill, Add, true);
      check_flag(&ierr, "Dblock, Global FillComplete", 1);
    }

    ZeroEblock(Eblock);
    FillEblock(Eblock);
    {
      TimeMonitor LocalTimer(*EBlockAssemblyTime);

      ierr = Eblock->GlobalAssemble(firstfill, Add, true);
      check_flag(&ierr, "Eblock, Global FillComplete", 1);
    }

    // ZeroFblock(Fblock);
    // FillFblock(Fblock);
    // {
    //   TimeMonitor LocalTimer(*FBlockAssemblyTime);

    //   ierr = Fblock->GlobalAssemble(*VelMap, *TMap, firstfill, Add, true);  
    //   check_flag(&ierr, "Fblock, Global FillComplete", 1);
    // }

    // ZeroGblock(Gblock);
    // FillGblock(Gblock);
    // {
    //   TimeMonitor LocalTimer(*GBlockAssemblyTime);

    //   ierr = Gblock->GlobalAssemble(*PsMap, *TMap, firstfill, Add, true);
    //   check_flag(&ierr, "Gblock, Global FillComplete", 1);
    // }
    
    ZeroHblock(Hblock);
    FillHblock(Hblock);
    {
      TimeMonitor LocalTimer(*HBlockAssemblyTime);

      ierr = Hblock->GlobalAssemble(firstfill, Add, true);
      check_flag(&ierr, "Hblock, Global FillComplete", 1);
    }

    if (firstfill) firstfill = false;

#if defined(CHECK_JAC)
    print_matrix(Ablock,"A");
    print_matrix(Bblock,"B");
    // print_matrix(Cblock,"C");
    print_matrix(Dblock,"D");
    print_matrix(Eblock,"E");
    // print_matrix(Fblock,"F");
    // print_matrix(Gblock,"G");
    print_matrix(Hblock,"H");
    cout << myid << " Check Jacobian Complete" << endl;
    haltrun();
#endif
    
    // approximate A block inverse hat{A}^{-1}
    {
      TimeMonitor LocalTimer(*BuildAhatInv);
      
#ifdef PREC_SIMPLEC_2D
      // approximate A^{-1} by hat{A}^{-1} SIMPLE C
      ierr = Ablock->InvRowSums(*Ahatinv);
      check_flag(&ierr, "Ahatinv, InvRowSums", 1);
#endif
      
#ifdef PREC_SIMPLE_2D
      // approximate A^{-1} by hat{A}^{-1} SIMPLE
      ierr = Ablock->ExtractDiagonalCopy(*Ahatinv);
      check_flag(&ierr, "Ahatinv, ExtractDiagonalCopy", 1);
      
      ierr = Ahatinv->Reciprocal(*Ahatinv);
      check_flag(&ierr, "Ahatinv, Reciprocal", 1);
#endif
    }

    // build hat{A}^{-1} * B
    {
      TimeMonitor LocalTimer(*BuildAhatInvB);
            
      // copy global B matrix
      // AhatinvB = rcp(new Epetra_FECrsMatrix(*Bblock));

      ierr = Bblock->LeftScale(*Ahatinv);
      check_flag(&ierr, "LeftScale: AhatinvB", 1);
    } 

    // build W = D * hat{A}^{-1} * B
    {
      TimeMonitor LocalTimer(*BuildDAhatInvB);
      
      ierr = EpetraExt::MatrixMatrix::Multiply(*Dblock, false, *Bblock, false, *Wmatrix, true);
      check_flag(&ierr, "Matrix Multiply: Wmatrix", 1);
    } 

#ifdef DEBUG_PRINT_ON      
      if (printproc) cout << "Computing Schur Complement S" << endl;
#endif

    // build Schur complement: S = E - D A^{-1} B
    {
      TimeMonitor LocalTimer(*BuildSchurS);
      
      Epetra_CrsMatrix* tempptr1 = &(*Smatrix);
      ierr = EpetraExt::MatrixMatrix::Add(*Eblock, false, 1.0, *Wmatrix, false, -1.0, tempptr1);
      check_flag(&ierr, "Matrix Add: Smatrix", 1);
    
      // Finish global matrix assembly (row map = domain map = range map)
      ierr = Smatrix->FillComplete();  
      check_flag(&ierr, "FillComplete: Smatrix", 1);   
    } 
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception: " << ierr << endl;
  }
    
  return ierr;
}


// ------------------------------------------------------------------------------
// Apply inverse of preconditioner, P^{-1} X = Y (Solve P Y = X)
//
// Identity Preconditioner
// ------------------------------------------------------------------------------
#ifdef PREC_IDENT
int BlockPreconditioner_Sparse::ApplyInverse(const Epetra_MultiVector& X, 
                                             Epetra_MultiVector& Y) const
{  
  TimeMonitor LocalTimer(*ApplyTime);

  // Identity Preconditioner
  Y = X;
  
  return 0;
}
#endif


// ------------------------------------------------------------------------------
// Apply inverse of preconditioner, P^{-1} X = Y (Solve P Y = X)
//
// SIMPLE and SIMPLE C Approximate Block Preconditioner on Fluid Variables
// ------------------------------------------------------------------------------
#if defined(PREC_SIMPLE_2D) || defined(PREC_SIMPLEC_2D) 
int BlockPreconditioner_Sparse::ApplyInverse(const Epetra_MultiVector& X, 
                                             Epetra_MultiVector& Y) const
{  
  // start timer
  TimeMonitor LocalTimer(*ApplyTime);

  Belos::ReturnType BelosRet;

  // error flag
  int ierr = 0;

  try{

#ifdef DEBUG_PRINT_ON
    double norm_X; X(0)->Norm1(&norm_X); 
    if(printproc) cout << " norm_X= " << norm_X << flush << endl;
#endif

    // ==========================================================================
    // build global rhs (combine duplicate unknowns)
    // ==========================================================================
  
    // NOTE: X and Y have lenghts nVel which includes locally repeated values but
    //       the vectors VelVec, PsVec, and TVec have length nVel_unique which does
    //       not include the locally repeated values. Thus we have to compensate
    //       by rescaling values added into *Vec_local from X and make sure to 
    //       double copy some values out of *Vec_local into Y.

    {
      TimeMonitor LocalTimer(*BuildRHSTime);

      // --------------------------------------------------------------------------
      // extract velocities from input vector and add into global vector
      // --------------------------------------------------------------------------
      double VelVals[nVel];

      for (int i = 0; i < nVel; i++) {
	VelVals[i] = VelWgt[i] * X[0][i]; 
      }
      
      VelVec_local->PutScalar(0.0);
      
      ierr = VelVec_local->SumIntoGlobalValues(nVel, VelVals, VelIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: VelVec_local",1);
      
      VelVec->PutScalar(0.0);
      
      ierr = VelVec->Export(*VelVec_local, *Vel_Exporter, Add);
      check_flag(&ierr,"Export: VelVec",1);
      
#ifdef DEBUG_PRINT_ON      
      double norm_VelVec; VelVec->Norm1(&norm_VelVec); 
      if(printproc) cout << " norm_VelVec= " << norm_VelVec << flush << endl;
#endif
      
      // --------------------------------------------------------------------------
      // extract surface pressures from input vector and add into global vector
      // --------------------------------------------------------------------------
      double PsVals[nPs];
      
      for (int i = 0; i < nPs; i++) {
	PsVals[i] = PsWgt[i] * X[0][i+nVel]; 
      }
      
      PsVec_local->PutScalar(0.0);
      
      ierr = PsVec_local->SumIntoGlobalValues(nPs, PsVals, PsIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: PsVec_local",1);
      
      PsVec->PutScalar(0.0);
      
      ierr = PsVec->Export(*PsVec_local, *Ps_Exporter, Add);
      check_flag(&ierr,"Export: PsVec",1);
      
#ifdef DEBUG_PRINT_ON      
      double norm_PsVec; PsVec->Norm1(&norm_PsVec); 
      if(printproc) cout << " norm_PsVec= " << norm_PsVec << flush << endl;
#endif
      
      // --------------------------------------------------------------------------
      // extract temperatures from input vector and add into global vector
      // --------------------------------------------------------------------------
      double TVals[nT];
      
      for (int i = 0; i < nT; i++) {
	TVals[i] = TWgt[i] * X[0][i+nVel+nPs];
      }
  
      TVec_local->PutScalar(0.0);
      
      ierr = TVec_local->SumIntoGlobalValues(nT, TVals, TIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: TVec_local",1);
      
      TVec->PutScalar(0.0);
      
      ierr = TVec->Export(*TVec_local, *T_Exporter, Add);
      check_flag(&ierr,"Export: TVec",1);
      
#ifdef DEBUG_PRINT_ON       
      double norm_TVec; TVec->Norm1(&norm_TVec); 
      if(printproc) cout << " norm_TVec= " << norm_TVec << flush << endl;
#endif
    } // close timer scope

    // ==========================================================================
    // apply preconditioner
    // ==========================================================================

    //---------------------------------------------------------------------------
    // Step 1: Solve A z1 = x1
    //---------------------------------------------------------------------------
#ifdef PRINT_SOLVE_INFO
    if(printproc) cout << "ASolve" << flush << endl;
#endif
    {
      TimeMonitor LocalTimer(*ASolveTime);     

      // initialize solution
      Ax->PutScalar(0.0);

      // A block solve
      AProblem->setOperator(Ablock);
      AProblem->setLHS(Ax);
      AProblem->setRHS(VelVec);
      ASolver->reset(Belos::Problem);

      BelosRet = ASolver->solve();
      *ATotalIters += ASolver->getNumIters();
    }

#ifdef DEBUG_PRINT_ON
    if (printproc) {
      if (BelosRet == Belos::Converged) {
	cout << " Belos A solve converged." << flush << endl;
      } else {
	cout << " Belos A solve did not converge." << flush << endl;
      }
    }
  
    double norm_Ax; Ax->Norm1(&norm_Ax);
    if(printproc) cout << " norm_Ax= " << norm_Ax << endl;  
    if(printproc) cout << " ATotalIters= "<< *ATotalIters << flush << endl;
#endif   

    //---------------------------------------------------------------------------
    // Step 2: Solve S z2 = x2 - D z1
    //---------------------------------------------------------------------------
#ifdef PRINT_SOLVE_INFO
    if(printproc) cout << "SSolve" << flush << endl;
#endif
    {
      TimeMonitor LocalTimer(*SSolveTime);     

      // D * z1
      Dblock->Multiply(false, *Ax, *DAx);

      // x2 = x2 - D * z1
      ierr = PsVec->Update(-1.0, *DAx, 1.0);

      // initialize solution
      Sx->PutScalar(0.0);

      // S block solve
      SProblem->setOperator(Smatrix);
      SProblem->setLHS(Sx);
      SProblem->setRHS(PsVec);
      SSolver->reset(Belos::Problem);

      BelosRet = SSolver->solve();
      *STotalIters += SSolver->getNumIters();
    }

#ifdef DEBUG_PRINT_ON
    if (printproc) {
      if (BelosRet == Belos::Converged) {
	cout << " Belos S solve converged." << flush << endl;
      } else {
	cout << " Belos S solve did not converge." << flush << endl;
      }
    }
  
    double norm_Sx; Sx->Norm1(&norm_Sx);
    if(printproc) cout << " norm_Sx= " << norm_Sx << endl;
    if(printproc) cout << " STotalIters= "<< *STotalIters << flush << endl;
#endif   

    //---------------------------------------------------------------------------
    // Step 3: Solve H z3 = x3
    //---------------------------------------------------------------------------
#ifdef PRINT_SOLVE_INFO
    if(printproc) cout << "PSolve" << flush << endl;
#endif
    {
      TimeMonitor LocalTimer(*PSolveTime);     

      // initialize solution
      Px->PutScalar(0.0);

      // P block solve
      PProblem->setOperator(Hblock);
      PProblem->setLHS(Px);
      PProblem->setRHS(TVec);
      PSolver->reset(Belos::Problem);

      BelosRet = PSolver->solve();
      *PTotalIters += PSolver->getNumIters();
    }

#ifdef DEBUG_PRINT_ON
    if (printproc) {
      if (BelosRet == Belos::Converged) {
	cout << " Belos P solve converged." << flush << endl;
      } else {
	cout << " Belos P solve did not converge." << flush << endl;
      }
    }
  
    double norm_Px; Px->Norm1(&norm_Px);
    if(printproc) cout << " norm_Px= " << norm_Px << endl;
    if(printproc) cout << " PTotalIters= "<< *PTotalIters << flush << endl;
#endif   

    //---------------------------------------------------------------------------
    // Step 4: y3 = z3
    //--------------------------------------------------------------------------- 
    // Do nothing

#ifdef DEBUG_PRINT_ON
    double norm_y3; Px->Norm1(&norm_y3); 
    if(printproc) cout << " norm_y3= " << norm_y3 << flush << endl;
#endif
  
    //---------------------------------------------------------------------------
    // Step 5: y2 = z2
    //--------------------------------------------------------------------------- 
    // Do nothing

#ifdef DEBUG_PRINT_ON
    double norm_y2; Sx->Norm1(&norm_y2); 
    if(printproc) cout << " norm_y2= " << norm_y2 << flush << endl;
#endif

    //---------------------------------------------------------------------------
    // Step 6: y1 = z1 - A^{-1}B z2
    //---------------------------------------------------------------------------

    Bblock->Multiply(false, *Sx, *w1);
  
    Ax->Update(-1.0, *w1, 1.0);

#ifdef DEBUG_PRINT_ON
    double norm_y1; Ax->Norm1(&norm_y1); 
    if(printproc) cout << " norm_y1= " << norm_y1 << flush << endl;
#endif

    // ==========================================================================
    // distribute global solution (copy solution to duplicated unknowns)
    // ==========================================================================

    {
      TimeMonitor LocalTimer(*BuildYvecTime);     
      
      // --------------------------------------------------------------------------
      // velocities
      // --------------------------------------------------------------------------
      ierr = VelVec_local->Import(*Ax, *Vel_Importer, Insert);
      check_flag(&ierr, "Import, VelVec_local", 1);
      
#ifdef DEBUG_PRINT_ON      
      double norm_VelVec_local; VelVec_local->Norm1(&norm_VelVec_local); 
      if(printproc) cout << " norm_VelVec_local= " << norm_VelVec_local << flush << endl;
#endif
      
      for (int i = 0; i < nVel_unique; i++) {
	for (int j = 0; j < nVelIdxMap[i]; j++) {
	  Y[0][VelIdxMap[i][j]] = (*VelVec_local)[i]; 
	}
      }
      
      // --------------------------------------------------------------------------
      // surface pressures
      // --------------------------------------------------------------------------
      ierr = PsVec_local->Import(*Sx, *Ps_Importer, Insert);
      check_flag(&ierr, "Import, PsVec_local", 1);

#ifdef DEBUG_PRINT_ON     
      double norm_PsVec_local; PsVec_local->Norm1(&norm_PsVec_local); 
      if(printproc) cout << " norm_PsVec_local= " << norm_PsVec_local << flush << endl;
#endif
      
      for (int i = 0; i < nPs_unique; i++) {
	for (int j = 0; j < nPsIdxMap[i]; j++) {
	  Y[0][PsIdxMap[i][j]] = (*PsVec_local)[i]; 
	}
      }

      // --------------------------------------------------------------------------
      // temperatures
      // --------------------------------------------------------------------------
      ierr = TVec_local->Import(*Px, *T_Importer, Insert);
      check_flag(&ierr, "Import, TVec_local", 1);

#ifdef DEBUG_PRINT_ON      
      double norm_TVec_local; TVec_local->Norm1(&norm_TVec_local); 
      if(printproc) cout << " norm_TVec_local= " << norm_TVec_local << flush << endl;
#endif
      
      for (int i = 0; i < nT_unique; i++) {
	for (int j = 0; j < nTIdxMap[i]; j++) {
	  Y[0][TIdxMap[i][j]] = (*TVec_local)[i];  
	}
      }
            
#ifdef DEBUG_PRINT_ON
      double norm_Y; Y(0)->Norm1(&norm_Y); 
      if(printproc) cout << " norm_Y= " << norm_Y << flush << endl;
#endif
    } // close timer scope
      
  }
  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception: " << ierr << endl;
  }

#ifdef PRINT_SOLVE_INFO
  if(printproc) cout << "PC Complete" << flush << endl;
#endif
      
  return ierr;
}
#endif
