#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosDefines.h" // Preconditioner Options
#include "trilinosUtils.hpp" // Helpful functions

#include "GMRESPreconditioner_Sparse.hpp"

//#define DEBUG_PRINT_ON

using std::cout;
using std::endl;
using std::flush;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// ==============================================================================
// Constructor for GMRES Preconditioner
// ==============================================================================
GMRESPreconditioner_Sparse::GMRESPreconditioner_Sparse
( const int np_, const int nlev_, const int nelemd_,
  RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, 
  void* StateData_, 
  void (*precUpdateFunction_)(double*, int, void*),
  const RCP<ParameterList>& JSolvePL_,
  int* JTotalIters_)
  : 
  SparsePreconditionerBase(np_, nlev_, nelemd_, StateData_, xVec_, xMap_),
  precUpdateFunction(precUpdateFunction_),
  JSolvePL(JSolvePL_),
  JTotalIters(JTotalIters_),
  firstfill(true)
{
  if (printproc) cout << ">>> Constructing GMRES Preconditioner <<<" << endl;

  // MPI communicator
  const Epetra_Comm& Comm = xVec_->Comm();
  
  // assembly timers
  ComputePCTime      = TimeMonitor::getNewCounter("PC: Compute PC");
  UpdateStateTime    = TimeMonitor::getNewCounter("PC: Update State");
  ZeroTime           = TimeMonitor::getNewCounter("PC: Zero PC");
  FillTime           = TimeMonitor::getNewCounter("PC: Fill PC");
  GlobalAssemblyTime = TimeMonitor::getNewCounter("PC: Global Assembly");

  // apply timers
  ApplyTime          = TimeMonitor::getNewCounter("PC: Apply PC");
  BuildRHSTime       = TimeMonitor::getNewCounter("PC: Build RHS Vec");
  JSolveTime         = TimeMonitor::getNewCounter("PC: Solve J");
  BuildYvecTime      = TimeMonitor::getNewCounter("PC: Build Y Vec");

  // combine indices for V, Ps, and T into one state vector
  int StateIdx[nState_unique];
  
  for (int i=0; i < nVel_unique; i++) {
    StateIdx[i] = VelIdx[i];
  }
  
  for (int i=0; i < nPs_unique; i++) {
    StateIdx[i+nVel_unique] = PsIdx[i];
  }
  
  for (int i=0; i < nT_unique; i++) {
    StateIdx[i+nVel_unique+nPs_unique] = TIdx[i];
  }
  
  // local and global Epetra Maps
  StateMap_local = rcp(new Epetra_Map(-1, nState_unique, StateIdx, 0, Comm)); 
  StateMap       = rcp(new Epetra_Map(Epetra_Util::Create_OneToOne_Map(*StateMap_local)));

  // exporter (source map = local, target map = global)
  State_Exporter = rcp(new Epetra_Export(*StateMap_local, *StateMap));

  // importer (target map = local, source map = global)
  State_Importer = rcp(new Epetra_Import(*StateMap_local, *StateMap));

  // local and global Epetra solution vectors
  bool zeroout=true;

  StateVec_local = rcp(new Epetra_Vector(*StateMap_local, zeroout));
  StateVec       = rcp(new Epetra_Vector(*StateMap,       zeroout));
  Jx             = rcp(new Epetra_Vector(*StateMap, zeroout));

  // Jacobian matrix
  int nnz_per_row = 2*(np+np-1)*nlev + np+np-1 + (np+np-1)*nlev;
  Jmatrix = rcp(new Epetra_FECrsMatrix(Copy, *StateMap, 2*nnz_per_row));

  // create Belos linear problem (Op, LHS, RHS) and solver manager
  JProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Jmatrix, Jx, StateVec));
  JSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(JProblem, JSolvePL));
}



// ==============================================================================
// Update preconditioner
// ==============================================================================
int GMRESPreconditioner_Sparse::computePreconditioner(RCP<const Epetra_Vector> xVec, 
                                                      void* StateData_)
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

    { 
      TimeMonitor LocalTimer(*ZeroTime);

      if (Jmatrix->Filled()) {
        Jmatrix->Scale(0.0);
      } else {
        ZeroAblock(Jmatrix);
        ZeroBblock(Jmatrix);
        ZeroCblock(Jmatrix);
        ZeroDblock(Jmatrix);
        ZeroEblock(Jmatrix);
        ZeroFblock(Jmatrix);
        ZeroGblock(Jmatrix);
        ZeroHblock(Jmatrix);
      }
    }

    {
      TimeMonitor LocalTimer(*FillTime);
      
      FillAblock(Jmatrix);
      FillBblock(Jmatrix);
      FillCblock(Jmatrix);
      FillDblock(Jmatrix);
      FillEblock(Jmatrix);
      FillFblock(Jmatrix);
      FillGblock(Jmatrix);
      FillHblock(Jmatrix);
    }
    
    {
      TimeMonitor LocalTimer(*GlobalAssemblyTime);
      
      // Finish global matrix assembly (row map = domain map = range map)
      ierr = Jmatrix->GlobalAssemble(firstfill, Add, true);
      check_flag(&ierr, "Jmatrix, Global FillComplete", 1);

      if (firstfill) firstfill = false;
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
// ------------------------------------------------------------------------------
int GMRESPreconditioner_Sparse::ApplyInverse(const Epetra_MultiVector& X,
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
    //       the vectors VelVec, PsVec, and TVec have length NVel which does 
    //       not include the locally repeated values. Thus we have to compensate
    //       by rescaling values added into *Vec_local from X and make sure to 
    //       double copy some values out of *Vec_local into Y.
    {
      TimeMonitor LocalTimer(*BuildRHSTime);

      StateVec_local->Scale(0.0);

      // Velocities
      double VelVals[nVel];

      for (int i=0; i < nVel; i++)
	VelVals[i] = VelWgt[i] * X[0][i]; 

      ierr = StateVec_local->SumIntoGlobalValues(nVel, VelVals, VelIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: VelVec_local Velocities",1);

      // Surface pressure
      double PsVals[nPs];

      for (int i = 0; i < nPs; i++)
	PsVals[i] = PsWgt[i] * X[0][i+nVel];

      ierr = StateVec_local->SumIntoGlobalValues(nPs, PsVals, PsIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: StateVec_local Surface Pressure",1);

      // Temperatures
      double TVals[nT];
      
      for (int i = 0; i < nT; i++)
	TVals[i] = TWgt[i] * X[0][i+nVel+nPs]; 
      
      ierr = StateVec_local->SumIntoGlobalValues(nT, TVals, TIdx_unsorted);
      check_flag(&ierr,"SumIntoGlobalValues: StateVec_local Temperature",1);

      // Export to local to global
      StateVec->Scale(0.0);

      ierr = StateVec->Export(*StateVec_local, *State_Exporter, Add);
      check_flag(&ierr,"Export: StateVec",1);
            
#ifdef DEBUG_PRINT_ON      
      double norm_StateVec; StateVec->Norm1(&norm_StateVec); 
      if(printproc) cout << " norm_StateVec= " << norm_StateVec << flush << endl;
#endif
    }

    // ==========================================================================
    // apply preconditioner
    // ==========================================================================
    //if(printproc) cout << "JSolve" << flush << endl;
    {
      TimeMonitor LocalTimer(*JSolveTime);     

      // initialize solution
      Jx->PutScalar(0.0);

      // A block solve
      JProblem->setOperator(Jmatrix);
      JProblem->setLHS(Jx);
      JProblem->setRHS(StateVec);
      JSolver->reset(Belos::Problem);

      BelosRet = JSolver->solve();
      *JTotalIters += JSolver->getNumIters();
    }

#ifdef DEBUG_PRINT_ON
    if (printproc) {
      if (BelosRet == Belos::Converged) {
	cout << " Belos A solve converged." << flush << endl;
      } else {
	cout << " Belos A solve did not converge." << flush << endl;
      }
    }
  
    double norm_Jx; Jx->Norm1(&norm_Jx);
    if(printproc) cout << " norm_Jx= " << norm_Jx << endl;
  
    if(printproc) cout << " JTotalIters= "<< *JTotalIters << flush << endl;
#endif   

    // ==========================================================================
    // distribute global solution (copy solution to duplicated unknowns)
    // ==========================================================================
    {
      TimeMonitor LocalTimer(*BuildYvecTime);     

      // initialize solution
      //StateVec_local->PutScalar(0.0);
      
      ierr = StateVec_local->Import(*Jx, *State_Importer, Insert);
      check_flag(&ierr, "Import: StateVec_local", 1);
      
      // velocities            
      for (int i = 0; i < nVel_unique; i++) {
	for (int j = 0; j < nVelIdxMap[i]; j++) {
	  Y[0][VelIdxMap[i][j]] = (*StateVec_local)[i]; 
	}
      }

      // surface pressures
      for (int i = 0; i < nPs_unique; i++) {
	for (int j = 0; j < nPsIdxMap[i]; j++) {
	  Y[0][PsIdxMap[i][j]] = (*StateVec_local)[i+nVel_unique]; 
	}
      }
      
      // temperatures
      for (int i = 0; i < nT_unique; i++) {
	for (int j = 0; j < nTIdxMap[i]; j++) {
	  Y[0][TIdxMap[i][j]] = (*StateVec_local)[i+nVel_unique+nPs_unique];  
	}
      }
      
#ifdef DEBUG_PRINT_ON
      double norm_Y; Y(0)->Norm1(&norm_Y); 
      if(printproc) cout << " norm_Y= " << norm_Y << flush << endl;
#endif
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

  //if(printproc) cout << "PC Complete" << flush << endl;
  
  return ierr;
}
