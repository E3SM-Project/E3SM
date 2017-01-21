#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosDefines.h" // Preconditioner Options
#include "trilinosUtils.hpp" // Helpful functions

#include "Precon_Operator.hpp"
#include "GMRESPreconditioner_Op.hpp"

//#define DEBUG_PRINT_ON

using std::cout;
using std::endl;
using std::flush;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// ------------------------------------------------------------------------------
// Prototypes for Fortran functions
// ------------------------------------------------------------------------------
extern "C" {

  // preconditioner functions for PRIM
  void prim_approx_jacobian_op(double*, int, double*, void*); 
  void prim_approx_jacobian_op_2(double*, int, double*, void*); 
  void prim_approx_jacobian_op_3(double*, int, double*, void*); 
}

// ------------------------------------------------------------------------------
// Constructor for GMRES Preconditioner
// ------------------------------------------------------------------------------
GMRESPreconditioner_Op::GMRESPreconditioner_Op
( int np_, int nlev_, int nelemd_, 
  RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, 
  void* StateData_, 
  void (*precUpdateFunction_)(double*, int, void*),
  const RCP<ParameterList>&  JSolvePL_,
  int* JTotalIters_)
  : 
  hommePreconditionerBase(xMap_), //Required Base Class construction
  np(np_), nlev(nlev_), nelemd(nelemd_),
  StateData(StateData_),
  precUpdateFunction(precUpdateFunction_),
  JSolvePL(JSolvePL_),
  JTotalIters(JTotalIters_)
{
  const Epetra_Comm& comm = xVec_->Comm();
  printproc = (comm.MyPID()==0);

  if (printproc) 
    cout << " >>> GMRESPreconditioner Op: Constructing Preconditioner <<< " << endl;

  // create timers
  AssemblyTime = TimeMonitor::getNewCounter("PC: Assemble PC");
  ApplyTime    = TimeMonitor::getNewCounter("PC: Apply PC");
  JSolveTime   = TimeMonitor::getNewCounter("PC: Solve A");

  // array lengths
  nState = (3*np*np*nlev + np*np)*(nelemd);
  nVel   = 2*np*np*nlev*nelemd;
  nPs    = np*np*nelemd;
  nT     = np*np*nlev*nelemd;

  StateMap = rcp(new Epetra_Map(-1, nState, 0, comm)); // state map

  bool zeroout=true;

  // J solve LHS and RHS
  Jx = rcp(new Epetra_Vector(*StateMap, zeroout));
  Jb = rcp(new Epetra_Vector(*StateMap, zeroout));

  // set which jacobian approximation to use
  Solve1 = prim_approx_jacobian_op;
  
  Jop      = rcp(new Precon_Operator(nState, StateMap, comm, StateData, Solve1));
  JProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Jop, Jx, Jb));
  JSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(JProblem, JSolvePL));
}



// ------------------------------------------------------------------------------
// Update preconditioner
// ------------------------------------------------------------------------------
int GMRESPreconditioner_Op::computePreconditioner(RCP<const Epetra_Vector> xVec, 
                                                  void* StateData_)
{
  TimeMonitor LocalTimer(*AssemblyTime);

  StateData = StateData_;

  // Update state in preconditioner code
  precUpdateFunction(xVec->Values(), nState, StateData);
  
  const Epetra_Comm& comm = xVec->Comm();
  Jop = rcp(new Precon_Operator(nState, StateMap, comm, StateData, Solve1));

  return 0;
}



// ------------------------------------------------------------------------------
// Apply inverse of preconditioner, P^{-1} X = Y
// ------------------------------------------------------------------------------
int GMRESPreconditioner_Op::ApplyInverse(const Epetra_MultiVector& X, 
                                         Epetra_MultiVector& Y) const
{  
  // start timer
  TimeMonitor LocalTimer(*ApplyTime);

  Belos::ReturnType BelosRet;

#ifdef DEBUG_PRINT_ON
  if (printproc) cout << "In ApplyInverse" << flush << endl;
#endif 
  
#ifdef DEBUG_PRINT_ON
  double norm_X; X(0)->Norm2(&norm_X); 
  if(printproc) cout << " norm_X= " << norm_X << flush << endl;
#endif
  
  //----------------------------------------------------------------------------
  // Step 1: Solve A x1 = b
  //----------------------------------------------------------------------------
  //if(printproc) cout << "ASolve" << flush << endl;
  {
    TimeMonitor LocalTimer(*JSolveTime);     
    
    // initialize RHS and solutin vectors
    Jx->PutScalar(0.0);
    
    // extract velocities from input vector
    for (int i = 0; i < nState; i++) 
      (*Jb)[i] = X[0][i]; 
    
#ifdef DEBUG_PRINT_ON
    double norm_Jb; Jb->Norm2(&norm_Jb);
    if(printproc) cout << " norm_Jb= "<< norm_Jb << flush << endl;
#endif 
    
    // Setup A block solve
    JProblem->setOperator( Jop );
    JProblem->setLHS( Jx );
    JProblem->setRHS( Jb );
    JSolver->reset( Belos::Problem );
    
    // Solve A block system
    BelosRet = JSolver->solve();
    *JTotalIters += JSolver->getNumIters();
  }
  
#ifdef DEBUG_PRINT_ON
  if (printproc) {
    if (BelosRet == Belos::Converged) {
      cout << " Belos J solve converged." << flush << endl;
    } else {
      cout << " Belos J solve did not converge." << flush << endl;
    }
  }
  
  double norm_Jx; Jx->Norm2(&norm_Jx);
  if(printproc) cout << " norm_Jx= " << norm_Jx << endl;
  if(printproc) cout << " JTotalIters= "<< *JTotalIters << flush << endl;
#endif   
  
  for (int i = 0; i < nState; i++) 
    Y[0][i] = (*Jx)[i];  
  
#ifdef DEBUG_PRINT_ON
  double norm_Y; Y(0)->Norm2(&norm_Y);
  if(printproc) cout << " norm_Y= " << norm_Y << endl;
#endif
  
  //if(printproc) cout << "PC Complete" << flush << endl;
  
  return 0;
}
