#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosDefines.h" // Preconditioner Options
#include "trilinosUtils.hpp" // Helpful functions

#include "Precon_Operator.hpp"
#include "BlockPreconditioner_Op.hpp"

//#define PREC_IDENT
#define PREC_SIMPLE

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

#ifdef _PRIM

  // preconditioner functions for PRIM
  void prim_Ablock_op(double*, int, double*, void*); 
  void prim_Bblock_op(double*, int, double*, int, void*); 
  void prim_Cblock_op(double*, int, double*, int, void*);  

  void prim_Dblock_op(double*, int, double*, int, void*);  
  void prim_Eblock_op(double*, int, double*, void*);

  // void prim_Fblock_op(double*, int, double*, int, void*);  
  // void prim_Gblock_op(double*, int, double*, int, void*);  
  void prim_Hblock_op(double*, int, double*, void*);

  void prim_SchurS_op(double*, int, double*, void*);

  void prim_AhatInvB_op(double*, int, double*, int, void*);
  void prim_DAhatInv_op(double*, int, double*, int, void*);

  void prim_ident(double*, int, double*, void*); 

  void prim_dti_ident(double*, int, double*, void*); 
  
#else //_SWIM

  // preconditioner functions for SWIM
  void sw_picard_block_11(double*, int, double*, void*); 
  void sw_picard_block_21(double*, int, double*, int, void*); 
  void sw_picard_DFinvBt(double*, int, double*, int, void*);
  void sw_picard_schur(double*, int, double*, void*); 

#endif
}



// ------------------------------------------------------------------------------
// Constructor for Block Preconditioner
// ------------------------------------------------------------------------------
BlockPreconditioner_Op::BlockPreconditioner_Op
( int np_, int nlev_, int nelemd_, 
  RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, 
  void* StateData_, 
  void (*precUpdateFunction_)(double*, int, void*),
  const RCP<ParameterList>&  ASolvePL_,
  const RCP<ParameterList>&  SSolvePL_,
  const RCP<ParameterList>&  PSolvePL_,
  int* ATotalIters_,
  int* STotalIters_,
  int* PTotalIters_)
  : 
  hommePreconditionerBase(xMap_), //Required Base Class construction
  np(np_), nlev(nlev_), nelemd(nelemd_),
  StateData(StateData_),
  precUpdateFunction(precUpdateFunction_),
  ASolvePL(ASolvePL_),
  SSolvePL(SSolvePL_),
  PSolvePL(PSolvePL_),
  ATotalIters(ATotalIters_),
  STotalIters(STotalIters_),
  PTotalIters(PTotalIters_)
{
  const Epetra_Comm& comm = xVec_->Comm();
  printproc = (comm.MyPID()==0);

  if (printproc) 
    cout << " >>> Constructing BlockPreconditioner Op <<< " << endl;

  // create timers
  AssemblyTime = TimeMonitor::getNewCounter("PC: Assemble PC");
  ApplyTime    = TimeMonitor::getNewCounter("PC: Apply PC");
  ASolveTime   = TimeMonitor::getNewCounter("PC: Solve A");
  SSolveTime   = TimeMonitor::getNewCounter("PC: Solve S");
  PSolveTime   = TimeMonitor::getNewCounter("PC: Solve P");

  // array lengths
#ifdef _PRIM
  nState = (3*np*np*nlev + np*np)*(nelemd);
  nVel   = 2*np*np*nlev*nelemd;
  nPs    = np*np*nelemd;
  nT     = np*np*nlev*nelemd;
#else //_SWIM
  nState = 3*np*np*nlev*nelemd;
  nVel   = 2*nState/3;
  nPs    = nState/3;
  nT     = 0;
#endif

  // array maps
  StateMap = rcp(new Epetra_Map(-1, nState, 0, comm)); // state
  VelMap   = rcp(new Epetra_Map(-1, nVel,   0, comm)); // velocity
  PsMap    = rcp(new Epetra_Map(-1, nPs,    0, comm)); // surface pressure
#ifdef _PRIM
  TMap     = rcp(new Epetra_Map(-1, nT,     0, comm)); // temperature
#endif

  bool zeroout=true;

  // temporary working vectors
  w1 = rcp(new Epetra_Vector(*PsMap, zeroout));
  w2 = rcp(new Epetra_Vector(*PsMap, zeroout));
  w3 = rcp(new Epetra_Vector(*VelMap, zeroout));

  // A solve LHS and RHS
  Ax = rcp(new Epetra_Vector(*VelMap, zeroout));
  Ab = rcp(new Epetra_Vector(*VelMap, zeroout));

  // S solve LHS and RHS
  Sx = rcp(new Epetra_Vector(*PsMap, zeroout));
  Sb = rcp(new Epetra_Vector(*PsMap, zeroout));

#ifdef _PRIM
  // P solve LHS and RHS
  Px = rcp(new Epetra_Vector(*TMap, zeroout));
  Pb = rcp(new Epetra_Vector(*TMap, zeroout));
#endif

#ifdef _PRIM
  Solve1 = prim_Ablock_op; // A Solve
  Solve2 = prim_SchurS_op; // S Solve
  Solve3 = prim_Hblock_op; // P Solve
#else // _SWIM
  Solve1 = sw_picard_block_11; // A solve
  Solve2 = sw_picard_schur;    // S solve
#endif
  
  Aop      = rcp(new Precon_Operator(nVel, VelMap, comm, StateData, Solve1));
  AProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Aop, Ax, Ab));
  ASolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(AProblem, ASolvePL));
  
  Sop      = rcp(new Precon_Operator(nPs, PsMap, comm, StateData, Solve2));
  SProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Sop, Sx, Sb));
  SSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(SProblem, SSolvePL));
  
#ifdef _PRIM
  Pop      = rcp(new Precon_Operator(nT, TMap, comm, StateData, Solve3));
  PProblem = rcp(new Belos::LinearProblem<ST,MV,OP>(Pop, Px, Pb));
  PSolver  = rcp(new Belos::BlockGmresSolMgr<ST,MV,OP>(PProblem, PSolvePL));
#endif

#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << "A Solver Parameters   " << endl;
  if (printproc) cout << "Block Size            " << ASolvePL->get<int>("Block Size")<<endl;
  if (printproc) cout << "Num Blocks            " << ASolvePL->get<int>("Num Blocks")<<endl;
  if (printproc) cout << "Maximum Iterations    " << ASolvePL->get<int>("Maximum Iterations")<<endl;
  if (printproc) cout << "Convergence Tolerance " << ASolvePL->get<double>("Convergence Tolerance")<<endl;
  if (printproc) cout << "Output Frequency      " << ASolvePL->get<int>("Output Frequency")<<endl;
  
  if (printproc) cout << "S Solver Parameters   " << endl;
  if (printproc) cout << "Block Size            " << SSolvePL->get<int>("Block Size")<<endl;
  if (printproc) cout << "Num Blocks            " << SSolvePL->get<int>("Num Blocks")<<endl;
  if (printproc) cout << "Maximum Iterations    " << SSolvePL->get<int>("Maximum Iterations")<<endl;
  if (printproc) cout << "Convergence Tolerance " << SSolvePL->get<double>("Convergence Tolerance")<<endl;
  if (printproc) cout << "Output Frequency      " << SSolvePL->get<int>("Output Frequency")<<endl;

#ifdef _PRIM
  if (printproc) cout << "P Solver Parameters   " << endl;
  if (printproc) cout << "Block Size            " << PSolvePL->get<int>("Block Size")<<endl;
  if (printproc) cout << "Num Blocks            " << PSolvePL->get<int>("Num Blocks")<<endl;
  if (printproc) cout << "Maximum Iterations    " << PSolvePL->get<int>("Maximum Iterations")<<endl;
  if (printproc) cout << "Convergence Tolerance " << PSolvePL->get<double>("Convergence Tolerance")<<endl;
  if (printproc) cout << "Output Frequency      " << PSolvePL->get<int>("Output Frequency")<<endl;
#endif
#endif

}



// ------------------------------------------------------------------------------
// Update preconditioner
// ------------------------------------------------------------------------------
int BlockPreconditioner_Op::computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_)
{
  TimeMonitor LocalTimer(*AssemblyTime);

  StateData = StateData_;

  // Update state in preconditioner code
  precUpdateFunction(xVec->Values(), nState, StateData);
  
  const Epetra_Comm& comm = xVec->Comm();
  Aop = rcp(new Precon_Operator(nVel, VelMap, comm, StateData, Solve1));
  Sop = rcp(new Precon_Operator(nPs,  PsMap,  comm, StateData, Solve2));
#ifdef _PRIM
  Pop = rcp(new Precon_Operator(nT,   TMap,   comm, StateData, Solve3));
#endif

  return 0;
}


// ------------------------------------------------------------------------------
// Identity preconditioner for testing
// ------------------------------------------------------------------------------
#ifdef PREC_IDENT
int BlockPreconditioner_Op::ApplyInverse(const Epetra_MultiVector& X, 
                                         Epetra_MultiVector& Y) const
{  
  // start timer
  TimeMonitor LocalTimer(*ApplyTime);

  Y = X;
  return 0;
}
#endif


// ------------------------------------------------------------------------------
// Apply inverse of preconditioner
//
// Applying the SIMPLE Block Preconditioner requries four steps:
// 
// 1. Solve1: A x1 = b1
// 2. Solve2: S x2 = -p B x1 + b2
// 3. Compute y1 = x1 - (1/alpha) A_hat^inv B^T x2, where A_hat is approx of A
// 4. Compute y2 = (1/alpha) x2
// If PRIM:
// 5. Solve3: P x3 = b3
// 
// A, S, and P are implemented as Epetra Operators because they are the LHS in a 
// a linear solve. p B and A_hat^inv B^T can be applied directly as functions. 
//
// Input:  X the RHS of MY = X where M is the preconditioner matirx 
// Output: Y preconditioned vector Y = M^{inv}*X 
// ------------------------------------------------------------------------------
#ifdef PREC_SIMPLE
int BlockPreconditioner_Op::ApplyInverse(const Epetra_MultiVector& X, 
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
    TimeMonitor LocalTimer(*ASolveTime);     
  
    // initialize RHS and solutin vectors
    Ax->PutScalar(0.0);
    
    // extract velocities from input vector
    for (int i = 0; i < nVel; i++) 
      (*Ab)[i] = X[0][i]; 
    
#ifdef DEBUG_PRINT_ON
    double norm_Ab; Ab->Norm2(&norm_Ab);
    if(printproc) cout << " norm_Ab= "<< norm_Ab << flush << endl;
#endif 
    
    // Setup A block solve
    AProblem->setOperator( Aop );
    AProblem->setLHS( Ax );
    AProblem->setRHS( Ab );
    ASolver->reset( Belos::Problem );
  
    
    // Solve A block system
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
  
  double norm_Ax; Ax->Norm2(&norm_Ax);
  if(printproc) cout << " norm_Ax= " << norm_Ax << endl;
  if(printproc) cout << " ATotalIters= "<< *ATotalIters << flush << endl;
#endif   
  
  //---------------------------------------------------------------------------
  // Step 2: Solve S z2 = x2 - D z1
  //---------------------------------------------------------------------------
  //if(printproc) cout << "SSolve" << flush << endl;
  {
    TimeMonitor LocalTimer(*SSolveTime);

    // initialize temporary vectors
    w1->PutScalar(0.0);
    w2->PutScalar(0.0); 
    
    // compute D x1 (block 2,1)
#ifdef _PRIM
    prim_Dblock_op((*Ax).Values(), nVel, (*w1).Values(), nPs, StateData);
#else //_SWIM
    sw_picard_block21((*Ax).Values(), nVel, (*w1).Values(), nPs, StateData);
#endif
    
#ifdef DEBUG_PRINT_ON
    double norm_w1; w1->Norm2(&norm_w1); 
    if(printproc) cout << " norm_w1= "<< norm_w1 << endl;
#endif
    
    // extract surface pressures from input vector
    for (int i = 0; i < nPs; i++) 
      (*w2)[i] = X[0][i+nVel];
    
#ifdef DEBUG_PRINT_ON
    double norm_w2; w2->Norm2(&norm_w2); 
    if(printproc) cout << " norm_w2= "<< norm_w2 << endl;
#endif
    
    // initialize RHS, and solution vectors
    Sb->PutScalar(0.0);
    Sx->PutScalar(0.0);
    
    // Set S Schur RHS, x2 - D x1
    for (int i = 0; i < nPs; i++) 
      (*Sb)[i] = (*w2)[i]-(*w1)[i];
    
#ifdef DEBUG_PRINT_ON
    double norm_Sb; Sb->Norm2(&norm_Sb); 
    if(printproc) cout << " norm_Sb= "<< norm_Sb << flush << endl;
#endif
    
    // Setup S block solve
    SProblem->setOperator( Sop );
    SProblem->setLHS( Sx );
    SProblem->setRHS( Sb );
    SSolver->reset( Belos::Problem ); 
    
    // Solve S system
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
  
  double norm_Sx; Sx->Norm2(&norm_Sx);
  if(printproc) cout << " norm_Sx= " << norm_Sx << endl;
  if(printproc) cout << " STotalIters= "<< *STotalIters << flush << endl;
#endif
  
  //----------------------------------------------------------------------------
  // Step 3: y1 = x1 - (1/alpha) A_hat^inv B^T x2, where A_hat = diag(1/dt)
  // Step 4: y2 = (1/alpha) x2
  //----------------------------------------------------------------------------

  // compute A_hat^inv B^T x2 and store in w1 
  w3->PutScalar(0.0);
  
#ifdef _PRIM
  prim_AhatInvB_op((*Sx).Values(), nPs, (*w3).Values(), nVel, StateData);
#else // _SWIM
  sw_picard_DFinvBt((*Sx).Values(), nPs, (*w3).Values(), nVel, StateData);
#endif
  
#ifdef DEBUG_PRINT_ON
  double norm_w3; w3->Norm2(&norm_w3);
  if(printproc) cout << " norm_w3=" << norm_w3 << flush << endl;
#endif
  
  // relaxation parameter
  double alpha    = 1.0;
  double alphainv = 1.0 / alpha;
  
  for (int i = 0; i < nVel; i++) 
    Y[0][i] = ((*Ax)[i] - alphainv * (*w3)[i]);
  
  for (int i = 0; i < nPs; i++) 
    Y[0][i+nVel] = (alphainv * (*Sx)[i]); 

#ifdef _PRIM
  //----------------------------------------------------------------------------
  // Step 5: Solve H z3 = x3
  //----------------------------------------------------------------------------
  //if(printproc) cout << "PSolve" << flush << endl;
  {
    TimeMonitor LocalTimer(*PSolveTime);     
    
    // initialize solution and RHS vectors
    Px->PutScalar(0.0);
    
    // extract temperature values
    for (int i = 0; i < nT; i++)
      (*Pb)[i] = X[0][i+nVel+nPs]; 
    
#ifdef DEBUG_PRINT_ON
    double norm_Pb; Pb->Norm2(&norm_Pb);
    if(printproc) cout << " norm_Pb= "<< norm_Pb << flush << endl;
#endif
    
    // setup linear system
    PProblem->setOperator( Pop );
    PProblem->setLHS( Px );
    PProblem->setRHS( Pb );
    PSolver->reset( Belos::Problem );
        
    // solve linear system
    BelosRet = PSolver->solve();
    *PTotalIters += PSolver->getNumIters();
  }    
  
#ifdef DEBUG_PRINT_ON
  if (printproc) {
    if (BelosRet == Belos::Converged) {
      cout << "Belos P solve converged." << flush << endl;
    } else {
      cout << "Belos P solve did not converge." << flush << endl;
    }
  }
  
  double norm_Px; Px->Norm2(&norm_Px);
  if(printproc) cout << " norm_Px= " << norm_Px << endl;
  if(printproc) cout << " PTotalIters= "<< *PTotalIters << flush << endl;
#endif   

  // pack output vector
  for (int i = 0; i < nT; i++)
    Y[0][i+nVel+nPs] = (*Px)[i];  

#endif //_PRIM

#ifdef DEBUG_PRINT_ON
  double norm_Y; Y(0)->Norm2(&norm_Y);
  if(printproc) cout << " norm_Y= " << norm_Y << endl;
#endif

  //if(printproc) cout << "PC Complete" << flush << endl;
  
  return 0;
}
#endif
