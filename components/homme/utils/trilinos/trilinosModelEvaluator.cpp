#include "trilinosModelEvaluator.hpp"

// Preconditioner Options
#include "trilinosDefines.h"

// Preconditioners
#include "Precon_Operator.hpp"
#include "IdentityPreconditioner.hpp"
#include "BlockPreconditioner_Op.hpp"
#include "BlockPreconditioner_Dense.hpp"
#include "BlockPreconditioner_Sparse.hpp"
#include "BlockPreconditioner_SparseML.hpp"
#include "GMRESPreconditioner_Sparse.hpp"
#include "GMRESPreconditioner_Op.hpp"
#include "MLPreconditioner.hpp"

// Trilinos
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

//#define DEBUG_PRINT_ON 

using std::cout;
using std::endl;

using Teuchos::RCP;
using Teuchos::rcp;

// ------------------------------------------------------------------------------
// Constructor 
// ------------------------------------------------------------------------------
trilinosModelEvaluator::trilinosModelEvaluator(const int np_, const int nlev_, const int nelemd_,
					       double* StateVector_, void* StateData_,
					       const Epetra_Comm& comm_,
					       void (*residualFunction_)(double*, double*, int, void*),
					       void (*jacFunction_)(double*, int, double*, void*),
					       void (*jacUpdateFunction_)(double*, int, void*),
					       void (*precUpdateFunction_)(double*, int, void*),
					       const RCP<ParameterList>&  ASolvePL_,
					       const RCP<ParameterList>&  SSolvePL_,
					       const RCP<ParameterList>&  PSolvePL_,
					       int* ATotalIt_,
					       int* STotalIt_,
					       int* PTotalIt_)
  : comm(comm_),
    StateData(StateData_),
    residualFunction(residualFunction_)
{ 
  printproc = (comm.MyPID()==0);

#ifdef _PRIM
  nState = (3*np_*np_*nlev_ + np_*np_)*nelemd_;
#else //_SWIM   
  nState = (3*np_*np_*nlev_*nelemd_);
#endif

  initialize(StateVector_);

#if defined(IDENT_PC)
  // Identity preconditioner
  precOp = rcp(new identityPreconditioner(nState, xVec, xMap));

#elif defined(BLOCK_PC_OP)
  // Block preconditioner using HOMME operators
  precOp = rcp(new BlockPreconditioner_Op(np_, nlev_, nelemd_,
                                          xVec, xMap, StateData,
                                          precUpdateFunction_, 
                                          ASolvePL_, SSolvePL_, PSolvePL_,
                                          ATotalIt_, STotalIt_, PTotalIt_)); 

#elif defined(BLOCK_PC_DENSE)
  // Block preconditioner using dense Jacobian blocks
  precOp = rcp(new BlockPreconditioner_Dense(np_, nlev_, nelemd_,
                                             xVec, xMap, StateData,
                                             precUpdateFunction_, 
                                             ASolvePL_, SSolvePL_, PSolvePL_,
                                             ATotalIt_, STotalIt_, PTotalIt_)); 

#elif defined(BLOCK_PC_SPARSE)
  // Block preconditioner using sparse Jacobian blocks
  precOp = rcp(new BlockPreconditioner_Sparse(np_, nlev_, nelemd_,
					      xVec, xMap, StateData,
					      precUpdateFunction_, 
					      ASolvePL_, SSolvePL_, PSolvePL_,
					      ATotalIt_, STotalIt_, PTotalIt_)); 
  
#elif defined(BLOCK_PC_SPARSE_ML)
  // Block preconditioner using sparse Jacobian blocks and ML solver
  precOp = rcp(new BlockPreconditioner_SparseML(np_, nlev_, nelemd_,
						xVec, xMap, StateData,
						precUpdateFunction_, 
						ASolvePL_, SSolvePL_, PSolvePL_,
						ATotalIt_, STotalIt_, PTotalIt_)); 

#elif defined(GMRES_PC_OP)
  // GMRES preconditioner using HOMME operators
  precOp = rcp(new GMRESPreconditioner_Op(np_, nlev_, nelemd_,
                                          xVec, xMap, StateData,
                                          precUpdateFunction_, 
                                          ASolvePL_,
                                          ATotalIt_));

#elif defined(GMRES_PC_SPARSE)
  // GMRES preconditioner using sparse Jacobian
  precOp = rcp(new GMRESPreconditioner_Sparse(np_, nlev_, nelemd_,
					      xVec, xMap, StateData,
					      precUpdateFunction_, 
					      ASolvePL_,
					      ATotalIt_));

#elif defined(ML_PC)
  // ML preconditioner
  precOp = rcp(new MLPreconditioner(np_, nlev_, nelemd_,
                                    xVec, xMap, StateData,
                                    precUpdateFunction_, 
                                    ASolvePL_,
                                    ATotalIt_));

#else
#error Unknown preconditioner option set in trilinosDefines.h
#endif

  // Create the operator for analytic Jacobian-Vector product
  jacOp = rcp(new hommeJacobianVectorProduct(nState, xVec, xMap, StateData, 
					     jacFunction_, jacUpdateFunction_));
}



// destructor
trilinosModelEvaluator::~trilinosModelEvaluator()
{ 
}



// Hooks to reset data from fortran to be passed back to fortran.
void trilinosModelEvaluator::resetInitialGuess(double* StateVector_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tME: Resetting Initial Guess" << endl;
#endif

  for (int i=0; i<nState; i++) (*xVec)[i] = StateVector_[i];
}


// reset pointer to homme derived data type
void trilinosModelEvaluator::resetStateData(void* StateData_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tME: Resetting State Data" << endl;
#endif

  StateData = StateData_; 
}



// ------------------------------------------------------------------------------
// Utility function called by all Constructors above
// ------------------------------------------------------------------------------
void trilinosModelEvaluator::initialize(double* StateVector_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tME: Initializing Maps and Vectors" << endl;
#endif

  bool succeeded=true;
  try {
    xMap = rcp(new Epetra_Map(-1, nState, 0, comm));
    xVec = rcp(new Epetra_Vector(Copy, *xMap, StateVector_));

    pMap = rcp(new Epetra_LocalMap(1, 0, comm));
    pVec = rcp(new Epetra_Vector(*pMap));
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);
}



// Return solution vector map
RCP<const Epetra_Map> trilinosModelEvaluator::get_x_map() const{
  return xMap;
}



// Return residual vector map
RCP<const Epetra_Map> trilinosModelEvaluator::get_f_map() const{
  return xMap;
}



// Return initial solution and x_dot init
RCP<const Epetra_Vector> trilinosModelEvaluator::get_x_init() const{
  return xVec;
}



RCP<EpetraExt::ModelEvaluator::Preconditioner>
trilinosModelEvaluator::create_WPrec() const
{
  // precOp is already constructed.
  //   bool is answer to: "Prec is already inverted?"
  return rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,true));
}



// When NOX asks for Jacobian operator, return jacOp
RCP<Epetra_Operator> trilinosModelEvaluator::create_W() const{
  return jacOp;
}



RCP<const Epetra_Map> trilinosModelEvaluator::get_p_map(int l) const{
  return pMap;
}



RCP<const Epetra_Vector> trilinosModelEvaluator::get_p_init(int l) const{
  return pVec;
}



RCP<const  Teuchos::Array<std::string> >  trilinosModelEvaluator::get_p_names(int l) const{
  RCP<Teuchos::Array<std::string> > p_names =
    rcp(new Teuchos::Array<std::string>(1) );
  (*p_names)[0] = "LOCAParameter";

  return p_names;
}



// Create InArgs
EpetraExt::ModelEvaluator::InArgs trilinosModelEvaluator::createInArgs() const{
  InArgsSetup inArgs;

  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);
  return inArgs;
}



// Create OutArgs
EpetraExt::ModelEvaluator::OutArgs trilinosModelEvaluator::createOutArgs() const{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_WPrec, true);
  // Declare that this interface can compute an analytic Jacobian, W
  outArgs.setSupports(OUT_ARG_W, true);

  return outArgs;
}



// Evaluate model on InArgs
void trilinosModelEvaluator::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const{

  // Get the solution vector x from inArgs and residual vector from outArgs
  RCP<const Epetra_Vector> x = inArgs.get_x();
  if (x == Teuchos::null) throw "trilinosModelEvaluator::evalModel: x was NOT specified!";

  // Check what is being computed in the evalModel call
  bool residualRequested=false, preconditionerRequested=false, jacRequested=false;

  // An evaluation of the residual has been requested by NOX
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f = outArgs.get_f();
  if (f != Teuchos::null) residualRequested=true;

  // A recompute of the preconditioner has been requested by NOX
  RCP<Epetra_Operator> WPrec = outArgs.get_WPrec();
  if (WPrec != Teuchos::null) preconditionerRequested=true;

  // A recompute of the Jacobian has been requested by NOX
  RCP<Epetra_Operator> W = outArgs.get_W();
  if (W != Teuchos::null) jacRequested=true;

  // Code for parameter continuation with LOCA; not set up yet
  RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  //  if (p_in.get()) set_parameter(&(*p_in)[0]);
  TEUCHOS_TEST_FOR_EXCEPTION(p_in.get(), std::logic_error,
			     "Parameter being set in Model Evaluator, but not implemented in code.");

  // Save the current solution, which makes it initial guess for next nonlinear solve
  *xVec = *x;

  if (residualRequested) {
    // Check if this is a perturbed eval. CESM only saves off matrices for unperturbed case.
    // int ispert =0; if  (f.getType() == EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV) ispert=1;

    f->PutScalar(0.0);
    //double nrm; x->Norm2(&nrm); cout << "EvalModel Norm x = " << nrm << endl;
    residualFunction(x->Values(), f->Values(), nState, StateData);
  }

  // recompute preconditioner with new vector, data
  if (preconditionerRequested) {
    precOp->computePreconditioner(x, StateData);
  }

  // recompute Jacobian now with new vector, data
  if (jacRequested) {
    jacOp->computeJacobian(x, StateData);
  }
}



// ==============================================================================
// Functions for analytic Jacobian vector product
// ==============================================================================



// ------------------------------------------------------------------------------
// constructor
// ------------------------------------------------------------------------------
hommeJacobianVectorProduct::hommeJacobianVectorProduct (int nState_, 
							RCP<Epetra_Vector> xVec_, 
							RCP<Epetra_Map> xMap_,
							void* StateData_,
							void (*jacFunction_)(double*, int, double*, void*),
							void (*jacUpdateFunction_)(double*, int, void*))
  : nState(nState_),
    xVec(xVec_),
    xMap(xMap_),
    StateData(StateData_),
    jacFunction(jacFunction_),
    jacUpdateFunction(jacUpdateFunction_)
{
  printproc = (xVec->Comm().MyPID() == 0);

#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tJVP: Constructing Analytic Jacobian <<<" << endl;
#endif
}



// ------------------------------------------------------------------------------
// update Jacobian matrix
// ------------------------------------------------------------------------------
int hommeJacobianVectorProduct::computeJacobian(RCP<const Epetra_Vector> xVecNew, void* StateData_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tJVP: Computing Analytic Jacobian <<<" << endl;
#endif

  // update pointer to state data
  StateData = StateData_;

  // update state in State Data 
  jacUpdateFunction(xVecNew->Values(), nState, StateData);

  return 0;
}



// ------------------------------------------------------------------------------
// compute Jacobian-vector product
// ------------------------------------------------------------------------------
int hommeJacobianVectorProduct::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> tJVP: Applying Analytic Jacobian <<<" << endl;
#endif

  jacFunction(X.Values(), nState, Y.Values(), StateData);

  return 0;
}
