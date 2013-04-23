#include "trilinosModelEvaluator.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// Need all these still?

#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_LinearProblem.h"

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <AztecOO.h>
#include "precon_interface.hpp"
#include "block_precon_interface.hpp"


#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>


#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

//#define DEBUG_PRINT_ON 
//#define SIMPLE_PREC
#define ID_PREC

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

//Identity preconditioner: 
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
    const Epetra_Comm& comm_,
    void* blackbox_res_, void* precdata_,
    void (*residualFunction_)(double *, double *, int, void *),
    void (*precUpdateFunction_)(double *, int, void *) ) :
  N(nelems),
  comm(comm_),
  blackbox_res(blackbox_res_),
  precdata(precdata_),
  residualFunction(residualFunction_),
  precUpdateFunction(precUpdateFunction_)
{ 
  initialize(statevector);
  precOp = Teuchos::rcp(new identityPreconditioner(N, xVec, xMap, blackbox_res,
                                                   precdata, precUpdateFunction) );
}

// Analytic Jacobian: AN_JAC_SCALAR_PREC_ON 
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
    const Epetra_Comm& comm_,
    void* blackbox_res_, void* precdata_, void* jacdata_,
    void (*residualFunction_)(double *, double *, int, void *),
    void (*precFunction_)(double *, int, double*, void *),
    void (*jacFunction_)(double *, int, double*, void *),
    void (*precUpdateFunction_)(double *, int, void *),
    void (*getJacVector_)(double *, int, void *)) :
  N(nelems),
  comm(comm_),
  blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
  residualFunction(residualFunction_),
  precFunction(precFunction_),
  jacFunction(jacFunction_),
  precUpdateFunction(precUpdateFunction_),
  getJacVector(getJacVector_)
{ 
  initialize(statevector);
  throw "No Preconditioner constructed in this modelEval AAA";
}


//SIMPLE preconditioner:  SIMPLE_PREC_ON
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
    const Epetra_Comm& comm_,
    void* blackbox_res_, void* precdata_, void* jacdata_,
    void (*residualFunction_)(double *, double *, int, void *),
    void (*precFunctionblock11_)(double *, int, double*, void *),
    void (*precFunctionblock12_)(double *, int, double*, void *),
    void (*precFunctionblock21_)(double *, int, double*, void *),
    void (*precFunctionblock22_)(double *, int, double*, void *),
    void (*jacFunction_)(double *, int, double*, void *),
    void (*precUpdateFunction_)(double *, int, void *),
    void (*getJacVector_)(double *, int, void *)) :
  N(nelems),
  comm(comm_),
  blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
  residualFunction(residualFunction_),
  precFunctionblock11(precFunctionblock11_),
  precFunctionblock12(precFunctionblock12_),
  precFunctionblock21(precFunctionblock21_),
  precFunctionblock22(precFunctionblock22_),
  jacFunction(jacFunction_),
  precUpdateFunction(precUpdateFunction_),
  getJacVector(getJacVector_)
{ 
  initialize(statevector);
  throw "No Preconditioner constructed in this modelEval BBB";
}

/* This interface is just for testing... comparing two block preconditioner formulations */

//two different formulations of SIMPLE: COMPARE_SIMPLE_BLOCK_VS_SEGGREGATED_ON
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
                const Epetra_Comm& comm_,
                void* blackbox_res_, void* precdata_, void* jacdata_,
                void (*residualFunction_)(double *, double *, int, void *),
                void (*precFunctionblock11_)(double *, int, double*, void *),
                void (*precFunctionblock12_)(double *, int, double*, void *),
                void (*precFunctionblock21_)(double *, int, double*, void *),
                void (*precFunctionblock22_)(double *, int, double*, void *),
                void (*auxprecFunctionblock11_)(double *, int, double*, void *),
                void (*auxprecFunctionblock12_)(double *, int, double*, void *),
                void (*auxprecFunctionblock21_)(double *, int, double*, void *),
                void (*auxprecFunctionblock22_)(double *, int, double*, void *),
                void (*jacFunction_)(double *, int, double*, void *),
                void (*precUpdateFunction_)(double *, int, void *),
                void (*getJacVector_)(double *, int, void *)) :
        N(nelems),
        comm(comm_),
        blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
        residualFunction(residualFunction_),
        precFunctionblock11(precFunctionblock11_),
        precFunctionblock12(precFunctionblock12_),
        precFunctionblock21(precFunctionblock21_),
        precFunctionblock22(precFunctionblock22_),
        auxprecFunctionblock11(auxprecFunctionblock11_),
        auxprecFunctionblock12(auxprecFunctionblock12_),
        auxprecFunctionblock21(auxprecFunctionblock21_),
        auxprecFunctionblock22(auxprecFunctionblock22_),
        jacFunction(jacFunction_),
        precUpdateFunction(precUpdateFunction_),
        getJacVector(getJacVector_)
{
  initialize(statevector);
  throw "No Preconditioner constructed in this modelEval CCC";
}

// Scalar Preconditioner: FD_JAC_SCALAR_PREC_ON 
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
    const Epetra_Comm& comm_,
    void* blackbox_res_, void* precdata_,
    void (*residualFunction_)(double *, double *, int, void *),
    void (*precFunction_)(double *, int, double*, void *),
    void (*precUpdateFunction_)(double *, int, void *)) :
  N(nelems),
  comm(comm_),
  blackbox_res(blackbox_res_),
        precdata(precdata_),
  residualFunction(residualFunction_),
  precFunction(precFunction_),
  precUpdateFunction(precUpdateFunction_)
{ 
  initialize(statevector);
  throw "No Preconditioner constructed in this modelEval DDD";
}

trilinosModelEvaluator::~trilinosModelEvaluator()
{ 
}


// Hooks to reset data from fortran to be passed back to fortran.
void trilinosModelEvaluator::resetInitialGuess(double* statevector_){
  for (int i=0; i<N; i++) (*xVec)[i] = statevector_[i];
}

void trilinosModelEvaluator::resetBlackbox(void* blackbox_res_,  void* precdata_,void* jacdata_){
  blackbox_res=blackbox_res_; 
  precdata=precdata_;
  jacdata=jacdata_;
}

void trilinosModelEvaluator::resetBlackbox(void* blackbox_res_,  void* precdata_){
  blackbox_res=blackbox_res_; 
  precdata=precdata_;
}

void trilinosModelEvaluator::resetBlackbox(void* blackbox_res_){
  blackbox_res=blackbox_res_; 
}


//-----------------------------------------------------------------------------

// Utility function called by all Constructors above
void trilinosModelEvaluator::initialize(double* statevector)
{
  bool succeeded=true;
  try {
    xMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));
    xVec = Teuchos::rcp(new Epetra_Vector(Copy, *xMap, statevector));

    pMap = Teuchos::rcp(new Epetra_LocalMap(1, 0, comm));
    pVec = Teuchos::rcp(new Epetra_Vector(*pMap));

    if (comm.MyPID()==0) printproc=true;
    else                 printproc=false;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);
}

/*******************************************************************************/
// Return solution vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_x_map() const{
  return xMap;
}

// Return residual vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_f_map() const{
  return xMap;
}

// Return initial solution and x_dot init
Teuchos::RCP<const Epetra_Vector> trilinosModelEvaluator::get_x_init() const{
  return xVec;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
trilinosModelEvaluator::create_WPrec() const
{
  // precOp is already constructed.
  //   bool is answer to: "Prec is already inverted?"
  return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,true));
}

Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_p_map(int l) const{
  return pMap;
}
Teuchos::RCP<const Epetra_Vector> trilinosModelEvaluator::get_p_init(int l) const{
  return pVec;
}
Teuchos::RCP<const  Teuchos::Array<std::string> >  trilinosModelEvaluator::get_p_names(int l) const{
    RCP<Teuchos::Array<std::string> > p_names =
      rcp(new Teuchos::Array<std::string>(1) );
    (*p_names)[0] = "LOCAParameter";

  return p_names;
}

/*******************************************************************************/
// Create InArgs
EpetraExt::ModelEvaluator::InArgs trilinosModelEvaluator::createInArgs() const{
  InArgsSetup inArgs;

  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1);
  return inArgs;
}

/*******************************************************************************/
// Create OutArgs
EpetraExt::ModelEvaluator::OutArgs trilinosModelEvaluator::createOutArgs() const{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_WPrec, true);

  return outArgs;
}

/*******************************************************************************/
// Evaluate model on InArgs
void trilinosModelEvaluator::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const{

  // Get the solution vector x from inArgs and residual vector from outArgs
  RCP<const Epetra_Vector> x = inArgs.get_x();
    if (x == Teuchos::null) throw "trilinosModelEvaluator::evalModel: x was NOT specified!";

  // Check what is being computed in the evalModel call
  bool residualRequested=false, preconditionerRequested=false;
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f = outArgs.get_f();
    if (f != Teuchos::null) residualRequested=true;
  RCP<Epetra_Operator> WPrec = outArgs.get_WPrec();
    if (WPrec != Teuchos::null) preconditionerRequested=true;
  

  // Code for parameter continuation with LOCA; not set up yet
  Teuchos::RCP<const Epetra_Vector> p_in = inArgs.get_p(0);
  //  if (p_in.get()) set_parameter(&(*p_in)[0]);
      TEUCHOS_TEST_FOR_EXCEPTION(p_in.get(), std::logic_error,
          "Parameter being set in Model Evaluator, but not implemented in code.");


  // Save the current solution, which makes it initial guess for next nonlinear solve
  *xVec = *x;

  if (residualRequested) {
    // Check if this is a perturbed eval. Glimmer only saves off matrices for unperturbed case.
    // int ispert =0; if  (f.getType() == EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV) ispert=1;

    f->PutScalar(0.0);
    //double nrm; x->Norm2(&nrm); cout << "EvalModel Norm x = " << nrm << endl;
    residualFunction(x->Values(), f->Values(), N, blackbox_res);
  }

  if (preconditionerRequested) {
    precOp->computePreconditioner(x);
  }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
identityPreconditioner::identityPreconditioner (
       int N_, RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_,
       void* blackbox_res_, void* precdata_,
       void (*precUpdateFunction_)(double *, int, void *) )
       : hommePreconditionerBase(xMap_), //Required Base Class construction
         N(N_), xVec(xVec_), xMap(xMap_),
         blackbox_res(blackbox_res_), precdata(precdata_),
         precUpdateFunction(precUpdateFunction_)

{
}

int identityPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Y = X;
  return 0;
}

int identityPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew)
{
  // Copy new values of current solution
  // Not Needed? xVec never used
//  *xVec = *xVecNew; 

  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);

  return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
// Add other preconditioners here. Only contructor and 2 methods need to be writen.
// All preconditioners must have xMap_ in constructor's argument list
//    and the following line to build the Base class.
//        : hommePreconditionerBase(xMap_), //Required Base Class construction
