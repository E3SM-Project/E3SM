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

#define DEBUG_PRINT_ON 

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
  precOp = Teuchos::rcp(new simplePreconditioner(N, xVec, xMap, blackbox_res,
                                                 precdata, 
                                                 precFunctionblock11_, precFunctionblock12_,
                                                 precFunctionblock21_, precFunctionblock22_,
                                                 precUpdateFunction) );
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
    // Check if this is a perturbed eval. CESM only saves off matrices for unperturbed case.
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
  bool printproc = false;
  if (xVec->Comm().MyPID() == 0) printproc=true;
  if (printproc) cout << "Constructing preconditioner:  identityPreconditioner" << endl;
}

int identityPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew)
{
  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);


  return 0;
}

int identityPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Y = X;
  return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
simplePreconditioner::simplePreconditioner (
       int N_, RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_,
       void* blackbox_res_, void* precdata_,
       void (*precFunctionblock11_)(double *, int, double*, void *),
       void (*precFunctionblock12_)(double *, int, double*, void *),
       void (*precFunctionblock21_)(double *, int, double*, void *),
       void (*precFunctionblock22_)(double *, int, double*, void *),
       void (*precUpdateFunction_)(double *, int, void *) )
       : hommePreconditionerBase(xMap_), //Required Base Class construction
         N(N_), xVec(xVec_), xMap(xMap_),
         blackbox_res(blackbox_res_), precdata(precdata_),
         precFunctionblock11(precFunctionblock11_),
         precFunctionblock12(precFunctionblock12_),
         precFunctionblock21(precFunctionblock21_),
         precFunctionblock22(precFunctionblock22_),
         precUpdateFunction(precUpdateFunction_)
{
  const Epetra_Comm& comm = xVec->Comm();
  if (comm.MyPID()==0) printproc=true;
  else                 printproc=false;

  if (printproc) cout << "Constructing preconditioner:  simplePreconditioner" << endl;

  F=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock11));
  S=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock22));
}

int simplePreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew)
{
  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);

  return 0;
}


/*use this applyinverse for applying the simple preconditioner for our tests*/
/* SIMPLE  */

//ApplyInverse routine for the SIMPLE preconditioner
//Convention: for SIMPLE algorithm based on Picard linearization
//11 block is F
//12 block is diag(F)^{-1}B'
//21 block is B
//22 block is S=G-Bdiag(F)^{-1}B'

// These operators receive and return complete state vectors. 
// Components where the operator is not defined will return a zero 
// Specifically,
// F is applied to velocities and returns velocities along with a height field of 0
// similarly S is applied to a height field and returns a height field along with a velocity field of 0
// B is applied to a velocity field and returns a scalar height feild along with a zero velocity field
// diag(F)^{-1}B' is applied to a height field and returns a velocity field with a zero height field.
// communication within each of these operators uses the appropriate edge1,2,3 data structures


// rhs vector is X
// we implicitly view X as partitioned into components b1 and b2, 
// return vector is Y
// we implicitly view Y as partitioned into components Y1 and Y2, 


//SIMPLE algorithm
//
// 1. Solve F x1=b1
// 2. Solve S x2=-Bx1+b2
// 3a. y1= x1-(1/alpha)*diag(F)^{-1}B'x2
// 3b. y2= (1/alpha)*x2
int simplePreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

#ifdef DEBUG_PRINT_ON
          if (printproc) cout << "In ApplyInverse" << flush<<endl;
#endif
//X  RHS
//Y=Ainv*X Solution

int numv= X.NumVectors();

	double n8; X(0)->Norm2(&n8); 
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Norm of RHS="<<n8<<endl;
#endif


       Y.PutScalar(0.0);

	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug+Belos::TimingDetails+Belos::IterationDetails;
//	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary;
	int verb = Belos::Warnings + Belos::Errors;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	myPL.set( "Convergence Tolerance", 1.0e-4 ); // It would be great to have this tolarance set as a parameter within the xml file
	myPL.set( "Output Frequency", frequency );
	myPL.set( "Timer Label","F Solve"  );


#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Set Preonditioner Parameterlist"<<flush<<endl;
#endif
//tmp rhs malloc
	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);
//CB Is there a more efficient way to zero out this section of the vector?
        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;

#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Initialized B1"<<flush<<endl;
#endif


//Creating RCP to RHS vector 
	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));
        double nfrhs; Fb->Norm2(&nfrhs);
#ifdef DEBUG_PRINT_ON
        if(printproc) cout << "normfrhs="<<nfrhs<<flush<<endl;
#endif

//temp soln malloc
        Epetra_MultiVector x1(Y);
        x1(0)->PutScalar(0.0);

//temp soln malloc
        Epetra_MultiVector y1(Y);
        y1(0)->PutScalar(0.0);

        double sum;
        B1.Norm1(&sum);

#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Set sum"<<flush<<endl;
#endif

        if(sum<1.e-8){//if rhs is zero then don't solve
#ifdef DEBUG_PRINT_ON
          if(printproc)cout<<"rhs sum="<<sum<<" returning 0 solution "<<flush<<endl;
#endif
         }
        else{
#ifdef DEBUG_PRINT_ON
          if(printproc)cout<<"rhs sum="<<sum<<"solving with GMRES"<<flush<<endl;
#endif
         
//temp soln rcp
	  Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          Fx->PutScalar(0.0);
          
	  
	  //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);


          Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );

	  bool Fret = FProblem->setProblem(); 


#ifdef DEBUG_PRINT_ON
          if (printproc) {
            if (Fret == true) {
             cout << "Belos F Linear Problem Set" << flush<<std::endl;
             } 
            else {
             cout << "Error setting Belos F Linear Problem" << std::endl;
            }
          }

#endif


	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

#ifdef DEBUG_PRINT_ON
          if(printproc) cout << "GMRES Solver Manager set"<<flush<<endl;
#endif
          Belos::ReturnType FsolverRet = FSolver.solve();

#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (FsolverRet == Belos::Converged) {
      cout << "Belos F converged."<<flush << std::endl;
    } else {
      cout << "Belos F did not converge." <<flush<< std::endl;
    }
  }
#endif


	  Teuchos::RCP<MV> FSol= FProblem->getLHS();
          x1=*FSol;
          y1=*Fb;
//temp malloc
Epetra_Vector tempx1a(*x1(0));
     double npva; tempx1a.Norm2(&npva);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "fsolnorm="<<npva<<endl;
#endif

//temp malloc
Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "frhsnorm="<<npya<<endl;
#endif

        }



	myPL.set( "Timer Label","Schur Solve"  );


// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

//temp rhs malloc
        Epetra_Vector bx1(*Y(0));
        bx1.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization

        precFunctionblock21(x1(0)->Values(),N, bx1.Values(), precdata);

     double nB21; bx1.Norm2(&nB21);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "normB21="<<nB21<<endl;
#endif

// Then set Sproblos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,x,b) );


// Insert loop here to add -bx1 and b2 and put in Schurb
// b2 the second component of the X vector


//temp malloc
        Epetra_MultiVector b2(X);
//CB more efficient zero out of this section of vector?
        for (int i=0; i<2*(N+1)/3;i++) b2[0][i] = 0.0;


	double nsa; b2(0)->Norm2(&nsa); 
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Norm of RHS Schur A="<<nsa<<endl;
#endif


//    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Norm of RHS Schur B="<<nsb<<flush<<endl;
#endif

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));

#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Schur rhs set="<<flush<<endl;
#endif

//temp rcp
        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Schur solution initialized "<<flush<<endl;
#endif
        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );


#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Schur Problem initialized"<<flush<<endl;
#endif


        bool Sret = SchurProblem->setProblem(); 


#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (Sret == true) {
      cout << "Belos S Linear Problem Set"<<flush<< std::endl;
    } else {
      cout << "Error setting Belos S Linear Problem" <<flush<< std::endl;
    }
  }
#endif

        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Schur GMRES Solver Set"<<flush<<endl;
#endif
        Belos::ReturnType SchursolverRet = SchurSolver.solve();



#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      cout << "Belos Schur converged." <<flush<< std::endl;
    } else {
      cout << "Belos Schur did not converge." << flush<<std::endl;
    }
  }
#endif


	Teuchos::RCP<MV> SchurSol= SchurProblem->getLHS();
        Epetra_MultiVector x2(*SchurSol);



//Next apply dDinvBt to x1 and store in Bx1

        Epetra_Vector dFinvBt(*Y(0));
        dFinvBt.PutScalar(0.0);

	precFunctionblock12(x2(0)->Values(),N, dFinvBt.Values(), precdata);

     double nBt; dFinvBt.Norm2(&nBt);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "normBt="<<nBt<<flush<<endl;
#endif



// Insert loop here to add x1 -(1/alpha)*dFinvBt*x2
// Since the vector is broken up into two components with zeroes padding the 
// 'non-active' entries we can build Y in one action
// We write down that formally this is being done:
// y1= x1 -(1/alpha)*dFinvBt*x2
// y2= (1/alpha)x2 
// 
// which is equivalentlly implemented as Y= x1 -(1/alpha)*dFinvBt*x2 + (1/alpha)x2

Epetra_Vector tempx1(*x1(0));
Epetra_Vector tempx2(*x2(0));
//double alpha=1.0;
double alpha=1.0;
double alphainv=1.0/alpha;


/*
         for (int i=0; i<2*(N+1)/3; i++) {
             tempx1[i]=x1[0][i]-alphainv*dFinvBt[i];
             //tempx1[i]=x1[0][i]-0.0*dFinvBt[i];
          }
         for (int i=2*(N+1)/3;i<N; i++) {
             tempx1[i]=(alphainv*tempx2[i]);
          }
*/
        //Epetra_Vector Ix1(*X(0));
         //for (int i=2*(N+1)/3;i<N; i++) {Ix1[i]=0.0;}
    //for (int i=0; i<N; i++) tempx1[i] = (Ix1[i]+tempx2[i]);
//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]+tempx2[i]);


     double npvf; tempx1.Norm2(&npvf);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "Diagonal components"<<flush<<endl;
if(printproc) cout << "factoredprecnormvec1="<<npvf<<flush<<endl;
#endif
double npvs; tempx2.Norm2(&npvs);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "factoredprecnormvec2="<<npvs<<flush<<endl;
#endif


// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU
//   for (int i=0; i<N; i++) tempx1[i] = tempx1[i]+(alphainv*tempx2[i]);//Diagonal and also DU


//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-dFinvBt[i]) +(tempx2[i]); //No G no alpha


     double npv; tempx1.Norm2(&npv);
#ifdef DEBUG_PRINT_ON
if(printproc) cout << "factoredprecnorm="<<npv<<endl;
#endif
Y=tempx1;


return 0;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
// Add other preconditioners here. Only contructor and 2 methods need to be writen.
// All preconditioners must have xMap_ in constructor's argument list
//    and the following line to build the Base class.
//        : hommePreconditionerBase(xMap_), //Required Base Class construction
