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
//#include <AztecOO.h>
#include "precon_interface.hpp"
#include "block_precon_interface.hpp"


#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>


#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

//#define DEBUG_PRINT_ON 

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



//SIMPLE preconditioner:  SIMPLE_PREC_CLIP_ON
trilinosModelEvaluator::trilinosModelEvaluator(int nelems, double* statevector,
		const Epetra_Comm& comm_,
		void* blackbox_res_, void* precdata_, void* jacdata_,
		void (*residualFunction_)(double *, double *, int, void *),
		void (*precFunctionblock11_)(double *, int, double*, void *),
		void (*precFunctionblock12_)(double *, int, double*, int, void *),
		void (*precFunctionblock21_)(double *, int, double*, int, void *),
		void (*precFunctionblock22_)(double *, int, double*, void *),
		void (*jacFunction_)(double *, int, double*, void *),
		void (*precUpdateFunction_)(double *, int, void *),
		void (*getJacVector_)(double *, int, void *),
                const RCP<ParameterList>&  FSolvePL_,
                const RCP<ParameterList>&  SchurSolvePL_,
                int* FTotalIt_,
                int* SchurTotalIt_
		):
	N(nelems),
	comm(comm_),
	blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
	residualFunction(residualFunction_),
	jacFunction(jacFunction_),
	precUpdateFunction(precUpdateFunction_),
	getJacVector(getJacVector_)
{ 

  initialize(statevector);

  precOp = Teuchos::rcp(new simpleClipPreconditioner(N, xVec, xMap, blackbox_res,
                                                 precdata, 
                                                 precFunctionblock11_, precFunctionblock12_,
                                                 precFunctionblock21_, precFunctionblock22_,
                                                 precUpdateFunction, FSolvePL_,SchurSolvePL_,FTotalIt_,SchurTotalIt_) );
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
    //double nrm; x->Norm2(&nrm); std::cout << "EvalModel Norm x = " << nrm << std::endl;
    residualFunction(x->Values(), f->Values(), N, blackbox_res);
  }

  if (preconditionerRequested) {
    precOp->computePreconditioner(x,precdata);
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
  if (printproc) std::cout << "Constructing preconditioner:  identityPreconditioner" << std::endl;
}

int identityPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew, void* precdata_)
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

  if (printproc) std::cout << "Constructing preconditioner:  simplePreconditioner" << std::endl;

  F=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock11));
  S=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock22));
}

int simplePreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew, void* precdata_)
{
    precdata = precdata_;
      // Update state in preconditioner code
     precUpdateFunction(xVecNew->Values(), N, precdata);
     const Epetra_Comm& comm = xVec->Comm();
     F=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock11));
     S=Teuchos::rcp ( new Precon_Interface(N,xMap,comm,precdata,precFunctionblock22));




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
          if (printproc) std::cout << "In ApplyInverse" << std::flush<<std::endl;
#endif
//X  RHS
//Y=Ainv*X Solution

int numv= X.NumVectors();

	double n8; X(0)->Norm2(&n8); 
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Norm of RHS="<<n8<<std::endl;
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
if(printproc) std::cout << "Set Preonditioner Parameterlist"<<std::flush<<std::endl;
#endif
//tmp rhs malloc
	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);
//CB Is there a more efficient way to zero out this section of the vector?
        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Initialized B1"<<std::flush<<std::endl;
#endif


//Creating RCP to RHS vector 
	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));
        double nfrhs; Fb->Norm2(&nfrhs);
#ifdef DEBUG_PRINT_ON
        if(printproc) std::cout << "normfrhs="<<nfrhs<<std::flush<<std::endl;
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
if(printproc) std::cout << "Set sum"<<std::flush<<std::endl;
#endif

        if(sum<1.e-8){//if rhs is zero then don't solve
#ifdef DEBUG_PRINT_ON
          if(printproc)std::cout<<"rhs sum="<<sum<<" returning 0 solution "<<std::flush<<std::endl;
#endif
         }
        else{
#ifdef DEBUG_PRINT_ON
          if(printproc)std::cout<<"rhs sum="<<sum<<"solving with GMRES"<<std::flush<<std::endl;
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
             std::cout << "Belos F Linear Problem Set" << std::flush<<std::endl;
             } 
            else {
             std::cout << "Error setting Belos F Linear Problem" << std::endl;
            }
          }

#endif


	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

#ifdef DEBUG_PRINT_ON
          if(printproc) std::cout << "GMRES Solver Manager set"<<std::flush<<std::endl;
#endif
          Belos::ReturnType FsolverRet = FSolver.solve();

#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (FsolverRet == Belos::Converged) {
      std::cout << "Belos F converged."<<std::flush << std::endl;
    } else {
      std::cout << "Belos F did not converge." <<std::flush<< std::endl;
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
if(printproc) std::cout << "fsolnorm="<<npva<<std::endl;
#endif

//temp malloc
Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "frhsnorm="<<npya<<std::endl;
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
if(printproc) std::cout << "normB21="<<nB21<<std::endl;
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
if(printproc) std::cout << "Norm of RHS Schur A="<<nsa<<std::endl;
#endif


//    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Norm of RHS Schur B="<<nsb<<std::flush<<std::endl;
#endif

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur rhs set="<<std::flush<<std::endl;
#endif

//temp rcp
        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur solution initialized "<<std::flush<<std::endl;
#endif
        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );


#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur Problem initialized"<<std::flush<<std::endl;
#endif


        bool Sret = SchurProblem->setProblem(); 


#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (Sret == true) {
      std::cout << "Belos S Linear Problem Set"<<std::flush<< std::endl;
    } else {
      std::cout << "Error setting Belos S Linear Problem" <<std::flush<< std::endl;
    }
  }
#endif

        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur GMRES Solver Set"<<std::flush<<std::endl;
#endif
        Belos::ReturnType SchursolverRet = SchurSolver.solve();



#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      std::cout << "Belos Schur converged." <<std::flush<< std::endl;
    } else {
      std::cout << "Belos Schur did not converge." << std::flush<<std::endl;
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
if(printproc) std::cout << "normBt="<<nBt<<std::flush<<std::endl;
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
if(printproc) std::cout << "Diagonal components"<<std::flush<<std::endl;
if(printproc) std::cout << "factoredprecnormvec1="<<npvf<<std::flush<<std::endl;
#endif
double npvs; tempx2.Norm2(&npvs);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "factoredprecnormvec2="<<npvs<<std::flush<<std::endl;
#endif


// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU
//   for (int i=0; i<N; i++) tempx1[i] = tempx1[i]+(alphainv*tempx2[i]);//Diagonal and also DU


//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-dFinvBt[i]) +(tempx2[i]); //No G no alpha


     double npv; tempx1.Norm2(&npv);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "factoredprecnorm="<<npv<<std::endl;
#endif
Y=tempx1;


return 0;
}






/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
simpleClipPreconditioner::simpleClipPreconditioner (
       int N_, RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_,
       void* blackbox_res_, void* precdata_,
       void (*precFunctionblock11_)(double *, int, double*, void *),
       void (*precFunctionblock12_)(double *,int,double*,int, void *),
       void (*precFunctionblock21_)(double *,int,double*,int, void *),
       void (*precFunctionblock22_)(double *, int, double*, void *),
       void (*precUpdateFunction_)(double *, int, void *),
       const RCP<ParameterList>&  FSolvePL_,
       const RCP<ParameterList>&  SchurSolvePL_,
       int* FTotalIt_,
       int* SchurTotalIt_
       )
       : hommePreconditionerBase(xMap_), //Required Base Class construction
         N(N_), xVec(xVec_), xMap(xMap_),
         blackbox_res(blackbox_res_), precdata(precdata_),
         precFunctionblock11(precFunctionblock11_),
         precFunctionblock12(precFunctionblock12_),
         precFunctionblock21(precFunctionblock21_),
         precFunctionblock22(precFunctionblock22_),
         precUpdateFunction(precUpdateFunction_),
	 FSolvePL(FSolvePL_),
	 SchurSolvePL(SchurSolvePL_),
	 FTotalIt(FTotalIt_),
	 SchurTotalIt(SchurTotalIt_)
{
  const Epetra_Comm& comm = xVec->Comm();
  if (comm.MyPID()==0) printproc=true;
  else                 printproc=false;

  if (printproc) std::cout << "Constructing preconditioner:  simpleClipPreconditioner" << std::endl;


  UVMap = Teuchos::rcp(new Epetra_Map(-1, 2*N/3, 0, comm));
  HMap = Teuchos::rcp(new Epetra_Map(-1, N/3, 0, comm));

  bool zeroout=true;
  workvector4 = Teuchos::rcp(new Epetra_Vector(*HMap,zeroout ));

  dFinvBt = Teuchos::rcp(new Epetra_Vector(*UVMap,zeroout ));

  bx1 = Teuchos::rcp(new Epetra_Vector(*HMap,zeroout ));

  Fb= Teuchos::rcp(new Epetra_Vector(*UVMap,zeroout )); //uv workvector
  Fx= Teuchos::rcp(new Epetra_Vector(*UVMap,zeroout )); //uv workvector

  Schurb= Teuchos::rcp(new Epetra_Vector(*HMap,zeroout )); //h workvector
  Schurx= Teuchos::rcp(new Epetra_Vector(*HMap,zeroout )); //h workvector
 
  F=Teuchos::rcp ( new Precon_Interface(2*N/3,UVMap,comm,precdata,precFunctionblock11_));
  S=Teuchos::rcp ( new Precon_Interface(N/3,HMap,comm,precdata,precFunctionblock22_));



#ifdef DEBUG_PRINT_ON      
       if(printproc) std::cout<<"F Solver Precon Parameters"<<std::endl;
       if(printproc) std::cout<<"Block Size "<<FSolvePL->get<int>("Block Size")<<std::endl;
       if(printproc) std::cout<<"Num Blocks "<<FSolvePL->get<int>("Num Blocks")<<std::endl;
       if(printproc) std::cout<<"Maximum Iterations "<<FSolvePL->get<int>("Maximum Iterations")<<std::endl;
       if(printproc) std::cout<<"Convergence Tolerance "<<FSolvePL->get<double>("Convergence Tolerance")<<std::endl;
       if(printproc) std::cout<<"Output Frequency "<<FSolvePL->get<int>("Output Frequency")<<std::endl;

       if(printproc) std::cout<<"Schur Solver Precon Parameters"<<std::endl;
       if(printproc) std::cout<<"Block Size "<<SchurSolvePL->get<int>("Block Size")<<std::endl;
       if(printproc) std::cout<<"Num Blocks "<<SchurSolvePL->get<int>("Num Blocks")<<std::endl;
       if(printproc) std::cout<<"Maximum Iterations "<<SchurSolvePL->get<int>("Maximum Iterations")<<std::endl;
       if(printproc) std::cout<<"Convergence Tolerance "<<SchurSolvePL->get<double>("Convergence Tolerance")<<std::endl;
       if(printproc) std::cout<<"Output Frequency "<<SchurSolvePL->get<int>("Output Frequency")<<std::endl;
#endif




  FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );
  FSolver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>( FProblem, FSolvePL ) );
  SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );
  SchurSolver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>( SchurProblem, SchurSolvePL ) );



}

int simpleClipPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew, void* precdata_)
{
  precdata = precdata_;
  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);
  const Epetra_Comm& comm = xVec->Comm();


  F=Teuchos::rcp ( new Precon_Interface(2*N/3,UVMap,comm,precdata,precFunctionblock11));
  S=Teuchos::rcp ( new Precon_Interface(N/3,HMap,comm,precdata,precFunctionblock22));

  return 0;
}


int simpleClipPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

// const Epetra_Vector & x_v = dynamic_cast<const Epetra_Vector&> (X);

#ifdef DEBUG_PRINT_ON
          if (printproc) std::cout << "In ApplyInverse" << std::flush<<std::endl;
#endif
//X  RHS
//Y=Ainv*X Solution

int numv= X.NumVectors();

	double n8; X(0)->Norm2(&n8); 
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Norm of RHS="<<n8<<std::endl;
#endif

       Y.PutScalar(0.0);

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Set Preconditioner Parameterlist"<<std::flush<<std::endl;
#endif



       Fb->PutScalar(0.0);


        for (int i=0;i<2*(N+1)/3; i++) (*Fb)[i] = X[0][i];


        double nfrhs; Fb->Norm2(&nfrhs);
#ifdef DEBUG_PRINT_ON
        if(printproc) std::cout << "normfrhs="<<nfrhs<<std::flush<<std::endl;
#endif

        double sum; Fb->Norm2(&sum);

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Set sum"<<std::flush<<std::endl;
#endif

        if(sum<1.e-8){//if rhs is zero then don't solve
#ifdef DEBUG_PRINT_ON
          if(printproc)std::cout<<"rhs sum="<<sum<<" returning 0 solution "<<std::flush<<std::endl;
#endif
         }
        else{
#ifdef DEBUG_PRINT_ON
          if(printproc)std::cout<<"rhs sum="<<sum<<"solving with GMRES"<<std::flush<<std::endl;
#endif
         
//temp soln rcp
	 // Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          Fx->PutScalar(0.0);
          
	  
	  //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);

// FProblem->reset( F,Fx,Fb );
	  FProblem->setOperator( F);
	  FProblem->setLHS( Fx);
	  FProblem->setRHS( Fb);


	  FSolver->reset( Belos::Problem );

          Belos::ReturnType FSolverRet = FSolver->solve();

#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (FSolverRet == Belos::Converged) {
      std::cout << "Belos F converged."<<std::flush << std::endl;
    } else {
      std::cout << "Belos F did not converge." <<std::flush<< std::endl;
    }
  }
#endif

	*FTotalIt+=FSolver->getNumIters();
     double npva; Fx->Norm2(&npva);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "fsolnorm="<<npva<<std::endl;
#endif

        }//end FSolve


// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

//temp rhs malloc
        bx1->PutScalar(0.0);

        precFunctionblock21((*Fx).Values(),2*(N+1)/3, (*bx1).Values(),(N+1)/3, precdata);

     double nB21; bx1->Norm2(&nB21);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "normB21="<<nB21<<std::endl;
#endif

        workvector4->PutScalar(0.0);
        
	for (int i=2*(N+1)/3;i<N;i++) (*workvector4)[i-2*(N+1)/3] = X[0][i];

	double nsa; workvector4->Norm2(&nsa); 
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Norm of RHS Schur A="<<nsa<<std::endl;
#endif
        Schurb->PutScalar(0.0);
        Schurx->PutScalar(0.0);

   //for (int i=2*(N+1)/3; i<N; i++) (*Schurb)[i] = (*workvector4)[i]-(*bx1)[i]; // Lower block

   for (int i=0; i<(N+1)/3; i++) (*Schurb)[i] = (*workvector4)[i]-(*bx1)[i]; // Lower block


	double nsb; Schurb->Norm2(&nsb); 
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Norm of RHS Schur B="<<nsb<<std::flush<<std::endl;
#endif


#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur rhs set="<<std::flush<<std::endl;
#endif


#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur solution initialized "<<std::flush<<std::endl;
#endif

	  SchurProblem->setOperator( S);
	  SchurProblem->setLHS( Schurx);
	  SchurProblem->setRHS( Schurb);

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur Problem initialized"<<std::flush<<std::endl;
#endif


	bool Sret = SchurProblem->setProblem(); 


#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (Sret == true) {
      std::cout << "Belos S Linear Problem Set"<<std::flush<< std::endl;
    } else {
      std::cout << "Error setting Belos S Linear Problem" <<std::flush<< std::endl;
    }
  }
#endif

	

	  SchurSolver->reset( Belos::Problem );
	  //SchurSolver->setProblem( SchurProblem );


#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur GMRES Solver Set"<<std::flush<<std::endl;
#endif
        Belos::ReturnType SchursolverRet = SchurSolver->solve();

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "Schur GMRES Solved "<<std::flush<<std::endl;
#endif

	*SchurTotalIt+=SchurSolver->getNumIters();

#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "SchurTotalIt="<<*SchurTotalIt<<std::flush<<std::endl;
#endif

#ifdef DEBUG_PRINT_ON
if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      std::cout << "Belos Schur converged." <<std::flush<< std::endl;
    } else {
      std::cout << "Belos Schur did not converge." << std::flush<<std::endl;
    }
  }
#endif


//Next apply dDinvBt to x1 and store in Bx1

     //   Epetra_Vector dFinvBt(*Y(0));
        dFinvBt->PutScalar(0.0);

	precFunctionblock12((*Schurx).Values(),(N+1)/3, (*dFinvBt).Values(),2*(N+1)/3, precdata);

     double nBt; dFinvBt->Norm2(&nBt);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "normBt="<<nBt<<std::flush<<std::endl;
#endif




double alpha=1.0;
double alphainv=1.0/alpha;


// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
  //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*(*dFinvBt)[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU
//   for (int i=0; i<N; i++) tempx1[i] = tempx1[i]+(alphainv*tempx2[i]);//Diagonal and also DU

    for (int i=0; i<2*(N+1)/3; i++) Y[0][i] = ((*Fx)[i]-alphainv*(*dFinvBt)[i]);
	    
    for (int i=2*(N+1)/3; i<N; i++) Y[0][i] = (alphainv*(*Schurx)[i-2*(N+1)/3]); 
	    



     double npv;Y.Norm2(&npv);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "factoredprecnorm="<<npv<<std::endl;
#endif

return 0;
}











/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
// Add other preconditioners here. Only contructor and 2 methods need to be writen.
// All preconditioners must have xMap_ in constructor's argument list
//    and the following line to build the Base class.
//        : hommePreconditionerBase(xMap_), //Required Base Class construction
