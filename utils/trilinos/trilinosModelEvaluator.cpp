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
#include "precon_interface.hpp"


#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>


#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"


#include "Epetra_InvOperator.h"
#include "Epetra_Util.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"


#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"


#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace std;


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


//SIMPLE preconditioner:  SIMPLE_PREC_ON
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

  precOp = Teuchos::rcp(new simplePreconditioner(N, xVec, xMap, blackbox_res,
                                                 precdata, 
                                                 precFunctionblock11_, precFunctionblock12_,
                                                 precFunctionblock21_, precFunctionblock22_,
                                                 precUpdateFunction, FSolvePL_,SchurSolvePL_,FTotalIt_,SchurTotalIt_) );
}




//SIMPLE_ML preconditioner:  SIMPLE_ML_PREC_ON
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
                int* SchurTotalIt_,
		void (*homme_globalIDs_)(int, int* ,void *),
		void (*helm_mat_)(int, int, double *, int *, void *),
		void (*helm_map_)(int, int, int *, void *),
                const RCP<ParameterList>&  HelmSolvePL_,
                int* HelmTotalIt_,
		int nets_,
		int nete_,
		int np_,
		int nlev_
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
  

  precOp = Teuchos::rcp(new simpleMLPreconditioner(N, xVec, xMap, blackbox_res,
                                                 precdata, 
                                                 precFunctionblock11_, precFunctionblock12_,
                                                 precFunctionblock21_, precFunctionblock22_,
                                                 precUpdateFunction, 
						 FSolvePL_,SchurSolvePL_,
						 FTotalIt_,SchurTotalIt_,
						 homme_globalIDs_,helm_mat_,helm_map_,
						 HelmSolvePL_,
						 HelmTotalIt_,
						 nets_,nete_,np_,nlev_ 
						 ) );
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

  if (printproc) std::cout << "Constructing preconditioner:  simplePreconditioner" << std::endl;


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







int simplePreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew, void* precdata_)
{
  precdata = precdata_;
  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);
  const Epetra_Comm& comm = xVec->Comm();


  F=Teuchos::rcp ( new Precon_Interface(2*N/3,UVMap,comm,precdata,precFunctionblock11));
  S=Teuchos::rcp ( new Precon_Interface(N/3,HMap,comm,precdata,precFunctionblock22));

  return 0;
}


int simplePreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
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


   for (int i=0; i<(N+1)/3; i++) (*Schurb)[i] = (*workvector4)[i]-(*bx1)[i]; // Schur RHS


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

        dFinvBt->PutScalar(0.0);

	precFunctionblock12((*Schurx).Values(),(N+1)/3, (*dFinvBt).Values(),2*(N+1)/3, precdata);

     double nBt; dFinvBt->Norm2(&nBt);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "normBt="<<nBt<<std::flush<<std::endl;
#endif




    double alpha=1.0;
    double alphainv=1.0/alpha;

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
simpleMLPreconditioner::simpleMLPreconditioner (
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
       int* SchurTotalIt_,
       void (*homme_globalIDs_)(int, int* ,void *),
       void (*helm_mat_)(int, int, double *, int *, void *),
       void (*helm_map_)(int, int, int *, void *),
       const RCP<ParameterList>&  HelmSolvePL_,
       int* HelmTotalIt_,
       int nets_, int nete_, int np_, int nlev_
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
	 SchurTotalIt(SchurTotalIt_),
	 homme_globalIDs(homme_globalIDs_),
	 helm_mat(helm_mat_),
	 helm_map(helm_map_),
	 HelmSolvePL(HelmSolvePL_),
	 HelmTotalIt(HelmTotalIt_),
	 nets(nets_), nete(nete_), np(np_),nlev(nlev_)
{
  const Epetra_Comm& comm = xVec->Comm();
  if (comm.MyPID()==0) printproc=true;
  else                 printproc=false;

  if (printproc) std::cout << "Constructing preconditioner:  simpleMLPreconditioner" << std::endl;
 
 


  initializeMatrixMaps(comm);

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


// Here is where we need to call the routines to build the maps for matrix assembly, build the helmholtz operator and associated rhs vector
// See the function templates below based on the belos_helm.cpp methods
// simpleMLPreconditioner::initializeMatrixMaps
// simpleMLPreconditioner::BuildMatrix
// simpleMLPreconditioner::BuildRHS
//
//
// After initially setting the problem and right hand side here, these will need to be called again to be updated appropriately within the simpleMLPreconditioner::ApplyInverse method (similar to what is currently done for the Schur solve)
// simpleMLPreconditioner::SetProblem
// simpleMLPreconditioner::HelmholtzSolve( double* statevector){

}







int simpleMLPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVecNew, void* precdata_)
{
  precdata = precdata_;
  // Update state in preconditioner code
  precUpdateFunction(xVecNew->Values(), N, precdata);
  const Epetra_Comm& comm = xVec->Comm();
// This is where we will need to update the preconditioner based on the current state of the Jacobian

  F=Teuchos::rcp ( new Precon_Interface(2*N/3,UVMap,comm,precdata,precFunctionblock11));
  S=Teuchos::rcp ( new Precon_Interface(N/3,HMap,comm,precdata,precFunctionblock22));

  return 0;
}


int simpleMLPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
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

        dFinvBt->PutScalar(0.0);

	precFunctionblock12((*Schurx).Values(),(N+1)/3, (*dFinvBt).Values(),2*(N+1)/3, precdata);

     double nBt; dFinvBt->Norm2(&nBt);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "normBt="<<nBt<<std::flush<<std::endl;
#endif




    double alpha=1.0;
    double alphainv=1.0/alpha;

    for (int i=0; i<2*(N+1)/3; i++) Y[0][i] = ((*Fx)[i]-alphainv*(*dFinvBt)[i]);
	    
    for (int i=2*(N+1)/3; i<N; i++) Y[0][i] = (alphainv*(*Schurx)[i-2*(N+1)/3]); 
	    



     double npv;Y.Norm2(&npv);
#ifdef DEBUG_PRINT_ON
if(printproc) std::cout << "factoredprecnorm="<<npv<<std::endl;
#endif

return 0;
}



// Utility function called by constructors above
void simpleMLPreconditioner::initializeMatrixMaps(const Epetra_Comm& comm)
{
  bool succeeded=true;
  void (*get_globalIDs)(int, int *, void *) = homme_globalIDs;

  try {
    
    my_GlobalIDs=new(int[N]);


    int globalids[N];
    get_globalIDs(N,my_GlobalIDs,precdata);

    for(int i=0;i<N;i++){
            globalids[i]=my_GlobalIDs[i];
    }

    //sorting ids
   
     std::sort(&my_GlobalIDs[0],&my_GlobalIDs[0]+ (N));

     int *end=unique(my_GlobalIDs,&my_GlobalIDs[0]+(N));

     NPPP=end-(&my_GlobalIDs[0]);

     bin=new int*[NPPP];
     for (int i=0;i<N;i++){
          bin[i]=new int[5];
     }

     for(int i=0;i<NPPP;i++){
         for (int j=0;j<5;j++){
                 bin[i][j]=N+1;//initialize to out index out of bounds
          }
     }


    for(int i=0;i<NPPP;i++){
         bin[i][0]=my_GlobalIDs[i];
    }

    int avail;
    bincount=new int[NPPP];

    for(int i=0;i<NPPP;i++){
       bincount[i]=0;
       avail=1;
       for(int j=0;j<N;j++){
           if(globalids[j]==bin[i][0]){
              bin[i][avail]=j;
              avail++;
              bincount[i]++;
           }
       }
    }

    MatrixMap = Teuchos::rcp(new Epetra_Map(-1, NPPP, my_GlobalIDs, 0, comm));

    bool high_rank_proc_owns_shared=false;

    ParMatrixMap = Teuchos::rcp(new Epetra_Map( Epetra_Util::Create_OneToOne_Map(*MatrixMap,high_rank_proc_owns_shared) ));

    exporter = Teuchos::rcp(new Epetra_Export(*MatrixMap, *ParMatrixMap));
    importer = Teuchos::rcp(new Epetra_Import(*MatrixMap, *ParMatrixMap));


  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);
}





void simpleMLPreconditioner::BuildMatrix( void* data ){

    void (*get_HelmElementMat)(int, int, double *,int *, void *)=helm_mat;

     try{

         int matsize=(np)*(np);
         double Aloc[matsize][matsize];
         double Zeroloc[matsize][matsize];
         int Index[matsize];
         int ColumnMajor=3;
         int ndx=my_GlobalIDs[NPPP-1];

         Epetra_DataAccess copy = ::Copy;
         const std::size_t approximate_indices_per_row = N;

         HelmMatrix = Teuchos::rcp(new Epetra_FECrsMatrix(copy, *MatrixMap,approximate_indices_per_row,StaticProfile));
         ParHelmMatrix = Teuchos::rcp(new Epetra_FECrsMatrix(copy, *ParMatrixMap,approximate_indices_per_row,StaticProfile));


         for(int i=0;i<matsize;i++){
             for(int j=0;j<matsize;j++){
                 Zeroloc[i][j]=0.0;
             }
         }

//Initialize Matrix to Zero
     for (int ie=nets;ie<nete+1;ie++){
           get_HelmElementMat(ie,matsize,&Aloc[0][0],&Index[0],data);
           HelmMatrix->InsertGlobalValues(matsize,&Index[0],matsize,&Index[0],&Zeroloc[0][0],ColumnMajor);
     }

//Sum Matrix entries Element by Element 
     for (int ie=nets;ie<nete+1;ie++){
          get_HelmElementMat(ie,matsize,&Aloc[0][0],&Index[0],data);
          HelmMatrix->SumIntoGlobalValues(matsize,&Index[0],matsize,&Index[0],&Aloc[0][0],ColumnMajor);
     }
    
     HelmMatrix->FillComplete();
     ParHelmMatrix->Export(*HelmMatrix, *exporter, Add);
     bool callFillComplete=true;
     ParHelmMatrix->GlobalAssemble(callFillComplete);

  } //end try

  catch (std::exception& e) {
  cout << e.what() << endl;
  }
  catch (const char *s) {
  cout << s << endl;
  }
  catch (...) {
  cout << "Caught unknown exception in BuildMatrix " << endl;
  }
}//end BuildMatrix




void simpleMLPreconditioner::BuildRHS(double* rhsVector, void* data ){
     void (*get_HelmMap)(int, int, int *, void *)=helm_map;

     try{
      int matsize=(np)*(np);
      double rhsLoc[matsize];
      int rhsIndex[matsize];
      int ndx=my_GlobalIDs[N-1];
      double rhsGlob[ndx];
      int rhsGIndex[ndx];

      for(int i=0;i<matsize;i++){
          rhsIndex[i]=0;
          rhsLoc[i]=0.0;
      }

      for(int i=0;i<ndx;i++){
	  rhsGlob[i]=0;
	  rhsGIndex[i]=0;
      }

      for (int ie=nets;ie<nete+1;ie++){
           get_HelmMap(ie,matsize,&rhsIndex[0],data);
	   //  for (int lev=0;lev<1;lev++){
	          int lndx=0;
		  for (int i=0;i<np;i++){
	               for (int j=0;j<np;j++){
		            rhsLoc[i*(np)+j]=rhsVector[ (ie-nets)*(np)*(np) +i*(np)+j] ;
		            rhsGIndex[rhsIndex[lndx]-1]=rhsIndex[lndx];
		            rhsGlob[rhsIndex[lndx]-1]+=rhsLoc[i*(np)+j];
		            lndx++;
		        }

                   }
            //}
      }

      HelmRHS = Teuchos::rcp(new Epetra_Vector(*MatrixMap,StaticProfile));
      ParHelmRHS = Teuchos::rcp(new Epetra_Vector(*ParMatrixMap,StaticProfile));
      HelmRHS->PutScalar(0.0);
      ParHelmRHS->PutScalar(0.0);

      HelmRHS->SumIntoGlobalValues(ndx,&rhsGlob[0],&rhsGIndex[0]);

      for (int i=0;i<N;i++){
           (*HelmRHS)[i]=(*HelmRHS)[i]/bincount[i];
      }

      ParHelmRHS->Export(*HelmRHS, *exporter,Insert);

     } //end try

     catch (std::exception& e) {
       cout << e.what() << endl;
     }
     catch (const char *s) {
       cout << s << endl;
     }
     catch (...) {
       cout << "Caught unknown exception in BuildRHS " << endl;
     }
}//end BuildRHS

	   

void simpleMLPreconditioner::SetProblem(){

   try{
    HelmSolution = Teuchos::rcp(new Epetra_Vector(*MatrixMap,StaticProfile ));
    ParHelmSolution = Teuchos::rcp(new Epetra_Vector(*ParMatrixMap,StaticProfile));

    HelmProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(ParHelmMatrix,ParHelmSolution,ParHelmRHS) );
    HelmProblem->setOperator( ParHelmMatrix);
    HelmProblem->setLHS( ParHelmSolution);
    HelmProblem->setRHS( ParHelmRHS);



    // Begin Preconditioner Set
    
    MLList = Teuchos::rcp(new Teuchos::ParameterList);
    ML_Epetra::SetDefaults("SA",*MLList);
    //ML_Epetra::SetDefaults("DD",*MLList);
    //ML_Epetra::SetDefaults("DD-ML",*MLList);

    MLPrec=Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*ParHelmMatrix,*MLList,true));
    bool Helmret = HelmProblem->setProblem();
    if (Helmret != true && printproc) cout<<"Error setting problem"<<flush<<endl;

    const bool CheckFiltering = false;
    int comp=MLPrec->ComputePreconditioner();

    MLPI=Teuchos::rcp(new Epetra_InvOperator(MLPrec.get()));

    HelmProblem->setRightPrec(MLPI);

    comp=MLPrec->IsPreconditionerComputed();

    const int cycles=1;
    comp= MLPrec->AnalyzeCycle(cycles);

    paramList = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::updateParametersFromXmlFile("helmbelos.xml", paramList.ptr());
    HelmSolvePL = Teuchos::rcp(new Teuchos::ParameterList);
    *HelmSolvePL = paramList->sublist("HelmSolvePL");
    HelmTotalIts=0;
    HelmTotalIt=&HelmTotalIts;

    HelmSolver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>( HelmProblem, HelmSolvePL ) );


   } //end try

   catch (std::exception& e) {
     cout << e.what() << endl;
   }
   catch (const char *s) {
     cout << s << endl;
   }
   catch (...) {
     cout << "Caught unknown exception in SetProblem  " << endl;
   }
}//end SetProblem


void simpleMLPreconditioner::HelmholtzSolve( double* statevector){

   try{

    HelmSolution->PutScalar(0.0);
    ParHelmSolution->PutScalar(0.0);

    HelmProblem->setOperator( ParHelmMatrix);
    HelmProblem->setLHS( ParHelmSolution);
    HelmProblem->setRHS( ParHelmRHS);
   
    HelmSolver->reset( Belos::Problem );

    Belos::ReturnType HelmSolveStatus = HelmSolver->solve();

    HelmTotalIts = HelmSolver->getNumIters();
    HelmTotalIt=&HelmTotalIts;

    HelmSolution->Import(*ParHelmSolution, *importer, Insert);

    if(printproc)cout<<"Trilinos Computed Helmholtz Solve in "<<*HelmTotalIt<<" iterations"<<flush<<endl;

    int ndx=my_GlobalIDs[NPPP-1];
    double solGlob[ndx];

    for (int j=1;j<5;j++){
         for(int i=0;i<NPPP;i++){
	     if(bin[i][j]!=N+1){
	        solGlob[bin[i][j]]=(*HelmSolution)[i];
	     }
	 }
    }

    for(int i=0;i<N;i++){
        statevector[i]=solGlob[i];
    }

    HelmRHS=Teuchos::null;
    ParHelmRHS=Teuchos::null;


   } //end try


   catch (std::exception& e) {
       cout << e.what() << endl;
   }
   catch (const char *s) {
       cout << s << endl;
   }
   catch (...) {
       cout << "Caught unknown exception in HelmholtzSolve " << endl;
   }
}//end HelmholtzSolve



void simpleMLPreconditioner::belosfinish(){
   delete bincount;
   delete bin;
   delete my_GlobalIDs;
   //comm=Teuchos::null;
   MatrixMap=Teuchos::null;
   ParMatrixMap=Teuchos::null;
   exporter=Teuchos::null;
   importer=Teuchos::null;
   HelmMatrix = Teuchos::null;
   ParHelmMatrix = Teuchos::null;
   HelmRHS = Teuchos::null;
   ParHelmRHS = Teuchos::null;
   HelmSolution = Teuchos::null;
   ParHelmSolution = Teuchos::null;
   HelmProblem = Teuchos::null;
   MLList = Teuchos::null;
   MLPrec=Teuchos::null;
   MLPI=Teuchos::null;
   paramList = Teuchos::null;
   HelmSolvePL = Teuchos::null;
   HelmSolver = Teuchos::null;
}







/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
// Add other preconditioners here. Only contructor and 2 methods need to be writen.
// All preconditioners must have xMap_ in constructor's argument list
//    and the following line to build the Base class.
//        : hommePreconditionerBase(xMap_), //Required Base Class construction





