// Trilinos Objects
#include "Piro_Epetra_NOXSolver.hpp"
//#include "Piro_Epetra_LOCASolver.hpp"
#include "trilinosModelEvaluator.hpp"

#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Teuchos_DefaultMpiComm.hpp"


//#define SIMPLE_PREC_ON
#define IDENTITY_PREC_ON

//#define SIMPLE_ML_PREC_ON

//#define PRINT_DEBUG

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

// Objects that are global to the file
static RCP<EpetraExt::ModelEvaluator> Nsolver;
static RCP<trilinosModelEvaluator> model;
static RCP<Teuchos::ParameterList> paramList;
static RCP<Epetra_MpiComm> Comm_;

static EpetraExt::ModelEvaluator::InArgs inArgs;
static EpetraExt::ModelEvaluator::OutArgs outArgs;
static bool printProc;
static int timeStep=1; // time step counter



static Teuchos::RCP<Teuchos::ParameterList> FSolvePL;
static Teuchos::RCP<Teuchos::ParameterList> SchurSolvePL;

static int* FTotalIt;
static int* SchurTotalIt;

static int FTotalIts;
static int SchurTotalIts;

static int OutputStep;



static Teuchos::RCP<Teuchos::ParameterList> HelmSolvePL;
static int* HelmTotalIt;
static int HelmTotalIts;


// Prototypes for function pointers
extern "C" {

/*Defining external functions for finite difference Jacobian*/
  void calc_f(double *, double *, int, void *);
  

  void test_id(double *, int, double *, void *);
  void update_prec_state(double *, int, void *); 

  void test_id(double *, int, double *, void *);
  void update_prec_state(double *, int, void *); 
  void get_jac_vector(double *, int, void *); 

  void sw_picard_block_11(double *, int, double *, void *); 
  void sw_picard_block_21(double *, int, double *, int, void *); 
  void sw_picard_DFinvBt(double *, int, double *, int, void *);
  void sw_picard_schur(double *, int, double *, void *); 

#ifdef SIMPLE_ML_PREC_ON

  void homme_globalIDs(int, int* ,void *);
  void helm_mat(int, int, double *, int *, void *);
  void helm_map(int, int, int *, void *);
  void get_discrete_params(int, int, int, int);

#endif

}

extern "C" {

  // Prototypes for functions in this file.
  RCP<Teuchos::ParameterList> readParameterListAndSetDefaults(MPI_Comm& mpi_comm_c);

  //void printNewtonKrylovFinalStats(Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos>& linsys);
  //void printNewtonKrylovFinalStats();


  void noxinit(int* nelems, double* statevector, int* mpi_comm_f,
               void* blackbox_res, void* blackbox_prec, void* jac_data)
  {

    int N=*nelems;
    void (*residualFunction)(double *, double *, int, void *) = calc_f;



    //void (*jacFunction)(double *,  int, double*,  void *) = sw_picard;
    //void (*jacFunction)(double *,  int, double*,  void *) = sw_jacobian;
    void (*jacFunction)(double *, int, double*,  void *) = test_id;
    
    void (*precFunction)(double *, int, double*,  void *) = test_id;
    
    void (*precUpdateFunction)(double *, int, void *)=update_prec_state;
    void (*getJacVector)(double *, int, void *)=get_jac_vector;

  /*Convention: for SIMPLE algorithm based on Picard linearization
    11 block is F
    12 block is diag(F)^{-1}B'
    21 block is B
    22 block is S=G-Bdiag(F)^{-1}B' */
    

    #ifdef SIMPLE_PREC_ON
    void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    void (*precFunctionblock12)(double *, int, double*,int,  void *) = sw_picard_DFinvBt;
    void (*precFunctionblock21)(double *, int, double*,int,  void *) = sw_picard_block_21;
    void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_schur;
    #endif

    #ifdef SIMPLE_ML_PREC_ON
    void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    void (*precFunctionblock12)(double *, int, double*,int,  void *) = sw_picard_DFinvBt;
    void (*precFunctionblock21)(double *, int, double*,int,  void *) = sw_picard_block_21;
    void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_schur;

    void (*get_globalIDs)(int, int *, void *) = homme_globalIDs;
    void (*get_HelmElementMat)(int, int, double *,int *, void *)=helm_mat;
    void (*get_HelmMap)(int, int, int *, void *)=helm_map;

    #endif


    bool succeeded=true;
    try {

      // Build the epetra communicator
      //  AGS: MPI_Comm_f2c did not work on my machine -- all procs were rank=0
      //       Switching to MPI_COMM_WORLD for now. 
      MPI_Comm mpi_comm_c = MPI_COMM_WORLD; //MPI_Comm_f2c(*mpi_comm_f);
      Comm_=rcp(new Epetra_MpiComm(mpi_comm_c));
      Epetra_Comm& Comm=*Comm_;
      printProc = (Comm_->MyPID() == 0);
    
      if (printProc) cout << "NOXINIT CALLED    for nelem=" << N 
                          << ". Looking for trilinosOptions.xml" << endl;

      try {
        paramList = readParameterListAndSetDefaults(mpi_comm_c);
      }
      catch (std::exception& e) {
        cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"
             << e.what() << "\nExiting: Invalid trilinosOptions.xml file."
             << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
        exit(1);
      }


      if (printProc) cout << "NOXInit: param list is: (delete this debug line)\n" << *paramList << endl;

      //Interface for scalar preconditioner
      /* Interface for IDENTITY  preconditioner */
      #ifdef IDENTITY_PREC_ON
        model = Teuchos::rcp(new trilinosModelEvaluator(N, statevector, Comm, 
                             blackbox_res, blackbox_prec, residualFunction, precUpdateFunction)); 
      #endif


      /* Interface for SIMPLE preconditioner*/
      #ifdef SIMPLE_PREC_ON
      model = Teuchos::rcp(new trilinosModelEvaluator(N, statevector, Comm, 
			      blackbox_res, blackbox_prec,jac_data, 
			      residualFunction, 
			      precFunctionblock11,precFunctionblock12, 
			      precFunctionblock21,precFunctionblock22, 
			      jacFunction,precUpdateFunction,getJacVector,
			      FSolvePL,SchurSolvePL,
			      FTotalIt,SchurTotalIt
			      ));

      #endif


      #ifdef SIMPLE_ML_PREC_ON

      int nets=1;
      int nete=N;
      int np=4;
      int nlev=1;
      get_discrete_params(nets,nete,np,nlev);

      model = Teuchos::rcp(new trilinosModelEvaluator(N, statevector, Comm, 
			      blackbox_res, blackbox_prec,jac_data, 
			      residualFunction, 
			      precFunctionblock11,precFunctionblock12, 
			      precFunctionblock21,precFunctionblock22, 
			      jacFunction,precUpdateFunction,getJacVector,
			      FSolvePL,SchurSolvePL,
			      FTotalIt,SchurTotalIt, 
			      get_globalIDs,get_HelmElementMat,get_HelmMap,
			      HelmSolvePL, HelmTotalIt, nets, nete, np, nlev 
			      ));

			      
      #endif



      Nsolver = rcp(new Piro::Epetra::NOXSolver(paramList, model));

      inArgs=Nsolver->createInArgs();
      outArgs=Nsolver->createOutArgs();

      // Ask the model for the converged solution from g(0)
      RCP<const Epetra_Map> xmap = Nsolver->get_g_map(0);
      RCP<Epetra_Vector> xout = rcp(new Epetra_Vector(*xmap));

      outArgs.set_g(0,xout);

      // Time step counter: just for deciding whether to use continuation on relaxatin param
      timeStep++;

    } //end try block
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
    if (!succeeded) exit(1);
  }

  void noxsolve(int* nelems, double* statevector, void* blackbox_res, void* blackbox_prec, void* jac_data)
  {
    try {
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(Nsolver), logic_error, 
          "Exception: noxsolve called with Nsolver=null: "
          << "either did call noxinit 1st, or called noxfinish already");

      // reset solver with new initial Guess
      model->resetInitialGuess(statevector);

      // reset interface with new blackbox
      //model->resetBlackbox(blackbox_res, blackbox_prec, jac_data);
      model->resetBlackbox(blackbox_res, blackbox_prec);

      // Solve    
      Nsolver->evalModel(inArgs,outArgs);

      // Copy out the solution
      RCP<Epetra_Vector> xout = outArgs.get_g(0); 
      if(xout == Teuchos::null) throw "evalModel is NOT returning a vector";

      for (int i=0; i<*nelems; i++) statevector[i] = (*xout)[i];

    }
    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception!" << endl;
    }
  }

  void noxfinish()
  {
      if (printProc) {
cout<<"FTotalIts="<< FTotalIts<<endl;
cout<<"SchurTotalIts="<< SchurTotalIts<<endl;
      }
    Nsolver   = Teuchos::null;
    model     = Teuchos::null;
    paramList = Teuchos::null;
    Comm_     = Teuchos::null;
  }

  RCP<Teuchos::ParameterList> readParameterListAndSetDefaults(MPI_Comm& mpi_comm_c)
  {
    RCP<Teuchos::ParameterList> 
    params = rcp(new Teuchos::ParameterList("Trilinos Options for Piro"));

    Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(mpi_comm_c));
    Teuchos::updateParametersFromXmlFileAndBroadcast(
                           "trilinosOptions.xml", params.ptr(),tcomm);

    // Set 2 defaults different then NOX defaults unless overwridden in input file
    if (!params->isParameter("Jacobian Operator"))
      params->set("Jacobian Operator","Matrix-Free");
    if (!params->isParameter("Lean Matrix Free"))
      params->set("Lean Matrix Free", true);

    Teuchos::ParameterList& nlPrintParams = params->sublist("NOX").sublist("Printing");
    nlPrintParams.set("MyPID", Comm_->MyPID());
    if (!nlPrintParams.isParameter("Output Information"))
      nlPrintParams.set("Output Information", 67);
    // Following is required for our usage of Piro::NOX.
    params->sublist("NOX").set("Reset Initial Guess",true);

    // Validate top level of parameter list
    Teuchos::ParameterList validPL("Valid List");;
    validPL.sublist("NOX"); validPL.sublist("LOCA");
    validPL.sublist("FSolvePL"); validPL.sublist("SchurSolvePL");
    validPL.set("Jacobian Operator","Matrix-Free");
    validPL.set("Lean Matrix Free", true);
    validPL.set("Matrix-Free Perturbation", 1.0e-6);
    params->validateParameters(validPL, 0);


    FSolvePL = Teuchos::rcp(new Teuchos::ParameterList);
    *FSolvePL = params->sublist("FSolvePL");
    FTotalIts=0;
    FTotalIt=&FTotalIts;


    SchurSolvePL = Teuchos::rcp(new Teuchos::ParameterList);
    *SchurSolvePL = params->sublist("SchurSolvePL");
    SchurTotalIts=0;
    SchurTotalIt=&SchurTotalIts;

    HelmSolvePL = Teuchos::rcp(new Teuchos::ParameterList);
    *HelmSolvePL = params->sublist("HelmSolvePL");
    HelmTotalIts=0;
    HelmTotalIt=&HelmTotalIts;

    return params;
  }

//Andy this is the function that I am trying to write to accumulate the total newton iteration and the total linear iteration

/*
  void printNewtonKrylovFinalStats()
  {
    static int totalNewtonIters=0;
    static int totalKrylovIters=0;
    static int stepNum=0;
    int NewtonIters = nlParams.sublist("Output").get("Nonlinear Iterations", -1000);
    int KrylovIters = linsys->getLinearItersTotal();
      
    totalNewtonIters += NewtonIters;
    stepNum++;
    if (stepNum > 1) {
    //make this 120 parameter something read in from file based on IO frequency
    //if (stepNum == 120) {
      globalData->locaUtils->out() << "Convergence Stats:"<<endl 
                             << "Linear Tolerance " << lintol << endl
                             << "Time Steps " << stepNum << endl
                             << "Total Newton Iterations " << totalNewtonIters <<endl
                             << "Total Krylov Iterations  " << KrylovIters <<endl
              << "Average Newton Steps per timestep" << (double) totalNewtonIters / (double) stepNum<<endl
              << "Average Linear Steps per Newton Step" << (double) KrylovIters / (double) totalNewtonIters <<endl
              << "Average Linear Steps per Time Step" << (double) KrylovIters / (double) stepNum<<endl;
     //}
    }
  }
  */
} //extern "C"
