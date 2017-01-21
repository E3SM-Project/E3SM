#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

// Preconditioner Options
#include "trilinosDefines.h"

// Trilinos header files
#include "Piro_Epetra_NOXSolver.hpp"

#include "trilinosModelEvaluator.hpp"

#include "Epetra_MpiComm.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TimeMonitor.hpp"

using std::cout;
using std::endl;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::TimeMonitor;

static RCP<EpetraExt::ModelEvaluator> Nsolver;
static RCP<trilinosModelEvaluator> model;
static RCP<Epetra_MpiComm> Comm_;

static EpetraExt::ModelEvaluator::InArgs inArgs;
static EpetraExt::ModelEvaluator::OutArgs outArgs;

static bool printproc;

static int timeStep=0; // time step counter

// Trilinos XML parameter list
static RCP<ParameterList> paramList;

// Parameter sub-lists in XML file for preconditioner solves
static RCP<ParameterList> ASolvePL;
static RCP<ParameterList> SSolvePL;
static RCP<ParameterList> PSolvePL;

// Iteration counters for preconditioner solves
static int* ATotalIt;
static int* STotalIt;
static int* PTotalIt;

static int ATotalIts;
static int STotalIts;
static int PTotalIts;


// ------------------------------------------------------------------------------
// Prototypes for Fortran functions
// ------------------------------------------------------------------------------
extern "C" {
  // residual function
  void calc_f(double*, double*, int, void*);

  // swim analytic Jacobian vector product
  void sw_jacobian(double*, int, double*, void*);

  // update Jacobian (when using analytic Jacobian)
  void update_jac_state(double*, int, void*);

  // update preconditioner
  void update_prec_state(double*, int, void*); 

  // output analytic dense, analytic sparse, and finite difference Jacobians
  // to disk for comparison
  void check_jacobian(double*, int, void*); 
}

extern "C" {
  // Prototypes for functions in this file
  RCP<ParameterList> readParameterListAndSetDefaults(MPI_Comm& mpi_comm_c);



  // ----------------------------------------------------------------------------
  // Initialize NOX solver
  // ----------------------------------------------------------------------------
  void noxinit(int* np_, int* nlev_, int* nelemd_, 
	       double* StateVector, void* StateData, int* mpi_comm_f)
  {
    // HOMME descretization parameters
    const int np     = *np_;
    const int nlev   = *nlev_;
    const int nelemd = *nelemd_;

    // local state vector length
#ifdef _PRIM
    int nState = (3*np*np*nlev + np*np)*nelemd;
#else //_SWIM
    int nState = 3*np*np*nlev*nelemd;
#endif

    // set residual function
    void (*residualFunction)(double*, double*, int, void*) = calc_f;

    // set Jacobian matvec and update functions
    void (*jacFunction)(double*, int, double*, void*) = NULL; 
    void (*jacUpdateFunction)(double*, int, void*)    = NULL;

    // set preconditioner update function
#ifdef CHECK_JAC
    void (*precUpdateFunction)(double*, int, void*) = check_jacobian;
#else
    void (*precUpdateFunction)(double*, int, void*) = update_prec_state;
#endif

    bool succeeded = true;
    try {
      // Build the epetra communicator
      //  AGS: MPI_Comm_f2c did not work on my machine -- all procs were rank=0
      //       Switching to MPI_COMM_WORLD for now. 
      MPI_Comm mpi_comm_c = MPI_COMM_WORLD; //MPI_Comm_f2c(*mpi_comm_f);
      Comm_=rcp(new Epetra_MpiComm(mpi_comm_c));
      Epetra_Comm& Comm=*Comm_;
      printproc = (Comm_->MyPID() == 0);
      
      if (printproc) cout << "NOXINIT CALLED    for nelem=" << nState 
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

      if (printproc) cout << "NOXInit: Parameter List: \n" << *paramList << endl;

      model = rcp(new trilinosModelEvaluator(np, nlev, nelemd, 
 					     StateVector, StateData, Comm, 
					     residualFunction, 
					     jacFunction, jacUpdateFunction,
					     precUpdateFunction,
					     ASolvePL, SSolvePL, PSolvePL,
					     ATotalIt, STotalIt, PTotalIt));
      
      Nsolver = rcp(new Piro::Epetra::NOXSolver(paramList, model));
      
      inArgs  = Nsolver->createInArgs();
      outArgs = Nsolver->createOutArgs();

      // Ask the model for the converged solution from g(0)
      RCP<const Epetra_Map> xmap = Nsolver->get_g_map(0);
      RCP<Epetra_Vector>    xout = rcp(new Epetra_Vector(*xmap));

      outArgs.set_g(0,xout);

    } //end try block
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
    if (!succeeded) exit(1);
  }
  

  // ----------------------------------------------------------------------------
  // Solve nonlinear system with NOX
  // ----------------------------------------------------------------------------
  void noxsolve(int* nelems, double* StateVector, void* StateData, int* ierr)
  {
    *ierr = 0; // solver success flag

    try {
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(Nsolver), std::logic_error, 
				 "Exception: noxsolve called with Nsolver=null: "
				 << "either did not call noxinit 1st, or called noxfinish already");

      // reset solver with new initial Guess
      model->resetInitialGuess(StateVector);

      // reset interface with new state data
      model->resetStateData(StateData);

      // solve    
      Nsolver->evalModel(inArgs, outArgs);

      // output timers
      // Teuchos::oblackholestream blackhole;
      // std::ostream &out = (printproc ? std::cout : blackhole);
      // TimeMonitor::summarize(out);

      // check for solver failure
      if (outArgs.isFailed()) {
	*ierr = 1;
	throw "Noxsolve Error: Newton Failed to Converge";    
      }

      // Copy out the solution
      RCP<Epetra_Vector> xout = outArgs.get_g(0); 
      if(xout == Teuchos::null) throw "evalModel is NOT returning a vector";

      for (int i=0; i<*nelems; i++) StateVector[i] = (*xout)[i];

      // Time step counter: just for deciding whether to use continuation on relaxatin param
      timeStep++;
    }
    catch (std::exception& e) {
      *ierr = -1;
      cout << e.what() << endl;
    }
    catch (const char *s) {
      *ierr = -1;
      cout << s << endl;
    }
    catch (...) {
      *ierr = -1;
      cout << "Caught unknown exception in noxsolve!" << endl;
    }
  }
  
  
  // ----------------------------------------------------------------------------
  // print solver data and clean up at end of run
  // ----------------------------------------------------------------------------
  void noxfinish()
  {
    // print timings
    Teuchos::oblackholestream blackhole;
    std::ostream &out = (printproc ? std::cout : blackhole);
    
    TimeMonitor::summarize(out);

    if (printproc) {
      cout << "TimeSteps= " << timeStep  << endl;
      cout << "ATotalIts= " << ATotalIts << endl;
      cout << "STotalIts= " << STotalIts << endl;
      cout << "PTotalIts= " << PTotalIts << endl;
    }
    Nsolver   = Teuchos::null;
    model     = Teuchos::null;
    paramList = Teuchos::null;
    Comm_     = Teuchos::null;
  }


  // ----------------------------------------------------------------------------
  // read solver xml parameter file
  // ----------------------------------------------------------------------------
  RCP<ParameterList> readParameterListAndSetDefaults(MPI_Comm& mpi_comm_c)
  {
    RCP<ParameterList> params = rcp(new ParameterList("Trilinos Options for Piro"));

    Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(mpi_comm_c));
    Teuchos::updateParametersFromXmlFileAndBroadcast("trilinosOptions.xml", 
						     params.ptr(), tcomm);

    // Set 2 defaults different then NOX defaults unless overwridden in input file
    if (!params->isParameter("Jacobian Operator"))
      params->set("Jacobian Operator","Matrix-Free");

    if (!params->isParameter("Lean Matrix Free"))
      params->set("Lean Matrix Free", true);

    ParameterList& nlPrintParams = params->sublist("NOX").sublist("Printing");
    nlPrintParams.set("MyPID", Comm_->MyPID());
    if (!nlPrintParams.isParameter("Output Information"))
      nlPrintParams.set("Output Information", 67);

    // Following is required for our usage of Piro::NOX.
    params->sublist("NOX").set("Reset Initial Guess", true);

    // Validate top level of parameter list (error catching)
    ParameterList validPL("Valid List");
    validPL.sublist("NOX"); 
    validPL.sublist("LOCA");
    validPL.sublist("ASolvePL"); 
    validPL.sublist("SSolvePL");
    validPL.sublist("PSolvePL");
    validPL.set("Jacobian Operator", "Matrix-Free");
    validPL.set("Lean Matrix Free", true);
    validPL.set("Matrix-Free Perturbation", 1.0e-6);
    params->validateParameters(validPL, 0);

    // Parameter list for A solve in preconditioner
    ASolvePL  = rcp(new ParameterList);
    *ASolvePL = params->sublist("ASolvePL");
    ATotalIts = 0;
    ATotalIt  = &ATotalIts;

    // Parameter list for S solve in preconditioner
    SSolvePL  = rcp(new ParameterList);
    *SSolvePL = params->sublist("SSolvePL");
    STotalIts = 0;
    STotalIt  = &STotalIts;

    // Parameter list for P solve in preconditioner
    PSolvePL  = rcp(new ParameterList);
    *PSolvePL = params->sublist("SSolvePL");
    PTotalIts = 0;
    PTotalIt  = &PTotalIts;

    return params;
  }


  // ----------------------------------------------------------------------------
  // Preconditioner type passed to Fortran:
  //  -1 = Check Jacobian test
  //   0 = Identity
  //   1 = PC with Spectral Element Operators
  //   2 = PC with Fotran Dense Matrices
  //   3 = PC with Trilinos Sparse Matrices
  // ----------------------------------------------------------------------------
  void get_ptype(int* ptype)
  {
#if defined(CHECK_JAC)
    *ptype = -1;
#elif defined(IDENT_PC)
    *ptype = 0;
#elif defined(BLOCK_PC_OP) || defined(GMRES_PC_OP)
    *ptype = 1;
#elif defined(BLOCK_PC_DENSE)
    *ptype = 2;    
#elif defined(BLOCK_PC_SPARSE) || \
  defined(BLOCK_PC_SPARSE_ML)  || \
  defined(GMRES_PC_SPARSE)     || \
  defined(ML_PC)
    *ptype = 3;
#else
#error Unknown preconditioner option set in trilinosDefines.h
#endif
  }

} //extern "C"
