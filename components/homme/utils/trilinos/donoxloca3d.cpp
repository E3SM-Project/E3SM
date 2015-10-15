#include "LOCA.H"
#include "LOCA_Epetra.H"
//STRAT1
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

// Trilinos Objects
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// User's application specific files 
#include "noxlocainterface.hpp" // Interface file to NOX
//#include "precon_interface.hpp" // Interface file to NOX
// Aaron! duplicate routine above for GMRES only preconditioner
//#include "precon_interface_gmres.hpp" // Interface file to NOX

// Required for reading and writing parameter lists from xml format
// Configure Trilinos with --enable-teuchos-extended
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

using namespace std;

// Define some variables that are Global to the file that
// allow NOX to be called in 3 stages (init, solve, finalize)
static Teuchos::RCP<NOX::Solver::Generic> solver;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess;
static Teuchos::RCP<Problem_Interface> interface;
static Teuchos::RCP<LOCA::GlobalData> globalData;
//  Define corresponding variables for NOX to be used by precon
//static Teuchos::RCP<NOX::Solver::Generic> solver_p;
//static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_p;
//static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_r;
//static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_v;
//static Teuchos::RCP<Precon_Interface> interface_p;
//static Teuchos::RCP<LOCA::GlobalData> globalData_p;
//static Teuchos::RCP<LOCA::Epetra::Group> grp_p;

extern "C" {
  void calc_f(double *, double *, int, void *);
//  void calc_f_lin(double *, double *, int, void *);
  void precon(double *, double *, int, double*, void *, void *);
// Aaron! add this to be the GMRES only preconditioner
//  void precon_wgmres(double *, double *, int, double*, void *, void *);
// when you create this, the bind_C name needs to be precon_wgmres

// once precon_gmres is ready, use this as dummy
//  void precon_si(double *, double *, int, double*, void *);
//  void precon(double *, double *, int, double*, void *);
}

extern "C" { 

 void printNewtonKrylovStats(Teuchos::ParameterList& nlParams);

// void printPreconStats(Teuchos::ParameterList& lsParams_p);

void noxinit(int* nelems, double* statevector, int* mpi_comm_ignored,
             void* blackbox_res, void* blackbox_prec)
{

  void (*residualFunction)(double *, double *, int, void *) = calc_f;
//  void (*precFunction)(double *, double *, int, double*, void *) = precon;
  void (*precFunction)(double *, double *, int, double*, void *, void *) = precon;

 try {

  // Create a communicator for Epetra objects
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  //int NumProc = Comm.NumProc();

  // Begin LOCA Solver ************************************
  //
  // Create and initialize the continuation/bifurcation parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("ContinuationParam", 0.0);

  // Create parameter (options) list
  Teuchos::RCP<Teuchos::ParameterList> paramList = 
    Teuchos::rcp(new Teuchos::ParameterList);

  // Read in the parameter list from a file
  Teuchos::updateParametersFromXmlFile("input.xml", paramList.ptr());

  // Set some default parameters
  Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Parameter", "ContinuationParam");
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
  Teuchos::ParameterList& nlDir = nlParams.sublist("Direction");
  Teuchos::ParameterList& nlNewton = nlDir.sublist("Newton");
  Teuchos::ParameterList& stratParams = nlNewton.sublist("Stratimikos");
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", MyPID);
   if (!nlPrintParams.isParameter("Output Information"))
//     nlPrintParams.set("Output Information",0
     nlPrintParams.set("Output Information",
                        NOX::Utils::OuterIteration  +
//                      NOX::Utils::OuterIterationStatusTest +
//                      NOX::Utils::InnerIteration +
//                      NOX::Utils::Details +
                        NOX::Utils::LinearSolverDetails 
//                      NOX::Utils::Warning 
//                      NOX::Utils::StepperIteration +
//                      NOX::Utils::StepperDetails +
//                      NOX::Utils::StepperParameters
                              );


  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  interface = Teuchos::rcp(new Problem_Interface(nelems, statevector, pVector, Comm,
                                       blackbox_res, blackbox_prec,
                                       residualFunction, precFunction));
//  interface = Teuchos::rcp(new Problem_Interface(nelems, statevector, pVector, Comm,
//                                       blackbox_res, blackbox_prec,
//                                       residualFunction));
  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

  // Interface is inherited from these 2 classes as well for user prec
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
  Teuchos::RCP<Epetra_Operator> precOperator = interface;
  
  // Create the Epetra_Vector for the  state vector
  Teuchos::RCP<Epetra_Vector> soln = interface->getVector();

  Teuchos::RCP<NOX::Epetra::MatrixFree> FD =
    Teuchos::rcp(new NOX::Epetra::MatrixFree(nlPrintParams, interface, soln));

  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = FD;

   //Create the linear systems
//  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys = 
//    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
//  				         iJac, FD, iPrec, precOperator, soln));

//STRAT3 giving Strat the matrix free Jacobian, precon and precon operator
  Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> linsysstrat = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(nlPrintParams, stratParams,
				         iJac, FD, iPrec, precOperator, soln, true));

  // Create the loca vector
  initialGuess = Teuchos::RCP<NOX::Epetra::Vector>(new NOX::Epetra::Vector(soln));

// LOCA specific, but here grp is one branch of eventual suite of parameter spaces
  // Create Epetra factory
  Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  globalData = LOCA::createGlobalData(paramList, epetraFactory);

  // Create the Group
//STRAT4
  Teuchos::RCP<LOCA::Epetra::Group> grp = 
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, 
			 iReq, *initialGuess, linsysstrat, pVector)); 

//grp->disableLinearResidualComputation(true);

  // Create the Solver convergence test
  Teuchos::RCP<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Convergence Tolerance",1.0e-8)));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

  Teuchos::RCP<Teuchos::ParameterList> nlParamsRCP = Teuchos::rcp(&nlParams, false);
  solver = NOX::Solver::buildSolver(grp, combo, nlParamsRCP);

//FD->setSolverForComputeJacobian(solver);

 } //end try block

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

void noxsolve(int* nelems, double* statevector,
             void* blackbox_res, void* blackbox_prec)
{

  try {
    TEUCHOS_TEST_FOR_EXCEPTION(is_null(solver), logic_error, 
                               "Exception: noxsolve called with solver=null: "
                               << "either did call noxinit 1st, or called noxfinish already");

    // reset solver with new initial Guess
    Epetra_Vector& initialGuessEV = initialGuess->getEpetraVector();
    for (int i=0; i<*nelems; i++) initialGuessEV[i] = statevector[i];
    solver->reset(*initialGuess);

    // reset interface with new blackbox
    interface->resetBlackbox(blackbox_res, blackbox_prec);

    // reset counter for iteration statistics
    Teuchos::ParameterList& nlParams =
         const_cast<Teuchos::ParameterList&>(solver->getList());
    nlParams.sublist("Direction").sublist("Newton").
                        sublist("Stratimikos").sublist("Output").
                        set("Total Number of Linear Iterations", 0);

    NOX::StatusTest::StatusType status = solver->solve();

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
        << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl;
      solver->getList().print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    printNewtonKrylovStats(nlParams);

    // copy out final solution
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>
          (solver->getSolutionGroup().getX())).getEpetraVector();
    for (int i=0; i<*nelems; i++) statevector[i] = finalSolution[i];
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
  LOCA::destroyGlobalData(globalData);
  initialGuess=Teuchos::null;
  interface=Teuchos::null;
  solver=Teuchos::null;
}

void printNewtonKrylovStats(Teuchos::ParameterList& nlParams)
{
  static int totalNewtonIters=0;
  static int totalKrylovIters=0;
  static int stepNum=0;
  int NewtonIters = nlParams.sublist("Output").get("Nonlinear Iterations", -1000);
  int KrylovIters = nlParams.sublist("Direction").sublist("Newton").
                    sublist("Linear Solver").sublist("Output").
                    get("Total Number of Linear Iterations", -1000);
  totalNewtonIters += NewtonIters;
  totalKrylovIters += KrylovIters;
  stepNum++;
  globalData->locaUtils->out() << "Convergence Stats: for step  #" << stepNum
       << " : Newton, Krylov, Kr/Ne: "
       << NewtonIters << "  " << KrylovIters << "  "
       << (double) KrylovIters / (double) (NewtonIters+1.0e-14) << endl;
  if (stepNum > 1)
  globalData->locaUtils->out() << "Convergence Stats: running total: Nn, Kr, Kr/Ne, Kr/Step: "
       << totalNewtonIters << "  " << totalKrylovIters << "  "
       << (double) totalKrylovIters / (double) totalNewtonIters
       << "  " << (double) totalKrylovIters / (double) stepNum << endl;
}

} //extern "C"
