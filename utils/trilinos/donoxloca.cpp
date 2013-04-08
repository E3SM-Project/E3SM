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
//! pal aztecoo no longer used
//#include "AztecOO.h"

// User's application specific files 
#include "noxlocainterface.hpp" // Interface file to NOX
#include "precon_interface.hpp" // Interface file to NOX
// Aaron! duplicate routine above for GMRES only preconditioner
//#include "precon_interface_gmres.hpp" // Interface file to NOX

// Required for reading and writing parameter lists from xml format
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace std;

// Define some variables that are Global to the file that
// allow NOX to be called in 3 stages (init, solve, finalize)
static Teuchos::RCP<NOX::Solver::Generic> solver;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess;
static Teuchos::RCP<Problem_Interface> interface;
static Teuchos::RCP<LOCA::GlobalData> globalData;
//pal Comm needs to be defined in this scope to be passed to preconditioner objects
static Teuchos::RCP<Epetra_MpiComm> Comm;
//pal print proc used for debug printing
static bool printproc;



// Define corresponding variables for NOX to be used by precon
// These are no longer used by preconditioner
/*
static Teuchos::RCP<NOX::Solver::Generic> solver_p;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_p;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_r;
static Teuchos::RCP<NOX::Epetra::Vector> initialGuess_v;
static Teuchos::RCP<Precon_Interface> interface_p;
static Teuchos::RCP<LOCA::GlobalData> globalData_p;
static Teuchos::RCP<LOCA::Epetra::Group> grp_p;
*/

extern "C" {

/*Defining external functions for finite difference Jacobian*/
  void calc_f(double *, double *, int, void *);
  void calc_f_lin(double *, double *, int, void *);
  
/*Defining external functions for analytic Jacobian*/
  void sw_jacobian(double *, int, double *, void *); 
/*Defining external functions for analytic Picard linearization*/
  void sw_picard(double *, int, double *, void *); 
/*Defining external functions for applying preconditioners and preconditioner blocks*/

  void simple_prec_op(double *, int, double *, void *);
  void test_id(double *, int, double *, void *);
  void picard_lin(double *, int, double *,  void *);
  void picard_diag(double *, int, double *,  void *);
  void test_picard_lin(double *, int, double *,  void *);
  void sw_picard_fdiag(double *, int, double *, void *); 
  void sw_picard_diag(double *, int, double *, void *); 
  void sw_picard_simple(double *, int, double *, void *); 
  void update_prec_state(double *, int, void *); 
  void get_jac_vector(double *, int, void *); 
  void sw_picard_block_FDiag(double *, int, double *, void *); 
  void sw_picard_block_FDiagInv(double *, int, double *, void *); 
  void sw_picard_block_ID11(double *, int, double *, void *); 
  void sw_picard_block_ID22(double *, int, double *, void *); 
  void sw_picard_block_11alt(double *, int, double *, void *); 
  void sw_picard_block_11(double *, int, double *, void *); 
  void sw_picard_block_12(double *, int, double *, void *);
  void sw_picard_block_21(double *, int, double *, void *); 
  void sw_picard_block_22(double *, int, double *, void *); 
  void sw_picard_FDFinvBt(double *, int, double *, void *); 
  void sw_picard_FDFinvBt2(double *, int, double *, void *); 
  void sw_picard_FBt(double *, int, double *, void *); 
  void sw_picard_DFinvBt(double *, int, double *, void *); 
  void sw_picard_schur(double *, int, double *, void *); 
  void sw_picard_alphaschur(double *, int, double *, void *);
// once precon_gmres is ready, use this as dummy
//  void precon_si(double *, double *, int, double*, void *);



/*pal no longer require the precon 'old' interface*/
//  void precon(double *, double *, int, double*, void *, void *);
  // Aaron! add this to be the GMRES only preconditioner
  //  void precon_wgmres(double *, double *, int, double*, void *, void *);
  // when you create this, the bind_C name needs to be precon_wgmres
}

extern "C" { 

  void printNewtonKrylovStats(Teuchos::ParameterList& nlParams);

  void printPreconStats(Teuchos::ParameterList& lsParams_p);

  void noxinit(int* nelems, double* statevector, int* mpi_comm_ignored,
               void* blackbox_res, void* blackbox_prec, void* jac_data)
  {

    void (*residualFunction)(double *, double *, int, void *) = calc_f;
    //pal the 'old' precon function is obsolete


    //void (*precFunction)(double *, double *, int, double*, void *, void *) = precon;

    //void (*jacFunction)(double *,  int, double*,  void *) = sw_picard;
    void (*jacFunction)(double *,  int, double*,  void *) = sw_jacobian;
    //void (*jacFunction)(double *,  int, double*,  void *) = picard_lin;
    //void (*jacFunction)(double *,  int, double*,  void *) = test_picard_lin;
    //void (*jacFunction)(double *, int, double*,  void *) = test_id;
    
    //void (*precFunction)(double *, int, double*,  void *) = simple_prec_op;
    void (*precFunction)(double *, int, double*,  void *) = test_id;
    
    //void (*precFunction)(double *, int, double*,  void *) = sw_picard;
    //void (*precFunction)(double *, int, double*,  void *) = sw_picard_fdiag;
    //void (*precFunction)(double *, int, double*,  void *) = sw_picard_diag;
    //void (*precFunction)(double *, int, double*,  void *) = sw_jacobian;
    
    
    //void (*precFunction)(double *, int, double*,  void *) = sw_picard_simple;
    //void (*precFunction)(double *,  int, double*,  void *) = picard_lin;
    //void (*precFunction)(double *,  int, double*,  void *) = picard_diag;
    //void (*precFunction)(double *,  int, double*,  void *) = test_picard_lin;
    
    //void (*precFunction)(double *,  int, double*,  void *) = sw_jacobian;
    void (*precUpdateFunction)(double *, int, void *)=update_prec_state;
    void (*getJacVector)(double *, int, void *)=get_jac_vector;
    
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_test_mass;
    
    /*Convention: for SIMPLE algorithm based on Picard linearization
    11 block is F
    12 block is diag(F)^{-1}B'
    21 block is B
    22 block is S=G-Bdiag(F)^{-1}B' */
    
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11alt;
    void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_FDiagInv;
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_FDiag;
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_ID11;
    void (*precFunctionblock12)(double *, int, double*,  void *) = sw_picard_DFinvBt;
    //void (*precFunctionblock12)(double *, int, double*,  void *) = sw_picard_block_12;
    void (*precFunctionblock21)(double *, int, double*,  void *) = sw_picard_block_21;
    void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_schur;
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_ID22;
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_12;//Testing only
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_22;
    
    //Convention: for Picard linearization preconditioner
    //11 block is F
    //12 block is B'
    //21 block is B
    //22 block is G
    
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    //void (*precFunctionblock12)(double *, int, double*,  void *) = sw_picard_block_12;
    //void (*precFunctionblock21)(double *, int, double*,  void *) = sw_picard_block_21;
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_22;
    
    
                    //Convention: for alpha=1 test for SIMPLE


    //11 block is F
    //12 block is FDiagFinvB'
    //21 block is B
    //22 block is G
    
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    //void (*precFunctionblock12)(double *, int, double*,  void *) = sw_picard_FDFinvBt;
    //void (*precFunctionblock21)(double *, int, double*,  void *) = sw_picard_block_21;
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_22;
    
    //Convention: Test for SIMPLE any alpha
    //11 block is F
    //12 block is FDiagFinvB'
    //21 block is B
    //22 block is G
    
    //void (*precFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    //void (*precFunctionblock12)(double *, int, double*,  void *) = sw_picard_FDFinvBt;
    //void (*precFunctionblock21)(double *, int, double*,  void *) = sw_picard_block_21;
    //void (*precFunctionblock22)(double *, int, double*,  void *) = sw_picard_alphaschur;
    
    //void (*auxprecFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_ID11;
    void (*auxprecFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_11;
    //void (*auxprecFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_FDiagInv;
    //void (*auxprecFunctionblock11)(double *, int, double*,  void *) = sw_picard_block_FDiag;
    //void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_FDFinvBt;
    void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_FDFinvBt2;
    //void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_DFinvBt;
    //void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_block_12;
    // void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_FBt;
    //void (*auxprecFunctionblock12)(double *, int, double*,  void *) = sw_picard_DFinvBt;
    void (*auxprecFunctionblock21)(double *, int, double*,  void *) = sw_picard_block_21;
    //void (*auxprecFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_22; 
    void (*auxprecFunctionblock22)(double *, int, double*,  void *) = sw_picard_alphaschur;
    //void (*auxprecFunctionblock22)(double *, int, double*,  void *) = sw_picard_schur;
    //void (*auxprecFunctionblock22)(double *, int, double*,  void *) = sw_picard_block_ID22;




    try {

      // Create a communicator for Epetra objects
//      Epetra_MpiComm Comm( MPI_COMM_WORLD );
//

      int N=*nelems;
      // Create a communicator for Epetra objects
      // Epetra_MpiComm Comm( MPI_COMM_WORLD );
      
//        cout<<"setting Comm"<<endl;
      
        Comm=Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      
 //       cout<<" Comm set"<<endl;
      
         if (Comm->MyPID()==0) printproc=true;
                 else   printproc=false;
      


      // Get the process ID and the total number of processors
      //
      int MyPID = Comm->MyPID();
      //if(printproc)cout<<"mpid="<<MyPID<<endl;
     
      // int MyPID = Comm.MyPID();
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
      //
      Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
      Teuchos::ParameterList& nlDir = nlParams.sublist("Direction");
      Teuchos::ParameterList& nlNewton = nlDir.sublist("Newton");
      //
     // Teuchos::ParameterList& stratParams = nlNewton.sublist("Stratimikos");
      Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
      nlPrintParams.set("MyPID", MyPID);
      if (!nlPrintParams.isParameter("Output Information"))
        //     nlPrintParams.set("Output Information",0
        nlPrintParams.set("Output Information",
            NOX::Utils::OuterIteration  +
            // NOX::Utils::OuterIterationStatusTest +
            // NOX::Utils::InnerIteration +
            // NOX::Utils::Details +
            NOX::Utils::LinearSolverDetails 
            // NOX::Utils::Warning 
            // NOX::Utils::StepperIteration +
            // NOX::Utils::StepperDetails +
            // NOX::Utils::StepperParameters
            );




/* Create the interface between the test problem and the nonlinear solver
   This is created by the user using inheritance of the abstract base class:
*/ 



      //FD Jacobian
      
      #if 1 
        interface = Teuchos::rcp(new Problem_Interface(N, statevector, pVector, *Comm, blackbox_res, blackbox_prec, residualFunction, precFunction,precUpdateFunction));                                                                                
      #endif
       	 


      //Analytic Jacobian
      #if 0 
      interface = Teuchos::rcp(new Problem_Interface(N, statevector, pVector, *Comm, blackbox_res, blackbox_prec,jac_data, residualFunction, precFunction,jacFunction,precUpdateFunction,getJacVector));
      #endif

      /* Interface for SIMPLE preconditioner */
      #if 0
      interface = Teuchos::rcp(new Problem_Interface(N, statevector, pVector, *Comm, 
			      blackbox_res, blackbox_prec,jac_data, 
			      residualFunction, 
			      precFunctionblock11,precFunctionblock12, 
			      precFunctionblock21,precFunctionblock22, 
			      jacFunction,precUpdateFunction,getJacVector)); 
      #endif
      //                                                                                                                                                                                                                                                                                                                                                                   /* interface for comparing two different formulations of SIMPLE  */
       #if 0
       interface = Teuchos::rcp(new Problem_Interface(N, statevector, pVector, *Comm, 
			       blackbox_res, blackbox_prec,jac_data, 
			       residualFunction, 
			       precFunctionblock11,precFunctionblock12, 
			       precFunctionblock21,precFunctionblock22, 
			       auxprecFunctionblock11,auxprecFunctionblock12, 
			       auxprecFunctionblock21,auxprecFunctionblock22, 
			       jacFunction,precUpdateFunction,getJacVector));
        #endif

       // interface->printCommID();
       
  //     if (printproc) cout<<"setup NOX interface components"<<endl<<flush;
       
       // Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
       // LOCA 2
       Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;
   //    if (printproc) cout<<"iReq"<<endl<<flush;
       
       // Interface is inherited from these 2 classes as well for user prec
       Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;
       Teuchos::RCP<Epetra_Operator> precOperator = interface;
    //   if (printproc) cout<<"iPrec"<<endl<<flush;
     //  if (printproc) cout<<"precOperator"<<endl<<flush;
       
       // Create the Epetra_Vector for the  state vector
       Teuchos::RCP<Epetra_Vector> soln = interface->getVector();



/*
    Epetra_Vector Ain(*soln);
        if (printproc){
        cout<<"Asoln=["<<endl;
        for(int i=0; i<N;i++) cout<<*(soln->Values()+i)<<endl;
        cout<<"];"<<flush;
        }
*/
       /*Finite Difference Jacobian*/
#if 1       
	   Teuchos::RCP<NOX::Epetra::MatrixFree> FD =
	       Teuchos::rcp(new NOX::Epetra::MatrixFree(nlPrintParams, interface, soln));
	           //Teuchos::rcp(new NOX::Epetra::MatrixFree(nlPrintParams, interface, soln));
	             Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = FD;
#endif	
        


       /*
	   Epetra_Vector B(*soln);
	//soln->PutScalar(0.0);
           if (printproc){
		   cout<<"B=["<<endl;
		   for(int i=0; i<N;i++) cout<<*(B(0)->Values()+i)<<endl;
		   cout<<"];"<<flush;
	   }
	   */


//Analytic Jacobian
#if 0
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<Epetra_Operator> FD = interface;
#endif


      //Create the linear systems & give Strat the matrix free Jacobian, precon and precon operator
      Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> linsys =
          Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(nlPrintParams, nlNewton,
          //Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(nlPrintParams, stratParams,
                                                   iJac, FD, iPrec, precOperator, soln, true));

      // Create the loca vector
      initialGuess = Teuchos::RCP<NOX::Epetra::Vector>(new NOX::Epetra::Vector(soln));

      // LOCA specific, but here grp is one branch of eventual suite of parameter spaces
      // Create Epetra factory
      Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory = Teuchos::rcp(new LOCA::Epetra::Factory);

      // Create global data object
      globalData = LOCA::createGlobalData(paramList, epetraFactory);

      // Create the Group
      Teuchos::RCP<LOCA::Epetra::Group> grp = 
        Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, 
              iReq, *initialGuess, linsys, pVector)); 

      // extra calc of F
      //grp->disableLinearResidualComputation(true);

      // Create the Solver convergence test
      Teuchos::RCP<NOX::StatusTest::NormF> wrms;
      const double wrms_tol = nlParams.get("Convergence Tolerance",1.0e-10);
      wrms = Teuchos::rcp(new NOX::StatusTest::NormF(wrms_tol));
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


  void noxsolve(int* nelems, double* statevector, void* blackbox_res, void* blackbox_prec, void* jac_data)
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
      interface->resetBlackbox(blackbox_res, blackbox_prec, jac_data);

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
      const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(
                                            solver->getSolutionGroup().getX())).getEpetraVector();
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
/*
  void precon_solve(int* nelems, double* statevector, double* rhs, 
      void* blackbox_res, void* blackbox_prec)
  {

    try {
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(solver_p), logic_error, 
          "Exception: precon_solve called with solver=null: "
          << "either did call init_prec 1st, or called finish_prec already");

      // reset solver with new initial Guess
      Epetra_Vector& initialGuessEVp = initialGuess_p->getEpetraVector();
      for (int i=0; i<*nelems; i++) initialGuessEVp[i] = statevector[i];
      // maybe a redundancy - check
      solver_p->reset(*initialGuess_p);

      // convert rhs from double * to epetra vector(just like for statevector)
      Epetra_Vector& initialGuessEVr = initialGuess_r->getEpetraVector();
      for (int i=0; i<*nelems; i++) initialGuessEVr[i] = rhs[i];

      // convert solution to appplyJI from double * to epetra vector
      Epetra_Vector& initialGuessEVv = initialGuess_v->getEpetraVector();
      for (int i=0; i<*nelems; i++) initialGuessEVv[i] = rhs[i];

      // reset interface with new blackbox 
      interface_p->resetBlackbox(blackbox_res, blackbox_prec);

      //    Teuchos::ParameterList& lsParams_p =
      //         const_cast<Teuchos::ParameterList&>(solver_p->getList());
      //    lsParams_p.sublist("Linear Solver").sublist("Output").
      //                        set("Total Number of Linear Iterations", 1);

      Teuchos::ParameterList& nlParams_p =
        const_cast<Teuchos::ParameterList&>(solver_p->getList());
      nlParams_p.sublist("Direction").sublist("Newton").
                 sublist("Linear Solver").sublist("Output").
                 set("Total Number of Linear Iterations", 0);

      Teuchos::ParameterList& lsParams_p =
        nlParams_p.sublist("Direction").sublist("Newton").sublist("Linear Solver");

      // call Apply Jacobian inverse method with z as the argument 
      grp_p->setX(*initialGuess_p);
      NOX::Abstract::Group::ReturnType stat1 = grp_p->computeJacobian();
      NOX::Abstract::Group::ReturnType status = 
        grp_p->applyJacobianInverse(lsParams_p,*initialGuess_r,*initialGuess_v); 

      printPreconStats(nlParams_p);

      // convert solution to appplyJI from double * to epetra vector
      for (int i=0; i<*nelems; i++) rhs[i] = initialGuessEVv[i];
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
*/

  void noxfinish()
  {
    LOCA::destroyGlobalData(globalData);
    initialGuess=Teuchos::null;
    interface=Teuchos::null;
    solver=Teuchos::null;
  }

/*
  void finish_prec()
  {
    LOCA::destroyGlobalData(globalData_p);
    initialGuess_p=Teuchos::null;
    initialGuess_r=Teuchos::null;
    initialGuess_v=Teuchos::null;
    interface_p=Teuchos::null;
    solver_p=Teuchos::null;
  }
*/

  void printNewtonKrylovStats(Teuchos::ParameterList& nlParams)
  {
    static int totalNewtonIters=0;
    static int totalKrylovIters=0;
    static int stepNum=0;
    int NewtonIters = nlParams.sublist("Output").get("Nonlinear Iterations", -1000);
// replace with thyra/belos version
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
    if (stepNum > 1) {
      globalData->locaUtils->out() << "Convergence Stats: running total: Nn, Kr, Kr/Ne, Kr/Step: "
                                   << totalNewtonIters << "  " << totalKrylovIters << "  "
                                   << (double) totalKrylovIters / (double) totalNewtonIters
                                   << "  " << (double) totalKrylovIters / (double) stepNum << endl;
    }
  }
/*
  void printPreconStats(Teuchos::ParameterList& lsParams_p)
  {
    static int totalKrylovIters=0;
    static int stepNum=0;
    int KrylovIters = lsParams_p.sublist("Direction").sublist("Newton").
                                 sublist("Linear Solver").sublist("Output").
                                 get("Total Number of Linear Iterations", -1000);
    totalKrylovIters += KrylovIters;
    stepNum++;
    globalData_p->locaUtils->out() << "Precon Convergence Stats: for step  #" << stepNum
                                   << " : Krylov " << KrylovIters << endl;
    if (stepNum > 1) {
      globalData_p->locaUtils->out() << "Precon Convergence Stats: running Krylov total: "
                                     << totalKrylovIters << "  " << endl;
    }
  }
*/

} //extern "C"
