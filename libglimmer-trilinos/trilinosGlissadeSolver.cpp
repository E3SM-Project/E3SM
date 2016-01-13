//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                             
//   trilinosGLissadeSolver.cpp - part of the Community Ice Sheet Model (CISM)  
//                                                              
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Copyright (C) 2005-2014
//   CISM contributors - see AUTHORS file for list of contributors
//
//   This file is part of CISM.
//
//   CISM is free software: you can redistribute it and/or modify it
//   under the terms of the Lesser GNU General Public License as published
//   by the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CISM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   Lesser GNU General Public License for more details.
//
//   You should have received a copy of the Lesser GNU General Public License
//   along with CISM. If not, see <http://www.gnu.org/licenses/>.
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"
//#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#ifdef GLIMMER_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Epetra_MpiComm.h"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#include "Epetra_SerialComm.h"
#endif

#include "config.inc"

// Uncomment this #define to write out linear system
//#define WRITE_OUT_LINEAR_SYSTEM
#ifdef WRITE_OUT_LINEAR_SYSTEM
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
int solvecount=0;
#endif

// Turn this on to check validity of sparse matrix entries
//#define CHECK_FOR_ROGUE_COLUMNS

// Define variables that are global to this file.
// If this were a C++ class, these would be member data.
Teuchos::RCP<Epetra_Vector> rhs;
Teuchos::RCP<Teuchos::ParameterList> paramList;
Teuchos::RCP<Teuchos::FancyOStream> tout;
Teuchos::RCP<Epetra_CrsMatrix> matrix;
Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > linOp;
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > linOpFactory;
Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraOp;
bool successFlag = true;

// Flag for operations done once per time step (e.g. define active unknowns)
bool firstMatrixAssemblyForTimeStep = true;

// Flag for operations done once per run (e.g. read in trilinosOPtions.xml)
bool firstCallToInitializeTGS = true;

int linSolveCount=0, linSolveSuccessCount=0, linSolveIters_last=0,  linSolveIters_total=0;
double linSolveAchievedTol;
bool printLinSolDetails=true; // Need to set in input file.

extern "C" {

  // Prototypes for locally called functions
  void linSolveDetails_tgs(Thyra::SolveStatus<double>& status);
  void check_for_rogue_columns_tgs( Epetra_CrsMatrix& mat);

  //================================================================
  // This needs to be called only once per time step in the beginning 
  // to set up the owned unknow map for the problem.
  //================================================================

  void FC_FUNC(initializetgs,initializetgs) 
      (int& mySize, int* myIndicies, int* mpi_comm_f) {
    // mySize: the number of active_owned_unknowns on this processor
    // myIndicies[]: global_active_owned_unknowns integer array in glissade-speak
    // mpi_comm_f: CISM's fortran mpi communicator

    // Define output stream that only prints on Proc 0
    tout = Teuchos::VerboseObjectBase::getDefaultOStream();

#ifdef GLIMMER_MPI
    // Make sure the MPI_Init in Fortran is recognized by C++.
    // We used to call an extra MPI_Init if (!flag), but the behavior of doing so is uncertain,
    // especially if CISM's MPI communicator is a subset of MPI_COMM_WORLD (as can be the case in CESM).
    // Thus, for now, we die with an error message if C++ perceives MPI to be uninitialized.
    // If this causes problems (e.g., if certain MPI implementations seem not to recognize 
    // that MPI has already been initialized), then we will revisit how to handle this.
       int flag;
       MPI_Initialized(&flag);
       if (!flag) {
	 *tout << "ERROR in initializetgs: MPI not initialized according to C++ code" << std::endl;
	 exit(1);
       }
    MPI_Comm mpi_comm_c = MPI_Comm_f2c(*mpi_comm_f);
    Epetra_MpiComm comm(mpi_comm_c);
    Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(mpi_comm_c));
#else
    Epetra_SerialComm comm;
    Teuchos::SerialComm<int> tcomm;
#endif


    // Read parameter list from XML file once per run
    if (firstCallToInitializeTGS) {
      // Set flag so following code is executed only once per code run
      firstCallToInitializeTGS = false;
      try { 
        paramList = Teuchos::rcp(new Teuchos::ParameterList("Trilinos Options"));
        Teuchos::updateParametersFromXmlFileAndBroadcast("trilinosOptions.xml", paramList.ptr(), tcomm);

        Teuchos::ParameterList validPL("Valid List");;
        validPL.sublist("Stratimikos"); validPL.sublist("Piro");
        paramList->validateParameters(validPL, 0);

        // Set the coordinate position of the nodes for ML for repartitioning (important for #procs > 100s)
        if (paramList->sublist("Stratimikos").isParameter("Preconditioner Type")) {
           if ("ML" == paramList->sublist("Stratimikos").get<std::string>("Preconditioner Type")) {
             *tout << "\nNOTE: ML preconditioner can work much better when interface is extended\n"
                  << "\tto include Nodal XYZ coordinates.\n" << std::endl;
             Teuchos::ParameterList& mlList =
                paramList->sublist("Stratimikos").sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");
             //mlList.set("x-coordinates",myX);
             //mlList.set("y-coordinates",myY);
             //mlList.set("z-coordinates",myZ);
             mlList.set("PDE equations", 2);
           }
        }

        // Set up solver (preconditioner, iterative method) based on XML file
        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
        linearSolverBuilder.setParameterList(Teuchos::sublist(paramList, "Stratimikos"));
        linOpFactory = linearSolverBuilder.createLinearSolveStrategy("");
        linOpFactory->setOStream(tout);
        linOpFactory->setVerbLevel(Teuchos::VERB_LOW);

        linOp=Teuchos::null;
        thyraOp=Teuchos::null;
      }
      catch (std::exception& e) {
        std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n" 
             << e.what() << "\nExiting: Invalid trilinosOptions.xml file."
             << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        exit(1);
      }
      catch (...) {
        std::cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n" 
             << "\nExiting: Invalid trilinosOptions.xml file."
             << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        exit(1);
      }
    }

    // Continue setup that needs to be redone every time step
    try { 
      // Flag to let subsequent functions know that a new matrix has been created
      firstMatrixAssemblyForTimeStep = true;

      Teuchos::RCP<const Epetra_Map> rowMap = 
        Teuchos::rcp(new Epetra_Map(-1, mySize, myIndicies, 1, comm) );

      TEUCHOS_TEST_FOR_EXCEPTION(!rowMap->UniqueGIDs(), std::logic_error,
         "Error: initializetgs, myIndicies array needs to have unique entries" 
          << " across all processors.");

      // Diagnostic output for partitioning
      int minSize, maxSize;
      comm.MinAll(&mySize, &minSize, 1);
      comm.MaxAll(&mySize, &maxSize, 1);
      if (comm.MyPID()==0) 
        *tout << "\nPartition Info in init_trilinos: Total nodes = " << rowMap->NumGlobalElements()
             << "  Max = " << maxSize << "  Min = " << minSize 
             << "  Ave = " << rowMap->NumGlobalElements() / comm.NumProc() << std::endl;

      // rhs is the b vector, rhs of linear system (owned, active)
      rhs  = Teuchos::rcp(new Epetra_Vector(*rowMap));

      // Reset counters every time step: can remove these lines to have averages over entire run
      linSolveIters_total = 0;
      linSolveCount=0;
      linSolveSuccessCount = 0;

      // Construct the CrsMatrix based on the row map and bandwidth estimate
      const int bandwidth = 54;
      matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, bandwidth));
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, successFlag);
    if (!successFlag) exit(1);

    //Teuchos::TimeMonitor::summarize(*tout,false,true,false/*zero timers*/);
  }

  //============================================================
  // Insert one row of entries into matrix and RHS
  //============================================================

  void FC_FUNC(insertrowtgs,INSERTROWTGS)
	       (int& rowInd, int& numColumns, int* columns,
                double* matrixValues, double& rhsValue ) {
   // rowInd: global row number
   // numColumns: number of columns in this row (typically 54, but less on boundaries)
   // columns[]: array with numColumns valid entries of global column numbers
   // matrixValues[]: array with corresponding matrix entries
   // rhsValue: entry into "b" vector for that same row.
   //
   //TEUCHOS_FUNC_TIME_MONITOR("> insertRowTGS");

  try {
    int ierr;
    const Epetra_Map& rowMap = matrix->RowMap();

    // If this row is not owned on this processor, then throw error
    TEUCHOS_TEST_FOR_EXCEPTION(!rowMap.MyGID(rowInd), std::logic_error,
       "Error: Trilinos matrix has detected an invalid row entry (row=" 
        << rowInd << ").\n");

    // Insert contribution to rhs a.k.a. b vector (as in  Au=b)
    rhs->ReplaceGlobalValues(1, &rhsValue, &rowInd);

    if (firstMatrixAssemblyForTimeStep) {

//#define ONE_PROC_DEBUG 
#ifdef ONE_PROC_DEBUG
      if (rowMap.Comm().NumProc()==1) 
        for (int col=0; col<numColumns; col++) {
          TEUCHOS_TEST_FOR_EXCEPTION(!rowMap.MyGID(columns[col]), std::logic_error,
            "Error: Trilinos matrix has column that is not in row map: entry A(" 
            << rowInd << ", " << columns[col] << ") = " << matrixValues[col]  
            << "\n\t(This is not an error for parallel runs.");
          }
#endif
      // The matrix has not been "FillComplete()"ed. First fill of time step.
      ierr = matrix->InsertGlobalValues(rowInd, numColumns, matrixValues, columns);

      if (ierr<0) {std::cout << "Error Code for " << rowInd << "  = ("<< ierr <<")"<<std::endl; exit(1);}
      else if (ierr>0) std::cout << "Warning Code for " << rowInd << "  = ("<< ierr <<")"<<std::endl;
    }
    else {
//#define ONE_ENTRY_DEBUG 
#ifdef ONE_ENTRY_DEBUG
      for (int col=0; col<numColumns; col++) {
      ierr = matrix->ReplaceGlobalValues(rowInd, 1, &matrixValues[col], &columns[col]);

      TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
        "Error: Trilinos matrix has detected a new column entry A(" 
        << rowInd << ", " << columns[col] << ") = " << matrixValues[col]  
        << "\n\t that did not exist before.");
      }
#else
      // Subsequent matrix fills of each time step.
      ierr = matrix->ReplaceGlobalValues(rowInd, numColumns, matrixValues, columns);

      TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
        "Error: Trilinos matrix has detected a new column entry in row (" 
        << rowInd << ")\n\t that did not exist before.");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, successFlag);
  if (!successFlag) exit(1);
  }

  //============================================================
  // Call to perform solve of previously assembled linear system
  //============================================================

  void FC_FUNC(solvevelocitytgs,SOLVEVELOCITYTGS)
	       (double* velocityResult) {
   // velocityResult[]: array of length mySize from initializetgs call, that
   //                   upon return will have the velocities from Au=b solve.
   //TEUCHOS_FUNC_TIME_MONITOR("> solveVelocityTGS");

   try {
    //Teuchos::Time linearTime("LinearTime"); linearTime.start();

    // Lock in sparsity pattern of CrsMatrix -- first solve only
    if (firstMatrixAssemblyForTimeStep) {
      firstMatrixAssemblyForTimeStep = false;

      matrix->FillComplete();
#ifdef CHECK_FOR_ROGUE_COLUMNS
      check_for_rogue_columns_tgs(*matrix);
#endif

      // Associate matrix with solver strategy layers
      thyraOp = Thyra::epetraLinearOp(matrix);
    }
    // Need to do this call to invoke fresh preconditioner
    linOp = Thyra::linearOpWithSolve(*linOpFactory, thyraOp);

    // Wrap velocity vector inside Epetra Vector data structure
    Teuchos::RCP<Epetra_Vector> solution 
      = Teuchos::rcp(new Epetra_Vector(View, matrix->RowMap(), velocityResult));

#ifdef WRITE_OUT_LINEAR_SYSTEM
    solvecount++; 
    if (solvecount==1) {
      EpetraExt::RowMatrixToMatrixMarketFile("matrix1", *matrix);
      EpetraExt::MultiVectorToMatrixMarketFile("vector1", *rhs);
    }
#endif

    // Wrap Epetra Vetors as Thyra vectors, as the solver requires
    Teuchos::RCP<Thyra::MultiVectorBase<double> >
      thyraRhs = Thyra::create_Vector(rhs, thyraOp->range() );
    Teuchos::RCP<Thyra::MultiVectorBase<double> >
      thyraSol = Thyra::create_Vector(solution, thyraOp->domain() );
    Thyra::SolveStatus<double>
      status = Thyra::solve(*linOp, Thyra::NOTRANS, *thyraRhs, thyraSol.ptr());

    if (printLinSolDetails) linSolveDetails_tgs(status);

    //elapsedTime = linearTime.stop(); 
    //*tout << "Total time elapsed for calling Solve(): " << elapsedTime << std::endl;
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, successFlag);
   if (!successFlag) exit(1);
  }

  //============================================================

  void linSolveDetails_tgs(Thyra::SolveStatus<double>& status) {
    ++linSolveCount;
    bool haveData=false;
    if (status.extraParameters != Teuchos::null) {
      if (status.extraParameters->isParameter("Belos/Iteration Count")) {
        linSolveIters_last = status.extraParameters->get<int>("Belos/Iteration Count");
        linSolveIters_total += linSolveIters_last;
        haveData=true;
      }
      if (status.extraParameters->isParameter("Belos/Achieved Tolerance"))
        linSolveAchievedTol = status.extraParameters->get<double>("Belos/Achieved Tolerance");
      if (status.extraParameters->isParameter("AztecOO/Iteration Count")) {
        linSolveIters_last = status.extraParameters->get<int>("AztecOO/Iteration Count");
        linSolveIters_total += linSolveIters_last;
        haveData=true;
      }
      if (status.extraParameters->isParameter("AztecOO/Achieved Tolerance"))
        linSolveAchievedTol = status.extraParameters->get<double>("AztecOO/Achieved Tolerance");

      if (haveData) {
        *tout <<  "Precon Linear Solve ";
        if (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED)
         {*tout <<  "Succeeded: "; ++linSolveSuccessCount;}
        else  *tout <<  "Failed: ";
        *tout << std::setprecision(3) 
             << linSolveAchievedTol << " drop in " 
             << linSolveIters_last << " its (avg: " 
             << linSolveIters_total / (double) linSolveCount << " its/slv, " 
             << 100.0* linSolveSuccessCount / (double) linSolveCount << "% success)"
             << std::endl;
      }
    }
  }

  /* Debugging utility to check if columns have been Inserted into the 
   * matrix that do not correspond to a row on any processor
   */
  void check_for_rogue_columns_tgs( Epetra_CrsMatrix& mat) {
    // Set up rowVector of 0s and column vector of 1s
    const Epetra_Map& rowMap = mat.RowMap();
    const Epetra_Map& colMap = mat.ColMap();
    Epetra_Vector rowVec(rowMap); rowVec.PutScalar(0.0);
    Epetra_Vector colVec(colMap); colVec.PutScalar(1.0);
    Epetra_Import importer(colMap, rowMap);

    // Overwrite colVec 1s with rowVec 0s 
    colVec.Import(rowVec, importer, Insert);

    // Check that all 1s have been overwritten
    double nrm=0.0;
    colVec.Norm1(&nrm); // nrm = number of columns not overwritten by rows

    // If any rogue columns, exit now (or just get nans later)
    if (nrm>=1.0) {
      *tout << "ERROR: Column map has " << nrm 
           << " rogue entries that are not associated with any row." << std::endl;
       rowMap.Comm().Barrier();
       exit(-3);
    }
    else {
      *tout << "Debugging check for rogue column indices passed." 
           << " Turn off for production runs.\n" << std::endl;
    }
  }

  //============================================================

} // extern"C"
