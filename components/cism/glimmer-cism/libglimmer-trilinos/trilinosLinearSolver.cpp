//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                             
//   trilinosLinearSolver.cpp - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
//                                                              
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   Copyright (C) 2005-2013
//   Glimmer-CISM contributors - see AUTHORS file for list of contributors
//
//   This file is part of Glimmer-CISM.
//
//   Glimmer-CISM is free software: you can redistribute it and/or modify it
//   under the terms of the Lesser GNU General Public License as published
//   by the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Glimmer-CISM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   Lesser GNU General Public License for more details.
//
//   You should have received a copy of the Lesser GNU General Public License
//   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CombineMode.h"
#include "matrixInterface.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#ifdef GLIMMER_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
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
#define CHECK_FOR_ROGUE_COLUMNS

// Define variables that are global to this file.
// If this were a C++ class, these would be member data.
Teuchos::RCP<TrilinosMatrix_Interface> interface;
Teuchos::RCP<Epetra_CrsMatrix> savedMatrix_A;
Teuchos::RCP<Epetra_CrsMatrix> savedMatrix_C;
Teuchos::RCP<Epetra_Vector> soln;
Teuchos::RCP<Teuchos::ParameterList> pl;
Teuchos::RCP<Teuchos::FancyOStream> out;
Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;
Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraOper;
bool success = true;

int linearSolveCount=0, linearSolveSuccessCount=0, linearSolveIters_last=0,  linearSolveIters_total=0;
double linearSolveAchievedTol;
bool printDetails=true; // Need to set in input file.

extern "C" {

  // Prototype for locally called function
  void linSolveDetails(Thyra::SolveStatus<double>& status);
  void check_for_rogue_columns( Epetra_CrsMatrix& mat);

  //================================================================
  //================================================================
  // RN_20091215: This needs to be called only once per time step
  // in the beginning to set up the problem.
  //================================================================
  void FC_FUNC(inittrilinos,INITTRILINOS) (int& bandwidth, int& mySize,
	       int* myIndicies, double* myX, double* myY, double* myZ,
	       int* mpi_comm_f) {
// mpi_comm_f: CISM's fortran mpi communicator

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
	 cout << "ERROR in inittrilinos: MPI not initialized according to C++ code" << endl;
	 exit(1);
       }
    MPI_Comm mpi_comm_c = MPI_Comm_f2c(*mpi_comm_f);
    Epetra_MpiComm comm(mpi_comm_c);
    Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(mpi_comm_c));
#else
    Epetra_SerialComm comm;
    Teuchos::SerialComm<int> tcomm;
#endif

    Teuchos::RCP<const Epetra_Map> rowMap = 
      Teuchos::rcp(new Epetra_Map(-1,mySize,myIndicies,1,comm) );

    TEUCHOS_TEST_FOR_EXCEPTION(!rowMap->UniqueGIDs(), std::logic_error,
       "Error: inittrilinos, myIndices array needs to have Unique entries" 
        << " across all processor.");

    // Diagnostic output for partitioning
    int minSize, maxSize;
    comm.MinAll(&mySize, &minSize, 1);
    comm.MaxAll(&mySize, &maxSize, 1);
    if (comm.MyPID()==0) 
      cout << "\nPartition Info in init_trilinos: Total nodes = " << rowMap->NumGlobalElements()
           << "  Max = " << maxSize << "  Min = " << minSize 
           << "  Ave = " << rowMap->NumGlobalElements() / comm.NumProc() << endl;

    soln = Teuchos::rcp(new Epetra_Vector(*rowMap));

    // Read parameter list once
    try { 
       pl = Teuchos::rcp(new Teuchos::ParameterList("Trilinos Options"));
       Teuchos::updateParametersFromXmlFileAndBroadcast("trilinosOptions.xml", pl.ptr(), tcomm);

       Teuchos::ParameterList validPL("Valid List");;
       validPL.sublist("Stratimikos"); validPL.sublist("Piro");
       pl->validateParameters(validPL, 0);
    }
    catch (std::exception& e) {
      cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n" 
           << e.what() << "\nExiting: Invalid trilinosOptions.xml file."
           << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
      exit(1);
    }
    catch (...) {
      cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n" 
           << "\nExiting: Invalid trilinosOptions.xml file."
           << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
      exit(1);
    }

    try { 
      // Set the coordinate position of the nodes for ML for repartitioning (important for #procs > 100s)
      if (pl->sublist("Stratimikos").isParameter("Preconditioner Type")) {
         if ("ML" == pl->sublist("Stratimikos").get<string>("Preconditioner Type")) {
           Teuchos::ParameterList& mlList =
              pl->sublist("Stratimikos").sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");
           mlList.set("x-coordinates",myX);
           mlList.set("y-coordinates",myY);
           mlList.set("z-coordinates",myZ);
           mlList.set("PDE equations", 1);
         }
      }

      out = Teuchos::VerboseObjectBase::getDefaultOStream();

      // Reset counters every time step: can remove these lines to have averages over entire run
      linearSolveIters_total = 0;
      linearSolveCount=0;
      linearSolveSuccessCount = 0;

      // Create an interface that holds a CrsMatrix instance and some useful methods.
      interface = Teuchos::rcp(new TrilinosMatrix_Interface(rowMap, bandwidth, comm));

      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      linearSolverBuilder.setParameterList(Teuchos::sublist(pl, "Stratimikos"));
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
      lowsFactory->setOStream(out);
      lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

      lows=Teuchos::null;
      thyraOper=Teuchos::null;
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) exit(1);
  }

  //============================================================
  // RN_20091118: This is to update the matrix with new entries.
  //============================================================

  void FC_FUNC(putintotrilinosmatrix,PUTINTOTRILINOSMATRIX)
	       (int& rowInd, int& colInd, double& val) {

  try {
    int ierr;
    const Epetra_Map& map = interface->getRowMap();
    // If this row is not owned on this processor, then throw error
    TEUCHOS_TEST_FOR_EXCEPTION(!map.MyGID(rowInd), std::logic_error,
       "Error: Trilinos matrix has detected an invalide row entry (row=" 
        << rowInd << ",col=" << colInd << ",val=" << val << ").\n");

    Epetra_CrsMatrix& matrix = *(interface->getOperator());

    if (!interface->isSparsitySet()) {

      // The matrix has not been "FillComplete()"ed. First fill of time step.
      ierr = matrix.InsertGlobalValues(rowInd, 1, &val, &colInd);
      if (ierr<0) {cout << "Error Code for " << rowInd << "  " << colInd << "  = ("<< ierr <<")"<<endl; exit(1);}
      else if (ierr>0) cout << "Warning Code for " << rowInd << "  " << colInd << "  = ("<< ierr <<")"<<endl;
    }
    else {
      // Subsequent matrix fills of each time step.
       ierr = matrix.ReplaceGlobalValues(rowInd, 1, &val, &colInd);

      TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
	 "Error: Trilinos matrix has detected a new entry (" 
             << rowInd << ", " << colInd << ", " << val
             << ")\n\t that did not exist before.");
    }
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }

  //========================================================
  // RN_20091118: This is to make calls to Trilinos solvers.
  //========================================================

  void FC_FUNC(solvewithtrilinos,SOLVEWITHTRILINOS)
	       (double* rhs, double* answer, double& elapsedTime) {
   try {
    //Teuchos::Time linearTime("LinearTime"); linearTime.start();

    // Lock in sparsity pattern
    if (!interface->isSparsitySet()) {
      interface->finalizeSparsity();
#ifdef CHECK_FOR_ROGUE_COLUMNS
      check_for_rogue_columns(*interface->getOperator());
#endif
    }

    const Epetra_Map& map = interface->getRowMap(); 
    Teuchos::RCP<Epetra_Vector> epetraSol = soln;
    Teuchos::RCP<Epetra_Vector> epetraRhs;
    epetraRhs = Teuchos::rcp(new Epetra_Vector(View, map, rhs));

    thyraOper = Thyra::epetraLinearOp(interface->getOperator());
    Teuchos::RCP<Thyra::MultiVectorBase<double> >
      thyraRhs = Thyra::create_Vector(epetraRhs, thyraOper->range() );
    Teuchos::RCP<Thyra::MultiVectorBase<double> >
      thyraSol = Thyra::create_Vector(epetraSol, thyraOper->domain() );

    lows = Thyra::linearOpWithSolve(*lowsFactory, thyraOper);

    // Uncomment following block to Dump out two matrices Avv, Auu. 
    // This function is called twice per Picard iter, which is twice
    // per outer GMRES step for Newton solves, so writing at 
    // solvecount==1 is first system, solvecount==51 is 26th Picard iter.
    
#ifdef WRITE_OUT_LINEAR_SYSTEM
    solvecount++; 
    if (solvecount==1) {
      EpetraExt::RowMatrixToMatrixMarketFile("matrix1", *interface->getOperator());
      EpetraExt::MultiVectorToMatrixMarketFile("vector1", *epetraRhs);
    }
#endif

    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows, Thyra::NOTRANS, *thyraRhs, thyraSol.ptr());

    if (printDetails) linSolveDetails(status);

    soln->ExtractCopy(answer);

    //elapsedTime = linearTime.stop(); *out << "Total time elapsed for calling Solve(): " << elapsedTime << endl;
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }


  void FC_FUNC(savetrilinosmatrix,SAVETRILINOSMATRIX) (int* i) {
   try {
    if (!interface->isSparsitySet()) {
      interface->finalizeSparsity();
#ifdef CHECK_FOR_ROGUE_COLUMNS
      check_for_rogue_columns(*interface->getOperator());
#endif
    }
    if (*i==0)
      savedMatrix_A = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
    else if (*i==1)
      savedMatrix_C = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
    else if (*i==2) {
      savedMatrix_A = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
      savedMatrix_C = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
    }
    else
      assert(false);
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }


  void FC_FUNC(restoretrilinosmatrix,RESTORTRILINOSMATRIX) (int* i) {
   try {
    if (*i==0)
      interface->updateOperator(savedMatrix_A);
    else if (*i==1)
      interface->updateOperator(savedMatrix_C);
    else
      assert(false);
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }

  void FC_FUNC(matvecwithtrilinos,MATVECWITHTRILINOS)
	       (double* x, double* answer) {
   try {
    const Epetra_Map& map = interface->getRowMap(); 

    Teuchos::RCP<Epetra_Vector> epetra_x;
    epetra_x  = Teuchos::rcp(new Epetra_Vector(View, map, x));

    Epetra_Vector y(map);
    interface->getOperator()->Multiply(false, *epetra_x, y);

    y.ExtractCopy(answer);
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }


  //============================================================
  // Functionality here is for FEM fills. These differ in that 
  // contributions to matrix entried can come in multiple parts,
  // so we need to ZeroOut and SumInto the matris, instead of 
  // Replace matrix entries.
  //
  // This first attempt will not work in parallel -- we need to
  // add functionality to deal with off-processor contributions.
  //============================================================

  void FC_FUNC(zeroouttrilinosmatrix,ZEROOUTTRILINOSMATRIX)() {
   try {
    // Zero out matrix. Don't do anything for first call, when matrix is empty.
    if (interface->isSparsitySet()) {
      Epetra_CrsMatrix& matrix = *(interface->getOperator());
      matrix.PutScalar(0.0);
    }
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }

  void FC_FUNC(sumintotrilinosmatrix,SUMINTOTRILINOSMATRIX)
	       (int& rowInd, int& numEntries, int* colInd, double* val) {

   try {
    const Epetra_Map& map = interface->getRowMap();

    Epetra_CrsMatrix& matrix = *(interface->getOperator());

    if (!interface->isSparsitySet()) {
      // The matrix has not been "FillComplete()"ed. First fill of time step.
      // Inserted values at this stage will be summed together later
      int ierr = matrix.InsertGlobalValues(rowInd, numEntries, val, colInd);
      if (ierr<0) {cout << "Error Code for " << rowInd << "  " << colInd[0] << "  = ("<< ierr <<")"<<endl; exit(1);}
      else if (ierr>0) cout << "Warning Code for " << rowInd << "  " << colInd[0] << "  = ("<< ierr <<")"<<endl;
    }
    else {
      // Subsequent matrix fills of each time step.
      int ierr = matrix.SumIntoGlobalValues(rowInd, numEntries, val, colInd);
    
      TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
	 "Error: Trilinos matrix has detected a new entry (" 
             << rowInd << ", " << colInd[0] << ", " << val[0] 
             << ")\n\t that did not exist before.");
    }
   }
   TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
   if (!success) exit(1);
  }

  void linSolveDetails(Thyra::SolveStatus<double>& status) {
    ++linearSolveCount;
    bool haveData=false;
    if (status.extraParameters != Teuchos::null) {
      if (status.extraParameters->isParameter("Belos/Iteration Count")) {
        linearSolveIters_last = status.extraParameters->get<int>("Belos/Iteration Count");
        linearSolveIters_total += linearSolveIters_last;
        haveData=true;
      }
      if (status.extraParameters->isParameter("Belos/Achieved Tolerance"))
        linearSolveAchievedTol = status.extraParameters->get<double>("Belos/Achieved Tolerance");
      if (status.extraParameters->isParameter("AztecOO/Iteration Count")) {
        linearSolveIters_last = status.extraParameters->get<int>("AztecOO/Iteration Count");
        linearSolveIters_total += linearSolveIters_last;
        haveData=true;
      }
      if (status.extraParameters->isParameter("AztecOO/Achieved Tolerance"))
        linearSolveAchievedTol = status.extraParameters->get<double>("AztecOO/Achieved Tolerance");

      if (haveData) {
        *out <<  "Precon Linear Solve ";
        if (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED)
         {*out <<  "Succeeded: "; ++linearSolveSuccessCount;}
        else  *out <<  "Failed: ";
        *out << std::setprecision(3) 
             << linearSolveAchievedTol << " drop in " 
             << linearSolveIters_last << " its (avg: " 
             << linearSolveIters_total / (double) linearSolveCount << " its/slv, " 
             << 100.0* linearSolveSuccessCount / (double) linearSolveCount << "% success)"
             << endl;
      }
    }
  }

  /* Debugging utility to check if columns have been Inserted into the 
   * matrix that do not correspond to a row on any processor
   */
  void check_for_rogue_columns( Epetra_CrsMatrix& mat) {
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
      *out << "ERROR: Column map has " << nrm 
           << " rogue entries that are not associated with any row." << endl;
       rowMap.Comm().Barrier();
       exit(-3);
    }
  }

  //============================================================

} // extern"C"
