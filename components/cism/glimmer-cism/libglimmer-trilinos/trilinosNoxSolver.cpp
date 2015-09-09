//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                             
//   trilinosNoxSolver.cpp - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

// Trilinos Objects
#include "Piro_Epetra_NOXSolver.hpp"
#include "Piro_Epetra_LOCASolver.hpp"
#include "trilinosModelEvaluator.hpp"

#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

#include "config.inc"

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
// Use continuation instead of straight Newton for this many time steps:

void setCismLocaDefaults(Teuchos::ParameterList& locaList) {
  Teuchos::ParameterList& predList = locaList.sublist("Predictor");
  Teuchos::ParameterList& stepperList = locaList.sublist("Stepper");
  Teuchos::ParameterList& stepSizeList = locaList.sublist("Step Size");

  // If not set in XML list, set these defaults instead
  (void) predList.get("Method","Constant");
  (void) stepperList.get("Continuation Method","Natural");
  (void) stepperList.get("Continuation Parameter","Effstrmin Factor");
  (void) stepperList.get("Initial Value",10.0);
  (void) stepperList.get("Max Steps",10);
  (void) stepperList.get("Max Value",100.0); // not used
  (void) stepperList.get("Min Value",0.0); // Important!!

  (void) stepSizeList.get("Initial Step Size",-3.0); // Important!!
  (void) stepSizeList.get("Aggressiveness",2.0); // Important!!
}


extern "C" {
void FC_FUNC(noxinit,NOXINIT) ( int* nelems, double* statevector,
               int* mpi_comm_f, void* blackbox_res)
// mpi_comm_f: CISM's fortran mpi communicator
{

 bool succeeded=true;
 try {

  // Build the epetra communicator
  MPI_Comm mpi_comm_c = MPI_Comm_f2c(*mpi_comm_f);
  Comm_=rcp(new Epetra_MpiComm(mpi_comm_c));
  Epetra_Comm& Comm=*Comm_;
  printProc = (Comm_->MyPID() == 0);
  Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(mpi_comm_c));
  
  if (printProc) cout << "NOXINIT CALLED    for nelem=" << *nelems << endl;

    try { // Check that the parameter list is valid at the top
      RCP<Teuchos::ParameterList> pl =
        rcp(new Teuchos::ParameterList("Trilinos Options for NOX"));
      Teuchos::updateParametersFromXmlFileAndBroadcast(
                             "trilinosOptions.xml", pl.ptr(),tcomm);
 
      Teuchos::ParameterList validPL("Valid List");;
      validPL.sublist("Stratimikos"); validPL.sublist("Piro");
      pl->validateParameters(validPL, 0);
      paramList = Teuchos::sublist(pl,"Piro",true);
    }
    catch (std::exception& e) {
      cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"
           << e.what() << "\nExiting: Invalid trilinosOptions.xml file."
           << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
      exit(1);
    }
  
  paramList->set("Lean Matrix Free",true); // Saves some GMRES steps
  if (printProc) cout << "NOXInit: param list is: (delete this debug line)\n" << *paramList << endl;

  model = rcp(new trilinosModelEvaluator(*nelems, statevector, Comm, blackbox_res));
    
  // Logic to see if we want to use LOCA continuation or NOX single steady solve
  // Turn on LOCA by having a LOCA sublist OR setting "CISM: Number of Time Steps To Use LOCA"
  bool useLoca=false;
  // If LOCA sublist exists, defaults to using it for 1 time step; but can be set in XML.
  int numStepsToUseLOCA = 0;
  if (paramList->isSublist("LOCA"))
    numStepsToUseLOCA = paramList->get("CISM: Number of Time Steps To Use LOCA",1);
  else 
    numStepsToUseLOCA = paramList->get("CISM: Number of Time Steps To Use LOCA",0);

  if (timeStep <= numStepsToUseLOCA) useLoca=true;

  if (useLoca) if (printProc)
    cout << "\nUsing LOCA continuation for first " << numStepsToUseLOCA << "  time steps." << endl;

  if (useLoca) {
    setCismLocaDefaults(paramList->sublist("LOCA"));
    Nsolver = rcp(new Piro::Epetra::LOCASolver(paramList, model));
  }
  else
    Nsolver = rcp(new Piro::Epetra::NOXSolver(paramList, model));

  inArgs=Nsolver->createInArgs();
  outArgs=Nsolver->createOutArgs();

  // Ask the model for the converged solution from g(0)
  RCP<const Epetra_Map> xmap = Nsolver->get_g_map(0);
  RCP<Epetra_Vector> xout = rcp(new Epetra_Vector(*xmap));

  outArgs.set_g(0,xout);

  // Set up parameter vector for continuation runs
  if (useLoca) {
    RCP<const Epetra_Map> pmap = Nsolver->get_p_map(0);
    RCP<Epetra_Vector> pvec = rcp(new Epetra_Vector(*pmap));
    inArgs.set_p(0, pvec);
  }

  // Time step counter: just for deciding whether to use continuation on relaxatin param
  timeStep++;

 } //end try block
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);
}

/****************************************************/
void FC_FUNC(noxsolve,NOXSOLVE) (int* nelems, double* statevector, void* blackbox_res)
{
  bool succeeded=true;
  try {
    TEUCHOS_TEST_FOR_EXCEPTION(Nsolver==Teuchos::null, logic_error, 
                          "Exception: noxsolve called with solver=null: \n"
       << "You either did not call noxinit first, or called noxfinish already");
    if (printProc) cout << "NOXSolve called" << endl;

    // Solve    
    Nsolver->evalModel(inArgs,outArgs);

    // Copy out the solution
    RCP<Epetra_Vector> xout = outArgs.get_g(0); 
    if(xout == Teuchos::null) throw "evalModel is NOT returning a vector";

    for (int i=0; i<*nelems; i++) statevector[i] = (*xout)[i];
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, succeeded);
  if (!succeeded) exit(1);

}

/****************************************************/ 
void FC_FUNC(noxfinish,NOXFINISH) (void)
{
 if (printProc) cout << "NOXFinish called" << endl;

 // Free memory
 Nsolver   = Teuchos::null;
 model     = Teuchos::null;
 paramList = Teuchos::null;
 Comm_     = Teuchos::null;
}

} //extern "C"
