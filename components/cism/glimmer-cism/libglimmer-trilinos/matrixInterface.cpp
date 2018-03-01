//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                             
//   matrixInterface.cpp - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
#include "Teuchos_TestForException.hpp"
#include "matrixInterface.hpp"

// Constructor
TrilinosMatrix_Interface::TrilinosMatrix_Interface
  (const Teuchos::RCP<const Epetra_Map>& rowMap,
   int bandwidth, const Epetra_Comm& comm)
  : rowMap_(rowMap), bandwidth_(bandwidth), matrixOrder_(-1), comm_(comm) {
  
  matrixOrder_ = rowMap->NumGlobalElements();

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, bandwidth) );
  isFillCompleted_ = false;
}

// Destructor
TrilinosMatrix_Interface::~TrilinosMatrix_Interface() {
}

// Accessor methods
bool TrilinosMatrix_Interface::isSparsitySet() const {return isFillCompleted_;}
int TrilinosMatrix_Interface::bandwidth() const {return bandwidth_;}
int TrilinosMatrix_Interface::matrixOrder() const {return matrixOrder_;}
const Epetra_Map& TrilinosMatrix_Interface::getRowMap() const {return *rowMap_;}
Teuchos::RCP<Epetra_CrsMatrix>& TrilinosMatrix_Interface::getOperator() {return operator_;}


// Fix the sparsity patter by calling FillComplete
void TrilinosMatrix_Interface::finalizeSparsity() {
  isFillCompleted_ = true;
  int ierr = operator_->FillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
     "Error: Trilinos Fill Complete  returned nozero error code ( " << ierr << " )\n");

}

// Update the operator and also the corresponding row map.
void TrilinosMatrix_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  isFillCompleted_ = operator_->Filled();
}
