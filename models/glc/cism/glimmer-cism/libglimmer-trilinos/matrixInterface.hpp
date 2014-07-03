#ifndef TRILINOSMATIX_INTERFACE_H
#define TRILINOSMATIX_INTERFACE_H

#include <iostream>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#ifdef GLIMMER_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

class TrilinosMatrix_Interface {
public:
  // Constructor
  TrilinosMatrix_Interface(const Teuchos::RCP<const Epetra_Map>& rowMap,
                           int bandwidth, const Epetra_Comm& comm);

  // Destructor
  ~TrilinosMatrix_Interface();

  // Accessors
  bool isSparsitySet() const;
  int bandwidth() const;
  int matrixOrder() const;
  const Epetra_Map& getRowMap() const;
  Teuchos::RCP<Epetra_CrsMatrix>& getOperator();

  // Mutators
  void finalizeSparsity(); // Call FillComplet to lock in sparsity pattern
  void updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator);

private:
  bool isFillCompleted_; // to indicate if operator_ is "FillComplete()"ed
  int bandwidth_;
  int matrixOrder_;
  const Epetra_Comm& comm_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<const Epetra_Map> rowMap_;
};
#endif
