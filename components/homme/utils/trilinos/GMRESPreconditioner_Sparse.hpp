#ifndef GMRESPREC_SPARSE_HPP
#define GMRESPREC_SPARSE_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

#include "SparsePreconditionerBase.hpp"

class GMRESPreconditioner_Sparse : public SparsePreconditionerBase
{
public:
  GMRESPreconditioner_Sparse(const int np_, const int nlev_, const int nelemd_,
			     RCP<Epetra_Vector> xVec, 
			     RCP<Epetra_Map> xMap,
			     void* StateData,
			     void (*precUpdateFunction)(double *,int,void *),
			     const RCP<ParameterList>&  JSolvePL_,
			     int* JTotalIters_);

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_);

private:
  bool firstfill;

  // preconditioner update
  void (*precUpdateFunction)(double *,int,void *);

  // Local row maps for vectors/matrices, not 1-to-1
  RCP<Epetra_Map> StateMap_local;

  // Global row maps for vectors/matrices, 1-to-1
  RCP<Epetra_Map> StateMap;

  // Exporters to go from local to global
  RCP<Epetra_Export> State_Exporter;

  // Importers to go from global to local
  RCP<Epetra_Import> State_Importer;

  // Local vectors
  RCP<Epetra_Vector> StateVec_local;

  // Global vectors
  RCP<Epetra_Vector> StateVec;

  // Jacobian blocks
  RCP<Epetra_FECrsMatrix> Jmatrix;
  
  RCP<Epetra_CrsGraph> Jgraph;
  
  // temporary vectors
  RCP<Epetra_Vector> Jx;

  RCP< Belos::LinearProblem<ST,MV,OP> > JProblem;
  RCP< Belos::SolverManager<ST,MV,OP> > JSolver;
  RCP< ParameterList > JSolvePL;
  int* JTotalIters;

  // Timers
  RCP<Teuchos::Time> ComputePCTime;

  RCP<Teuchos::Time> UpdateStateTime;
  RCP<Teuchos::Time> ZeroTime;
  RCP<Teuchos::Time> FillTime;
  RCP<Teuchos::Time> GlobalAssemblyTime;

  RCP<Teuchos::Time> BuildRHSTime;
  RCP<Teuchos::Time> BuildYvecTime;
  RCP<Teuchos::Time> ApplyTime;
  RCP<Teuchos::Time> JSolveTime;
};
#endif
