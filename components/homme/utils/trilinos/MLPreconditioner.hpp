#ifndef MLPREC_HPP
#define MLPREC_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"

#include "SparsePreconditionerBase.hpp"

class MLPreconditioner : public SparsePreconditionerBase
{
public:
  MLPreconditioner(const int np_, const int nlev_, const int nelemd_,
                   RCP<Epetra_Vector> xVec, 
                   RCP<Epetra_Map> xMap,
                   void* StateData,
                   void (*precUpdateFunction)(double *,int,void *),
                   const RCP<ParameterList>&  JSolvePL_,
                   int* JTotalIters_);
    
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_);

private:

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

  // temporary vectors
  RCP<Epetra_Vector> Jx;

  RCP<ML_Epetra::MultiLevelPreconditioner> MLJsolve;
  RCP< ParameterList > JSolvePL;
  int* JTotalIters;

  ML* ml_object;

  // Timers
  RCP<Teuchos::Time> AssemblyTime;
  RCP<Teuchos::Time> BuildRHSTime;
  RCP<Teuchos::Time> BuildYvecTime;
  RCP<Teuchos::Time> ApplyTime;
  RCP<Teuchos::Time> JSolveTime;

  RCP<Teuchos::Time> MLSetupTime;
};
#endif
