#ifndef BLOCKPREC_SPARSE_HPP
#define BLOCKPREC_SPARSE_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "EpetraExt_MatrixMatrix.h"   // Matrix-Matrix operations
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

#include "SparsePreconditionerBase.hpp"

class BlockPreconditioner_Sparse : public SparsePreconditionerBase
{
public:
  BlockPreconditioner_Sparse(const int np_, 
                             const int nlev_,
                             const int nelemd_, 
                             RCP<Epetra_Vector> xVec, 
                             RCP<Epetra_Map> xMap,
                             void* StateData,
                             void (*precUpdateFunction)(double *,int,void *),
                             const RCP<ParameterList>&  ASolvePL_,
                             const RCP<ParameterList>&  SSolvePL_,
                             const RCP<ParameterList>&  PSolvePL_,
                             int* ATotalIters_,
                             int* STotalIters_,
                             int* PTotalIters_);
  
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  int computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_);
  
private:
  bool firstfill;

  // preconditioner update 
  void (*precUpdateFunction)(double *,int,void *);

  // Local row maps for vectors/matrices, not 1-to-1
  RCP<Epetra_Map> VelMap_local;
  RCP<Epetra_Map> PsMap_local;
  RCP<Epetra_Map> TMap_local;

  // Global row maps for vectors/matrices, 1-to-1
  RCP<Epetra_Map> VelMap;
  RCP<Epetra_Map> PsMap;
  RCP<Epetra_Map> TMap;

  // Exporters to go from local to global
  RCP<Epetra_Export> Vel_Exporter;
  RCP<Epetra_Export> Ps_Exporter;
  RCP<Epetra_Export> T_Exporter;
  
  // Importers to go from global to local
  RCP<Epetra_Import> Vel_Importer;
  RCP<Epetra_Import> Ps_Importer;
  RCP<Epetra_Import> T_Importer;

  // Local vectors 
  RCP<Epetra_Vector> VelVec_local;
  RCP<Epetra_Vector> PsVec_local;
  RCP<Epetra_Vector> TVec_local;

  // Global vectors
  RCP<Epetra_Vector> VelVec;
  RCP<Epetra_Vector> PsVec;
  RCP<Epetra_Vector> TVec;

  // Jacobian blocks 
  RCP<Epetra_FECrsMatrix> Ablock;
  RCP<Epetra_FECrsMatrix> Bblock;
  RCP<Epetra_FECrsMatrix> Cblock;

  RCP<Epetra_FECrsMatrix> Dblock;
  RCP<Epetra_FECrsMatrix> Eblock;

  RCP<Epetra_FECrsMatrix> Fblock;
  RCP<Epetra_FECrsMatrix> Gblock;
  RCP<Epetra_FECrsMatrix> Hblock;

  // Diagonal approximations of A^{-1}
  RCP<Epetra_Vector> Ahatinv; 
  
  RCP<Epetra_FECrsMatrix> Wmatrix; // W = D hat{A}^{-1} * B

  // Schur Complement (Note Crs Matrix to work with Add)
  RCP<Epetra_CrsMatrix> Smatrix;

  // temporary vectors
  RCP<Epetra_Vector> DAx;
  RCP<Epetra_Vector> w1;

  // A matrix solve
  RCP<Epetra_Vector> Ax; 

  RCP< Belos::LinearProblem<ST,MV,OP> > AProblem;
  RCP< Belos::SolverManager<ST,MV,OP> > ASolver;
  RCP< ParameterList > ASolvePL;
  int* ATotalIters;

  // S Schur matrix solve
  RCP<Epetra_Vector> Sx; 

  RCP< Belos::LinearProblem<ST,MV,OP> > SProblem;
  RCP< Belos::SolverManager<ST,MV,OP> > SSolver;
  RCP< ParameterList >  SSolvePL;
  int* STotalIters; 

  // P Schur matrix solve
  RCP<Epetra_Vector> Px; 

  RCP< Belos::LinearProblem<ST,MV,OP> > PProblem;
  RCP< Belos::SolverManager<ST,MV,OP> > PSolver;
  RCP<ParameterList>  PSolvePL;
  int* PTotalIters; 

  // Timers
  RCP<Teuchos::Time> ComputePCTime;

  RCP<Teuchos::Time> UpdateStateTime;

  RCP<Teuchos::Time> ABlockAssemblyTime;
  RCP<Teuchos::Time> BBlockAssemblyTime;
  RCP<Teuchos::Time> CBlockAssemblyTime;
  RCP<Teuchos::Time> DBlockAssemblyTime;
  RCP<Teuchos::Time> EBlockAssemblyTime;
  RCP<Teuchos::Time> FBlockAssemblyTime;
  RCP<Teuchos::Time> GBlockAssemblyTime;
  RCP<Teuchos::Time> HBlockAssemblyTime;

  RCP<Teuchos::Time> BuildAhatInv;
  RCP<Teuchos::Time> BuildAhatInvB;
  RCP<Teuchos::Time> BuildAhatInvC;
  RCP<Teuchos::Time> BuildDAhatInvB; // W = DAB
  RCP<Teuchos::Time> BuildDAhatInvC; // X = DAC
  RCP<Teuchos::Time> BuildFAhatInvB; // Y = FAB
  RCP<Teuchos::Time> BuildFAhatInvC; // Y = FAC
  RCP<Teuchos::Time> BuildSchurS;
  RCP<Teuchos::Time> BuildShatInv;
  RCP<Teuchos::Time> BuildX; // X
  RCP<Teuchos::Time> BuildR;
  RCP<Teuchos::Time> BuildRX;
  RCP<Teuchos::Time> BuildSchurP;

  RCP<Teuchos::Time> ApplyTime;

  RCP<Teuchos::Time> BuildRHSTime;
  RCP<Teuchos::Time> BuildYvecTime;

  RCP<Teuchos::Time> ASolveTime;
  RCP<Teuchos::Time> SSolveTime;
  RCP<Teuchos::Time> PSolveTime;
};
#endif
