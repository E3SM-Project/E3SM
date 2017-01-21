#ifndef HOMME_BLOCK_PREC_DENSE_HPP
#define HOMME_BLOCK_PREC_DENSE_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

class BlockPreconditioner_Dense : public hommePreconditionerBase 
{
public:
  BlockPreconditioner_Dense(int np_, int nlev_, int nelemd_, 
                            RCP<Epetra_Vector> xVec, RCP<Epetra_Map> xMap,
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
  bool printproc; // output process flag

  int nState; // length of full state vector

  int np;     // number of nodes along an element edge
  int nlev;   // number of vertical levels
  int nelemd; // number of elements on this process

  int nVel;   // length of velocity vector
  int nPs;    // length of surface pressure vector
  int nT;     // length of temperature vector

  // residual data
  void* StateData;

  // preconditioner update 
  void (*precUpdateFunction)(double *,int,void *);

  // Fortran functions for preconditioner block operators
  void (*Solve1)(double*, int, double*, void*);
  void (*Solve2)(double*, int, double*, void*);
  void (*Solve3)(double*, int, double*, void*);

  // vector maps
  RCP<Epetra_Map> VelMap;
  RCP<Epetra_Map> PsMap;
  RCP<Epetra_Map> TMap;

  // temporary work vectors
  RCP<Epetra_Vector> w1;
  RCP<Epetra_Vector> w2;
  RCP<Epetra_Vector> w3;

  // A matrix solve
  RCP<Epetra_Operator> Aop;
  RCP<Epetra_Vector> Ab; 
  RCP<Epetra_Vector> Ax; 

  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > AProblem;
  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > ASolver;
  Teuchos::RCP< ParameterList >  ASolvePL;
  int* ATotalIters;

  // S Schur matrix solve
  RCP<Epetra_Operator> Sop;
  RCP<Epetra_Vector> Sb; 
  RCP<Epetra_Vector> Sx; 

  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SProblem;
  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > SSolver;
  Teuchos::RCP< ParameterList >  SSolvePL;
  int* STotalIters; 

  // P Schur matrix solve
  RCP<Epetra_Operator> Pop;
  RCP<Epetra_Vector> Pb; 
  RCP<Epetra_Vector> Px; 

  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > PProblem;
  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > PSolver;
  Teuchos::RCP<ParameterList>  PSolvePL;
  int* PTotalIters; 

  // Timers
  RCP<Teuchos::Time> AssemblyTime;
  RCP<Teuchos::Time> ApplyTime;
  RCP<Teuchos::Time> ASolveTime;
  RCP<Teuchos::Time> SSolveTime;
  RCP<Teuchos::Time> PSolveTime;
};

#endif
