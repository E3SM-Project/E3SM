#ifndef HOMME_GMRES_PREC_OP_HPP
#define HOMME_GMRES_PREC_OP_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

class GMRESPreconditioner_Op : public hommePreconditionerBase 
{
public:
  GMRESPreconditioner_Op(int np_, int nlev_, int nelemd_, 
                         RCP<Epetra_Vector> xVec, RCP<Epetra_Map> xMap,
                         void* StateData,
                         void (*precUpdateFunction)(double *,int,void *),
                         const RCP<ParameterList>&  JSolvePL_,
                         int* JTotalIters_);
  
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
  void (*precUpdateFunction)(double*, int, void*);

  // function for preconditioner blocks
  void (*Solve1)(double*, int, double*, void*);

  // vector maps
  RCP<Epetra_Map> StateMap;

  // J matrix solve
  RCP<Epetra_Operator> Jop;
  RCP<Epetra_Vector> Jb; 
  RCP<Epetra_Vector> Jx; 

  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > JProblem;
  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > JSolver;
  Teuchos::RCP< ParameterList >  JSolvePL;
  int* JTotalIters;

  // Timers
  RCP<Teuchos::Time> AssemblyTime;
  RCP<Teuchos::Time> ApplyTime;
  RCP<Teuchos::Time> JSolveTime;
};

#endif
