// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "noxlocainterface.hpp"

#include "EpetraExt_RowMatrixOut.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(int* nelems, double* statevector,
         const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
         void* blackbox_res_, void* blackbox_prec_,
         void (*residualFunction_)(double *, double *, int, void *),
//         void (*precFunction_)(double *, double *, int, double*, void *)) :
         void (*precFunction_)(double *, double *, int, double*, void *, void *)) :
  N(*nelems),
  comm(comm_),
  pVector(pVector_),
  blackbox_res(blackbox_res_),
  blackbox_prec(blackbox_prec_),
  residualFunction(residualFunction_),
  precFunction(precFunction_)
{ 

  globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

  solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));
}

Problem_Interface::~Problem_Interface()
{ 
}

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  F.PutScalar(0.0);
//  static int icount=0; cout << "Residual called" << icount++ <<  endl;
  residualFunction(x.Values(), F.Values(), N, blackbox_res);

  return true;
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  pVector = params;
}

void Problem_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
//    cout << setprecision(5) << "Solution at: " << conParam << "  is:  " 
//    << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " "
//    << x[4] << " " << x[5] << " ... " << x[N-1] << endl;
}

Teuchos::RCP<Epetra_Vector> Problem_Interface::getVector() const
{ return solution;}

// Preconditioner is 2 steps. Only computePreconditioner is given the state,
//  which we store in solution.
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                                              Epetra_Operator& Prec,
                                              Teuchos::ParameterList* p)
{
//cout << "XXX Inside Problem_Interface::computePreconditioner. Paramlist is:" << endl; p->print(cout);

  (*solution) = x;
  return true;
}

// ApplyInverse is the action of the preconditioner. The state was already sent above.
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

double n2; X(0)->Norm2(&n2);
//cout << "Before: Norm of z="<<n2<<endl;

//   static int icount=0; cout << "Preconditioner called" << icount++ <<  endl;
   precFunction(Y(0)->Values(), X(0)->Values(), N, solution->Values(), 
                                             blackbox_res, blackbox_prec);
// SI preconditioner calls derived type with all the cg solve stuff
//   precFunction(Y(0)->Values(), X(0)->Values(), N, solution->Values(), blackbox_prec);

double n1; Y(0)->Norm2(&n1); 
//cout << "After: Norm of vv="<<n1<<", Norm of z="<<n2<<endl;

   return 0;
}

//-----------------------------------------------------------------------------

