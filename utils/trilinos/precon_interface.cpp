// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "precon_interface.hpp"

#include "EpetraExt_RowMatrixOut.h"

//-----------------------------------------------------------------------------
Precon_Interface::Precon_Interface(int* nelems, double* statevector, double* rhs,
         const LOCA::ParameterVector& pVector_p_, const Epetra_Comm& comm_,
         void* blackbox_res_, void* blackbox_prec_,
         void (*residualFunction_lin_)(double *, double *, int, void *),
         void (*precFunction_)(double *, double *, int, double*, void *, void *)) :
//         void (*residualFunction_nd_)(double *, double *, int, void *),
//         void (*precFunction_gmres_)(double *, double *, int, double*, void *, void *)) :
  N(*nelems),
  comm(comm_),
  pVector_p(pVector_p_),
  blackbox_res(blackbox_res_),
  blackbox_prec(blackbox_prec_),
  residualFunction_lin(residualFunction_lin_),
  precFunction(precFunction_)
//  residualFunction_nd(residualFunction_nd_),
//  precFunction_gmres(precFunction_gmres_)
{ 

  globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

  solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));
  solution_rhs = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, rhs));

}

Precon_Interface::~Precon_Interface()
{ 
}

bool Precon_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  F.PutScalar(0.0);
//  static int icount=0; cout << "Residual called" << icount++ <<  endl;
  residualFunction_lin(x.Values(), F.Values(), N, blackbox_res);
//  residualFunction_nd(x.Values(), F.Values(), N, blackbox_res);

  return true;
}

void Precon_Interface::setParameters(const LOCA::ParameterVector& params)
{
  pVector_p = params;
}

void Precon_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
//    cout << setprecision(5) << "Solution at: " << conParam << "  is:  " 
//    << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " "
//    << x[4] << " " << x[5] << " ... " << x[N-1] << endl;
}

Teuchos::RCP<Epetra_Vector> Precon_Interface::getVector_p() const
{ return solution;}

Teuchos::RCP<Epetra_Vector> Precon_Interface::getVector_rhs() const
{ return solution_rhs;}

// Preconditioner is 2 steps. Only computePreconditioner is given the state,
//  which we store in solution.
bool Precon_Interface::computePreconditioner(const Epetra_Vector& x,
                                              Epetra_Operator& Prec,
                                              Teuchos::ParameterList* p)
{
//cout << "XXX Inside Precon_Interface::computePreconditioner. Paramlist is:" << endl; p->print(cout);

  (*solution) = x;
  return true;
}

// ApplyInverse is the action of the preconditioner. The state was already sent above.
int Precon_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

double n2; X(0)->Norm2(&n2);

//   static int icount=0; cout << "Preconditioner called" << icount++ <<  endl;
   precFunction(Y(0)->Values(), X(0)->Values(), N, solution->Values(), 
                                             blackbox_res, blackbox_prec);
//   precFunction_gmres(Y(0)->Values(), X(0)->Values(), N, solution->Values(), 
//                                             blackbox_res, blackbox_prec);
// SI preconditioner calls derived type with all the cg solve stuff
//   precFunction(Y(0)->Values(), X(0)->Values(), N, solution->Values(), blackbox_prec);

double n1; Y(0)->Norm2(&n1); 
//cout << "After: Norm of vv="<<n1<<", Norm of z="<<n2<<endl;

   return 0;
}

//-----------------------------------------------------------------------------

