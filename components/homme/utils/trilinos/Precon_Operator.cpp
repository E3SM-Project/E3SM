#include "Precon_Operator.hpp"

using Teuchos::RCP;


Precon_Operator::Precon_Operator(int nelem_, 
				 RCP<Epetra_Map> globalMap_, 
				 const Epetra_Comm& comm_, 
				 void* StateData_, 
				 void (*precFunction_)(double*, int, double*, void*))
  : nelem(nelem_),
    globalMap(globalMap_),
    comm(comm_),
    StateData(StateData_),
    precFunction(precFunction_)
{ 
  printproc = (comm.MyPID()==0);
}


Precon_Operator::~Precon_Operator() { }


int Precon_Operator::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y)const
{
  Y.PutScalar(0.0);

  // Apply preconditioner function P to X and store result in X, StateData is the necessary p.c. data
  precFunction(X(0)->Values(), nelem, Y(0)->Values(), StateData);

  // Uncomment the line below to make an identity preconditioner 
  //Y=X;

  return 0;
}
