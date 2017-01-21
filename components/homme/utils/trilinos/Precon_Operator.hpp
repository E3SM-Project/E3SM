#ifndef PRECON_OPERATOR_HPP
#define PRECON_OPERATOR_HPP

#include "Epetra_Operator.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include <iostream>

using std::cout;
using std::endl;

using Teuchos::RCP;


class  Precon_Operator : public Epetra_Operator
{
public:
  Precon_Operator(int nelem_, 
		  RCP<Epetra_Map> globalMap_, 
		  const Epetra_Comm& comm_,
		  void* StateData_,
		  void (*precFunction_)(double *,int,double*,void *));

  ~Precon_Operator();

  // 10 Methods for inheritance from Epetra_Operator
  // Only ApplyInverse is non-trivial -- first 9 satisfied here in header

  int SetUseTranspose(bool UseTranspose)
  { std::cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<endl; return -1;};
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
  { std::cerr<<"ERROR: ApplyInverse"<<endl; return -1;};
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  double NormInf() const { std::cerr<<"norminf"<<endl; return 1.0;};
  const char* Label() const { return "noxlocainterface::user preconditioner";};
  bool UseTranspose() const { return false;};
  bool HasNormInf() const { return false;};
  const Epetra_Comm& Comm() const {return comm;};
  const Epetra_Map& OperatorDomainMap() const {cout<<"returning domain map"<<endl; return *globalMap;};
  const Epetra_Map& OperatorRangeMap() const {cout<<"returning range map"<<endl; return *globalMap;};

private:
  bool printproc;

  int nelem;

  RCP<Epetra_Map> globalMap;

  const Epetra_Comm& comm;

  void* StateData;

  void (*precFunction)(double*, int, double*, void *);
};

#endif

