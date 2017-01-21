#ifndef HOMME_IDENTITY_PREC_HPP
#define HOMME_IDENTITY_PREC_HPP

#include "trilinosModelEvaluator.hpp" // include preconditioner base class

class identityPreconditioner : public hommePreconditionerBase 
{
public:
  identityPreconditioner(int nState_,
			 RCP<Epetra_Vector> xVec_,
			 RCP<Epetra_Map> xMap_);

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int computePreconditioner(RCP<const Epetra_Vector> xVec,  void* StateData_);

private:
  // flag for output processor
  bool printproc;

  // local state vector length
  int nState;
};

#endif
