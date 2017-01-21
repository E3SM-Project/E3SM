#include "IdentityPreconditioner.hpp"

//#define DEBUG_PRINT_ON

using std::cout;
using std::endl;

// ------------------------------------------------------------------------------
// constructor for identity preconditioner
// ------------------------------------------------------------------------------
identityPreconditioner::identityPreconditioner(int nState_, 
					       RCP<Epetra_Vector> xVec_, 
					       RCP<Epetra_Map> xMap_)
  : hommePreconditionerBase(xMap_), //Required Base Class construction
    nState(nState_)
{
  printproc = (xVec_->Comm().MyPID() == 0);
  if (printproc) cout << endl << ">>> IP: Constructing Identity Preconditioner <<<" << endl << endl;
}


// ------------------------------------------------------------------------------
// update preconditioner, does nothing for identity preconditioner
// ------------------------------------------------------------------------------
int identityPreconditioner::computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_)
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> IP: Computing Identity Preconditioner <<<" << endl;
#endif

  // do nothing for identity preconditioner
  return 0;
}


// ------------------------------------------------------------------------------
// apply inverse of preconditioner matrix
// ------------------------------------------------------------------------------
int identityPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
#ifdef DEBUG_PRINT_ON      
  if (printproc) cout << ">>> IP: Applying Inverse of Identity Preconditioner <<<" << endl;
#endif

  Y = X;

  return 0;
}
