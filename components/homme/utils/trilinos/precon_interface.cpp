// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "precon_interface.hpp"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"




#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"



//-----------------------------------------------------------------------------
//
Precon_Interface::Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_, void* precdata_, void (*precFunction_)(double *,int,double*,void *)):
	N(nelems),
	globalMap(gmap),
 	comm(comm_),
	precdata(precdata_),
	precFunction(precFunction_)
{ 
  if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
}



Precon_Interface::~Precon_Interface() { }

int Precon_Interface::Apply(const Epetra_MultiVector &X,Epetra_MultiVector &Y)const
{

       Y.PutScalar(0.0);

	// Apply preconditioner function P to X and store result in X, prec_data is the necessary p.c. data
			 precFunction(X(0)->Values(),N,Y(0)->Values(),precdata);
        // Uncomment the line below to make an identity preconditioner 
	//Y=X;

	return 0;
}






//-----------------------------------------------------------------------------

