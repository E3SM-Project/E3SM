// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "block_precon_interface.hpp"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"




#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"






//-----------------------------------------------------------------------------
//
Block_Precon_Interface::Block_Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_, void* precdata_, void (*precFunction_)(double *,int,double*,void *)):
	N(nelems),
	globalMap(gmap),
 	comm(comm_),
	precdata(precdata_),
	precFunction(precFunction_)
{ 
  if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
//if(printproc)std::cout<<"mypid_precon="<<comm.MyPID()<<std::endl;
}

Block_Precon_Interface::Block_Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_, void* precdata_, 
void (*precFunctionblock11_)(double *,int,double*,void *),
void (*precFunctionblock12_)(double *,int,double*,void *),
void (*precFunctionblock21_)(double *,int,double*,void *),
void (*precFunctionblock22_)(double *,int,double*,void *)):
	N(nelems),
	globalMap(gmap),
 	comm(comm_),
	precdata(precdata_),
	precFunctionblock11(precFunctionblock11_),
	precFunctionblock12(precFunctionblock12_),
	precFunctionblock21(precFunctionblock21_),
	precFunctionblock22(precFunctionblock22_)
{ 
  if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
//if(printproc)std::cout<<"mypid_precon="<<comm.MyPID()<<std::endl;
}


Block_Precon_Interface::~Block_Precon_Interface() { }

/*
int Precon_Interface::Apply(const Epetra_MultiVector &X,Epetra_MultiVector &Y)const
{
	//	int n=X.MyLength();

       Y.PutScalar(0.0);

//if(printproc)std::cout<<"mypid_precon="<<comm.MyPID()<<std::endl;


//double n0; X(0)->Norm2(&n0);
//if(printproc) std::cout << "ApplyPrec Norm of x in="<<n0<<std::endl;

// std::cout<<"Applying Preconditioning Operator"<<std::endl;
//if (printproc) std::cout << "PX address in:  " << &X << "   PY address in:  " << &Y << std::endl;
//Y=X;
//if (printproc) std::cout << "PX address out:  " << &X << "   PY address out:  " << &Y << std::endl;

//double n1; Y(0)->Norm2(&n1);
//if(printproc) std::cout << "ApplyPrec Norm of y="<<n1<<std::endl;

//double n2; X(0)->Norm2(&n2);
//if(printproc) std::cout << "ApplyPrec Norm of x out="<<n2<<std::endl;

//       Y.PutScalar(0.0);



	//	Epetra_MultiVector tmp=Y;
	// Apply preconditioner function P to X and store result in X, prec_data is the necessary p.c. data
			 precFunction(X(0)->Values(),N,Y(0)->Values(),precdata);

//Y=X;

	return 0;
}
*/
int Block_Precon_Interface::Apply(const Epetra_MultiVector &X,Epetra_MultiVector &Y)const
{
        //      int n=X.MyLength();

       Y.PutScalar(0.0);
Epetra_Vector P11(*Y(0));
Epetra_Vector P12(*Y(0));
Epetra_Vector P21(*Y(0));
Epetra_Vector P22(*Y(0));


P11.PutScalar(0.0);
P12.PutScalar(0.0);
P21.PutScalar(0.0);
P22.PutScalar(0.0);

        // Apply preconditioner function P to X and store result in X, prec_data is the necessary p.c. data
      precFunctionblock11(X(0)->Values(),N,P11.Values(),precdata);



//check norm of x-p11 to see the effect of the mass matrix

//then create p22 that does mass inverse and check norm of x-p22 to see if we get zero 



      precFunctionblock12(X(0)->Values(),N,P12.Values(),precdata);
      precFunctionblock21(X(0)->Values(),N,P21.Values(),precdata);
      precFunctionblock22(X(0)->Values(),N,P22.Values(),precdata);

Epetra_Vector Yt(*Y(0));
Yt.PutScalar(0.0);


// Epetra_Vector Ix1(*X(0));
 //        for (int i=2*(N+1)/3;i<N; i++) {Ix1[i]=0.0;}





    //for (int i=0; i<N; i++) Yt[i] = P11[i]+P21[i]+P22[i];
   // for (int i=0; i<N; i++) Yt[i] = P11[i]+P12[i]+P21[i];// Saddle point
    for (int i=0; i<N; i++) Yt[i] = P11[i]+P12[i]+P21[i]+P22[i];
    //for (int i=0; i<N; i++) Yt[i] = P11[i]+P22[i]; //Block Diagonal Precon
  // for (int i=0; i<N; i++) Yt[i] = P11[i]+P12[i]+P22[i]; //DU
   // for (int i=0; i<N; i++) Yt[i] = P11[i]+P21[i]+P22[i]; //LD
    //for (int i=0; i<N; i++) Yt[i] = Ix1[i]+P22[i];
Y=Yt;


/* ... or Apply Identity preonditioner */
//Y=X;
        return 0;
}







//-----------------------------------------------------------------------------

