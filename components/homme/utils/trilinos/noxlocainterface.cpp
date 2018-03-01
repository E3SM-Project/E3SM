// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "noxlocainterface.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Operator.h"

#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "Epetra_LinearProblem.h"
//#include "AztecOO.h"
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <AztecOO.h>
#include "precon_interface.hpp"
#include "block_precon_interface.hpp"


#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>


#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"




//-----------------------------------------------------------------------------


Problem_Interface::Problem_Interface(int nelems, double* statevector_,
		const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
		void* blackbox_res_, void* precdata_, void* jacdata_,
		void (*residualFunction_)(double *, double *, int, void *),
		void (*precFunction_)(double *, int, double*, void *),
		void (*jacFunction_)(double *, int, double*, void *),
		void (*precUpdateFunction_)(double *, int, void *),
		void (*getJacVector_)(double *, int, void *)) :
	N(nelems),
	statevector(statevector_),
	comm(comm_),
	pVector(pVector_),
	blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
	residualFunction(residualFunction_),
	precFunction(precFunction_),
	jacFunction(jacFunction_),
	precUpdateFunction(precUpdateFunction_),
	getJacVector(getJacVector_)
{ 

	globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

	solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));

        if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
}



Problem_Interface::Problem_Interface(int nelems, double* statevector_,
		const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
		void* blackbox_res_, void* precdata_, void* jacdata_,
		void (*residualFunction_)(double *, double *, int, void *),
		void (*precFunctionblock11_)(double *, int, double*, void *),
		void (*precFunctionblock12_)(double *, int, double*, void *),
		void (*precFunctionblock21_)(double *, int, double*, void *),
		void (*precFunctionblock22_)(double *, int, double*, void *),
		void (*jacFunction_)(double *, int, double*, void *),
		void (*precUpdateFunction_)(double *, int, void *),
		void (*getJacVector_)(double *, int, void *)) :
	N(nelems),
	statevector(statevector_),
	comm(comm_),
	pVector(pVector_),
	blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
	residualFunction(residualFunction_),
	precFunctionblock11(precFunctionblock11_),
	precFunctionblock12(precFunctionblock12_),
	precFunctionblock21(precFunctionblock21_),
	precFunctionblock22(precFunctionblock22_),
	jacFunction(jacFunction_),
	precUpdateFunction(precUpdateFunction_),
	getJacVector(getJacVector_)
{ 

	globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

	solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));

        if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
}

/* This interface is just for testing... comparing two block preconditioner formulations */

Problem_Interface::Problem_Interface(int nelems, double* statevector_,
                const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
                void* blackbox_res_, void* precdata_, void* jacdata_,
                void (*residualFunction_)(double *, double *, int, void *),
                void (*precFunctionblock11_)(double *, int, double*, void *),
                void (*precFunctionblock12_)(double *, int, double*, void *),
                void (*precFunctionblock21_)(double *, int, double*, void *),
                void (*precFunctionblock22_)(double *, int, double*, void *),
                void (*auxprecFunctionblock11_)(double *, int, double*, void *),
                void (*auxprecFunctionblock12_)(double *, int, double*, void *),
                void (*auxprecFunctionblock21_)(double *, int, double*, void *),
                void (*auxprecFunctionblock22_)(double *, int, double*, void *),
                void (*jacFunction_)(double *, int, double*, void *),
                void (*precUpdateFunction_)(double *, int, void *),
                void (*getJacVector_)(double *, int, void *)) :
        N(nelems),
        statevector(statevector_),
        comm(comm_),
        pVector(pVector_),
        blackbox_res(blackbox_res_),
        precdata(precdata_),
        jacdata(jacdata_),
        residualFunction(residualFunction_),
        precFunctionblock11(precFunctionblock11_),
        precFunctionblock12(precFunctionblock12_),
        precFunctionblock21(precFunctionblock21_),
        precFunctionblock22(precFunctionblock22_),
        auxprecFunctionblock11(auxprecFunctionblock11_),
        auxprecFunctionblock12(auxprecFunctionblock12_),
        auxprecFunctionblock21(auxprecFunctionblock21_),
        auxprecFunctionblock22(auxprecFunctionblock22_),
        jacFunction(jacFunction_),
        precUpdateFunction(precUpdateFunction_),
        getJacVector(getJacVector_)
{

        globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

        solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));

        if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
}







Problem_Interface::Problem_Interface(int nelems, double* statevector_,
		const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
		void* blackbox_res_, void* precdata_,
		void (*residualFunction_)(double *, double *, int, void *),
		void (*precFunction_)(double *, int, double*, void *),
		void (*precUpdateFunction_)(double *, int, void *)) :
	N(nelems),
	statevector(statevector_),
	comm(comm_),
	pVector(pVector_),
	blackbox_res(blackbox_res_),
        precdata(precdata_),
	residualFunction(residualFunction_),
	precFunction(precFunction_),
	precUpdateFunction(precUpdateFunction_)
{ 

	globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

	solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));

        if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
}

/*
Problem_Interface::Problem_Interface(int nelems, double* statevector_,
		const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
		void* blackbox_res_, void* precdata_,
		void (*residualFunction_)(double *, double *, int, void *),
		void (*precFunction_)(double *, int, double*, void *),
		void (*precUpdateFunction_)(double *, int, void *)) :
	N(nelems),
	statevector(statevector_),
	comm(comm_),
	pVector(pVector_),
	blackbox_res(blackbox_res_),
        precdata(precdata_),
	residualFunction(residualFunction_),
	precFunction(precFunction_),
	precUpdateFunction(precUpdateFunction_)
{ 

	globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

	solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));

        if (comm.MyPID()==0) printproc=true;
        else   printproc=false;
//if(printproc)cout<<"mypid_interface="<<comm.MyPID()<<endl;


}
*/

Problem_Interface::~Problem_Interface()
{ 
}


void Problem_Interface::printCommID(){
        if (printproc)   cout<<"Printing CommID"<<endl<<flush;
        if (comm.MyPID()==0) cout<<"MyPID_interface="<<comm.MyPID()<<endl<<flush;
}



void Problem_Interface::resetBlackbox(void* blackbox_res_,  void* precdata_,void* jacdata_){
	blackbox_res=blackbox_res_; 
	precdata=precdata_;
	jacdata=jacdata_;
}



void Problem_Interface::resetBlackbox(void* blackbox_res_,  void* precdata_){
	blackbox_res=blackbox_res_; 
	precdata=precdata_;
}



void Problem_Interface::resetBlackbox(void* blackbox_res_){
	blackbox_res=blackbox_res_; 
}



bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
cout.precision(16);
//     double n1; x.Norm2(&n1);
//if(printproc) cout << "computeF Norm of x="<<n1<<endl;

	F.PutScalar(0.0);
	//  static int icount=0; cout << "Residual called" << icount++ <<  endl;
	residualFunction(x.Values(), F.Values(), N, blackbox_res);
	//residualFunction(x.Values(), F.Values(), N, jacdata);

//     double n2; F.Norm2(&n2);
//if(printproc) cout << "computeF Norm of F="<<n2<<endl;
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


void Problem_Interface::getJacVec(Epetra_Vector &y)const
{ 
        getJacVector(y.Values(),N,jacdata);
        //getJacVector(y.Values(),N,blackbox_res);

}


void Problem_Interface::getbbVec(Epetra_Vector &y)const
{ 
        getJacVector(y.Values(),N,blackbox_res);

}




bool Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
//        if (printproc)   cout<<"Updating Jacobian data before nonlinear solve"<<endl;
        // cout<<"Updating Jacobian data before nonlinear solve"<<endl;

//precUpdateFunction should be more appropriately named UpdateStateVector
//test 
//Epetra_Vector a(x);
//a.PutScalar(0.0);
//precUpdateFunction(a.Values(),N,jacdata);


/*
cout.precision(16);



Epetra_Vector Jacv(x);
Epetra_Vector F(x); 
Epetra_Vector step(x);
step.PutScalar(0.0);

Epetra_Vector a(x);
a.PutScalar(0.0);
getJacVec(a); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a

     double na1; a.Norm2(&na1);
if(printproc) cout << "computeJacobian norm of na1="<<na1<<endl;

     double n0; x.Norm2(&n0);
if(printproc) cout << "computeJacobian norm of x="<<n0<<endl;


	residualFunction(x.Values(), F.Values(), N, blackbox_res);
     double n1; F.Norm2(&n1);
if(printproc) cout << "computeF Norm of F="<<n1<<endl;


        jacFunction(step.Values(),N,Jacv.Values(),jacdata);

double n2; Jacv.Norm2(&n2);
if(printproc) cout << "computeJacobian Norm of preJac="<<n2<<endl;
*/

//only necessary line
      precUpdateFunction(x.Values(),N,jacdata);
      precUpdateFunction(x.Values(),N,precdata);


/*
Epetra_Vector b(x);
b.PutScalar(0.0);
getJacVec(b); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a
     double nb2; b.Norm2(&nb2);
if(printproc) cout << "computeJacobian norm of nb2="<<nb2<<endl;
	
residualFunction(x.Values(), F.Values(), N, blackbox_res);
     double n3; F.Norm2(&n3);
if(printproc) cout << "computeF Norm of F="<<n3<<endl;

        jacFunction(step.Values(),N,Jacv.Values(),jacdata);
double n4; Jacv.Norm2(&n4);
if(printproc) cout << "computeJacobian Norm of postJac="<<n4<<endl;



//revert back for test
      precUpdateFunction(a.Values(),N,jacdata);
      precUpdateFunction(a.Values(),N,precdata);


Epetra_Vector c(x);
c.PutScalar(0.0);
getJacVec(c); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a

     double nc1; c.Norm2(&nc1);
if(printproc) cout << "computeJacobian norm of nc1="<<nc1<<endl;
	
residualFunction(c.Values(), F.Values(), N, blackbox_res);
     double nc3; F.Norm2(&nc3);
if(printproc) cout << "computeF Norm of F="<<nc3<<endl;

        jacFunction(step.Values(),N,Jacv.Values(),jacdata);
double nc4; Jacv.Norm2(&nc4);
if(printproc) cout << "computeJacobian Norm of postJac="<<nc4<<endl;



        precUpdateFunction(x.Values(),N,jacdata);
        precUpdateFunction(x.Values(),N,precdata);
//        precUpdateFunction(x.Values(),N,precdata);
//        precUpdateFunction(x.Values(),N,blackbox_res);

*/
/*
Epetra_Vector a(x);
a.PutScalar(10.0);
getJacVec(a);


if (printproc) {
cout<<"Jacx=["<<endl;
for(int i=0; i<N;i++) cout<<*(a.Values()+i)<<endl;
cout<<"];"<<endl;
}

*/


//if(printproc) cout << "computejac x="<<x<<endl<<flush;


//        if (printproc)   cout<<x;
        //precUpdateFunction(x.Values(),N,blackbox_res);

	return true;
}


// Preconditioner is 2 steps. Only computePreconditioner is given the state,
//  which we store in solution.
/*
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
		Epetra_Operator& Prec,
		Teuchos::ParameterList* p)
{
//        if (printproc)   cout<<"Updating Preconditioner before nonlinear solve"<<endl;
// precUpdateFunction(x.Values(),N,precdata);

 //       if (printproc)   cout<<"Updating Preconditioner before nonlinear solve"<<endl;

	//A=Teuchos::rcp ( new Precon_Interface(comm,N,globalMap,precdata,precFunction));

     //   if (comm.MyPID()==0) printproc=true;

  //      if (printproc)   cout<<"N="<<N<<endl;

// if (comm.MyPID()==0) cout<<"comm.mypid="<<comm.MyPID()<<endl<<flush;


	A=Teuchos::rcp ( new Precon_Interface(N,globalMap,comm,precdata,precFunction));

   //     if (printproc)   cout<<"Updated Preconditioner before nonlinear solve"<<endl;
	return true;
}
*/



// Preconditioner is 2 steps. Only computePreconditioner is given the state,
//  which we store in solution.
bool Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                Epetra_Operator& Prec,
                Teuchos::ParameterList* p)
{
//        if (printproc)   cout<<"Updating Preconditioner before nonlinear solve"<<endl;


//This is how we define the non-seggregated Picard operator as a preconditioning operator
// precUpdateFunction(x.Values(),N,precdata);


//This is how we define the seggregated block Picard operator as a preconditioning operator
        //A=Teuchos::rcp ( new Precon_Interface(N,globalMap,comm,precdata,precFunctionblock11,precFunctionblock12,precFunctionblock21,precFunctionblock22));
        //A=Teuchos::rcp ( new Block_Precon_Interface(N,globalMap,comm,precdata,precFunctionblock11,precFunctionblock12,precFunctionblock21,precFunctionblock22));
        A=Teuchos::rcp ( new Block_Precon_Interface(N,globalMap,comm,precdata,auxprecFunctionblock11,auxprecFunctionblock12,auxprecFunctionblock21,auxprecFunctionblock22));


//This is how we setup the SIMPLE algorithm
//Convention: for SIMPLE algorithm based on Picard linearization
//11 block is F
//12 block is diag(F)^{-1}B'
//21 block is B
//22 block is S=G-Bdiag(F)^{-1}B'
        F=Teuchos::rcp ( new Precon_Interface(N,globalMap,comm,precdata,precFunctionblock11));
        S=Teuchos::rcp ( new Precon_Interface(N,globalMap,comm,precdata,precFunctionblock22));

        return true;
}









/* Analytic Jacobian */
#if 1
int Problem_Interface::Apply(const Epetra_MultiVector &X,Epetra_MultiVector &Y)const {

cout.precision(16);
       Y.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization
			 //jacFunction(X(0)->Values(),N,Y(0)->Values(),jacdata);
			 //Test with precdata here instead of jacdata, I believe they should be the same
//Print out vector we are applying J to before and after the apply
        //if (printproc)  cout<<"a="<<X<<endl;


     double n0; X(0)->Norm2(&n0);
if(printproc) cout << "Apply Norm of x="<<n0<<endl;

Epetra_Vector Jac(*Y(0)); //y=J_a*x
Epetra_Vector tempx(*X(0)); //x
Epetra_Vector a(*Y(0));
a.PutScalar(10.0);
getJacVec(a); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a


//double n1; a.Norm2(&n1);
//if(printproc) cout << "Analytic Jac Norm of a="<<n1<<endl;
     
//double n2; tempx.Norm2(&n2);
//if(printproc) cout << "Analytic Jac Norm of x="<<n2<<endl;

			 //jacFunction(X(0)->Values(),N,Y(0)->Values(),jacdata);
//			 jacFunction(tempx.Values(),N,Jac.Values(),blackbox_res);
			 jacFunction(tempx.Values(),N,Jac.Values(),jacdata);

//			 jacFunction(X(0)->Values(),N,Y(0)->Values(),blackbox_res);
//			 jacFunction(X(0)->Values(),N,Y(0)->Values(),blackbox_res);
        //if (printproc)  cout<<"b="<<Y<<endl;

//double n3; Jac.Norm2(&n3);
//if(printproc) cout << "Analytic Jac Norm of Jac="<<n3<<endl;


Y=Jac;
	return 0;
}
#endif

// Compare FD Jacobian with Analytic Jacobian
#if 0
int Problem_Interface::Apply(const Epetra_MultiVector &X,Epetra_MultiVector &Y)const {
//Finite Difference Jacobian calculation, can be used as alternative to trilinos FD calculation
// Should make a vector to store F(x) in order to just have 1 function eval per application
cout.precision(16);



       Y.PutScalar(0.0);

Epetra_Vector a(*Y(0));
a.PutScalar(10.0);
getJacVec(a); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a


//J_x (y)=(F(x+delta y)-F(x))/delta

double normJacx;a.Norm2(&normJacx);     // here x is the vector that we are linearizing about, this is passed in via jacdata J_x, we get this vector data using getJacVec(a);
double normJacy;X(0)->Norm2(&normJacy); //here y is the vector that we are applying J_x to (y is the const input vector X in Apply)
double lambda=1.e-6;
//double lambda=5.e-8;
double delta;
// Formulate norm of x and norm of y for a finite differnce jacobian calculation as done in trilinos


//     double na; a.Norm2(&na);
//if(printproc) cout << "FD jac norm a="<<na<<endl;

//     double nb; bb.Norm2(&nb);
//if(printproc) cout << "FD jac norm bb="<<nb<<endl;

if(abs(normJacy)<1.e-12){
delta=lambda*lambda;
return 0;
}
else
//delta=lambda*(lambda+normJacx/normJacy);
delta=1.e-8*(lambda+normJacx/normJacy);


Epetra_Vector tempx(a) ; //x
Epetra_Vector tempy(*X(0)); //y
Epetra_Vector tempa(*X(0));// y

     double n0; X(0)->Norm2(&n0);
//if(printproc) cout << "Apply Norm of x0="<<n0<<endl;


//     double ny; tempy.Norm2(&ny);
//if(printproc) cout << "FD jac norm of y="<<ny<<endl;


tempa.PutScalar(0.0);

Epetra_Vector tempb(*X(0));
tempb.PutScalar(0.0);
    for (int i=0; i<N; i++) tempx[i] += delta*tempy[i];

const Epetra_Vector b(tempx);


        //residualFunction(b.Values(), tempa.Values(), N,jacdata); //F(x+delta y)
        residualFunction(b.Values(), tempa.Values(), N,blackbox_res); //F(x+delta y)


 //    double n2; tempa.Norm2(&n2);
//if(printproc) cout << "Apply Norm of F="<<n2<<endl;


const Epetra_Vector c(a);
        //residualFunction(c.Values(), tempb.Values(), N, jacdata); //F(x)
        residualFunction(c.Values(), tempb.Values(), N, blackbox_res); //F(x)
     
//double n4; tempb.Norm2(&n4);
//if(printproc) cout << "Apply Norm of F="<<n4<<endl;

Epetra_Vector d(tempy);
d.PutScalar(0.0);

// i.e. J_x y= (F(x+delta y)-F(x))/delta, where delta= lambda(lambda +norm x/norm y)
    for (int i=0; i<N; i++) d[i] = (tempa[i]-tempb[i])/delta;

//double n5; d.Norm2(&n5);
//if(printproc) cout << "FD Jac Norm="<<n5<<endl;
Y=d;


Epetra_Vector Jac(*Y(0)); //y=J_a*x
Epetra_Vector xval(*X(0)); //x

			 jacFunction(xval.Values(),N,Jac.Values(),jacdata);
//double nc4; Jac.Norm2(&nc4);
//if(printproc) cout << "FD Analytic Jacobian Norm="<<nc4<<endl;


Epetra_Vector diff(*Y(0));
    for (int i=0; i<N; i++) diff[i] = (d[i]-Jac[i]);
double ndiff; diff.Norm2(&ndiff);
if(printproc) cout << "FD - Analytic Jacobian diff Norm="<<ndiff<<endl;


//Test Analytic Jacobian
 Y=Jac;

      //  return 0;


//Epetra_Vector Schur(*Y(0)); //y=J_a*x
//Epetra_Vector pval(*X(0)); //

//pval.PutScalar(1.0);
//precFunctionblock22(pval.Values(),N,Schur.Values(),jacdata);

//double nschur; Schur.Norm2(&nschur);
//if(printproc) cout << "Schur Norm="<<nschur<<endl;

return 0;

}

#endif



#if 1
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

//    if (printproc) {
//cout<<"rhs=["<<endl;
//    for (int i=0; i<N; i++) cout<<" "<<X[0][i];
//cout<<"];"<<endl;
//}

	Y=X;
	return 0;

}
#endif


#if 0
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

//	double n1; X(0)->Norm2(&n1); 
//if(printproc) cout << "ApplyInverse Norm of x="<<n1<<endl;

int numv= X.NumVectors();

//if(printproc) cout << "numv="<<numv<<endl;

	double n8; Y(0)->Norm2(&n8); 


       Y.PutScalar(0.0);

	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary + Belos::TimingDetails;
	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug;
	//int verb = Belos::TimingDetails;
	//int verb = 33;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	myPL.set( "Convergence Tolerance", 1.0e-10 );
	//myPL.set( "Convergence Tolerance", 1.0e-8 );
	//myPL.set( "Convergence Tolerance", 1.0e-11 );
	myPL.set( "Output Frequency", frequency );


	Teuchos::RCP<const MV>b=Teuchos::rcp ( new MV(X));
        Teuchos::RCP<MV>x=Teuchos::rcp ( new MV(Y));
	



 Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > myProblem = Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(A,x,b) );

	bool ret = myProblem->setProblem(); 
	if (ret != true) {
		cout<<"Error setting problem"<<endl;
	}
	Belos::BlockGmresSolMgr<ST,MV,OP> mySolver( myProblem, rcp(&myPL,false) );
        Belos::ReturnType solverRet = mySolver.solve();

	Teuchos::RCP<MV> sol = myProblem->getLHS();

	

	Y=*sol;
	double n0; Y(0)->Norm2(&n0); 
if(printproc) cout << "normprecvec="<<n0<<endl;


return 0;
} 
#endif 

/*
//Test FD
Epetra_Vector a(*Y(0));
a.PutScalar(10.0);
getJacVec(a); //copies the data from np1 entries of jacdata and put them into an Epetra_Vector a


//J_x (y)=(F(x+delta y)-F(x))/delta

double normJacx;a.Norm2(&normJacx);     // here x is the vector that we are linearizing about, this is passed in via jacdata J_x, we get this vector data using getJacVec(a);
double normJacy;Y(0)->Norm2(&normJacy); //here y is the vector that we are applying J_x to (y is the const input vector X in Apply)
double lambda=1.e-6;
//double lambda=5.e-8;
double delta;

if(abs(normJacy)<1.e-12){
delta=lambda*lambda;
return 0;
}
else
//delta=lambda*(lambda+normJacx/normJacy);
delta=1.e-6*(lambda+normJacx/normJacy);


Epetra_Vector tempx(a) ; //x
Epetra_Vector tempy(*Y(0)); //y
Epetra_Vector tempa(*Y(0));// y

     double n0; X(0)->Norm2(&n0);
tempa.PutScalar(0.0);

Epetra_Vector tempb(*X(0));
tempb.PutScalar(0.0);

    for (int i=0; i<N; i++) tempx[i] += delta*tempy[i];

const Epetra_Vector FDb(tempx);

//     double n1; b.Norm2(&n1);
//if(printproc) cout << "Apply Norm of x="<<n1<<endl;

        //residualFunction(b.Values(), tempa.Values(), N,jacdata); //F(x+delta y)
        residualFunction(FDb.Values(), tempa.Values(), N,blackbox_res); //F(x+delta y)


 //    double n2; tempa.Norm2(&n2);
//if(printproc) cout << "Apply Norm of F="<<n2<<endl;


const Epetra_Vector c(a);
//    double n3; c.Norm2(&n3);
//if(printproc) cout << "FD Jac Norm of x="<<n3<<endl;
        //residualFunction(c.Values(), tempb.Values(), N, jacdata); //F(x)
        residualFunction(c.Values(), tempb.Values(), N, blackbox_res); //F(x)

//double n4; tempb.Norm2(&n4);
//if(printproc) cout << "Apply Norm of F="<<n4<<endl;

Epetra_Vector FDJac(tempy);
FDJac.PutScalar(0.0);

// i.e. J_x y= (F(x+delta y)-F(x))/delta, where delta= lambda(lambda +norm x/norm y)
    for (int i=0; i<N; i++) FDJac[i] = (tempa[i]-tempb[i])/delta;

Epetra_Vector Xin(*X(0)); //X
Epetra_Vector FDdiff(*Y(0));
    for (int i=0; i<N; i++) FDdiff[i] = (Xin[i]-FDJac[i]);
double nFDdiff; FDdiff.Norm2(&nFDdiff);
if(printproc) cout << "X - FD*(Pinv*X)"<<nFDdiff<<endl;





//Test Jac

Epetra_Vector d(*X(0)); //X
Epetra_Vector Jac(*X(0)); //A*(Pinv*X)
Epetra_Vector xval(*Y(0)); //Pinv*X

jacFunction(xval.Values(),N,Jac.Values(),jacdata);


Epetra_Vector diff(*Y(0));
    for (int i=0; i<N; i++) diff[i] = (d[i]-Jac[i]);
double ndiff; diff.Norm2(&ndiff);
if(printproc) cout << "X - A*(Pinv*X)"<<ndiff<<endl;
return 0;
} 
*/




/* seggregated SIMPLE */

//ApplyInverse routine for the SIMPLE preconditioner
//Convention: for SIMPLE algorithm based on Picard linearization
//11 block is F
//12 block is diag(F)^{-1}B'
//21 block is B
//22 block is S=G-Bdiag(F)^{-1}B'

// These operators receive and return complete state vectors. 
// Components where the operator is not defined will return a zero 
// Specifically,
// F is applied to velocities and returns velocities along with a height field of 0
// similarly S is applied to a height field and returns a height field along with a velocity field of 0
// B is applied to a velocity field and returns a scalar height feild along with a zero velocity field
// diag(F)^{-1}B' is applied to a height field and returns a velocity field with a zero height field.
// communication within each of these operators uses the appropriate edge1,2,3 data structures


// rhs vector is X
// we implicitly view X as partitioned into components b1 and b2, 
// return vector is Y
// we implicitly view Y as partitioned into components Y1 and Y2, 


//SIMPLE algorithm
//
// 1. Solve F x1=b1
// 2. Solve S x2=-Bx1+b2
// 3a. y1= x1-(1/alpha)*diag(F)^{-1}B'x2
// 3b. y2= (1/alpha)*x2
#if 0
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

//X  RHS
//Y=Ainv*X Solution


int numv= X.NumVectors();


	double n8; X(0)->Norm2(&n8); 
if(printproc) cout << "Norm of RHS="<<n8<<endl;


       Y.PutScalar(0.0);


	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug;
	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	//myPL.set( "Convergence Tolerance", 1.0e-10);
	//myPL.set( "Convergence Tolerance", 1.0e-12 );
	myPL.set( "Convergence Tolerance", 1.0e-4 );
	myPL.set( "Output Frequency", frequency );


	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);


        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;


	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));

        Epetra_MultiVector x1(Y);
        x1(0)->PutScalar(0.0);

        Epetra_MultiVector y1(Y);
        y1(0)->PutScalar(0.0);



         double sum;
         B1.Norm1(&sum);





        if(sum<1.e-10){//if rhs is zero then don't solve
         if(printproc) cout<<"sum="<<sum<<endl;
         }
        else{
          Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          Fx->PutScalar(0.0);
          //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);

     double nfrhs; Fb->Norm2(&nfrhs);
if(printproc) cout << "normfrhs="<<nfrhs<<endl;


          Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );

	  bool Fret = FProblem->setProblem(); 
	  if (Fret != true) {
  		cout<<"Error setting problem"<<endl;
           }


if (printproc) {
    if (Fret == true) {
      cout << "Belos F Linear Problem Set" << std::endl;
    } else {
      cout << "Error setting Belos F Linear Problem" << std::endl;
    }
  }



	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

          Belos::ReturnType FsolverRet = FSolver.solve();

if (printproc) {
    if (FsolverRet == Belos::Converged) {
      cout << "Belos F converged." << std::endl;
    } else {
      cout << "Belos F did not converge." << std::endl;
    }
  }



	  Teuchos::RCP<MV> FSol= FProblem->getLHS();
          x1=*FSol;
          y1=*Fb;

Epetra_Vector tempx1a(*x1(0));
     double npva; tempx1a.Norm2(&npva);
if(printproc) cout << "fsolnorm="<<npva<<endl;


Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
if(printproc) cout << "frhsnorm="<<npya<<endl;

        }

// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

        Epetra_Vector bx1(*Y(0));
        bx1.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization

        precFunctionblock21(x1(0)->Values(),N, bx1.Values(), precdata);

     double nB21; bx1.Norm2(&nB21);
if(printproc) cout << "normB21="<<nB21<<endl;

// Then set Sproblos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,x,b) );


// Insert loop here to add -bx1 and b2 and put in Schurb
// b2 the second component of the X vector

        Epetra_MultiVector b2(X);
        for (int i=0; i<2*(N+1)/3;i++) b2[0][i] = 0.0;


	double nsa; b2(0)->Norm2(&nsa); 
if(printproc) cout << "Norm of RHS Schur A="<<nsa<<endl;


//    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
if(printproc) cout << "Norm of RHS Schur B="<<nsb<<endl;

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));



        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );






        bool Sret = SchurProblem->setProblem(); 


if (printproc) {
    if (Sret == true) {
      cout << "Belos S Linear Problem Set" << std::endl;
    } else {
      cout << "Error setting Belos S Linear Problem" << std::endl;
    }
  }


        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
        Belos::ReturnType SchursolverRet = SchurSolver.solve();

if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      cout << "Belos Schur converged." << std::endl;
    } else {
      cout << "Belos Schur did not converge." << std::endl;
    }
  }



	Teuchos::RCP<MV> SchurSol= SchurProblem->getLHS();
        Epetra_MultiVector x2(*SchurSol);
        //Epetra_MultiVector x2(b2);



//Next apply dDinvBt to x1 and store in Bx1

        Epetra_Vector dFinvBt(*Y(0));
        dFinvBt.PutScalar(0.0);

	precFunctionblock12(x2(0)->Values(),N, dFinvBt.Values(), precdata);

     double nBt; dFinvBt.Norm2(&nBt);
if(printproc) cout << "normBt="<<nBt<<endl;



// Insert loop here to add x1 -(1/alpha)*dFinvBt*x2
// Since the vector is broken up into two components with zeroes padding the 
// 'non-active' entries we can build Y in one action
// We write down that formally this is being done:
// y1= x1 -(1/alpha)*dFinvBt*x2
// y2= (1/alpha)x2 
// 
// which is equivalentlly implemented as Y= x1 -(1/alpha)*dFinvBt*x2 + (1/alpha)x2

Epetra_Vector tempx1(*x1(0));
Epetra_Vector tempx2(*x2(0));
//double alpha=1.0;
double alpha=1.0;
double alphainv=1.0/alpha;





// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU




Y=tempx1;

return 0;
}
#endif

























/* Testing seggregated simple vs unseggregated simple */

//ApplyInverse routine for the SIMPLE preconditioner
//Convention: for SIMPLE algorithm based on Picard linearization
//11 block is F
//12 block is diag(F)^{-1}B'
//21 block is B
//22 block is S=G-Bdiag(F)^{-1}B'

// These operators receive and return complete state vectors. 
// Components where the operator is not defined will return a zero 
// Specifically,
// F is applied to velocities and returns velocities along with a height field of 0
// similarly S is applied to a height field and returns a height field along with a velocity field of 0
// B is applied to a velocity field and returns a scalar height feild along with a zero velocity field
// diag(F)^{-1}B' is applied to a height field and returns a velocity field with a zero height field.
// communication within each of these operators uses the appropriate edge1,2,3 data structures


// rhs vector is X
// we implicitly view X as partitioned into components b1 and b2, 
// return vector is Y
// we implicitly view Y as partitioned into components Y1 and Y2, 


//SIMPLE algorithm
//
// 1. Solve F x1=b1
// 2. Solve S x2=-Bx1+b2
// 3a. y1= x1-(1/alpha)*diag(F)^{-1}B'x2
// 3b. y2= (1/alpha)*x2
#if 0
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

//X  RHS
//Y=Ainv*X Solution

//if(printproc) cout << "ApplyInverse Norm of x="<<n1<<endl;

int numv= X.NumVectors();

//if(printproc) cout << "numv="<<numv<<endl;

	double n8; X(0)->Norm2(&n8); 
if(printproc) cout << "Norm of RHS="<<n8<<endl;


       Y.PutScalar(0.0);
//if (printproc)   cout<<"X1"<<X<<endl;

// precFunction(X(0)->Values(),N,Y(0)->Values(),precdata);
//if (printproc)  cout << "X address in:  " << &X << "   Y address in:  " << &Y << endl;
    
 //   if (printproc)   cout<<"In ApplyInverse"<<endl;
  //  if (printproc)   cout<<"N="<<N<<endl;
   // if (comm.MyPID()==0) cout<<"comm.mypid="<<comm.MyPID()<<endl<<flush;

	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug;
	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	//myPL.set( "Convergence Tolerance", 1.0e-10);
	//myPL.set( "Convergence Tolerance", 1.0e-12 );
	myPL.set( "Convergence Tolerance", 1.0e-4 );
	myPL.set( "Output Frequency", frequency );


	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);


        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;


	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));

        Epetra_MultiVector x1(Y);
        x1(0)->PutScalar(0.0);

        Epetra_MultiVector y1(Y);
        y1(0)->PutScalar(0.0);






#if 0
       Epetra_MultiVector b1(X);

    //for (int i=0; i<N; i++) b2[0][i] = b2[0][i]-bx1[i];

        Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(b1));

        Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));

        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );

        bool Fret = FProblem->setProblem();
        if (Sret != true) {
                cout<<"Error setting problem"<<endl;
        }
        Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );
        Belos::ReturnType FsolverRet = FSolver.solve();


 if (FsolverRet!=Belos::Converged ) {
       if(printproc) cout << "Belos Fsolve Failed"<<endl;
  }


        Teuchos::RCP<MV> FSol= FProblem->getLHS();
        Epetra_MultiVector x1(*FSol);
        //Epetra_MultiVector x2(b2);

#endif




         double sum;
         B1.Norm1(&sum);





     //    double sum;
      //   for (int i=0; i<2*(N+1)/3; i++) {
       //      sum+=X[0][i]*X[0][i];
        //  }
//Test M and Minv
/*
        if(1){//if rhs is zero then don't solve


     double nx1;B1(0)->Norm2(&nx1);
if(printproc) cout << "normb1="<<nx1<<endl;
        //auxprecFunctionblock11(B1(0)->Values(),N, x1(0)->Values(), precdata);
        //precFunctionblock11(x1(0)->Values(),N, y1(0)->Values(), precdata);

        precFunctionblock11(B1(0)->Values(),N, x1(0)->Values(), precdata);
        auxprecFunctionblock11(x1(0)->Values(),N, y1(0)->Values(), precdata);


     double ny1; y1(0)->Norm2(&ny1);
if(printproc) cout << "normy1="<<ny1<<endl;


Epetra_Vector diagdiff(*Y(0));
    for (int i=0; i<N; i++) diagdiff[i] = (B1[0][i]-y1[0][i]);

double diagndiff; diagdiff.Norm2(&diagndiff);
if(printproc) cout << "diag vec diff="<<diagndiff/nx1<<endl;
}
*/
//end M and Minv


#if 0
//Test diag(F)*inv(Diag(F))
Epetra_Vector ei1(*Y(0));
        ei1.PutScalar(1.0);
for (int i=2*(N+1)/3; i<N; i++){
ei1[i]=0.0;
}

Epetra_Vector fd1(*Y(0));
        fd1.PutScalar(0.0);

Epetra_Vector f1(*Y(0));
        f1.PutScalar(0.0);

//cout << "ei1="<<ei1<<endl<<flush;
        precFunctionblock11(ei1.Values(),N, fd1.Values(), precdata);
//cout << "fd="<<fd1<<endl<<flush;
        auxprecFunctionblock11(fd1.Values(),N, f1.Values(), precdata);
//cout << "f1="<<f1<<endl<<flush;


Epetra_Vector diagdiff1(*Y(0));
        diagdiff1.PutScalar(0.0);
for (int i=0; i<N; i++){
diagdiff1[i]=ei1[i]-f1[i];
}

double diagndiff1; diagdiff1.Norm2(&diagndiff1);
if(printproc) cout << "diag vec diff norm="<<diagndiff1<<endl<<flush;
// end test diag(F)*inv(Diag(F))
#endif






#if 0
//Test inv(Diag(F))*B' -(inv(Diag(F))B')
Epetra_Vector ei1(*Y(0));
        ei1.PutScalar(0.0);



for (int i=2*(N+1)/3; i<N; i++){
//ei1[i]=0.0;
ei1[i]=rand()*acos(-1.0);
}

Epetra_Vector fd1(*Y(0));
        fd1.PutScalar(0.0);

Epetra_Vector f1(*Y(0));
        f1.PutScalar(0.0);

Epetra_Vector fb1(*Y(0));
        fb1.PutScalar(0.0);

//cout << "ei1="<<ei1<<endl<<flush;
        precFunctionblock12(ei1.Values(),N, fd1.Values(), precdata);
        precFunctionblock11(fd1.Values(),N, fb1.Values(), precdata);
        auxprecFunctionblock12(ei1.Values(),N, f1.Values(), precdata);


Epetra_Vector diagdiff1(*Y(0));
        diagdiff1.PutScalar(0.0);
for (int i=0; i<N; i++){
diagdiff1[i]=fb1[i]-f1[i];
if(printproc) cout<<fd1[i]<<f1[i]<<endl;
}

double diagndiff1; diagdiff1.Norm2(&diagndiff1);
if(printproc) cout << "FDi Bt-FDiBt vec diff norm="<<diagndiff1<<endl<<flush;
// end test diag(F)*inv(Diag(F))
#endif





#if 0
//Test F*(inv(Diag(F))*B') -(F*inv(diag(f))*B')
Epetra_Vector ei1(*Y(0));
        ei1.PutScalar(0.0);



for (int i=2*(N+1)/3; i<N; i++){
//ei1[i]=0.0;
ei1[i]=rand()*acos(-1.0);
}

Epetra_Vector fd1(*Y(0));
        fd1.PutScalar(0.0);
Epetra_Vector fd2(*Y(0));
        fd2.PutScalar(0.0);

Epetra_Vector f1(*Y(0));
        f1.PutScalar(0.0);

Epetra_Vector fb1(*Y(0));
        fb1.PutScalar(0.0);

//cout << "ei1="<<ei1<<endl<<flush;
        precFunctionblock12(ei1.Values(),N, fd1.Values(), precdata);
        precFunctionblock11(fd1.Values(),N, fd2.Values(), precdata);
        auxprecFunctionblock12(ei1.Values(),N, f1.Values(), precdata);
        precFunctionblock22(ei1.Values(),N, fb1.Values(), precdata);


Epetra_Vector diagdiff1(*Y(0));
        diagdiff1.PutScalar(0.0);

Epetra_Vector diagdiff2(*Y(0));
        diagdiff2.PutScalar(0.0);

Epetra_Vector diagdiff3(*Y(0));
        diagdiff3.PutScalar(0.0);
for (int i=0; i<N; i++){
diagdiff1[i]=fd2[i]-f1[i];
diagdiff2[i]=fd2[i]-fb1[i];
diagdiff3[i]=fb1[i]-f1[i];
//if(printproc) cout<<fd1[i]<<f1[i]<<endl;
}


double nfd2; fd2.Norm2(&nfd2);
if(printproc) cout << "nfd2="<<nfd2<<endl<<flush;

double nf1; f1.Norm2(&nf1);
if(printproc) cout << "nf1="<<nf1<<endl<<flush;


double diagndiff1; diagdiff1.Norm2(&diagndiff1);
if(printproc) cout << "F( FDi Bt)-(F FDi FBt) vec diff norm="<<diagndiff1<<endl<<flush;

double diagndiff2; diagdiff2.Norm2(&diagndiff2);
if(printproc) cout << "FD*(FDi Bt)-Bt vec diff norm="<<diagndiff2<<endl<<flush;

double diagndiff3; diagdiff3.Norm2(&diagndiff3);
if(printproc) cout << "Bt- (FD* FDi FBt) vec diff norm="<<diagndiff3<<endl<<flush;

// end test diag(F)*inv(Diag(F))
#endif







//Test diag(F) and Fdiag
        if(0){//if rhs is zero then don't solve

     //double nfd;B1(0)->Norm2(&nx1);
//if(printproc) cout << "normb1="<<nx1<<endl;



Epetra_Vector diagdiff(*Y(0));
        diagdiff.PutScalar(0.0);
Epetra_Vector ei(*Y(0));
        ei.PutScalar(0.0);
Epetra_Vector fd(*Y(0));
        fd.PutScalar(0.0);
Epetra_Vector f(*Y(0));
        f.PutScalar(0.0);

Epetra_Vector fdd(*Y(0));
        fdd.PutScalar(0.0);
Epetra_Vector DDf(*Y(0));
        DDf.PutScalar(0.0);



    for (int i=0; i<N; i++){
        ei.PutScalar(0.0);
        fd.PutScalar(0.0);
        f.PutScalar(0.0);
        ei[i]=1.0;

        precFunctionblock11(ei.Values(),N, fd.Values(), precdata);
        auxprecFunctionblock11(ei.Values(),N, f.Values(), precdata);
DDf[i]=f[i];
fdd[i]=fd[i];
        diagdiff[i] = (f[i]-fd[i]);
//if(printproc) cout << "DDf("<<i+1<<")=["<<DDf[i]<<"];"<<endl;
//if(printproc) cout << "fdd("<<i+1<<")=["<<fdd[i]<<"];"<<endl;
      }

double diagndiff; diagdiff.Norm2(&diagndiff);
if(printproc) cout << "diag vec diff="<<diagndiff<<endl;


double DDfnorm; DDf.Norm2(&DDfnorm);
if(printproc) cout << "Diag(F)norm="<<DDfnorm<<endl;

double fddnorm; fdd.Norm2(&fddnorm);
if(printproc) cout << "FDiagnorm="<<fddnorm<<endl;

if(printproc) cout << "N="<<N<<endl;





}
//end test diag(F) FDiag 



//Test G -B*DFInv*Bt 
/*
        if(1){//if rhs is zero then don't solve

        Epetra_MultiVector xb2(X);

        for (int i=0; i<2*(N+1)/3;i++) xb2[0][i] = 0.0;

     Epetra_Vector xdFinvBt(*Y(0));
        xdFinvBt.PutScalar(0.0);
     
    Epetra_Vector ydFinvBt(*Y(0));
        ydFinvBt.PutScalar(0.0);

    Epetra_Vector zdFinvBt(*Y(0));
        zdFinvBt.PutScalar(0.0);


     double nxb2; xb2.Norm2(&nxb2);
if(printproc) cout << "normxb2="<<nxb2<<endl;

       // DFinv*Bt
        precFunctionblock12(xb2(0)->Values(),N, xdFinvBt.Values(), precdata);
        //precFunctionblock12(xb2(0)->Values(),N, ydFinvBt.Values(), precdata);
        //B*DFinv*Bt
        auxprecFunctionblock21(xdFinvBt.Values(),N, ydFinvBt.Values(), precdata);


        //G
        auxprecFunctionblock22(xb2(0)->Values(),N, zdFinvBt.Values(), precdata);


     double cb1; ydFinvBt.Norm2(&cb1);
if(printproc) cout << "S1="<<cb1<<endl;

        //G-B*DFinv*Bt
 xdFinvBt.PutScalar(0.0);
       precFunctionblock22(xb2(0)->Values(),N, xdFinvBt.Values(), precdata);

     double cb2; xdFinvBt.Norm2(&cb2);
if(printproc) cout << "S2="<<cb2<<endl;

        Epetra_Vector xbdiff(*X(0));
        xbdiff.PutScalar(0.0);
        for (int i=0; i<N;i++) xbdiff[i] = xdFinvBt[i]-(zdFinvBt[i]-ydFinvBt[i]); //+ to compare -BDFinvBt - -BDFinvBt

     double nyb2; xbdiff.Norm2(&nyb2);
if(printproc) cout << "normyb2="<<nyb2<<endl;
         }
*/
//End Test  DFInvBt 


        if(sum<1.e-10){//if rhs is zero then don't solve
         if(printproc) cout<<"sum="<<sum<<endl;
         }
        else{
          Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          Fx->PutScalar(0.0);
          //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);

     double nfrhs; Fb->Norm2(&nfrhs);
if(printproc) cout << "normfrhs="<<nfrhs<<endl;


          Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );

	  bool Fret = FProblem->setProblem(); 
//	  if (Fret != true) {
 // 		cout<<"Error setting problem"<<endl;
  //         }


if (printproc) {
    if (Fret == true) {
      cout << "Belos F Linear Problem Set" << std::endl;
    } else {
      cout << "Error setting Belos F Linear Problem" << std::endl;
    }
  }



	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

          Belos::ReturnType FsolverRet = FSolver.solve();

if (printproc) {
    if (FsolverRet == Belos::Converged) {
      cout << "Belos F converged." << std::endl;
    } else {
      cout << "Belos F did not converge." << std::endl;
    }
  }



	  Teuchos::RCP<MV> FSol= FProblem->getLHS();
          x1=*FSol;
          y1=*Fb;

Epetra_Vector tempx1a(*x1(0));
     double npva; tempx1a.Norm2(&npva);
if(printproc) cout << "fsolnorm="<<npva<<endl;


Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
if(printproc) cout << "frhsnorm="<<npya<<endl;

        }

// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

        Epetra_Vector bx1(*Y(0));
        bx1.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization

        precFunctionblock21(x1(0)->Values(),N, bx1.Values(), precdata);

     double nB21; bx1.Norm2(&nB21);
if(printproc) cout << "normB21="<<nB21<<endl;

// Then set Sproblos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,x,b) );


// Insert loop here to add -bx1 and b2 and put in Schurb
// b2 the second component of the X vector

        Epetra_MultiVector b2(X);
        for (int i=0; i<2*(N+1)/3;i++) b2[0][i] = 0.0;


	double nsa; b2(0)->Norm2(&nsa); 
if(printproc) cout << "Norm of RHS Schur A="<<nsa<<endl;


//    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
if(printproc) cout << "Norm of RHS Schur B="<<nsb<<endl;

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));



        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );






        bool Sret = SchurProblem->setProblem(); 


if (printproc) {
    if (Sret == true) {
      cout << "Belos S Linear Problem Set" << std::endl;
    } else {
      cout << "Error setting Belos S Linear Problem" << std::endl;
    }
  }


        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
        Belos::ReturnType SchursolverRet = SchurSolver.solve();

if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      cout << "Belos Schur converged." << std::endl;
    } else {
      cout << "Belos Schur did not converge." << std::endl;
    }
  }



	Teuchos::RCP<MV> SchurSol= SchurProblem->getLHS();
        Epetra_MultiVector x2(*SchurSol);
        //Epetra_MultiVector x2(b2);



//Next apply dDinvBt to x1 and store in Bx1

        Epetra_Vector dFinvBt(*Y(0));
        dFinvBt.PutScalar(0.0);

	precFunctionblock12(x2(0)->Values(),N, dFinvBt.Values(), precdata);

     double nBt; dFinvBt.Norm2(&nBt);
if(printproc) cout << "normBt="<<nBt<<endl;



// Insert loop here to add x1 -(1/alpha)*dFinvBt*x2
// Since the vector is broken up into two components with zeroes padding the 
// 'non-active' entries we can build Y in one action
// We write down that formally this is being done:
// y1= x1 -(1/alpha)*dFinvBt*x2
// y2= (1/alpha)x2 
// 
// which is equivalentlly implemented as Y= x1 -(1/alpha)*dFinvBt*x2 + (1/alpha)x2

Epetra_Vector tempx1(*x1(0));
Epetra_Vector tempx2(*x2(0));
//double alpha=1.0;
double alpha=.99;
double alphainv=1.0/alpha;


/*
         for (int i=0; i<2*(N+1)/3; i++) {
             tempx1[i]=x1[0][i]-alphainv*dFinvBt[i];
             //tempx1[i]=x1[0][i]-0.0*dFinvBt[i];
          }
         for (int i=2*(N+1)/3;i<N; i++) {
             tempx1[i]=(alphainv*tempx2[i]);
          }
*/
        //Epetra_Vector Ix1(*X(0));
         //for (int i=2*(N+1)/3;i<N; i++) {Ix1[i]=0.0;}
    //for (int i=0; i<N; i++) tempx1[i] = (Ix1[i]+tempx2[i]);
//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]+tempx2[i]);


     double npvf; tempx1.Norm2(&npvf);
if(printproc) cout << "Diagonal components"<<endl;
if(printproc) cout << "factoredprecnormvec1="<<npvf<<endl;
double npvs; tempx2.Norm2(&npvs);
if(printproc) cout << "factoredprecnormvec2="<<npvs<<endl;


// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU





     double npv; tempx1.Norm2(&npv);
if(printproc) cout << "factoredprecnorm="<<npv<<endl;

Y=tempx1;


       Y.PutScalar(0.0);

        typedef double ST;
        typedef Epetra_MultiVector MV;
        typedef Epetra_Operator OP;
        typedef Belos::MultiVecTraits<ST,MV>    MVT;
        typedef Belos::OperatorTraits<ST,MV,OP> OPT;


        Teuchos::RCP<const MV>b=Teuchos::rcp ( new MV(X));
        Teuchos::RCP<MV>x=Teuchos::rcp ( new MV(Y));


        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > myProblem = Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(A,x,b) );

        bool ret = myProblem->setProblem();


if (printproc) {
    if (ret == true) {
      cout << "Belos Full Linear Problem Set" << std::endl<<flush;
    } else {
      cout << "Error setting Belos Full Linear Problem" << std::endl<<flush;
    }
  }


        Belos::BlockGmresSolMgr<ST,MV,OP> mySolver( myProblem, rcp(&myPL,false) );
        Belos::ReturnType solverRet = mySolver.solve();

if (printproc) {
    if (solverRet == Belos::Converged) {
      cout << "Belos Full converged." << std::endl<<flush;
    } else {
      cout << "Belos Full did not converge." << std::endl<<flush;
    }
  }


        Teuchos::RCP<MV> sol = myProblem->getLHS();

        Y=*sol;
        double n0; Y(0)->Norm2(&n0);
        if(printproc) cout << "fullprecnorm="<<n0<<endl<<flush;





Epetra_Vector PC(*Y(0)); //y=PC*x


Epetra_Vector PC11(*Y(0)); //y=PC*x
Epetra_Vector PC22(*Y(0)); //y=PC*x

Epetra_Vector FPC1(tempx1); //y=PC*x
Epetra_Vector FPC2(tempx1); //y=PC*x


    for (int i=0; i<2*(N+1)/3; i++){
 FPC2[i] = 0.0;
 PC22[i] = 0.0;
}

    for (int i=2*(N+1)/3;i<N; i++){
 FPC1[i] = 0.0;
 PC11[i] = 0.0;
}

double npvf11; PC11.Norm2(&npvf11);
if(printproc) cout << "normprecvecfull1="<<npvf11<<endl;
double npvs22; PC22.Norm2(&npvs22);
if(printproc) cout << "normprecvecfull2="<<npvs22<<endl;



//Factored Solution Norms
double fnpvf11; FPC1.Norm2(&fnpvf11);
if(printproc) cout << "fnormprecvecfull1="<<fnpvf11<<endl;
double fnpvs22; FPC2.Norm2(&fnpvs22);
if(printproc) cout << "fnormprecvecfull2="<<fnpvs22<<endl;





Epetra_Vector diff(*Y(0));
    for (int i=0; i<N; i++) diff[i] = (PC[i]-tempx1[i]);

double ndiff; diff.Norm2(&ndiff);
if(printproc) cout << "precon vec diff="<<ndiff<<endl;

if(printproc) cout << "relative precon vec diff="<<ndiff/n0<<endl;


Y=tempx1;

return 0;
}
#endif






#if 0
/* SIMPLE  */

//ApplyInverse routine for the SIMPLE preconditioner
//Convention: for SIMPLE algorithm based on Picard linearization
//11 block is F
//12 block is diag(F)^{-1}B'
//21 block is B
//22 block is S=G-Bdiag(F)^{-1}B'

// These operators receive and return complete state vectors. 
// Components where the operator is not defined will return a zero 
// Specifically,
// F is applied to velocities and returns velocities along with a height field of 0
// similarly S is applied to a height field and returns a height field along with a velocity field of 0
// B is applied to a velocity field and returns a scalar height feild along with a zero velocity field
// diag(F)^{-1}B' is applied to a height field and returns a velocity field with a zero height field.
// communication within each of these operators uses the appropriate edge1,2,3 data structures


// rhs vector is X
// we implicitly view X as partitioned into components b1 and b2, 
// return vector is Y
// we implicitly view Y as partitioned into components Y1 and Y2, 


//SIMPLE algorithm
//
// 1. Solve F x1=b1
// 2. Solve S x2=-Bx1+b2
// 3a. y1= x1-(1/alpha)*diag(F)^{-1}B'x2
// 3b. y2= (1/alpha)*x2
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

          if (printproc) cout << "In ApplyInverse" << flush<<endl;
          if (printproc) system("date");
//X  RHS
//Y=Ainv*X Solution

int numv= X.NumVectors();

	double n8; X(0)->Norm2(&n8); 
if(printproc) cout << "Norm of RHS="<<n8<<endl;


       Y.PutScalar(0.0);

	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug+Belos::TimingDetails+Belos::IterationDetails;
	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	//myPL.set( "Convergence Tolerance", 1.0e-10);
	myPL.set( "Convergence Tolerance", 1.0e-4 );
	myPL.set( "Output Frequency", frequency );
	myPL.set( "Timer Label","F Solve"  );

if(printproc) cout << "Set Preonditioner Parameterlist"<<flush<<endl;

	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);


        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;

if(printproc) cout << "Initialized B1"<<flush<<endl;

	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));
          double nfrhs; Fb->Norm2(&nfrhs);
          if(printproc) cout << "normfrhs="<<nfrhs<<flush<<endl;

        Epetra_MultiVector x1(Y);
        x1(0)->PutScalar(0.0);

        Epetra_MultiVector y1(Y);
        y1(0)->PutScalar(0.0);

          double sum;
         B1.Norm1(&sum);


/*         double sum;
         for (int i=0; i<2*(N+1)/3; i++) {
             sum+=X[0][i]*X[0][i];
          }
*/

if(printproc) cout << "Set sum"<<flush<<endl;

        if(sum<1.e-8){//if rhs is zero then don't solve
          if(printproc)cout<<"rhs sum="<<sum<<" returning 0 solution "<<flush<<endl;
         }
        else{
          if(printproc)cout<<"rhs sum="<<sum<<"solving with GMRES"<<flush<<endl;
          Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          if(printproc) cout << "created rcp fx"<<flush<<endl;
          Fx->PutScalar(0.0);
          if(printproc) cout << "assigned fx"<<flush<<endl;
          //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);


          Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );

	  bool Fret = FProblem->setProblem(); 


          if (printproc) {
            if (Fret == true) {
             cout << "Belos F Linear Problem Set" << flush<<std::endl;
             } 
            else {
             cout << "Error setting Belos F Linear Problem" << std::endl;
            }
          }



	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

          if(printproc) cout << "GMRES Solver Manager set"<<flush<<endl;
          Belos::ReturnType FsolverRet = FSolver.solve();

if (printproc) {
    if (FsolverRet == Belos::Converged) {
      cout << "Belos F converged."<<flush << std::endl;
    } else {
      cout << "Belos F did not converge." <<flush<< std::endl;
    }
  }



	  Teuchos::RCP<MV> FSol= FProblem->getLHS();
          x1=*FSol;
          y1=*Fb;

Epetra_Vector tempx1a(*x1(0));
     double npva; tempx1a.Norm2(&npva);
if(printproc) cout << "fsolnorm="<<npva<<endl;


Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
if(printproc) cout << "frhsnorm="<<npya<<endl;

        }



	myPL.set( "Timer Label","Schur Solve"  );


// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

        Epetra_Vector bx1(*Y(0));
        bx1.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization

        precFunctionblock21(x1(0)->Values(),N, bx1.Values(), precdata);

     double nB21; bx1.Norm2(&nB21);
if(printproc) cout << "normB21="<<nB21<<endl;

// Then set Sproblos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,x,b) );


// Insert loop here to add -bx1 and b2 and put in Schurb
// b2 the second component of the X vector

        Epetra_MultiVector b2(X);
        for (int i=0; i<2*(N+1)/3;i++) b2[0][i] = 0.0;


	double nsa; b2(0)->Norm2(&nsa); 
if(printproc) cout << "Norm of RHS Schur A="<<nsa<<endl;


//    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
if(printproc) cout << "Norm of RHS Schur B="<<nsb<<flush<<endl;

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));

if(printproc) cout << "Schur rhs set="<<flush<<endl;


        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

if(printproc) cout << "Schur solution initialized "<<flush<<endl;
        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );


if(printproc) cout << "Schur Problem initialized"<<flush<<endl;


        bool Sret = SchurProblem->setProblem(); 


if (printproc) {
    if (Sret == true) {
      cout << "Belos S Linear Problem Set"<<flush<< std::endl;
    } else {
      cout << "Error setting Belos S Linear Problem" <<flush<< std::endl;
    }
  }


        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
if(printproc) cout << "Schur GMRES Solver Set"<<flush<<endl;
        Belos::ReturnType SchursolverRet = SchurSolver.solve();



if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      cout << "Belos Schur converged." <<flush<< std::endl;
    } else {
      cout << "Belos Schur did not converge." << flush<<std::endl;
    }
  }



	Teuchos::RCP<MV> SchurSol= SchurProblem->getLHS();
        Epetra_MultiVector x2(*SchurSol);



//Next apply dDinvBt to x1 and store in Bx1

        Epetra_Vector dFinvBt(*Y(0));
        dFinvBt.PutScalar(0.0);

	precFunctionblock12(x2(0)->Values(),N, dFinvBt.Values(), precdata);

     double nBt; dFinvBt.Norm2(&nBt);
if(printproc) cout << "normBt="<<nBt<<flush<<endl;



// Insert loop here to add x1 -(1/alpha)*dFinvBt*x2
// Since the vector is broken up into two components with zeroes padding the 
// 'non-active' entries we can build Y in one action
// We write down that formally this is being done:
// y1= x1 -(1/alpha)*dFinvBt*x2
// y2= (1/alpha)x2 
// 
// which is equivalentlly implemented as Y= x1 -(1/alpha)*dFinvBt*x2 + (1/alpha)x2

Epetra_Vector tempx1(*x1(0));
Epetra_Vector tempx2(*x2(0));
double alpha=1.0;
double alphainv=1.0/alpha;


/*
         for (int i=0; i<2*(N+1)/3; i++) {
             tempx1[i]=x1[0][i]-alphainv*dFinvBt[i];
             //tempx1[i]=x1[0][i]-0.0*dFinvBt[i];
          }
         for (int i=2*(N+1)/3;i<N; i++) {
             tempx1[i]=(alphainv*tempx2[i]);
          }
*/
        //Epetra_Vector Ix1(*X(0));
         //for (int i=2*(N+1)/3;i<N; i++) {Ix1[i]=0.0;}
    //for (int i=0; i<N; i++) tempx1[i] = (Ix1[i]+tempx2[i]);
//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]+tempx2[i]);


     double npvf; tempx1.Norm2(&npvf);
if(printproc) cout << "Diagonal components"<<flush<<endl;
if(printproc) cout << "factoredprecnormvec1="<<npvf<<flush<<endl;
double npvs; tempx2.Norm2(&npvs);
if(printproc) cout << "factoredprecnormvec2="<<npvs<<flush<<endl;


// for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU
//   for (int i=0; i<N; i++) tempx1[i] = tempx1[i]+(alphainv*tempx2[i]);//Diagonal and also DU


//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-dFinvBt[i]) +(tempx2[i]); //No G no alpha


     double npv; tempx1.Norm2(&npv);
if(printproc) cout << "factoredprecnorm="<<npv<<endl;

Y=tempx1;

          if (printproc) system("date");
          if (printproc) cout << "Leaving ApplyInverse" << flush<<endl;

return 0;
}
#endif




#if 0
int Problem_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
{

          if (printproc) cout << "In ApplyInverse" << flush<<endl;
//X  RHS
//Y=Ainv*X Solution

int numv= X.NumVectors();

	double n8; X(0)->Norm2(&n8); 
if(printproc) cout << "Norm of RHS="<<n8<<endl;


       Y.PutScalar(0.0);

	using Teuchos::ParameterList;
	using Teuchos::RCP; 
	using Teuchos::rcp;

	typedef double ST;
	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Belos::MultiVecTraits<ST,MV>	MVT; 
	typedef Belos::OperatorTraits<ST,MV,OP> OPT;

	//int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary+Belos::Debug;
	int verb = Belos::Warnings + Belos::Errors + Belos::FinalSummary;
	int frequency=1;

	ParameterList myPL;
	myPL.set( "Verbosity", verb ); 
	myPL.set( "Block Size", 1 ); 
	myPL.set( "Num Blocks", 100 ); 
	myPL.set( "Maximum Iterations", 500 ); 
	//myPL.set( "Convergence Tolerance", 1.0e-10);
	myPL.set( "Convergence Tolerance", 1.0e-4 );
	myPL.set( "Output Frequency", frequency );

if(printproc) cout << "Set Preonditioner Parameterlist"<<flush<<endl;

	Epetra_MultiVector B1(X);
	Epetra_MultiVector B2(X);


        for (int i=2*(N+1)/3;i<N; i++) B1[0][i] = 0.0;

if(printproc) cout << "Initialized B1"<<flush<<endl;

	Teuchos::RCP<const MV>Fb=Teuchos::rcp ( new MV(B1));
          double nfrhs; Fb->Norm2(&nfrhs);
          if(printproc) cout << "normfrhs="<<nfrhs<<flush<<endl;

        Epetra_MultiVector x1(Y);
        x1(0)->PutScalar(0.0);

        Epetra_MultiVector y1(Y);
        y1(0)->PutScalar(0.0);


         double sum;
         B1.Norm1(&sum);
//         for (int i=0; i<2*(N+1)/3; i++) {
 //            sum+=X[0][i]*X[0][i];
  //        }

// sum this across processors

if(printproc) cout << "Set sum"<<flush<<endl;
cout<<"comm.mypid="<<comm.MyPID()<<"   normfrhs"<<sum<<endl<<flush;

        if(sum<1.e-10){//if rhs is zero then don't solve
          if(printproc)cout<<"rhs sum="<<sum<<" returning 0 solution "<<flush<<endl;
         }
        else{
          if(printproc)cout<<"rhs sum="<<sum<<" solving with GMRES"<<flush<<endl;
          Teuchos::RCP<MV>Fx=Teuchos::rcp ( new MV(Y));
          if(printproc) cout << "created rcp fx"<<flush<<endl;
          Fx->PutScalar(0.0);
          if(printproc) cout << "assigned fx"<<flush<<endl;
          //We initialized Y to zero and now have based Fx on Y //Fx.PutScalar(0.0);


          Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,Fx,Fb) );
          if(printproc) cout << "Problem declared"<<flush<<endl;

	  bool Fret = FProblem->setProblem(); 
          if(printproc) cout << "Problem set"<<flush<<endl;


          if (printproc) {
            if (Fret == true) {
             cout << "Belos F Linear Problem Set" << flush<<std::endl;
             } 
            else {
             cout << "Error setting Belos F Linear Problem" << std::endl;
            }
          }



	  Belos::BlockGmresSolMgr<ST,MV,OP> FSolver( FProblem, rcp(&myPL,false) );

          if(printproc) cout << "GMRES Solver Manager set"<<flush<<endl;
          Belos::ReturnType FsolverRet = FSolver.solve();

if (printproc) {
    if (FsolverRet == Belos::Converged) {
      cout << "Belos F converged."<<flush << std::endl;
    } else {
      cout << "Belos F did not converge." <<flush<< std::endl;
    }
  }



	  Teuchos::RCP<MV> FSol= FProblem->getLHS();
          x1=*FSol;
          y1=*Fb;

Epetra_Vector tempx1a(*x1(0));
     double npva; tempx1a.Norm2(&npva);
if(printproc) cout << "fsolnorm="<<npva<<endl;


Epetra_Vector tempy1a(*y1(0));
     double npya; tempy1a.Norm2(&npya);
 cout << "frhsnorm="<<npya<<endl;
//if(printproc) cout << "frhsnorm="<<npya<<endl;

        }

// Next apply B to x1 and store in Bx1
// We don't need to make B and DFinvBt Epetra Operators, only F and S, these other two can be applied directly as functions 

        Epetra_Vector bx1(*Y(0));
        bx1.PutScalar(0.0);
// Apply Analytic Jacobian function J to X and store result in Y, jacdata (@np1) holds data for point of linearization

//        precFunctionblock21(x1(0)->Values(),N, bx1.Values(), precdata);

//     double nB21; bx1.Norm2(&nB21);
//if(printproc) cout << "normB21="<<nB21<<endl;

// Then set Sproblos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(F,x,b) );


// Insert loop here to add -bx1 and b2 and put in Schurb
// b2 the second component of the X vector

        Epetra_MultiVector b2(X);
        for (int i=0; i<2*(N+1)/3;i++) b2[0][i] = 0.0;


	double nsa; b2(0)->Norm2(&nsa); 
if(printproc) cout << "Norm of RHS Schur A="<<nsa<<endl;


    for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]; // No L block
//   for (int i=2*(N+1)/3; i<N; i++) b2[0][i] = b2[0][i]-bx1[i]; // Lower block

	double nsb; b2(0)->Norm2(&nsb); 
if(printproc) cout << "Norm of RHS Schur B="<<nsb<<flush<<endl;

	Teuchos::RCP<const MV>Schurb=Teuchos::rcp ( new MV(b2));

if(printproc) cout << "Schur rhs set="<<flush<<endl;


        Teuchos::RCP<MV>Schurx=Teuchos::rcp ( new MV(Y));
        //We initialized Y to zero and now have based Schurx on Y //Schurx.PutScalar(0.0);

if(printproc) cout << "Schur solution initialized "<<flush<<endl;
        Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(S,Schurx,Schurb) );


if(printproc) cout << "Schur Problem initialized"<<flush<<endl;

cout<<"Schur initialized comm.mypid="<<comm.MyPID()<<endl<<flush;

        bool Sret = SchurProblem->setProblem(); 


if (printproc) {
    if (Sret == true) {
      cout << "Belos S Linear Problem Set"<<flush<< std::endl;
    } else {
      cout << "Error setting Belos S Linear Problem" <<flush<< std::endl;
    }
  }


        Belos::BlockGmresSolMgr<ST,MV,OP> SchurSolver( SchurProblem, rcp(&myPL,false) );
if(printproc) cout << "Schur GMRES Solver Set"<<flush<<endl;
        Belos::ReturnType SchursolverRet = SchurSolver.solve();



if (printproc) {
    if (SchursolverRet == Belos::Converged) {
      cout << "Belos Schur converged." <<flush<< std::endl;
    } else {
      cout << "Belos Schur did not converge." << flush<<std::endl;
    }
  }



	Teuchos::RCP<MV> SchurSol= SchurProblem->getLHS();
        Epetra_MultiVector x2(*SchurSol);



//Next apply dDinvBt to x1 and store in Bx1

//        Epetra_Vector dFinvBt(*Y(0));
//        dFinvBt.PutScalar(0.0);

//	precFunctionblock12(x2(0)->Values(),N, dFinvBt.Values(), precdata);

 //    double nBt; dFinvBt.Norm2(&nBt);
//if(printproc) cout << "normBt="<<nBt<<flush<<endl;



// Insert loop here to add x1 -(1/alpha)*dFinvBt*x2
// Since the vector is broken up into two components with zeroes padding the 
// 'non-active' entries we can build Y in one action
// We write down that formally this is being done:
// y1= x1 -(1/alpha)*dFinvBt*x2
// y2= (1/alpha)x2 
// 
// which is equivalentlly implemented as Y= x1 -(1/alpha)*dFinvBt*x2 + (1/alpha)x2

Epetra_Vector tempx1(*x1(0));
Epetra_Vector tempx2(*x2(0));
double alpha=1.0;
double alphainv=1.0/alpha;


/*
         for (int i=0; i<2*(N+1)/3; i++) {
             tempx1[i]=x1[0][i]-alphainv*dFinvBt[i];
             //tempx1[i]=x1[0][i]-0.0*dFinvBt[i];
          }
         for (int i=2*(N+1)/3;i<N; i++) {
             tempx1[i]=(alphainv*tempx2[i]);
          }
*/
        //Epetra_Vector Ix1(*X(0));
         //for (int i=2*(N+1)/3;i<N; i++) {Ix1[i]=0.0;}
    //for (int i=0; i<N; i++) tempx1[i] = (Ix1[i]+tempx2[i]);
//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]+tempx2[i]);


     double npvf; tempx1.Norm2(&npvf);
if(printproc) cout << "Diagonal components"<<flush<<endl;
if(printproc) cout << "factoredprecnormvec1="<<npvf<<flush<<endl;
double npvs; tempx2.Norm2(&npvs);
if(printproc) cout << "factoredprecnormvec2="<<npvs<<flush<<endl;


 for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]) +(alphainv*tempx2[i]); //Diagonal Block Preconditioner
//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-alphainv*dFinvBt[i]) +(alphainv*tempx2[i]); //Upper
 //  for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-0.0*dFinvBt[i]) +(alphainv*tempx2[i]);//Diagonal and also DU
//   for (int i=0; i<N; i++) tempx1[i] = tempx1[i]+(alphainv*tempx2[i]);//Diagonal and also DU


//    for (int i=0; i<N; i++) tempx1[i] = (tempx1[i]-dFinvBt[i]) +(tempx2[i]); //No G no alpha


     double npv; tempx1.Norm2(&npv);
if(printproc) cout << "factoredprecnorm="<<npv<<endl;

Y=tempx1;


return 0;
}
#endif








void Problem_Interface::runPreIterate (const NOX::Solver::Generic &solver)
{

//    if (printproc)   cout<<"Pre Iterate"<<endl;
 //   if (printproc)   cout<<"N="<<N<<endl;
//    if (comm.MyPID()==0) cout<<"comm.mypid="<<comm.MyPID()<<endl<<flush;


//  cout<<"Updating Preconditioner before nonlinear solve"<<endl;
//  precUpdateFunction(statevector,N,precdata);
//  precUpdateFunction(statevector,N,blackbox_res);
}


//-----------------------------------------------------------------------------

