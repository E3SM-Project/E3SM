//STRAT1
#include "NOX_Epetra_LinearSystem_Stratimikos.H"

// Trilinos Objects
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"


#include "Epetra_LinearProblem.h"
#include <BelosLinearProblem.hpp>
#include "BelosConfigDefs.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "BelosLinearProblem.hpp"

#include "Epetra_InvOperator.h"
#include "Epetra_Util.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"

// Required for reading and writing parameter lists from xml format
#include "Teuchos_XMLParameterListHelpers.hpp"

//#define DEBUG
using namespace std;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

typedef double                            ST;
typedef Teuchos::ScalarTraits<ST>        SCT;
typedef SCT::magnitudeType                MT;
typedef Epetra_MultiVector                MV;
typedef Epetra_Operator                   OP;
typedef Belos::MultiVecTraits<ST,MV>     MVT;
typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

// Define some variables that are Global to the file that
// allow Belos to be called in 3 stages (init, solve, finalize)
//static Teuchos::RCP<Epetra_MpiComm> Comm;
static Teuchos::RCP<Epetra_Comm> Comm;
static bool printproc;
static Teuchos::RCP<Teuchos::ParameterList> paramList;
static Teuchos::RCP<Teuchos::ParameterList> HelmSolvePL;
static Teuchos::RCP<Teuchos::ParameterList> MLList;
static int* HelmTotalIt;
static int HelmTotalIts;
static int OutputStep;


static int *my_GlobalIDs;
static Teuchos::RCP<Epetra_Map> ParMatrixMap;

static Teuchos::RCP<Epetra_Map> MatrixMap; 
//static Teuchos::RCP<Epetra_Export::Epetra_Export> exporter;
//static Teuchos::RCP<Epetra_Import::Epetra_Import> importer;

static Teuchos::RCP<Epetra_Export> exporter;
static Teuchos::RCP<Epetra_Import> importer;
static Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;
static Teuchos::RCP<Epetra_InvOperator> MLPI;

static int *bincount;
static int **bin;

static Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > HelmProblem;
static Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > HelmSolver;
static Teuchos::RCP<Epetra_Vector>   HelmSolution; 
static Teuchos::RCP<Epetra_Vector>   ParHelmSolution; 

static Teuchos::RCP<Epetra_Vector>   ParHelmRHS; 
static Teuchos::RCP<Epetra_Vector>   HelmRHS; 
static Teuchos::RCP<Epetra_FECrsMatrix> HelmMatrix;
static Teuchos::RCP<Epetra_FECrsMatrix> ParHelmMatrix;

static bool StaticProfile=false;
static int N; //number of points on each processor (unique to processor but not across processors) 

extern "C" {
/*Defining external functions for finite difference Jacobian*/
  void homme_globalIDs(int, int* ,void *);
  void helm_mat(int, int, double *, int *, void *);
  void helm_map(int, int, int *, void *);
}

extern "C" { 

/*
This routine builds the maps for the Helmholtz Matrix assembly and rhs assembly

the data from HOMME is configured so that DOF's can be represented multiple times
we first remove these multiplities and record the indices of the vector corresponding
to the unique index in the "bin" structure. A point can be represented at most on 4 non-unique indices
the first index of bin[i] stores the unique index on the local processor, the second index stores the
non-index indices. bincount[i] counts the number of bins (multiplicties) for a given in index. The local rhs vector
will need to be scaled by this number in the assembly process.

The second level assembly is also performed to remove non-unique entries across processors. This is done using the ParMatrixMap and the inporter/exporter between the ParMatrixMap and the local MatrixMap 

*/
void BuildMaps(int* nelems, int* nets, int* nete, int* np, int* nlev, void* data){

     void (*get_globalIDs)(int, int *, void *) = homme_globalIDs;
     try {

     my_GlobalIDs=new(int[*nelems]); 
           
     Comm=Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

     if (Comm->MyPID()==0) printproc=true;
     else   printproc=false;
	   

     int globalids[*nelems];

     get_globalIDs(*nelems,my_GlobalIDs,data);
	   
     for(int i=0;i<*nelems;i++){
         globalids[i]=my_GlobalIDs[i];
     }
	   

     //sorting ids

     std::sort(&my_GlobalIDs[0],&my_GlobalIDs[0]+ (*nelems));

     int *end=unique(my_GlobalIDs,&my_GlobalIDs[0]+(*nelems));
	   
     N=end-(&my_GlobalIDs[0]);

     bin=new int*[N];
     for (int i=0;i<N;i++){
	  bin[i]=new int[5];
     }

     for(int i=0;i<N;i++){
         for (int j=0;j<5;j++){
   	      bin[i][j]=*nelems+1;//initialize to out index out of bounds
	 }
     }

     for(int i=0;i<N;i++){
         bin[i][0]=my_GlobalIDs[i];
     }

     int avail;
     bincount=new int[N];

     for(int i=0;i<N;i++){
         bincount[i]=0;
         avail=1;
	 for(int j=0;j<*nelems;j++){
	     if(globalids[j]==bin[i][0]){
	        bin[i][avail]=j;
	        avail++;
		bincount[i]++;
	     }
	 }
     }


     MatrixMap = Teuchos::rcp(new Epetra_Map(-1, N, my_GlobalIDs, 0, *Comm));

     bool high_rank_proc_owns_shared=false;

     ParMatrixMap = Teuchos::rcp(new Epetra_Map(
		     Epetra_Util::Create_OneToOne_Map(*MatrixMap,high_rank_proc_owns_shared) ));

     exporter = Teuchos::rcp(new Epetra_Export(*MatrixMap, *ParMatrixMap));
     importer = Teuchos::rcp(new Epetra_Import(*MatrixMap, *ParMatrixMap));

     } //end try

    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception in BuildMaps " << endl;
    }
}//end BuildMaps



void BuildMatrix(int* nelems, int* nets, int* nete, int* np, int* nlev,  void* data ){
     
    void (*get_HelmElementMat)(int, int, double *,int *, void *)=helm_mat;

 try{

     int matsize=(*np)*(*np);
     double Aloc[matsize][matsize];
     double Zeroloc[matsize][matsize];
     int Index[matsize];
     int ColumnMajor=3;
     int ndx=my_GlobalIDs[N-1];

     Epetra_DataAccess copy = ::Copy;
     const std::size_t approximate_indices_per_row = (*nelems);
    // HelmMatrix=Teuchos::null;
    // ParHelmMatrix=Teuchos::null;

     HelmMatrix = Teuchos::rcp(new Epetra_FECrsMatrix(copy, *MatrixMap,approximate_indices_per_row,StaticProfile));
     ParHelmMatrix = Teuchos::rcp(new Epetra_FECrsMatrix(copy, *ParMatrixMap,approximate_indices_per_row,StaticProfile));


     for(int i=0;i<matsize;i++){
 	 for(int j=0;j<matsize;j++){
  	     Zeroloc[i][j]=0.0;
	 }
     }




//Initialize Matrix to Zero
     for (int ie=*nets;ie<*nete+1;ie++){
          get_HelmElementMat(ie,matsize,&Aloc[0][0],&Index[0],data);
          HelmMatrix->InsertGlobalValues(matsize,&Index[0],matsize,&Index[0],&Zeroloc[0][0],ColumnMajor);
     }

//Sum Matrix entries Element by Element 
     for (int ie=*nets;ie<*nete+1;ie++){
          get_HelmElementMat(ie,matsize,&Aloc[0][0],&Index[0],data);
          HelmMatrix->SumIntoGlobalValues(matsize,&Index[0],matsize,&Index[0],&Aloc[0][0],ColumnMajor);
     }

     HelmMatrix->FillComplete();
     ParHelmMatrix->Export(*HelmMatrix, *exporter, Add);
     bool callFillComplete=true;
     ParHelmMatrix->GlobalAssemble(callFillComplete);

  } //end try

    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception in BuildMatrix " << endl;
    }
}//end BuildMatrix





void BuildRHS(int* nets, int* nete, int* np, int* nlev, double* rhsVector, void* data ){

    void (*get_HelmMap)(int, int, int *, void *)=helm_map;

    try{

     int matsize=(*np)*(*np);
     double rhsLoc[matsize];
     int rhsIndex[matsize];
     int ndx=my_GlobalIDs[N-1];
     double rhsGlob[ndx];
     int rhsGIndex[ndx];

     for(int i=0;i<matsize;i++){
 	 rhsIndex[i]=0;
 	 rhsLoc[i]=0.0;
     }

     for(int i=0;i<ndx;i++){
	 rhsGlob[i]=0;
	 rhsGIndex[i]=0;
     }

     for (int ie=*nets;ie<*nete+1;ie++){
          get_HelmMap(ie,matsize,&rhsIndex[0],data);
      //  for (int lev=0;lev<1;lev++){
          int lndx=0;

   	  for (int i=0;i<*np;i++){
	       for (int j=0;j<*np;j++){
	     	    rhsLoc[i*(*np)+j]=rhsVector[ (ie-*nets)*(*np)*(*np) +i*(*np)+j] ;
		    rhsGIndex[rhsIndex[lndx]-1]=rhsIndex[lndx];
		    rhsGlob[rhsIndex[lndx]-1]+=rhsLoc[i*(*np)+j];
		    lndx++;
	       }

          }
        //}

     }

     HelmRHS = Teuchos::rcp(new Epetra_Vector(*MatrixMap,StaticProfile));
     ParHelmRHS = Teuchos::rcp(new Epetra_Vector(*ParMatrixMap,StaticProfile));
     HelmRHS->PutScalar(0.0);
     ParHelmRHS->PutScalar(0.0);

     HelmRHS->SumIntoGlobalValues(ndx,&rhsGlob[0],&rhsGIndex[0]);

     for (int i=0;i<N;i++){
 	  (*HelmRHS)[i]=(*HelmRHS)[i]/bincount[i];
     }

     ParHelmRHS->Export(*HelmRHS, *exporter,Insert);

  } //end try

    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception in BuildRHS " << endl;
    }
}//end BuildRHS



void SetProblem(){

   try{
     HelmSolution = Teuchos::rcp(new Epetra_Vector(*MatrixMap,StaticProfile ));
     ParHelmSolution = Teuchos::rcp(new Epetra_Vector(*ParMatrixMap,StaticProfile));

     HelmProblem= Teuchos::rcp( new Belos::LinearProblem<ST,MV,OP>(ParHelmMatrix,ParHelmSolution,ParHelmRHS) );
     HelmProblem->setOperator( ParHelmMatrix);
     HelmProblem->setLHS( ParHelmSolution);
     HelmProblem->setRHS( ParHelmRHS);

// Begin Preconditioner Set
    double factor=4.0/3.0;

     MLList = Teuchos::rcp(new Teuchos::ParameterList);
     ML_Epetra::SetDefaults("SA",*MLList);
     //ML_Epetra::SetDefaults("DD",*MLList);
     //ML_Epetra::SetDefaults("DD-ML",*MLList);
      
     MLPrec=Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*ParHelmMatrix,*MLList,true));

     bool Helmret = HelmProblem->setProblem();
     if (Helmret != true && printproc) cout<<"Error setting problem"<<flush<<endl;

     const bool CheckFiltering = false;
     int comp=MLPrec->ComputePreconditioner();

     MLPI=Teuchos::rcp(new Epetra_InvOperator(MLPrec.get()));

     HelmProblem->setRightPrec(MLPI);

     comp=MLPrec->IsPreconditionerComputed();

     const int cycles=1;
     comp= MLPrec->AnalyzeCycle(cycles);

     // MLPrec->PrintUnused(0);
     //
     //
     //
    


     paramList = Teuchos::rcp(new Teuchos::ParameterList);
     Teuchos::updateParametersFromXmlFile("helmbelos.xml", paramList.ptr());
     HelmSolvePL = Teuchos::rcp(new Teuchos::ParameterList);
     *HelmSolvePL = paramList->sublist("HelmSolvePL");
     HelmTotalIts=0;
     HelmTotalIt=&HelmTotalIts;	   

     HelmSolver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>( HelmProblem, HelmSolvePL ) );
     
 
     } //end try

    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception in SetProblem  " << endl;
    }
}//end SetProblem



void HelmholtzSolve(int *nelems, double* statevector){

    try{

     HelmSolution->PutScalar(0.0);
     ParHelmSolution->PutScalar(0.0);

     HelmProblem->setOperator( ParHelmMatrix);
     HelmProblem->setLHS( ParHelmSolution);
     HelmProblem->setRHS( ParHelmRHS);

//     HelmProblem->Reset( ParHelmSolution, ParHelmRHS );

     HelmSolver->reset( Belos::Problem );

     Belos::ReturnType HelmSolveStatus = HelmSolver->solve();

     HelmTotalIts = HelmSolver->getNumIters();
     HelmTotalIt=&HelmTotalIts;

     HelmSolution->Import(*ParHelmSolution, *importer, Insert);

     if(printproc)cout<<"Trilinos Computed Helmholtz Solve in "<<*HelmTotalIt<<" iterations"<<flush<<endl;

     int ndx=my_GlobalIDs[N-1];
     double solGlob[ndx];

     for (int j=1;j<5;j++){
     	for( int i=0;i<N;i++){
	      if(bin[i][j]!=*nelems+1){
		      solGlob[bin[i][j]]=(*HelmSolution)[i];
		      }
	      }
      }

      for(int i=0;i<*nelems;i++){
          statevector[i]=solGlob[i];
      }

    
      HelmRHS=Teuchos::null;
      ParHelmRHS=Teuchos::null;


   } //end try


    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (const char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception in HelmholtzSolve " << endl;
    }
  }//end HelmholtzSolve


  void belosfinish()
  {
    delete bincount;
    delete bin;
    delete my_GlobalIDs;
    Comm=Teuchos::null;
    MatrixMap=Teuchos::null;
    ParMatrixMap=Teuchos::null;
    exporter=Teuchos::null;
    importer=Teuchos::null;
    HelmMatrix = Teuchos::null;
    ParHelmMatrix = Teuchos::null;
    HelmRHS = Teuchos::null;
    ParHelmRHS = Teuchos::null;
    HelmSolution = Teuchos::null;
    ParHelmSolution = Teuchos::null;
    HelmProblem = Teuchos::null;
    MLList = Teuchos::null;
    MLPrec=Teuchos::null;
    MLPI=Teuchos::null;
    paramList = Teuchos::null;
    HelmSolvePL = Teuchos::null;
    HelmSolver = Teuchos::null;
  }

} //extern "C"
