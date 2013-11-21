#ifndef HOMME_MODELEVALUATOR_HPP
#define HOMME_MODELEVALUATOR_HPP

#include "Teuchos_RCP.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"

#include <iostream>
#include "LOCA_Parameter_Vector.H"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "LOCA_Epetra.H"
#include <NOX_Abstract_PrePostOperator.H>


#include "BelosLinearProblem.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp" 


using Teuchos::RCP;
using Teuchos::ParameterList;
typedef double ST;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;
typedef void (*residFnPtr)(double *, double *, int, void *);

class hommePreconditionerBase;



class trilinosModelEvaluator : public EpetraExt::ModelEvaluator {
public:
  
  // 1. Constructors. Different interfaces for different preconditioners

  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precUpdateFunction)(double *,int,void *));

  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,void * jacdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precFunction)(double *,int,double*,void *),
        void (*jacFunction)(double *,int,double*,void *),
        void (*precUpdateFunction)(double *,int,void *),
        void (*getJacVector)(double *, int, void *));

  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,void * jacdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precFunctionblock11)(double *,int,double*,void *),
        void (*precFunctionblock12)(double *,int,double*,void *),
        void (*precFunctionblock21)(double *,int,double*,void *),
        void (*precFunctionblock22)(double *,int,double*,void *),
        void (*jacFunction)(double *,int,double*,void *),
        void (*precUpdateFunction)(double *,int,void *),
        void (*getJacVector)(double *, int, void *));


  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,void * jacdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precFunctionblock11)(double *,int,double*,void *),
        void (*precFunctionblock12)(double *,int,double*,int, void *),
        void (*precFunctionblock21)(double *,int,double*,int, void *),
        void (*precFunctionblock22)(double *,int,double*,void *),
        void (*jacFunction)(double *,int,double*,void *),
        void (*precUpdateFunction)(double *,int,void *),
        void (*getJacVector)(double *, int, void *),
        const RCP<ParameterList>&  FSolvePL_, 
        const RCP<ParameterList>&  SchurSolvePL_, 
        int *FTotalIt, 
        int *SchurTotalIt);


  /* interface for comparing two block preconditioner formulations */
  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,void * jacdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precFunctionblock11)(double *,int,double*,void *),
        void (*precFunctionblock12)(double *,int,double*,void *),
        void (*precFunctionblock21)(double *,int,double*,void *),
        void (*precFunctionblock22)(double *,int,double*,void *),
        void (*auxprecFunctionblock11)(double *,int,double*,void *),
        void (*auxprecFunctionblock12)(double *,int,double*,void *),
        void (*auxprecFunctionblock21)(double *,int,double*,void *),
        void (*auxprecFunctionblock22)(double *,int,double*,void *),
        void (*jacFunction)(double *,int,double*,void *),
        void (*precUpdateFunction)(double *,int,void *),
        void (*getJacVector)(double *, int, void *));

  trilinosModelEvaluator(int nelems, double* statevector_,
        const Epetra_Comm& comm_,
        void* blackbox_res, void* precdata,
        void (*residualFunction)(double *, double *, int, void *),
        void (*precFunction)(double *,int,double*,void *),
        void (*precUpdateFunction)(double *,int,void *));

    ~trilinosModelEvaluator();
 
  // 2. Required Methods from ModelEvaluator Base Class
  //@{
  
  //! Return solution vector map
  RCP<const Epetra_Map> get_x_map() const;
  
  //! Return residual vector map
  RCP<const Epetra_Map> get_f_map() const;
  
  //! Return initial solution and x_dot init
  RCP<const Epetra_Vector> get_x_init() const;

  RCP<EpetraExt::ModelEvaluator::Preconditioner> create_WPrec() const;

  //! Parameter setting functions for LOCA continuation
  RCP<const Epetra_Map> get_p_map(int l) const;
  RCP<const Epetra_Vector> get_p_init(int l) const;
  RCP<const  Teuchos::Array<std::string> >  get_p_names(int l) const;
  
  //! Create InArgs
  InArgs createInArgs() const;
  
  //! Create OutArgs
  OutArgs createOutArgs() const;
  
  //! Reset State
  //    void ResetState(double *statevector,void* blackbox_res_);
  
  //! Evaluate model on InArgs
  void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;
  //@}
 
  // 3. Extra interface functions specific to this implementation
  void resetInitialGuess(double* statevector_);
  void resetBlackbox(void* blackbox_res_);
  void resetBlackbox(void* blackbox_res_,  void* precdata_);
  void resetBlackbox(void* blackbox_res_,  void* precdata_, void* jacdata_);


  // 4. Utility functions not accessible from outside this class.
private:
  void initialize(double* statevector);


private:
  // Solution vector and map
  int N;
  const Epetra_Comm& comm;
  RCP<Epetra_Map> xMap;
  RCP<Epetra_Vector> xVec;
  RCP<Epetra_Map> pMap;
  RCP<Epetra_Vector> pVec;
  bool printproc;

  void* blackbox_res;
  RCP<hommePreconditionerBase> precOp;

  void *precdata;
  void *jacdata;
  void (*residualFunction)(double *, double *, int, void *);
  void (*precFunction)(double *,int,double*,void *);
  void (*precFunctionblock11)(double *,int,double*,void *);
  void (*precFunctionblock12)(double *,int,double*,void *);
  void (*precFunctionblock21)(double *,int,double*,void *);
  void (*precFunctionblock22)(double *,int,double*,void *);


  void (*auxprecFunctionblock11)(double *,int,double*,void *);
  void (*auxprecFunctionblock12)(double *,int,double*,void *);
  void (*auxprecFunctionblock21)(double *,int,double*,void *);
  void (*auxprecFunctionblock22)(double *,int,double*,void *);
  void (*jacFunction)(double *,int,double*,void *);
  void (*precUpdateFunction)(double *,int,void *);
  void (*getJacVector)(double *,int,void *);
  Teuchos::RCP<Epetra_Operator>A;
  Teuchos::RCP<Epetra_Operator>F;
  Teuchos::RCP<Epetra_Operator>DFinvBt;
  Teuchos::RCP<Epetra_Operator>B;
  Teuchos::RCP<Epetra_Operator>S;


};


/*******************************************************************************/
/***** Base Class under all Preconditioners, sets defaults *********************/
/*******************************************************************************/
class hommePreconditionerBase : public Epetra_Operator {
  
public:
  // Preconditioner as  Epetra_Operator required methods
  
  hommePreconditionerBase(RCP<Epetra_Map> xMap) : xMap_(xMap) {};

  // Methods that MUST be implemented for preconditioner
  virtual int computePreconditioner(RCP<const Epetra_Vector> xVec, void* precdata_) = 0;
  virtual int ApplyInverse(const Epetra_MultiVector& V, Epetra_MultiVector& Y) const = 0;

  // Trivial implemetations set in this base class
  int SetUseTranspose(bool UseTranspose) { TEUCHOS_TEST_FOR_EXCEPT(UseTranspose); return 0;};
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
    { throw "No Apply() in TrilinosPreconditioner";};
  double NormInf() const { throw "NO NormInf Implemented in trilinosPrecon";};
  const char* Label () const { return "trilinosPrec"; };
  bool UseTranspose() const { return false; };
  bool HasNormInf() const { return false; };
  const Epetra_Comm & Comm() const { return xMap_->Comm();};
  const Epetra_Map& OperatorDomainMap () const { return *xMap_;};
  const Epetra_Map& OperatorRangeMap  () const { return *xMap_;};
  
private:
  hommePreconditionerBase(); //NOT ACCESSIBLE, NEED TO USE OTHER CONSTRUCTOR};
  RCP<Epetra_Map> xMap_;
};

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
   // Add Preconditioners, inheriting from hommePreconditionerBase to get
   // default implementations of most required methods of an Epetra_Operator
   // above. Only need to implement constructor, computePreconditioner,
   // and ApplyInverse. 
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
class identityPreconditioner : public hommePreconditionerBase {
public:
  identityPreconditioner(int N, RCP<Epetra_Vector> xVec, RCP<Epetra_Map> xMap,
                         void* blackbox_res, void* precdata,
                         void (*precUpdateFunction)(double *,int,void *) );
  int ApplyInverse(const Epetra_MultiVector& V, Epetra_MultiVector& Y) const;
  int computePreconditioner(RCP<const Epetra_Vector> xVec,  void* precdata_);

private:
  int N;
  RCP<Epetra_Vector> xVec;
  RCP<Epetra_Map> xMap;
  void* blackbox_res;
  void* precdata;
  void (*precUpdateFunction)(double *,int,void *);
};
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
class simplePreconditioner : public hommePreconditionerBase {
public:
  simplePreconditioner(int N, RCP<Epetra_Vector> xVec, RCP<Epetra_Map> xMap,
                       void* blackbox_res, void* precdata,
                       void (*precFunctionblock11)(double *,int,double*,void *),
                       void (*precFunctionblock12)(double *,int,double*,int, void *),
                       void (*precFunctionblock21)(double *,int,double*,int, void *),
                       void (*precFunctionblock22)(double *,int,double*,void *),
                       void (*precUpdateFunction)(double *,int,void *),
		       const RCP<ParameterList>&  FSolvePL_,
		       const RCP<ParameterList>&  SchurSolvePL_,
		       int* FTotalIt_,
		       int* SchurTotalIt_ );

  int ApplyInverse(const Epetra_MultiVector& V, Epetra_MultiVector& Y) const;
  int computePreconditioner(RCP<const Epetra_Vector> xVec, void* precdata_);

private:
  int N;
  bool printproc;
  RCP<Epetra_Vector> xVec;
  RCP<Epetra_Map> xMap;
  void* blackbox_res;
  void* precdata;
  void (*precFunctionblock11)(double *,int,double*,void *);
  void (*precFunctionblock12)(double *,int,double*,int, void *);
  void (*precFunctionblock21)(double *,int,double*,int, void *);
  void (*precFunctionblock22)(double *,int,double*,void *);
  void (*precUpdateFunction)(double *,int,void *);
  RCP<Epetra_Operator>F;
  RCP<Epetra_Operator>S;

  int *FTotalIt;
  int *SchurTotalIt;


  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > FSolver;
  Teuchos::RCP< Belos::SolverManager<ST,MV,OP> > SchurSolver;

  Teuchos::RCP<Epetra_Map> UVMap;
  Teuchos::RCP<Epetra_Map> HMap;

  Teuchos::RCP<ParameterList>  FSolvePL;
  Teuchos::RCP<ParameterList>  SchurSolvePL;

  Teuchos::RCP<Epetra_Vector> workvector4;
  Teuchos::RCP<Epetra_Vector> dFinvBt;
  Teuchos::RCP<Epetra_Vector> bx1;

  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > FProblem;
  Teuchos::RCP< Belos::LinearProblem<ST,MV,OP> > SchurProblem;

  Teuchos::RCP<Epetra_Vector> Fb; //uv workvector
  Teuchos::RCP<Epetra_Vector> Fx; //uv workvector

  Teuchos::RCP<Epetra_Vector> Schurb; //height workvector
  Teuchos::RCP<Epetra_Vector> Schurx; //height workvector

};

#endif
