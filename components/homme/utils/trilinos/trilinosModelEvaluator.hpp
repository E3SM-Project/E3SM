#ifndef HOMME_MODELEVALUATOR_HPP
#define HOMME_MODELEVALUATOR_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Operator.h"
#include "Epetra_Util.h"

#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"

#include "EpetraExt_ModelEvaluator.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp" 

#include <iostream>


using Teuchos::RCP;
using Teuchos::ParameterList;

typedef double ST;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;

class hommePreconditionerBase;
class hommeJacobianVectorProduct;


class trilinosModelEvaluator : public EpetraExt::ModelEvaluator {
public:
  
  // 1. Constructor
  trilinosModelEvaluator(int np_, int nlev_, int nelemd_,
			 double* StateVector_, void* StateData_, 
			 const Epetra_Comm& comm_,
			 void (*residualFunction)(double*, double*, int, void*),
			 void (*jacFunction)(double*, int, double*, void*),
			 void (*jacUpdateFunction)(double*, int, void*),
			 void (*precUpdateFunction)(double*, int, void*),
			 const RCP<ParameterList>&  ASolvePL_, 
			 const RCP<ParameterList>&  SSolvePL_, 
			 const RCP<ParameterList>&  PSolvePL_, 
			 int* ATotalIt, 
			 int* STotalIt,
			 int* PTotalIt);

  /* Destructor */
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

  // analytic Jacobian vector product
  RCP<Epetra_Operator> create_W() const;

  //! Parameter setting functions for LOCA continuation
  RCP<const Epetra_Map> get_p_map(int l) const;
  RCP<const Epetra_Vector> get_p_init(int l) const;
  RCP<const  Teuchos::Array<std::string> >  get_p_names(int l) const;
  
  //! Create InArgs
  InArgs createInArgs() const;
  
  //! Create OutArgs
  OutArgs createOutArgs() const;
  
  //! Reset State
  //    void ResetState(double *StateVector, void* StateData_);
  
  //! Evaluate model on InArgs
  void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;
  //@}
 
  // 3. Extra interface functions specific to this implementation
  void resetInitialGuess(double* StateVector_);
  void resetStateData(void* StateData_);

  // 4. Utility functions not accessible from outside this class.
private:
  void initialize(double* StateVector_);


private:
  int nState; // local state vector length

  const Epetra_Comm& comm; // MPI communicator

  RCP<Epetra_Map> xMap;    // state vector map
  RCP<Epetra_Vector> xVec; // state vector

  RCP<Epetra_Map> pMap;
  RCP<Epetra_Vector> pVec;

  bool printproc; // flag for output process

  void* StateData; // homme element derived data type

  // function pointer to residual function
  void (*residualFunction)(double*, double*, int, void*);

  // preconditioner operator
  RCP<hommePreconditionerBase> precOp;

  // Jacobian operator
  RCP<hommeJacobianVectorProduct> jacOp;
};



/*******************************************************************************/
/************ Epetra Operator for Jacobian-Vector Product  *********************/
/*******************************************************************************/
class hommeJacobianVectorProduct : public Epetra_Operator {
  
public:
  // Epetra_Operator required methods
  
  // Constructor, called once per time step (or maybe once per run)
  hommeJacobianVectorProduct(int nState_, 
			     RCP<Epetra_Vector> xVec_, 
			     RCP<Epetra_Map> xMap_,
			     void* StateData_, 
			     void (*jacFunction)(double*, int, double*, void*),                
			     void (*jacUpdateFunction)(double*, int, void*));

  // compute Jacobian called once per Newton step
  virtual int computeJacobian(RCP<const Epetra_Vector> xVec, void* StateData_);

  // Apply Jacobian called once per Krylov iter
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  // Trivial implemetations of other required functions set here
  int SetUseTranspose(bool UseTranspose) { TEUCHOS_TEST_FOR_EXCEPT(UseTranspose); return 0;};
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { throw "No ApplyInverse() in JacobianVectorProduct";};
  double NormInf() const { throw "NO NormInf Implemented in JacobianVectorProduct";};
  const char* Label () const { return "JacobianVectorProduct"; };
  bool UseTranspose() const { return false; };
  bool HasNormInf() const { return false; };
  const Epetra_Comm & Comm() const { return xMap->Comm();};
  const Epetra_Map& OperatorDomainMap () const { return *xMap;};
  const Epetra_Map& OperatorRangeMap  () const { return *xMap;};
  
private:
  bool printproc; // flag for output process

  int nState; // local state vector length

  RCP<Epetra_Map> xMap;    // state vector map
  RCP<Epetra_Vector> xVec; // state vector

  void* StateData;

  void (*jacFunction)(double*, int, double*, void*);
  void (*jacUpdateFunction)(double*, int, void*);
};



/*******************************************************************************/
// Base Class under all Preconditioners, sets defaults
//
// Add Preconditioners, inheriting from hommePreconditionerBase to get
// default implementations of most required methods of an Epetra_Operator.
//
// Only need to implement constructor, computePreconditioner, and ApplyInverse.
//
// All preconditioners must have xMap_ in constructor's argument list and the
// following line to build the Base class.
//     : hommePreconditionerBase(xMap_), //Required Base Class construction
/*******************************************************************************/
class hommePreconditionerBase : public Epetra_Operator {
  
public:
  // Preconditioner as Epetra_Operator required methods
  
  hommePreconditionerBase(RCP<Epetra_Map> xMap) : xMap_(xMap) {};

  // Methods that MUST be implemented for preconditioner
  virtual int computePreconditioner(RCP<const Epetra_Vector> xVec, void* StateData_) = 0;

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

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
  hommePreconditionerBase(); // NOT ACCESSIBLE, NEED TO USE OTHER CONSTRUCTOR
  RCP<Epetra_Map> xMap_;     // state vector map
};

#endif
