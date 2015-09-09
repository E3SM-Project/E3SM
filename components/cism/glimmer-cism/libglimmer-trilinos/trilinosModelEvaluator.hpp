#ifndef GLIMMER_MODELEVALUATOR_HPP
#define GLIMMER_MODELEVALUATOR_HPP

#include "Teuchos_RCP.hpp"
#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Operator.h"

using Teuchos::RCP;


class trilinosModelEvaluator : public EpetraExt::ModelEvaluator {
public:
  
  
  trilinosModelEvaluator(int N_,
                         double* statevector,
                         const Epetra_Comm& comm_,
                         void* blackbox_res);
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

private:
  // Solution vector and map
  int N;
  RCP<Epetra_Map> xMap;
  RCP<Epetra_Vector> xVec;
  const Epetra_Comm& comm;
  void* blackbox_res;
  RCP<Epetra_Operator> precOp;

  RCP<Epetra_Map> pMap;
  RCP<Epetra_Vector> pVec;
};


class trilinosPreconditioner : public Epetra_Operator {
  
public:
  // Preconditioner as  Epetra_Operator required methods
  
  trilinosPreconditioner(int N, RCP<Epetra_Vector> xVec, RCP<Epetra_Map> xMap,
                         void* blackbox_res);

  int ApplyInverse(const Epetra_MultiVector& V, Epetra_MultiVector& Y) const;

  // Trivial implemetations
  int SetUseTranspose(bool UseTranspose) { TEUCHOS_TEST_FOR_EXCEPT(UseTranspose); return 0;};
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
    { throw "No Apply() in TrilinosPreconditioner";};
  double NormInf() const { throw "NO NormInf Implemented in trilinosPrecon";};
  const char* Label () const { return "trilinosPrec"; };
  bool UseTranspose() const { return false; };
  bool HasNormInf() const { return false; };
  const Epetra_Comm & Comm() const { return xMap->Comm();};
  const Epetra_Map& OperatorDomainMap () const { return *xMap;};
  const Epetra_Map& OperatorRangeMap  () const { return *xMap;};
  
private:
  int N;
  RCP<Epetra_Vector> xVec;
  RCP<Epetra_Map> xMap;
  void* blackbox_res;
};

#endif

