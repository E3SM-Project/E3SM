//-----------------------------------------------------------------------------
#ifndef noxlocainterface_H
#define noxlocainterface_H

// Interface to the NLS_PetraGroup to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>

#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "LOCA_Epetra.H"
#include "LOCA_Parameter_Vector.H"

typedef void (*residFnPtr)(double *, double *, int, void *);


class  Problem_Interface : 
  public LOCA::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Preconditioner,
  public Epetra_Operator
{
public:
  Problem_Interface(int* nelems, double* statevector,
                    const LOCA::ParameterVector& pVector_,
                    const Epetra_Comm& comm_,
                    void* blackbox_res, void* blackbox_prec,
                    void (*residualFunction)(double *, double *, int, void *),
//                    void (*precFunction)(double *,double *,int,double*,void *));
                    void (*precFunction)(double *,double *,int,double*,void *,void *));
  ~Problem_Interface();
  
  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag);

  //! Set a parameter in the user's code.
  void setParameters(const LOCA::ParameterVector& params);

  //! Print solution to output file
  virtual void printSolution(const Epetra_Vector& x, double conParam);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 

  Teuchos::RCP<Epetra_Vector> getVector() const;

  // 1 Method for inheritance from NOX::Epetra::Interface::Preconditioner
  // Compute preconditioner \f$M\f$.
  virtual bool computePreconditioner(const Epetra_Vector& x,
                                     Epetra_Operator& Prec,
                                     Teuchos::ParameterList* p = 0);



  // 10 Methods for inheritance from Epetra_Operator
  // Only ApplyInverse is non-trivial -- first 9 satisfied here in header

//STRAT5
    int SetUseTranspose(bool UseTranspose) {
       if (UseTranspose) {cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<endl; return -1;}
       else return 0;
    };
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
        { cerr<<"ERROR: No noxlocainterface::Apply"<<endl; return -1;};
    double NormInf() const
        { cerr<<"ERROR: No noxlocainterface::Apply"<<endl; return 1.0;};
    const char* Label() const { return "noxlocainterface::user preconditioner";};
    bool UseTranspose() const { return false;};
    bool HasNormInf() const { return false;};
    const Epetra_Comm& Comm() const {return comm;};
    const Epetra_Map& OperatorDomainMap() const {return *globalMap;};
    const Epetra_Map& OperatorRangeMap() const {return *globalMap;};

    void resetBlackbox(void* blackbox_res_,  void* blackbox_prec_) {
       blackbox_res=blackbox_res_; blackbox_prec=blackbox_prec_; }

    //! Apply the preconditioner
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


  private:

    int N;
    const Epetra_Comm& comm;
    Teuchos::RCP<Epetra_Vector> solution;
    //Teuchos::RCP<Epetra_Vector> global_solution;
    //Teuchos::RCP<Epetra_Import> global_importer;
    Teuchos::RCP<Epetra_Map> globalMap;
    LOCA::ParameterVector pVector;
    void* blackbox_res;
    void* blackbox_prec;
    void (*residualFunction)(double *, double *, int, void *);
    void (*precFunction)(double *, double *, int, double*, void *, void *);
//    void (*precFunction)(double *, double *, int, double*, void *);  // SI precon

};

#endif

