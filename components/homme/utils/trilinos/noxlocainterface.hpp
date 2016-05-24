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
#include "precon_interface.hpp"
#include <NOX_Abstract_PrePostOperator.H>






typedef void (*residFnPtr)(double *, double *, int, void *);


class  Problem_Interface : 
	//public NOX::Epetra::Interface::Required,
	//LOCA 1
	public LOCA::Epetra::Interface::Required,
	public NOX::Epetra::Interface::Preconditioner,
	public NOX::Abstract::PrePostOperator,
	public Epetra_Operator,
        public NOX::Epetra::Interface::Jacobian

{
	public:


		Problem_Interface(int nelems, double* statevector_,
				const LOCA::ParameterVector& pVector_,
				const Epetra_Comm& comm_,
				void* blackbox_res, void* precdata,void * jacdata,
				void (*residualFunction)(double *, double *, int, void *),
				void (*precFunction)(double *,int,double*,void *),
				void (*jacFunction)(double *,int,double*,void *),
				void (*precUpdateFunction)(double *,int,void *),
                   		void (*getJacVector)(double *, int, void *));

		Problem_Interface(int nelems, double* statevector_,
				const LOCA::ParameterVector& pVector_,
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

/* interface for comparing two block preconditioner formulations */
                Problem_Interface(int nelems, double* statevector_,
                                const LOCA::ParameterVector& pVector_,
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



        	Problem_Interface(int nelems, double* statevector_,
				const LOCA::ParameterVector& pVector_,
				const Epetra_Comm& comm_,
				void* blackbox_res, void* precdata,
				void (*residualFunction)(double *, double *, int, void *),
				void (*precFunction)(double *,int,double*,void *),
				void (*precUpdateFunction)(double *,int,void *));

/*

		Problem_Interface(int nelems, double* statevector_,
				const LOCA::ParameterVector& pVector_,
				const Epetra_Comm& comm_,
				void* blackbox_res, void* blackbox_prec,
				void (*residualFunction)(double *, double *, int, void *),
				void (*precFunction)(double *,int,double*,void *),
				void (*precUpdateFunction)(double *,int,void *));
*/
		~Problem_Interface();

		//! Compute and return F
		bool computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag);
		bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac);


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


void printCommID();


		// 10 Methods for inheritance from Epetra_Operator
		// Only ApplyInverse is non-trivial -- first 9 satisfied here in header

		//STRAT5
		int SetUseTranspose(bool UseTranspose) {
			if (UseTranspose) {cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<endl; return -1;}
			else return 0;
		};
		double NormInf() const
		{ cerr<<"ERROR: No noxlocainterface::Apply"<<endl; return 1.0;};
		//const char* Label() const { return "noxlocainterface::user preconditioner";};
		const char* Label() const { return label.c_str();};
		bool UseTranspose() const { return false;};
		bool HasNormInf() const { return false;};
		const Epetra_Comm& Comm() const {return comm;};
		const Epetra_Map& OperatorDomainMap() const {return *globalMap;};
		const Epetra_Map& OperatorRangeMap() const {return *globalMap;};

		void resetBlackbox(void* blackbox_res_);
		void resetBlackbox(void* blackbox_res_,  void* precdata_);
		//  { blackbox_res=blackbox_res_; blackbox_prec=blackbox_prec_; }
		void resetBlackbox(void* blackbox_res_,  void* precdata_, void* jacdata_);
		void getJacVec(Epetra_Vector& y) const;
		void getbbVec(Epetra_Vector& y) const;

		//! Apply the preconditioner
		int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

	//Apply the Jacobian
		int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const;

//                 PrePostOperator (const Teuchos::RCP< NOX::Utils > &utils, Teuchos::ParameterList &solverOptionsSubList);
//  virtual ~PrePostOperator ();

		virtual void runPreIterate (const NOX::Solver::Generic &solver);

	private:


		int N;
		double* statevector;
		const Epetra_Comm& comm;
		LOCA::ParameterVector pVector;
		void* blackbox_res;
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
		Teuchos::RCP<Epetra_Map> globalMap;
		Teuchos::RCP<Epetra_Vector> solution;
                bool printproc;
     		Teuchos::RCP<Epetra_Operator>A;
     		Teuchos::RCP<Epetra_Operator>F;
     		Teuchos::RCP<Epetra_Operator>DFinvBt;
     		Teuchos::RCP<Epetra_Operator>B;
     		Teuchos::RCP<Epetra_Operator>S;




		std::string label;

};

#endif

