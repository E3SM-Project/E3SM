//-----------------------------------------------------------------------------
#ifndef block_precon_interface_H
#define block_precon_interface_H

// Interface to the NLS_PetraGroup to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>

#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "LOCA_Epetra.H"
#include "LOCA_Parameter_Vector.H"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"













class  Block_Precon_Interface : 
	public Epetra_Operator
{
	public:
		Block_Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_,
				void* prec_data_,
				void (*precFunction_)(double *,int,double*,void *));

		Block_Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_,
				void* prec_data_,
				void (*precFunctionblock11_)(double *,int,double*,void *),
				void (*precFunctionblock12_)(double *,int,double*,void *),
				void (*precFunctionblock21_)(double *,int,double*,void *),
				void (*precFunctionblock22_)(double *,int,double*,void *));

	//	Precon_Interface(const Epetra_Comm& comm_,int nelems,Teuchos::RCP<Epetra_Map> gmap);


		~Block_Precon_Interface();

		// 10 Methods for inheritance from Epetra_Operator
		// Only ApplyInverse is non-trivial -- first 9 satisfied here in header

		int SetUseTranspose(bool UseTranspose)
		{ cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<endl; return -1;};
		int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
		{ cerr<<"ERROR: ApplyInverse"<<endl; return -1;};
		int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
		double NormInf() const { cerr<<"norminf"<<endl; return 1.0;};
		const char* Label() const { return "noxlocainterface::user preconditioner";};
		bool UseTranspose() const { return false;};
		bool HasNormInf() const { return false;};
		const Epetra_Comm& Comm() const {return comm;};
		const Epetra_Map& OperatorDomainMap() const {cout<<"returning domain map"<<endl; return *globalMap;};
		const Epetra_Map& OperatorRangeMap() const {cout<<"returning range map"<<endl; return *globalMap;};

	private:
		int N;
		Teuchos::RCP<Epetra_Map> globalMap;
		const Epetra_Comm& comm;
		void* precdata;
		void (*precFunction)(double *, int, double*, void *);
		void (*precFunctionblock11)(double *, int, double*, void *);
		void (*precFunctionblock12)(double *, int, double*, void *);
		void (*precFunctionblock21)(double *, int, double*, void *);
		void (*precFunctionblock22)(double *, int, double*, void *);
		std::string label;
                bool printproc;

};

#endif

