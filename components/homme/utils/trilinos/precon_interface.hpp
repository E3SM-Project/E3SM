//-----------------------------------------------------------------------------
#ifndef precon_interface_H
#define precon_interface_H

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













class  Precon_Interface : 
	public Epetra_Operator
{
	public:
		Precon_Interface(int nelems,Teuchos::RCP<Epetra_Map> gmap,const Epetra_Comm& comm_,
				void* prec_data_,
				void (*precFunction_)(double *,int,double*,void *));


	//	Precon_Interface(const Epetra_Comm& comm_,int nelems,Teuchos::RCP<Epetra_Map> gmap);


		~Precon_Interface();

		// 10 Methods for inheritance from Epetra_Operator
		// Only ApplyInverse is non-trivial -- first 9 satisfied here in header

		int SetUseTranspose(bool UseTranspose)
		{ std::cerr<<"ERROR: No noxlocainterface::SetUseTranspose"<<std::endl; return -1;};
		int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)const
		{ std::cerr<<"ERROR: ApplyInverse"<<std::endl; return -1;};
		int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
		double NormInf() const { std::cerr<<"norminf"<<std::endl; return 1.0;};
		const char* Label() const { return "noxlocainterface::user preconditioner";};
		bool UseTranspose() const { return false;};
		bool HasNormInf() const { return false;};
		const Epetra_Comm& Comm() const {return comm;};
		const Epetra_Map& OperatorDomainMap() const {std::cout<<"returning domain map"<<std::endl; return *globalMap;};
		const Epetra_Map& OperatorRangeMap() const {std::cout<<"returning range map"<<std::endl; return *globalMap;};

	private:
		int N;
		Teuchos::RCP<Epetra_Map> globalMap;
		const Epetra_Comm& comm;
		void* precdata;
		void (*precFunction)(double *, int, double*, void *);
		std::string label;
                bool printproc;

};

#endif

