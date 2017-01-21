#ifndef SPARSE_PREC_BASE_HPP
#define SPARSE_PREC_BASE_HPP

#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "trilinosModelEvaluator.hpp" // include preconditioner base class
#include "EpetraExt_MatrixMatrix.h"   // Matrix-Matrix operations
#include "Teuchos_TimeMonitor.hpp"    // Trilinos timers
#include "Teuchos_oblackholestream.hpp"

#include "Epetra_FECrsMatrix.h"

#include <string>

class SparsePreconditionerBase : public hommePreconditionerBase
{
public:
  SparsePreconditionerBase(const int np_,
                           const int nlev_,
                           const int nelemd_,
                           void* StateData_,
                           RCP<Epetra_Vector> xVec_,
                           RCP<Epetra_Map> xMap_);

  ~SparsePreconditionerBase() {

    delete[] VelWgt;
    delete[] PsWgt;
    delete[] TWgt;

    delete[] VelIdx;
    delete[] PsIdx;
    delete[] TIdx;

    delete[] VelIdx_unsorted;
    delete[] PsIdx_unsorted;
    delete[] TIdx_unsorted;

    delete[] nVelIdxMap;
    delete[] nPsIdxMap;
    delete[] nTIdxMap;

    for (int i=0; i < nVel_unique; i++) {
      delete[] VelIdxMap[i];
    }
    delete[] VelIdxMap;

    for (int i=0; i < nPs_unique; i++) {
      delete[] PsIdxMap[i];
    }
    delete[] PsIdxMap;

    for (int i=0; i < nT_unique; i++) {
      delete[] TIdxMap[i];
    }
    delete[] TIdxMap;

    delete[] row_nnz_varies;

    delete[] Urows;
    delete[] Ucols;
    delete[] Uvals;

    delete[] Vrows;
    delete[] Vcols;
    delete[] Vvals;
  }
  

  // methods that MUST be implemented by the preconditioner
  virtual int computePreconditioner(RCP<const Epetra_Vector> xVec, 
                                    void* StateData_) = 0;

  virtual int ApplyInverse(const Epetra_MultiVector& X, 
                           Epetra_MultiVector& Y) const = 0;

protected:
  int ZeroAblock(RCP<Epetra_FECrsMatrix> Ablock);
  int ZeroBblock(RCP<Epetra_FECrsMatrix> Bblock);
  int ZeroCblock(RCP<Epetra_FECrsMatrix> Cblock);
  int ZeroDblock(RCP<Epetra_FECrsMatrix> Dblock);
  int ZeroEblock(RCP<Epetra_FECrsMatrix> Eblock);
  int ZeroFblock(RCP<Epetra_FECrsMatrix> Fblock);
  int ZeroGblock(RCP<Epetra_FECrsMatrix> Gblock);
  int ZeroHblock(RCP<Epetra_FECrsMatrix> Hblock);

  int FillAblock(RCP<Epetra_FECrsMatrix> Ablock);
  int FillBblock(RCP<Epetra_FECrsMatrix> Bblock);
  int FillCblock(RCP<Epetra_FECrsMatrix> Cblock);
  int FillDblock(RCP<Epetra_FECrsMatrix> Dblock);
  int FillEblock(RCP<Epetra_FECrsMatrix> Eblock);
  int FillFblock(RCP<Epetra_FECrsMatrix> Fblock);
  int FillGblock(RCP<Epetra_FECrsMatrix> Gblock);
  int FillHblock(RCP<Epetra_FECrsMatrix> Hblock);

  void print_matrix(RCP<Epetra_FECrsMatrix> Matrix, 
                    const std::string &MatrixName);

protected:
  int myid;
  bool printproc;

  int np;     // number of nodes along an element edge
  int nlev;   // number of vertical levels
  int nelemd; // number of elements on this process

  void* StateData;

  // vector lengths WITH duplicate values
  int nState;
  int nVel;
  int nPs;
  int nT;

  // vector lengths WITHOUT duplicate values
  int nState_unique;
  int nVel_unique;
  int nPs_unique;
  int nT_unique;

  // sorted node IDs
  int* VelIdx;
  int* PsIdx;
  int* TIdx;

  // unsorted node IDs
  int* VelIdx_unsorted;
  int* PsIdx_unsorted;
  int* TIdx_unsorted;

  // number of times a node ID occurs
  int* nVelIdxMap;
  int* nPsIdxMap;
  int* nTIdxMap;

  // locations of sorted unique node IDs in the unsorted array with repeated nodes
  int** VelIdxMap;
  int** PsIdxMap;
  int** TIdxMap;

  // weights for scaling duplicate vector entries
  double* VelWgt;
  double* PsWgt;
  double* TWgt;

  // number of non-zeros
  int nnz_2D;
  int nnz_3D;
  int nnz_max;
  
private:  
  int nRows_2D;
  int nRows_3D;

  int row_nnz_fixed;
  int* row_nnz_varies;

  int* Urows;
  int* Ucols;
  double* Uvals;
  
  int* Vrows;
  int* Vcols;
  double* Vvals;

  // Timers
  RCP<Teuchos::Time> ZeroATime;
  RCP<Teuchos::Time> ZeroBTime;
  RCP<Teuchos::Time> ZeroCTime;
  RCP<Teuchos::Time> ZeroDTime;
  RCP<Teuchos::Time> ZeroETime;
  RCP<Teuchos::Time> ZeroFTime;
  RCP<Teuchos::Time> ZeroGTime;
  RCP<Teuchos::Time> ZeroHTime;

  RCP<Teuchos::Time> FillATime;
  RCP<Teuchos::Time> FillBTime;
  RCP<Teuchos::Time> FillCTime;
  RCP<Teuchos::Time> FillDTime;
  RCP<Teuchos::Time> FillETime;
  RCP<Teuchos::Time> FillFTime;
  RCP<Teuchos::Time> FillGTime;
  RCP<Teuchos::Time> FillHTime;
};

#endif
