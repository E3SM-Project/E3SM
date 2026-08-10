#ifndef OMEGA_MOABINTERFACE_H
#define OMEGA_MOABINTERFACE_H
//===-- drivers/coupled/MoabInterface.h - MOAB coupling bridge -*- C++ -*-===//
//
// \file
// \brief
//
//
//===---------------------------------------------------------------------===//

#include "DataTypes.h"

#include <mpi.h>
#include <string>

namespace OMEGA {

int moabInit(MPI_Comm Comm, int OcnID);

void moabDefineTagStoragee(int Pid, const std::string &Cpl2OcnFieldNames,
                           const std::string &Ocn2CplFieldNames)

    // Read x2o tag storage (NFields, NCellsOwned) into buffer
    void moabImportTagStorgae(Real *Buffer, int Len);

// Write Buffer (NFields, NCellsOwned) into the o2x tag storage
void moabExportTagStorage(const Real *Buffer, int Len)

} // namespace OMEGA
#endif // !OMEGA_MOABINTERFACE_H
