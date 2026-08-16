#ifndef OMEGA_MOABINTERFACE_H
#define OMEGA_MOABINTERFACE_H
//===-- drivers/coupled/MoabInterface.h - MOAB coupling bridge -*- C++ -*-===//
//
// \file
// \brief Bridge between Omega's ocean core and the MOAB coupling driver.
//
// MoabInterface.h/.cpp are the only Omega source files that include
// moab/iMOAB.h.
//
//===---------------------------------------------------------------------===//

#include "DataTypes.h"

#include <mpi.h>
#include <string>

namespace OMEGA {

// Registers this ocean instance as a MOAB application, builds its mesh from
// the default Decomp/HorzMesh (owned cells only), and defines the domain
// tags (GLOBAL_ID, lon, lat, area, mask, frac). Returns the MOAB
// application id.
int moabInit(MPI_Comm Comm, int OcnID);

// Defines dense, double, single-component tags for the coupler's full
// x2o/o2x field lists (colon-separated, null-terminated) on the ocean's
// MOAB application.
void moabDefineTagStorage(int Pid, const std::string &Cpl2OcnFieldNames,
                           const std::string &Ocn2CplFieldNames);

// Reads the x2o tag storage (NFields, NCellsOwned, unrolled by tag) into
// Buffer.
void moabImportTagStorage(Real *Buffer, int Len);

// Writes Buffer (NFields, NCellsOwned, unrolled by tag) into the o2x tag
// storage.
void moabExportTagStorage(const Real *Buffer, int Len);

} // namespace OMEGA
#endif // !OMEGA_MOABINTERFACE_H
