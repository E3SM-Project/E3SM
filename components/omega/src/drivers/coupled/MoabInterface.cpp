//===-- drivers/coupled/MoabInterface.cpp - MOAB coupling bridge
//-----------===//
//
//===---------------------------------------------------------------------===//

#include "MoabInterface.h"

#include "Decomp.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"

#include "moab/iMOAB.h"

#include <vector>

namespace OMEGA {

namespace {

int Pid = -1;

// Registers this ocean instance as a MOAB application and returns its pid
int registerApplication(MPI_Comm Comm, int OcnID) {

   ErrCode Err;
   int LocalPid;
   int CompID = OcnID;

   std::string AppName = "OMEGA_MB_" + std::format("{:02}", OcnID);

   Err = iMOAB_RegisterApplication(AppName.c_str(), &Comm, &CompID, &LocalPid);

   if (Err != moab::MB_SUCCESS)
      ABORT_ERROR("iMOAB_RegisterApplication: FAIL");

   return LocalPid;
}
} // namespace

int moabInit(MPI_Comm Comm, int OcnID) {
   Pid = registerApplication(Comm, OcnID);

   return Pid;
}
} // namespace OMEGA
