//===-- infra/Broacast.cpp - Implements Broadcasting functions --*- C++ -*-===//
//
/// \file
/// \brief implements blocking and non-blocking broadcasting functions
///
/// This implements blocking and non-blocking broadcasting functions. Two
/// function names are overloaded accross Omega data types: I4, I8, R4, R8,
/// Real, bool, and std::string: 1) Broadcast for blocking mode and 2) mpiIbcast
/// for non- blocking mode.
//
//===----------------------------------------------------------------------===//

#include "Broadcast.h"
#include "Error.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Broadcast I4 scalar Value
void Broadcast(I4 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast(&Value, 1, MPI_INT32_T, Root, InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast I4: Error in MPI Broadcast");
   }

   return;
} // end Broadcast for I4 data type

void Broadcast(I4 &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast for I4 data type

//------------------------------------------------------------------------------
// Broadcast I8 scalar Value
void Broadcast(I8 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast(&Value, 1, MPI_INT64_T, Root, InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast I8: Error in MPI Broadcast");
   }

   return;
} // end Broadcast for I8 data type

void Broadcast(I8 &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast for I8 data type

//------------------------------------------------------------------------------
// Broadcast R4 scalar Value
void Broadcast(R4 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast(&Value, 1, MPI_FLOAT, Root, InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast R4: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(R4 &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 scalar Value
void Broadcast(R8 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast(&Value, 1, MPI_DOUBLE, Root, InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast R8: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(R8 &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast bool scalar Value
void Broadcast(bool &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast(&Value, 1, MPI_C_BOOL, Root, InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast bool: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(bool &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast std::string Value
void Broadcast(std::string &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root, StrSize;
   MPI_Comm Comm = InEnv->getComm();
   int MyTask    = InEnv->getMyTask();

   if (InEnv->isMember()) {
      Root = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      // For strings, we must ensure the destination has the correct size
      if (MyTask == Root)
         StrSize = Value.length();
      RetVal = MPI_Bcast(&StrSize, 1, MPI_INT, Root, Comm);
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast string: Error in Broadcast of string size");
      if (MyTask != Root)
         Value.resize(StrSize);

      // Now broadcast the string
      RetVal = MPI_Bcast(&Value[0], Value.size(), MPI_CHAR, Root, Comm);
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast string: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(std::string &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I4 array
void Broadcast(std::vector<I4> &Value, const MachEnv *InEnv,
               const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_INT32_T, Root,
                         InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast I4 vector: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(std::vector<I4> &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I8 array
void Broadcast(std::vector<I8> &Value, const MachEnv *InEnv,
               const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_INT64_T, Root,
                         InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast I8 vector: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(std::vector<I8> &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R4 array
void Broadcast(std::vector<R4> &Value, const MachEnv *InEnv,
               const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_FLOAT, Root,
                         InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast R4 vector: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(std::vector<R4> &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 array
void Broadcast(std::vector<R8> &Value, const MachEnv *InEnv,
               const int RankBcast) {
   int RetVal, Root;

   if (InEnv->isMember()) {
      Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
      RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_DOUBLE, Root,
                         InEnv->getComm());
      if (RetVal != MPI_SUCCESS)
         ABORT_ERROR("Broadcast R8 vector: Error in MPI Broadcast");
   }

   return;
} // end Broadcast

void Broadcast(std::vector<R8> &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
}

//------------------------------------------------------------------------------
// Broadcast bool array
void Broadcast(std::vector<bool> &Value, const MachEnv *InEnv,
               const int RankBcast) {

   // Due to special packing of std::vector<bool> we convert to an integer
   I4 VecSize = Value.size();
   std::vector<I4> TmpValue(VecSize);
   for (int I = 0; I < VecSize; ++I) {
      if (Value[I])
         TmpValue[I] = 1;
      else
         TmpValue[I] = 0;
   }
   Broadcast(TmpValue, InEnv, RankBcast);

   // Convert back to bool
   for (int I = 0; I < VecSize; ++I) {
      if (TmpValue[I] == 1)
         Value[I] = true;
      else
         Value[I] = false;
   }

   return;
}

void Broadcast(std::vector<bool> &Value, const int RankBcast) {
   Broadcast(Value, MachEnv::getDefault(), RankBcast);
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//
