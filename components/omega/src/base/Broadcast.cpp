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

namespace OMEGA {

//------------------------------------------------------------------------------
// Broadcast I4 scalar Value
int Broadcast(I4 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast(&Value, 1, MPI_INT32_T, Root, InEnv->getComm());

   return RetVal;
} // end Broadcast for I4 data type

int Broadcast(I4 &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast for I4 data type

//------------------------------------------------------------------------------
// Broadcast I8 scalar Value
int Broadcast(I8 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast(&Value, 1, MPI_INT64_T, Root, InEnv->getComm());

   return RetVal;
} // end Broadcast for I8 data type

int Broadcast(I8 &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast for I8 data type

//------------------------------------------------------------------------------
// Broadcast R4 scalar Value
int Broadcast(R4 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast(&Value, 1, MPI_FLOAT, Root, InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(R4 &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 scalar Value
int Broadcast(R8 &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast(&Value, 1, MPI_DOUBLE, Root, InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(R8 &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast bool scalar Value
int Broadcast(bool &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast(&Value, 1, MPI_C_BOOL, Root, InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(bool &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast std::string Value
int Broadcast(std::string &Value, const MachEnv *InEnv, const int RankBcast) {
   int RetVal, Root, StrSize;
   MPI_Comm Comm = InEnv->getComm();
   int MyTask    = InEnv->getMyTask();

   Root = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   // For strings, we must ensure the destination has the correct size
   if (MyTask == Root)
      StrSize = Value.length();
   RetVal = MPI_Bcast(&StrSize, 1, MPI_INT, Root, Comm);
   if (MyTask != Root)
      Value.resize(StrSize);

   // Now broadcast the string
   RetVal =
       MPI_Bcast((void *)Value.c_str(), Value.size(), MPI_CHAR, Root, Comm);

   return RetVal;
} // end Broadcast

int Broadcast(std::string &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I4 array
int Broadcast(std::vector<I4> &Value, const MachEnv *InEnv,
              const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_INT32_T, Root,
                      InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(std::vector<I4> &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast I8 array
int Broadcast(std::vector<I8> &Value, const MachEnv *InEnv,
              const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_INT64_T, Root,
                      InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(std::vector<I8> &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R4 array
int Broadcast(std::vector<R4> &Value, const MachEnv *InEnv,
              const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_FLOAT, Root,
                      InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(std::vector<R4> &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
} // end Broadcast

//------------------------------------------------------------------------------
// Broadcast R8 array
int Broadcast(std::vector<R8> &Value, const MachEnv *InEnv,
              const int RankBcast) {
   int RetVal, Root;

   Root   = (RankBcast < 0) ? InEnv->getMasterTask() : RankBcast;
   RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_DOUBLE, Root,
                      InEnv->getComm());

   return RetVal;
} // end Broadcast

int Broadcast(std::vector<R8> &Value, const int RankBcast) {
   return Broadcast(Value, MachEnv::getDefault(), RankBcast);
}

//------------------------------------------------------------------------------
// Broadcast bool array
// Elements of vector<bool> seem to be non-addressable
// int Broadcast(std::vector<bool> &Value, const MachEnv *InEnv, const int
// RankBcast ) {
//    int RetVal;
//
//    RetVal = MPI_Bcast((void *)Value.data(), Value.size(), MPI_C_BOOL,
//    RankBcast, InEnv->getComm());
//
//    return RetVal;
//} // end Broadcast
//
// int Broadcast(std::vector<bool> &Value, const int RankBcast
//) {
//    return Broadcast(Value, MachEnv::getDefault(), RankBcast);
//}

} // namespace OMEGA
