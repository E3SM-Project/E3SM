#ifndef OMEGA_BROADCAST_H
#define OMEGA_BROADCAST_H
//===-- infra/Broadcast.h - Broadcasting Values -----------------*- C++ -*-===//
//
/// \file
/// \brief Defines MPI broadcasting functions
///
/// This header defines functions to broadcast scalars and vectors from one MPI
/// rank to all others.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "mpi.h"

namespace OMEGA {

// blocking broadcast scalar
void Broadcast(I4 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(I4 &Value, const int RankBcast);

void Broadcast(I8 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(I8 &Value, const int RankBcast);

void Broadcast(R4 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(R4 &Value, const int RankBcast);

void Broadcast(R8 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(R8 &Value, const int RankBcast);

void Broadcast(bool &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(bool &Value, const int RankBcast);

void Broadcast(std::string &Value, const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast = -1);
void Broadcast(std::string &Value, const int RankBcast);

// blocking broadcast vector
void Broadcast(std::vector<I4> &Value,
               const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast  = -1);
void Broadcast(std::vector<I4> &Value, const int RankBcast);

void Broadcast(std::vector<I8> &Value,
               const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast  = -1);
void Broadcast(std::vector<I8> &Value, const int RankBcast);

void Broadcast(std::vector<R4> &Value,
               const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast  = -1);
void Broadcast(std::vector<R4> &Value, const int RankBcast);

void Broadcast(std::vector<R8> &Value,
               const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast  = -1);
void Broadcast(std::vector<R8> &Value, const int RankBcast);

void Broadcast(std::vector<bool> &Value,
               const MachEnv *InEnv = MachEnv::getDefault(),
               const int RankBcast  = -1);
void Broadcast(std::vector<bool> &Value, const int RankBcast);

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_BROADCAST_H
