#ifndef OMEGA_BROADCAST_H
#define OMEGA_BROADCAST_H
//===-- infra/Broadcast.h - Broadcasting Values --*- C++ -*-===//
//
/// \file
/// \brief Defines MPI broadcasting functions
///
/// This header defines functions to broadcast Values from one MPI rank
/// to another.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "mpi.h"

namespace OMEGA {

// blocking broadcast scalar
int Broadcast(I4 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(I4 &Value, const int RankBcast);

int Broadcast(I8 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(I8 &Value, const int RankBcast);

int Broadcast(R4 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(R4 &Value, const int RankBcast);

int Broadcast(R8 &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(R8 &Value, const int RankBcast);

int Broadcast(bool &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(bool &Value, const int RankBcast);

int Broadcast(std::string &Value, const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast = -1);
int Broadcast(std::string &Value, const int RankBcast);

// blocking broadcast array
int Broadcast(std::vector<I4> &Value,
              const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast  = -1);
int Broadcast(std::vector<I4> &Value, const int RankBcast);

int Broadcast(std::vector<I8> &Value,
              const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast  = -1);
int Broadcast(std::vector<I8> &Value, const int RankBcast);

int Broadcast(std::vector<R4> &Value,
              const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast  = -1);
int Broadcast(std::vector<R4> &Value, const int RankBcast);

int Broadcast(std::vector<R8> &Value,
              const MachEnv *InEnv = MachEnv::getDefault(),
              const int RankBcast  = -1);
int Broadcast(std::vector<R8> &Value, const int RankBcast);

// NOTE: Elements of vector<bool> seem to be non-addressable
// int Broadcast(std::vector<bool> &Value,
//              const MachEnv *InEnv = MachEnv::getDefault(),
//              const int RankBcast = -1);
// int Broadcast(std::vector<bool> &Value, const int RankBcast);

// non-blocking broadcast scalar:

// void IBroadcast(const I4 Value, const MachEnv *InEnv, const int RankIBcast =
// -1);

// void IBroadcast(const I8 Value, const MachEnv *InEnv, const int RankIBcast =
// -1);

// void IBroadcast(const R4 Value, const MachEnv *InEnv, const int RankIBcast =
// -1);

// void IBroadcast(const R8 Value, const MachEnv *InEnv, const int RankIBcast =
// -1);

// void IBroadcast(const Real Value, const MachEnv *InEnv, const int RankIBcast
// = -1);

// void IBroadcast(const bool Value, const MachEnv *InEnv, const int RankIBcast
// = -1);

// void IBroadcast(const std::string Value, const MachEnv *InEnv, const int
// RankIBcast = -1);

// non-blocking broadcast array

// void IBroadcast(const std::vector<I4> Value, const MachEnv *InEnv, const int
// RankBcast = -1);

// void IBroadcast(const std::vector<I8> Value, const MachEnv *InEnv, const int
// RankBcast = -1);

// void IBroadcast(const std::vector<R4> Value, const MachEnv *InEnv, const int
// RankBcast = -1);

// void IBroadcast(const std::vector<R8> Value, const MachEnv *InEnv, const int
// RankBcast = -1);

// void IBroadcast(const std::vector<Real> Value, const MachEnv *InEnv, const
// int RankBcast = -1);

// void IBroadcast(const std::vector<bool> Value, const MachEnv *InEnv, const
// int RankBcast = -1);

} // namespace OMEGA

#endif // OMEGA_BROADCAST_H
