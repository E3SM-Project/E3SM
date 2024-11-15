#ifndef PACER_H
#define PACER_H

//===-- Pacer.h - Pacer timing interface -------------*- C++
//-*-===//
//
// \file
// \brief Provides timer functionality for E3SM
//
// The Pacer class provides an interface to timers for
// E3SM components. 
//
//===------------------------------------------------===//

#include <mpi.h>
#include <gptl.h>
#include <unordered_map>
#include <string>

namespace Pacer {
   /// Flag to determine if the timing infrastructure is initialized 
    static bool IsInitialized;

    static MPI_Comm InternalComm;

    static int MyRank;

    static std::unordered_map<std::string,int> OpenTimers;

    // Initialize Pacer timing. 
    // InComm: overall MPI communicator used by application.
    bool initialize(MPI_Comm InComm);

    bool start(const std::string &TimerName);

    bool stop(const std::string &TimerName);

    bool setPrefix(const std::string &Prefix);

    bool unsetPrefix();

    bool print(const std::string &TimerFilePrefix, bool PrintAllRanks = false);

    bool finalize();

};

#endif
