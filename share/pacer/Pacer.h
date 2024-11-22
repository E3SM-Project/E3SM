//===-- Pacer.h - Pacer timing interface -------------*- C++
//-*-===//
//
// \file
// \brief Provides timer functionality for E3SM
//
// The Pacer class provides an interface to timers for
// E3SM components.
//
//===--------------------------------------------------===//
#ifndef PACER_H
#define PACER_H

#include <mpi.h>
#include <gptl.h>
#include <unordered_map>
#include <string>

namespace Pacer {

    /// Flag to determine if the timing infrastructure is initialized
    static bool IsInitialized = false;

    /// MPI communicator used within Pacer
    static MPI_Comm InternalComm;

    /// MPI rank of process
    static int MyRank;

    /// hash table of open timers with count (for multiple parents)
    static std::unordered_map<std::string,int> OpenTimers;

    enum PacerModeType { PACER_STANDALONE, PACER_INTEGRATED };

    /// Pacer Mode: standalone or within CIME
    static PacerModeType PacerMode;

    /// Initialize Pacer timing
    /// InComm: overall MPI communicator used by application.
    /// InMode: Pacer Mode: standalone (default) or within CIME
    bool initialize(MPI_Comm InComm, PacerModeType InMode = PACER_STANDALONE);

    /// Check if Pacer is initialized
    /// Returns true if initialized
    bool isInitialized(void);

    /// IntegratedMode: Pacer standalone (default:false) or within CIME (true)
    // bool initialize(MPI_Comm InComm, bool IntegratedMode = false);

    /// Start the time named TimerName
    bool start(const std::string &TimerName);

    /// Stop the time named TimerName
    /// Issues warning if timer hasn't been started yet
    bool stop(const std::string &TimerName);

    /// Sets named prefix for all subsequent timers
    bool setPrefix(const std::string &Prefix);

    /// Unsets prefix for all subsequent timers
    bool unsetPrefix();

    /// Prints timing statistics and global summary files
    /// Output Files: TimerFilePrefix.timing.<MyRank>
    /// TimerFilePrefix.summary
    /// PrintAllRanks: flag to control if per rank timing files are printed
    bool print(const std::string &TimerFilePrefix, bool PrintAllRanks = false);

    /// Cleans up Pacer
    /// Issues warning if any timers are still open
    bool finalize();

};

#endif
