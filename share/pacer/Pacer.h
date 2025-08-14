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
#include <string>

namespace Pacer {

    enum PacerModeType { PACER_STANDALONE, PACER_INTEGRATED };

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

}

#endif
