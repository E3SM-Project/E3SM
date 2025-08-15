//===-- Pacer.cpp - Timers for E3SM --*- C++ -*-===//
//
// \file
// \brief Timer infrastructure for E3SM C++ components
//
//
////===------------------------------------------===//

#include "Pacer.h"
#include <mpi.h>
#include <gptl.h>
#include <iostream>

#define PACER_CHECK_INIT() {\
    if (!IsInitialized) { \
        std::cerr << "[ERROR] Pacer: Not initialized." << std::endl; \
        return false; \
    } \
}

#define PACER_CHECK_ERROR(x) {\
    if ( (x) != 0 ) { \
        std::cerr << "[ERROR] Pacer: Failure calling GPTL function: " << #x << std::endl; \
        return false; \
    } \
}

/// Helper function to check if GPTL is initialized
/// Function declaration is missing in gptl.h
/// Hence the declaration here
extern "C" {
    extern int GPTLis_initialized(void);
}

/// Check if Pacer is initialized
/// Returns true if initialized
inline bool Pacer::isInitialized(void){
    if (!IsInitialized) {
        std::cerr << "[ERROR] Pacer: Not initialized." << std::endl;
        return false;
    }
    return true;
}

/// Initialize Pacer timing
/// InComm: overall MPI communicator used by application.
/// InMode: Pacer standalone (default) or within CIME
bool Pacer::initialize(MPI_Comm InComm, PacerModeType InMode /* = PACER_STANDALONE */) {

    int errCode;

    PacerMode = InMode;

    // Check if already initialized and return
    if (IsInitialized)
        return true;

    // Duplicate comm and get MPI rank for both standalone and integrated modes
    if (MPI_Comm_dup(InComm, &InternalComm) != MPI_SUCCESS)
        std::cerr << "Pacer: Error duplicating MPI communicator" << std::endl;
    MPI_Comm_rank(InternalComm, &MyRank);

    if (PacerMode == PACER_STANDALONE ) {
        // GPTL set default options
        PACER_CHECK_ERROR(GPTLsetoption(GPTLdepthlimit, 20));
        PACER_CHECK_ERROR(GPTLsetoption(GPTLdopr_quotes, 1));
        PACER_CHECK_ERROR(GPTLsetoption(GPTLprofile_ovhd, 1));
        // GPTL default is set to 52
        // Presently setting to 64
        PACER_CHECK_ERROR(GPTLsetoption(GPTLmaxthreads, 64));

        PACER_CHECK_ERROR(GPTLsetutr(GPTLmpiwtime));

        errCode = GPTLinitialize();

        if (errCode) {
            std::cerr << "[ERROR] Pacer: Unable to initialize GPTL." << std::endl;
            return false;
        }
        else
            IsInitialized = true;
    }
    else if (PacerMode == PACER_INTEGRATED) {
        // GPTL is assumed to be initialized by E3SM driver in coupled modeling context

        // Check if GPTL is initialized
        errCode = GPTLis_initialized();

        if (errCode == true) {
            IsInitialized = true;
        }
        else {
            IsInitialized = false;
            std::cerr << "[ERROR] Pacer: GPTL is not initialized in PACER_INTEGRATED mode." << std::endl;
            return false;
        }
    }

    return true;
}

/// Start the time named TimerName
bool Pacer::start(const std::string &TimerName)
{
    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLstart(TimerName.c_str()));

    auto it = OpenTimers.find(TimerName);
    if (it != OpenTimers.end() )
        OpenTimers[TimerName]++;
    else
        OpenTimers[TimerName] = 1;
    return true;
}

/// Stop the time named TimerName
/// Issues warning if timer hasn't been started yet
bool Pacer::stop(const std::string &TimerName)
{
    PACER_CHECK_INIT();

    auto it = OpenTimers.find(TimerName);

    if (it != OpenTimers.end() ) {
        PACER_CHECK_ERROR(GPTLstop(TimerName.c_str()));

        if ( OpenTimers[TimerName] == 1 )
            OpenTimers.erase(TimerName);
        else
            OpenTimers[TimerName]--;
    }
    else {
        std::cerr << "[WARNING] Pacer: Trying to stop timer: \""
            << TimerName << "\" before starting it." << std::endl;

        return false;
    }

    return true;
}

/// Sets named prefix for all subsequent timers
bool Pacer::setPrefix(const std::string &Prefix)
{
    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLprefix_set(Prefix.c_str()));

    return true;
}

/// Unsets prefix for all subsequent timers
bool Pacer::unsetPrefix()
{
    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLprefix_unset());

    return true;
}

/// Prints timing statistics and global summary files
/// Output Files: TimerFilePrefix.timing.<MyRank>
/// TimerFilePrefix.summary
/// PrintAllRanks: flag to control if per rank timing files are printed
bool Pacer::print(const std::string &TimerFilePrefix, bool PrintAllRanks /*= = false */)
{
    PACER_CHECK_INIT();

    std::string TimerFileName = TimerFilePrefix + ".timing." + std::to_string(MyRank);
    std::string SummaryFileName = TimerFilePrefix + ".summary";

    PACER_CHECK_ERROR(GPTLpr_summary_file(InternalComm, SummaryFileName.c_str()));

    if ( PrintAllRanks == false ) {
        if (MyRank == 0) {
            PACER_CHECK_ERROR(GPTLpr_file(TimerFileName.c_str()));
        }
    }
    else
        PACER_CHECK_ERROR(GPTLpr_file(TimerFileName.c_str()));

    return true;
}

/// Cleans up Pacer
/// Issues warning if any timers are still open
bool Pacer::finalize()
{
    PACER_CHECK_INIT();

    if ( PacerMode == PACER_STANDALONE )
        PACER_CHECK_ERROR(GPTLfinalize());

    if ( (MyRank == 0) && ( OpenTimers.size() > 0) ){
        std::cerr << "[WARNING] Pacer: Following " << OpenTimers.size() << " timer(s) is/are still open." << std::endl;
        for (auto i = OpenTimers.begin(); i != OpenTimers.end(); i++)
            std::cerr << '\t' << i->first << std::endl;
    }
    OpenTimers.clear();

    // Clear Pacer state and free communicator
    IsInitialized = false;
    MPI_Comm_free(&InternalComm);

    return true;
}


//===----------------------------------------------------------------------------===//
