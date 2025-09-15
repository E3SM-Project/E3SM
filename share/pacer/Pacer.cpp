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
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef PACER_HAVE_KOKKOS
#include <Kokkos_Core.hpp>

#if defined(PACER_ADD_RANGES) && defined(KOKKOS_ENABLE_CUDA)
#include <nvtx3/nvToolsExt.h>
#endif

#endif

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

namespace Pacer {

/// Flag to determine if the timing infrastructure is initialized
static bool IsInitialized = false;

/// MPI communicator used within Pacer
static MPI_Comm InternalComm;

/// MPI rank of process
static int MyRank;

/// Vector-based stack of open timers
static std::vector<std::string> OpenTimers;

/// GPTL doesn't seem to provide a function to obtain the current prefix
/// so we track it ourselves
static std::string CurrentPrefix = "";

/// Pacer Mode: standalone or within CIME
static PacerModeType PacerMode;

/// Global timing level
static int TimingLevel = 0;

/// Flag to determine if timing barriers are enabled
static bool TimingBarriersEnabled = false;

/// Flag to determine if automatic Kokkos fences are enabled
static bool AutoFenceEnabled = false;

/// Check if Pacer is initialized
/// Returns true if initialized
inline bool isInitialized(){
    if (!IsInitialized) {
        std::cerr << "[ERROR] Pacer: Not initialized." << std::endl;
        return false;
    }
    return true;
}

/// Initialize Pacer timing
/// InComm: overall MPI communicator used by application.
/// InMode: Pacer standalone (default) or within CIME
bool initialize(MPI_Comm InComm, PacerModeType InMode /* = PACER_STANDALONE */) {

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

/// Start the timer named TimerName active when TimingLevel >= Level 
bool start(const std::string &TimerName, int Level)
{
    // Return immediately if this timer level is above the global timing level
    if (Level > TimingLevel) {
        return true;
    }

#ifdef PACER_HAVE_KOKKOS
    if (AutoFenceEnabled) {
        Kokkos::fence();
    }

    // If CUDA is enabled start an NVTX range
#if defined(PACER_ADD_RANGES) && defined(KOKKOS_ENABLE_CUDA)
    nvtxRangePush(TimerName.c_str());
#endif

#endif

    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLstart(TimerName.c_str()));

    // Push this timer onto the stack
    OpenTimers.push_back(TimerName);

    return true;
}

/// Stop the timer named TimerName active when TimingLevel >= Level 
/// Issues warning if timer hasn't been started yet
bool stop(const std::string &TimerName, int Level)
{
    // Return immediately if this timer level is above the global timing level
    if (Level > TimingLevel) {
        return true;
    }

    PACER_CHECK_INIT();

    auto it = std::find(OpenTimers.begin(), OpenTimers.end(), TimerName);

    if (it != OpenTimers.end() ) {

#ifdef PACER_HAVE_KOKKOS

        // If CUDA is enabled end the NVTX range
#if defined(PACER_ADD_RANGES) && defined(KOKKOS_ENABLE_CUDA)
        nvtxRangePop();
#endif

        if (AutoFenceEnabled) {
            Kokkos::fence();
        }
#endif

        PACER_CHECK_ERROR(GPTLstop(TimerName.c_str()));

        // Pop this timer from the stack
        OpenTimers.pop_back();
    }
    else {
        std::cerr << "[WARNING] Pacer: Trying to stop timer: \""
            << TimerName << "\" before starting it." << std::endl;

        return false;
    }

    return true;
}

/// Sets named prefix for all subsequent timers
bool setPrefix(const std::string &Prefix)
{
    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLprefix_set(Prefix.c_str()));

    CurrentPrefix = Prefix;

    return true;
}

/// Unsets prefix for all subsequent timers
bool unsetPrefix()
{
    PACER_CHECK_INIT();

    PACER_CHECK_ERROR(GPTLprefix_unset());

    CurrentPrefix = "";

    return true;
}

/// Adds the current enclosing timer name to the prefix for all subsequent timers
bool addParentPrefix()
{
    PACER_CHECK_INIT();

    const std::string NewPrefix = CurrentPrefix + OpenTimers.back() + ":";
    
    PACER_CHECK_ERROR(GPTLprefix_set(NewPrefix.c_str()));
    
    CurrentPrefix = NewPrefix;

    return true;
}

/// Removed the current enclosing timer name to the prefix for all subsequent timers
bool removeParentPrefix()
{
    PACER_CHECK_INIT();
    
    std::string NewPrefix = CurrentPrefix;
    const int NCharsToErase = OpenTimers.back().length() + 1;
    NewPrefix.erase(NewPrefix.length() - NCharsToErase);
    
    PACER_CHECK_ERROR(GPTLprefix_set(NewPrefix.c_str()));
    
    CurrentPrefix = NewPrefix;
    
    return true;
}

void setTimingLevel(int Level)
{
    TimingLevel = Level;
}

bool disableTiming()
{
    PACER_CHECK_ERROR(GPTLdisable());
    
    return true;
}

bool enableTiming()
{
    PACER_CHECK_ERROR(GPTLenable());
    
    return true;
}

void enableTimingBarriers()
{
    TimingBarriersEnabled = true;
}

void disableTimingBarriers()
{
    TimingBarriersEnabled = false;
}

void enableAutoFence()
{
    AutoFenceEnabled = true;
}

void disableAutoFence()
{
    AutoFenceEnabled = false;
}

bool timingBarrier(const std::string &TimerName, int Level, MPI_Comm Comm)
{
    bool Ok = true;
    if (TimingBarriersEnabled && Level <= TimingLevel) {

        Ok = Ok && start(TimerName, Level);
        MPI_Barrier(Comm);
        Ok = Ok && stop(TimerName, Level);
    }
    return Ok;
}

/// Prints timing statistics and global summary files
/// Output Files: TimerFilePrefix.timing.<MyRank>
/// TimerFilePrefix.summary
/// PrintAllRanks: flag to control if per rank timing files are printed
bool print(const std::string &TimerFilePrefix, bool PrintAllRanks /*= = false */)
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
bool finalize()
{
    PACER_CHECK_INIT();

    if ( PacerMode == PACER_STANDALONE )
        PACER_CHECK_ERROR(GPTLfinalize());

    if ( (MyRank == 0) && ( OpenTimers.size() > 0) ){
        std::cerr << "[WARNING] Pacer: Following " << OpenTimers.size() << " timer(s) is/are still open." << std::endl;
        for (const auto& Timer : OpenTimers)
            std::cerr << '\t' << Timer << std::endl;
    }
    OpenTimers.clear();

    // Clear Pacer state and free communicator
    IsInitialized = false;
    MPI_Comm_free(&InternalComm);

    return true;
}

}


//===----------------------------------------------------------------------------===//
