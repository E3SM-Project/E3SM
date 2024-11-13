#ifndef PACER_H
#define PACER_H

//===-- Pacer.h - time stepper -----------------------*- C++
////-*-===//
////
///// \file
///// \brief Provides timer functionality for E3SM
/////
///// The Pacer class provides an interface to timers for
///// E3SM components. 
////
////===------------------------------------------------===//

#include <mpi.h>
#include <gptl.h>

#define STANDALONE_OMEGA

namespace Pacer
{
 //private:
   /// Flag to determine if the timing infrastructure is initialized 
    static bool IsInitialized;

    /// Timers will be output with this filename or the
    /// constructed filename based on this template
    //static std::string TimerFilePrefix;

    static MPI_Comm InternalComm;

    static std::map<std::string,int> openTimers;

    // public:
    bool initialize(MPI_Comm InComm);

    bool start(const std::string &TimerName);

    bool stop(const std::string &TimerName);

    bool setPrefix(const std::string &Prefix);

    bool unsetPrefix();

    bool print(const std::string &TimerFilePrefix);

    bool finalize();

};



#endif
