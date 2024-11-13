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

#define STANDALONE_OMEGA
    
    bool Pacer::initialize(MPI_Comm InComm) {

#ifdef STANDALONE_OMEGA
        // GPTL set default options
        GPTLsetoption(GPTLdepthlimit, 20);
        GPTLsetoption(GPTLdopr_quotes, 1);
        GPTLsetoption(GPTLprofile_ovhd, 1);
        // default is set to 52
        // GPTLsetoption(GPTLmaxthreads)
        
        GPTLsetutr(GPTLmpiwtime);

        if (MPI_Comm_dup(InComm, &InternalComm) != MPI_SUCCESS)
            std::cerr << "Pacer: Error duplicating MPI communicator" << std::endl;
        MPI_Comm_rank(InternalComm, &MyRank);

        GPTLinitialize();

        IsInitialized = true;
#endif
        return true;
    }


    bool Pacer::start(const std::string &TimerName)
    {
        GPTLstart(TimerName.c_str());
        auto it = OpenTimers.find(TimerName); 
        if (it != OpenTimers.end() )
           OpenTimers[TimerName]++;
        else
           OpenTimers[TimerName] = 1;
        return true;
    }

    bool Pacer::stop(const std::string &TimerName)
    {
        GPTLstop(TimerName.c_str());
        if ( OpenTimers[TimerName] == 1 )
            OpenTimers.erase(TimerName);
        else 
            OpenTimers[TimerName]--;

        return true;
    }

    bool Pacer::setPrefix(const std::string &Prefix)
    {
        GPTLprefix_set(Prefix.c_str());

        return true;
    }

    bool Pacer::unsetPrefix()
    {
        GPTLprefix_unset();
        return true;
    }

    bool Pacer::print(const std::string &TimerFilePrefix)
    {
        // https://github.com/E3SM-Project/E3SM/blob/master/share/timing/perf_mod.F90
        //GPTLpr(0);
        if (MyRank == 0) {
            std::string TimerFileName = TimerFilePrefix + ".timing";
            std::string SummaryFileName = TimerFilePrefix + ".summary";
            GPTLpr_file(TimerFileName.c_str());
            GPTLpr_summary_file(InternalComm, SummaryFileName.c_str());
        }
        // https://github.com/jmrosinski/GPTL/blob/master/tests/global.c
        return true;
    }

    bool Pacer::finalize()
    {
#ifdef STANDALONE_OMEGA
        GPTLfinalize();
#endif

        if ( (MyRank == 0) && ( OpenTimers.size() > 0) ){ 
            std::cerr << "PACER Warning: Following timers are not closed." << std::endl;
            for (auto i = OpenTimers.begin(); i != OpenTimers.end(); i++)
                std::cerr << i->first << std::endl;

        }

        return true;
    }


//===----------------------------------------------------------------------===//
