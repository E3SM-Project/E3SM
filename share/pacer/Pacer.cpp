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

#define PACER_STANDALONE_MODE

    bool Pacer::initialize(MPI_Comm InComm) {

#ifdef PACER_STANDALONE_MODE
        // GPTL set default options
        GPTLsetoption(GPTLdepthlimit, 20);
        GPTLsetoption(GPTLdopr_quotes, 1);
        GPTLsetoption(GPTLprofile_ovhd, 1);
        // default is set to 52
        // GPTLsetoption(GPTLmaxthreads)

        GPTLsetutr(GPTLmpiwtime);

        GPTLinitialize();

        IsInitialized = true;
#endif

        if (MPI_Comm_dup(InComm, &InternalComm) != MPI_SUCCESS)
            std::cerr << "Pacer: Error duplicating MPI communicator" << std::endl;
        MPI_Comm_rank(InternalComm, &MyRank);

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
        auto it = OpenTimers.find(TimerName);

        if (it != OpenTimers.end() ) {
            GPTLstop(TimerName.c_str());

            if ( OpenTimers[TimerName] == 1 )
                OpenTimers.erase(TimerName);
            else
                OpenTimers[TimerName]--;
        }
        else {
            std::cerr << "Pacer Warning: Trying to stop timer:"
                << TimerName << "before starting it." << std::endl;

            return false;
        }

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

    bool Pacer::print(const std::string &TimerFilePrefix, bool PrintAllRanks /*= = false */)
    {
        std::string TimerFileName = TimerFilePrefix + ".timing." + std::to_string(MyRank);
        std::string SummaryFileName = TimerFilePrefix + ".summary";
        GPTLpr_summary_file(InternalComm, SummaryFileName.c_str());

        if ( PrintAllRanks == false ) {
            if (MyRank == 0) {
                GPTLpr_file(TimerFileName.c_str());
            }
        }
        else
            GPTLpr_file(TimerFileName.c_str());

        return true;
    }

    bool Pacer::finalize()
    {
#ifdef PACER_STANDALONE_MODE
        GPTLfinalize();
#endif
        if ( (MyRank == 0) && ( OpenTimers.size() > 0) ){
            std::cerr << "PACER Warning: Following timers are not closed." << std::endl;
            for (auto i = OpenTimers.begin(); i != OpenTimers.end(); i++)
                std::cerr << '\t' << i->first << std::endl;

        }

        return true;
    }


//===----------------------------------------------------------------------===//
