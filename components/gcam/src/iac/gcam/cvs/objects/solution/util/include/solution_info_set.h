#ifndef _SOLUTION_INFO_SET_H_
#define _SOLUTION_INFO_SET_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */


/*! 
* \file solution_info_set.h
* \ingroup Solution
* \brief A class which contains a set of SolutionInfo objects.
* \author Josh Lurz, Sonny Kim
*/
#include "solution/util/include/solution_info.h" // Maybe use pointer instead.
class Marketplace;
class ILogger;
class World;
class ISolutionInfoFilter;
class SolutionInfoParamParser;

class SolutionInfoSet {
    friend std::ostream& operator<<( std::ostream& os, const SolutionInfoSet& aSolutionInfoSet ){
        aSolutionInfoSet.print( os );
        return os;
    }
public:
    enum UpdateCode {
        UNCHANGED,
        ADDED,
        REMOVED,
        ADDED_AND_REMOVED
    };

    typedef std::vector<SolutionInfo>::iterator SetIterator;
    typedef std::vector<SolutionInfo>::const_iterator ConstSetIterator;
    SolutionInfoSet( Marketplace* marketplace );
    SolutionInfoSet( const std::vector<SolutionInfo> aSolutionSet );
    void init( const unsigned int aPeriod, const double aDefaultSolutionTolerance, const double aDefaultSolutionFloor,
               const SolutionInfoParamParser* aSolutionInfoParamParser );
    void merge( const std::vector<SolutionInfo> aSolutionSet );
    void updateFromMarkets();
    void updateToMarkets();
    UpdateCode updateSolvable( const ISolutionInfoFilter* aSolutionInfoFilter );
    void updateElasticities();
    void adjustBrackets();
    void storeValues();
    void restoreValues();
    void resetBrackets();
    bool checkAndResetBrackets();
    SolutionInfo* getWorstSolutionInfo( const bool aIgnoreBisected = false );
    SolutionInfo* getWorstSolutionInfoReverse( const bool aIgnoreBisected = false );
    SolutionInfo* getPolicyOrWorstSolutionInfo();
    SolutionInfo* getPolicySolutionInfo();
    double getMaxRelativeExcessDemand() const;
    double getMaxAbsoluteExcessDemand() const;
    bool isAllBracketed() const;
    const std::vector<double> getDemands() const; // move derivatives and make me private!
    const std::vector<double> getSupplies() const;
    unsigned int getNumSolvable() const;
    unsigned int getNumTotal() const;
    const SolutionInfo& getSolvable( unsigned int index ) const;
    SolutionInfo& getUnsolved( unsigned int index );
    SolutionInfo& getSolvable( unsigned int index );
    const SolutionInfo& getAny( unsigned int index ) const;
    SolutionInfo& getAny( unsigned int index );
    std::vector<SolutionInfo> getSolvableSet() const;
    std::vector<SolutionInfo> getSolvedSet() const;
    std::vector<SolutionInfo> getUnsolvedSet() const;
    bool isAllSolved();
    bool hasSingularUnsolved();
    void unsetBisectedFlag();
    void printUnsolved( std::ostream& out );
    void findAndPrintSD( World* aWorld, Marketplace* aMarketplace, const int aPeriod, ILogger& aLogger );
    void printMarketInfo( const std::string& comment, const double worldCalcCount, std::ostream& out ) const;
    void printDerivatives( std::ostream& aOut ) const;
private:
    unsigned int period;
    Marketplace* marketplace;
    std::vector<SolutionInfo> solvable;
    std::vector<SolutionInfo> unsolved; // solvable markets that are not currently solved
    std::vector<SolutionInfo> unsolvable;
    void print( std::ostream& out ) const;
};

#endif // _SOLUTION_INFO_SET_H_
