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
 * \file sector_cycle_breaker.cpp
 * \ingroup objects
 * \brief The SectorCycleBreaker source file.
 * \author Jim Naslund
 */

#include "util/base/include/definitions.h"
#include <list>
#include <stack>
#include <algorithm>
#include "containers/include/sector_cycle_breaker.h"
#include "containers/include/dependency_finder.h"
#include "util/logger/include/ilogger.h"
#include "marketplace/include/marketplace.h"

using namespace std;

SectorCycleBreaker::SectorCycleBreaker( Marketplace* aMarketPlace,
                                        const string& aRegionName):
ICycleBreaker(),
mRegionName( aRegionName),
mMarketplace( aMarketPlace )
{
}

//! Virtual destructor.
SectorCycleBreaker::~SectorCycleBreaker(){
}

void SectorCycleBreaker::breakCycle( DependencyFinder &aDependencyFinder,
                                     const size_t aFirstSector,
                                     const size_t aSecondSector ){
    // Notify that we are removing a cycle.
    ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );
    depFinderLog.setLevel( ILogger::DEBUG );
    depFinderLog << "Breaking cycle between " 
                 << aDependencyFinder.getNameFromIndex( aFirstSector ) << " and "
                 << aDependencyFinder.getNameFromIndex( aSecondSector ) << "." << endl;

    // Add simul markets to remove the dependency. Note that one of these
    // sectors may already have a simul market setup for it, the marketplace
    // will ignore the request to convert the market in that case.
    mMarketplace->resetToPriceMarket( aDependencyFinder.getNameFromIndex( aFirstSector ),
                                      mRegionName );
    mMarketplace->resetToPriceMarket( aDependencyFinder.getNameFromIndex( aSecondSector ),
                                      mRegionName );
    
    // Remove the cycle from the graph by removing both edges.
    aDependencyFinder.removeDependency( aFirstSector, aSecondSector );
    aDependencyFinder.removeDependency( aSecondSector, aFirstSector );
}


