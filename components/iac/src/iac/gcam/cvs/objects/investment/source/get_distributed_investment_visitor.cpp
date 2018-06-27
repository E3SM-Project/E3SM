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
 * \file get_distributed_investment.cpp
 * \ingroup Objects
 * \brief The GetDistributedInvestmentVisitor class source file.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "investment/include/get_distributed_investment_visitor.h"
#include "technologies/include/base_technology.h"
#include "util/base/include/util.h"
#include "sectors/include/subsector.h"
#include "sectors/include/sector.h"
#include "marketplace/include/marketplace.h"
 
using namespace std;
extern Scenario* scenario;

/*!
 *\brief Constructor
 * \param aRegionName Name of the region since we starts the visiting below the
 *        region level.
 */
GetDistributedInvestmentVisitor::GetDistributedInvestmentVisitor( const string& aRegionName )
:mCurrentRegionName( aRegionName )
{
}

void GetDistributedInvestmentVisitor::startVisitSector( const Sector* aSector,
                                                 const int aPeriod )
{
    mCurrentSectorName = aSector->getName();
}

void GetDistributedInvestmentVisitor::endVisitSector( const Sector* aSector,
                                               const int aPeriod )
{
    mCurrentSectorName.clear();
}

void GetDistributedInvestmentVisitor::startVisitSubsector( const Subsector* aSubsector,
                                                    const int aPeriod )
{
    // check the hasCalibrationMarket flag which has preivously determined if
    // this investable needs calibrating
    mCurrentSubSectorName = aSubsector->getName();
    if ( aSubsector->hasCalibrationMarket() ){
        // get the amount of investment and set it into the calibration market
        // as the supply, note this includes fixed investment
        /*!
         * \warning If a subsector is all fixed it's calibration market should not be solvable.
         */
        double subsectorInvestment = aSubsector->getAnnualInvestment( aPeriod );
        Marketplace* marketplace = scenario->getMarketplace();
        string calibrationMrkName = "SubSec-"+mCurrentSubSectorName;
        marketplace->addToSupply( calibrationMrkName, mCurrentRegionName, subsectorInvestment, aPeriod, true );
    }
}

void GetDistributedInvestmentVisitor::endVisitSubsector( const Subsector* aSubsector,
                                                  const int aPeriod )
{
    mCurrentSubSectorName.clear();
}

void GetDistributedInvestmentVisitor::startVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                         const int aPeriod )
{
    // check the hasCalibrationMarket flag which has preivously determined if
    // this investable needs calibrating
    if ( aTechnology->hasCalibrationMarket() ){
        // get the amount of investment and set it into the calibration market
        // as the supply, note this includes fixed investment
        /*!
         * \warning If a technology is all fixed it's calibration market should not be solvable.
         */
        double technologyInvestment = aTechnology->getAnnualInvestment( aPeriod );
        Marketplace* marketplace = scenario->getMarketplace();
        string calibrationMrkName = "Tech-"+aTechnology->getName();
        marketplace->addToSupply( calibrationMrkName, mCurrentRegionName, technologyInvestment, aPeriod, true );
    }
}

void GetDistributedInvestmentVisitor::endVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                       const int aPeriod )
{
}
