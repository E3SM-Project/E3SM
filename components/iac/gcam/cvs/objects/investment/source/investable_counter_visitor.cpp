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
 * \file investable_counter_visitor.cpp
 * \ingroup Objects
 * \brief The InvestableCounterVisitor class source file.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "investment/include/investable_counter_visitor.h"
#include "technologies/include/base_technology.h"
#include "util/base/include/util.h"
#include "sectors/include/subsector.h"
#include "sectors/include/sector.h"
#include "marketplace/include/marketplace.h"
 
using namespace std;
extern Scenario* scenario;

/*!
 *\brief Constructor
 * \param aRegionName Name of the region since it will start the visiting below the
 *        region level.
 * \param aInvestableCount The count of investables below the sector level
 *          which will be used to determine if a calibration market is necessary
 *          at the subsector level.
 */
InvestableCounterVisitor::InvestableCounterVisitor( const string& aRegionName,
                                                    const int aInvestableCount )
:mCurrentRegionName( aRegionName ),
mMultipleTechnologies( false ),
mSubSectorCount( 0 ),
mTechCount( 0 ),
mMultipleSubsectors( aInvestableCount > 1 )
{
}

void InvestableCounterVisitor::startVisitSubsector( const Subsector* aSubsector,
                                                    const int aPeriod )
{
    mCurrentSubSectorName = aSubsector->getName();
    // Check if there are more than one technologies
    int countTechs = 0;
    for( unsigned int i = 0; i < aSubsector->baseTechs.size(); ++i ){
        if( aSubsector->baseTechs[ i ]->isNewInvestment( aPeriod ) ){
            ++countTechs;
        };
    }
    // if there are more than one then set the flag so that when we are visiting
    // technologies they can know whether they should create calibration markets
    mMultipleTechnologies = countTechs > 1;

    // only create calibration markets if there are multiple subsectors
    if ( mMultipleSubsectors ){
        Marketplace* marketplace = scenario->getMarketplace();
        double subsectorInvestment = aSubsector->getAnnualInvestment( aPeriod );
        double shareWeight = aSubsector->getShareWeight( aPeriod );
        string calibrationMrkName = "SubSec-"+mCurrentSubSectorName;

        // create the calibration market of type subsidy since we will only have access to the required investment
        // now and so need to make sure it does not get cleared from the market
        marketplace->createMarket( mCurrentRegionName, mCurrentRegionName, calibrationMrkName, IMarketType::SUBSIDY );
        marketplace->addToDemand( calibrationMrkName, mCurrentRegionName, subsectorInvestment, aPeriod, true );
        // aSubsector is const, but we need to set this flag to simplify the other calibration visitor routines
        const_cast<Subsector*>(aSubsector)->doCalibration = true;
        // For the first subsector
        if( mSubSectorCount == 0 || 
            ( aSubsector->getFixedInvestment( aPeriod) - subsectorInvestment ) == 0 ){
            // Share weights are relative so set first to 1,
            // and do not set market to solve.
            // we also need to make sure the fixed subsectors do not get set to solve
            shareWeight = 1.0;
            marketplace->setPrice( calibrationMrkName, mCurrentRegionName, shareWeight, aPeriod, true );
        }
        else{
            // we will need to solve this one and can start with the initial share weight set
            marketplace->setMarketToSolve( calibrationMrkName, mCurrentRegionName, aPeriod );
            marketplace->setPrice( calibrationMrkName, mCurrentRegionName, shareWeight, aPeriod, true );
        }

        // increase the subsector count so that we know we have already have created a subsector
        // that will be help at a constant share weight
        ++mSubSectorCount;
    }
}

void InvestableCounterVisitor::endVisitSubsector( const Subsector* aSubsector,
                                                  const int aPeriod )
{
    // make sure to clear our temporary flags and counters so we have
    // a fresh count for the next subsector
    mCurrentSubSectorName.clear();
    mMultipleTechnologies = false;
    mTechCount = 0;
}

void InvestableCounterVisitor::startVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                         const int aPeriod )
{
    // only create the market if there are multiple technologies and only when we are visiting
    // the new investment for that technology
    if ( mMultipleTechnologies && aTechnology->isNewInvestment( aPeriod ) ){
        //aTechnology is const, but we need to set this flag to simplify the other calibration visitor routines
        const_cast<BaseTechnology*>(aTechnology)->doCalibration = true;

        double technologyInvestment = aTechnology->getAnnualInvestment( aPeriod );
        double shareWeight = aTechnology->getShareWeight( aPeriod );
        Marketplace* marketplace = scenario->getMarketplace();
        string calibrationMrkName = "Tech-"+aTechnology->getName();

        // create the calibration market of type subsidy since we will only have access to the required investment
        // now and so need to make sure it does not get cleared from the market
        marketplace->createMarket( mCurrentRegionName, mCurrentRegionName, calibrationMrkName, IMarketType::SUBSIDY );
        marketplace->addToDemand( calibrationMrkName, mCurrentRegionName, technologyInvestment, aPeriod, true );
        
        // set the market to solve if this is not the first one which we want to hold at a constant
        // share weight and the technology is fixed investment
        if( mTechCount == 0 || aTechnology->getFixedInvestment( aPeriod ) > 0 ) {
            shareWeight = 1.0;
        }
        else {
            marketplace->setMarketToSolve( calibrationMrkName, mCurrentRegionName, aPeriod );
        }
        marketplace->setPrice( calibrationMrkName, mCurrentRegionName, shareWeight, aPeriod, true );
        ++mTechCount;
    }
}

void InvestableCounterVisitor::endVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                       const int aPeriod )
{
    // we don't have to bother clearing the flags and counters since they will get
    // overwritten by the subsector visit anyways
}
