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
 * \file set_share_weight_visitor.cpp
 * \ingroup Objects
 * \brief The SetShareWeightVisitor class source file.
 * \author Sonny Kim
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "investment/include/set_share_weight_visitor.h"
#include "technologies/include/base_technology.h"
#include "util/base/include/util.h"
#include "sectors/include/subsector.h"
#include "sectors/include/sector.h"
#include "marketplace/include/marketplace.h"
 
using namespace std;
extern Scenario* scenario;

/*!
 *\brief Constructor
 * \param aRegionName Name of the region if starting the visiting below the
 *        region level.
 */
SetShareWeightVisitor::SetShareWeightVisitor( const string& aRegionName )
:mCurrentRegionName( aRegionName )
{
}

void SetShareWeightVisitor::startVisitSector( const Sector* aSector,
                                                 const int aPeriod )
{
    mCurrentSectorName = aSector->getName();
}

void SetShareWeightVisitor::endVisitSector( const Sector* aSector,
                                               const int aPeriod )
{
    mCurrentSectorName.clear();
}

void SetShareWeightVisitor::startVisitSubsector( const Subsector* aSubsector,
                                                    const int aPeriod )
{
    // if the subsector has a calibration market get the price from it and
    // set it into the subsector
    mCurrentSubSectorName = aSubsector->getName();
    if ( aSubsector->hasCalibrationMarket() ){
        Marketplace* marketplace = scenario->getMarketplace();
        string calibrationMrkName = "SubSec-"+mCurrentSubSectorName;
        double trialShareWeight = marketplace->getPrice( calibrationMrkName, mCurrentRegionName, aPeriod, true );
        // aSubsector is const
        const_cast<Subsector*>(aSubsector)->mShareWeights[ aPeriod ] = trialShareWeight;
    }
}

void SetShareWeightVisitor::endVisitSubsector( const Subsector* aSubsector,
                                                  const int aPeriod )
{
    mCurrentSubSectorName.clear();
}

void SetShareWeightVisitor::startVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                         const int aPeriod )
{
    // if the technology has a calibration market get the price from it and
    // set it into the technology
    if ( aTechnology->hasCalibrationMarket() ){
        Marketplace* marketplace = scenario->getMarketplace();
        string calibrationMrkName = "Tech-"+aTechnology->getName();
        double trialShareWeight = marketplace->getPrice( calibrationMrkName, mCurrentRegionName, aPeriod, true );
        //aTechnology is const
        const_cast<BaseTechnology*>(aTechnology)->mShareWeight = trialShareWeight;
    }
}

void SetShareWeightVisitor::endVisitBaseTechnology( const BaseTechnology* aTechnology,
                                                       const int aPeriod )
{
}
