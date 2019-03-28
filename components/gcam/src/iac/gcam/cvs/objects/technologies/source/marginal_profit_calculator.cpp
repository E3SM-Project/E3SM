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
 * \file vintage_production_state.cpp
 * \ingroup Objects
 * \brief VintageProductionState class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include "technologies/include/marginal_profit_calculator.h"
#include "technologies/include/technology.h"
#include "util/base/include/util.h"

using namespace std;

/*!
 * \brief Constructor.
 * \param aTechnology Technology for which to calculate marginal profits.
 */
MarginalProfitCalculator::MarginalProfitCalculator( const Technology* aTechnology )
: mTechnology( aTechnology )
{}

/*!
 * \brief Calculate the short term marginal profit as a proportion of variable
 *        costs.
 * \param aRegionName Region name.
 * \param aPeriod Model period.
 * \return Short term marginal profit as a proportion of variable costs.
 * \todo This calculation will have to be improved when a Technology has
 *       multiple and correctly differentiated fixed and variable costs. The
 *       code currently assumes that all non-energy costs are fixed, which would
 *       not be true for O&M for example.
 */
double MarginalProfitCalculator::calcShortTermMarginalProfit( const string& aRegionName,
                                                              const string& aSectorName,
                                                              const int aPeriod ) const
{
    // TODO: Marginal revenue has already deducted the ghg value.  If we could avoid
    // recalculating that value it could be a performance increase.
    double ghgValue = mTechnology->getTotalGHGCost( aRegionName, aSectorName, aPeriod );
    // The GHG cost should be added as a cost however it has already been included
    // in the revenue so we need to remove it from there first.
    double marginalRevenue = mTechnology->getMarginalRevenue( aRegionName,
                                                              aSectorName,
                                                              aPeriod ) + ghgValue;

    double variableCosts = mTechnology->getEnergyCost( aRegionName,
                                                       aSectorName,
                                                       aPeriod ) + ghgValue;
    

    /* Marginal profit is defined here as the percentage that price exceeds
    variable costs. Absolute value is used in the denominator since there is at least
    one situation (i.e., bio+CCS where the carbon credit exceeds the costs of using
    biomass) where, at least as the model defines it, variable costs can be negative.
    In this case, marginal profit will always be very positive, so there should be no
    problem introduced by using the absolute value in the denominator. In the vast
    majority of situations, variable costs will be positive and using absolute value
    will have no impact. A Very Small Number is also added to the denominator in case
	that the variable cost is zero, as it could be with renewables, which should
	not be shut down based on marginal profit either (MAW 8-6-09) */

    double marginalProfit = ( marginalRevenue - variableCosts ) / 
		                    ( fabs(variableCosts) + util::getVerySmallNumber() );
    return marginalProfit;
}
