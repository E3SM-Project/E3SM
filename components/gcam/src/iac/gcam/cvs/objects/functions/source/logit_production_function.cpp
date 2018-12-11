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
* \file logit_production_function.cpp
* \ingroup Objects
* \brief The LogitProductionFunction class source file.
* \author Pralit Patel
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/logit_production_function.h"

using namespace std;

double LogitProductionFunction::calcCoefficient( InputSet& input, double consumption, const std::string& regionName,
                            const std::string& sectorName, int period, double sigma, double IBT,
                            double capitalStock ) const
{
    // period is really the period which we are setting the coef in
    // we want to use the base period prices and demands
    const int basePeriod = 0;
    // the first input is the parent node input
    InputSet::iterator it = input.begin();
    double totalDemand = 0;
    for( ++it; it != input.end(); ++it ) {
        // get physical is really giving us currency at this point so divide
        // by the price to get physical
        totalDemand += (*it)->getPhysicalDemand( basePeriod ) / (*it)->getPricePaid( regionName, basePeriod );
    }

    // we need to set this price back into the node because the node price does matter
    // when mixing logit and CES so we cannot leave this at 1
    // TODO: make sure this is true
    it = input.begin();
    (*it)->setPricePaid( (*it)->getPhysicalDemand( basePeriod ) / totalDemand, basePeriod );

    // get the first input which we will make the rest relative to
    // TODO: we could get the biggest input and make the rest relative
    // to it that way the coefs will have meaning
    ++it;
    double relativeFromPrice = (*it)->getPricePaid( regionName, basePeriod );
    double relativeFromShare = ( (*it)->getPhysicalDemand( basePeriod ) / (*it)->getPricePaid( regionName, basePeriod ) )
        / totalDemand;
    for( ; it != input.end(); ++it ) {
        double currShare = ( (*it)->getPhysicalDemand( basePeriod ) / (*it)->getPricePaid( regionName, basePeriod ) )
            / totalDemand;
        double shareWeight = ( currShare / relativeFromShare )
            * pow( relativeFromPrice / (*it)->getPricePaid( regionName, basePeriod ), sigma );
        (*it)->setCoefficient( shareWeight, period );
        (*it)->setCoefficient( shareWeight, basePeriod );
    }
    return 1;
}

double LogitProductionFunction::changeElasticity( InputSet& input, const std::string& aRegionName, double priceReceived,
                             double aProfits, double capitalStock, const int aPeriod, double alphaZero,
                             double sigmaNew, double sigmaOld ) const
{
    // TODO: this probably should be implemented
    return 0;
}

double LogitProductionFunction::calcDemand( InputSet& input, double consumption, const std::string& regionName,
                       const std::string& sectorName, const double aShutdownCoef, int period,
                       double capitalStock, double alphaZero, double sigma, double IBT ) const
{
    double totalSum = 0;
    for( InputSet::const_iterator it = input.begin(); it != input.end(); ++it ) {
        totalSum += (*it)->getCoefficient( period ) * pow( (*it)->getPricePaid( regionName, period ), sigma );
    }

    double totalDemand = 0.0;
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        // consumption is really output at the current node
        double shareDemand = (*it)->getCoefficient( period ) * pow( (*it)->getPricePaid( regionName, period ), sigma )
            / totalSum * consumption;
        (*it)->setPhysicalDemand( shareDemand, regionName, period );
        totalDemand += shareDemand;
    }

    // assert everything got distributed
    // TODO: this was off by a very small number so use util::isEqual with a reasonable
    // tolerance instead
    //assert( totalDemand == consumption );
    return totalDemand;
}

double LogitProductionFunction::calcLevelizedCost( const InputSet& aInputs, const std::string& aRegionName,
                         const std::string& aSectorName, int aPeriod, double aAlphaZero, double aSigma ) const
{
    double totalSum = 0;
    for( InputSet::const_iterator it = aInputs.begin(); it != aInputs.end(); ++it ) {
        totalSum += (*it)->getCoefficient( aPeriod ) * pow( (*it)->getPricePaid( aRegionName, aPeriod ), aSigma );
    }

    double totalCost = 0.0;
    for( InputSet::const_iterator it = aInputs.begin(); it != aInputs.end(); ++it ) {
        // consumption is really output at the current node
        double share = (*it)->getCoefficient( aPeriod ) * pow( (*it)->getPricePaid( aRegionName, aPeriod ), aSigma )
            / totalSum;
        totalCost += share * (*it)->getPricePaid( aRegionName, aPeriod );
    }
    return totalCost;
}
