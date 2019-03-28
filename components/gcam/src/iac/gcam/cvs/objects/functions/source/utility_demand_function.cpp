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
* \file utility_demand_function.cpp
* \ingroup Objects
* \brief The UtilityDemandFunction class source file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/utility_demand_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/model_time.h"

using namespace std;

extern Scenario* scenario;

//! Calculate Demand
double UtilityDemandFunction::calcDemand( InputSet& input, double personalIncome, 
                                            const string& regionName, const string& sectorName,
                                            const double aShutdownCoef,
                                            int period, double capitalStock, double alphaZero, 
                                            double sigma, double IBT ) const 
{
    // KVC: capitalStock is the mUtilityParameterA from the consumer
    double A = capitalStock;

    const string utilityMarketName = sectorName+"-utility";
    Marketplace* marketplace = scenario->getMarketplace();
    double totalUtility = 0;
    const double trialUtility = marketplace->getPrice( utilityMarketName, regionName, period, true );
    double totalDemand = 0; // total demand used for scaling
    double availablePersonalIncome = personalIncome;
    
    // subtract off minimum demand to get the available income
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        availablePersonalIncome -= (*it)->getPricePaid( regionName, period ) *
            (*it)->getCoefficient( period );
    }
    
    // calc the demands and utilites for each input
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        // TODO: income elasticity is really alpha, price elasticity is really beta,
        //       and getCoefficient is really gamma
        assert( (*it)->getPricePaid( regionName, period ) > 0 );
        assert( (*it)->getPriceElasticity() >= 0 );
        assert( (*it)->getIncomeElasticity() >= 0 );
        
        double demand = (*it)->getCoefficient( period ) + ( 1 / (*it)->getPricePaid( regionName, period ) )
            * phi( (*it), trialUtility ) * availablePersonalIncome;
        assert( util::isValidNumber( demand ) );
        (*it)->setPhysicalDemand( demand, regionName, period );
        
        double inputUtility = 0;
        if ( demand - (*it)->getCoefficient( period ) > 0 ) {
            // TODO: consider storing this somewhere
            inputUtility = phi( (*it), trialUtility ) * log( ( demand - (*it)->getCoefficient( period ) )
                / ( A * /*exp(*/ trialUtility /*)*/ ) ); // the solver is working in e^u see the g function
        }
            
        totalDemand += demand * (*it)->getPricePaid( regionName, period );
        totalUtility += inputUtility;
    }
    
    /*!
     * \pre The market supply and demand is zero as this is the only place it is added
     */
    assert( marketplace->getSupply( utilityMarketName, regionName, period ) 
        < util::getVerySmallNumber() );
    assert( marketplace->getDemand( utilityMarketName, regionName, period ) 
        < util::getVerySmallNumber() );
    // the constraint for input utility is that it should sum to 1
    marketplace->addToSupply( utilityMarketName, regionName, 1, period, true );
    marketplace->addToDemand( utilityMarketName, regionName, totalUtility, period, true );
    return totalDemand;
}

double UtilityDemandFunction::calcCoefficient( InputSet& input, double consumption,
                                                 const string& regionName, const string& sectorName,
                                                 int period, double sigma, double IBT,
                                                 double capitalStock ) const
{
    // TODO: this does not actually calculate coefficients however that could be something
    // we do if they are not read in and just assume LES

    // back out our initial utility, to do this we can use our values to back out
    // a phi for an input and then back out the utility from that
    const int basePeriod = 0;
    // the first input is the parent node of the other inputs
    InputSet::iterator it = input.begin();
    // first calculate our available income past the subsitance level
    // physical demand is really currency at this point
    double availablePersonalIncome = (*it)->getPhysicalDemand( basePeriod );
    // subtract off minimum demand to get the available income
    for( ++it; it != input.end(); ++it ) {
        // getCoefficient is really gamma here
        availablePersonalIncome -= (*it)->getPricePaid( regionName, basePeriod ) *
            (*it)->getCoefficient( basePeriod );
    }

    // pick an input to back out phi and u
    it = input.begin() + 1;
    assert( it != input.end() );
    // physical demand is really currency at this point
    const double phi = ( (*it)->getPhysicalDemand( basePeriod ) - ( (*it)->getPricePaid( regionName, basePeriod ) * 
        (*it)->getCoefficient( basePeriod ) ) ) / availablePersonalIncome;
    // price elasticity is really alpha and income elasticity is really beta
    const double initialUtility = ( (*it)->getPriceElasticity() - phi ) / ( phi - (*it)->getIncomeElasticity() );
    // note that initialUtility is really g( u ) however solve is currently working in e^u so this is what
    // we want to set as the price
    const string utilityMarketName = sectorName+"-utility";
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->setPrice( utilityMarketName, regionName, initialUtility, basePeriod );

    return 1; // return null for CES AlphaZero only
}

double UtilityDemandFunction::applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                              const string& regionName, const string& sectorName,
                                              const int aPeriod, double alphaZero, double sigma ) const 
{
    // do nothing
    return 0;
}

/*!
 * \brief Calculate the phi function at the given utility level.
 * \details This method is really just calculating a share at a given utility
 *          that is between the lower bound alpha, and the upper bound beta also
 *          g() is just some smooth function.
 * \param aInput The input for which we are calculating a share for.
 * \param aUtility The utility level at which we are calculating the share.
 * \return The input's share of the demand at a given utility level.
 */
double UtilityDemandFunction::phi( const InputSet::value_type& aInput, const double aUtility ) const {
    // TODO: income elasticity is really beta and price elasticity is really alpha
    return ( aInput->getPriceElasticity() + aInput->getIncomeElasticity() * g( aUtility ) )
        / ( 1 + g( aUtility ) );
}

/*!
 * \brief The g() function to be used when calculating phi.
 * \details This function needs to be a smooth function to define the behavoir
 *          between alpha and beta.  I think it must be twice differentiable but
 *          I should double check.
 * \param aUtility The utility level to calculate the function at.
 * \note Since the solver needs to work in only positive numbers aUtility is really
 *       e^u so g will just return that.
 */
double UtilityDemandFunction::g( const double aUtility ) const {
    // currently the solver is working in e^u since it needs to work in
    // only positive values
    return aUtility; //exp( aUtility );
}
