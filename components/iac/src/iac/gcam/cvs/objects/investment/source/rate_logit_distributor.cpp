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
 * \file rate_logit_distributor.cpp
 * \ingroup Objects
 * \brief RateLogitDistributor class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <cassert>
#include <cmath>
#include <numeric>
#include "investment/include/rate_logit_distributor.h"
#include "investment/include/iinvestable.h"
#include "investment/include/iexpected_profit_calculator.h"
#include "investment/include/investment_utils.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"

using namespace std;

/*! \brief Constructor
*/
RateLogitDistributor::RateLogitDistributor() {
}

/*! \brief Distribute investment based on the expected profit rates of the
*          children and a logit.
* \details Distributes a level of investment passed from the level above between
*          the aInvestable objects. The following steps are performed:
*          <ol><li>Calculate the expected profit rate for the object containing the
*          investable children, the sector or subsector</li>
*          <li>Use each child's expect profit level, the overall expected profit level,
*          and a logit to determine an unnormalized share of investment</li>
*          <li>Normalize the investment shares to 1.</li>
*          <li>Sum the fixed investment of all investable children</li>
*          <li>Calculate the amount of variable investment which will be distributed by
*          subtracting the fixed amount from the quantity passed in.</li> <li>Distribute
*          a quantity of investment to each investable child equal to the investable
*          childs share multiplied by the total variable investment.</li>
*          <li>Check that the total amount distributed is equal to the original amount passed in.</li>
* \param aRateCalc An expected profit rate calculation object.
* \param aInvestables A vector of children which will receive investment.
* \param aNationalAccount The national account for this region.
* \param aInvestmentExp Logit exp used to distribute. TODO: remove if
*        can think of a better way.
* \param aRegionName The name of the containing region.
* \param aSectorName The name of the containing sector.
* \param aAmount The amount of investment to distribute.
* \param aPeriod The period in which to distribute investment.
* \return The total amount of investment actually distributed.
* \author Josh Lurz
*/
double RateLogitDistributor::distribute( const IExpectedProfitRateCalculator* aRateCalc,
                                         vector<IInvestable*>& aInvestables,
                                         NationalAccount& aNationalAccount,
                                         const double aInvestmentExp,
                                         const string& aRegionName,
                                         const string& aSectorName,
                                         const double aAmount,
                                         const int aPeriod ) const 
{
    // Calculate the normalized investment shares.
    // Create a vector of investment shares, one per subsector.
    vector<double> shares = calcInvestmentShares( aInvestables, aRateCalc, 
                                                  aNationalAccount, aInvestmentExp,
                                                  aRegionName, aSectorName, aPeriod );
    // Should be one share per investable.
    assert( shares.size() == aInvestables.size() );

    // Calculate the sum of the shares for error checking.
    const double sumShares = accumulate( shares.begin(), shares.end(), 0.0 );
    
    // Calculate the variable amount of investment.
    const double totalFixedAmount = InvestmentUtils::sumFixedInvestment( aInvestables, aPeriod );
    const double varInvestment = aAmount - totalFixedAmount;

    // Check if the shares could not be normalized, which means they have
    // expected profit rates of zero, and there is variable investment which
    // should be distributed.
    if( !util::isEqual( sumShares, 1.0 ) && varInvestment > 0 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Variable investment will not be distributed as there are no children"
                << " with positive expected profit rates in sector " 
                << aSectorName << " in region " << aRegionName << "." << endl;
    }

    // Distribute new investment.
    double sumDistributed = 0;
    for( unsigned int i = 0; i < aInvestables.size(); ++i ){
        // If the subsector only has fixed investment the share will be zero.
        // Fixed investment will be used at the child level. This will pass zero
        // investment to the child for fixed investment.
        sumDistributed += aInvestables[ i ]->distributeInvestment( this,
                                                                   aNationalAccount,
                                                                   aRateCalc,
                                                                   aRegionName,
                                                                   aSectorName,
                                                                   varInvestment * shares[ i ],
                                                                   aPeriod );
    }
    // Check if the amount distributed is equal to that was passed in.
    if( !util::isEqual( sumDistributed, aAmount ) ){
      /*  ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Investment Problem!  Region: "<< aRegionName << "  Sector: " << aSectorName 
                << "  Requested Amount: " << aAmount << "  Distributed Amount: " << sumDistributed 
                << "  Difference: " << aAmount - sumDistributed << endl
                << "  Variable Amount: " << varInvestment << "  Sum of Fixed Amount: " << totalFixedAmount
                << endl; 
      */
    }
    return sumDistributed;
}

/*! \brief Calculate the amount of capital required to produce a unit of output.
* \details Determines the aggregate capital to output ratio of a group of
*          investables given that the shares used to calculate the ratio are the
*          same as those used to distribute investment. For this reason
*          investment distribution and calculation of capital to output ratios
*          should always use the same instance of this class. To calculate the
*          output ratio, the function first determines the set of shares that
*          will be used eventually to distribute investment. It then creates a
*          weighted average capital to output ratio using the shares as weights
*          on the capital to output ratios on the investable children.
* \param aInvestables A vector of children to calculate the capital to output
*        ratio for.
* \param aRateCalc The object responsible for calculating expected profit rates.
* \param aNationalAccount The national account for this region.
* \param aRegionName The name of the containing region.
* \param aSectorName The name of the containing sector.
* \param aPeriod The period in which to calculate the ratio.
* \return The capital to output ratio.
* \todo Check that this is correct for multiple technologies.
*/
double RateLogitDistributor::calcCapitalOutputRatio( const vector<IInvestable*>& aInvestables,
                                                     const IExpectedProfitRateCalculator* aRateCalc,  
                                                     const NationalAccount& aNationalAccount,
                                                     const string& aRegionName,
                                                     const string& aSectorName,
                                                     const int aPeriod ) const
{   
    // Create a vector of investment shares, one per subsector.
    // TODO: not correct although this should probably never get called.
    vector<double> shares = calcInvestmentShares( aInvestables, aRateCalc, 
                                                  aNationalAccount, 0, aRegionName,
                                                  aSectorName, aPeriod );
    // Should be one share per investable.
    assert( shares.size() == aInvestables.size() );

    // Now calculate the sector level average capital output ratio.
    double outputRatio = 0;
    for( unsigned int i = 0; i < aInvestables.size(); ++i ){
        const double capOutputRatio = aInvestables[ i ]->getCapitalOutputRatio( this, 
                                                                                aRateCalc,
                                                                                aNationalAccount,
                                                                                aRegionName,
                                                                                aSectorName,
                                                                                aPeriod );
        if( capOutputRatio > 0 ){
            outputRatio += shares[ i ] * 1 / capOutputRatio;
        }
    }
    // If the output ratio is positive, return the inverse. Otherwise return 0.
    return ( outputRatio > 0 ) ? 1 / outputRatio : 0;
}

/*! \brief Calculate a set of shares to distribute investment to a set of
*          investable children.
* \details Calculates a set of investment shares given a logit to shape the
*          distribution and an object which calculates expected profits. This is
*          accomplished by first calculating a parent level average expected
*          profit rate, and then calculating a share for each investable child
*          relative to that rate. These shares are then normalized to one.
* \param aInvestables A vector of children for which to calculate the investment
*        shares.
* \param aRateCalc The object responsible for calculating expected profit rates.
* \param aNationalAccount The national account for this region.
* \param aInvestmentExp An investment logit exponent.
* \param aRegionName The name of the containing region.
* \param aSectorName The name of the containing sector.
* \param aPeriod The period in which to calculate shares.
* \return A set of investment shares.
*/
const vector<double> RateLogitDistributor::calcInvestmentShares( const vector<IInvestable*>& aInvestables,
                                                   const IExpectedProfitRateCalculator* aRateCalc,
                                                   const NationalAccount& aNationalAccount,
                                                   const double aInvestmentExp,
                                                   const string& aRegionName,
                                                   const string& aSectorName,
                                                   const int aPeriod ) const
{
    assert( aRateCalc );
    // Calculate sector level expected profit.
    double expProfitRateTotal = aRateCalc->calcSectorExpectedProfitRate( aInvestables,
                                                                         aNationalAccount,
                                                                         aRegionName,
                                                                         aSectorName,
                                                                         aInvestmentExp,
                                                                         true,
                                                                         true,
                                                                         aPeriod );
    // Create a vector of investment shares, one per subsector.
    vector<double> shares( aInvestables.size(), 0 );
    
    // If there is no expected profit then the shares will stay at zero.
    if( expProfitRateTotal > 0 ){
        // Loop through subsectors first to calculate shares.
        for( unsigned int i = 0; i < aInvestables.size(); ++i ){
            double currExpProfitRate = aInvestables[ i ]->getExpectedProfitRate( aNationalAccount,
                                                                                 aRegionName,
                                                                                 aSectorName,
                                                                                 aRateCalc,
                                                                                 aInvestmentExp,
                                                                                 false,
                                                                                 true,
                                                                                 aPeriod );

            // Store the subsector investment share. This will be zero for negative profit subsectors.
            if( currExpProfitRate > 0 ) {
                shares[ i ] = aInvestables[ i ]->getShareWeight( aPeriod ) * pow( currExpProfitRate, aInvestmentExp ) 
                    / expProfitRateTotal;
            }
        }
    }
    // Normalize the investment shares.
    InvestmentUtils::normalizeShares( shares );

    return shares;
}
