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
* \file nested_ces_production_function.cpp
* \ingroup Objects
* \brief The NestedCESProductionFunction class source file.
* \author Pralit Patel
* \author Ron Sands
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/nested_ces_production_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

double NestedCESProductionFunction::calcCoefficient( InputSet& input, double consumption, const std::string& regionName,
                            const std::string& sectorName, int period, double sigma, double IBT,
                            double capitalStock ) const
{
    // the first input is the parent node input
    InputSet::iterator it = input.begin();
    const double r = 1 - sigma;
    // this is really currency
    const double totalDemand = (*it)->getPhysicalDemand( period );
    const double outputPrice = (*it)->getPrice( regionName, period );
    // period is really the period which we are setting the coef in
    // we want to use the base period prices and demands
    const int basePeriod = 0;
    for( ++it; it != input.end(); ++it ) {
        double priceRatio = (*it)->getPricePaid( regionName, basePeriod ) / outputPrice;
        // this will really give currency demand
        double tempCoef = pow1( ( (*it)->getPhysicalDemand( basePeriod ) / totalDemand ), -1 / r ) *
            priceRatio;
        (*it)->setCoefficient( tempCoef, period );
        (*it)->setCoefficient( tempCoef, basePeriod );
    }
    return 1.0;
}

double NestedCESProductionFunction::applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                 const std::string& regionName,const std::string& sectorName, const int aPeriod, 
                                 double alphaZero, double sigma ) const
{
    // TODO: should I just use FunctionUtils::applyTechnicalChangeInternal
    int timeStep = scenario->getModeltime()->gettimestep( aPeriod );
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        double techChange = FunctionUtils::getTechChangeForInput( *it,
                                                                  aTechChange,
                                                                  aPeriod );
                                                                  
        if( techChange != 0 ) {
            double scaleFactor = pow1( 1 + techChange, timeStep );
            double newCoef = (*it)->getCoefficient( aPeriod ) * scaleFactor;
            (*it)->setCoefficient( newCoef, aPeriod );
        }
    }
    // do I need to do anything for alpha zero?
    return 0;
}

double NestedCESProductionFunction::changeElasticity( InputSet& input, const std::string& aRegionName, double priceReceived,
                             double aProfits, double capitalStock, const int aPeriod, double alphaZero,
                             double sigmaNew, double sigmaOld ) const
{
    // sigmaNew means sigma new capital conversly old is old capital and we change from
    // new to old which may seem counter intuitive if you were thinking new is the one
    // we will be using and old was the one we used to use
    const double alphaExp = ( sigmaNew - 1 ) / ( sigmaOld - 1 );
    const double priceRatioExp = ( sigmaNew - sigmaOld ) / ( sigmaOld - 1 );
    const int basePeriod = scenario->getModeltime()->getBasePeriod();
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        // Note aProfits is a hack to indicate wether we were intending on changing sigmas from new capital
        // to old capital (in which case aProfits would be zero and we are changing the coef from the previous
        // period), or chaging the sigmas over time (in which case aProfits would be 1 and we are changing the coef
        // which was set in the current period after being passed forward).
        double newCoef = pow1( alphaZero * (*it)->getCoefficient( aProfits == 0 ? aPeriod - 1 : aPeriod ), alphaExp ) *
            pow1( priceReceived / (*it)->getPricePaid( aRegionName, basePeriod  ), priceRatioExp );
        (*it)->setCoefficient( newCoef, aPeriod );
    }
    return 0;
}

double NestedCESProductionFunction::calcDemand( InputSet& input, double consumption, const std::string& regionName,
                       const std::string& sectorName, const double aShutdownCoef, int period,
                       double capitalStock, double alphaZero, double sigma, double IBT ) const
{   
    // IBT is really output price
    const double outputPrice = IBT;
    double totalDemand = 0.0;
    for( InputSet::iterator it = input.begin(); it != input.end(); ++it ) {
        double ioRatio = calcIORatio( *it, regionName, period, alphaZero, outputPrice, sigma );
        // consumption is really output at the current node (the path io ratio has already been applied)
        double tempDemand = ioRatio * consumption * aShutdownCoef;
        (*it)->setPhysicalDemand( tempDemand, regionName, period );
        totalDemand += tempDemand;
    }
    return totalDemand;
}

double NestedCESProductionFunction::calcOutput( InputSet& input, const std::string& regionName,
                       const std::string& sectorName, const double aShutdownCoef,
                       int period, double capitalStock, double alphaZero, double sigma ) const
{
    // TODO: everything
    assert( false );
    return 0;
}

double NestedCESProductionFunction::getCapitalOutputRatio( const InputSet& aInputs, const std::string& aRegionName,
                                  const std::string& aSectorName, double aLifeTimeYears, int aPeriod,
                                  double aAlphaZero, double aSigma ) const
{
    // we are only expecting a single capital input here
    assert( aInputs.size() == 1 );
    IInput* capitalInput = aInputs.front();
    //assert( capitalInput->hasTypeFlag( IInput::CAPITAL ) );

    // aLifeTimeYears is really the node price
    return calcIORatio( capitalInput, aRegionName, aPeriod, aAlphaZero, aLifeTimeYears, aSigma );
}

double NestedCESProductionFunction::calcExpProfitRate( const InputSet& input, const std::string& regionName,
                              const std::string& sectorName, double aLifeTimeYears, int period, double alphaZero,
                              double sigma ) const
{
    // TODO: figure out what to do about all of this since it is no longer used
    // TODO: this is sort of unecessary going to through to find the capital twice.. might just have to live with it
    IInput* capInput = FunctionUtils::getCapitalInput( input );
    if( !capInput ) {
        return 0;
    }
    double profitRate = capInput->getPricePaid( regionName, period );
    // this is the real discount rate
    capInput->setPricePaid( capInput->getPrice( regionName, period ) + capInput->getPriceAdjustment(),
        period );
    // Calculate the net present value multiplier to determine expected prices.
    const double netPresentValueMult = FunctionUtils::getNetPresentValueMult( input, regionName, aLifeTimeYears, period );
    // set the profit rate back in capital
    capInput->setPricePaid( profitRate, period );
    return profitRate * netPresentValueMult;
}

double NestedCESProductionFunction::calcLevelizedCost( const InputSet& aInputs, const std::string& aRegionName,
                         const std::string& aSectorName, int aPeriod, double aAlphaZero, double aSigma ) const
{
    // TODO: shouldn't aPeriod, aAlphaZero, and aSigma all be const?

    double r = 1 - aSigma;
    double currLevelizedCost = 0.0;
    for( InputSet::const_iterator inputIter = aInputs.begin(); inputIter != aInputs.end(); ++inputIter ) {
        double currPrice = (*inputIter)->getPricePaid( aRegionName, aPeriod );
        currLevelizedCost += pow1( currPrice / (*inputIter)->getCoefficient( aPeriod ), r );
    }
    currLevelizedCost = pow1( currLevelizedCost, 1 / r ) / aAlphaZero;
    
    return currLevelizedCost;
}

double NestedCESProductionFunction::calcUnscaledProfits( const InputSet& aInputs, 
                                const std::string& aRegionName,
                                const std::string& aSectorName,
                                const int aPeriod,
                                const double aCapitalStock,
                                const double aAlphaZero,
                                const double aSigma ) const
{
    // TODO: this is a hack aSigma is really the output carbon tax adjustment
    const double priceRecieved = FunctionUtils::getPriceReceived( aRegionName, aSectorName, aPeriod ) + aSigma;
    const double cost = calcCosts( aInputs, aRegionName, aAlphaZero, aPeriod );
    // note profits can be negative
    return ( aCapitalStock * priceRecieved ) - cost; // aCapitalStock is really output
}

double NestedCESProductionFunction::calcCapitalScaler( const InputSet& input, double aAlphaZero, double sigma,
                              double capitalStock, const int aPeriod ) const
{
    // TODO: do we need this?
    assert( false );

    return 0;
}

double NestedCESProductionFunction::calcIORatio( const InputSet::value_type aInput, const string& aRegionName, const int aPeriod,
                                const double aAlphaZero, const double aParentPrice, const double aSigma ) const
{
    // optimazation to speed up calculations when sigma is zero
    if( aSigma == 0 ) {
        return 1 / aInput->getCoefficient( aPeriod );
    }
    // TODO: price paid should be > 0 put the assert back in
    double pricePaid = aInput->getPricePaid( aRegionName, aPeriod );
    //assert( pricePaid >= 0 );
    double ret = pow( aAlphaZero * aInput->getCoefficient( aPeriod ), aSigma - 1 ) * 
        pow( aParentPrice / pricePaid, aSigma );
    return pricePaid != 0 ? ret : 0;
}

/*!
 * \brief A wrapper around to pow function intended to provide a performance boost.
 * \details  When using elasticities of 0 for the CES function, which can happen often,
 *           we end up trying to raise somthing to the power of 1 of which the result
 *           will of course be the same value.  In order to avoid the time costly call
 *           to pow we will just check for the exponent of 1 explicitly before hand.
 * \param base A value to be raised to a power.
 * \param exp The power to raise by.
 * \return base ^ exp.
 */
inline double NestedCESProductionFunction::pow1( double base, double exp ) const {
    return exp != 1 ? pow( base, exp ) : base;
}
