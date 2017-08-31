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
* \file leontief_production_function.cpp
* \ingroup Objects
* \brief The LeontiefProductionFunction class source file.
*
*  Detailed description.
*
* \author Pralit Patel
* \author Sonny Kim
* \author Josh Lurz
* \date $Date: 2005/06/01 21:23:59 $
* \version $Revision: 1.2 $
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/leontief_production_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

//! Calculate the capital scaler.
double LeontiefProductionFunction::calcCapitalScaler( const InputSet& input, double aAlphaZero, 
                                                      double sigma, double capitalStock, const int aPeriod ) const 
{
    double capitalCoef = FunctionUtils::getCapitalInput( input )->getCoefficient( aPeriod );
    double qCapital = 0;
    if( capitalCoef != 0 && aAlphaZero != 0 ){
        qCapital = 1 / capitalCoef / aAlphaZero;
    }
    assert( util::isValidNumber( qCapital ) );
    return qCapital;
}

/*! \brief Calculate Leontief coefficients
* \return alpha zero.
*/
double LeontiefProductionFunction::calcCoefficient( InputSet& input, double consumption, 
												    const string& regionName, const string& sectorName, 
                                                    int period, double sigma, double indBusTax, 
                                                    double capitalStock ) const 
{
    assert( sigma < 0.05 ); // Should be a CES.
    // If there is no capital stock, initialize technical coefficients to zero.
    if( capitalStock < 1 ){
        for( unsigned int i = 0; i < input.size(); ++i ){
            input[ i ]->setCoefficient( 0, period );
        }
    }
    
    // Calculate output quantity.
    double demandCurrencySum = 0; // find real name.
    for( unsigned int i = 0; i < input.size(); ++i ){
        demandCurrencySum += input[ i ]->getCurrencyDemand( period );
    }
    
    // use price received(psubj) for the good in the next equation
    double priceReceived = FunctionUtils::getPriceReceived( regionName, sectorName, period );
    assert( priceReceived > 0 );
    double qSubJ = demandCurrencySum / priceReceived;
    assert( qSubJ > 0 );

    // Set the coefficients. 
    for( unsigned int i = 0; i < input.size(); ++i ){
		if( input[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            input[ i ]->setCoefficient( capitalStock / qSubJ, period );
        }
        else {
            if( input[ i ]->getPricePaid( regionName, period ) > 0 ){
                double currToPrice = input[ i ]->getCurrencyDemand( period )
					                 / input[ i ]->getPricePaid( regionName, period );
                input[ i ]->setCoefficient( currToPrice / qSubJ, period );
            }
            else {
                input[ i ]->setCoefficient( 0, period );
            }
        }
    }
    // Initialize alpha zero to unity.
    return 1;
}

/*! \brief Transform production coefficients based on long-term elasticity to
*          short-term elasticity.
* \author Sonny Kim
*/
double LeontiefProductionFunction::changeElasticity( InputSet& input, const string& aRegionName, 
													 double priceReceived, double aProfits, double capitalStock,
													 const int aPeriod, double alphaZero, double sigmaNew,
													 double sigmaOld ) const 
{
    assert( sigmaNew > 0 );
    assert( aProfits > 0 );
    const double rhoNew = FunctionUtils::getRho( sigmaNew );
    const double alphaZeroExp = rhoNew / ( 1 - rhoNew );

    for	( unsigned int i = 0; i	< input.size();	++i	) {
		// all inputs other than capital
        double priceRatio = 0;
        if(	!input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
			priceRatio = priceReceived / input[i]->getPricePaid( aRegionName, aPeriod );
        }
		// for capital or OVA input
		else {
            priceRatio = priceReceived / (aProfits/capitalStock);
		}
        double newCoef = pow( alphaZero, alphaZeroExp )
                         * pow( input[i]->getCoefficient( aPeriod ), sigmaNew )
                         * pow( priceRatio, sigmaNew );

		// set new transformed coefficient in input object
        input[i]->setCoefficient( newCoef, aPeriod );
    }
    
    // Set alpha zero to 1.
	return 1;
}

/*! \brief Calculate Demand.
* \return Returns total demand
*/
double LeontiefProductionFunction::calcDemand( InputSet& input, double personalIncome, 
											   const string& regionName, const string& sectorName,
                                               const double aShutdownCoef,
                                               int period, double capitalStock, double alphaZero, 
                                               double sigma, double IBT ) const 
{
    double qCapital = calcCapitalScaler( input, alphaZero, sigma, capitalStock, period );

    double qTemp = aShutdownCoef * qCapital * capitalStock;
    
    // Calculate direct demands. Not sure why this doesn't use scaleFactor.
    double totalDemand = 0;
    for( unsigned int i = 0; i < input.size(); ++i ){
        if( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            double currDemand = qTemp * input[ i ]->getCoefficient( period ) / alphaZero;
            totalDemand += currDemand;
            if( currDemand < 0 ){
				ILogger& mainLog = ILogger::getLogger( "main_log" );
				mainLog.setLevel( ILogger::WARNING );
                mainLog << "Trying to add negative demand currency for " << input[ i ]->getName() 
					    << " in " << sectorName << endl;
            }
            input[ i ]->setCurrencyDemand( currDemand, regionName, period );
        }
    }
    assert( totalDemand >= 0 );
    assert( util::isValidNumber( totalDemand ) );
    return totalDemand;
}

//! Calculate Expected Profit Rate.
double LeontiefProductionFunction::calcExpProfitRate( const InputSet& input, const string& regionName, 
													  const string& sectorName, double aLifetimeYears,
													  int period, double alphaZero, double sigma ) const 
{
    double qCapital = calcCapitalScaler( input, alphaZero, sigma, 0, period );
    double valQ = qCapital * FunctionUtils::getExpectedPriceReceived( input, regionName, sectorName,
                                                                      aLifetimeYears, period );
    assert( valQ >= 0 );
    assert( util::isValidNumber( valQ ) );
    const double netPresentValueMult = FunctionUtils::getNetPresentValueMult( input, regionName,
																			  aLifetimeYears, period );
    double sum = 0;
    for( unsigned int i = 0; i < input.size(); ++i ){
        if( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            sum += qCapital / alphaZero * input[ i ]->getCoefficient( period ) * netPresentValueMult;
        }
    }
    assert( sum >= 0 );
    assert( util::isValidNumber( sum ) );
    return max( valQ - sum, 0.0 );
}

/*! \brief Calculate the levelized cost for this technology.
* \details 
* 
* \author Josh Lurz
* \param aInputs Vector of inputs for the technology.
* \param aRegionName Name of the region containing the production function.
* \param aSectorName Name of the sector containing the production function.
* \param aPeriod Period in which to calculate levelized cost.
* \param aAlphaZero Out front scalar.
* \param aSigma Sigma coefficient.
* \return The levelized cost.
*/
double LeontiefProductionFunction::calcLevelizedCost( const InputSet& aInputs, const string& aRegionName,
													  const string& aSectorName, int aPeriod,
													  double aAlphaZero, double aSigma ) const
{
    // Loop through the inputs and add on their costs.
    double levelizedCost = 0;
    for( unsigned int i = 0; i < aInputs.size(); ++i ){
        // Capital is done specially.
        if( aInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            levelizedCost += aInputs[ i ]->getCoefficient( aPeriod ) 
                             * aInputs[ i ]->getPrice( aRegionName, aPeriod )
							 // I think this should be + not *.
                             * aInputs[ i ]->getPriceAdjustment();
        }
        else {
            levelizedCost += aInputs[ i ]->getCoefficient( aPeriod ) 
				             * aInputs[ i ]->getPricePaid( aRegionName, aPeriod );
        }
    }
    return levelizedCost / aAlphaZero;
}

//! Calculate profits.
double LeontiefProductionFunction::calcUnscaledProfits( const InputSet& input, const string& regionName,
														const string& sectorName, const int period, 
													    const double capitalStock, const double alphaZero, 
                                                        const double sigma ) const 
{
    double qCapital = calcCapitalScaler( input, alphaZero, sigma, capitalStock, period );

    double valQ = qCapital * FunctionUtils::getPriceReceived( regionName, sectorName, period );
    
    assert( valQ >= 0 );
    assert( util::isValidNumber( valQ ) );

    double sum = 0;
    for( unsigned int i = 0; i < input.size(); ++i ){
        if( !input[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            sum += qCapital * input[ i ]->getCoefficient( period ) / alphaZero
				   * input[ i ]->getPricePaid( regionName, period );
        }
    }
    assert( sum >= 0 );
    assert( util::isValidNumber( sum ) );
    return max( ( valQ - sum ) * capitalStock, 0.0 );
}

//! Calculate output.
double LeontiefProductionFunction::calcOutput( InputSet& input, const string& regionName, 
											   const string& sectorName, const double aShutdownCoef,
                                               int period, double capitalStock, 
                                               double alphaZero, double sigma ) const 
{
    double qCapital = calcCapitalScaler( input, alphaZero, sigma, capitalStock, period );
    return aShutdownCoef * qCapital * capitalStock;
}

//! Return the amount of output produced by one unit of capital.
double LeontiefProductionFunction::getCapitalOutputRatio( const InputSet& aInputs, const string& aRegionName,
														  const string& aSectorName, double aLifetimeYears, 
                                                          int aPeriod, double aAlphaZero,
                                                          double aSigma ) const
{
    // this is actually 1 / ratio, might be better to return the ratio.
	return aAlphaZero * FunctionUtils::getCapitalInput( aInputs )->getCoefficient( aPeriod );
}

/*! \brief Apply technical change to production functions.
* \details 
* \note This function currently makes a call to
*       FunctionUtils::applyTechChangeInternal so that it can share code with
*       ADemandFunction::applyTechnicalChange. In the future these
*       implementations may vary so that function is not called directly.
* \param input Vector of inputs for the demand function.
* \param aTechChange A structure containing the various possible types of
*        technical change.
* \param regionName Name of the region containing the function.
* \param sectorName Nmae of the sector containing the function.
* \param alphaZero The up-front scaler.
* \param sigma Sigma coefficient.
* \return The new alpha zero.
*/
double LeontiefProductionFunction::applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                                         const string& regionName, const string& sectorName,
                                                         const int aPeriod, double alphaZero, double sigma ) const 
{
    return FunctionUtils::applyTechnicalChangeInternal( input, aTechChange, regionName, sectorName,
                                                        aPeriod, alphaZero, sigma );
}
