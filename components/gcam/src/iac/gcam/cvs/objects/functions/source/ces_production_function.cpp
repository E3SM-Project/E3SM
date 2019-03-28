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
* \file ces_production_function.cpp
* \ingroup Objects
* \brief The CESProductionFunction class source file.
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

#include "functions/include/ces_production_function.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "util/base/include/model_time.h"
#include "functions/include/function_utils.h"
#include "util/logger/include/ilogger.h"

using namespace std;

extern Scenario* scenario;

// Calculate CES coefficients
double CESProductionFunction::calcCoefficient( InputSet& input, double consumption, 
                                               const string& regionName, const string& sectorName, 
                                               int period, double sigma, double indBusTax, 
                                               double capitalStock ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );
    assert( capitalStock > 0 );

    if( input.empty() ) {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Error no inputs while trying to calculate CES Coefficient" << endl;
    }
    // first part of calculating coefficient, divide IO value by price
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            double tempCoefficient = input[i]->getCurrencyDemand( period ) 
                                     / input[ i ]->getPrice( regionName, period );
            input[i]->setCoefficient( tempCoefficient, period );
        }
    }
    // get price of the numeraireIndex.
    const IInput* numeraireInput = FunctionUtils::getNumeraireInput( input );
    assert( numeraireInput );
    const double priceNumeraire = numeraireInput->getPrice( regionName, period );
    const double coefNumeraire = numeraireInput->getCoefficient( period );

    double total = 0;
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            total += input[i]->getCoefficient( period ) / coefNumeraire * input[ i ]->getPrice( regionName, period );
        }
    }

    // calculate alpha 0
    double priceOutput = FunctionUtils::getCurrencyDemandSum( input, period ) 
                         / ( FunctionUtils::getCurrencyDemandSum( input, period ) + indBusTax);
    double rho = FunctionUtils::getRho( sigma );
    double mu = rho/(1-rho);
    double alphaZero = pow( priceOutput/priceNumeraire,  -(1/rho) )
        * pow( numeraireInput->getCoefficient( period )
        / ( FunctionUtils::getCurrencyDemandSum( input, period )+indBusTax), (1/mu) );

    double Z = 1 - pow( (alphaZero * priceOutput), mu ) * total;

    // second part of calculation coefficient, normalize to numeraireIndex and apply elasticities
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            double tempCoefficient = pow( (input[i]->getCoefficient( period ) / coefNumeraire ),
                ( 1 / sigma ) ) * input[ i ]->getPrice( regionName, period ) / priceNumeraire;
            input[i]->setCoefficient( tempCoefficient, period );
        }
    }

    // calculate alpha for capital
    double capCoef = pow( (FunctionUtils::getCurrencyDemandSum(input, period )+indBusTax)/alphaZero/capitalStock, rho) * Z;
    IInput* capInput = FunctionUtils::getCapitalInput( input );
    assert( capInput );
    capInput->setCoefficient( capCoef, period );

    // third part of calculation coefficient, normalize to total and then to 100
    // normalize alphaZero
    alphaZero *= pow( ( FunctionUtils::getCoefSum(input, period )/100), (1/rho) );
    // normalize and set new coefficients
    FunctionUtils::scaleCoefficientInputs(input, 100/ FunctionUtils::getCoefSum(input, period ), period );
    // normalize alphaZero to 1
    alphaZero = normalizeAlphaZero( input, alphaZero, sigma, period );
    return alphaZero;
}

/*! \brief Normalize alpha zero scaler to 1 and readjust inputs
* \param 
* \author Sonny Kim
* \return alpha zero.
*/
double CESProductionFunction::normalizeAlphaZero( InputSet& input, double aAlphaZero,
                                                  double sigma, const int aPeriod ) const {
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    const double rho = FunctionUtils::getRho( sigma );
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        double newCoef = input[i]->getCoefficient( aPeriod ) * pow( aAlphaZero, rho );
        input[i]->setCoefficient( newCoef, aPeriod );
    }
    return 1;
}

/*! \brief Transform production coefficients based on long-term elasticity to
*          short-term elasticity.
* \param 
* \author Sonny Kim
* \return alpha zero.
*/
double CESProductionFunction::changeElasticity( InputSet& input, const string& aRegionName,
                                                double priceReceived, double aProfits, double capitalStock,
                                                const int aPeriod, double alphaZero, double sigmaNew, double sigmaOld ) const 
{
    // Note: This could actually happen and we should handle it by shutting down
    // the technology. I think we could just set alphaZero to zero.
    assert( aProfits > 0 );
    // Calculate parameters.
    const double rhoNew = FunctionUtils::getRho( sigmaNew );
    const double rhoOld = FunctionUtils::getRho( sigmaOld );
    const double sigmaRatio1 = sigmaNew / sigmaOld;
    const double sigmaRatio2 = ( rhoNew - rhoOld ) * sigmaNew;

    for ( unsigned int i = 0; i < input.size(); ++i ) {
        // all inputs other than capital
        double newCoef = 0;
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            double priceRatio = priceReceived / input[i]->getPricePaid( aRegionName, aPeriod );
            newCoef = pow( input[i]->getCoefficient( aPeriod ), sigmaRatio1 ) * pow( priceRatio, sigmaRatio2 );
        }
        // for capital or OVA input
        else {
            double priceRatio = priceReceived / (aProfits/capitalStock);
            newCoef = pow( input[i]->getCoefficient( aPeriod ), sigmaRatio1 ) * pow( priceRatio, sigmaRatio2 );
        }
        // set new transformed coefficient in input object
        input[i]->setCoefficient( newCoef, aPeriod );
    }

    // new alpha zero coefficient
    const double sigmaRatio = rhoNew / rhoOld * sigmaNew / sigmaOld; 
    return pow( alphaZero, sigmaRatio );
}

/*! \brief Calculate capital scaler for production technology, includes capital stock */
double CESProductionFunction::calcCapitalScaler( const InputSet& input, double aAlphaZero, double sigma, 
                                                 double capitalStock, const int aPeriod ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    const double rho = FunctionUtils::getRho( sigma );

    // find capital coefficient and multiply by stock to the power of rho
    IInput* capitalInput = FunctionUtils::getCapitalInput( input );
    assert( capitalInput );
    assert( capitalStock > 0 );

    double temp = capitalInput->getCoefficient( aPeriod ) * pow( capitalStock, rho );
    assert( util::isValidNumber( temp ) );
    double tempCapitalScaler = pow( temp, ( 1 / rho ) );

    /*! \post The capital scaler is a valid number. */
    assert( util::isValidNumber( tempCapitalScaler ) );
    return tempCapitalScaler;
}

/*! \brief Calculate capital rate scaler for expected profit rate calculation,
*          does not include actual capital stock
* \return The capital rate scaler.
*/
double CESProductionFunction::calcCapitalRateScaler( const InputSet& input, double sigma, const int aPeriod ) const {
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief production function. */
    assert( sigma >= 0.05 );
    const double rho = FunctionUtils::getRho( sigma );
    const IInput* capInput = FunctionUtils::getCapitalInput( input );
    assert( capInput );
    // do only capital to get capital contribution to production level
    assert( capInput->getCoefficient( aPeriod ) > 0 );
    return pow( capInput->getCoefficient( aPeriod ), ( 1 / rho ) );
}

/*! \brief Calculate profit scaler for production technology.
*/
double CESProductionFunction::calcFinalProfitScaler( const InputSet& input, const string& regionName, 
                                                     const string& sectorName, int period, 
                                                     double alphaZero, double sigma ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    const double rho = FunctionUtils::getRho( sigma );
    // calculate tempZ for all inputs except capital
    double tempZ = 0; // temporary scaler
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        // capital contribution calculated separately, see calcCapitalScaler()
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            double tempCoef = pow( input[i]->getCoefficient( period ), sigma );
            tempZ += tempCoef * pow( input[i]->getPricePaid( regionName, period ), (-rho * sigma) );
        }
    }
    // use price received for the good in the next equation
    double priceReceived = FunctionUtils::getPriceReceived( regionName, sectorName, period );
    double Z = 1 - tempZ * pow( ( priceReceived * alphaZero), (rho * sigma) );
    return ( Z > 0 ) ? pow( Z, -(1/rho) ) : 0;
}

//! Calculate Demand
double CESProductionFunction::calcDemand( InputSet& input, double personalIncome, 
                                          const string& regionName, const string& sectorName,
                                          const double aShutdownCoef, int period, double capitalStock,
                                          double alphaZero, 
                                          double sigma, double IBT) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    const double Z1 = calcCapitalScaler( input, alphaZero, sigma, capitalStock, period );
    const double Z2 = calcFinalProfitScaler( input, regionName, sectorName, period, alphaZero, sigma );
    const double Z = Z1 * Z2 * 
        pow( ( FunctionUtils::getPriceReceived( regionName, sectorName, period ) * alphaZero), sigma );
    
    
    double totalDemand = 0; // total demand used for scaling
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        // capital input name should be changed to OtherValueAdded
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            assert( input[i]->getPricePaid( regionName, period ) >= 0 );
            double pricePaid = max( input[i]->getPricePaid( regionName, period ), util::getSmallNumber() );
            double demand = pow( ( input[i]->getCoefficient( period ) / pricePaid ), sigma ) * aShutdownCoef * Z;
            if( demand < 0 ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::NOTICE );
                mainLog << "Demand less than zero for region " << regionName << " sector " << sectorName 
                        << " for input " << input[i]->getName() << " of " << demand << "." << endl;
            }
            input[i]->setCurrencyDemand( demand, regionName, period );
            totalDemand += demand;
        }
    }
    return totalDemand;
}

/*! \brief Calculate expected profit scaler for investment decision.
*/
double CESProductionFunction::calcExpProfitScaler( const InputSet& input, double aLifetimeYears, 
                                                   const string& regionName, const string& sectorName, 
                                                   int period, double alphaZero, double sigma ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );
    const double rho = FunctionUtils::getRho( sigma );
    
    // Calculate the net present value multiplier to determine expected prices.
    const double netPresentValueMult = FunctionUtils::getNetPresentValueMult( input, regionName, aLifetimeYears, period );

    // calculate tempZ for all inputs except capital
    double tempZ = 0; // temporary scaler
    for ( unsigned int i = 0; i < input.size(); ++i ) {
        // capital contribution calculated separately, see calcCapitalScaler()
        if( !input[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            tempZ += pow(input[i]->getCoefficient( period ), sigma) 
                     * pow( input[i]->getPricePaid( regionName, period ) * netPresentValueMult , ( -rho * sigma ) );
        }
    }
    // Calculate the expected price received.
    const double expPriceReceived = FunctionUtils::getExpectedPriceReceived( input, regionName,
                                                                             sectorName,
                                                                             aLifetimeYears, period );

    // use expected price received for the good in the next equation
    double Z = 1 - tempZ * pow( ( expPriceReceived * alphaZero ), ( rho * sigma ) );
    // okay to have negative Z, it means that technology is not profitable set Z
    // = 0 so that power does not blow up
    Z = max( Z, 0.0 );
    assert( util::isValidNumber( Z ) );
    return Z;
}

//! Calculate Expected Profit Rate
double CESProductionFunction::calcExpProfitRate( const InputSet& input, const string& regionName, 
                                                 const string& sectorName, double aLifetimeYears, int period,
                                                 double alphaZero, double sigma ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    const double rho = FunctionUtils::getRho( sigma );
    const double Z1 = calcCapitalRateScaler( input, sigma, period );
    const double Z2 = calcExpProfitScaler( input, aLifetimeYears, regionName, 
                                           sectorName, period, alphaZero, sigma );
    // Note: Dividing by the numeraire price at all here may not be correct.
    const IInput* numInput = FunctionUtils::getNumeraireInput( input );
    assert( numInput );
    const double pricePaidNumeraire = numInput->getPricePaid( regionName, period );
    assert( pricePaidNumeraire > 0 );
    // returns a rate, using price ratios
    const double expPriceReceived = FunctionUtils::getExpectedPriceReceived( input, regionName, sectorName,
                                                              aLifetimeYears, period );
    double expectedProfitRate = alphaZero * ( expPriceReceived / pricePaidNumeraire ) * 
        Z1 *  pow( Z2, (1/(-rho*sigma)) );

    assert( expectedProfitRate >= 0 );
    assert( util::isValidNumber( expectedProfitRate ) );
    return expectedProfitRate;
}

/*! \brief Calculate the levelized cost for this technology.
* \details Performs a levelized cost calculation, or cost per unit of output,
*          for the current period. This is performed by summing the coefficient,
*          or quantity used per unit of output, multiplied by the price paid for
*          that good. Capital is treated specially as the price paid is replaced
*          by the discount factor, which is equal to the interest rate
*          multiplied by the input's price adjustment. The rho and exp exponents
*          are used to allow substitution between inputs to a limited degree.
*          This function does not correct for the non-sensical levelized cost of
*          zero, it only avoids calculating invalid exponents. It is assumed the
*          caller will determine how to handle the zero. 
* 
* \author Josh Lurz
* \param aInputs Vector of inputs for the technology.
* \param aRegionName Name of the region containing the production function.
* \param aSectorName Name of the sector containing the production function.
* \param aPeriod Period in which to calculate levelized cost.
* \param aAlphaZero Out front scalar.
* \param aSigma Sigma coefficient.
* \warning This function does not correct for zero levelized cost.
* \return The levelized cost.
*/
double CESProductionFunction::calcLevelizedCost( const InputSet& aInputs,
                                                 const string& aRegionName,
                                                 const string& aSectorName,
                                                 int aPeriod,
                                                 double aAlphaZero,
                                                 double aSigma ) const
{
    // Calculate the exponents.
    const double rho = FunctionUtils::getRho( aSigma );
    const double exp = ( rho - 1 ) / rho;

    // Loop through the inputs and add on their costs.
    double levelizedCost = 0;
    for( unsigned int i = 0; i < aInputs.size(); ++i ){
        // Capital is done specially.
        if( aInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            // Store the capital price because getting it requires a call to the marketplace.
            const double capPrice = aInputs[ i ]->getPrice( aRegionName, aPeriod );
            // Check that the coefficient and adjusted price can both be raised
            // to exponents without overflowing or underflowing. Note: The price
            // adjustment is not included in price paid as the adjustment is not
            // correct. There should be a discount adjustment however.
            if( aInputs[ i ]->getCoefficient( aPeriod ) > util::getVerySmallNumber() &&
                ( capPrice > util::getVerySmallNumber() ) )
            {
                levelizedCost += pow( aInputs[ i ]->getCoefficient( aPeriod ), aSigma )
                    * pow( capPrice, exp );
                assert( util::isValidNumber( levelizedCost ) );
            }
        }
        // All other inputs.
        else {
            // Check that the coefficient and adjusted price can both be raised
            // to exponents without overflowing or underflowing.
            if( aInputs[ i ]->getCoefficient( aPeriod ) > util::getVerySmallNumber() &&
                ( aInputs[ i ]->getPricePaid( aRegionName, aPeriod ) > util::getVerySmallNumber() ) )
            {
                levelizedCost += pow( aInputs[ i ]->getCoefficient( aPeriod ), aSigma )
                    * pow( aInputs[ i ]->getPricePaid( aRegionName, aPeriod ), exp );
                assert( util::isValidNumber( levelizedCost ) );
            }
        }
    }
    /*! \post Levelized cost is greater than zero, as output must not be free. */
    assert( levelizedCost > 0 );
    // Make sure this calculation succeeds even if levelized cost is
    // non-sensical.
    return ( levelizedCost > 0 ) ? pow( levelizedCost, aSigma ) / aAlphaZero: 0;
}

//! Calculate profits
double CESProductionFunction::calcUnscaledProfits( const InputSet& input, const string& regionName,
                                                   const string& sectorName, const int period, 
                                                   const double capitalStock, const double alphaZero, 
                                                   const double sigma ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    double rho = FunctionUtils::getRho( sigma );
    double Z1 = calcCapitalScaler( input, alphaZero, sigma, capitalStock, period );
    double Z2 = calcFinalProfitScaler( input, regionName, sectorName, period, alphaZero, sigma );
    double priceReceived = FunctionUtils::getPriceReceived( regionName, sectorName, period );
    return alphaZero * priceReceived * Z1 * pow( Z2, (1 - rho) );
}

//! Calculate output
double CESProductionFunction::calcOutput( InputSet& input, const string& regionName, 
                                          const string& sectorName, const double aShutdownCoef,
                                          int period, double capitalStock, 
                                          double alphaZero, double sigma ) const 
{
    /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief
    *        production function. 
    */
    assert( sigma >= 0.05 );

    return alphaZero * 
           aShutdownCoef * 
           calcCapitalScaler( input, alphaZero, sigma, capitalStock, period ) * 
           calcFinalProfitScaler( input, regionName, sectorName, period, alphaZero, sigma );
}

//! Return the amount of output produced by one unit of capital.
double CESProductionFunction::getCapitalOutputRatio( const InputSet& aInputs, const string& aRegionName,
                                                     const string& aSectorName, double aLifetimeYears, int aPeriod,
                                                     double aAlphaZero, double aSigma ) const
{
    // Calculate the expected profit rate.
    const double profitRate = calcExpProfitRate( aInputs, aRegionName, aSectorName, aLifetimeYears, 
                                                 aPeriod, aAlphaZero, aSigma );
    
    // Check for negative profit rates.
    if( profitRate > 0 ){
        const double rho = FunctionUtils::getRho( aSigma );
        const double expA0 = rho / ( 1 - rho ); // var name is bad.
        const double capCoef = FunctionUtils::getCapitalInput( aInputs )->getCoefficient( aPeriod );
        const double expPriceReceived = FunctionUtils::getExpectedPriceReceived( aInputs, aRegionName, aSectorName,
                                                                                 aLifetimeYears, aPeriod );
        assert( capCoef > 0 );
        const double capQ = pow( aAlphaZero, expA0 ) * pow( capCoef, aSigma ) *
            pow( expPriceReceived / profitRate, aSigma );

        return capQ; // this is actually 1 / ratio, might be better to return
                     // the ratio.
    }
    // If there is no profit rate.
    return 0;
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
double CESProductionFunction::applyTechnicalChange( InputSet& input, const TechChange& aTechChange,
                                                    const string& regionName, const string& sectorName,
                                                    const int aPeriod, double alphaZero, double sigma ) const 
{
     /*! \pre sigma is greater than 0.05, otherwise we should be using a Leontief production function. */
    assert( sigma >= 0.05 );
    // Need to apply to a0 for hicks neutral, also material and energy.
    const double rho = FunctionUtils::getRho( sigma );
    for( unsigned int i = 0; i < input.size(); ++i ){
        double techChange = FunctionUtils::getTechChangeForInput( input[ i ],
                                                                  aTechChange,
                                                                  aPeriod );
        if( techChange > 0 ) {
            double scaleFactor = pow( 1 + techChange, scenario->getModeltime()->gettimestep( aPeriod ) );
            double newCoef = input[i]->getCoefficient( aPeriod ) * pow( scaleFactor, rho );
            input[i]->setCoefficient( newCoef, aPeriod );
        }
    }
    // Apply hicks tech change. Check to make sure this works correctly as it is untested.
    if( aTechChange.mHicksTechChange > 0 ){
        double scaleFactor = pow( 1 + aTechChange.mHicksTechChange,
                                  scenario->getModeltime()->gettimestep( aPeriod ) );
        alphaZero *= scaleFactor;
    }
    return alphaZero;
}
