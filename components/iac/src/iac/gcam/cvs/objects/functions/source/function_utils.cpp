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
* \file function_utils.cpp
* \ingroup Objects
* \brief The FunctionUtils class source file.
* \author Sonny Kim
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "functions/include/function_utils.h"
#include "functions/include/iinput.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/util.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "functions/include/ifunction.h" // for TechChange.
#include "containers/include/iinfo.h"
#include "util/logger/include/ilogger.h"
#include "functions/include/inested_input.h"
#include "functions/include/leaf_input_finder.h"

using namespace std;

extern Scenario* scenario; // for marketplace.

typedef InputSet::const_iterator CInputIterator;
typedef InputSet::iterator InputIterator;

//! Scale Input Coefficients
void FunctionUtils::scaleCoefficientInputs( InputSet& input,
                                            double scaler,
                                            const int aPeriod )
{
    for( unsigned int i = 0; i < input.size(); ++i ) {
        input[i]->setCoefficient( input[i]->getCoefficient( aPeriod ) * scaler, aPeriod );
    }
}

//! Function to return sum of all currency demand inputs.
double FunctionUtils::getCurrencyDemandSum( const InputSet& aInputs,
                                            const int aPeriod )
{
    double sum = 0;
    for( unsigned int i = 0; i < aInputs.size(); ++i ) {
        sum += aInputs[i]->getCurrencyDemand( aPeriod );
    }
    return sum;
}

//! Function to calculate the periodic payment amount of the loan with the first payment due a 
// full payment period into the loan.
// Future value is assumed to be 0.
// TODO: change aPV to Value class to incorporate units calculation.
double FunctionUtils::PMT( double aRate, double aNper, double aPV )
{
    // Formula from http://en.wikipedia.org/wiki/Amortization_calculator
    assert( aRate > 0 );
    assert( aNper > 0 );
    if( ( aNper > util::getVerySmallNumber() ) | 
        ( aRate > util::getVerySmallNumber() ) ) {
            return ( aPV * aRate ) / ( 1 - pow( ( 1 + aRate ), -aNper ) );
        }
    return -1;
}

//! Function to return average number of hours in a year.
double FunctionUtils::HOURS_PER_YEAR( void )
{
    return 8766;
}

//! Function to return hours per day.
double FunctionUtils::HOURS_PER_DAY( void )
{
    return 24;
}

//! Function to return 2003 to 1975 US $ deflator.
double FunctionUtils::DEFLATOR_1975_PER_DEFLATOR_2003( void )
{
    return 0.357; // Ratio of 1975/2003 deflators
}

//! Function to return 2005 to 1975 US $ deflator.
double FunctionUtils::DEFLATOR_1975_PER_DEFLATOR_2005( void )
{
    return 0.337; // Ratio of 1975/2005 deflators
}

//! Function to return 1990 to 1975 US $ deflator.
double FunctionUtils::DEFLATOR_1990_PER_DEFLATOR_1975( void )
{
    return 2.212; // Ratio of 1990/1975 deflators
}

//! Function to return kWh to GJ conversion.
double FunctionUtils::GJ_PER_KWH( void )
{
    return 0.003600457; // GJ/kWh
}

//! Function to return kWh to GJ conversion.
double FunctionUtils::GJ_PER_MWH( void )
{
    return 3.600457; // GJ/MWh
}

//! Function to return EJ to GJ conversion.
double FunctionUtils::GJ_PER_EJ( void )
{
    return 1.0E+9; // GJ/EJ
}

//! Function to return EJ to GJ conversion.
double FunctionUtils::MWH_PER_GWH( void )
{
    return 1.0E+3; // MWh/GWh
}

//! Function to return GWh to EJ conversion.
double FunctionUtils::EJ_PER_GWH( void )
{
    return GJ_PER_MWH() / GJ_PER_EJ() * MWH_PER_GWH(); // EJ/GWh
}

//! Function to return GJ to GWh conversion.
double FunctionUtils::GWH_PER_GJ( void )
{
    return 2.778E-4; // GWh/GJ
}

//! Function to return GJ to MWh conversion.
double FunctionUtils::MWH_PER_GJ( void )
{
    return 2.778E-1; // MWh/GJ
}

//! Function to return metric ton to kg conversion.
double FunctionUtils::KG_PER_METRIC_TON( void )
{
    return 1000; // kg/tonne
}

//! Function to return sum of all physical demand inputs
double FunctionUtils::getPhysicalDemandSum( const InputSet& aInputs,
                                            const int aPeriod )
{
    double sum = 0;
    for( unsigned int i = 0; i < aInputs.size(); ++i ) {
        sum += aInputs[i]->getPhysicalDemand( aPeriod );
    }
    return sum;
}

//! Return sum of Inputs
double FunctionUtils::getCoefSum( const InputSet& input,
                                  const int aPeriod )
{
    double sum = 0;
    for( unsigned int i = 0; i < input.size(); ++i ) {
        sum += input[i]->getCoefficient( aPeriod );
    }
    return sum;
}

/*! \brief Return a specific input given its name.
* \param aInputs The vector of inputs to search.
* \param aInputName The name of the input to find.
* \return A constant pointer to the input, null if it does not exist.
*/
IInput* FunctionUtils::getInput( const InputSet& aInputs,
                                 const string& aInputName )
{
    for( InputSet::const_iterator input = aInputs.begin(); input != aInputs.end(); ++input ) {
        if( ( *input )->getName() == aInputName ) {
            return *input;
        }
    }
    return 0;
}

//! Helper function to find the index of capital.
IInput* FunctionUtils::getCapitalInput( const InputSet& aInputs )
{
    for( unsigned int i = 0; i < aInputs.size(); ++i ) {
        if( aInputs[i]->hasTypeFlag( IInput::CAPITAL ) ) {
            return aInputs[i];
        }
    }
    return 0;
}

//! Helper function to get an input by type.
IInput* FunctionUtils::getInputByType( const InputSet& aInputs, const int aTypeFlag )
{
    for( unsigned int i = 0; i < aInputs.size(); ++i ) {
        if( aInputs[i]->hasTypeFlag( aTypeFlag ) ) {
            return aInputs[i];
        }
    }
    return 0;
}

/*! \brief Helper function to find the index of the numeraire input.
* \param aInputs The vector of inputs to the production function.
* \return A pointer to the numeraire input, null if there is not one(an error).
* \author Josh Lurz
*/
IInput* FunctionUtils::getNumeraireInput( const InputSet& aInputs )
{
    for( unsigned int i = 0; i < aInputs.size(); ++i ) {
        if( aInputs[i]->hasTypeFlag( IInput::NUMERAIRE ) ) {
            return aInputs[i];
        }
    }
    return 0;
}

//! Calculate Rho from Sigma and return
double FunctionUtils::getRho( const double aSigma )
{
    assert( aSigma > 0 );
    return ( aSigma - 1 ) / aSigma;
}

/*! \brief Calculate the net present value multiplier using the current price of
*          the capital input.
* \details 
* \param aInputs The vector of inputs to the production function.
* \param aRegionName Region name.
* \param aLifetimeYears The nameplate lifetime of the vintage.
* \param aPeriod Model period.
* \return The net present value multiplier.
* \author Josh Lurz
*/
double FunctionUtils::getNetPresentValueMult( const InputSet& aInputs,
                                              const string& aRegionName,
                                              const double aLifetimeYears,
                                              const int aPeriod )
{
    // Find the capital input.
    const IInput* capInput = getCapitalInput( aInputs );
    assert( capInput );
    double discountRate = capInput->getPricePaid( aRegionName, aPeriod ); // already includes adjustments

    // calculate net present value multiplier
    double netPresentValueMult = FunctionUtils::calcNetPresentValueMult( discountRate, aLifetimeYears );
    assert( util::isValidNumber( netPresentValueMult ) );
    assert( netPresentValueMult >= 0 );
    return netPresentValueMult;
}

/*! \brief Function that calculates the net present value based on discount rate
*          and lifetime.
* \details This function separates the net present value calculation from the
*          model as a utility function that can be called and used from
*          anywhere. More specifically this function is used to calculate the
*          expected future price received and price paid. The
*          netPresentValueMult is multiplied by the price paid and received.
* \param aDiscountRate The discount rate
* \param aLifetime The lifetime over which to calculate the multiplier.
* \return The net present value multiplier.
* \author Sonny Kim
*/
double FunctionUtils::calcNetPresentValueMult( const double aDiscountRate,
                                               const double aLifetime )
{
    double netPresentValueMult = pow( ( 1 / ( 1 + aDiscountRate ) ), 0.5 )
                                * ( 1 + ( 1 / aDiscountRate ) )
                                * ( 1 - pow( ( 1 / ( 1 + aDiscountRate ) ), aLifetime ) );
    return netPresentValueMult;
}

/*! \brief Set the price paid for a good into the marketplace.
* \details Helper function which gets the information object from the market and
*          sets the price paid.
* \param aRegionName Name of the region for which to set the price paid.
* \param aGoodName Name of the good for which to set the price paid.
* \param aPeriod Model period in which to set the price paid.
* \param aPricePaid The price paid.
*/
void FunctionUtils::setPricePaid( const string& aRegionName,
                                  const string& aGoodName,
                                  const int aPeriod,
                                  const double aPricePaid )
{
    assert( aGoodName != "USA" );
    IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName, aPeriod, true );

    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    marketInfo->setDouble( "pricePaid", aPricePaid );
}

/*! \brief Gets the price paid for the good by querying the marketplace.
* \details Helper function which gets the information object from the market,
*          and queries it for the price paid.
* \param aRegionName Name of the region for which to get the price paid.
* \param aGoodName Name of the good for which to get the price paid.
* \param aPeriod Model period in which to get the price paid.
* \return The price paid, 0 if there is not one set.
*/
double FunctionUtils::getPricePaid( const string& aRegionName,
                                    const string& aGoodName,
                                    const int aPeriod )
{
    assert( aGoodName != "USA" );
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* marketInfo = marketplace->getMarketInfo( aGoodName, aRegionName, aPeriod, true );

    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    return marketInfo->getDouble( "pricePaid", true );
}

/*! \brief Set the price received for a good into the marketplace.
* \details Helper function which gets the information object from the market and
*          sets the price received.
* \param aRegionName Name of the region for which to set the price received.
* \param aGoodName Name of the good for which to set the price received.
* \param aPeriod Model period in which to set the price received.
* \param aPriceReceived Price received to set.
*/
void FunctionUtils::setPriceReceived( const string& aRegionName,
                                      const string& aGoodName,
                                      const int aPeriod,
                                      const double aPriceReceived )
{
    assert( aGoodName != "USA" );
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( aGoodName, aRegionName, aPeriod, true );

    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    marketInfo->setDouble( "priceReceived", aPriceReceived );
}

/*! \brief Gets the price received for the good by querying the marketplace.
* \details Helper function which gets the information object from the market,
*          and queries it for the price received.
* \param aRegionName Name of the region for which to get the price received.
* \param aGoodName Name of the good for which to get the price received.
* \param aPeriod Model period in which to get the price received.
* \return The price received, 0 if there is not one set.
*/
double FunctionUtils::getPriceReceived( const string& aRegionName,
                                        const string& aGoodName,
                                        const int aPeriod )
{
    assert( aGoodName != "USA" );
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* marketInfo = marketplace->getMarketInfo( aGoodName, aRegionName, aPeriod, true );
    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    return marketInfo->getDouble( "priceReceived", true );
}

/*! \brief Calculate the expected price received for the good produced by the
*          production function.
* \details Calculate the expected price received for a good over the lifetime.
*          This is done by performing a net present value calculation using the
*          interest rate as determined from the capital input. 
* \param aInputs The vector of inputs to the production function.
* \param aRegionName The name of the region containing the production function.
* \param aGoodName The name of the sector containing the production function.
* \param aLifetimeYears The nameplate lifetime of vintage.
* \param aPeriod Period
* \return The expected price received for the good.
* \author Josh Lurz
*/
double FunctionUtils::getExpectedPriceReceived( const InputSet& aInputs,
                                                const string& aRegionName,
                                                const string& aGoodName,
                                                const double aLifetimeYears,
                                                const int aPeriod )
{
    // calculate expected price received for the produced good
    return getPriceReceived( aRegionName, aGoodName, aPeriod )
           * getNetPresentValueMult( aInputs, aRegionName, aLifetimeYears, aPeriod );
}

/*! \brief Apply technical change to a vector of inputs.
* \details 
* \param input Vector of inputs for the demand function.
* \param aTechChange A structure containing the various possible types of
*        technical change.
* \param regionName Name of the region containing the function.
* \param sectorName Name of the sector containing the function.
* \param alphaZero The up-front scaler.
* \param sigma Sigma coefficient.
* \return The new alpha zero.
*/
double FunctionUtils::applyTechnicalChangeInternal( InputSet& aInputs,
                                                    const TechChange& aTechChange,
                                                    const string& aRegionName,
                                                    const string& aSectorName,
                                                    const int aPeriod,
                                                    double aAlphaZero,
                                                    double aSigma )
{
    const int timeStep = scenario->getModeltime()->gettimestep( aPeriod );
    for( InputIterator i = aInputs.begin(); i != aInputs.end(); ++i ) {
        double techChange = FunctionUtils::getTechChangeForInput( *i,
                                                                  aTechChange,
                                                                  aPeriod );
        if( techChange > 0 ) {
            double scaleFactor = pow( 1 + techChange, timeStep );
            double newCoef = (*i)->getCoefficient( aPeriod ) / scaleFactor;
            (*i)->setCoefficient( newCoef, aPeriod );
        }
    }
    // Apply hicks tech change.
    if( aTechChange.mHicksTechChange > 0 ) {
        double scaleFactor = pow( 1 + aTechChange.mHicksTechChange,
                                  timeStep );
        aAlphaZero *= scaleFactor;
    }

    return aAlphaZero;
}

/*!
 * \brief Determine the technical change coefficient to apply to the input.
 * \details The technical change coefficient applied to an input may be one of
 *          three items. If the input has a read-in technical change
 *          coefficient(currently greater than zero), then that technical change
 *          will be used. Otherwise, the Technology's energy technical change
 *          will be used if it is an energy input, and the Technology's material
 *          technical change input will be if it is a material input. An input
 *          cannot both be a material and energy input. Note that Hicks neutral
 *          technical change will always be applied in addition.
 * \param aInput Input for which to get the technical change.
 * \param aTechChange Technical change containing structure.
 * \param aPeriod Model period.
 * \return The technical change to apply to the input.
 */
double FunctionUtils::getTechChangeForInput( const IInput* aInput,
                                             const TechChange& aTechChange,
                                             const int aPeriod )
{
    // Get the input specific technical change.
    double techChange = aInput->getTechChange( aPeriod );

    // If an input specific technical change was not read-in, use an energy
    // or material technical change. TODO: Check for not-initialized vs.
    // zero.
    if( techChange == 0 ){
        if( aInput->hasTypeFlag( IInput::ENERGY ) ){
            techChange = aTechChange.mEnergyTechChange;
        }
        else if( aInput->hasTypeFlag( IInput::MATERIAL ) ){
            techChange = aTechChange.mMaterialTechChange;
        }
    }
    return techChange;
}

/*! \brief Return whether a good is a fixed price good.
* \param aRegionName Region name.
* \param aGoodName Good name.
* \param aPeriod Model period.
* \return Whether the good is a fixed price good.
*/
bool FunctionUtils::isFixedPrice( const string& aRegionName,
                                  const string& aGoodName,
                                  const int aPeriod )
{
    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName, aPeriod, false );
    return marketInfo && marketInfo->getBoolean( "IsFixedPrice", false );
}

/*! \brief Static function which returns the conversion factor for the good
*          which is stored in the marketplace.
* \param aRegionName The region name.
* \param aGoodName The good name.
* \param aMustExist If the conversion factor must exist.  It would be alright
*                   for it to not exist for consumers.
* \return The conversion factor for the good.
*/
double FunctionUtils::getMarketConversionFactor( const string& aRegionName,
                                                 const string& aGoodName,
                                                 const bool aMustExist )
{
    assert( !aGoodName.empty() && !aRegionName.empty() );
    assert( aGoodName != "USA" );

    const IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( aGoodName, aRegionName, 0, aMustExist );

    return marketInfo ? marketInfo->getDouble( "ConversionFactor", aMustExist ) : 0;
}

/*!
 * \brief Copy information about input paramaters forward to the next period.
 * \param aPrevInputs Previous period's set of inputs.
 * \param aCurrInputs Current period's set of inputs.
 * \param aPeriod Period in which the copying is occurring.
 */
void FunctionUtils::copyInputParamsForward( const InputSet& aPrevInputs,
                                            InputSet& aCurrInputs,
                                            const int aPeriod )
{
    for( CInputIterator prev = aPrevInputs.begin(); prev != aPrevInputs.end(); ++prev ){
        // Find the matching input in the current set.
        IInput* currInput = FunctionUtils::getInput( aCurrInputs, (*prev)->getName() );
        if( !currInput ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Input " << (*prev)->getName()
                    << " does not have a matching input in the next period." << endl;
            continue;
        }
        currInput->copyParam( *prev, aPeriod );
    }
}

/*!
 * \brief Static function which returns the CO2 content for a good
 *        if it were to be emitted.  The coefficient is stored in the marketplace.
 * \param aRegionName The region name.
 * \param aGoodName The good name.
 * \param aPeriod The period for which to look.
 * \param aMustExist If the CO2 coefficient must exist.  It would be alright
 *                   for it to not exist for consumers.
 * \return The CO2 emissions coefficient for the good.
 */
double FunctionUtils::getCO2Coef( const string& aRegionName,
                                  const string& aGoodName,
                                  const int aPeriod,
                                  const bool aMustExist )
{
    // Initialize the cached CO2 coefficient.
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* productInfo = marketplace->getMarketInfo( aGoodName, aRegionName, aPeriod, aMustExist );

    // Output ratio is determined by the CO2 coefficient and the ratio of output
    // to the primary good. The info should not be null except in cases of
    // improperly constructed input files. This function will have already
    // warned in that case.
    return productInfo ? productInfo->getDouble( "CO2coefficient", false ) : 0;
}

/*!
 * \brief Get the ratio of an input's price in the current period to the base
 *        period.
 * \details Returns the ratio of the input's price in the given region in the
 *          current period to the base period. If the input's price in the base
 *          period is zero or if the base period is greater than or equal to the
 *          current period this will return 1.
 * \param aRegionName Name of the region in which to find the ratio.
 * \param aInput Input for which to calculate the price ratio.
 * \param aBasePeriod Base period for the ratio.
 * \param aCurrentPeriod Current period for the ratio.
 * \return Price ratio for the input.
 * \author Josh Lurz
 */
double FunctionUtils::calcPriceRatio( const string& aRegionName,
                                      const IInput* aInput,
                                      const int aBasePeriod,
                                      const int aCurrentPeriod )
{
    // The price ratio is always 1 in the base period.
    double priceRatio = 1;
    if( aCurrentPeriod > aBasePeriod ) {
        double basePrice = aInput->getPrice( aRegionName, aBasePeriod );
        if( basePrice > util::getSmallNumber() ){
            priceRatio = aInput->getPrice( aRegionName, aCurrentPeriod ) / basePrice;
        }
    }
    return priceRatio;
}

/*!
 * \brief Gets the leaf inputs from a nested input structure.
 * \details This method will use the LeafInputFinder visitor to create a 
 *          list of any leaf inputs that are a descendant (or self) of the passed 
 *          in INestedInput.  See LeafInputFinder to determine which inputs
 *          are considered leaves.
 * \param aNestedInput The nested input to search for leaves.
 * \return An InputSet of the leaves that are withing the passed in nested
 *         input.
 * \author Pralit Patel
 */
InputSet FunctionUtils::getLeafInputs( const INestedInput* aNestedInput ) {
    assert( aNestedInput );
    LeafInputFinder inputFinder;
    aNestedInput->accept( &inputFinder, -1 );
    return inputFinder.getInputs();
}

/*!
 * \breif Sets the price of the capital good.
 * \details The price of the capital good is the CES aggregate of
 *          goods that make up capital usually determined by inputs
 *          to the investment consumer.
 * \param aRegionName The name of the region in which to set the price.
 * \param aPeriod The period in which we want to set the price.
 * \param aCapitalGoodPrice The price to set.
 * \author Pralit Patel
 */
void FunctionUtils::setCapitalGoodPrice( const string& aRegionName,
                                         const int aPeriod,
                                         const double aCapitalGoodPrice )
{
    static const string capitalGoodName = "Capital";
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( capitalGoodName, aRegionName, aPeriod, true );

    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    marketInfo->setDouble( "CapitalGoodPrice", aCapitalGoodPrice );
}

/*!
 * \breif Gets the price of the capital good.
 * \details The price of the capital good is the CES aggregate of
 *          goods that make up capital usually determined by inputs
 *          to the investment consumer.  This price is used to convert
 *          a quantity amount of capital to a dollar amount which can
 *          be used for investment.
 * \param aRegionName The name of the region in which to get the price.
 * \param aPeriod The period in which we want to get the price.
 * \return The price of the capital good.
 * \author Pralit Patel
 */
double FunctionUtils::getCapitalGoodPrice( const string& aRegionName,
                                           const int aPeriod )
{
    static const string capitalGoodName = "Capital";
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* marketInfo = marketplace->getMarketInfo( capitalGoodName, aRegionName, aPeriod, true );
    /*! \invariant The market and market info must exist. */
    assert( marketInfo );
    return marketInfo->getDouble( "CapitalGoodPrice", true );
}
