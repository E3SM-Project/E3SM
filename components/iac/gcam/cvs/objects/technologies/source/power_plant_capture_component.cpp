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
 * \file power_plant_capture_component.cpp
 * \ingroup Objects
 * \brief PowerPlantCaptureComponent source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "technologies/include/power_plant_capture_component.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/iinfo.h"
#include "containers/include/scenario.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/dependency_finder.h"
#include "functions/include/iinput.h"
#include "containers/include/scenario.h"

using namespace std;

extern Scenario* scenario;

/*!
 * \brief Constructor.
 * \details Protected constructor which prevents the capture component from
 *          being created without using the CaptureComponentFactory.
 */
PowerPlantCaptureComponent::PowerPlantCaptureComponent():
mSequesteredAmount( scenario->getModeltime()->getmaxper() ),
mRemoveFraction( 0 ),
mCaptureEnergy( 0 ),
mNonEnergyCostPenalty( 0 )
{
}

PowerPlantCaptureComponent* PowerPlantCaptureComponent::clone() const {
    return new PowerPlantCaptureComponent( *this );
}

bool PowerPlantCaptureComponent::isSameType( const std::string& aType ) const {
    return aType == getXMLNameStatic();
}

/*!
 * \brief Get the XML node name in static form for comparison when parsing XML.
 * \details This public function accesses the private constant string, XML_NAME.
 *          This way the tag is always consistent for both read-in and output
 *          and can be easily changed. The "==" operator that is used when
 *          parsing, required this second function to return static.
 * \note A function cannot be static and virtual.
 * \author Josh Lurz, James Blackwood
 * \return The constant XML_NAME as a static.
 */
const string& PowerPlantCaptureComponent::getXMLNameStatic() {
    const static string XML_NAME = "power-plant-capture-component";
    return XML_NAME;
}

const string& PowerPlantCaptureComponent::getName() const {
    return getXMLNameStatic();
}

// Documentation inherits.
bool PowerPlantCaptureComponent::XMLParse( const xercesc::DOMNode* node ){
    /*! \pre Assume we are passed a valid node. */
    assert( node );

    const xercesc::DOMNodeList* nodeList = node->getChildNodes();
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
        const xercesc::DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() != xercesc::DOMNode::ELEMENT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "storage-market" ){
            mStorageMarket = XMLHelper<string>::getValue( curr );
        }
        else if( nodeName == "remove-fraction" ){
            mRemoveFraction = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "capture-energy" ){
            mCaptureEnergy = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "non-energy-penalty" ){
            mNonEnergyCostPenalty = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "target-gas" ){
            mTargetGas = XMLHelper<string>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Unknown tag " << nodeName << " encountered while processing "
                    << getXMLNameStatic() << endl;
        }
    }
    // TODO: Handle success and failure better.
    return true;
}

void PowerPlantCaptureComponent::toInputXML( ostream& aOut,
                                             Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElementCheckDefault( mStorageMarket, "storage-market", aOut, aTabs, string( "" ) );
    XMLWriteElementCheckDefault( mRemoveFraction, "remove-fraction", aOut, aTabs, 0.0 );
    XMLWriteElementCheckDefault( mCaptureEnergy, "capture-energy", aOut, aTabs, 0.0 );
    XMLWriteElementCheckDefault( mNonEnergyCostPenalty, "non-energy-penalty", aOut, aTabs, 0.0 );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void PowerPlantCaptureComponent::toDebugXML( const int aPeriod,
                                             ostream& aOut,
                                             Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mStorageMarket, "storage-market", aOut, aTabs );
    XMLWriteElement( mRemoveFraction, "remove-fraction", aOut, aTabs );
    XMLWriteElement( mCaptureEnergy, "capture-energy", aOut, aTabs );
    XMLWriteElement( mNonEnergyCostPenalty, "non-energy-penalty", aOut, aTabs );
    XMLWriteElement( mSequesteredAmount[ aPeriod ], "sequestered-amount", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void PowerPlantCaptureComponent::completeInit( const string& aRegionName,
                                               const string& aSectorName,
                                               DependencyFinder* aDependencyFinder )
{
    // Add the storage market as a dependency of the sector. This is because
    // this sector will have to be ordered first so that the total demand and
    // price for storage are known.
    aDependencyFinder->addDependency( aSectorName, mStorageMarket );
    
    // Check that the remove fraction is valid.
    if( mRemoveFraction < 0 || mRemoveFraction > 1 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid removal fraction of " << mRemoveFraction << "." << endl;
        mRemoveFraction = 0;
    }
    // Default the target gas to CO2.
    if( mTargetGas.empty() ){
        mTargetGas = "CO2";
    }
}

void PowerPlantCaptureComponent::initCalc( const string& aRegionName,
                                           const string& aSectorName,
                                           const string& aFuelName,
                                           const int aPeriod )
{
    // Calculate the emissions coefficient of the fuel.
    const IInfo* fuelInfo = scenario->getMarketplace()->getMarketInfo( aFuelName, aRegionName, aPeriod, true );
    mCachedFuelCoef = fuelInfo ? fuelInfo->getDouble( "CO2Coef", true ) : 0;
}

/**
 * \details Storage cost is only valid for mTargetGas ("CO2")
 * \param aRegionName 
 * \param aGHGName 
 * \param aPeriod 
 * \return storage cost
 */
double PowerPlantCaptureComponent::getStorageCost( const string& aRegionName,
                                                   const string& aGHGName,
                                                   const int aPeriod ) const
{
    if( aGHGName != mTargetGas ){
        return 0;
    }

    // Check if there is a market for storage.
    double storageMarketPrice = scenario->getMarketplace()->getPrice( mStorageMarket,
                                                                      aRegionName,
                                                                      aPeriod, true );
    // Check if there is a carbon market.
    double carbonMarketPrice = scenario->getMarketplace()->getPrice( mTargetGas,
                                                                     aRegionName,
                                                                     aPeriod, false );

    // Retrieve proportional tax rate.
    const Marketplace* marketplace = scenario->getMarketplace();
    const IInfo* marketInfo = marketplace->getMarketInfo( "CO2", aRegionName, aPeriod, false );
    // Note: the key includes the region name.
    const double proportionalTaxRate = 
        ( marketInfo && marketInfo->hasValue( "proportional-tax-rate" + aRegionName ) ) 
        ? marketInfo->getDouble( "proportional-tax-rate" + aRegionName, true )
        : 1.0;
    // Adjust greenhouse gas tax with the proportional tax rate.
    carbonMarketPrice *= proportionalTaxRate;

    // If there is no carbon market, return a large number to disable the
    // capture technology.
    if( carbonMarketPrice == Marketplace::NO_MARKET_PRICE || carbonMarketPrice < util::getSmallNumber() ){
        return util::getLargeNumber();
    }
    // If carbon and storage markets exists use the storage market price.
    return ( storageMarketPrice == Marketplace::NO_MARKET_PRICE ) ? util::getLargeNumber()
                                                                  : storageMarketPrice;
}

/**
 * \details Has a valid remove fraction for the target gas only (currently CO2).
 * \param aGHGName 
 * \return remove fraction
 */
double PowerPlantCaptureComponent::getRemoveFraction( const string& aGHGName ) const {
    return aGHGName == mTargetGas ? mRemoveFraction : 0;
}

/**
 * \details Calculate sequestered amount for all gases, but do not add
 *  to market demand if gas is not mTargetGas.
 * \param aRegionName 
 * \param aGHGName 
 * \param aTotalEmissions 
 * \param aPeriod 
 * \return emissions sequestered
 */
double PowerPlantCaptureComponent::calcSequesteredAmount( const string& aRegionName,
                                                          const string& aGHGName,
                                                          const double aTotalEmissions,
                                                          const int aPeriod )
{
    // Calculate the amount.
    mSequesteredAmount[ aPeriod ] = mRemoveFraction * aTotalEmissions;
    
    // Add the demand to the marketplace.
    if( mSequesteredAmount[ aPeriod ] > 0 ){
        // set sequestered amount as demand side of carbon storage market
        Marketplace* marketplace = scenario->getMarketplace();
        if( aGHGName == mTargetGas ){
            marketplace->addToDemand( mStorageMarket, aRegionName, mSequesteredAmount[ aPeriod ],
                aPeriod, false );
        }
    }
    return mSequesteredAmount[ aPeriod ];
}

/**
 * \param aGHGName 
 * \param aGetGeologic 
 * \param aPeriod 
 * \return sequestered amount
 */
double PowerPlantCaptureComponent::getSequesteredAmount( const string& aGHGName,
                                                         const bool aGetGeologic,
                                                         const int aPeriod ) const 
{
    // Only return emissions if the type of the sequestration equals is geologic.
    if( aGetGeologic ){
        return mSequesteredAmount[ aPeriod ];
    }
    return 0;
}

/*! \brief Adjust the set of inputs for a technology for the costs and
*          efficiency losses due to the capture component.
* \param aRegionName Name of the region.
* \param aInputs Vector of technology inputs.
* \param aPeriod Model period.
*/
void PowerPlantCaptureComponent::adjustInputs( const string& aRegionName,
                                               vector<IInput*>& aInputs,
                                               const int aPeriod ) const
{
    // TODO: Improve this code!
    double baseEnergyIntensity = 0;
    double effectiveEnergyIntensity = 0;
    double fuelEmissCoef = 0;

    // Loop through the inputs and search for energy inputs.
    for( unsigned int i = 0; i < aInputs.size(); ++i ){
        if( aInputs[ i ]->hasTypeFlag( IInput::ENERGY ) ){
            // Store the unadjusted coefficient.
            // TODO: Handle multiple energy inputs.
            baseEnergyIntensity = aInputs[ i ]->getCoefficient( aPeriod );
            // TODO: Unhardcode has name.
            fuelEmissCoef = aInputs[ i ]->getCO2EmissionsCoefficient( mTargetGas, aPeriod );
            adjustEnergyInput( aInputs[ i ], aPeriod );
            effectiveEnergyIntensity = aInputs[ i ]->getCoefficient( aPeriod );
        }
    }

    // Now adjust the non-energy input.
    // TODO: What if an energy input wasn't found?
    for( unsigned int i = 0; i < aInputs.size(); ++i ){
        if( aInputs[ i ]->hasTypeFlag( IInput::CAPITAL ) ){
            adjustNonEnergyInput( aInputs[ i ], aRegionName, baseEnergyIntensity,
                                  effectiveEnergyIntensity, fuelEmissCoef, aPeriod );
        }
    }
}

void PowerPlantCaptureComponent::adjustEnergyInput( IInput* aEnergyInput,
                                                    const int aPeriod ) const
{
    assert( aEnergyInput > 0 );

    // Calculate effective intensity: This increases the intensity by first converting
	// to an efficiency then subtracting by the product of capture energy, CO2 coefficient
	// and removal fraction, and finally converting back to an intensity.
 
	double adjustedIntensity = 1/(1/aEnergyInput->getCoefficient( aPeriod ) 
	                          - mCaptureEnergy
	                          * aEnergyInput->getCO2EmissionsCoefficient( mTargetGas, aPeriod )
	                          * mRemoveFraction);

    aEnergyInput->setCoefficient( adjustedIntensity, aPeriod );
}

void PowerPlantCaptureComponent::adjustNonEnergyInput( IInput* aNonEnergyInput,
                                                       const string& aRegionName,
                                                       const double aBaseEnergyIntensity,
                                                       const double aEffectiveEnergyIntensity,
                                                       const double aFuelEmissCoef,
                                                       const int aPeriod ) const
{
    assert( aBaseEnergyIntensity >= 0 );
    assert( aEffectiveEnergyIntensity >= 0 );
    assert( aFuelEmissCoef >= 0 );
    assert( aEffectiveEnergyIntensity >= aBaseEnergyIntensity );

    // Calculate the "a" term.
    const double a =  aBaseEnergyIntensity * aFuelEmissCoef * mRemoveFraction;
    
    // A must be positive.
    assert( a >= 0 );   
    
    // Calculate the total non-energy cost.
    const double totalNonEnergyCost = ( aNonEnergyInput->getPrice( aRegionName, aPeriod )
                                      + a * mNonEnergyCostPenalty ) * aEffectiveEnergyIntensity
                                      / aBaseEnergyIntensity;

    // Total non-energy cost is greater or equal to zero.
    assert( totalNonEnergyCost >= 0 );
    aNonEnergyInput->setPrice( aRegionName, totalNonEnergyCost, aPeriod );
}
