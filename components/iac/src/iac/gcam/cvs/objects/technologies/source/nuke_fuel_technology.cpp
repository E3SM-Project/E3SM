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
* \file nuke_fuel_technology.cpp
* \ingroup Objects
* \brief Nuclear fuel technology class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <cmath>
#include "technologies/include/nuke_fuel_technology.h"
#include "emissions/include/by_product.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "util/logger/include/ilogger.h"
#include "technologies/include/iproduction_state.h"
#include "technologies/include/marginal_profit_calculator.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

typedef vector<IInput*>::iterator InputIterator;
typedef vector<IOutput*>::iterator OutputIterator;

//! Constructor.
NukeFuelTechnology::NukeFuelTechnology( const string& aName, const int aYear ): Technology( aName, aYear ) {
    // default values for nuclear fuel
    fertileFuelName = "none"; // name of secondary fertile material used for making nuclear fuel
    blanketFuelName = "none"; // name of secondary fertile material used for breeding fissile material
    blanketFuelRatio = 0.0; // Ratio of blanket to fuel materials (kgBlanket/kgFuel)
    burnup = 50; // designed burnup of fuel associated with nuclear plant (MWd/kgHM)
    conversionCost = 5; // uranium ore conversion cost ($/kgU)
    enrichmentProd = 0.045; // fissile material enrichment 
    enrichmentFeed = 0.0071; // feed material enrichment 
    enrichmentTail = 0.003; // tails enrichment 
    enrichmentCost = 100; // uranium enrichment cost ($/SWU)
    fabricationCost = 200; // enriched uranium fuel fabrication cost ($/kgHM)
    blanketFabCost = 0.0; // blanket material fabrication cost ($/kgHM)
    interimStorageCost = 200; // interim storage cost of spent fuel ($/kgHM)
    geologicWasteDisposalCost = 400; // cost of permenant waste disposal ($/kgHM)
    reprocessingCost = 0; // reprocessing cost of spent fuel ($/kgHM)
    mConversionFactor = 1;
}

NukeFuelTechnology* NukeFuelTechnology::clone() const {
    return new NukeFuelTechnology( *this );
}

const string& NukeFuelTechnology::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& NukeFuelTechnology::getXMLNameStatic() {
    const static string XML_NAME = "nuclearFuelTechnology";
    return XML_NAME;
}

bool NukeFuelTechnology::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if( nodeName == "burnup" ){
        burnup = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "fertileFuelName" ){
        fertileFuelName = XMLHelper<string>::getValue( curr );
    }
    else if( nodeName == "blanketFuelName" ){
        blanketFuelName = XMLHelper<string>::getValue( curr );
    }
    else if( nodeName == "blanketFuelRatio" ){
        blanketFuelRatio = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "conversionCost" ){
        conversionCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "enrichmentCost" ){
        enrichmentCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "enrichmentProd" ){
        enrichmentProd = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "enrichmentFeed" ){
        enrichmentFeed = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "enrichmentTail" ){
        enrichmentTail = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "fabricationCost" ){
        fabricationCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "blanketFabCost" ){
        blanketFabCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "interimStorageCost" ){
        interimStorageCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "geologicWasteDisposalCost" ){
        geologicWasteDisposalCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "reprocessingCost" ){
        reprocessingCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "fMultiplier" ){
        mConversionFactor = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

void NukeFuelTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {  
    XMLWriteElementCheckDefault( fertileFuelName, "fertileFuelName", out, tabs, string( "none" ) );
    XMLWriteElementCheckDefault( blanketFuelName, "blanketFuelName", out, tabs, string( "none" ) );
    XMLWriteElementCheckDefault( burnup, "burnup", out, tabs, 50.0 );
    XMLWriteElementCheckDefault( blanketFuelRatio, "blanketFuelRatio", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( conversionCost, "conversionCost", out, tabs, 5.0 );
    XMLWriteElementCheckDefault( enrichmentProd, "enrichmentProd", out, tabs, 0.045 );
    XMLWriteElementCheckDefault( enrichmentFeed, "enrichmentFeed", out, tabs, 0.0071 );
    XMLWriteElementCheckDefault( enrichmentTail, "enrichmentTail", out, tabs, 0.003 );
    XMLWriteElementCheckDefault( enrichmentCost, "enrichmentCost", out, tabs, 100.0 );
    XMLWriteElementCheckDefault( fabricationCost, "fabricationCost", out, tabs, 200.0 );
    XMLWriteElementCheckDefault( blanketFabCost, "blanketFabCost", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( interimStorageCost, "interimStorageCost", out, tabs, 200.0 );
    XMLWriteElementCheckDefault( geologicWasteDisposalCost, "geologicWasteDisposalCost", out, tabs, 400.0 );
    XMLWriteElementCheckDefault( reprocessingCost, "reprocessingCost", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mConversionFactor, "fMultiplier", out, tabs, 1.0 );
}	

void NukeFuelTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const { 
    XMLWriteElement( fertileFuelName, "fertileFuelName", out, tabs );
    XMLWriteElement( blanketFuelName, "blanketFuelName", out, tabs );
    XMLWriteElement( burnup, "burnup", out, tabs );
    XMLWriteElement( blanketFuelRatio, "blanketFuelRatio", out, tabs );
    XMLWriteElement( conversionCost, "conversionCost", out, tabs );
    XMLWriteElement( enrichmentProd, "enrichmentProd", out, tabs );
    XMLWriteElement( enrichmentFeed, "enrichmentFeed", out, tabs );
    XMLWriteElement( enrichmentTail, "enrichmentTail", out, tabs );
    XMLWriteElement( enrichmentCost, "enrichmentCost", out, tabs );
    XMLWriteElement( fabricationCost, "fabricationCost", out, tabs );
    XMLWriteElement( blanketFabCost, "blanketFabCost", out, tabs );
    XMLWriteElement( interimStorageCost, "interimStorageCost", out, tabs );
    XMLWriteElement( geologicWasteDisposalCost, "geologicWasteDisposalCost", out, tabs );
    XMLWriteElement( reprocessingCost, "reprocessingCost", out, tabs );
    XMLWriteElement( mConversionFactor, "fMultiplier", out, tabs);
}	

void NukeFuelTechnology::initCalc( const string& aRegionName,
                                  const string& aSectorName,
                                  const IInfo* aSubsectorInfo,
                                  const Demographic* aDemographics,
                                  PreviousPeriodInfo& aPrevPeriodInfo,
                                  const int aPeriod )
{
    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
        aDemographics, aPrevPeriodInfo, aPeriod );

    // TODO: This is a hack. Adjust the coeffient for the energy input
    // to be the special effiency.
    // if representing pure fissile material, getFeedProductRatio() should return
    // the fissile enrichment percentage
    double kgFissilePerGJ = getInitialMass() * getFeedProductRatio(); //fissile or ore mass

    // Find the energy input. There should only be one.
    for( InputIterator i = mInputs.begin(); i != mInputs.end(); ++i ){
        (*i)->setCoefficient( kgFissilePerGJ, aPeriod );
    }
}
void NukeFuelTechnology::completeInit( const string& aRegionName,
                                      const string& aSectorName,
                                      const string& aSubsectorName,
                                      DependencyFinder* aDepFinder,
                                      const IInfo* aSubsectorInfo,
                                      ILandAllocator* aLandAllocator )
{
    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder,
        aSubsectorInfo, aLandAllocator );

    for( OutputIterator i = mOutputs.begin(); i != mOutputs.end(); ++i ){
        ( *i )->scaleCoefficient( getInitialMass() );
    }

    if( enrichmentFeed < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Feed enrichment must be greater than zero." << endl;
    }

    if( enrichmentTail < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Tail enrichment must be greater than zero." << endl;
    }

    if( enrichmentProd < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Product enrichment must be greater than zero." << endl;
    }

    if( enrichmentFeed - enrichmentTail < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Feed enrichment must be greater than tail enrichment." << endl;
    }
    if( enrichmentProd - enrichmentTail < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Product enrichment must be greater than tail enrichment." << endl;
    }

    if( burnup < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Burnup ratio must be positive." << endl;
    }
}

void NukeFuelTechnology::postCalc( const string& aRegionName, const int aPeriod ) {
    Technology::postCalc( aRegionName, aPeriod );
}

void NukeFuelTechnology::calcCost( const string& aRegionName,
                                  const string& aSectorName,
                                  const int aPeriod )
{
    Technology::calcCost( aRegionName, aSectorName, aPeriod );
}

double NukeFuelTechnology::getTotalInputCost( const string& aRegionName,
                                             const string& aSectorName,
                                             const int aPeriod ) const
{
    double inputCost = Technology::getTotalInputCost( aRegionName, 
        aSectorName,
        aPeriod );

    // add enrichment, fabrication, interim storage, permanent disposal,
    // and reprocessing costs
    inputCost += ( enrichmentCost * getSWUperProduct() + fabricationCost + interimStorageCost
        + geologicWasteDisposalCost + reprocessingCost ) * getInitialMass();
    // if blanket material is necessary to breed fissile material,
    // add blanket cost to necost
    if (blanketFuelName != "none") {
        inputCost += ( blanketFabCost + interimStorageCost
            + geologicWasteDisposalCost + reprocessingCost ) * getInitialMass() * blanketFuelRatio;
    }
    return inputCost;
}

//! Calculates mutiple inputs for nuclear fuel technology and output.
/*! Adds demands for fuels and ghg emissions to markets in the marketplace
* \author Sonny Kim
*/
void NukeFuelTechnology::production( const string& aRegionName,
                                    const string& aSectorName, 
                                    double aVariableDemand,
                                    double aFixedOutputScaleFactor,
                                    const GDP* aGDP,
                                    const int aPeriod )
{
    // Can't have a scale factor and positive demand.
    assert( aFixedOutputScaleFactor == 1 || aVariableDemand == 0 );

    // Can't have negative variable demand.
    assert( aVariableDemand >= 0 && util::isValidNumber( aVariableDemand ) );

    // Check for positive variable demand and positive fixed output.
    assert( mFixedOutput == getFixedOutputDefault() || util::isEqual( aVariableDemand, 0.0 ) );

    // Check that a state has been created for the period.
    assert( aPeriod < static_cast<int>( mProductionState.size() ) );

    // Construct a marginal profit calculator. This allows the calculation of 
    // marginal profits to be lazy.
    MarginalProfitCalculator marginalProfitCalc( this );

    double primaryOutput =
        mProductionState[ aPeriod ]->calcProduction( aRegionName,
        aSectorName,
        aVariableDemand,
        &marginalProfitCalc,
        aFixedOutputScaleFactor,
        mShutdownDeciders,
        aPeriod );

    // Calculate input demand.
    mProductionFunction->calcDemand( mInputs, primaryOutput, aRegionName, aSectorName,
        1, aPeriod, 0, mAlphaZero );

    // add demand for fertile material
    Marketplace* marketplace = scenario->getMarketplace();
    if( fertileFuelName != "none" ) {
        double inputFertile = primaryOutput / getFertileEfficiency( aPeriod );
        marketplace->addToDemand( fertileFuelName, aRegionName, inputFertile, aPeriod );
    }
    // add demand for blanket material
    if( blanketFuelName != "none" ) {
        double inputBlanket = primaryOutput / getBlanketEfficiency( aPeriod );
        marketplace->addToDemand( blanketFuelName, aRegionName, inputBlanket, aPeriod );
    }

    // calculate by-products from technology (shk 10/11/04) mass of initial
    // heavy metal (mass of nuclear fuel) output unit (EJ) is thermal energy
    // content of nuclear fuel kg/GJ equal to MMT/EJ. byProduct coeffcient is
    // composition of waste in percentage.
    calcEmissionsAndOutputs( aRegionName, primaryOutput, aGDP, aPeriod );
}

double NukeFuelTechnology::calcShare( const string& aRegionName,
                                     const string& aSectorName,
                                     const GDP* aGDP,
                                     const double aLogitExp,
                                     const int aPeriod ) const
{
    return Technology::calcShare( aRegionName, aSectorName, aGDP, aLogitExp, aPeriod );
}


/*! \brief Value function for calculating SWU
*
* Function represents value function for SWU 
*
* \author Sonny Kim
* \param aWeightFraction enrichment weight fraction
* \return The result of the value function.
*/
double NukeFuelTechnology::getSWValue( const double aWeightFraction ) {
    // TODO: Some input files seem to have this...?
    if( aWeightFraction >= 1 ){
        return 0;
    }

    return ( 2 * aWeightFraction - 1 ) * log( aWeightFraction / ( 1 - aWeightFraction ) );
}
/*! \brief Separative Work Unit per unit product
*
* Function calculates separative work unit per unit product using value function
*
* \author Sonny Kim
* \return The separtive work unit per unit product.
*/
double NukeFuelTechnology::getSWUperProduct() const {
    double S = 0;
    // do only if enrichment has a cost and is therefore necessary 
    if( enrichmentCost > 0) {
        S = getSWValue( enrichmentProd ) - getSWValue( enrichmentTail )
            - getFeedProductRatio()
            * ( getSWValue( enrichmentFeed ) - getSWValue( enrichmentTail ) );
        assert( S > 0 );
    }
    return S;
}

/*! \brief Feed to product ratio for nuclear fuel.
*
* Function calculate the feed to product ratio of nuclear fuel.
*
* \author Sonny Kim
* \return The feed to product ratio.
*/
double NukeFuelTechnology::getFeedProductRatio() const {
    /*! \pre Input enrichment is positive. */
    assert( enrichmentFeed - enrichmentTail > 0 );

    /*! \pre Output enrichment is positive. */
    assert( enrichmentProd - enrichmentTail > 0 );

    return ( enrichmentProd - enrichmentTail )
        /( enrichmentFeed - enrichmentTail );
}

/*! \brief Get the calculated mass of initial heavy metal per unit thermal
*          energy of nuclear fuel.
* \details Calculate amount of fuel per heat output in (kgHM/GJ) or inverse of
*          burnup assuming burnup unit in GWd/MetricTonHM.
* \return The initial mass.
*/
double NukeFuelTechnology::getInitialMass() const {
    const unsigned int HOURS_PER_DAY = 24;
    const double GWH_TO_GJ = 2.78E-4;
    const double KG_PER_METRIC_TON = 1000;

    /*! \pre Burnup must be positive. */
    assert( burnup > 0 );

    double initialMass = 1 / ( burnup * HOURS_PER_DAY / GWH_TO_GJ / KG_PER_METRIC_TON );

    /*! \post Initial mass is positive. */
    assert( initialMass > 0 );
    return initialMass;
}

/*! \brief Get the fertile efficiency.
* \param aPeriod Model period.
* \return The fertile efficiency.
*/
double NukeFuelTechnology::getFertileEfficiency( const int aPeriod ) const {
    // if reprocessed fuel, fuel is combination of fissile and fertile material
    double fertileEff = 1;
    if ( fertileFuelName != "none" ) {
        double kgFertilePerGJ = getInitialMass() * ( 1 - enrichmentProd );

        /*! \pre Fertile mass per GJ is positive.*/
        assert( kgFertilePerGJ > 0 );

        fertileEff = 1 / kgFertilePerGJ;
    }
    /*! \post Fertile efficiency is positive. */
    return fertileEff;
}

/*! \brief Get the blanket efficiency.
* \param aPeriod Model period.
* \return The fertile efficiency.
*/
double NukeFuelTechnology::getBlanketEfficiency( const int aPeriod ) const {
    // if blanket material is necessary to breed fissile material.
    double blanketEff = 1;
    if (blanketFuelName != "none") {
        // multiple of fuel amount
        double kgBlanketPerGJ = getInitialMass() * blanketFuelRatio;
        blanketEff = 1 / kgBlanketPerGJ;
    }
    /*! \post Blanket efficiency must be greater than zero. */
    assert( blanketEff > 0 );
    return blanketEff;
}
