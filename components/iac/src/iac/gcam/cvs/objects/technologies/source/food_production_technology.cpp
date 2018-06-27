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
* \file food_production_technology.cpp
* \ingroup Objects
* \brief FoodProductionTechnology class source file.
* \author James Blackwood
*/

#include "util/base/include/definitions.h"
#include "technologies/include/food_production_technology.h"
#include "land_allocator/include/iland_allocator.h"
#include "emissions/include/aghg.h"
#include "containers/include/scenario.h"
#include "util/base/include/xml_helper.h"
#include "marketplace/include/marketplace.h"
#include "containers/include/iinfo.h"
#include "technologies/include/ical_data.h"
#include "technologies/include/iproduction_state.h"
#include "technologies/include/ioutput.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! 
 * \brief Constructor.
 * \param aName Technology name.
 * \param aYear Technology year.
 */
FoodProductionTechnology::FoodProductionTechnology( const string& aName, const int aYear )
:Technology( aName, aYear ){
    mLandAllocator = 0;
    variableCost = -1;
    mNonLandCostTechChange = 0;
    calLandUsed = -1;
    calYield = -1;
    mMaxYield = -1;
    calObservedYield = -1;
    agProdChange = 0;
    mAboveGroundCarbon = 0;
    mBelowGroundCarbon = 0;
    mMatureAge = 1;        
    mHarvestedToCroppedLandRatio = 1;
    mLUCPenaltyMultiplier = 0;
}

// ! Destructor
FoodProductionTechnology::~FoodProductionTechnology() {
}

//! Parses any input variables specific to derived classes
bool FoodProductionTechnology::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if ( nodeName == "variableCost" ) {
        variableCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "NonLandCostTechChange" ) {
        mNonLandCostTechChange = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "landType" ) {
        landType = XMLHelper<string>::getValue( curr );
    }
    else if( nodeName == "calLandUsed" ) {
        calLandUsed = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "calYield" ) {
        calYield = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "maxYield" ) {
        mMaxYield = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "agProdChange" ) {
        agProdChange = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "above-ground-carbon" ){
        mAboveGroundCarbon = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "harvested-to-cropped-land-ratio" ){
        mHarvestedToCroppedLandRatio = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "luc-carbon-penalty-multiplier" ){
        mLUCPenaltyMultiplier = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "below-ground-carbon" ){
        mBelowGroundCarbon = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "mature-age" ) {
        mMatureAge = XMLHelper<int>::getValue( curr );        
    }
    else {
        return false;
    }
    return true;
}

/*! \brief Derived class visitor.
*
* This may have already been implemented in multi-inputs version, so replace with that if so.
* \author Steve Smith
*/
void FoodProductionTechnology::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
    // Derived visit.
    aVisitor->startVisitFoodProductionTechnology( this, aPeriod );
    // End the derived class visit.
    aVisitor->endVisitFoodProductionTechnology( this, aPeriod );
}

//! write object to xml output stream
void FoodProductionTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( landType, "landType", out, tabs );
    XMLWriteElement( variableCost, "variableCost", out, tabs );
    XMLWriteElementCheckDefault( mNonLandCostTechChange, "NonLandCostTechChange", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( calYield, "calYield", out, tabs, -1.0 );
    if ( getXMLName1D() != "ForestProductionTechnology" ) {
        XMLWriteElementCheckDefault( calLandUsed, "calLandUsed", out, tabs, -1.0 );
    }
    XMLWriteElementCheckDefault( agProdChange, "agProdChange", out, tabs, 0.0 );   
    XMLWriteElementCheckDefault( mHarvestedToCroppedLandRatio, "harvested-to-cropped-land-ratio", out, tabs, 1.0 );
    XMLWriteElementCheckDefault( mLUCPenaltyMultiplier, "luc-carbon-penalty-multiplier", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mAboveGroundCarbon, "above-ground-carbon", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mBelowGroundCarbon, "below-ground-carbon", out, tabs, 0.0 );
    XMLWriteElement( mMatureAge, "mature-age", out, tabs );
    XMLWriteElement( mMaxYield, "maxYield", out, tabs );
}

//! write object to xml output stream
void FoodProductionTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( landType, "landType", out, tabs );
    XMLWriteElement( variableCost, "variableCost", out, tabs );
    XMLWriteElement( calYield, "calYield", out, tabs );
    XMLWriteElement( calObservedYield, "calObservedYield", out, tabs );
    XMLWriteElement( calLandUsed, "calLandUsed", out, tabs );
    XMLWriteElement( agProdChange, "agProdChange", out, tabs );
    XMLWriteElement( mHarvestedToCroppedLandRatio, "harvested-to-cropped-land-ratio", out, tabs );
    XMLWriteElement( mLUCPenaltyMultiplier, "luc-carbon-penalty-multiplier", out, tabs );
    XMLWriteElement( mAboveGroundCarbon, "above-ground-carbon", out, tabs );
    XMLWriteElement( mBelowGroundCarbon, "below-ground-carbon", out, tabs );
    XMLWriteElement( mMatureAge, "mature-age", out, tabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& FoodProductionTechnology::getXMLName1D() const {
    return getXMLNameStatic1D();
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
const string& FoodProductionTechnology::getXMLNameStatic1D() {
    const static string XML_NAME1D = "FoodProductionTechnology";
    return XML_NAME1D;
}

//! Clone Function. Returns a deep copy of the current technology.
FoodProductionTechnology* FoodProductionTechnology::clone() const {
    return new FoodProductionTechnology( *this );
}

void FoodProductionTechnology::initCalc( const string& aRegionName,
                                         const string& aSectorName,
                                         const IInfo* aSubsectorInfo,
                                         const Demographic* aDemographics,
                                         PreviousPeriodInfo& aPrevPeriodInfo,
                                         const int aPeriod )
{
    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
                          aDemographics, aPrevPeriodInfo, aPeriod );
  
    const Modeltime* modeltime = scenario->getModeltime();

    // Only setup calibration information if this is the initial year of the
    // technology.
    if( !mProductionState[ aPeriod ]->isNewInvestment() ){
        return;
    }

    // If calibration data is present for this year, then zero out ag prod change
    // for previous periods 
    if ( mCalValue.get() ) {
        for ( int pastPeriod = 0; pastPeriod <= aPeriod; ++pastPeriod ) {
            mLandAllocator->applyAgProdChange( landType, mName, 0.0, pastPeriod, pastPeriod );
        }
    }
    
    // Apply technical change.
    mLandAllocator->applyAgProdChange( landType, mName, agProdChange, aPeriod, aPeriod );
    
    // TODO: Use a better method of passing forward calibration information.
    // Since the market may be global, create a unique regional string for the
    // calibrated variable cost. This is hacked to work for multiple technologies too
    // This needs to be fixed ASAP.
    const string preVarCostName = "preVarCost-" + mName + "-" + aRegionName;
    double calVarCost;
    double preVarCost = 0.0;

    // Get the information object for this market.
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( aSectorName, aRegionName, aPeriod, true );
    assert( marketInfo );
    if( calObservedYield != -1 ) {
        double calPrice = marketInfo->getDouble( "calPrice", true );

        // Calculate the calibrated variable cost. (Yield is adjusted from agronomic to economic yield per acre)
        // TODO: This is the only access of this variable from outside the AgLU. Change this function
        // to getCalNumeraireAveObservedRate.
        calVarCost = calPrice - mLandAllocator->getAverageLandProfitRate( aPeriod, landType )
                                / calcDiscountFactor()
                                / ( calObservedYield * mHarvestedToCroppedLandRatio );
       // assert( util::isValidNumber( calVarCost ) );

        // Set the variable cost for the technology to the calibrated variable cost.
        if ( calVarCost > util::getSmallNumber() ) {
            // TODO: Add warning if there was a read-in variable cost.
            variableCost = calVarCost;
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Read in value for calPrice in " << aRegionName << " "
                    << mName << " is too low by " << fabs( calVarCost ) << endl;
        }
    }
    else {
        // Get the calibrated variable cost from the market info.
        if ( variableCost == -1 ) {
            preVarCost = marketInfo->getDouble( preVarCostName, aPeriod != 0 ? true : false );
            if ( preVarCost > util::getSmallNumber() ) {
                // Adjust last period's variable cost by tech change
                int timestep = modeltime->gettimestep( aPeriod );
                variableCost = preVarCost / pow(1 + mNonLandCostTechChange , timestep);
            }
            else if( aPeriod > 0 ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::DEBUG );
                mainLog << "No variable cost or caldata Read " << aRegionName << " "
                    << mName << endl;
            }
        }
    }

    // If this is not the end period of the model, set the market info
    // variable for the next period.
    if( aPeriod + 1 < modeltime->getmaxper() ){
        IInfo* nextPerMarketInfo = marketplace->getMarketInfo( aSectorName, aRegionName, aPeriod + 1, true );
        assert( nextPerMarketInfo );
        nextPerMarketInfo->setDouble( preVarCostName, variableCost );
    }

    // Set the above and below ground carbon for this technology.
    // TODO: This may need to be moved if the carbon content is calculated dynamically.
    mLandAllocator->setCarbonContent( landType, mName, mAboveGroundCarbon, mBelowGroundCarbon, mMatureAge, aPeriod );

    mLandAllocator->setMaxYield( landType, mName, mMaxYield, aPeriod );
    
}

void FoodProductionTechnology::postCalc( const string& aRegionName,
                                         const int aPeriod )
{
    Technology::postCalc( aRegionName, aPeriod );
}

void FoodProductionTechnology::completeInit( const std::string& aRegionName,
                                             const std::string& aSectorName,
                                             const std::string& aSubsectorName,
                                             DependencyFinder* aDepFinder,
                                             const IInfo* aSubsectorInfo,
                                             ILandAllocator* aLandAllocator,
                                             const GlobalTechnologyDatabase* aGlobalTechDB )
{
    // Store away the land allocator.
    mLandAllocator = aLandAllocator;

    // Setup the land allocators for the secondary outputs
    if ( mOutputs.size() ) {
        // Technology::completeInit() will add the primary output.
        // At this point, all are secondary outputs
        for ( vector<IOutput*>::iterator outputIter = mOutputs.begin(); outputIter != mOutputs.end(); ++outputIter ) {
           ( *outputIter )->setLandAllocator( aLandAllocator, mName, landType );
        }
    }

    // Note: Technology::completeInit() loops through the outputs.
    //       Therefore, if any of the outputs need the land allocator,
    //       the call to Technology::completeInit() must come afterwards
    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder, aSubsectorInfo,
                              aLandAllocator, aGlobalTechDB );

    // Setup the land usage for this production.
    int techPeriod = scenario->getModeltime()->getyr_to_per( year );
    
    // Check for existance of land allocator so that can gracefully exit if does not exist
    if ( mLandAllocator ) {
        mLandAllocator->addLandUsage( landType, mName, ILandAllocator::eCrop, techPeriod );
    }
    else {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::SEVERE );
        mainLog << "No land allocator is present for Ag Production technology. Exiting model." << endl;
        exit( -1 );
    }
    
    // Technical change may only be applied after the base period.
    if( agProdChange > util::getSmallNumber() && mCalValue.get() )
    {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Food production technologies may not have technical change"
                << " in a calibration period." << endl;
        agProdChange = 0;
    }

    if( mHarvestedToCroppedLandRatio < util::getSmallNumber() )
    {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid value of harvested-to-cropped-land-ratio. Reset to 1." << endl;
        mHarvestedToCroppedLandRatio = 1.0;
    }

    setCalLandValues();
}

/*! \brief Sets calibrated land values to land allocator.
*
* This utility function is called once for food technologies and twice for forest 
* technologies (see forest version of this function). Call in completeInit sets initial
* land-use and calibration values in the land allocator.
*
* \author Steve Smith
*/
void FoodProductionTechnology::setCalLandValues() {

    // TODO: We should change the input so that only one of these options can be read
    //       as they are mutually exclusive.
    if ( mCalValue.get() && calLandUsed != -1 ) {

        // KVC: prevent divide by zero errors.
        if ( calLandUsed > 0 ){
            calObservedYield = mCalValue->getCalOutput() / calLandUsed;
        }
        else {
            calObservedYield = -1.0;
        }

        // Setup the land usage for this production.
        int techPeriod = scenario->getModeltime()->getyr_to_per( year );

        // Warn the user that the calibrated yield will not be used since an
        // observed yield can be calculated.
        if( calYield != -1 ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Calibrated yield will be overridden by the observed yield." << endl;
        }

        // Want to pass in yield in units of GCal/kHa
        mLandAllocator->setCalLandAllocation( landType, mName, calLandUsed / mHarvestedToCroppedLandRatio, techPeriod, techPeriod );
        mLandAllocator->setCalObservedYield( landType, mName, calObservedYield * mHarvestedToCroppedLandRatio, techPeriod );
    } 
    else if ( calYield != -1 ) {
        int techPeriod = scenario->getModeltime()->getyr_to_per( year );
        mLandAllocator->setCalObservedYield( landType, mName, calYield * mHarvestedToCroppedLandRatio, techPeriod );
    }
}

/*!
* \brief Calculate unnormalized technology unnormalized shares.
* \details Since food and forestry technologies are profit based, they do not
*          directly calculate a share. Instead, their share of total supply is
*          determined by the sharing which occurs in the land allocator. To
*          facilitate this the technology sets the intrinsic rate for the land
*          use into the land allocator. The technology share itself is set to 1.
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
* \return The technology share, always 1 for FoodProductionTechnologies.
* \author James Blackwood, Steve Smith
*/
double FoodProductionTechnology::calcShare( const string& aRegionName,
                                            const string& aSectorName,
                                            const GDP* aGDP,
                                            const int aPeriod ) const
{
    assert( mProductionState[ aPeriod ]->isNewInvestment() );
    
    // Food production technologies are profit based, so the amount of output
    // they produce is independent of the share.
    return 1;
}

void FoodProductionTechnology::calcCost( const string& aRegionName,
                                         const string& aSectorName,
                                         const int aPeriod )
{
    if( !mProductionState[ aPeriod ]->isOperating() ){
        return;
    }

    // Calculate the technology's profit rate. This rate is in $/GCal
    double profitRate = calcProfitRate( aRegionName, aSectorName, aPeriod );

    // If yield is GCal/Ha and prices are $/GCal, then rental rate is $/Ha
    // Passing in rate as $/GCal and setIntrinsicRate will transform it to  $/Ha inside the land leaf.
    mLandAllocator->setIntrinsicRate( aRegionName, landType, mName, profitRate, aPeriod );

    // Override costs to a non-zero value as the cost for a food production
    // technology is not used for the shares.
    mCosts[ aPeriod ] = 1;
}

double FoodProductionTechnology::getNonEnergyCost( const string& aRegionName,
                                                   const int aPeriod ) const
{
    return 0;
}

/*! \brief Calculates the output of the technology.
* \details Calculates the amount of current forestry output based on the amount
*          of planted forestry land and it's yield. Forestry production
*          technologies are profit based and determine their supply
*          independently of the passed in subsector demand. However, since this
*          is a solved market, in equilibrium the sum of the production of
*          technologies within a sector will equal the demand for the sector.
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aVariableDemand Subsector demand for output.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
* \todo Need better way to deal with biomass units than the hard coded
*       conversion below. 
*/
void FoodProductionTechnology::production( const string& aRegionName,
                                           const string& aSectorName,
                                           const double aVariableDemand,
                                           const double aFixedOutputScaleFactor,
                                           const GDP* aGDP,
                                           const int aPeriod )
{
    // TODO: Food production technologies are not currently vintaged.
    if( !mProductionState[ aPeriod ]->isOperating() ){
        // Set physical output to zero.
        mOutputs[ 0 ]->setPhysicalOutput( 0, aRegionName,
                                          mCaptureComponent.get(),
                                          aPeriod );
        return;
    }

    // Calculate the profit rate.
    double profitRate = calcProfitRate( aRegionName, aSectorName, aPeriod );

    // Calculate the yield.
    mLandAllocator->calcYield( landType, mName, aRegionName, 
                               profitRate, aPeriod, aPeriod );

    // Calculate the output of the technology.
    double primaryOutput = calcSupply( aRegionName, aSectorName, aPeriod );

    // This output needs to be in EJ instead of GJ.
    // TODO: Fix this once units framework is complete.
    if( mName == "biomass" ) {
        primaryOutput /= 1e9;
    }

    // Set the input to be the land used.
    // TODO: Add a land input.
    //mInput[ aPeriod ] = mLandAllocator->getLandAllocation( landType, mName, aPeriod );
    calcEmissionsAndOutputs( aRegionName, primaryOutput, aGDP, aPeriod );

    // TODO -- need to set the input to be the land used for multi-inputs version.
    // This will be necessary before GHG's will work in multi-input for ag production sectors.
}

/*! \brief Calculate the profit rate for the technology.
* \details Calculates the profit rate which is equal to the market price minus
*          the variable cost.
*          Profit rate can be negative.
* \param aRegionName Region name.
* \param aProductName Name of the product for which to calculate the profit
*        rate. Must be an output of the technology.
* \return The profit rate.
*/
double FoodProductionTechnology::calcProfitRate( const string& aRegionName,
                                                 const string& aProductName,
                                                 const int aPeriod ) const
{
    // Conversion from 1990 to 1975 dollars
    const double CVRT_75_TO_90 = 2.212;
    
    // Calculate profit rate.
    const Marketplace* marketplace = scenario->getMarketplace();

    // TODO: Units here will be wrong for anything other than biomass because prices will be in $/Gcal
    // as GHG costs/profits are always in $/GJ.
    double secondaryValue = calcSecondaryValue( aRegionName, aPeriod );
    double profitRate = ( marketplace->getPrice( aProductName, aRegionName, aPeriod ) + secondaryValue ) 
                         * CVRT_75_TO_90 - variableCost;

    // Add carbon penalty
    // Check if a carbon market exists and has a non-zero price.
    double carbonPrice = marketplace->getPrice( "CO2", aRegionName, aPeriod, false );
    if( carbonPrice != Marketplace::NO_MARKET_PRICE ){
       profitRate -= mLUCPenaltyMultiplier * carbonPrice;
    }
    
    return profitRate;
}

/*! \brief Calculate the factor to discount between the present period and the
*          harvest period.
* \return The discount factor.
*/
double FoodProductionTechnology::calcDiscountFactor() const {
    // Food products are produced in a single year, so they do not have a
    // discount factor.
    return 1;
}

/*! \brief Calculate the supply for the technology.
* \details Calculates the food produced which is equal to the yield multiplied
*          by the land allocated.
* \param aProductName Product name.
* \param aPeriod Period.
* \return The supply produced by the technology.
*/
double FoodProductionTechnology::calcSupply( const string& aRegionName,
                                             const string& aProductName,
                                             const int aPeriod ) const
{
    // Get yield per acre of land
    double economicYield = mLandAllocator->getYield( landType, mName, aPeriod );
    assert( economicYield >= 0 );
    
    // Convert to agronomic yield, which is per harvest
    // TODO: Should we write-out agronomic yield to XML DB? Might be useful to have.
    double agronomicYield = economicYield / mHarvestedToCroppedLandRatio;
    
    double landAllocation = mLandAllocator->getLandAllocation( landType, mName, aPeriod );
    
    // Convert from physical acres of land to land allocated
    double harvestedLand = landAllocation * mHarvestedToCroppedLandRatio; 
    
    // Check that if yield is zero the land allocation is zero.
    // TODO: Determine why a small number is too large.
    // TODO: Checking the variable cost of zero is a hack so that this works
    //       for the unmanaged sector.
    if( agronomicYield < util::getSmallNumber() && landAllocation > 0.1 && variableCost > util::getTinyNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Zero production of " << aProductName << " by technology " << mName
                << " in region " << aRegionName << " with a positive land allocation of "
                << landAllocation << "." << endl;
    }

    assert( agronomicYield * harvestedLand >= 0 );

    // Set output to yield times amount of land.
    return agronomicYield * harvestedLand;
}

double FoodProductionTechnology::getTotalInputCost( const string& aRegionName, const string& aSectorName,
                                       const int aPeriod ) const
{
    return Technology::getTotalInputCost( aRegionName, aSectorName, aPeriod );
}

