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
* \file ag_production_technology.cpp
* \ingroup Objects
* \brief AgProductionTechnology class source file.
* \author Marshall Wise, Kate Calvin
*/

#include "util/base/include/definitions.h"
#include "technologies/include/ag_production_technology.h"
#include "land_allocator/include/iland_allocator.h"
#include "land_allocator/include/aland_allocator_item.h"
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
AgProductionTechnology::AgProductionTechnology( const string& aName, const int aYear )
:Technology( aName, aYear ),
   mLandAllocator( 0 ),
   mProductLeaf( 0 ),
   mNonLandVariableCost( 0 ), 
   mNonLandCostTechChange( 0),
   mAgProdChange ( 0 ),
   mHarvestsPerYear( 1 ),
   mYield ( 0 )
{
}

/*!
 * \brief Copy constructor.
 * \details Does not copy variables which should get initialized through normal
 *          model operations.
 * \param aAgTech Tech AgProductionTechnology to copy.
 */
AgProductionTechnology::AgProductionTechnology( const AgProductionTechnology& aAgTech )
:Technology( aAgTech ),
mNonLandVariableCost( aAgTech.mNonLandVariableCost ),
mNonLandCostTechChange( aAgTech.mNonLandCostTechChange ),
mAgProdChange( aAgTech.mAgProdChange ),
mHarvestsPerYear( aAgTech.mHarvestsPerYear ),
// The following do not get copied as they are initialized through other means
mLandAllocator( 0 ),
mProductLeaf( 0 ),
mYield( 0 ),
mYieldAdjForDensity( 1.0 )
{
}

// ! Destructor
AgProductionTechnology::~AgProductionTechnology() {
}

//! Parses any input variables specific to derived classes
bool AgProductionTechnology::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    if ( nodeName == "nonLandVariableCost" ) {
        mNonLandVariableCost = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "nonLandCostTechChange" ) {
        mNonLandCostTechChange = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "yield" ) {
        mYield = XMLHelper<double>::getValue( curr );
    }    
    else if( nodeName == "agProdChange" ) {
        mAgProdChange = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "harvests-per-year" ){
        mHarvestsPerYear = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

/*! \brief Derived class visitor.
*
* 
*/
void AgProductionTechnology::acceptDerived( IVisitor* aVisitor, const int aPeriod ) const {
    // Derived visit.
    aVisitor->startVisitAgProductionTechnology( this, aPeriod );
    // End the derived class visit.
    aVisitor->endVisitAgProductionTechnology( this, aPeriod );
}

//! write object to xml output stream
void AgProductionTechnology::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mNonLandVariableCost, "nonLandVariableCost", out, tabs );
    XMLWriteElementCheckDefault( mNonLandCostTechChange, "nonLandCostTechChange", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mYield, "yield", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( mAgProdChange, "agProdChange", out, tabs, 0.0 );   
    XMLWriteElementCheckDefault( mHarvestsPerYear, "harvests-per-year", out, tabs, 1.0 );
}

//! write object to xml output stream
void AgProductionTechnology::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteElement( mNonLandVariableCost, "nonLandVariableCost", out, tabs );
    XMLWriteElement( mYield, "yield", out, tabs );
    XMLWriteElement( mAgProdChange, "agProdChange", out, tabs );
    XMLWriteElement( mHarvestsPerYear, "harvests-per-year", out, tabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& AgProductionTechnology::getXMLName() const {
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
const string& AgProductionTechnology::getXMLNameStatic() {
    const static string XML_NAME = "AgProductionTechnology";
    return XML_NAME;
}

//! Clone Function. Returns a deep copy of the current technology.
AgProductionTechnology* AgProductionTechnology::clone() const {
    return new AgProductionTechnology( *this );
}

void AgProductionTechnology::initCalc( const string& aRegionName,
                                         const string& aSectorName,
                                         const IInfo* aSubsectorInfo,
                                         const Demographic* aDemographics,
                                         PreviousPeriodInfo& aPrevPeriodInfo,
                                         const int aPeriod )
{
    Technology::initCalc( aRegionName, aSectorName, aSubsectorInfo,
                          aDemographics, aPrevPeriodInfo, aPeriod );
  
    const Modeltime* modeltime = scenario->getModeltime();

    // Only do tech changes if this is the initial year of the
    // technology.
    if( !mProductionState[ aPeriod ]->isNewInvestment() ){
        return;
    }

    // Compute tech change values for this period for both ag productivity and
    // the nonLandVariableCost.  Since technologies are distinct by vintage and don't
    // previous period technologies, need to save a previous period compounded cumulative
    // change in the MarketInfo
 
    // Create a unique regional string for the yield and for the
    // variable cost. 

    const string preVarCostName = "preVarCost-" + mName + "-" + aRegionName;
    double preVarCost = 0.0;

    const string preYieldName = "preYield-" + mName + "-" + aRegionName;
    double preYield = 0.0;


    // Get the information object for this market.
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( aSectorName, aRegionName, aPeriod, true );
    assert( marketInfo );

    // If no nonLandVariableCost is read in, get the previous period cost from the market info.
    // Note: you can never overwrite a positive yield with a zero yield. If the model sees a
    // zero non-land cost, it will copy from the previous period.
    if ( mNonLandVariableCost == 0 && aPeriod != 0 ) {
         preVarCost = marketInfo->getDouble( preVarCostName, true );
         // Adjust last period's variable cost by tech change
         int timestep = modeltime->gettimestep( aPeriod );
         mNonLandVariableCost = preVarCost / pow(1 + mNonLandCostTechChange , timestep);
    }

    // Only do the ag productivity change calc if a calibration value is not read in that period
    if( !mCalValue.get() ){       
        // Unless a yield is read in for this period, get the previous period yield from the market info.
        // Note: you can never overwrite a positive yield with a zero yield. If the model sees a
        // zero yield, it will copy from the previous period.
        if ( mYield == 0 && aPeriod != 0 ) {
            preYield = marketInfo->getDouble( preYieldName, true );
            // Adjust last period's variable cost by tech change
            int timestep = modeltime->gettimestep( aPeriod );
            mYield = preYield * pow(1 + mAgProdChange , timestep);
        }
    }

    // If this is not the end period of the model, set the market info
    // variable cost and Yield for the next period.
    if( aPeriod + 1 < modeltime->getmaxper() ){
        IInfo* nextPerMarketInfo = marketplace->getMarketInfo( aSectorName, aRegionName, aPeriod + 1, true );
        assert( nextPerMarketInfo );
        nextPerMarketInfo->setDouble( preVarCostName, mNonLandVariableCost );
        nextPerMarketInfo->setDouble( preYieldName, mYield );
    }

    // Adjust yield for changes in carbon density.
    // Note that we want to make sure this adjustment is for this year's yield only
    mYield *= mYieldAdjForDensity;

    // We need to and do set the profit rate in the land allocator
    //           so that it can calculate share weights.
    // Calculate the technology's profit rate. This rate is in $/GCal
    // Note that we are hard coding period 0 so that the cal price is used
    // during calibration.  Another approach may be worth consideration.
    double profitRate = calcProfitRate( aRegionName, aSectorName, 0 );

    // If yield is GCal/kHa and prices are $/GCal, then rental rate is $/kHa
    // And this is what is now passed in ($/kHa)
    mLandAllocator->setProfitRate( aRegionName, mName, profitRate, aPeriod );
}

void AgProductionTechnology::postCalc( const string& aRegionName,
                                         const int aPeriod )
{
    Technology::postCalc( aRegionName, aPeriod );
}

void AgProductionTechnology::completeInit( const std::string& aRegionName,
                                             const std::string& aSectorName,
                                             const std::string& aSubsectorName,
                                             DependencyFinder* aDepFinder,
                                             const IInfo* aSubsectorInfo,
                                             ILandAllocator* aLandAllocator )
{
    // Store away the land allocator.
    mLandAllocator = aLandAllocator;
    mProductLeaf = aLandAllocator->findProductLeaf( mName );
 
    // Send "pointer to the land allocator" to each of the secondary outputs, e.g, residue biomass
    if ( mOutputs.size() ) {
        // Note: the only outputs in the mOutputs vector at this point are secondary outputs.
        //       The primary output isn't initialized until Technology::completeInit() is called, 
        //       which happens after this step.
        for ( vector<IOutput*>::iterator outputIter = mOutputs.begin(); outputIter != mOutputs.end(); ++outputIter ) {
            ( *outputIter )->sendLandAllocator( aLandAllocator, mName );
        }
    }

    // Note: Technology::completeInit() loops through the outputs.
    //       Therefore, if any of the outputs need the land allocator,
    //       the call to Technology::completeInit() must come afterwards
    Technology::completeInit( aRegionName, aSectorName, aSubsectorName, aDepFinder, aSubsectorInfo,
                              aLandAllocator );

    // Make some tests for bad inputs

    // Technical change may only be applied after the base period.
    if( mAgProdChange > 0.0 && mCalValue.get() )
    {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Ag technologies may not have technical change"
                << " in a calibration period."
                << aRegionName << " " << mName << endl;
        mAgProdChange = 0;
    }

    if( mHarvestsPerYear < 0.0 )
    {
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid value of harvests-per-year. Reset to 1."
                << aRegionName << " " << mName << endl;
        mHarvestsPerYear = 1.0;
    }

    setCalYields( aRegionName );
}

/*! \brief Sets the calibrated yields
*
* Call in completeInit sets initial
* land-use and calibration values in the land allocator.
*
* \author Marshall Wise
*/
void AgProductionTechnology::setCalYields(const std::string& aRegionName) {
    // if a calibrated output is read in for this period, use it to compute yield	
    if ( mCalValue.get() ) {
        // technology knows the year it started in member variable "year"
        int techPeriod = scenario->getModeltime()->getyr_to_per( year );
        double calLandUsed = mLandAllocator->getLandAllocation( mName, techPeriod );
        // if land is also read in, compute yield, else write a warning and set
        // yield to 0
        if ( calLandUsed > 0 ) {
            mYield = mCalValue->getCalOutput() / calLandUsed;
        }
        else if ( mCalValue->getCalOutput() > 0 ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Caloutput read in but no land read in for technology"
                << aRegionName << " " << mName << endl;
            mYield=0;
        }
    }
}

/*!
* \brief Calculate unnormalized technology unnormalized shares.
* \details Since ag technologies compute output based on land, they do not
*          directly calculate a share. Instead, their total supply is
*          determined by the sharing which occurs in the land allocator. To
*          facilitate this the technology sets the profit rate for the land
*          use into the land allocator. The technology share itself is set to 1
           but is not used.
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
* \return The technology share, always 1 for AgProductionTechnologies.
* \author James Blackwood, Steve Smith
*/
double AgProductionTechnology::calcShare( const string& aRegionName,
                                            const string& aSectorName,
                                            const GDP* aGDP,
                                            const double aLogitExp,
                                            const int aPeriod ) const
{
    assert( mProductionState[ aPeriod ]->isNewInvestment() );
    
    // Ag production technologies of output is determined by land amount
    // and yield, so the share among technologies is not used.
    return 1;
}


/* agTechnologies are not shared on cost, so this calCost method is overwritten
   by a calculation of technology profit which is passed to the land allocator
   where it is used for sharing land.  */

void AgProductionTechnology::calcCost( const string& aRegionName,
                                         const string& aSectorName,
                                         const int aPeriod )
{
    if( !mProductionState[ aPeriod ]->isOperating() ){
        return;
    }

    // Calculate the technology's profit rate. This rate is in $/m2
    double profitRate = calcProfitRate( aRegionName, aSectorName, aPeriod );

    mProductLeaf->setProfitRate( aRegionName, mName, profitRate, aPeriod );

    // Override costs to a non-zero value as the cost for a ag production
    // technology is not used for the shares.  
    mCosts[ aPeriod ] = 1;
}

double AgProductionTechnology::getNonEnergyCost( const string& aRegionName,
                                                   const int aPeriod ) const
{
    return 0;
}

/*! \brief Calculates the output of the technology.
* \details Calculates the amount of current ag output based on the amount
*          land and it's yield. 
* \param aRegionName Region name.
* \param aSectorName Sector name, also the name of the product.
* \param aVariableDemand Subsector demand for output.
* \param aGDP Regional GDP container.
* \param aPeriod Model period.
*/
void AgProductionTechnology::production( const string& aRegionName,
                                           const string& aSectorName,
                                           const double aVariableDemand,
                                           const double aFixedOutputScaleFactor,
                                           const GDP* aGDP,
                                           const int aPeriod )
{
    // This code is also in the technology class, but we do not
    // call the parent function so it must happen here too.
    if( !mProductionState[ aPeriod ]->isOperating() ){
        // Set physical output to zero.
        mOutputs[ 0 ]->setPhysicalOutput( 0, aRegionName,
                                          mCaptureComponent.get(),
                                          aPeriod );
        return;
    }

    // Calculate the profit rate.  
    // KVC_AGLU:: NOT SURE IF THIS IS NEEDED, SINCE IT ISN'T PASSED TO LAND ALLOCATOR
    // AT THIS POINT...THAT HAPPENS IN calcCost()
 //   double profitRate = calcProfitRate( aRegionName, aSectorName, aPeriod );

    // Calculate the output of the technology.
    double primaryOutput = calcSupply( aRegionName, aSectorName, aPeriod );

    // Calculate input demand.
    mProductionFunction->calcDemand( mInputs, primaryOutput, aRegionName, aSectorName,
        1, aPeriod, 0, mAlphaZero );
	
    // This call to the technology::calcEmissionsAndOutputs() is where the physical output
    // of the ag technology is set.
    calcEmissionsAndOutputs( aRegionName, primaryOutput, aGDP, aPeriod );

}

/*! \brief Calculate the profit rate for the technology.
* \details Calculates the profit rate which is equal to the market price minus
*          the variable cost.
*          Profit rate can be negative.
* \param aRegionName Region name.
* \param aProductName Name of the product for which to calculate the profit
*        rate. Must be an output of the technology.
* \return The profit rate.
\\ Profit rate is now in 1975$ per billion m2, so computation includes yield
*/
double AgProductionTechnology::calcProfitRate( const string& aRegionName,
                                                 const string& aProductName,
                                                 const int aPeriod ) const
{
    
    // Calculate profit rate.
    const Marketplace* marketplace = scenario->getMarketplace();

    // TODO: consider adding the residue biomass value to crop value
    // First, need to change residue biomass output as per unit of land
    // in order to prevent a simultaneity.  Then, we can include this value.
    double secondaryValue = calcSecondaryValue( aRegionName, aPeriod );

    // nonlandvariable cost units are now assumed to be in $/kg
    double price = marketplace->getPrice( aProductName, aRegionName, aPeriod );

    // Compute cost of variable inputs (such as water and fertilizer)
	// TODO: modify profit to reflect these costs.
    double inputCosts = getTotalInputCost( aRegionName, aProductName, aPeriod );

    // Price in model is 1975$/kg.  land and ag costs are now assumed to be in 1975$ also
	// We are assuming that secondary values will be in 1975$/kg
    double profitRate = ( price - mNonLandVariableCost + secondaryValue) * mYield; 

    // We multiply by 1e9 since profitRate above is in $/m2
    // and the land allocator needs it in $/billion m2. This assumes yield is in kg/m2
    profitRate *= 1e9;

    return profitRate;
}

/*! \brief Calculate the supply for the technology.
* \details Calculates the ag product produced which is equal to the yield multiplied
*          by the land allocated.
* \param aProductName Product name.
* \param aPeriod Period.
* \return The supply produced by the technology.
*/
double AgProductionTechnology::calcSupply( const string& aRegionName,
                                             const string& aProductName,
                                             const int aPeriod ) const
{
    double landAllocation = mProductLeaf->getLandAllocation( mName, aPeriod );

    // Set output to yield times amount of land.
    return mYield * landAllocation;
}

void AgProductionTechnology::doInterpolations( const Technology* aPrevTech, const Technology* aNextTech ) {
    Technology::doInterpolations( aPrevTech, aNextTech );

    const AgProductionTechnology* prevAgTech = static_cast<const AgProductionTechnology*> ( aPrevTech );
    const AgProductionTechnology* nextAgTech = static_cast<const AgProductionTechnology*> ( aNextTech );
    
    /*!
     * \pre We were given a valid previous ag production technology.
     */
    assert( prevAgTech );
    
    /*!
     * \pre We were given a valid next ag production technology.
     */
    assert( nextAgTech );
    
    // productivity change is a tech change and should be held constant at the next
    // technologies values to retain the same behavior
    mAgProdChange = nextAgTech->mAgProdChange;
    
    mNonLandCostTechChange = nextAgTech->mNonLandCostTechChange;
    
    // Non land variable costs will be interpoalted.
    mNonLandVariableCost = util::linearInterpolateY( year, prevAgTech->year, nextAgTech->year,
                                                     prevAgTech->mNonLandVariableCost,
                                                     nextAgTech->mNonLandVariableCost );
}

Value AgProductionTechnology::getParsedShareWeight() const {
    // Ag production technologies do not have shares and thus no share
    // weights.  We are going to return an intialized share-weight here
    // anyways so that consistency checks don't fail.
    const Value defaultShareWeight( 1.0 );
    return defaultShareWeight;
}
