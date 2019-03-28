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
* \file region_minicam.cpp
* \ingroup Objects
* \brief The RegionMiniCAM class source file.
* \author Sonny Kim
*/

//TODO: Clean up #includes

#include "util/base/include/definitions.h"
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <algorithm>
#include <memory>
#include <stack>

#include "containers/include/region_minicam.h"
#include "containers/include/gdp.h"
#include "containers/include/scenario.h"
#include "containers/include/info_factory.h"
#include "containers/include/iinfo.h"
#include "containers/include/dependency_finder.h"
#include "sectors/include/cal_quantity_tabulator.h"

// TODO: This needs a factory.
#include "sectors/include/sector.h"
#include "sectors/include/supply_sector.h"
#include "sectors/include/production_sector.h"
#include "sectors/include/demand_sector.h"
#include "sectors/include/ag_supply_sector.h"
#include "sectors/include/building_final_demand.h"
#include "sectors/include/export_sector.h"

#include "resources/include/resource.h"

#include "demographics/include/demographic.h"

#include "marketplace/include/marketplace.h"

#include "land_allocator/include/land_allocator.h"
#include "emissions/include/emissions_summer.h"
#include "emissions/include/luc_emissions_summer.h"
#include "policy/include/policy_ghg.h"
#include "emissions/include/total_sector_emissions.h"

#include "util/base/include/input_finder.h"
#include "util/base/include/summary.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "util/base/include/configuration.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"

#include "reporting/include/indirect_emissions_calculator.h"

#include "containers/include/sector_cycle_breaker.h"
#include "util/base/include/calibrate_share_weight_visitor.h"
#include "util/base/include/calibrate_resource_visitor.h"

using namespace std;
using namespace xercesc;

typedef std::vector<AFinalDemand*>::iterator FinalDemandIterator;
typedef std::vector<AFinalDemand*>::const_iterator CFinalDemandIterator;
typedef std::vector<AResource*>::iterator ResourceIterator;
typedef std::vector<GHGPolicy*>::iterator GHGPolicyIterator;
typedef std::vector<GHGPolicy*>::const_iterator CGHGPolicyIterator;
typedef std::vector<Sector*>::iterator SectorIterator;
typedef std::vector<Sector*>::reverse_iterator SectorReverseIterator;
typedef std::vector<Sector*>::const_iterator CSectorIterator;

extern Scenario* scenario;

//! Default constructor
RegionMiniCAM::RegionMiniCAM() {
    /*! \pre The modeltime object must be read-in before the Region can be
    *        parsed.
    */
    assert( scenario->getModeltime() );

    // Resize all vectors to maximum period
    const int maxper = scenario->getModeltime()->getmaxper();
    summary.resize( maxper );
    calibrationGDPs.resize( maxper );
    GDPcalPerCapita.resize( maxper );

    mRotationPeriod = 0;
    mInterestRate = 0;
}

//! Default destructor destroys sector, demsector, Resource, and
//! population objects.
RegionMiniCAM::~RegionMiniCAM() {
    clear();
}

//! Clear member variables and initialize elemental members.
void RegionMiniCAM::clear(){

    for ( FinalDemandIterator demIter = mFinalDemands.begin(); demIter != mFinalDemands.end(); ++demIter ) {
        delete *demIter;
    }
}

// for parsing derived region classes
bool RegionMiniCAM::XMLDerivedClassParse( const std::string& nodeName, const xercesc::DOMNode* curr ) {
    if( nodeName == "PrimaryFuelCO2Coef" ) {
        primaryFuelCO2Coef[ XMLHelper<string>::getAttr( curr, "name" ) ] = XMLHelper<double>::getValue( curr );
    }

    else if( nodeName == GDP::getXMLNameStatic() ){
        parseSingleNode( curr, gdp, new GDP );
    }
    else if( nodeName == "coolingDegreeDays" ){
        mCoolingDegreeDays = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "heatingDegreeDays" ) {
        mHeatingDegreeDays = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "rotationPeriod" ) {
        mRotationPeriod = XMLHelper<int>::getValue( curr );
    }
    else if( nodeName == "interest-rate" ){
        mInterestRate = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "cooling-fraction-of-year-active" ){
        mCoolingFractionOfYearActive = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "heating-fraction-of-year-active" ){
        mHeatingFractionOfYearActive = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == SupplySector::getXMLNameStatic() ){
        parseContainerNode( curr, supplySector, supplySectorNameMap, new SupplySector( name ) );
    }
    else if( nodeName == ExportSector::getXMLNameStatic() ){
        parseContainerNode( curr, supplySector, supplySectorNameMap, new ExportSector( name ) );
    }
    else if( nodeName == AgSupplySector::getXMLNameStatic() ) {
        parseContainerNode( curr, supplySector, supplySectorNameMap, new AgSupplySector( name ) );
    }
    else if( nodeName == EnergyFinalDemand::getXMLNameStatic() ){
        parseContainerNode( curr, mFinalDemands, new EnergyFinalDemand );
    }
    else if( nodeName == BuildingFinalDemand::getXMLNameStatic() ){
        parseContainerNode( curr, mFinalDemands, new BuildingFinalDemand );
    }
    else if( nodeName == TotalSectorEmissions::getXMLNameStatic() ){
        parseContainerNode( curr, mAggEmissionsCalculators, new TotalSectorEmissions );
    }
    else if ( nodeName == LandAllocator::getXMLNameStatic() ) {
        parseSingleNode( curr, mLandAllocator, new LandAllocator );
    }
    // regional economic data
    else if( nodeName == "calibrationdata" ){
        // get all child nodes.
        DOMNodeList* nodeListChild = curr->getChildNodes();
        // loop through the child nodes.
        string nodeNameChild;
        for( unsigned int j = 0; j < nodeListChild->getLength(); j++ ){
            DOMNode* currChild = nodeListChild->item( j );
            nodeNameChild = XMLHelper<string>::safeTranscode( currChild->getNodeName() );
            const Modeltime* modeltime = scenario->getModeltime();

            if( nodeNameChild == "#text" ) {
                continue;
            }
            else if( nodeNameChild == "GDPcal" ) { // TODO: MOVE TO GDP
                XMLHelper<double>::insertValueIntoVector( currChild, calibrationGDPs, modeltime );
            }
            // Per-capita value -- is converted to total GDP using population
            else if( nodeNameChild == "GDPcalPerCapita" ) { // TODO: MOVE TO GDP
                XMLHelper<double>::insertValueIntoVector( currChild, GDPcalPerCapita, modeltime );
            }
            else {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Unrecognized text string: " << nodeNameChild << " found while parsing region->calibrationdata." << endl;
            }
        }
    }
#if SORT_TESTING
    // A list representing the correct order in which to calculate the sectors.
    else if( nodeName == "SectorOrderList" ){
        // get all child nodes.
        string nodeNameChild;
        DOMNodeList* nodeListChild = curr->getChildNodes();
        // loop through the child nodes.
        for( unsigned int j = 0; j < nodeListChild->getLength(); j++ ){
            DOMNode* currChild = nodeListChild->item( j );
            nodeNameChild = XMLHelper<string>::safeTranscode( currChild->getNodeName() );

            if( nodeNameChild == "#text" ) {
                continue;
            }
            else if( nodeNameChild == "SectorName" ){
                sectorOrderList.push_back( XMLHelper<string>::getValue( currChild ) );
            }
            else {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Unrecognized text string: " << nodeNameChild << " found while parsing region->SectorOrderList." << endl;
            }
        }
    }
#endif
    else {
        return false;
    }
    return true;
}


/*! Complete the initialization. Get the size of vectors, initialize AGLU,
*   create all markets, call complete initialization
*  functions for nested objects, update the fuel map, and find simultaneities.
* \todo I think since there is one indirect ghg object for each sector, it might
*       be better in sector. This may require deriving supply sector.
*/
void RegionMiniCAM::completeInit() {
    Region::completeInit();

    // Region info has no parent Info.
    mRegionInfo.reset( InfoFactory::constructInfo( 0, name ) );
    mRegionInfo->setInteger( "rotationPeriod", mRotationPeriod );
    mRegionInfo->setDouble( "heating-fraction-of-year-active",
                            mHeatingFractionOfYearActive );
    mRegionInfo->setDouble( "cooling-fraction-of-year-active",
                            mCoolingFractionOfYearActive );
    mRegionInfo->setDouble( "heatingDegreeDays", mHeatingDegreeDays );
    mRegionInfo->setDouble( "coolingDegreeDays", mCoolingDegreeDays );

    if( mInterestRate == 0 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "No interest rate was read-in for region " << name << "." << endl;
    }

    // Add the interest rate to the region info.
    mRegionInfo->setDouble( "interest-rate", mInterestRate );

    // initialize demographic
    if( demographic.get() ){
        demographic->completeInit();
    }

    // Initialize the GDP
    if( gdp.get() ){
        gdp->initData( demographic.get() );
    }

    // Initialize the dependency finder.
    // This may need to move to init calc when markets change by period.
    // Note that the ICycleBreaker construction is on a seperate line because
    // C++ is given liberty when evaluating expressions across statements.
    SectorCycleBreaker sectorCycleBreaker( scenario->getMarketplace(), name );
    DependencyFinder depFinder( &sectorCycleBreaker );
    for( SectorIterator sectorIter = supplySector.begin(); sectorIter != supplySector.end(); ++sectorIter ) {
        ( *sectorIter )->completeInit( mRegionInfo.get(), &depFinder, mLandAllocator.get() );
    }

    if ( mLandAllocator.get() ) {
        mLandAllocator->completeInit( name, mRegionInfo.get() );
    }

    for( FinalDemandIterator demandSectorIter = mFinalDemands.begin();
        demandSectorIter != mFinalDemands.end(); ++demandSectorIter )
    {
        // Pass null in for the dependency finder so that the demand sector
        // technologies will not add their dependencies. This is because demand
        // sectors are not currently included in the ordering. Pass null for the land allocator
        // since demand sectors cannot allocate land.
        ( *demandSectorIter )->completeInit( name, mRegionInfo.get() );
    }

#if SORT_TESTING
    // Reorder sectors before sorting to test it.
    if( !sectorOrderList.empty() ){
        reorderSectors( sectorOrderList );
    }
#endif
    // Sort the sectors according to the order given by the dependency finder.

    // This needs to be at the end of completeInit so that all objects have
    // completed their calls to aDepFinder->addDependency. 
    depFinder.createOrdering();
    const vector<string> orderList = depFinder.getOrdering();
    
    // Check for an empty order list in which case the reordering shoud be
    // skipped. This occurs for SGM.
    if( !orderList.empty() ){
        util::reorderContainer<vector<Sector*>::iterator, Sector*>(
                                supplySector.begin(), supplySector.end(),
                                depFinder.getOrdering() );
    }
    else{
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Skipping sector reordering due to an empty order list."
                << "  This is bad." << endl;
    }

    setupCalibrationMarkets();
}

/*! \brief Initialize the GDP calibration markets.

* \todo Move GDP calibration parameters into GDP object
*/
void RegionMiniCAM::setupCalibrationMarkets() {
    if( !ensureGDP() ){
        return;
    }

    const Modeltime* modeltime = scenario->getModeltime();
    // Change per capita calibration values into absolute values
    for ( int period = 0; period < modeltime->getmaxper(); period++ ) {
        if ( GDPcalPerCapita [ period ] != 0 && demographic.get() ) {
            if ( calibrationGDPs[ period ] != 0 ) {
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::WARNING );
                mainLog << "Both GDPcal and GDPcalPerCapita read in region " << name << ". GDPcalPerCapita used." << endl;
            }
            // Convert from GDP/cap ($) to millions of dollars GDP total. Pop is in 1000's, so need an additional 1e3 scale to convert to millions
            // TODO: This will eventually write back out, so if this input file were used again the above warning
            // would be triggered.
            calibrationGDPs [ period ]  = GDPcalPerCapita [ period ] * demographic->getTotal( period )/1e3;
        }
    }

    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        gdp->setupCalibrationMarkets( name, calibrationGDPs );
    }
}

/*! \brief Reorder the sectors based on a list of sector names.
* \details This function is used to reorder a list of sectors using a supplied
*          ordering. In the partial-equilibrium model, the ordering of the
*          supply sectors is critical to the correct solving of the model. This
*          function handles errors as follows:
* <ul><li> If a sector is not specified in the list, the function will issue a
*          warning and remove the sector.</li>
*     <li> If a sector name in the ordering list is not an existing sector in
*          the model, a debugging warning will be issued and the sector will be
*          skipped. This can occur for resources and other non-sector
*          fuels.</li></ul>
* \param aOrderList A list of sector names in the order in which the sectors
*        should be put.
* \return Whether all sectors had orderings in the passed in order list.
*/
bool RegionMiniCAM::reorderSectors( const vector<string>& aOrderList ){
    // Check for an empty order list in which case the reordering should be
    // skipped. This occurs for SGM.
    if( aOrderList.empty() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Skipping sector reordering due to an empty order list. This is normal for SGM." << endl;
        return false;
    }

    // Create temporary copy of the sectors and the sector name map.
    vector<Sector*> originalOrder = supplySector;
    map<string,int> originalNameMap = supplySectorNameMap;

    // Clear the list of sectors and sectorNames.
    supplySector.clear();
    supplySectorNameMap.clear();

    // Get the dependency finder logger.
    ILogger& depFinderLog = ILogger::getLogger( "dependency_finder_log" );

    // Loop through the sector order vector.
    typedef vector<string>::const_iterator NameIterator;
    typedef map<string,int>::iterator NameMapIterator;

    for( NameIterator currSectorName = aOrderList.begin(); currSectorName != aOrderList.end(); ++currSectorName ){
        NameMapIterator origSectorPosition = originalNameMap.find( *currSectorName );

        // Check if the sector name in the sector ordering list exists
        // currently.
        if( origSectorPosition != originalNameMap.end() ){
            // Assign the sector.
            supplySector.push_back( originalOrder[ origSectorPosition->second ] );
            supplySectorNameMap[ *currSectorName ] = static_cast<int>( supplySector.size() - 1 );

            // Remove the original sector mapping. This will allow us to clean
            // up sectors that were not assigned a new ordering. This also is
            // more efficient as there are less entries to search.
            originalNameMap.erase( origSectorPosition );
        }
        else {
            depFinderLog.setLevel( ILogger::DEBUG );
            depFinderLog << *currSectorName << " is not the name of an existing sector. "
                << "It will not be included in the sector ordering." << endl;
        } // end else
    } // end for.

    // Check if there are any unassigned sectors and remove them.
    for( NameMapIterator currSecName = originalNameMap.begin(); currSecName != originalNameMap.end(); ++currSecName ){
        ILogger& mainLog = ILogger::getLogger( "main_log");
        mainLog.setLevel( ILogger::ERROR );
        mainLog << currSecName->first << " was not assigned a position in the explicit sector ordering list." << endl
            << "This sector will be removed from the model." << endl;

        depFinderLog.setLevel( ILogger::ERROR );
        depFinderLog << currSecName->first << " was not assigned a position in the explicit sector ordering list." << endl
            << "This sector will be removed from the model." << endl;


        // This sector is not in the new list, so free its memory.
        delete originalOrder[ currSecName->second ];
    }
    // This was successful if the sector name map is empty, so that all sectors
    // were ordered.

    // Print out sector ordering for debugging purposes
    depFinderLog.setLevel( ILogger::DEBUG );
    depFinderLog << "SECTOR ORDERING FOR REGION " << name << endl;
    depFinderLog << "  ";
    for( NameIterator currSectorName = aOrderList.begin(); currSectorName != aOrderList.end(); ++currSectorName ){
        depFinderLog << *currSectorName << ", ";
    }
    depFinderLog << endl;

    return originalNameMap.empty();
}

void RegionMiniCAM::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    // Write out the Co2 Coefficients.
    for( map<string,double>::const_iterator coefAllIter = primaryFuelCO2Coef.begin(); coefAllIter != primaryFuelCO2Coef.end(); coefAllIter++ ) {
        XMLWriteElement( coefAllIter->second, "PrimaryFuelCO2Coef", out, tabs, 0, coefAllIter->first );
    }

    XMLWriteElementCheckDefault( mHeatingDegreeDays, "heatingDegreeDays", out, tabs, Value( 0.0 ) );
    XMLWriteElementCheckDefault( mCoolingDegreeDays, "coolingDegreeDays", out, tabs, Value( 0.0 ) );
    XMLWriteElementCheckDefault( mCoolingFractionOfYearActive, "cooling-fraction-of-year-active", out, tabs, Value( 0.0 ) );
    XMLWriteElementCheckDefault( mHeatingFractionOfYearActive, "heating-fraction-of-year-active", out, tabs, Value( 0.0 ) );
    XMLWriteElementCheckDefault( mRotationPeriod, "rotationPeriod", out, tabs, 0 );
    XMLWriteElementCheckDefault( mInterestRate, "interest-rate", out, tabs, 0.0 );

    // write the xml for the class members.

    // write out data for land allocator
    if( mLandAllocator.get() ){
        mLandAllocator->toInputXML( out, tabs);
    }

    if( gdp.get() ){ // Check if gdp object exists
        gdp->toInputXML( out, tabs );
    }

    // write out demand sector objects.
    for( CFinalDemandIterator k = mFinalDemands.begin(); k != mFinalDemands.end(); k++ ){
        ( *k )->toInputXML( out, tabs );
    }

    for( unsigned int i = 0; i < mAggEmissionsCalculators.size(); ++i ){
        mAggEmissionsCalculators[ i ]->toInputXML( out, tabs );
    }

    // Note: The count function is an STL algorithm that counts the number of
    // times a value occurs within the a range of a container. The first two
    // arguments to the function are the range of the container to search, the
    // third is the value to search for.
    if( ( count( calibrationGDPs.begin(), calibrationGDPs.end(), 0 )
        != static_cast<int>( calibrationGDPs.size() ) ) ){ // makes sure tags aren't printed if no real data

            // Write out regional economic data
            XMLWriteOpeningTag( "calibrationdata", out, tabs );

            // write out calibration GDP
            const Modeltime* modeltime = scenario->getModeltime();
            XMLWriteVector( calibrationGDPs, "GDPcal", out, tabs, modeltime, 0.0 );

            XMLWriteClosingTag( "calibrationdata", out, tabs );
            // End write out regional economic data
        } // close calibration IF
}

void RegionMiniCAM::toDebugXMLDerived( const int period, std::ostream& out, Tabs* tabs ) const {
    // write out basic datamembers
    XMLWriteElement( mHeatingDegreeDays, "heatingDegreeDays", out, tabs );
    XMLWriteElement( mCoolingDegreeDays, "coolingDegreeDays", out, tabs );
    XMLWriteElement( mCoolingFractionOfYearActive, "cooling-fraction-of-year-active", out, tabs  );
    XMLWriteElement( mHeatingFractionOfYearActive, "heating-fraction-of-year-active", out, tabs );
    XMLWriteElement( mRotationPeriod, "rotationPeriod", out, tabs );
    XMLWriteElement( mInterestRate, "interest-rate", out, tabs );

    XMLWriteElement( calibrationGDPs[ period ], "calibrationGDPs", out, tabs );
    XMLWriteElement( getEndUseServicePrice( period ), "priceSer", out, tabs );
    
    // Write out the Co2 Coefficients.
    for( map<string,double>::const_iterator coefAllIter = primaryFuelCO2Coef.begin(); coefAllIter != primaryFuelCO2Coef.end(); coefAllIter++ ) {
        XMLWriteElement( coefAllIter->second, "PrimaryFuelCO2Coef", out, tabs, 0, coefAllIter->first );
    }

    // write the xml for the class members.
    // write out the single population object.
    if( gdp.get() ){
        gdp->toDebugXML( period, out, tabs );
    }

    // Write out the land allocator.
    if ( mLandAllocator.get() ) {
        mLandAllocator->toDebugXML( period, out, tabs );
    }

    // write out demand sector objects.
    for( CFinalDemandIterator currSector = mFinalDemands.begin(); currSector != mFinalDemands.end(); ++currSector ){
        (*currSector)->toDebugXML( period, out, tabs );
    }

    for( unsigned int i = 0; i < mAggEmissionsCalculators.size(); ++i ){
        mAggEmissionsCalculators[ i ]->toDebugXML( period, out, tabs );
    }
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. This function may be virtual to be overridden by derived class
* pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const std::string& RegionMiniCAM::getXMLName() const {
    return getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME. This way
* the tag is always consistent for both read-in and output and can be easily
* changed. The "==" operator that is used when parsing, required this second
* function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const std::string& RegionMiniCAM::getXMLNameStatic() {
    // TODO: Change to region-minicam, this will require changing input files
    static const string XML_NAME = "region";
    return XML_NAME;
}

//! Calculate the region.
void RegionMiniCAM::calc( const int period ) {
    // Store configuration variables locally as statics.
    static const Configuration* conf = Configuration::getInstance();
    static const bool calibrationActive = conf->getBool( "CalibrationActive" );

    // Write back calibrated values to the member variables.
    // These are still trial values.
    if( calibrationActive && gdp.get() ) {
        gdp->writeBackCalibratedValues( name, period );
    }
    // calculate regional GDP
    calcGDP( period );

    // determine prices for all supply sectors (e.g., refined fuels, electricity, etc.)
    calcFinalSupplyPrice( period );
    
    // adjust GDP for energy cost changes
    adjustGDP( period );

    // determine end-use demand for energy and other goods
    calcEndUseDemand( period );

    // determine supply of final energy and other goods based on demand
    setFinalSupply( period );

    // determine supply of primary resources
    calcResourceSupply( period );

    // Perform calibrations
    if( calibrationActive ) {
        calibrateRegion( period );
    }
}

/*! \brief Set regional ghg constraint from input data to market supply.
* \param period Model time period
*/
void RegionMiniCAM::calcResourceSupply( const int period ){
    if( !ensureGDP() ){
        return;
    }
    
    for( ResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
        (*currResource)->calcSupply( name, gdp.get(), period );
    }
}

/*!
 * \brief Calculate prices of all sectors.
 * \details If calibration is active we must have the CalibrateShareWeightVisitor
 *          calculate share weights at the current prices before calculating the
 *          price for that sector.  Doing this ensures consistent prices and share
 *          weights.
 * \todo Move this functionality into sector.
 * \param period Model time period
 */
void RegionMiniCAM::calcFinalSupplyPrice( const int period ) {
    if( !ensureGDP() ){
        return;
    }

    // determine if we should calibrate our share weights
    const bool calibrationPeriod = period > 0 && period <= scenario->getModeltime()->getFinalCalibrationPeriod();
    static const bool calibrationActive = Configuration::getInstance()->getBool( "CalibrationActive" );

    CalibrateShareWeightVisitor calibrator( name, gdp.get() );
    for( SectorIterator currSupply = supplySector.begin(); currSupply != supplySector.end(); ++currSupply ){
        // if calibration is active and we are in a calibration period we must calibrate our share weights
        // before calculating prices
        if( calibrationActive && calibrationPeriod ) {
                (*currSupply)->accept( &calibrator, period );
        }
        (*currSupply)->calcFinalSupplyPrice( gdp.get(), period );
    }

    if( calibrationActive && calibrationPeriod ) {
        CalibrateResourceVisitor resourceCalibrator( name );
        for( ResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
            (*currResource)->accept( &resourceCalibrator, period );
        }
    }

    if ( mLandAllocator.get() ) {
        mLandAllocator->calcFinalLandAllocation( name, period );
    }
}

/*! Calculates supply of final energy and other goods.
* \todo Move functionality into sector
* \param period Model time period
*/
void RegionMiniCAM::setFinalSupply( const int period ) {
    if( !ensureGDP() ){
        return;
    }

    // loop through all sectors in reverse once to get total output.
    typedef vector<Sector*>::reverse_iterator ReverseSectorIterator;
    for ( ReverseSectorIterator currSupply = supplySector.rbegin(); currSupply != supplySector.rend(); ++currSupply ) {
        (*currSupply)->supply( gdp.get(), period );
    }
}

/*! Calculate initial gdp value (without feedbacks)
*
* \param period Model time period
*/
void RegionMiniCAM::calcGDP( const int period ){
    if( !ensureGDP() || !ensureDemographics() ){
        return;
    }
    gdp->initialGDPcalc( period, demographic->getTotal( period ) );
}

/*! \brief Calculate forward-looking gdp (without feedbacks) for AgLU use.
* \details It is necessary to have a gdp without feedbacks so that all values
*          are known and AgLU can calibrate. This routine runs through each
*          period and calculates a series of gdp values without use of the
*          energy price feedback.
* \author Steve Smith, Josh Lurz
* \warning This will interfere with the normal gdp calculation if this is used
*          after model calc starts.
* \todo check to see if this works with AgLU. Not sure about conversions.
*/
const vector<double> RegionMiniCAM::calcFutureGDP() const {
    if( !ensureGDP() || !ensureDemographics() ){
        return vector<double>( 0 );
    }

    const Modeltime* modeltime = scenario->getModeltime();
    vector<double> gdps( modeltime->getmaxper() );

    for ( int period = 0; period < modeltime->getmaxper(); period++ ) {
        gdp->initialGDPcalc( period, demographic->getTotal( period ) );
        gdps[ period ] = gdp->getApproxScaledGDPperCap( period );
    }
    return gdps;
}

double RegionMiniCAM::getEndUseServicePrice( const int period ) const {
    double servicePrice = 0;
    for ( CFinalDemandIterator currDemSector = mFinalDemands.begin(); currDemSector != mFinalDemands.end(); ++currDemSector ) {
        servicePrice += (*currDemSector)->getWeightedEnergyPrice( name, period );
    }
    return servicePrice;
}

/*! Adjust regional gdp for energy.
*
* \param period Model time period
* \todo Move this calculation down to GDP
*/
void RegionMiniCAM::adjustGDP( const int period ){
    if( !ensureGDP() ){
        return;
    }

    const Modeltime* modeltime = scenario->getModeltime();

    double tempratio = 1;
    if ( period > modeltime->getFinalCalibrationPeriod() ){
        tempratio = getEndUseServicePrice( period ) / getEndUseServicePrice( period - 1 );
    }
    gdp->adjustGDP( period, tempratio );
}

/*! \brief Perform regional calibration of GDP.
* \details This method handles calibration of the GDP.  For Calibrations
*          of subsectors and technologies see RegionMiniCAM::calcFinalSupplyPrice
* \param period Model time period
*/
void RegionMiniCAM::calibrateRegion( const int period ) {
    if( !ensureGDP() ){
        return;
    }
    // TODO: Move this into GDP object.
    if( calibrationGDPs[ period ] > 0 ){
        const string goodName = "GDP";
        Marketplace* marketplace = scenario->getMarketplace();
        marketplace->addToSupply( goodName, name, gdp->getGDP( period ), period );
    }
}

/*! \brief Test to see if calibration worked for all sectors in this region
*
* Compares the sum of calibrated + fixed values to output of each sector.
*
* If calibrations are not on then will only printout out diagnostics (if
*
* \author Steve Smith
* \param period Model period
* \param calAccuracy value to which calibrations must match in order to pass
*        calibration test.
* \param printWarnings flag to turn on logging of warnings if calibrations are
*        not accurate.
* \return Boolean true if calibration is ok.
*/
bool RegionMiniCAM::isAllCalibrated( const int period, double calAccuracy, const bool printWarnings ) const {
    const static bool calOn = Configuration::getInstance()->getBool( "CalibrationActive" );
    
    // Don't check calibration in the base period or if calibration is off.
    if( !calOn || period == 0 ){
        return true;
    }

    bool returnVal = true;
    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        if ( !supplySector[ i ]->isAllCalibrated( period, calAccuracy, printWarnings ) ) {
            returnVal = false;
        }
    }

    // always return true if calibrations are not on.
    return returnVal;
}

/*! \brief Call any initializations that are only done once per period.
* \author Steve Smith, Josh Lurz, Sonny Kim
* \param period Model period
* \todo Once postCalc is present, add a check to verify that aggregate emissions worked properly
*/
void RegionMiniCAM::initCalc( const int period )
{

    // Set aggregate emissions factor to region info if needed
    for ( unsigned int i = 0; i < mAggEmissionsCalculators.size(); i++ ) {
        mAggEmissionsCalculators[ i ]->setAggregateEmissionFactor( name,
                                                                   supplySector,
                                                                   mRegionInfo.get(),
                                                                   period );
    }

    mRegionInfo->setDouble( "heatingDegreeDays", mHeatingDegreeDays );
    mRegionInfo->setDouble( "coolingDegreeDays", mCoolingDegreeDays );

    // Set the CO2 coefficients into the Marketplace before the Technologies and
    // GHGs are initialized so they can be accessed.
    setCO2CoefsIntoMarketplace( period );


    for( SectorIterator currSector = supplySector.begin(); currSector != supplySector.end(); ++currSector ){
        (*currSector)->initCalc( /*nationalAccount[ period ]*/0, demographic.get(), period );
    }
    for ( FinalDemandIterator currSector = mFinalDemands.begin(); currSector != mFinalDemands.end(); ++currSector ) {
        (*currSector)->initCalc( name, gdp.get(), demographic.get(), period  );
    }
    for( ResourceIterator currResource = mResources.begin(); currResource != mResources.end(); ++currResource ){
        (*currResource)->initCalc( name, period );
    }

    // Call initCalc for land allocator last. It needs profit from the ag sectors
    // before it can calculate share weights
    if ( mLandAllocator.get() ) {
        mLandAllocator->initCalc( name, period );
    }

}

/*
* \brief Initialize the CO2 coefficients read in by the Region into the
*        marketplace.
* \details In each period the Region must set the CO2 coefficients for all goods
*          into the Marketplace, so that Technologies and GHGs can access them.
* \param aPeriod Period.
*/
void RegionMiniCAM::setCO2CoefsIntoMarketplace( const int aPeriod ){
    const static string CO2COEF = "CO2coefficient";
    Marketplace* marketplace = scenario->getMarketplace();
    for( map<string, double>::const_iterator coef = primaryFuelCO2Coef.begin();
        coef != primaryFuelCO2Coef.end(); ++coef )
    {
        // Markets may not exist for incorrect fuel names.
        IInfo* fuelInfo = marketplace->getMarketInfo( coef->first, name, aPeriod, false );
        if( fuelInfo ){
            fuelInfo->setDouble( CO2COEF, coef->second );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::DEBUG );
            mainLog << "Cannot set emissions factor of zero for fuel " << coef->first
                    << " because the name does not match the name of a market." << endl;
        }
    }
}


/*!
 * \brief Call any calculations that are only done once per period after
 *        solution is found.
 * \author Sonny Kim
 * \param aPeriod Model period
 */
void RegionMiniCAM::postCalc( const int aPeriod ) {
    Region::postCalc( aPeriod );
    
    mLandAllocator->postCalc( name, aPeriod );
 }

/*! \brief Tabulates all calibrated supplies and demands with in the region.
*
* The result of this routine, once called for all regions, is a set of calSupply
* and calDemand values in marketInfo. The calSupply variable has a non-negative
* value for each supply sector that is completely calibrated. Similarly, the
* calDemand variable has a non-negative value for each demand good that is
* completely calibrated. Only those markets where all supplies or demands are
* calibrated will have non-zero values for these variables.
* Note that currently only the detailed buildings rely on these calibration values
* being set.  See the calculation in
* CalibrateShareWeightVisitor::startVisitBuildingGenericDmdTechnology.  We could
* remove this code if we just read in the proper input coefficients into
* BuildingDemandInput however that tends to be a bit of work due to all the 
* internal gains to keep track of.
*
* \author Steve Smith
* \param period Model period
*/

void RegionMiniCAM::setCalSuppliesAndDemands( const int period ) {
    // Check for fully calibrated/fixed supplies. This will search sectors and
    // resources for calibrated and fixed supples.
    CalQuantityTabulator calSupplyTabulator( "" );
    accept( &calSupplyTabulator, period );

    Marketplace* marketplace = scenario->getMarketplace();

    // Flag that the market is not all fixed.
    const double MARKET_NOT_ALL_FIXED = -1;

    // Name of the cal supply info string.
    const string CAL_SUPPLY_NAME = "calSupply";

    // Iterate over the entire map of goods.
    const CalQuantityTabulator::CalInfoMap& calSupplies = calSupplyTabulator.getSupplyInfo();

    for( CalQuantityTabulator::CalInfoMap::const_iterator iter = calSupplies.begin();
         iter != calSupplies.end(); ++iter )
    {
        // Get the market info for the good.
        // TODO: Require market info to exist once demand sectors are separated.
        IInfo* marketInfo = marketplace->getMarketInfo( iter->first, name, period, false );

        // Incorrect input files may cause the market to not exist. This will be
        // flagged elsewhere.
        if( !marketInfo ){
            continue;
        }

        double existingCalSupply = marketInfo->getDouble( CAL_SUPPLY_NAME, false );
        // Check is to see if some other region has flagged this as having all
        // inputs not fixed.
        if ( existingCalSupply != MARKET_NOT_ALL_FIXED && iter->second.mAllFixed ) {
            // If supply of this good has not been elimiated from the search and
            // output is fixed then add to fixed supply value.
            marketInfo->setDouble( CAL_SUPPLY_NAME, max( existingCalSupply, 0.0 )
                                   + iter->second.mCalQuantity + iter->second.mFixedQuantity );
        }
        else {
            // If supply of this good is not fixed then set flag to eliminate
            // this from other searches
            marketInfo->setDouble( CAL_SUPPLY_NAME, MARKET_NOT_ALL_FIXED );
        }
    }

    // Check for fully calibrated/fixed demands. Search through all demand
    // technologies
    for ( unsigned int i = 0; i < mFinalDemands.size(); i++ ) {
        mFinalDemands[ i ]->tabulateFixedDemands(name, demographic.get(), gdp.get(), period );
    }
}

/*! \brief Initializes calDemand and calSupply values for all sectors
*
* All calSupply and calDemand values are reset to -1 This erases any information
* from a prior period and eliminates warning messages for value not found in
* information store
*
* \author Steve Smith
* \param period Model period
*/
void RegionMiniCAM::initializeCalValues( const int period ) {
    Marketplace* marketplace = scenario->getMarketplace();
    // This will get set to -1 if this market is not fixed and a value >= 0
    // otherwise.
    const double DEFAULT_CAL_VALUE = -2;

    typedef list<string>::const_iterator CI;
    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        // Set default values for supply sector
        IInfo* marketInfo = marketplace->getMarketInfo( supplySector[ i ]->getName(), name, period, true );
        marketInfo->setDouble( "calSupply", DEFAULT_CAL_VALUE );
        marketInfo->setDouble( "calDemand", DEFAULT_CAL_VALUE );
    }

    for ( unsigned int i = 0; i < mResources.size(); i++ ) {
        // Now prepare to loop through all the fuels used by this sector
        IInfo* resourceMarketInfo = marketplace->getMarketInfo( mResources[ i ]->getName(), name, period, true );
        resourceMarketInfo->setDouble( "calSupply", DEFAULT_CAL_VALUE );
        resourceMarketInfo->setDouble( "calDemand", DEFAULT_CAL_VALUE );
    }
}

//! Calculate regional demand for energy and other goods for all sectors.
void RegionMiniCAM::calcEndUseDemand( const int period ) {
    if( !ensureGDP() ){
        return;
    }

    for ( FinalDemandIterator currDemSector = mFinalDemands.begin(); currDemSector != mFinalDemands.end(); ++currDemSector ){
        // calculate aggregate demand for end-use sector services
        // set fuel demand from aggregate demand for services
        (*currDemSector)->setFinalDemand( name, demographic.get(), gdp.get(), period );
    }

}

//! Calculate regional emissions from resources.
void RegionMiniCAM::calcEmissions( const int period ) {
    summary[period].clearemiss(); // clear emissions map

    // need to call emissions function but sum is not needed
    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        supplySector[i]->emission(period);
        summary[period].updateemiss(supplySector[i]->getemission(period));
    }

    // Determine land use change emissions and add to the summary.
	if( period > 0 ) {
		LUCEmissionsSummer co2LandUseSummer( "CO2NetLandUse" );


		const int year = scenario->getModeltime()->getper_to_yr( period );
		accept( &co2LandUseSummer, period );
		map<string,double> agEmissions;
		if ( co2LandUseSummer.areEmissionsSet( year ) ) {
			agEmissions[ "CO2NetLandUse" ] = co2LandUseSummer.getEmissions( year );
			summary[ period ].updateemiss( agEmissions );
		}
	}
}

/*! \brief Calculate regional emissions by fuel for reporting.
* \warning This function assumes emission has already been called, as this
*          function cannot clear the summary emissions.
* \param The global list of primary fuels.
* \param period The model period.
*/
void RegionMiniCAM::calcEmissFuel( const list<string>& aPrimaryFuelList, const int period )
{
    map<string, double> fuelemiss; // tempory emissions by fuel

    for( list<string>::const_iterator fuelIter = aPrimaryFuelList.begin();
        fuelIter != aPrimaryFuelList.end(); ++fuelIter )
    {
        fuelemiss[ *fuelIter ] = summary[period].get_pemap_second( *fuelIter )
                                 * util::searchForValue( primaryFuelCO2Coef, *fuelIter );
    }

    summary[period].updateemiss(fuelemiss); // add CO2 emissions by fuel
}

//! Write all outputs to file.
void RegionMiniCAM::csvOutputFile() const {
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    vector<double> temp(maxper);
    // function protocol
    void fileoutput3(string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    // write population results to database
    if( demographic.get() ){
        demographic->csvOutputFile( name );
    }
    if( gdp.get() ){
        gdp->csvOutputFile( name );
    }

    // write total emissions for region
    for (int m= 0;m<maxper;m++) {
        temp[m] = summary[m].get_emissmap_second("CO2");
    }
    fileoutput3(name," "," "," ","CO2 emiss","MTC",temp);

    // write ag emissions for region
    for (int m= 0;m<maxper;m++) {
        temp[m] = summary[m].get_emissmap_second( "CO2NetLandUse" );
    }
    fileoutput3(name," "," "," ","Net Land Use CO2 emiss","MTC",temp);

    // region primary energy consumption by fuel type
    typedef map<string,double>:: const_iterator CI;
    map<string,double> tpemap = summary[0].getpecons();
    for (CI pmap=tpemap.begin(); pmap!=tpemap.end(); ++pmap) {
        for (int m= 0;m<maxper;m++) {
            temp[m] = summary[m].get_pemap_second(pmap->first);
        }
        fileoutput3(name,"Pri Energy","Consumption","",pmap->first,"EJ",temp);
    }


    // regional Pri Energy Production Total
    for (int m= 0;m<maxper;m++) {
        temp[m] = summary[m].get_peprodmap_second("zTotal");
    }
    fileoutput3(name,"Pri Energy","total","","zTotal","EJ",temp);

    // Calculate indirect emissions so they can be written to the database.
    IndirectEmissionsCalculator indEmissCalc;
    for( int i = 0; i < modeltime->getmaxper(); ++i ){
        accept( &indEmissCalc, i );
    }

    // write resource results to file
    for ( unsigned int i = 0; i < mResources.size(); i++ )
        mResources[i]->csvOutputFile( name );

    // write supply sector results to file
    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        supplySector[i]->csvOutputFile( gdp.get(), &indEmissCalc );
    }

    // write end-use sector demand results to file
    for ( unsigned int i = 0; i < mFinalDemands.size(); i++ ){
        mFinalDemands[i]->csvOutputFile( name );
    }
}

//! Write MiniCAM style outputs to file.
void RegionMiniCAM::dbOutput( const list<string>& aPrimaryFuelList ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    vector<double> temp(maxper),temptot(maxper);
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // write population results to database
    if( demographic.get() ){
        demographic->dbOutput( name );
    }
    if( gdp.get() ){
        gdp->dbOutput( name );
    }

    // CO2 emissions by fuel
    for( list<string>::const_iterator fuelIter = aPrimaryFuelList.begin();
        fuelIter != aPrimaryFuelList.end(); fuelIter++ )
    {
        for ( int m= 0;m<maxper;m++) {
            temp[m] = summary[m].get_emissfuelmap_second( *fuelIter );
            temptot[m] += temp[m];
        }
        dboutput4(name,"CO2 Emiss","by Fuel",*fuelIter,"MTC",temp);
    }
    // add amount of geologic sequestration to emissions by fuel
    // todo change hard-coded category name
    for ( int m= 0;m<maxper;m++) {
        // note the negative value for sequestered amount
        temp[m] = - summary[m].get_emissmap_second( "CO2sequestGeologic" );
        temptot[m] += temp[m];
    }
    dboutput4(name,"CO2 Emiss","by Fuel","geologic sequestration","MTC",temp);

    // add amount of sequestration from non-energy use to emissions by fuel
    // todo change hard-coded category name
    for ( int m= 0;m<maxper;m++) {
        // note the negative value for sequestered amount
        temp[m] = - summary[m].get_emissmap_second( "CO2sequestNonEngy" );
        temptot[m] += temp[m];
    }
    dboutput4(name,"CO2 Emiss","by Fuel","non-energy use","MTC",temp);

    // total emissions by sector for region
    for ( int m= 0;m<maxper;m++) {
        temp[m] = summary[m].get_emissmap_second("CO2");
    }
    // CO2 emissions by fuel and sector totals use same value
    dboutput4(name,"CO2 Emiss","by Fuel","zTotal","MTC",temptot);
    dboutput4(name,"CO2 Emiss","by Sector","zTotal","MTC",temp);

    // total ag sector emissions by sector for region
    for ( int m= 0;m<maxper;m++) {
        temp[m] = summary[m].get_emissmap_second( "CO2NetLandUse" );
    }
    dboutput4(name, "CO2 Emiss", "Net Land Use", "total", "MTC", temp );

    // regional emissions for all greenhouse gases
    typedef map<string,double>:: const_iterator CI;
    map<string,double> temissmap = summary[0].getemission(); // get gases for period 0
    for (CI gmap=temissmap.begin(); gmap!=temissmap.end(); ++gmap) {
        for ( int m= 0;m<maxper;m++) {
            temp[m] = summary[m].get_emissmap_second(gmap->first);
        }
        dboutput4(name,"Emissions","by gas",gmap->first,"MTC",temp);
    }

    // regional fuel consumption (primary and secondary) by fuel type
    map<string,double> tfuelmap = summary[0].getfuelcons();
    for (CI fmap=tfuelmap.begin(); fmap!=tfuelmap.end(); ++fmap) {
        for ( int m= 0;m<maxper;m++) {
            temp[m] = summary[m].get_fmap_second(fmap->first);
        }
        if( fmap->first == "" ){
            dboutput4(name,"Fuel Consumption"," Region by fuel","No Fuelname","EJ",temp);
        }
        else {
            dboutput4(name,"Fuel Consumption"," Region by fuel",fmap->first,"EJ",temp);
        }
    }

    // region primary energy consumption by fuel type
    map<string,double> tpemap = summary[0].getpecons();
    for (CI pmap = tpemap.begin(); pmap != tpemap.end(); ++pmap) {
        for ( int m = 0; m < maxper; ++m ) {
            temp[m] = summary[m].get_pemap_second(pmap->first);
        }
        dboutput4(name,"Primary Energy","Consumption by fuel(renew error)",pmap->first,"EJ",temp);
    }

    // regional Pri Energy Production by fuel type
    tpemap = summary[0].getpeprod();
    for (CI pmap = tpemap.begin(); pmap != tpemap.end(); ++pmap) {
        for ( int m = 0; m < maxper; ++m ) {
            temp[m] = summary[m].get_peprodmap_second(pmap->first);
        }
        dboutput4(name,"Primary Energy","Production by fuel(renew error)",pmap->first,"EJ",temp);
    }

    // region primary energy trade by fuel type
    tpemap = summary[0].getpetrade();
    for (CI pmap = tpemap.begin(); pmap != tpemap.end(); ++pmap) {
        for ( int m = 0; m < maxper; ++m ) {
            temp[m] = summary[m].get_petrmap_second(pmap->first);
        }
        dboutput4(name,"Primary Energy","Trade by fuel",pmap->first,"EJ",temp);
    }

    // write resource results to database
    for ( unsigned int i = 0; i < mResources.size(); i++ ) {
        mResources[i]->dbOutput( name );
    }

    // Calculate indirect emissions so they can be written to the database.
        IndirectEmissionsCalculator indEmissCalc;
        for( int i = 0; i < modeltime->getmaxper(); ++i ){
            accept( &indEmissCalc, i );
    }

    // write supply sector results to database
    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        supplySector[i]->dbOutput( gdp.get(), &indEmissCalc );
    }
    // write end-use sector demand results to database
    for ( unsigned int i = 0; i < mFinalDemands.size(); i++ ) {
        mFinalDemands[i]->dbOutput( name );
    }
}

//! update regional summaries for reporting
void RegionMiniCAM::updateSummary( const list<string>& aPrimaryFuelList, const int period ) {
    calcEmissions( period );

    summary[period].clearpeprod();
    summary[period].clearfuelcons();
    summary[period].clearemfuelmap();

    for ( unsigned int i = 0; i < mResources.size(); i++ ) {
        summary[period].initpeprod( aPrimaryFuelList, mResources[i]->getName(),
            mResources[i]->getAnnualProd( name, period) );
    }

    for ( unsigned int i = 0; i < supplySector.size(); i++ ) {
        // call update for supply sector
        supplySector[i]->updateSummary( aPrimaryFuelList, period );
        // update regional fuel consumption (primary and secondary) for supply sector
        summary[period].updatefuelcons( aPrimaryFuelList, supplySector[i]->getfuelcons(period));
        summary[ period ].updateemfuelmap( supplySector[ i ]->getemfuelmap( period ) );
    }
    // update primary energy trade from consumption and production amounts
    summary[period].updatepetrade();

    calcEmissFuel( aPrimaryFuelList, period );
}

//! Return the summary object for the given period.
/*! \todo This is a temporary fix to get the global CO2. This should be
*         restructured.
* \param period Model period to return the summary for.
* \return The summary object.
*/
const Summary& RegionMiniCAM::getSummary( const int period ) const {
    return summary[ period ];
}

//! update regional output tables for reporting
void RegionMiniCAM::updateAllOutputContainers( const int period ) {
}

/*! \brief Check whether the GDP object exists and report a warning if it does not.
* \return Whether the GDP object exists.
* \author Josh Lurz
*/
bool RegionMiniCAM::ensureGDP() const {
    if( !gdp.get() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "GDP object has not been created and is required. "
            << "Check for a region name mismatch." << endl;
        return false;
    }
    return true;
}

/*! \brief Check whether the Demographics object exists and report a warning if it does not.
* \return Whether the Demographics object exists.
* \author Josh Lurz
*/
bool RegionMiniCAM::ensureDemographics() const {
    if( !demographic.get() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "Population object has not been created and is required. "
            << "Check for a region name mismatch." << endl;
        return false;
    }
    return true;
}

/*! \brief Update a visitor for a Region.
* \param aVisitor Visitor to update.
* \param aPeriod Period to update.
*/
void RegionMiniCAM::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitRegionMiniCAM( this, aPeriod );
    Region::accept( aVisitor, aPeriod );

    // Visit GDP object
    if ( gdp.get() ){
        gdp->accept( aVisitor, aPeriod );
    }

    // Visit LandAllocator object
    if ( mLandAllocator.get() ){
        mLandAllocator->accept( aVisitor, aPeriod );
    }

    // loop for final demand sectors.
    for( CFinalDemandIterator currDem = mFinalDemands.begin(); currDem != mFinalDemands.end(); ++currDem ){
        (*currDem)->accept( aVisitor, aPeriod );
    }

    aVisitor->endVisitRegionMiniCAM( this, aPeriod );
}
