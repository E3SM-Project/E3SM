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
* \file world.cpp
* \ingroup Objects
* \brief world class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <string>
#include <cassert>
#include <vector>
#include <map>
#include <algorithm>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/xml_helper.h"
#include "containers/include/world.h"
#include "containers/include/region_minicam.h"
#include "containers/include/region_cge.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/configuration.h"
#include "util/base/include/util.h"
#include "util/base/include/summary.h"
#include "util/curves/include/curve.h"
#include "util/curves/include/point_set_curve.h"
#include "util/curves/include/point_set.h"
#include "util/curves/include/explicit_point_set.h"
#include "util/curves/include/xy_data_point.h"
#include "solution/util/include/calc_counter.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/ivisitor.h"
#include "climate/include/iclimate_model.h"
// Could hide with a factory method.
#include "climate/include/magicc_model.h"
#include "util/base/include/hash_map.h"
#include "util/base/include/atom_registry.h"
#include "emissions/include/emissions_summer.h"
#include "emissions/include/luc_emissions_summer.h"
#include "technologies/include/global_technology_database.h"
#include "reporting/include/energy_balance_table.h"

#include "reporting/include/rcp_emissions_visitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor.
World::World():
mCalcCounter( new CalcCounter() )
{
}

//! World destructor. 
World::~World(){
    clear();
}

//! Helper member function for the destructor. Performs memory deallocation. 
void World::clear(){
    for ( RegionIterator regionIter = regions.begin(); regionIter != regions.end(); regionIter++ ) {
        delete *regionIter;
    }
}

//! parses World xml object
void World::XMLParse( const DOMNode* node ){
    // assume we are passed a valid node.
    assert( node );

    // get all the children.
    DOMNodeList* nodeList = node->getChildNodes();

    for( unsigned int i = 0;  i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == GlobalTechnologyDatabase::getXMLNameStatic() ) {
            GlobalTechnologyDatabase::getInstance()->XMLParse( curr );
        }
        // MiniCAM regions
        else if( nodeName == RegionMiniCAM::getXMLNameStatic() ){
            parseContainerNode( curr, regions, regionNamesToNumbers, new RegionMiniCAM() );
        }
		// Read in parameters for climate model
        else if( nodeName == MagiccModel::getXMLNameStatic() ){
            parseSingleNode( curr, mClimateModel, new MagiccModel( scenario->getModeltime() ) );
        }
		// SGM regions
        else if( nodeName == RegionCGE::getXMLNameStatic() ){
            parseContainerNode( curr, regions, regionNamesToNumbers, new RegionCGE() );
        }
        else if( nodeName == "set-carbon-density" ) {
            mAdjDensityVisitor.XMLParse( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing World." << endl;
        }
    }
}

/*! \brief Complete the initialization
*
* This routine is only called once per model run
*
* \author Josh Lurz
*/
void World::completeInit() {

    // Initialize the region lookup hashmap.
    createFastLookupMap();

    // Finish initializing all the regions.
    for( RegionIterator regionIter = regions.begin(); regionIter != regions.end(); regionIter++ ) {
        ( *regionIter )->completeInit();
    }

    //If none has been read in, instantiate the default climate model
    if ( !mClimateModel.get() ) {
        mClimateModel.reset( new MagiccModel( scenario->getModeltime() ) );
    }
    
    // Initialize Climate Model
    mClimateModel->completeInit( scenario->getName() );

    // Push carbon density adjustments into the model
    //    mAdjDensityVisitor.doInterpolations();
    accept( &mAdjDensityVisitor, -1 );
}

/*! \brief Initialize the region partial derivative calculation hash map.
* \details Constructs a mapping of region Atom to index within the region
*          vector. This allows for rapid lookups of region numbers during
*          derivative calculations where not all regions are calculated at once.
*          This must be done before markets are created since the markets will
*          store the atoms of the regions that they contain.
*/
void World::createFastLookupMap(){
    // Construct the hashmap as 20 percent full for performance.
    mRegionLookupMap.reset( new FastRegionLookupMap( regions.size() * 5 ) );

    // Add each region id to number mapping to the hashmap.
    for( unsigned int i = 0; i < regions.size(); ++i ){
        // An atom does not need to be created if this is the second or later run
        // of a batch of scenarios.
        const objects::Atom* regionID = objects::AtomRegistry::getInstance()->findAtom( regions[ i ]->getName() );
        
        if( !regionID ){
            // Construct an atom for the region name. The Atom will automatically be
            // registered.
            regionID = new objects::Atom( regions[ i ]->getName() );
        }

        // Add an entry to the hashmap for the id.
        mRegionLookupMap->insert( make_pair( regionID, i ) );
    }
}

//! Write out datamembers to XML output stream.
void World::toInputXML( ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag ( getXMLNameStatic(), out, tabs );

    GlobalTechnologyDatabase::getInstance()->toInputXML( out, tabs );
	
    // Climate model parameters
    // note that due to a dependency in the carbon cycle model this
    // must be written out before any of the carbon cycle historical
    // year data is written which is contained in the regions
    if ( mClimateModel.get() ) {
        mClimateModel->toInputXML( out, tabs );
    }

    mAdjDensityVisitor.toInputXML( out, tabs );

    for( CRegionIterator i = regions.begin(); i != regions.end(); i++ ){
        ( *i )->toInputXML( out, tabs );
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

//! Write out XML for debugging purposes.
/*! \warning This only call Region::toInputXML for the US. */
void World::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag ( getXMLNameStatic(), out, tabs, "", period );

    // write the xml for the class members.

    scenario->getMarketplace()->toDebugXML( period, out, tabs );

    // Only print debug XML information for the specified region to avoid
    // unmanagably large XML files.
    const static string debugRegion = Configuration::getInstance()->getString( "debug-region", "USA" );
    for( CRegionIterator i = regions.begin(); i != regions.end(); i++ ) {
        if( ( *i )->getName() == debugRegion ){
            ( *i )->toDebugXML( period, out, tabs );
        }
    }

    // Climate model parameters
    if ( !mClimateModel.get() ) {
        mClimateModel->toDebugXML( period, out, tabs );
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
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
const std::string& World::getXMLNameStatic() {
    const static string XML_NAME = "world";
    return XML_NAME;
}

/*! \brief Returns the name of the World.
* \details Although there is only one World in the model so a name is
*          unnecessary, this is needed for compatibility with the IParsable
*          interface. The function returns instead the XML name.
* \note In the future World object may have read-in names.
* \return The name of the world.
*/
const string& World::getName() const {
    return getXMLNameStatic();
}

//! initialize anything that won't change during the calculation
/*! Examples: share weight scaling due to previous calibration, 
* cumulative technology change, etc.
*/
void World::initCalc( const int period ) {

    for( vector<Region*>::iterator i = regions.begin(); i != regions.end(); i++ ){
        // Add supplies and demands to the marketplace in the base year for checking data consistency
        // and for getting demand and supply totals.
        // Need to update markets here after markets have been null by scenario.
        // TODO: This should be combined with check data.
        if( period == 0 ){
            ( *i )->updateMarketplace( period );
        }
        ( *i )->initCalc( period );
    }
    
    Configuration* conf = Configuration::getInstance();
    if( conf->getBool( "CalibrationActive" ) ){
        // print an I/O table for debuging before we do any calibration
        ILogger& calLog = ILogger::getLogger( "calibration_log" );
        calLog.setLevel( ILogger::DEBUG );
        for( CRegionIterator reigonIt = regions.begin(); reigonIt != regions.end(); ++reigonIt ){
            // for this table we will want a condensed table without non-calibrated values
            // so the user can get an easy to see view of what they put in
            EnergyBalanceTable table( (*reigonIt)->getName(), calLog, true, false );
            (*reigonIt)->accept( &table, period );
            table.finish();
        }
        
        // TODO: could rename although we could also just get rid of it
        // see the comments for RegionMiniCAM::setCalSuppliesAndDemands
        checkCalConsistancy( period );
    }
    
    // Reset the calc counter.
    mCalcCounter->startNewPeriod();
}

/*! \brief Tabulates calibrated supplies and demands and stores them in the
 *         corresponding market info.
 * \details Has each region set their calibrated supplies and demands.  This method
 *          no longer checks for consistency however.
 * \param period The current model period.
 * \author Steve Smith
 * \todo Shouldn't this all be at the regional level?
 * \todo This method could be eliminated if it were not for detailed buildings.
 * \return It will always return false currently.
*/
bool World::checkCalConsistancy( const int period ) {

    // Don't check for this unless calibration is active
    Configuration* conf = Configuration::getInstance();
    if( !conf->getBool( "CalibrationActive" ) ){
        return false;
    }

    //Setup for checking by initializing fixed supplies and demands counter to null value
    for( RegionIterator i = regions.begin(); i != regions.end(); ++i ){
        ( *i )->initializeCalValues( period );
    }

    //Setup for checking by adding up all fixed supplies and demands
    for( RegionIterator i = regions.begin(); i != regions.end(); ++i ){
        ( *i )->setCalSuppliesAndDemands( period );
    }

    return false;
}

/*! \brief Calculate supply and demand and emissions for all regions.
* \details Loops through regions so that regions can define their own
*          calculation ordering.
* \param aPeriod Period to calculate.
* \param aRegionsToCalc Optional subset of the regions to calculate, if excluded
*        all regions will be calculated.
*/
void World::calc( const int aPeriod, const AtomVector& aRegionsToCalc ) {   
    // Get the list of valid region numbers to solve.
    const vector<unsigned int> regionNumbersToCalc = getRegionIndexesToCalculate( aRegionsToCalc );

    /*! \invariant The number of regions to calculate must be between 0 and the
    *              number of regions inclusive. 
    */
    assert( regionNumbersToCalc.size() <= regions.size() );

    // Increment the world.calc count based on the number of regions to solve. 
    mCalcCounter->incrementCount( static_cast<double>( regionNumbersToCalc.size() ) / static_cast<double>( regions.size() ) );
    
    // Perform calculation loop on each region to calculate. 
    for ( vector<unsigned int>::const_iterator currIndex = regionNumbersToCalc.begin(); currIndex != regionNumbersToCalc.end(); ++currIndex ) {
        regions[ *currIndex ]->calc( aPeriod );
    }
}

//! Update all summary information for reporting
// Orginally in world.calc, removed to call only once after solved
void World::updateSummary( const list<string> aPrimaryFuelList, const int period ) {
    for( RegionIterator i = regions.begin(); i != regions.end(); i++ ){
        ( *i )->updateSummary( aPrimaryFuelList, period );
        ( *i )->updateAllOutputContainers( period );
    }
}

/*! Calculates the global emissions.
*/
void World::runClimateModel() {
    // Declare visitors which will aggregate emissions by period.
    EmissionsSummer co2Summer( "CO2" );
    LUCEmissionsSummer co2LandUseSummer( "CO2NetLandUse" );
    EmissionsSummer ch4Summer( "CH4" );
    EmissionsSummer ch4agrSummer( "CH4_AGR" );
    EmissionsSummer ch4awbSummer( "CH4_AWB" );
    EmissionsSummer coSummer( "CO" );
    EmissionsSummer coagrSummer( "CO_AGR" );
    EmissionsSummer coawbSummer( "CO_AWB" );
    EmissionsSummer n2oSummer( "N2O" );
    EmissionsSummer n2oagrSummer( "N2O_AGR" );
    EmissionsSummer n2oawbSummer( "N2O_AWB" );
    EmissionsSummer noxSummer( "NOx" );
    EmissionsSummer noxagrSummer( "NOx_AGR" );
    EmissionsSummer noxawbSummer( "NOx_AWB" );
    EmissionsSummer so21Summer( "SO2_1" );
    EmissionsSummer so22Summer( "SO2_2" );
    EmissionsSummer so23Summer( "SO2_3" );
    EmissionsSummer so24Summer( "SO2_4" );
    EmissionsSummer so21awbSummer( "SO2_1_AWB" );
    EmissionsSummer so22awbSummer( "SO2_2_AWB" );
    EmissionsSummer so23awbSummer( "SO2_3_AWB" );
    EmissionsSummer so24awbSummer( "SO2_4_AWB" );
    EmissionsSummer cf4Summer( "CF4" );
    EmissionsSummer c2f6Summer( "C2F6" );
    EmissionsSummer sf6Summer( "SF6" );
    EmissionsSummer hfc125Summer( "HFC125" );
    EmissionsSummer hfc134aSummer( "HFC134a" );
    EmissionsSummer hfc245faSummer( "HFC245fa" );
    EmissionsSummer vocSummer( "NMVOC" );
    EmissionsSummer vocagrSummer( "NMVOC_AGR" );
    EmissionsSummer vocawbSummer( "NMVOC_AWB" );
    EmissionsSummer bcSummer( "BC" );
    EmissionsSummer ocSummer( "OC" );
    EmissionsSummer bcawbSummer( "BC_AWB" );
    EmissionsSummer ocawbSummer( "OC_AWB" );
    
    // Group the EmissionsSummer together for improved performance.
    GroupedEmissionsSummer allSummer;
    allSummer.addEmissionsSummer( &co2Summer );
    allSummer.addEmissionsSummer( &ch4Summer );
    allSummer.addEmissionsSummer( &ch4agrSummer );
    allSummer.addEmissionsSummer( &ch4awbSummer );
    allSummer.addEmissionsSummer( &coSummer );
    allSummer.addEmissionsSummer( &coagrSummer );
    allSummer.addEmissionsSummer( &coawbSummer );
    allSummer.addEmissionsSummer( &n2oSummer );
    allSummer.addEmissionsSummer( &n2oagrSummer );
    allSummer.addEmissionsSummer( &n2oawbSummer );
    allSummer.addEmissionsSummer( &noxSummer );
    allSummer.addEmissionsSummer( &noxagrSummer );
    allSummer.addEmissionsSummer( &noxawbSummer );
    allSummer.addEmissionsSummer( &so21Summer );
    allSummer.addEmissionsSummer( &so22Summer );
    allSummer.addEmissionsSummer( &so23Summer );
    allSummer.addEmissionsSummer( &so24Summer );
    allSummer.addEmissionsSummer( &so21awbSummer );
    allSummer.addEmissionsSummer( &so22awbSummer );
    allSummer.addEmissionsSummer( &so23awbSummer );
    allSummer.addEmissionsSummer( &so24awbSummer );
    allSummer.addEmissionsSummer( &cf4Summer );
    allSummer.addEmissionsSummer( &c2f6Summer );
    allSummer.addEmissionsSummer( &sf6Summer );
    allSummer.addEmissionsSummer( &hfc125Summer );
    allSummer.addEmissionsSummer( &hfc134aSummer );
    allSummer.addEmissionsSummer( &hfc245faSummer );
    allSummer.addEmissionsSummer( &vocSummer );
    allSummer.addEmissionsSummer( &vocagrSummer );
    allSummer.addEmissionsSummer( &vocawbSummer );
    allSummer.addEmissionsSummer( &bcSummer );
    allSummer.addEmissionsSummer( &ocSummer );
    allSummer.addEmissionsSummer( &bcawbSummer );
    allSummer.addEmissionsSummer( &ocawbSummer );

   const double TG_TO_PG = 1000;
   const double N_TO_N2O = 1.571132; 
   const double N_TO_NO2 = 3.2857;
   const double S_TO_SO2 = 2.0; 
   const double HFC_CA_TO_FA = ( 950 / 640 ); 
    
    // Update all emissions values.
    accept( &allSummer, -1 );
	
    // The Climate model reads in data for the base period, so skip passing it in.
    for( int period = 1; period < scenario->getModeltime()->getmaxper(); ++period){
        
        accept( &co2LandUseSummer, period );
        
        // Only set emissions if they are valid. If these are not set
        // MAGICC will use the default values.
        if( co2Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "CO2", period,
                                          co2Summer.getEmissions( period )
                                          / TG_TO_PG );
        }
        
        const int currYear = scenario->getModeltime()->getper_to_yr( period );
        const int startYear = currYear - scenario->getModeltime()->gettimestep( period ) + 1;
        for ( int i = startYear; i <= currYear; i++ ) {
            if( co2LandUseSummer.areEmissionsSet( i ) ){
                mClimateModel->setLUCEmissions( "CO2NetLandUse", i,
                                                co2LandUseSummer.getEmissions( i )
                                                / TG_TO_PG );
            }
        }
        
        if( ch4Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "CH4", period,
                                          ch4Summer.getEmissions( period ) +
                                          ch4agrSummer.getEmissions( period ) + 
                                          ch4awbSummer.getEmissions( period ));
        }

        if( coSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "CO", period,
                                          coSummer.getEmissions( period ) +
                                          coagrSummer.getEmissions( period ) +
                                          coawbSummer.getEmissions( period ));
        }

        // MAGICC wants N2O emissions in Tg N, but miniCAM calculates Tg N2O
        if( n2oSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "N2O", period,
                                          ( n2oSummer.getEmissions( period ) +
                                          n2oawbSummer.getEmissions( period ) +
                                          n2oagrSummer.getEmissions( period )  )
                                        / N_TO_N2O );
        }

        // MAGICC wants NOx emissions in Tg N, but miniCAM calculates Tg NOx
        // FORTRAN code uses the conversion for NO2
        if( noxSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "NOx", period,
                                          ( noxSummer.getEmissions( period ) +
                                            noxagrSummer.getEmissions( period ) +
                                            noxawbSummer.getEmissions( period ))
                                          / N_TO_NO2 );
        }

        // MAGICC wants SO2 emissions in Tg S, but miniCAM calculates Tg SO2
        // Region 1 includes SO21 and 60% of SO24 (FSU)
        if( so21Summer.areEmissionsSet( period ) && so24Summer.areEmissionsSet( period )){
            mClimateModel->setEmissions( "SOXreg1", period,
                                          ( so21Summer.getEmissions( period ) +
                                            so21awbSummer.getEmissions( period )
                                            + 0.6*so24Summer.getEmissions( period ) 
                                            + 0.6*so24awbSummer.getEmissions( period ))
                                          / S_TO_SO2 );
        }

        // MAGICC wants SO2 emissions in Tg S, but miniCAM calculates Tg SO2
        // Region 2 includes SO22 and 40% of SO24 (FSU)
        if( so22Summer.areEmissionsSet( period ) && so24Summer.areEmissionsSet( period )){
            mClimateModel->setEmissions( "SOXreg2", period,
                                          ( so22Summer.getEmissions( period ) +
                                            so22awbSummer.getEmissions( period )
                                            + 0.4*so24Summer.getEmissions( period ) 
                                            + 0.4*so24awbSummer.getEmissions( period ))
                                          / S_TO_SO2 );
        }

        // MAGICC wants SO2 emissions in Tg S, but miniCAM calculates Tg SO2
        if( so23Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "SOXreg3", period,
                                          ( so23Summer.getEmissions( period ) +
                                          so23awbSummer.getEmissions( period ) )
                                          / S_TO_SO2 );
        }

        if( cf4Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "CF4", period,
                                          cf4Summer.getEmissions( period ) );
        }

        if( c2f6Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "C2F6", period,
                                          c2f6Summer.getEmissions( period ) );
        }

        if( sf6Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "SF6", period,
                                          sf6Summer.getEmissions( period ) );
        }

        if( hfc125Summer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "HFC125", period,
                                          hfc125Summer.getEmissions( period ) );
        }

        if( hfc134aSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "HFC134a", period,
                                          hfc134aSummer.getEmissions( period ) );
        }

        // MAGICC needs HFC245fa in kton of HFC245ca
        if( hfc245faSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "HFC245ca", period,
                                          hfc245faSummer.getEmissions( period ) 
                                          / HFC_CA_TO_FA );
        }

        // MAGICC needs this in tons of VOC. Input is in TgC
        if( vocSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "NMVOCs", period,
                                          ( vocSummer.getEmissions( period ) +
                                          vocagrSummer.getEmissions( period ) +
                                          vocawbSummer.getEmissions( period ) ));
        }
        
        // MAGICC needs this in GgC. Model output is in TgC
        if( bcSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "BC", period,
                                         ( bcSummer.getEmissions( period ) +
                                         bcawbSummer.getEmissions( period ) )
                                         * TG_TO_PG );
        }
        
        // MAGICC needs this in GgC. Model output is in TgC
        if( ocSummer.areEmissionsSet( period ) ){
            mClimateModel->setEmissions( "OC", period,
                                         ( ocSummer.getEmissions( period ) +
                                         ocawbSummer.getEmissions( period ) )
                                         * TG_TO_PG );
        }

    }
	
    // Run the model.
    mClimateModel->runModel();
}

//! write results for all regions to file
void World::csvOutputFile() const {

    // Write global data
    csvGlobalDataFile();
    
    for( CRegionIterator i = regions.begin(); i != regions.end(); i++ ){
        ( *i )->csvOutputFile();
    }
}

//! write global results to file
void World::csvGlobalDataFile() const {
    const int maxper = scenario->getModeltime()->getmaxper();
    vector<double> temp(maxper);
    // function protocol
    void fileoutput3(string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    // write total emissions for World
    for ( int m = 0; m < maxper; m++ ){
        // Sum emissions by period.
        for( CRegionIterator iter = regions.begin(); iter != regions.end(); ++iter ) {
            // This interface needs to be fixed.
            temp[ m ] += ( *iter )->getSummary( m ).get_emissmap_second( "CO2" );
        }
    }
    fileoutput3( "global"," "," "," ","CO2 emiss","MTC",temp);

    // Write out concentrations.
    mClimateModel->printFileOutput();
}

//! MiniCAM style output to database
void World::dbOutput( const list<string>& aPrimaryFuelList ) const {
    // Write out concentrations
    mClimateModel->printDBOutput();

    // call regional output
    for( CRegionIterator i = regions.begin(); i != regions.end(); i++ ){
        ( *i )->dbOutput( aPrimaryFuelList );
    }
}

/*! \brief Test to see if calibration worked for all regions
*
* Compares the sum of calibrated + fixed values to output of each sector.
*
* \author Steve Smith
* \param period Model period
* \param calAccuracy Calibration tolerance.
* \param printWarnings flag to turn on logging of warnings if calibrations are
*        not accurate
* \return Boolean true if calibration is ok.
*/
bool World::isAllCalibrated( const int period, double calAccuracy, const bool printWarnings ) const {
    bool isAllCalibrated = true;
    ILogger& calLog = ILogger::getLogger( "calibration_log" );
    calLog.setLevel( ILogger::DEBUG );
    for( CRegionIterator i = regions.begin(); i != regions.end(); i++ ){
        bool currRegionCalibrated = ( *i )->isAllCalibrated( period, calAccuracy, printWarnings );
        isAllCalibrated &= currRegionCalibrated;
        // if we did not calibrate this region correctly and we are printing warnings then give the
        // user some I/O tables to help them understand what was inconsistent
        if( !currRegionCalibrated && printWarnings ) {
            // we want to give the user two tables to use one with a condensed view
            // with all inputs and outputs so they can see what didn't calibrate
            // and another table fully expanded with just the calibrated values
            calLog << "Energy balance table where inputs and outputs have been replaced by a"
                << " calibrated value if it exists:" << endl;
            EnergyBalanceTable condensedTable( (*i)->getName(), calLog, true, true );
            (*i)->accept( &condensedTable, period );
            condensedTable.finish();
            
            calLog << "Full energy balalce table with just cal values:" << endl;
            EnergyBalanceTable fullTable( (*i)->getName(), calLog, false, false );
            (*i)->accept( &fullTable, period );
            fullTable.finish();
        }
    }
	
    return isAllCalibrated;
}

/*! \brief This function returns a special mapping of strings to ints for use in
*          the outputs. 
* \details This map is created such that global maps to zero, region 0 maps to
*          1, etc. It is similiar to the regionNamesToNumbers map but has the
*          global element and each region number in the regionMap is 1 + the
*          number in the regionNamesToNumbers map.
* \warning This function should only be used by the database output functions. 
* \return The map of region names to numbers.
*/
const map<string,int> World::getOutputRegionMap() const {
    map<string,int> regionMap;

    for ( unsigned int i = 0; i < regions.size(); i++ ) {
        regionMap[regions[i]->getName()] = i+1; // start index from 1
    }
    // hardcode for now
    regionMap["global"] = 0;
    return regionMap;
}

/*! \brief This function returns a vector of IDs for all regions that exist in
*          the world.
* \details This function creates a vector of region IDs in the same order in
*          which they exist in the world.
* \return A constant vector of region IDs.
*/
const World::AtomVector World::getRegionIDs() const {
    AtomVector regionIDs;
    
    // Iterate over the regions and lookup the ID for each region.
    for( CRegionIterator i = regions.begin(); i != regions.end(); ++i ) {
        // Find the atom for the region name.
        const objects::Atom* regionID = objects::AtomRegistry::getInstance()->findAtom( (*i)->getName() );
        /*! \invariant The region atom has already been created, otherwise the
        *              atom was not registered when the fast lookup map was
        *              created. 
        */
        assert( regionID );
        
        /*! \invariant The ID of region atom found is equal to the region name. */
        assert( regionID->getID() == (*i)->getName() );
        regionIDs.push_back( regionID );
    }
    /*! \post The size of the region ID vector is equal to the size of the
    *         region vector. 
    */
    assert( regionIDs.size() == regions.size() );
    return regionIDs;
}

/*! \brief Set a fixed tax for all regions.
* \param aTax Tax.
*/
void World::setTax( const GHGPolicy* aTax ){
    for( RegionIterator iter = regions.begin(); iter != regions.end(); ++iter ){
        (*iter)->setTax( aTax );
    }
}

/*! \brief Get the climate model.
* \return The climate model.
*/
const IClimateModel* World::getClimateModel() const {
    return mClimateModel.get();
}

/*! \brief A function to generate a series of ghg emissions quantity curves based on an already performed model run.
* \details This function used the information stored in it to create a series of curves, one for each region,
* with each datapoint containing a time period and an amount of gas emissions.
* \note The user is responsible for deallocating the memory in the returned Curves.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a set of curves for.
* \return A map with keys as region names and Curves as values representing the quantity of ghg emissions by time period.
*/
const map<const string,const Curve*> World::getEmissionsQuantityCurves( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    const string GLOBAL_NAME = "global";

    map<const string,const Curve*> emissionsQCurves;

    for( CRegionIterator rIter = regions.begin(); rIter != regions.end(); rIter++ ){
        emissionsQCurves[ (*rIter)->getName() ] = (*rIter)->getEmissionsQuantityCurve( ghgName );
    }

    // Add an entry for the global emissions. Should do this better. 
    ExplicitPointSet* globalQs = new ExplicitPointSet();
    const Marketplace* marketplace = scenario->getMarketplace();
    const Modeltime* modeltime = scenario->getModeltime();

    for( int per = 0; per < modeltime->getmaxper(); per++ ){
        globalQs->addPoint( new XYDataPoint( modeltime->getper_to_yr( per ), marketplace->getDemand( ghgName, "USA", per ) ) );
    }
    emissionsQCurves[ GLOBAL_NAME ] = new PointSetCurve( globalQs );
    return emissionsQCurves;
}

/*! \brief A function to generate a series of ghg emissions price curves based on an already performed model run.
* \details This function used the information stored in it to create a series of curves, one for each period,
* with each datapoint containing a time period and the price gas emissions. 
* \note The user is responsible for deallocating the memory in the returned Curves.
* \author Josh Lurz
* \param ghgName The name of the ghg to create a set of Curves for.
* \return A map with keys as region names and Curves as values representing the price of ghg emissions by time period. 
*/
const map<const string,const Curve*> World::getEmissionsPriceCurves( const string& ghgName ) const {
    /*! \pre The run has been completed. */
    map<const string,const Curve*> emissionsPCurves;
    const string GLOBAL_NAME = "global";
    
    for( CRegionIterator rIter = regions.begin(); rIter != regions.end(); rIter++ ){
        emissionsPCurves[ (*rIter)->getName() ] = (*rIter)->getEmissionsPriceCurve( ghgName );
    }

    // Add an entry for the global emissions. Should do this better. 
    ExplicitPointSet* globalQs = new ExplicitPointSet();
    const Marketplace* marketplace = scenario->getMarketplace();
    const Modeltime* modeltime = scenario->getModeltime();
    for( int per = 0; per < modeltime->getmaxper(); per++ ){
        globalQs->addPoint( new XYDataPoint( modeltime->getper_to_yr( per ), marketplace->getPrice( ghgName, "USA", per ) ) );
    }
    emissionsPCurves[ GLOBAL_NAME ] = new PointSetCurve( globalQs );
    return emissionsPCurves;
}

/*!
 * \brief Gets the reference to the calc counter.
 * \details The calc counter would be needed by solvers however since there
 *          could be many solvers the world object will contain it and each
 *          solver will be able to get it from here.
 * \return A reference to the Calc Counter
 */
CalcCounter* World::getCalcCounter() const {
    return mCalcCounter.get();
}

/*! \brief Protected function which takes a listing of region names to calculate
*          and returns a list of region indexes to calculate. 
* \details This function translates a passed in list of region names to solve to
*          a vector of valid region numbers to solve. If passed an empty list,
*          it will return the full list of region numbers to solve, this is the
*          default. It checks whether the region name is valid, but does not
*          check for duplicates.
* \param aRegionsToCalc A vector of region names to calculate.
* \return A vector of region numbers to solve.
*/
const vector<unsigned int> World::getRegionIndexesToCalculate( const AtomVector& aRegionsToCalc ){
    vector<unsigned int> regionNumbersToCalc;
    // Check for the empty list of names, return the full list of region numbers.
    if ( aRegionsToCalc.empty() ) {
        for( unsigned int regionNumber = 0; regionNumber < regions.size(); ++regionNumber ) {
            regionNumbersToCalc.push_back( regionNumber );
        }
    }
    // Check if each name is valid and add its region number to the vector if it is. 
    else {
        for( AtomVector::const_iterator regionID = aRegionsToCalc.begin();
             regionID != aRegionsToCalc.end(); ++regionID )
        {
            // Lookup the region atom in the world's fast lookup map.
            FastRegionLookupMap::const_iterator iter = mRegionLookupMap->find( *regionID );

            /*! \invariant The region ID was found otherwise the solution
            *              mechanism is instructing the world to calculate a
            *              non-existant region. 
            */
            assert( iter != mRegionLookupMap->end() );

            /*! \invariant The ID of the atom is equal to the region name of the
            *              region at the index. 
            */
            assert( (*regionID)->getID() == regions[ iter->second ]->getName() );
            
            // Add the index to the list of indexes to calculate.
            regionNumbersToCalc.push_back( iter->second );
        }
    }

    /*! \post The number of regions to calculate is between 0 and the total
    *         number of regions inclusive.
    */
    assert( regionNumbersToCalc.size() <= regions.size() );

    return regionNumbersToCalc;
}

/*! \brief Call any calculations that are only done once per period after
*          solution is found.
* \details This function is used to calculate and store variables which are only
*          needed after the current period is complete. 
* \param aPeriod The period to finalize.
* \author Sonny Kim, Josh Lurz
*/
void World::postCalc( const int aPeriod ){
    // Finalize sectors.
    for( RegionIterator region = regions.begin(); region != regions.end(); ++region ){
        (*region)->postCalc( aPeriod );
    }
    
    const int numGases = 8;
    const int numCategories = 12;
    const string gases[] = { "BC", "CH4", "CO", "NH3", "NMVOC", "NOx", "OC", "SO2" };
    const string categories[] = { "AGR", "AIR", "AWB", "DOM", "ENE", "IND", "LCF", "SAV", "SHIP", "SLV", "TRA", "WST" };
    
    ILogger& rcpLog = ILogger::getLogger( "rcp_log" );
    RCPEmissionsVisitor emissionsVisitor;
    const Modeltime* modeltime = scenario->getModeltime();
    const int year = modeltime->getper_to_yr( aPeriod );
    accept( &emissionsVisitor, aPeriod );
    for( RegionIterator region = regions.begin(); region != regions.end(); ++region ){
        for( int gasI = 0; gasI < numGases; ++gasI ) {
            for( int catI = 0; catI < numCategories; ++catI ) {
                rcpLog << ( region - regions.begin() + 1 ) << ", " << gases[ gasI ] << ", " << categories[ catI ] << ", "
                       << year << ", "
                       << emissionsVisitor.getEmissions( (*region)->getName(), categories[ catI ], gases[ gasI ], year ) << endl;
            }
        }
    }
}

void World::csvSGMOutputFile( ostream& aFile, const int period ) const {
    for( CRegionIterator rIter = regions.begin(); rIter != regions.end(); ++rIter ){
        ( *rIter )->csvSGMOutputFile( aFile, period );
    }
}

void World::csvSGMGenFile( ostream& aFile ) const {
    for( CRegionIterator rIter = regions.begin(); rIter != regions.end(); ++rIter ){
        ( *rIter )->csvSGMGenFile( aFile );
    }
}

/*! \brief Update a visitor for the World.
* \param aVisitor Visitor to update.
* \param aPeriod Period to update.
*/
void World::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitWorld( this, aPeriod );

    // Visit the marketplace
    scenario->getMarketplace()->accept( aVisitor, aPeriod );

    // Visit the climate model.
    mClimateModel->accept( aVisitor, aPeriod );

    // loop for regions
    for( CRegionIterator currRegion = regions.begin(); currRegion != regions.end(); ++currRegion ){
        (*currRegion)->accept( aVisitor, aPeriod );
    }

    aVisitor->endVisitWorld( this, aPeriod );
}

/*!
 * \brief Get a reference to the visitor which is used to adjust carbon densities.
 * \details Adjusting carbon densities can be used to incorparate feedbacks and will
 *          affect land use change emissions as well as yields.
 * \return A reference to the member visitor object which can be used to set carbon
 *         densities.
 */
SetCarbonDensity& World::getAdjustCarbonDensityVisitor() {
    return mAdjDensityVisitor;
}
std::vector<Region*> World::getRegions() {
    return regions;
}

