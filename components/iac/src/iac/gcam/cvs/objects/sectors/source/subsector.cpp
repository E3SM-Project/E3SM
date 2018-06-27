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
* \file subsector.cpp
* \ingroup Objects
* \brief Subsector class source file.
* \author Sonny Kim, Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "util/base/include/configuration.h"
#include "sectors/include/subsector.h"
#include "technologies/include/itechnology_container.h"
#include "technologies/include/technology_container.h"
#include "technologies/include/stub_technology_container.h"
#include "technologies/include/itechnology.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "util/base/include/xml_helper.h"
#include "marketplace/include/marketplace.h"
#include "util/base/include/summary.h"
#include "containers/include/gdp.h"
#include "containers/include/info_factory.h"
#include "containers/include/iinfo.h"
#include "technologies/include/base_technology.h"
#include "consumers/include/consumer.h"
#include "consumers/include/household_consumer.h"
#include "consumers/include/govt_consumer.h"
#include "consumers/include/trade_consumer.h"
#include "consumers/include/invest_consumer.h"
#include "technologies/include/production_technology.h"
#include "util/base/include/ivisitor.h"
#include "technologies/include/technology_type.h"
#include "investment/include/idistributor.h"
#include "investment/include/iexpected_profit_calculator.h"
#include "sectors/include/sector_utils.h"
#include "investment/include/investment_utils.h"
#include "reporting/include/indirect_emissions_calculator.h"
#include "util/base/include/interpolation_rule.h"

using namespace std;
using namespace xercesc;
using namespace objects;

extern Scenario* scenario;
// static initialize.
const string Subsector::XML_NAME = "subsector";

/*! \brief Default constructor.
*
* Constructor initializes member variables with default values, sets vector sizes, etc.
*
* \author Sonny Kim, Steve Smith, Josh Lurz
*/
const Value LOGIT_EXP_DEFAULT = -6;

Subsector::Subsector( const string& aRegionName, const string& aSectorName ):
regionName( aRegionName ),
sectorName( aSectorName ),
doCalibration( false ),
mTechLogitExp( LOGIT_EXP_DEFAULT )
{
    // resize vectors.
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    summary.resize(maxper); // object containing summaries
    fuelPrefElasticity.resize( maxper );
    summary.resize( maxper );
    mInvestments.resize( maxper );
    mFixedInvestments.resize( maxper, -1 );
}

/*! \brief Default destructor.
*
* deletes all Technology objects associated  with this sector.
*
* \author Josh Lurz
*/
Subsector::~Subsector() {
    clear();
}

//! Deallocate the subsector memory.
void Subsector::clear(){
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        delete *techIter;
    }
    
    for( BaseTechIterator delTech = baseTechs.begin(); delTech != baseTechs.end(); ++delTech ){
        delete *delTech;
    }
    for( map<string, TechnologyType*>::iterator techType = mTechTypes.begin(); techType != mTechTypes.end();
        ++techType )
    {
        delete techType->second;
    }
    clearInterpolationRules();
}

/*!
 * \brief A helper method to delete and clear the vector of interpolation rules.
 */
void Subsector::clearInterpolationRules() {
    for( CInterpRuleIterator ruleIter = mShareWeightInterpRules.begin(); ruleIter != mShareWeightInterpRules.end(); ++ruleIter ) {
        delete *ruleIter;
    }
    mShareWeightInterpRules.clear();
}

/*! \brief Returns sector name
*
* \author Sonny Kim
* \return sector name as a string
*/
const string& Subsector::getName() const {
    return name;
}

//! Initialize Subsector with xml data
void Subsector::XMLParse( const DOMNode* node ) {   

    /*! \pre Make sure we were passed a valid node. */
    assert( node );

    // get the name attribute.
    name = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    const Modeltime* modeltime = scenario->getModeltime();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "share-weight" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mParsedShareWeights, modeltime );
        }
        else if( nodeName == "logit-exponent" ){
            XMLHelper<Value>::insertValueIntoVector( curr, mTechLogitExp, modeltime );
        }
        else if( nodeName == "fuelprefElasticity" ){
            XMLHelper<double>::insertValueIntoVector( curr, fuelPrefElasticity, modeltime );  
        }
        // Fixed investment
        else if( nodeName == "FixedInvestment" ){
            XMLHelper<double>::insertValueIntoVector( curr, mFixedInvestments, scenario->getModeltime() );
        }
        // household consumer object for final demands
        else if( nodeName == HouseholdConsumer::getXMLNameStatic() ) {
            parseBaseTechHelper( curr, new HouseholdConsumer() );
        }
        // government consumer object for final demands
        else if( nodeName == GovtConsumer::getXMLNameStatic() ) {
            parseBaseTechHelper( curr, new GovtConsumer() );
        }
        // Trade consumer object for final demands
        else if( nodeName == TradeConsumer::getXMLNameStatic() ) {
            parseBaseTechHelper( curr, new TradeConsumer() );
        }
        // government consumer object for final demands
        else if( nodeName == InvestConsumer::getXMLNameStatic() ) {
            parseBaseTechHelper( curr, new InvestConsumer() );
        }
        // production technology object for production sectors
        else if( nodeName == ProductionTechnology::getXMLNameStatic() ) {
            parseBaseTechHelper( curr, new ProductionTechnology() );
        }
        else if( nodeName == InterpolationRule::getXMLNameStatic() && XMLHelper<string>::getAttr( curr, "apply-to" ) == "share-weight" ) {
            // if the delete flag is set then for interpolation rules that means to clear
            // out any previously parsed rules
            if( XMLHelper<bool>::getAttr( curr, "delete" ) ) {
                clearInterpolationRules();
            }

            InterpolationRule* tempRule = new InterpolationRule();
            tempRule->XMLParse( curr );
            mShareWeightInterpRules.push_back( tempRule );
        }
        else if( TechnologyContainer::hasTechnologyType( nodeName ) ) {
            parseContainerNode( curr, mTechContainers, new TechnologyContainer );
        }
        else if( nodeName == StubTechnologyContainer::getXMLNameStatic() ) {
            parseContainerNode( curr, mTechContainers, new StubTechnologyContainer );
        }
        // parsed derived classes
        else if( !XMLDerivedClassParse( nodeName, curr ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Unknown element " << nodeName << " encountered while parsing " << getXMLName() << endl;
        }
    }
}

//! Helper function which parses any type of base Technology correctly.
void Subsector::parseBaseTechHelper( const DOMNode* aCurr, BaseTechnology* aNewTech ){
    // Ensure a valid technology was passed.
    assert( aNewTech );

    // Use an auto_ptr to take responsibility for the memory.
    auto_ptr<BaseTechnology> newTech( aNewTech );
    
    // Check if the base technology already exists.
    const string id = BaseTechnology::createIdentifier( XMLHelper<string>::getAttr( aCurr, "name" ),
                      XMLHelper<int>::getAttr( aCurr, "year" ) );

    map<string,int>::const_iterator baseTechMapIter = baseTechNameMap.find( id );
    if( baseTechMapIter != baseTechNameMap.end() ) { 
        // already exists, so tell the existing one to parse
        baseTechs[ baseTechMapIter->second ]->XMLParse( aCurr );
    }
    else { 
        // doesn't exist so use the new passed in base Technology type.
        newTech->XMLParse( aCurr );

        // Add the new Technology to the vector and the map.
        baseTechs.push_back( newTech.release() ); // Releases ownership of the memory.
        baseTechNameMap[ baseTechs.back()->getIdentifier() ] = static_cast<int>( baseTechs.size() ) - 1;

        // the Technology type may not exist yet.
        map<string,TechnologyType*>::iterator typePos = mTechTypes.find( baseTechs.back()->getName() );
        if( typePos == mTechTypes.end() ){
            // create the tech type, set the iterator to the new item.
            // Insert returns the pair of the iterator position the item was inserted in and whether 
            // the item was inserted, so set the iterator to the first spot in the pair.
            typePos = mTechTypes.insert( make_pair( baseTechs.back()->getName(), new TechnologyType ) ).first;
        }
        typePos->second->addVintage( baseTechs.back() );

        // Set the Technology type helper object to the Technology. This may be moved to the constructor
        // or removed if Technology type is made to inherit from IInvestable.
        baseTechs.back()->setTypeHelper( typePos->second );
    }
}

//! Parses any input variables specific to derived classes
bool Subsector::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ) {
    // do nothing
    // defining method here even though it does nothing so that we do not
    // create an abstract class.
    return false;
}

//! Output the Subsector member variables in XML format.
void Subsector::toInputXML( ostream& out, Tabs* tabs ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteOpeningTag( getXMLName(), out, tabs, name );

    // TODO: create a XMLWriteVector that skips !Value.isInited() rather than a default value.
    for( unsigned int period = 0; period < mParsedShareWeights.size(); period++ ){
        // Determine the correct year.
        unsigned int year = modeltime->getper_to_yr( period );

        if( mParsedShareWeights[ period ].isInited() ) {
            XMLWriteElement( mParsedShareWeights[ period ], "share-weight", out, tabs, year, "", false );
        }
    }

    for( CInterpRuleIterator ruleIt = mShareWeightInterpRules.begin(); ruleIt != mShareWeightInterpRules.end(); ++ruleIt ) {
        (*ruleIt)->toInputXML( out, tabs );
    }

    XMLWriteVector( mTechLogitExp, "logit-exponent", out, tabs, modeltime, LOGIT_EXP_DEFAULT );
    
    XMLWriteVector( fuelPrefElasticity, "fuelprefElasticity", out, tabs, modeltime, 0.0 );

    toInputXMLDerived( out, tabs );

    for ( unsigned int i = 0; i < baseTechs.size(); i++ ){
        baseTechs[i]->toInputXML( out, tabs );
    }
    
    XMLWriteVector( mFixedInvestments, "FixedInvestment", out, tabs, modeltime, -1.0 );

    // write out the technology objects.
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->toInputXML( out, tabs );
    }    
    
    // finished writing xml for the class members.
    
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Write information useful for debugging to XML output stream
*
* Function writes market and other useful info to XML. Useful for debugging.
*
* \author Josh Lurz
* \param period model period
* \param out reference to the output stream
*/
void Subsector::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    
    XMLWriteOpeningTag( getXMLName(), out, tabs, name );
    
    // Write the data for the current period within the vector.
    XMLWriteElement( mShareWeights[ period ], "share-weight", out, tabs );
    XMLWriteElement( mTechLogitExp[ period ], "logit-exponent", out, tabs );
    XMLWriteElement( fuelPrefElasticity[ period ], "fuelprefElasticity", out, tabs );
    XMLWriteElement( getEnergyInput( period ), "input", out, tabs );
    XMLWriteElement( getOutput( period ), "output", out, tabs );
    XMLWriteElement( getTotalCalOutputs( period ), "total-cal-outputs", out, tabs );
    XMLWriteElement( mInvestments[ period ], "investment", out, tabs );
    XMLWriteElement( mFixedInvestments[ period ], "FixedInvestment", out, tabs );
    XMLWriteElement( getCalibrationStatus( period ), "calibration-status", out, tabs );

    toDebugXMLDerived( period, out, tabs );
    // Write out the summary object.
    // summary[ period ].toDebugXML( period, out );
    // write out the Technology objects.

    for ( unsigned int j = 0; j < baseTechs.size(); j++ ) {
        // This isn't right, techs with years other the current year could change output.
        if (baseTechs[j]->getYear() == scenario->getModeltime()->getper_to_yr( period ) ) {
            baseTechs[j]->toDebugXML( period, out, tabs );
        }
    }
    
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->toDebugXML( period, out, tabs );
    }    
    
    // finished writing xml for the class members.
    
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& Subsector::getXMLName() const {
    return XML_NAME;
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
const string& Subsector::getXMLNameStatic() {
    return XML_NAME;
}

/*! \brief Complete the initialization
*
* This routine is only called once per model run
* \param aSectorInfo The parent sector info object.
* \param aDependencyFinder The regional dependency finder.
* \param aLandAllocator Regional land allocator.
* \author Josh Lurz
* \warning markets are not necessarily set when completeInit is called
*/
void Subsector::completeInit( const IInfo* aSectorInfo,
                              DependencyFinder* aDependencyFinder,
                              ILandAllocator* aLandAllocator )
{
    mSubsectorInfo.reset( InfoFactory::constructInfo( aSectorInfo, regionName + "-" + sectorName + "-" + name ) );
    
    for( unsigned int i = 0; i < baseTechs.size(); i++) {
        baseTechs[i]->completeInit( regionName, sectorName, name );
    }
    // TODO: make sure that isInitialYear() flag on the baseTechs is consistent
    // They would be consisitent if exactly one technology per tech name had the 
    // flag set to true.

    for( unsigned int j = 0; j < baseTechs.size(); ++j ){
        // Remove empty inputs for the initial year of the tech only.
        // removing unecessary inputs into the future will be handled
        // by the nested input structure in the technology
        if( baseTechs[ j ]->isInitialYear() ){
            baseTechs[ j ]->removeEmptyInputs();
        }
    }
    
    for ( TechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->completeInit( regionName, sectorName, name, aDependencyFinder, mSubsectorInfo.get(),
                                   aLandAllocator );
    }    
    
    const Modeltime* modeltime = scenario->getModeltime();
    // Initialize working share weights with parsed share weights up to and including
    // final calibration period.
    for( int per = 0; per <= modeltime->getFinalCalibrationPeriod(); ++per ) {
       if( mParsedShareWeights[ per ].isInited() ) {
           mShareWeights[ per ] = mParsedShareWeights[ per ];
       }
    }

    // Set missing technology logit exponent using the next available parsed
    // technology logit exponent.
    SectorUtils::fillMissingPeriodVectorNextAvailable( mTechLogitExp );
}

/*!
* \brief Perform any initializations needed for each period.
* \details Perform any initializations or calculations that only need to be done
*          once per period (instead of every iteration) should be placed in this
*          function.
* \warning The ghg part of this routine assumes the existence of technologies in
*          the previous and future periods
* \author Steve Smith, Sonny Kim
* \param aNationalAccount National accounts container.
* \param aDemographics Regional demographics container.
* \param aMoreSectorInfo SGM sector info object.
* \param aPeriod Model period
*/
void Subsector::initCalc( NationalAccount* aNationalAccount,
                          const Demographic* aDemographics,
                          const MoreSectorInfo* aMoreSectorInfo,
                          const int aPeriod )
{
    // Initialize all technologies.
    for( TechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->initCalc( regionName, sectorName, mSubsectorInfo.get(), aDemographics, aPeriod );
    }

    // Initialize the baseTechs. This might be better as a loop over tech types. 
    const Modeltime* modeltime = scenario->getModeltime();
    for( unsigned int j = 0; j < baseTechs.size(); j++ ){
        if( aPeriod == 0 && baseTechs[ j ]->isInitialYear() ){
            // TODO: remove capital stock as it is no longer used
            double totalCapital = 0;
            baseTechs[ j ]->initCalc( aMoreSectorInfo, regionName, sectorName, *aNationalAccount,
                                      aDemographics, totalCapital, aPeriod );
        
            // copy base year tech to old vintages
            mTechTypes[ baseTechs[ j ]->getName() ]->initializeTechsFromBase( baseTechs[ j ]->getYear(),
                                      aMoreSectorInfo, regionName, sectorName, *aNationalAccount,
                                      aDemographics, totalCapital );
        }
        // If the current tech is from the previous period, initialize the current tech with its parameters.
        else if ( aPeriod > 0 && baseTechs[ j ]->getYear() == modeltime->getper_to_yr( aPeriod - 1 ) ) {
            BaseTechnology* newTech = mTechTypes[ baseTechs[ j ]->getName() ]->initOrCreateTech( modeltime->getper_to_yr( aPeriod ), baseTechs[ j ]->getYear() );
            // Check if initOrCreate created a Technology which needs to be added to the base tech vector and map.
            if( newTech ){
                // If the tech already existed, it will get initCalc called on it later in this loop. 
                baseTechs.push_back( newTech );
                baseTechNameMap[ baseTechs.back()->getName() + util::toString( baseTechs.back()->getYear() ) ] = static_cast<int>( baseTechs.size() ) - 1;
            }
        } 
    }

    if(aPeriod > 0){
        for( unsigned int j = 0; j < baseTechs.size(); j++ ){
            if( baseTechs[ j ]->getYear() <= modeltime->getper_to_yr( aPeriod ) ){
                baseTechs[ j ]->initCalc( aMoreSectorInfo, regionName, sectorName, *aNationalAccount, aDemographics, 0, aPeriod );
            }
        }
    }

    // If calibration is active, reinitialize share weights.
    if( Configuration::getInstance()->getBool( "CalibrationActive" ) ){
        // Check for zero shareweight for subsector with calibration values
        if( getTotalCalOutputs( aPeriod ) > util::getSmallNumber() 
            && mShareWeights[ aPeriod ] == 0 ) 
        {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::NOTICE );
            mainLog << "Resetting zero shareweight for Subsector " << getName()
                << " in sector " << sectorName << " in region " << regionName
                << " since calibration values are present." << endl;
            mShareWeights[ aPeriod ] = 1.0;
        }

        // Reinitialize share weights to 1 for competing subsector with non-zero read-in share weight
        // for calibration periods only, but only if calibration values are read in any of the technologies
        // in this subsector (as there are cases where you want to fix these shares on a pass-through sector).
        else if( mShareWeights[ aPeriod ] != 0 && aPeriod <= modeltime->getFinalCalibrationPeriod() 
                && getTotalCalOutputs( aPeriod ) > 0.0 ){
            // Reinitialize to 1 to remove bias, calculate new share weights and
            // normalize in postCalc to anchor to dominant subsector.
            mShareWeights[ aPeriod ] = 1.0;
        }
    }

    // For subsectors with only fixed output technologies or all null technology share weights.
    if( containsOnlyFixedOutputTechnologies( aPeriod ) ){
        // Reset share weight for all periods to 0 to indicate that it does not have an 
        // impact on the results.
        mParsedShareWeights[ aPeriod ].set( 0 );
        mShareWeights[ aPeriod ].set( 0 );
    }

    // Apply share weight interpolation rules and fill in missing share weights.
    interpolateShareWeights( aPeriod );
}

/*! \brief Returns the subsector price.
* \details Calculates and returns share-weighted total price (subsectorprice)
*          and cost of fuel (fuelprice). 
* \author Sonny Kim
* \param aGDP Regional GDP object.
* \param aPeriod Model period
*/
double Subsector::getPrice( const GDP* aGDP, const int aPeriod ) const {
    double subsectorPrice = 0; // initialize to 0 for summing
    const vector<double> techShares = calcTechShares( aGDP, aPeriod );
    for ( unsigned int i = 0; i < mTechContainers.size(); ++i ) {
        // Technologies with zero share cannot affect the marginal price.
        if( techShares[ i ] > util::getSmallNumber() ){
            double currCost = mTechContainers[i]->getNewVintageTechnology(aPeriod)->getCost( aPeriod );
            // calculate weighted average price for Subsector
            if( currCost >= util::getSmallNumber() ) {
                subsectorPrice += techShares[ i ] * currCost;
            }
            // We want to allow regional and delivered biomass prices to be negative
            else if ( ( sectorName == "regional biomass"  || sectorName == "delivered biomass" ) ){
                subsectorPrice += techShares[ i ] * currCost;
            }

        }
    }

    // We want to allow regional and delivered biomass prices to be negative
    if ( ( sectorName == "regional biomass"  || sectorName == "delivered biomass" ) 
        && ( subsectorPrice < 0 ) ){
        return subsectorPrice;
    }


    // Check for the condition where all technologies were fixed.
    return ( subsectorPrice >= util::getSmallNumber() ) ? subsectorPrice : -1;
}

/*! \brief Returns whether the subsector should be calibrated.
* \details If either the Subsector output, or the output of all the technologies
*          under this Subsector (not including those with zero output) are
*          calibrated, then the Subsector should calibrate.
* \author Steve Smith
* \param aPeriod Model period
*/
bool Subsector::getCalibrationStatus( const int aPeriod ) const {
    
    // Check all the technologies for the period.
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        // Check whether there is any calibration input, not one for a specific fuel.
        if ( mTechContainers[ i ]->getNewVintageTechnology( aPeriod )->hasCalibratedValue( aPeriod) ) {
            return true;
        }
    }
    return false;
}


/*! \brief returns Subsector fuel price times share
* \details Returns the share-weighted fuel price, which is later summed to get
*          the sector-weighted fuel price.
* \author Sonny Kim
* \param aGDP GDP container.
* \param aPeriod Model period.
* \return share-weighted fuel price
*/
double Subsector::getAverageFuelPrice( const GDP* aGDP, const int aPeriod ) const {
    // Determine the average fuel price.
    double fuelPrice = 0;

    // The base period is not solved so the current shares can be calculated and
    // used. In future periods the previous period's shares must be used as the
    // current period's are unknown.
    const int sharePeriod = ( aPeriod == 0 ) ? aPeriod : aPeriod - 1;

    const vector<double> techShares = calcTechShares( aGDP, sharePeriod );
    for ( unsigned int i = 0; i < mTechContainers.size(); ++i) {
        // calculate weighted average price of fuel only
        // Technology shares are based on total cost
        fuelPrice += techShares[ i ] * mTechContainers[i]->getNewVintageTechnology( aPeriod )->getEnergyCost( regionName, sectorName, aPeriod );
    }
    /*! \post Fuel price must be positive. */
    assert( fuelPrice >= 0 );
    return fuelPrice;
}

/*! \brief calculate Technology shares within Subsector
*
* Calls Technology objects to first calculate cost, then their share. Follows this by normalizing shares. 
*
* \author Marshall Wise, Josh Lurz
* \param regionName region name
* \param period model period
* \return A vector of technology shares.
*/
const vector<double> Subsector::calcTechShares( const GDP* aGDP, const int aPeriod ) const {
    vector<double> techShares( mTechContainers.size() );
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        // determine shares based on Technology costs
        techShares[ i ] = mTechContainers[i]->getNewVintageTechnology(aPeriod)->calcShare( regionName, sectorName, aGDP,
            mTechLogitExp[ aPeriod ], aPeriod );

        // Check that Technology shares are valid.
        assert( util::isValidNumber( techShares[ i ] ) );
    }
    // Normalize technology shares.
    SectorUtils::normalizeShares( techShares );
    return techShares;
}

/*!
* \brief Calculate the cost of the Subsector.
* \details Instructs all technologies to calculate their costs. The subsector
*          can calculate it's costs dynamically once all Technologies have
*          calculated their costs, so the Subsector cost is not stored.
* \param aPeriod Model period.
*/
void Subsector::calcCost( const int aPeriod ){
    // Instruct all technologies up to and including the current period to
    // calculate their costs. Future Technologies cannot have a cost as they do
    // not yet exist.
    for( TechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        for( ITechnologyContainer::TechRangeIterator vintageIter = (*techIter)->getVintageBegin( aPeriod ); vintageIter != (*techIter)->getVintageEnd(); ++vintageIter ) {
            (*vintageIter).second->calcCost( regionName, sectorName, aPeriod );
        }
    }
}

/*! \brief calculate Subsector unnormalized shares
* \details Calculates the unnormalized share for this sector. Also calculates
*          the sector aggregate price (or cost)
* \author Sonny Kim, Josh Lurz
* \param aPeriod model period
* \param aGdp gdp object
* \param aLogitExp The logit exponent which controls the behavoir of this
*                  subsector nest.
* \warning There is no difference between demand and supply technologies.
*          Control behavior with value of parameter fuelPrefElasticity
* \return The subsector share.
*/
double Subsector::calcShare(const int aPeriod, const GDP* aGdp, const double aLogitExp ) const {
    double subsectorPrice = getPrice( aGdp, aPeriod );

    if( subsectorPrice >= util::getSmallNumber() ){
        double scaledGdpPerCapita = aGdp->getBestScaledGDPperCap( aPeriod );
        double share = mShareWeights[ aPeriod ] * pow( subsectorPrice, aLogitExp )
                        * pow( scaledGdpPerCapita, fuelPrefElasticity[ aPeriod ] );
        /*! \post Share is zero or positive. */
        // Check for invalid shares.
        if( share < -util::getSmallNumber() || !util::isValidNumber( share ) ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Invalid share for " << name << " in " << regionName 
                << " of " << share << endl;
        }
        return share;
    }

    // We want to allow regional and delivered biomass prices to be negative
    if( subsectorPrice < 0 && ( sectorName == "regional biomass"  || sectorName == "delivered biomass" ) ) {
        // This assumes that there is only one subsector in regional biomass
        return 1;
    }

    return 0;
}

/*! \brief Return the total exogenously fixed Technology output for this sector.
* \author Steve Smith
* \param period model period
*/
double Subsector::getFixedOutput( const int aPeriod ) const {
    double fixedOutput = 0;
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::CTechRangeIterator endIter = (*techIter)->getVintageEnd();
        for( ITechnologyContainer::CTechRangeIterator vintageIter = (*techIter)->getVintageBegin( aPeriod ); vintageIter != endIter ; ++vintageIter ) {
            double currFixedOutput = (*vintageIter).second->getFixedOutput( regionName, sectorName, false, "", aPeriod );
            /*! \invariant Fixed output for each Technology must be -1 or
            *              positive. 
            */
            assert( fixedOutput == -1 || fixedOutput >= 0 );
            if( currFixedOutput > 0 ){
                fixedOutput += currFixedOutput;
            }
        }
    }
    /*! \post Fixed output total must be positive. */
    assert( fixedOutput >= 0 );
    return fixedOutput;
}

/*!
 * \brief Apply share weight interpolation rules for the subsector
 * \details Rules will only apply in the first period after calibration.  It
 *          is an error to have uninitialized share weights after rules have
 *          been applied.
 * \param aPeriod Current model period.
 */
void Subsector::interpolateShareWeights( const int aPeriod ) {
    // Do not apply rules if this is not the first period after calibration.
    const Modeltime* modeltime = scenario->getModeltime();
    if( aPeriod != ( modeltime->getFinalCalibrationPeriod() + 1 ) ) {
        return;
    }

    // Make sure that calibrated values get stored back into the parsed share weights vector so that
    // they get written out.  All other parsed values will initialize the working share weights
    for( int per = 0; per < mParsedShareWeights.size(); ++per ) {
        if( per <= modeltime->getFinalCalibrationPeriod() ) {
            mParsedShareWeights[ per ].set( mShareWeights[ per ] );
        }
        else if( mParsedShareWeights[ per ].isInited() )  {
            mShareWeights[ per ] = mParsedShareWeights[ per ];
        }
    }
    for( CInterpRuleIterator ruleIt = mShareWeightInterpRules.begin(); ruleIt != mShareWeightInterpRules.end(); ++ruleIt ) {
        (*ruleIt)->applyInterpolations( mShareWeights, mParsedShareWeights );
    }

    // Fill in missing period vectors such as from time-step functionality.
    // All missing values are filled, whether intended or not.
    SectorUtils::fillMissingPeriodVectorInterpolated( mShareWeights );

    // All periods must have set a share weight value at this point, not having one is an error.
    for( int per = modeltime->getFinalCalibrationPeriod() + 1; per < mShareWeights.size(); ++per ) {
        if( !mShareWeights[ per ].isInited() ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Found uninitialized share weight in subsector: " << name << " in period " << per << endl;
            exit( 1 );
        }
    }
}

/*! \brief The demand passed to this function is shared out at the Technology
*          level.
* \details Variable demand (could be energy or energy service) is passed to
*          technologies and then shared out at the Technology level.
* \author Sonny Kim, Josh Lurz
* \param aSubsectorVariableDemand Total variable demand for this subsector.
* \param aFixedOutputScaleFactor Scale factor to scale down fixed output
*        technologies.
* \param aPeriod Model period
* \param aGDP Regional GDP container.
*/
void Subsector::setOutput( const double aSubsectorVariableDemand, 
                           const double aFixedOutputScaleFactor,
                           const GDP* aGDP,
                           const int aPeriod )
{
    assert( util::isValidNumber( aSubsectorVariableDemand ) && aSubsectorVariableDemand >= 0 );
    
    // Calculate the technology shares.
    const vector<double> shares = calcTechShares( aGDP, aPeriod );
    for( TechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::TechRangeIterator vintageIter = (*techIter)->getVintageBegin( aPeriod );
        
        // The first year is the current vintage, only pass variable output to current vintage.
        // Make sure that a new vintage technology exists for production.
        if( vintageIter != (*techIter)->getVintageEnd() ) {
            (*vintageIter).second->production( regionName, sectorName,
                                            aSubsectorVariableDemand * shares[ techIter - mTechContainers.begin() ],
                                            aFixedOutputScaleFactor, aGDP, aPeriod );
        }
        
        // Loop over old vintages which do not get variable demand.
        for( ++vintageIter; vintageIter != (*techIter)->getVintageEnd(); ++vintageIter ) {
            // calculate Technology output and fuel input for past vintages
            (*vintageIter).second->production( regionName, sectorName, 0,
                                            aFixedOutputScaleFactor, aGDP, aPeriod );
        }
    }
}

/*! \brief Test to see if calibration worked for this subsector
* \author Josh Lurz
* \param aPeriod The model period.
* \param aCalAccuracy Accuracy (fraction) to check if calibrations are within.
* \param aPrintWarnings Whether to print a warning.
* \return Whether calibration was successful.
*/
bool Subsector::isAllCalibrated( const int aPeriod, double aCalAccuracy, const bool aPrintWarnings ) const {
    // Check if each technology is calibrated.
    bool isAllCalibrated = true;
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        bool isThisTechCalibrated = (*techIter)->getNewVintageTechnology( aPeriod )->isAllCalibrated( aPeriod, aCalAccuracy,
                                               regionName, sectorName, name, aPrintWarnings );
 
        // if this (or any) technology not calibrated, indicate that subsector is not calibrated
        if (!isThisTechCalibrated) {
                isAllCalibrated = false;
        }
    }

    // if all technologies are calibrated, return true for the subsector, otherwise return false
    return isAllCalibrated;

}

/*! \brief returns the total calibrated output from this sector.
*
* Routine adds up calibrated values from both the sub-sector and (if not
* calibrated at Subsector), Technology levels. This returns only calibrated
* outputs, not values otherwise fixed (as fixed or zero share weights)
*
* \author Steve Smith
* \param period Model period
* \return Total calibrated output for this Subsector
*/
double Subsector::getTotalCalOutputs( const int period ) const {
    double sumCalValues = 0;
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        double currCalOutput = (*techIter)->getNewVintageTechnology( period )->getCalibrationOutput( false, "", period );
        if( currCalOutput > 0 ){
            sumCalValues += currCalOutput;
        }
    }
    /*! \post Total calibrated output is greater than or equal to zero. */
    assert( sumCalValues >= 0 );

    return sumCalValues;
}

/*! \brief returns true if all output is either fixed or calibrated.
*
* If output is is calibrated, fixed, or share weight is zero for this Subsector or all technologies in this subsector returns true.
*
* \author Steve Smith
* \param period Model period
* \return Total calibrated output for this Subsector
*/
bool Subsector::allOutputFixed( const int period ) const {
    // If there is no shareweight for this subsector than it cannot produce any
    // output, and so the output must be fixed.
    if( util::isEqual<Value>( mShareWeights[ period ], 0.0 ) ){
        return true;
    }

    // if not fixed at sub-sector level, then check at the Technology level
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        if ( !(*techIter)->getNewVintageTechnology( period )->isOutputFixed( false, "", period ) ) {
            return false;
        }
    }
    return true;
}

/*!\brief Returns a boolean for whether the subsector contains only fixed output technologies
* or at least one technology that competes on the margin.
*\author Sonny Kim
*\return Boolean for determining whether subsector contains only fixed output technologies.
*/
bool Subsector::containsOnlyFixedOutputTechnologies( const int aPeriod ) const {
    // Returns true if all technologies in the subsector has fixed exogenous outputs
    // for new investments.
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        // If at least one technology does not have fixed output return false.
        // Do not consider a technology that has a zero share weight
        const ITechnology* newVintageTech = (*techIter)->getNewVintageTechnology( aPeriod );
        if ( newVintageTech->getShareWeight() != 0 
            && !newVintageTech->isFixedOutputTechnology( aPeriod ) ) {
            return false;
        }
    }
    // Otherwise subsector contains only fixed output technologies or all technology share
    // weights are null.
    return true;
}

/*! \brief returns share weight for this Subsector
*
* Needed so that share weights can be scaled by sector
*
* \author Steve Smith
* \param period Model period
* \return share weight
*/
double Subsector::getShareWeight( const int period ) const {
    /*! \post Shareweight is valid and greater than or equal to zero. */
    assert( util::isValidNumber( mShareWeights[ period ] ) && mShareWeights[ period ] >= 0 );
    return mShareWeights[ period ];
}

//! write Subsector output to database
// TODO: Fix up this output to handle multiple vintages correctly.
void Subsector::csvOutputFile( const GDP* aGDP, 
                               const IndirectEmissionsCalculator* aIndirectEmissCalc ) const {

    // function protocol
    void fileoutput3( string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);
    
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    const string& outputUnit = mSubsectorInfo->getString( "output-unit", true );
    const string& priceUnit = mSubsectorInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    
    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    // total Subsector output
    for( int per = 0; per < maxper; ++per ){
        temp[ per ] = getOutput( per );
    }
    fileoutput3( regionName,sectorName,name," ","production",outputUnit,temp);
    // Subsector price
    for( int m = 0; m < maxper; m++ ){
        temp[ m ] = getPrice( aGDP, m );
    }
    fileoutput3( regionName,sectorName,name," ","price",priceUnit,temp);

    for ( int m= 0;m<maxper;m++){
        temp[m] = summary[m].get_emissmap_second("CO2");
    }
    fileoutput3( regionName,sectorName,name," ","CO2 emiss","MTC",temp);

    // do for all technologies in the Subsector
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        // sjs -- bad coding here, hard-wired period. But is difficult to do
        // something different with current output structure. This is just for
        // csv file. This should just use the emissions map.
        vector<string> ghgNames;
        ghgNames = mTechContainers[i]->getNewVintageTechnology( 2 )->getGHGNames();        
        for ( unsigned int ghgN =0; ghgN < ghgNames.size(); ghgN++ ) {
            if ( ghgNames[ ghgN ] != "CO2" ) {
                for ( int m=0;m<maxper;m++) {
                    temp[m] = mTechContainers[i]->getNewVintageTechnology( m )->getEmissionsByGas( ghgNames[ ghgN ], m );
                }
                fileoutput3( regionName,sectorName,name,mTechContainers[i]->getName(), ghgNames[ ghgN ] + " emiss","Tg",temp);
            }
        }

        // output or demand for each technology
        for ( int m= 0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getOutput( m );
        }
        fileoutput3( regionName,sectorName,name,mTechContainers[i]->getName(),"production",outputUnit,temp);
        // Technology share
        if( mTechContainers.size() > 1 ) {
            for ( int m=0;m<maxper;m++) {
                double subsecOutput = getOutput( m );
                if( subsecOutput > 0 ){
                    temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getOutput( m ) / subsecOutput;
                }
                else {
                    temp[ m ] = 0;
                }
            }
            fileoutput3( regionName,sectorName,name,mTechContainers[i]->getName(),"tech share","%",temp);
        }
        // Technology cost
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getCost( m );
        }
        fileoutput3( regionName,sectorName,name,mTechContainers[i]->getName(),"price",priceUnit,temp);

        // Technology CO2 emission
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getEmissionsByGas("CO2", m );
        }
        fileoutput3( regionName,sectorName,name,mTechContainers[i]->getName(),"CO2 emiss","MTC",temp);
    }
}

/*! \brief Write supply sector MiniCAM style Subsector output to database.
*
* Writes outputs with titles and units appropriate to supply sectors.
*
* \author Sonny Kim
*/
void Subsector::MCoutputSupplySector( const GDP* aGDP ) const {
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);
    
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    const double CVRT_90 = 2.212; //  convert '75 price to '90 price
    const string& outputUnit = mSubsectorInfo->getString( "output-unit", true );
    const string& priceUnit = mSubsectorInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    
    // total Subsector output
    for( int per = 0; per < maxper; ++per ){
        temp[ per ] = getOutput( per );
    }
    dboutput4(regionName,"Secondary Energy Prod",sectorName,name,outputUnit,temp);
    // Subsector price
    for( int m = 0; m < maxper; m++ ){
        temp[ m ] = getPrice( aGDP, m );
    }
    dboutput4(regionName,"Price",sectorName,name,priceUnit,temp);
    // for electricity sector only
    if (sectorName == "electricity") {
        for ( int m=0;m<maxper;m++) {
            temp[m] = getPrice( aGDP, m ) * CVRT_90 * 0.36;
        }
        dboutput4(regionName,"Price",sectorName+" C/kWh",name,"90C/kWh",temp);
    }
    
    // do for all technologies in the Subsector
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){

        // secondary energy and price output by tech
        // output or demand for each Technology
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getOutput( m );
        }

        dboutput4( regionName, "Secondary Energy Prod", sectorName + "_tech-new-investment", 
            mTechContainers[i]->getName(), outputUnit, temp );
        
        // Output for all vintages.
        for ( int m=0; m < maxper;m++) {
            temp[ m ] = 0;
            // Only sum output to the current period.
            for( int j = 0; j <= m; ++j ){
                temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getOutput( m );
            }
        }
        dboutput4( regionName, "Secondary Energy Prod", sectorName + "_tech-total", mTechContainers[i]->getName(), outputUnit, temp );
        // Technology cost
        for ( int m=0;m<maxper;m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getCost( m ) * CVRT_90;
        }
        dboutput4( regionName, "Price", sectorName + "_tech", mTechContainers[i]->getName(), "90$/GJ", temp );
    }
}

/*! \brief Write demand sector MiniCAM style Subsector output to database.
*
* Writes outputs with titles and units appropriate to demand sectors.
* Part B is for demand sector, titles and units are different from Part A
*
* \author Sonny Kim
*/
void Subsector::MCoutputDemandSector( const GDP* aGDP ) const {
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    const string& outputUnit = mSubsectorInfo->getString( "output-unit", true );
    const string& priceUnit = mSubsectorInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    
    // total Subsector output
    for( int per = 0; per < maxper; ++per ){
        temp[ per ] = getOutput( per );
    }
    dboutput4(regionName,"End-Use Service",sectorName+" by Subsec",name,outputUnit,temp);
    dboutput4(regionName,"End-Use Service",sectorName+" "+name,"zTotal",outputUnit,temp);
    // Subsector price
    for( int m = 0; m < maxper; m++ ){
        temp[ m ] = getPrice( aGDP, m );
    }
    dboutput4(regionName,"Price",sectorName,name+" Tot Cost",priceUnit,temp);
    
    // do for all technologies in the Subsector
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        if( mTechContainers.size() > 1 ) {  // write out if more than one Technology
            // output or demand for each Technology
            for ( int m=0;m<maxper;m++) {
                temp[ m ] = 0;
                for( unsigned int j = 0; j <= m; ++j ){
                    temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getOutput( m );
                }
            }
            dboutput4(regionName,"End-Use Service",sectorName+" "+name,mTechContainers[i]->getName(),outputUnit,temp);
            // total Technology cost
            for ( int m=0;m<maxper;m++) {
                temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getCost( m );
            }

            dboutput4(regionName,"Price",sectorName+" "+name,mTechContainers[i]->getName(),priceUnit,temp);
        }
    }
}

/*! \brief Write common MiniCAM style Subsector output to database.
*
* Writes outputs that are common to both supply and demand sectors.
*
* \author Sonny Kim
*/
void Subsector::MCoutputAllSectors( const GDP* aGDP,
                                    const IndirectEmissionsCalculator* aIndirectEmissCalc, 
                                    const vector<double> aSectorOutput ) const {
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    const string outputUnit = mSubsectorInfo->getString( "output-unit", true );
    const string inputUnit = mSubsectorInfo->getString( "input-unit", true );
    const string priceUnit = mSubsectorInfo->getString( "price-unit", true );
    vector<double> temp(maxper);
    
    // Subsector share
    for( int m = 0; m < maxper; m++ ){
        if( aSectorOutput[ m ] > 0 ){
            temp[ m ] = getOutput( m ) / aSectorOutput[ m ];
        }
        else {
            temp[ m ] = 0;
        }
    }
    dboutput4(regionName,"Subsec Share",sectorName,name,"100%", temp );
    
    dboutput4( regionName, "Subsec Share Wts", sectorName, name, "NoUnits", convertToVector(mShareWeights) );
    // Technology share weight by Subsector-technology for each sector
    for ( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        for ( unsigned int m = 0; m < maxper; ++m ) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getShareWeight();
        }
        dboutput4( regionName, "Tech Share Wts", sectorName, name+"-"+mTechContainers[i]->getName(),
            "NoUnits", temp );
    }

    // Subsector emissions for all greenhouse gases

    // subsector CO2 emission. How is this different then below?
    for ( int m = 0; m < maxper; m++ ) {
        temp[ m ] = 0;
        for ( unsigned int i = 0; i < mTechContainers.size(); ++i ){
            for ( unsigned int j = 0; j <  maxper; ++j ) {
                // this gives Subsector total CO2 emissions
                // get CO2 emissions for each Technology
                temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getEmissionsByGas("CO2", m );
            }
        }
    }
    dboutput4( regionName, "CO2 Emiss", sectorName, name, "MTC", temp );

    typedef map<string,double>::const_iterator CI;
    map<string,double> temissmap = summary[0].getemission(); // get gas names for period 0
    for (CI gmap=temissmap.begin(); gmap!=temissmap.end(); ++gmap) {
        for ( int m= 0;m<maxper;m++) {
            temp[m] = summary[m].get_emissmap_second(gmap->first);
        }
        dboutput4( regionName, "Emissions",  "Subsec-" + sectorName + "_" + name, gmap->first, "Tg", temp );
    }

    // do for all technologies in the Subsector
    for( unsigned int i = 0; i < mTechContainers.size(); ++i ){
        const string subsecTechName = name + mTechContainers[i]->getName();
        dboutput4(regionName,"CO2 Emiss(ind)",sectorName, subsecTechName,"MTC",temp);

        // Technology share
        for ( int m = 0; m < maxper; m++) {
            temp[ m ] = 0;
            double subsecOutput = getOutput( m );
            if( subsecOutput > 0 ){
                // sums all periods
                // does this exclude non operating vintages?
                for( unsigned int j = 0; j < maxper; ++j ){
                    temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getOutput( m ) / subsecOutput;
                }
            }
        }
        dboutput4(regionName,"Total Tech Share",sectorName, subsecTechName,"%",temp);

        // New technology share
        for ( int m = 0; m < maxper; m++) {
            temp[ m ] = 0;
            double subsecOutput = getOutput( m );
            if( subsecOutput > 0 ){
                temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getOutput( m ) / subsecOutput;
            }
        }
        dboutput4(regionName,"Tech Share (New)",sectorName, subsecTechName,"%",temp);

        // New technology share of investment.
        for ( int m=0;m<maxper;m++) {
            const vector<double> shares = calcTechShares( aGDP, m );
            temp[m] = shares[ i ];
        }

        dboutput4(regionName,"Tech Inv Share",sectorName, subsecTechName,"%",temp);

        // Old technology share
        for ( int m = 0; m < maxper; m++) {
            temp[ m ] = 0;
            double subsecOutput = getOutput( m );
            if( subsecOutput > 0 ){
                for( int j = 0; j < m; ++j ){
                    // does this exclude non-operating vintages
                    temp[m] += mTechContainers[i]->getNewVintageTechnology(j)->getOutput( m ) / subsecOutput;
                }
            }
        }
        dboutput4(regionName,"Tech Share (Old)",sectorName, subsecTechName,"%",temp);

        // New technology share of investment.
        for ( int m = 0; m < maxper; m++) {
            const vector<double> shares = calcTechShares( aGDP, m );
            temp[m] = shares[ i ];
        }
        dboutput4(regionName,"Tech Invest Share",sectorName, subsecTechName,"%",temp);

        // ghg tax and storage cost applied to Technology if any
        for ( int m = 0; m < maxper; m++) {
            temp[m] = mTechContainers[i]->getNewVintageTechnology(m)->getTotalGHGCost( regionName, sectorName, m );
        }
        dboutput4( regionName, "Total GHG Cost", sectorName, subsecTechName, priceUnit, temp);
    }
}

//! calculate GHG emissions from annual production of each Technology
void Subsector::emission( const int period ){
    /*! \pre period is less than max period. */
    assert( period < scenario->getModeltime()->getmaxper() );
    summary[period].clearemiss(); // clear emissions map
    summary[period].clearemfuelmap(); // clear emissions map
    
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::CTechRangeIterator endIter = (*techIter)->getVintageEnd();
        for( ITechnologyContainer::CTechRangeIterator vintageIter = (*techIter)->getVintageBegin( period ); vintageIter != endIter; ++vintageIter ) {
            summary[period].updateemiss( (*vintageIter).second->getEmissions( sectorName, period ) );
            summary[period].updateemfuelmap( (*vintageIter).second->getEmissionsByFuel( sectorName, period ) );
        }
    }
}

/*! \brief returns Subsector output
*
* output summed every time to ensure consistency
* this is never called for demand sectors!
*
* \author Sonny Kim, Josh Lurz
* \param period Model period
* \return sector output
*/
double Subsector::getOutput( const int period ) const {
    /*! \pre period is less than max period. */
    assert( period < scenario->getModeltime()->getmaxper() );
    double outputSum = 0;
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::CTechRangeIterator endIter = (*techIter)->getVintageEnd();
        for( ITechnologyContainer::CTechRangeIterator vintageIter = (*techIter)->getVintageBegin( period ); vintageIter != endIter; ++vintageIter ) {
            outputSum += (*vintageIter).second->getOutput( period );
        }
    }

    // Add on the base techs output too.
    for( CBaseTechIterator currTech = baseTechs.begin(); currTech != baseTechs.end(); ++currTech ){
        outputSum += (*currTech)->getOutput( period );
    }
    /*! \post Total subsector output is positive. */
    assert( outputSum >= 0 );
    return outputSum;
}

/*!
 * \brief Get the energy input for the Subsector.
 * \param aPeriod Period.
 * \return Total energy input.
 */
double Subsector::getEnergyInput( const int aPeriod ) const {
    double totalEnergy = 0;
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::CTechRangeIterator endIter = (*techIter)->getVintageEnd();
        for( ITechnologyContainer::CTechRangeIterator vintageIter = (*techIter)->getVintageBegin( aPeriod ); vintageIter != endIter; ++vintageIter ) {
            totalEnergy += (*vintageIter).second->getEnergyInput( aPeriod );
        }
    }
    return totalEnergy;
}


/*! \brief Return the total annual investment in all technologies for a given period.
* \param aPeriod Period in which to determine total investmetn.
* \return Total investment for a given period.
* \author Josh Lurz
*/
double Subsector::getAnnualInvestment( const int aPeriod ) const {
    double totalInvestment = 0;
    for( CBaseTechIterator tech = baseTechs.begin(); tech != baseTechs.end(); ++tech ){
        totalInvestment += (*tech)->getAnnualInvestment( aPeriod );
    }
    return totalInvestment;
}

/*! \brief Distribute new investment determined by the SectorInvestment object.
* \param aDistributor An object which contains the algorithm for distributing
*        investment.
* \param aNationalAccount National account object needed to calculate share.
* \param aExpProfitRateCalc An object which contains the algorithm for
*        calculating expected profits.
* \param aRegionName The name of the region containing this subsector.
* \param aSectorName The name of the sector containing this subsector.
* \param aNewInvestment The new subsector investment to be distributed.
* \param aPeriod The period in which to add investment. 
* \return The actual amount of annual investment distributed.
*/
double Subsector::distributeInvestment( const IDistributor* aDistributor,
                                        NationalAccount& aNationalAccount,
                                        const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                        const string& aRegionName,
                                        const string& aSectorName,
                                        const double aNewInvestment,
                                        const int aPeriod )
{
    // If investment is fixed used that instead of the passed in investment.
    double actInvestment = aNewInvestment;
    if( mFixedInvestments[ aPeriod ] != -1 ){
        // Check that zero investment was passed in for this case.
        if( aNewInvestment > 0 ){
            cout << "Warning: Passed in positive investment to a fixed investment subsector." << endl;
        }
        actInvestment = mFixedInvestments[ aPeriod ];
    }

    // Set the investment amount for the subsector to the quantity actually distributed.
    // Ensure that distribute() for current period is not affected by looping thru future technologies.
    vector<IInvestable*> investables = InvestmentUtils::getTechInvestables( baseTechs, aPeriod );
    mInvestments[ aPeriod ] = aDistributor->distribute( aExpProfitRateCalc,
                                                        investables,
                                                        aNationalAccount,
                                                        mTechLogitExp[ aPeriod ],
                                                        aRegionName,
                                                        aSectorName,
                                                        actInvestment,
                                                        aPeriod );

    // Check that the full amount of investment was distributed.
    if( !util::isEqual( mInvestments[ aPeriod ], actInvestment ) ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Warning: full amount of the investment was not distributed.  Actual: " << actInvestment
            << "  Distributed: " << mInvestments[ aPeriod ] << "  Difference: " 
            << fabs( mInvestments[ aPeriod ] - actInvestment ) 
            << "  in " << name << " " << aSectorName << " " << aRegionName << endl;
    }
    // Return the actual amount of investment that occurred.
    return mInvestments[ aPeriod ];
}

/*! \brief returns gets fuel consumption map for this subsector
*
* \author Sonny Kim, Josh Lurz
* \param period Model period
* \pre updateSummary
* \return fuel consumption map
*/
map<string, double> Subsector::getfuelcons( const int period ) const {
    /*! \pre period is less than max period. */
    assert( period < scenario->getModeltime()->getmaxper() );
    
    return summary[period].getfuelcons();
}

/*! \brief returns GHG emissions map for this subsector
*
* \author Sonny Kim, Josh Lurz
* \param period Model period
* \return GHG emissions map
*/
map<string, double> Subsector::getemission( const int period ) const {
    return summary[ period ].getemission();
}

/*! \brief returns map of GHG emissions by fuel for this subsector
*
* \author Sonny Kim, Josh Lurz
* \param period Model period
* \return map of GHG emissions by fuel
*/
map<string, double> Subsector::getemfuelmap( const int period ) const {
    return summary[ period ].getemfuelmap();
}

/*! \brief update summaries for reporting
*
* \author Sonny Kim, Josh Lurz
* \param aPrimaryFuelList List of primary fuels.
* \param period Model period
*/
void Subsector::updateSummary( const list<string>& aPrimaryFuelList,
                               const int period )
{
    // clears Subsector fuel consumption map
    summary[period].clearfuelcons();
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        ITechnologyContainer::CTechRangeIterator endIter = (*techIter)->getVintageEnd();
        for( ITechnologyContainer::CTechRangeIterator vintageIter = (*techIter)->getVintageBegin( period ); vintageIter != endIter; ++vintageIter ) {
            summary[period].updatefuelcons( aPrimaryFuelList, (*vintageIter).second->getFuelMap( period ) );
        }
    }
}

/*! \brief Return the  expected profit rate.
* \param aNationalAccount The regional accounting object.
* \param aRegionName The name of the region containing this subsector.
* \param aSectorName The name of the sector containing this subsector.
* \param aExpProfitRateCalc The calculator of expected profit rates.
* \param aInvestmentLogitExp The investment logit exponential.
* \param aIsShareCalc Whether this expected profit rate is being used to
*        calculate shares. Not great.
* \param aIsDistributing Whether this expected profit rate is being used
*        to distribute investment.
* \param aPeriod The period for which to get the expected profit rate.
* \return The expected profit rate for the subsector.
* \author Josh Lurz
*/
double Subsector::getExpectedProfitRate( const NationalAccount& aNationalAccount,
                                         const string& aRegionName,
                                         const string& aSectorName,
                                         const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                         const double aInvestmentLogitExp,
                                         const bool aIsShareCalc,
                                         const bool aIsDistributing,
                                         const int aPeriod ) const
{   
    assert( aExpProfitRateCalc );

    // Check for fixed investment.
    if( mFixedInvestments[ aPeriod ] != -1 && aIsDistributing ){
        return 0;
    }

    // Use the passed in expected profit calculator to determine the rate.
    return aExpProfitRateCalc->calcSectorExpectedProfitRate( 
                                         InvestmentUtils::getTechInvestables( baseTechs, aPeriod ),
                                         aNationalAccount,
                                         aRegionName,
                                         aSectorName,
                                         mTechLogitExp[ aPeriod ],
                                         aIsShareCalc,
                                         aIsDistributing,
                                         aPeriod );
}

/*! \brief Get the capital output ratio.
* \param aRegionName The name of the region containing this subsector.
* \param aSectorName The name of the sector containing this subsector.
* \param aPeriod The period.
* \return The capital output ratio.
* \author Josh Lurz
*/
double Subsector::getCapitalOutputRatio( const IDistributor* aDistributor,
                                         const IExpectedProfitRateCalculator* aExpProfitRateCalc,
                                         const NationalAccount& aNationalAccount,
                                         const string& aRegionName,
                                         const string& aSectorName, 
                                         const int aPeriod ) const
{
    assert( aDistributor );
    // Use the passed in investment distributor to calculate the share weighted average
    // capital to output ratio for the subsector.
    const double capOutputRatio = aDistributor->calcCapitalOutputRatio( 
                                    InvestmentUtils::getTechInvestables( baseTechs, aPeriod ),
                                    aExpProfitRateCalc,
                                    aNationalAccount,
                                    aRegionName,
                                    aSectorName,
                                    aPeriod );
    assert( capOutputRatio >= 0 );
    return capOutputRatio;
}
/*! \brief Operate the capital in the base technologies for this subsector. 
* \author Josh Lurz
* \param aMode Whether or not to operate all capital.
* \param aPeriod Period to operate in.
*/
void Subsector::operate( NationalAccount& aNationalAccount, const Demographic* aDemographic, const MoreSectorInfo* aMoreSectorInfo, const bool isNewVintageMode, const int aPeriod ){
    const Modeltime* modeltime = scenario->getModeltime();
    typedef vector<BaseTechnology*>::iterator BaseTechIterator;
    for( BaseTechIterator currTech = baseTechs.begin(); currTech != baseTechs.end(); ++currTech ){
        // SHK only operate for technology vintages up to current period
        if ( (*currTech)->getYear() <= modeltime->getper_to_yr( aPeriod ) ){
            (*currTech)->operate( aNationalAccount, aDemographic, aMoreSectorInfo, regionName, sectorName, isNewVintageMode, aPeriod );
        }
    }
}

/*! \brief Initialize the marketplaces in the base year to get initial demands from each Technology
 * \author Pralit Patel
 * \param period The period will most likely be the base period
 */
void Subsector::updateMarketplace( const int period ) {
    const Modeltime* modeltime = scenario->getModeltime();
    for( unsigned int j = 0; j < baseTechs.size(); j++ ) {
        if( baseTechs[ j ]->getYear() == modeltime->getper_to_yr( period ) ){ 
            baseTechs[ j ]->updateMarketplace( sectorName, regionName, period );
        }
    }
}

/*! \brief Function to finalize objects after a period is solved.
* \details This function is used to calculate and store variables which are only needed after the current
* period is complete. 
* \param aPeriod The period to finalize.
* \author Josh Lurz, Sonny Kim
*/
void Subsector::postCalc( const int aPeriod ){

    // Finalize base technologies.
    for( BaseTechIterator baseTech = baseTechs.begin(); baseTech != baseTechs.end(); ++baseTech ){
        (*baseTech)->postCalc( regionName, sectorName, aPeriod );
    }

    // Finalize all technologies in all periods.
    for( TechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->postCalc( regionName, aPeriod );
    }
}

/*! \brief For outputing SGM data to a flat csv File
 * \author Pralit Patel
 * \param period The period which we are outputing for
 */
void Subsector::csvSGMOutputFile( ostream& aFile, const int period ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    for( unsigned int j = 0; j < baseTechs.size(); j++ ) {
        if( baseTechs[ j ]->getYear() <= modeltime->getper_to_yr( period ) ){ 
            baseTechs[ j ]->csvSGMOutputFile( aFile, period );
        }
    }
}

void Subsector::accept( IVisitor* aVisitor, const int period ) const {
    aVisitor->startVisitSubsector( this, period );
    const Modeltime* modeltime = scenario->getModeltime();
    if( period == -1 ){
        // Output all techs.
        for( unsigned int j = 0; j < baseTechs.size(); j++ ) {
            baseTechs[ j ]->accept( aVisitor, period );
        }
    }
    else {
        for( unsigned int j = 0; j < baseTechs.size(); j++ ) {
            if( baseTechs[ j ]->getYear() <= modeltime->getper_to_yr( period ) ){ // should be unneeded.
                baseTechs[ j ]->accept( aVisitor, period );
            }
        }
    }
    for( CTechIterator techIter = mTechContainers.begin(); techIter != mTechContainers.end(); ++techIter ) {
        (*techIter)->accept( aVisitor, period );
    }
            
    aVisitor->endVisitSubsector( this, period );
}

/*! \brief Return fixed investment.
* \param aPeriod Period to return fixed investment for.
* \return Subsector level fixed investment.
* \todo Can't have a zero investment subsector.
*/
double Subsector::getFixedInvestment( const int aPeriod ) const {
    // Return amount of subsector fixed investment.
    if( mFixedInvestments[ aPeriod ] != -1 ){
        return mFixedInvestments[ aPeriod ];
    }
    // Sum any fixed investment at the vintage level.
    double totalFixedTechInvestment = InvestmentUtils::sumFixedInvestment( 
                                      InvestmentUtils::getTechInvestables( baseTechs, aPeriod ), aPeriod );
    
    /*! \post Fixed investment must be positive */
    assert( totalFixedTechInvestment >= 0 );
    return totalFixedTechInvestment;
}

bool Subsector::hasCalibrationMarket() const {
    return doCalibration;
}
