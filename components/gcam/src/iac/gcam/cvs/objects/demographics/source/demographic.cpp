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
* \file demographic.cpp
* \ingroup Objects-SGM
* \brief Demographic class source file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include "util/base/include/definitions.h"
#include <vector>
#include <map>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "demographics/include/demographic.h"
#include "demographics/include/population.h"
#include "demographics/include/population_sgm_fixed.h"
#include "demographics/include/population_sgm_rate.h"
#include "demographics/include/population_mini_cam.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"
// class for reporting
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default constructor.
Demographic::Demographic() {
}

//! Demographic destructor. 
Demographic::~Demographic(){
    clear();
}

//! Helper member function for the destructor. Performs memory deallocation. 
void Demographic::clear(){
    for( PopulationIterator popIter = population.begin(); popIter != population.end(); ++popIter ){
        delete *popIter;
    }
}

//! parses Demographics xml object
void Demographic::XMLParse( const xercesc::DOMNode* node ){
    // make sure we were passed a valid node.
    assert( node );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == PopulationMiniCAM::getXMLNameStatic() ) {
            parseContainerNode( curr, population, yearToMapIndex, new PopulationMiniCAM(), "year" );
        }
        else if( nodeName == PopulationSGMFixed::getXMLNameStatic() ){
            parseContainerNode( curr, population, yearToMapIndex, new PopulationSGMFixed(), "year" );
        }
        else if( nodeName == PopulationSGMRate::getXMLNameStatic() ){
            parseContainerNode( curr, population, yearToMapIndex, new PopulationSGMRate(), "year" );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing demographics." << endl;
        }
    }
}

//! Write out data members to XML output stream.
void Demographic::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag ( getXMLName(), out, tabs );

    for( CPopulationIterator i = population.begin(); i != population.end(); ++i ){
        if( ( *i )->mIsParsed ){
            ( *i )->toInputXML( out, tabs );
        }
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Write out XML for debugging purposes.
void Demographic::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag ( getXMLName(), out, tabs );
    // Convert the period into an index into the demographics object.
    int index = convertPeriodToPopulationIndex( period );
    // Check if the conversion into an index failed.
    if( index != -1 ){
        population[ index ]->toDebugXML( out, tabs );
    }
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Complete the initialization.
void Demographic::completeInit(){
    // First make sure we have a population for each model period.  If we do not have one in a given
    // period interpolate between ones we do have.
    const Modeltime* modeltime = scenario->getModeltime();
    Population* prevPopulation = 0;
    for( int period = 0; period < modeltime->getmaxper(); ++period ) {
        int modelYear = modeltime->getper_to_yr( period );
        CYearMapIterator iter = yearToMapIndex.find( util::toString( modelYear ) );
        if( iter == yearToMapIndex.end() ) {
            // the user must have read in the population for the first model period
            assert( period > 0 );
            
            iter = ++yearToMapIndex.find( util::toString( prevPopulation->getYear() ) );
            // the user must have read in the population for the final model period
            assert( iter != yearToMapIndex.end() );
            prevPopulation = prevPopulation->cloneAndInterpolate( modelYear, population[ (*iter).second ] );
            population.push_back( prevPopulation );
            yearToMapIndex[ util::toString( modelYear ) ] = population.size() - 1;
        }
        else {
            prevPopulation = population[ (*iter).second ];
        }
    }

    // Sort the population vector by year after cloning and interpolating.
    // If population is not interpolated sort is not required.
    YearComparator comp;
    sort( population.begin(), population.end(), comp );

    // remap year to map index after sorting
    int index = 0;
    for( CPopulationIterator i = population.begin(); i != population.end(); ++i ){
        yearToMapIndex[ util::toString( ( *i )->getYear() ) ] = index;
        ++index;
    }
    
    for( PopulationIterator popIter = population.begin(); popIter != population.end(); ++popIter ) {
        if( popIter == population.begin() ){
            ( *popIter )->completeInit();
        }
        else {
            // Initialize each population with the surviving population from the previous period.
            ( *popIter )->completeInit( (*(popIter - 1))->getSurvFemalePop(), (*( popIter - 1 ))->getSurvMalePop()  );
        }
    }
}

//! initialize anything that won't change during the calcuation
void Demographic::initCalc(){
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overriden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& Demographic::getXMLName() const {
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
const string& Demographic::getXMLNameStatic(){
    const static string XML_NAME = "demographics";
    return XML_NAME;
}

//! return total population
double Demographic::getTotal( const int per ) const {
    // Convert the period into an index into the demographics object.
    int index = convertPeriodToPopulationIndex( per );
    // Check if the conversion into an index failed.
    if( index == -1 ){
        return 0;
    }
    return population[ index ]->getTotal();
}

//! return the male working age population
double Demographic::getWorkingAgePopulationMales( const int per ) const {
    // Convert the period into an index into the demographics object.
    int index = convertPeriodToPopulationIndex( per );
    // Check if the conversion into an index failed.
    if( index == -1 ){
        return 0;
    }
    return population[ index ]->getWorkingAgePopMale();
}

//! return the female working age population
double Demographic::getWorkingAgePopulationFemales( const int per ) const {
    // Convert the period into an index into the demographics object.
    int index = convertPeriodToPopulationIndex( per );
    // Check if the conversion into an index failed.
    if( index == -1 ){
        return 0;
    }
    return population[ index ]->getWorkingAgePopFemale();
}

//! return total working age population (male and female)
double Demographic::getWorkingAgePopulation( const int per ) const {
    // Convert the period into an index into the demographics object.
    int index = convertPeriodToPopulationIndex( per );
    // Check if the conversion into an index failed.
    if( index == -1 ){
        return 0;
    }
    return population[ index  ]->getWorkingAgePop();
}

//! Translate a period into the index within the demographic object of the population.
int Demographic::convertPeriodToPopulationIndex( int aPeriod ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    assert( aPeriod >= 0 && aPeriod < modeltime->getmaxper() );

    int year = modeltime->getper_to_yr( aPeriod ); // get year from model period

    // print out error if year doesn't exist
    CYearMapIterator iter = yearToMapIndex.find( util::toString( year ) );
    if( iter == yearToMapIndex.end() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "In convertPeriodToPopulationIndex, Year " << year 
                << " corresponding to period " << aPeriod << " doesn't exist." << endl;
        return -1;
    }
    return iter->second;
}

/*! \brief Return a vector of total population.
* \return A vector of population data.
* \note This is for the Fortran AgLU only and should be removed when that
*       component is no longer used.
*/
const vector<double> Demographic::getTotalPopVec() const {
    // Create a vector with one slot per model period. Add an extra
    // slot for the population period before the base period.
    const Modeltime* modeltime = scenario->getModeltime();
    vector<double> newTotalVector( modeltime->getmaxper() + 1 );
    
    // Add the extra earlier period's value. The year for this will be base year
    // minus timestep. This can't be done with the helper function because 
    // modeltime will not handle period -1.
    int prevYear = modeltime->getStartYear() - modeltime->gettimestep( 0 );
    CYearMapIterator iter = yearToMapIndex.find( util::toString( prevYear ) );

    // If the year was not found print a warning and leave the value as zero.
    if( iter == yearToMapIndex.end() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::ERROR );
        mainLog << "There is no population object corresponding to " << prevYear << endl;
    }
    else {
        // Set the initial spot.
        newTotalVector[ 0 ] = population[ iter->second ]->getTotal();
    }

    // Fill in the vector with the total population for each period.
    for ( int i = 0; i < modeltime->getmaxper(); ++i ){
        int index = convertPeriodToPopulationIndex( i );
        // Check for invalid indices.
        if( index != -1 ){
            newTotalVector[ i + 1 ] = population[ index ]->getTotal();
        }
    }
    return newTotalVector;
}


//! MiniCAM output to file
void Demographic::dbOutput( const string& regionName ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxPeriod = modeltime->getmaxper();
    vector<double> temp( maxPeriod );

    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // write population to temporary array since not all will be sent to output
    for ( int i = 0; i < maxPeriod; i++ ){
        int index = convertPeriodToPopulationIndex( i );
        // Check for invalid indices.
        if( index != -1 ){
            temp[ i ] = population[ index ]->getTotal();
        }
    }
    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    dboutput4( regionName, "General", "Population", "Total", "thous", temp );
}

//! outputing population info to file
void Demographic::csvOutputFile( const string& regionName ) const {
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxPeriod = modeltime->getmaxper();
    vector<double> temp( maxPeriod );

    // function protocol
    void fileoutput3( string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    // write population to temporary array since not all will be sent to output
    for ( int i = 0; i < maxPeriod; i++ ){
        int index = convertPeriodToPopulationIndex( i );
        // Check for invalid indices.
        if( index != -1 ){
            temp[i] = population[ index ]->getTotal();
        }
    }

    // function arguments are variable name, double array, db name, table name
    // the function writes all years
    fileoutput3( regionName," "," "," ","population","1000s",temp);
}

void Demographic::csvSGMOutputFile( ostream& aFile, const int period ) const {
    aFile << "Demographic Data for Labor Force and Government Transfers" << endl << endl;
    aFile << getTotal( period ) << ',' << "Total Population" << endl;
    aFile << getWorkingAgePopulationMales( period ) << ',' << "Working Age Pop. Male" << endl;
    aFile << getWorkingAgePopulationFemales( period ) << ',' << "Working Age Pop. Females" << endl;
    aFile << endl;

    int index = convertPeriodToPopulationIndex( period );
    // Check for invalid indices.
    if( index != -1 ){
        population[ index ]->csvSGMOutputFile( aFile, period );
    }
}

// for reporting
void Demographic::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitDemographic( this, aPeriod );

    // If the period is -1 update the output container with information about all the populations.
    if( aPeriod == -1 ){
        for( unsigned int i = 0; i < population.size(); ++i ){
            population[ i ]->accept( aVisitor, aPeriod );
        }
    }
    // Otherwise only update the one for the given period.
    else {
    int index = convertPeriodToPopulationIndex( aPeriod );
        // Check for invalid indices.
        if( index != -1 ){
            population[ index ]->accept( aVisitor, aPeriod );
        }
    }
    aVisitor->endVisitDemographic( this, aPeriod );
}
