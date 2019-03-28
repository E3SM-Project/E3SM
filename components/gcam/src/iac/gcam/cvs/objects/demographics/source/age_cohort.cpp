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
* \file age_cohort.cpp
* \ingroup Objects-SGM
* \brief AgeCohort class source file.
* \author Sonny Kim
* \author Katherine Chung
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <string>

// Needed for int->string conversions.
#include <sstream>
#include <algorithm>

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "demographics/include/age_cohort.h"
#include "demographics/include/gender.h"
#include "demographics/include/male.h"
#include "demographics/include/female.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern ofstream outputFile;

//! default constructor
AgeCohort::AgeCohort() {
    mLowerAgeBound = -1;
    mUpperAgeBound = -1;
}

//! Default Destructor. Needed because we are using auto_ptr so that we avoid
//! destroying incomplete types.
AgeCohort::~AgeCohort(){
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
const std::string& AgeCohort::getXMLNameStatic(){
    const static string XML_NAME = "ageCohort";
    return XML_NAME;
}

//! parses AgeCohort xml object
void AgeCohort::XMLParse( const xercesc::DOMNode* node ){
    // make sure we were passed a valid node.
    assert( node );

    ageGroup = XMLHelper<string>::getAttr( node, "ageGroup" );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        // Note: Male and female can either be parsed as <gender ="male"> or <male>
        else if( nodeName == "gender" ){
            if( !parseGender( curr ) ){
                cout << "Failed to parse a gender element." << endl;
            }
        }
        else if( parseGender( curr ) ){
            // do nothing but don't print a warning.
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }
}

//! Attempt to parse a gender.
bool AgeCohort::parseGender( DOMNode* aNode ) {
    // We could either be parsing out of the "type" attribute
    // or the element name itself.
    string type = XMLHelper<string>::getAttr( aNode, "type" );
    if( type == "" ){
        type = XMLHelper<string>::safeTranscode( aNode->getNodeName() );
    }

    if( type == Male::getXMLNameStatic() ){
        if( !male.get() ){
            male.reset( new Male() );
        }
        male->XMLParse( aNode );
    }
    else if( type == Female::getXMLNameStatic() ) {
        if( !female.get() ){
            female.reset( new Female() );
        }
        female->XMLParse( aNode );
    }
    else {
        return false;
    }
    return true;
}
//! Write out data members to XML output stream.
void AgeCohort::toInputXML( std::ostream& out, Tabs* tabs ) const {
    // write the beginning tag.
    tabs->writeTabs( out );
    out << "<" << getXMLNameStatic() << " ageGroup=\"" << ageGroup << "\">"<< endl;

    // increase the indent.
    tabs->increaseIndent();

    male->toInputXML( out, tabs );
    female->toInputXML( out, tabs );
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

//! Write out XML for debugging purposes.
void AgeCohort::toDebugXML( ostream& out, Tabs* tabs ) const {
    map<string, string> attrMap;
    attrMap[ "ageGroup" ] = ageGroup;
    // write the beginning tag.
    tabs->writeTabs( out );
    out << "<" << getXMLNameStatic() << " ageGroup=\"" << ageGroup << "\">"<< endl;

    // increase the indent.
    tabs->increaseIndent();

    male->toDebugXML( out, tabs );
    female->toDebugXML( out, tabs );

    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

//! Complete the initialization.
void AgeCohort::completeInit(){
    // Check for a not completely initialized cohort.
    if( !male.get() ){
        // Warn the user and create an empty cohort.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Creating empty male group for cohort " << getAgeGroup() << "." << endl;
        male.reset( new Male() );
        male->setPopulation( 0 );
    }
    if( !female.get() ){
        // Warn the user and create an empty cohort.
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::NOTICE );
        mainLog << "Creating empty female group for cohort " << getAgeGroup() << "." << endl;
        female.reset( new Female() );
        female->setPopulation( 0 );
    }
}

//! initialize anything that won't change during the calcuation
void AgeCohort::initCalc(){
    male->calcSurvivingPop();
    female->calcSurvivingPop();
}

//! sets the male population to the param passed in
void AgeCohort::setMalePop( double aMalePopulation ){
    male->setPopulation( aMalePopulation );
}

//! sets the female population to the param passed in
void AgeCohort::setFemalePop(double aFemalePopulation ){
    female->setPopulation( aFemalePopulation );
}

//! returns the male population for this agecohort
double AgeCohort::getMalePop() const {
    return male->getPopulation();
}

//! returns the female population for this agecohort
double AgeCohort::getFemalePop() const {
    return female->getPopulation();
}

// calculate the male births from this age group
// and return the value
// only females give birth
double AgeCohort::calcMaleBirth() {
    return female->calcMaleBirth();
}

// calculate the female births from this age group
// and return the value
// only females give birth
double AgeCohort::calcFemaleBirth() {
    return female->calcFemaleBirth();
}

// calculate the surviving male population of this age group
// and return the value
double AgeCohort::calcSurvMalePop() {
    return male->calcSurvivingPop();
}

// calculate the surviving male population of this age group
// and return the value
double AgeCohort::calcSurvFemalePop() {
    return female->calcSurvivingPop();
}

const string& AgeCohort::getAgeGroup() const {
    return ageGroup;
}

/*! \brief Get the lower limit to the ages contained in this cohort.
* \details This calculates and returns the lower limit of the ages stored
* in this cohort. The function also caches the result locally so that the 
* calculation does not have to be repeated each iteration. If the age limit 
* cannot be parsed the function returns -1.
* \return The lower age limit for this cohort and -1 if there is an error.
* \author Josh Lurz
*/
int AgeCohort::getLowerAgeBound() const {
    // Check if we cached the lower age bound.
    if( mLowerAgeBound != -1 ){
        return mLowerAgeBound;
    }

    // Check if this is the last age cohort, specified as +x
    string lowerRangeString;
    if( ageGroup.find_first_of( '+' ) != ageGroup.npos ){
        lowerRangeString = ageGroup.substr( 1, ageGroup.size() - 1 );
    }
    else {
        // Check and make sure there is one and only one hyphen.
        if( count( ageGroup.begin(), ageGroup.end(), '-' ) != 1 ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Invalid age range " << ageGroup << " found while calculating getLowerAgeBound " << endl;
            return -1;
        }
        // Find the hyphen we now know is there.
        const string::size_type hyphenPos = ageGroup.find_first_of( '-' );
        // Get the lower range string.
        lowerRangeString = ageGroup.substr( 0, hyphenPos );
    }

    // Convert it to an int.
    istringstream converter( lowerRangeString );
    int lowerRange = -1;
    converter >> lowerRange;

    // Cache the value since this function is slow.
    mLowerAgeBound = lowerRange;
    return lowerRange;
}

/*! \brief Get the upper limit to the ages contained in this cohort.
* \details This calculates and returns the upper limit of the ages stored
* in this cohort. The function also caches the result locally so that the 
* calculation does not have to be repeated each iteration. If the age limit 
* cannot be parsed the function returns -1. If the upper age range is unbounded,
* in the case where this is the last cohort, it will return AGE_MAX.
* \return The upper age limit for this cohort, AGE_MAX if it is unbounded and -1 if
* there is an error.
* \author Josh Lurz
*/
int AgeCohort::getUpperAgeBound() const {
    const int AGE_MAX = 500;

    // Check if we cached the upper age bound.
    if( mUpperAgeBound != -1 ){
        return mUpperAgeBound;
    }
    // Check if this is the last age cohort, specified as +x
    if( ageGroup.find_first_of( '+' ) != ageGroup.npos ){
        // This age group has no upper bound.
        return AGE_MAX;
    }

    // Check and make sure there is one and only one hyphen.
    if( count( ageGroup.begin(), ageGroup.end(), '-' ) != 1 ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Invalid age range " << ageGroup << " found while calculating getUpperAgeBound " << endl;
        return -1;
    }
    // Find the hyphen we now know is there.
    const string::size_type hyphenPos = ageGroup.find_first_of( '-' );
    // Get the upper range string.
    const string upperRangeString = ageGroup.substr( hyphenPos + 1, ageGroup.size() -1 );
    // Convert it to an int.
    istringstream converter( upperRangeString );
    int upperRange = -1;
    converter >> upperRange;

    // Cache the value since this function is slow.
    mUpperAgeBound = upperRange;
    return upperRange;
}

/*! \brief Update a visitor with information about a AgeCohort.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to update.
*/
void AgeCohort::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitAgeCohort( this, aPeriod );
    male->accept( aVisitor, aPeriod );
    female->accept( aVisitor, aPeriod );
    aVisitor->endVisitAgeCohort( this, aPeriod );
}
