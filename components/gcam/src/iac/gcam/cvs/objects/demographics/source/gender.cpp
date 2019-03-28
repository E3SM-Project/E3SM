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
* \file gender.cpp
* \ingroup Objects-SGM
* \brief world class source file.
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "demographics/include/gender.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

// static initialization
const string Gender::XML_NAME = "gender";

//! Default constructor
Gender::Gender(){
    mPopulation = -1;
    mSurvivingPop = 0;
    mSurvivalRate = 0;
}

//! parse SOME xml data for Male or Female objects
void Gender::XMLParse( const DOMNode* node ) {
    /*! \pre make sure we were passed a valid node. */
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
        else if ( nodeName == "population" ) {
            mPopulation = XMLHelper<double>::getValue( curr );
        }
        else if (nodeName == "survivalRate" ) {
            mSurvivalRate = XMLHelper<double>::getValue( curr );
        }
        else if( !XMLDerivedClassParse( nodeName, curr ) ){
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string " << nodeName << " encountered while parsing " << getXMLName() << endl;
        }
    }
}

//! Output to XML data
void Gender::toInputXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs, "", 0, getXMLName() );

    XMLWriteElementCheckDefault( mPopulation, "population", out, tabs );
    XMLWriteElementCheckDefault( mSurvivingPop, "survivingPop", out, tabs );
    XMLWriteElementCheckDefault( mSurvivalRate, "survivalRate", out, tabs );

    // write out variables for derived classes
    toInputXMLDerived( out, tabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

//! Output debug info to XML data
void Gender::toDebugXML( ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag ( getXMLName(), out, tabs );

    XMLWriteElement( mPopulation, "population", out, tabs );
    XMLWriteElement( mSurvivingPop, "survivingPop", out, tabs );
    XMLWriteElement( mSurvivalRate, "survivalRate", out, tabs );

    // write out variables for derived classes
    toDebugXMLDerived( out, tabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

// calculate the suviving population which is 
// passed on to the next age group and time
double Gender::calcSurvivingPop() {
    assert( mPopulation != -1 );
    mSurvivingPop = mPopulation * mSurvivalRate;
    return mSurvivingPop;
}

//! Returns the population
double Gender::getPopulation() const {
    assert( mPopulation != -1 );
    return mPopulation;
}

//! Sets the population to the param passed in
void Gender::setPopulation( double aPopulation ){
    mPopulation = aPopulation;
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
const std::string& Gender::getXMLNameStatic(){
    return XML_NAME;
}

/*! \brief Update a visitor with information about a Gender.
* \param aVisitor Visitor to update.
* \param aPeriod Period for which to update.
*/
void Gender::accept( IVisitor* aVisitor, const int aPeriod ) const {
	aVisitor->startVisitGender( this, aPeriod );
	aVisitor->endVisitGender( this, aPeriod );
}
