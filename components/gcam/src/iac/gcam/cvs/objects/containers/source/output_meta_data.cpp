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
* \file output_meta_data.cpp
* \ingroup Objects
* \brief The OutputMetaData class source file
* \author Josh Lurz
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/output_meta_data.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

//! Constructor
OutputMetaData::OutputMetaData() {
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
const std::string& OutputMetaData::getXMLNameStatic() {
    const static string XML_NAME = "output-meta-data";
    return XML_NAME;
}

/*! \brief Write out XML data for input.
* \param aOut Output stream.
* \param aTabs Tabs object.
*/
void OutputMetaData::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag ( getXMLNameStatic(), aOut, aTabs );
    typedef list<string>::const_iterator CListIterator;
    map<string, string> attrs;
    for( CListIterator iter = mPrimaryFuels.begin(); iter != mPrimaryFuels.end(); ++iter ) {
        attrs[ "var" ] = *iter;
        XMLWriteElementWithAttributes( "", "primary-fuel", aOut, aTabs, attrs );
    }
    for( CListIterator iter = mSummableVariables.begin(); iter != mSummableVariables.end(); ++iter ) {
        attrs[ "var" ] = *iter;
        XMLWriteElementWithAttributes( "", "summable", aOut, aTabs, attrs );
    }
    for( CListIterator iter = mHasYearVariables.begin(); iter != mHasYearVariables.end(); ++iter ) {
        attrs[ "var" ] = *iter;
        XMLWriteElementWithAttributes( "", "has-year", aOut, aTabs, attrs );
    }
    XMLWriteElement( mScenarioSummary, "summary", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Parse the meta-data from XML.
* \details The model does not use this meta-data internally but reads it from
*          the input XML file here and passes that information along to both the
*          database output xml and the derived xml input file.
* \param aNode Root node of the object's DOM subtree.
*/
void OutputMetaData::XMLParse( const DOMNode* aNode ) {
    /*! \pre make sure we were passed a valid node. */
    assert( aNode );

    // get all child nodes.
    const DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        if( curr->getNodeType() != DOMNode::ELEMENT_NODE ){
            continue;
        }
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
        if( nodeName == "primary-fuel" ){
            mPrimaryFuels.push_back( XMLHelper<string>::getAttr( curr, "var" ) );
        }
        else if( nodeName == "summable" ){
            mSummableVariables.push_back( XMLHelper<string>::getAttr( curr, "var" ) );
        }
        else if( nodeName == "has-year" ){
            mHasYearVariables.push_back( XMLHelper<string>::getAttr( curr, "var" ) );
        }
        else if ( nodeName == "summary" ){
            mScenarioSummary = XMLHelper<string>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }
}

void OutputMetaData::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitOutputMetaData( this, aPeriod );
    aVisitor->endVisitOutputMetaData( this, aPeriod );
}

//! Get the primary fuel list. Remove this function once the output database is removed.
const list<string>& OutputMetaData::getPrimaryFuelList() const {
    return mPrimaryFuels;
}
