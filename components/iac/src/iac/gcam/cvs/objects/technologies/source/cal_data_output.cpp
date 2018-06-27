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
 * \file cal_data_output.cpp
 * \ingroup Objects
 * \brief CalDataOutput class source file.
 * \author James Blackwood
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNodeList.hpp>
#include "technologies/include/cal_data_output.h"

#include "util/base/include/xml_helper.h"

using namespace std;
using namespace xercesc;

/*! \brief Constructor.
* \author James Blackwood
*/
CalDataOutput::CalDataOutput() {
    mCalOutputValue = 0;
}

/*! \brief Clone the current object.
* \return A clone of the object.
*/
CalDataOutput* CalDataOutput::clone() const {
    return new CalDataOutput( *this );
}

/*! \brief Parses XML for the object.
* \author James Blackwood
* \param aNode pointer to the current node in the XML input tree
*/
void CalDataOutput::XMLParse( const DOMNode* aNode ){
	// assume we are passed a valid node.
	assert( aNode );

	// get all the children.
	DOMNodeList* nodeList = aNode->getChildNodes();

	for( unsigned int i = 0;  i <  nodeList->getLength(); ++i ){
		const DOMNode* curr = nodeList->item( i );
		const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

		if( nodeName == "#text" ) {
			continue;
		}
        else if( nodeName == "calOutputValue" ) {
            mCalOutputValue = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
	        mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << endl;
		}
	}
}

//! write object to xml output stream
void CalDataOutput::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mCalOutputValue, "calOutputValue", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

//! Write object to debugging xml output stream.
void CalDataOutput::toDebugXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mCalOutputValue, "calOutputValue", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
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
const std::string& CalDataOutput::getXMLNameStatic() {
    const static string XML_NAME = "CalDataOutput";
    return XML_NAME;
}

void CalDataOutput::initCalc( const Demographic* aDemographics, const int aPeriod ) {
}

void CalDataOutput::completeInit(){
}

double CalDataOutput::getCalOutput() {
    return mCalOutputValue;
}
