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
 * \file standard_technical_change_calc.cpp
 * \ingroup Objects
 * \brief StandardTechnicalChangeCalc source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "technologies/include/standard_technical_change_calc.h"
#include "util/logger/include/ilogger.h"
#include "functions/include/iinput.h"
#include "functions/include/ifunction.h"
#include "technologies/include/technology.h" // for PrevPeriodInfo.

using namespace std;

extern Scenario* scenario;

//! Constructor
StandardTechnicalChangeCalc::StandardTechnicalChangeCalc()
{
}


StandardTechnicalChangeCalc* StandardTechnicalChangeCalc::clone() const {
    return new StandardTechnicalChangeCalc( *this );
}

bool StandardTechnicalChangeCalc::isSameType( const string& aType ) const {
	return aType == getXMLNameStatic();
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
* \details This public function accesses the private constant string, XML_NAME.
*          This way the tag is always consistent for both read-in and output and
*          can be easily changed. The "==" operator that is used when parsing,
*          required this second function to return static.
* \note A function cannot be static and virtual.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& StandardTechnicalChangeCalc::getXMLNameStatic() {
    const static string XML_NAME = "tech-change";
    return XML_NAME;
}

const string& StandardTechnicalChangeCalc::getName() const {
    return getXMLNameStatic();
}

// Documentation inherits.
bool StandardTechnicalChangeCalc::XMLParse( const xercesc::DOMNode* node ){
	/*! \pre Assume we are passed a valid node. */
	assert( node );

	const xercesc::DOMNodeList* nodeList = node->getChildNodes();
	for( unsigned int i = 0; i < nodeList->getLength(); i++ ) {
		const xercesc::DOMNode* curr = nodeList->item( i );
		if( curr->getNodeType() != xercesc::DOMNode::ELEMENT_NODE ){
			continue;
		}
		const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );
		if( nodeName == "hicks-neutral" ){
			mHicksNeutralTechChange = XMLHelper<double>::getValue( curr );
		}
	    else if( nodeName == "energy-only" ){
            mEnergyTechChange = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "material-only" ){
            mMaterialTechChange = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::ERROR );
            mainLog << "Unknown tag " << nodeName << " encountered while processing " << getXMLNameStatic() << endl;
        }
	}

    // TODO: Handle success and failure better.
    return true;
}

void StandardTechnicalChangeCalc::toInputXML( ostream& aOut,
                                              Tabs* aTabs ) const
{
	XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
	XMLWriteElementCheckDefault( mHicksNeutralTechChange, "hicks-neutral", aOut, aTabs, Value( 0 ) );
    XMLWriteElementCheckDefault( mEnergyTechChange, "energy-only", aOut, aTabs, Value( 0 ) );
    XMLWriteElementCheckDefault( mMaterialTechChange, "material-only", aOut, aTabs, Value( 0 ) );
	XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void StandardTechnicalChangeCalc::toDebugXML( const int aPeriod,
                                              ostream& aOut,
                                              Tabs* aTabs ) const
{
	XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
	XMLWriteElement( mHicksNeutralTechChange, "hicks-neutral", aOut, aTabs );
    XMLWriteElement( mEnergyTechChange, "energy-only", aOut, aTabs );
    XMLWriteElement( mMaterialTechChange, "material-only", aOut, aTabs );
	XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void StandardTechnicalChangeCalc::completeInit()
{
}

double StandardTechnicalChangeCalc::calcAndAdjustForTechChange( vector<IInput*>& aInputs,
                                                                PreviousPeriodInfo& aPreviousPeriodInfo,
                                                                const IFunction* aProductionFunc,
                                                                const string& aRegionName,
                                                                const string& aSectorName,
                                                                const int aPeriod ) const
{
    // Allow the production function to adjust for technical change given the
    // read-in energy and material technical change. Note that the inputs have
    // already been adjusted for previous technical change by copying the
    // coeffients forward.
    TechChange techChange( mMaterialTechChange, mEnergyTechChange, mHicksNeutralTechChange );
    return aProductionFunc->applyTechnicalChange( aInputs, techChange, aRegionName,
                                                  aSectorName, aPeriod,
                                                  aPreviousPeriodInfo.mCumulativeHicksNeutralTechChange );
}
