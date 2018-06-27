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
* \file land_carbon_densities.cpp
* \ingroup Objects
* \brief LandCarbonDensities class source file.
* \author James Blackwood
*/

#include "util/base/include/definitions.h"
#include "util/base/include/xml_helper.h"
#include <xercesc/dom/DOMNodeList.hpp>
#include "ccarbon_model/include/land_carbon_densities.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "ccarbon_model/include/carbon_model_utils.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! \brief Constructor.
* \author James Blackwood
*/
LandCarbonDensities::LandCarbonDensities():
mAboveGroundCarbon( scenario->getModeltime()->getStartYear(), CarbonModelUtils::getEndYear(), 0.0 ),
mBelowGroundCarbon( scenario->getModeltime()->getStartYear(), CarbonModelUtils::getEndYear(), 0.0 ),
mAvgAboveGroundCarbon( 0.0 ),
mAvgBelowGroundCarbon( 0.0 ),
mMatureAge( 1 ),
mAboveGroundAdj( scenario->getModeltime()->getStartYear(), CarbonModelUtils::getEndYear(), 1.0 ),
mBelowGroundAdj( scenario->getModeltime()->getStartYear(), CarbonModelUtils::getEndYear(), 1.0 )
{
}

//! Default destructor
LandCarbonDensities::~LandCarbonDensities() {
}

LandCarbonDensities* LandCarbonDensities::clone() const{
    return 0;
}

/*! \brief Parses all XML data for the class.
* \author James Blackwood
* \param aCurr Pointer to the current node in the XML input tree
*/
bool LandCarbonDensities::XMLParse( const DOMNode* aCurr ) {
    // Assume we are passed a valid node.
    assert( aCurr );

    // get all the children.
    DOMNodeList* nodeList = aCurr->getChildNodes();

    //    ILogger& mainLog = ILogger::getLogger( "main_log" );
    //    mainLog.setLevel( ILogger::WARNING );
    for( unsigned int i = 0;  i < nodeList->getLength(); ++i ){
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "above-ground-carbon-density" ) {
            mAvgAboveGroundCarbon = XMLHelper<double>::getValue( curr );
	    //            mainLog << "xmlparse getting mAvgAboveGroundCarbon " <<  mAvgAboveGroundCarbon << endl;
        }
        else if( nodeName == "below-ground-carbon-density" ) {
            mAvgBelowGroundCarbon = XMLHelper<double>::getValue( curr );
	    //            mainLog << "xmlparse getting mAvgBelowGroundCarbon " <<  mAvgBelowGroundCarbon << endl;
        }
        else if( nodeName == "mature-age" ) {    
			// Because the mature age for some land areas is read in 
			// as a double, we are converting it to an integer
			// during XML parse.
            double temp = XMLHelper<double>::getValue( curr );
			
			mMatureAge = util::round( temp );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string " << nodeName
                    << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }

    // TODO: Improve error handling.
    return true;
}

void LandCarbonDensities::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    const Modeltime* modeltime = scenario->getModeltime();
    const int year = modeltime->getper_to_yr( aPeriod );
    XMLWriteElement( mAvgAboveGroundCarbon, "above-ground-carbon-density", aOut, aTabs );
    XMLWriteElement( mAvgBelowGroundCarbon, "below-ground-carbon-density", aOut, aTabs );
    XMLWriteElement( mLandUse[ aPeriod ], "land-use", aOut, aTabs );
    XMLWriteElement( mTotalEmissions[ year ], "total-emissions", aOut, aTabs );
    XMLWriteElement( mMatureAge, "mature-age", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void LandCarbonDensities::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs );
    XMLWriteElement( mAvgAboveGroundCarbon, "above-ground-carbon-density", aOut, aTabs, 0.0 );
    XMLWriteElement( mAvgBelowGroundCarbon, "below-ground-carbon-density", aOut, aTabs, 0.0 );
    XMLWriteElement( mMatureAge, "mature-age", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author James Blackwood
* \return The constant XML_NAME as a static.
*/
const string& LandCarbonDensities::getXMLNameStatic() {
    const static string XML_NAME = "land-carbon-densities";
    return XML_NAME;
}

/*! 
* \brief Perform initializations that only need to be done once.
* \author Kate Calvin
*/
void LandCarbonDensities::completeInit() {
    for ( int i = scenario->getModeltime()->getStartYear(); i <= CarbonModelUtils::getEndYear(); ++i ){
        mAboveGroundCarbon[ i ] = mAvgAboveGroundCarbon;
        mBelowGroundCarbon[ i ] = mAvgBelowGroundCarbon;
    }
    
    // force the sigmoid to get precalculated.
    setMatureAge( mMatureAge );
}

void LandCarbonDensities::setActualAboveGroundCarbonDensity( const double aAboveGroundCarbonDensity,
                                                     const int aYear )
{
  if (getActualAboveGroundCarbonDensity( aYear ) > 0)
    mAboveGroundAdj[ aYear ] = aAboveGroundCarbonDensity / getActualAboveGroundCarbonDensity( aYear );
}

void LandCarbonDensities::setActualBelowGroundCarbonDensity( const double aBelowGroundCarbonDensity,
                                                     const int aYear )
{
  if (getActualBelowGroundCarbonDensity( aYear ) > 0)
    mBelowGroundAdj[ aYear ] = aBelowGroundCarbonDensity / getActualBelowGroundCarbonDensity( aYear );
}

void LandCarbonDensities::setMatureAge( const int aMatureAge )    
{
    //assert( mMatureAge > 0 );
    mMatureAge = aMatureAge;
    
    // Precompute the sigmoid curve differnce to avoid doing it during calc.
    // Note this is only necessary when the mature age is greater than 1.
    if( mMatureAge > 1 ) {
        precalc_sigmoid_diff.resize( CarbonModelUtils::getEndYear() - CarbonModelUtils::getStartYear() + 1 );
        double prevSigmoid = pow( 1 - exp( ( -3.0 * 0 ) / mMatureAge ), 2.0 );
        for ( int i = CarbonModelUtils::getStartYear(); i <= CarbonModelUtils::getEndYear(); ++i ){
            const int offestYear = i - CarbonModelUtils::getStartYear();
            double currSigmoid = pow( 1 - exp( ( -3.0 * ( offestYear + 1 ) ) / mMatureAge ), 2.0 );
            precalc_sigmoid_diff[ offestYear ] = currSigmoid - prevSigmoid;
            prevSigmoid = currSigmoid;
        }
    }
}

double LandCarbonDensities::getActualAboveGroundCarbonDensity( const int aYear ) const {
    return mAboveGroundCarbon[ aYear ] * mAboveGroundAdj[ aYear ];
}

double LandCarbonDensities::getActualBelowGroundCarbonDensity( const int aYear ) const {
    return mBelowGroundCarbon[ aYear ] * mBelowGroundAdj[ aYear ];
}

int LandCarbonDensities::getMatureAge( ) const {
	return mMatureAge;
}
