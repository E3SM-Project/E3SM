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
 * \file so2_emissions.cpp
 * \ingroup objects
 * \brief SO2 class source file.
 * \author Nick Fernandez
 */
#include "util/base/include/definitions.h"
#include <string>
#include <iostream>
#include <xercesc/dom/DOMNode.hpp>
#include "emissions/include/so2_emissions.h"
#include "util/base/include/xml_helper.h"
#include "util/logger/include/ilogger.h"
#include "emissions/include/read_emissions_coef.h"

using namespace std;
using namespace xercesc;

//! Default Constructor.
SO2Emissions::SO2Emissions():
AComplexEmissions(),
ashRetention( 0 ),
percentSulfur( 0 ),
gjPerTonne( 1000 ),
finalSulfur( -1 )
{
}

//! Clone operator.
SO2Emissions* SO2Emissions::clone() const {
    return new SO2Emissions( *this );
}

/*!
 * \brief Get the XML node name for output to XML.
 * \details This public function accesses the private constant string, XML_NAME.
 *          This way the tag is always consistent for both read-in and output and can be easily changed.
 *          This function may be virtual to be overridden by derived class pointers.
 * \author Josh Lurz, James Blackwood
 * \return The constant XML_NAME.
 */
const string& SO2Emissions::getXMLName() const {
    return getXMLNameStatic();
}

const string& SO2Emissions::getXMLNameStatic() {
    static const string XML_NAME = "SO2";
    return XML_NAME;
}
const string& SO2Emissions::getName() const {
    return getXMLNameStatic();
}

void SO2Emissions::initCalc( const string& aRegionName,
                             const string& aFuelName,
                             const IInfo* aLocalInfo,
                             const int aPeriod )
{
    // TODO: Determine why this doesn't call the base class method.
    mEmissionsCoef.reset( new ReadEmissionsCoef( 1 ) );
}

/*! \brief Parses any child nodes specific to derived classes.
* \details Method parses any input data from child nodes that are specific to the classes derived from this class.
* \author Nick Fernandez
* \param nodeName name of current node
* \param curr pointer to the current node in the XML input tree
*/
bool SO2Emissions::XMLDerivedClassParse( const string& nodeName, const DOMNode* curr ){
    if( AComplexEmissions::XMLDerivedClassParse( nodeName, curr ) ){
    }
    else if( nodeName == "ashRetention" ){
        ashRetention = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "percentSulfur" ){
        percentSulfur = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "gjTonne" ){
        gjPerTonne = XMLHelper<double>::getValue( curr );
    }
    else if( nodeName == "finalSulfur" ){
        finalSulfur = XMLHelper<double>::getValue( curr );
    }
    else {
        return false;
    }
    return true;
}

// Do nothing because the name is always SO2.
void SO2Emissions::parseName( const string& aNameAttr ){
}

void SO2Emissions::toInputXMLDerived( ostream& out, Tabs* tabs ) const {
    AComplexEmissions::toInputXMLDerived( out, tabs );
    XMLWriteElementCheckDefault( ashRetention, "ashRetention", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( percentSulfur, "percentSulfur", out, tabs, 0.0 );
    XMLWriteElementCheckDefault( gjPerTonne, "gjTonne", out, tabs, 1000.0 );
    XMLWriteElementCheckDefault( finalSulfur, "finalSulfur", out, tabs, -1.0 );
}

void SO2Emissions::toDebugXMLDerived( const int period, ostream& out, Tabs* tabs ) const {
    AComplexEmissions::toDebugXMLDerived( period, out, tabs );
    XMLWriteElement( ashRetention, "ashRetention", out, tabs );
    XMLWriteElement( percentSulfur, "percentSulfur", out, tabs );
    XMLWriteElement( gjPerTonne, "gjTonne", out, tabs );
    XMLWriteElement( finalSulfur, "finalSulfur", out, tabs );
}

/*!
 * \brief Returns the emissions driver for SO2
 * \detailed SO2 emissions are converted from EJ to Tg by multiplying by 1000 
 *           and dividing by the ratio of gj per metric ton. This input is then multiplied by the percent
 *           sulfur content to give the weight of the sulfur, and then multiplied by the amount of that sulfur
 *           that escapes (is not retained in the form of ash)
 * \author Nick Fernandez
 * \author Jim Naslund
 * \param inputIn energy input
 * \param outputIn energy output
 * \return The emissions driver
 */
double SO2Emissions::emissionsDriver( const double inputIn, const double outputIn ) const {
    return AGHG::emissionsDriver( inputIn, outputIn ) *
           ( 1000 / gjPerTonne ) * ( percentSulfur / 100 ) * ( 1 - ( ashRetention / 100 ) );
}

/*!
 * \brief sets the variable adjMaxCntrl using finalSulfur (percent).
 * \detailed This is specific to sulfur emissions.
 * \author Nick Fernandez
 */
void SO2Emissions::setAdjMaxCntrl(){
    if ( ( finalSulfur <= percentSulfur ) && ( finalSulfur >= 0 ) ){
        double newMaxCntrl = ( 1 - ( finalSulfur / percentSulfur ) ) * 100;
        adjMaxCntrl = newMaxCntrl/maxCntrl;
    }
    else if (finalSulfur != -1){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog <<"finalSulfur is not in Valid range, percentSulfur:"<<percentSulfur<<" finalSulfur: "<<finalSulfur<<endl;
    }
}

void SO2Emissions::adjustMaxCntrl(const double GDPcap){
    // Set max control for final SO2 content
    setAdjMaxCntrl();

    AComplexEmissions::adjustMaxCntrl( GDPcap );
}
