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
* \file total_sector_emissions.cpp
* \ingroup Objects
* \brief TotalSectorEmissions class source file.
* \author James Blackwood
*/

#include "util/base/include/definitions.h"
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/iinfo.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "emissions/include/total_sector_emissions.h"
#include "sectors/include/cal_quantity_tabulator.h"
#include "sectors/include/sector.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*! \brief Constructor
*/
TotalSectorEmissions::TotalSectorEmissions()
// TODO: Storing the none string is a waste.
: mType( "none" ),
mApplicableYear( 0 ){
}

void TotalSectorEmissions::XMLParse( const DOMNode* aNode ){
	/*! \pre make sure we were passed a valid node. */
    assert( aNode );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( aNode, "name" );

    // get all child nodes.
    const DOMNodeList* nodeList = aNode->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
		else if( nodeName == "type" ){
            mType = XMLHelper<string>::getValue( curr );
        }
		else if( nodeName == "value" ){
            mAggregateEmissions = XMLHelper<double>::getValue( curr );
        }
		else if( nodeName == "year" ){
            mApplicableYear = XMLHelper<int>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing "
                    << getXMLNameStatic() << "." << endl;
        }
    }
}

void TotalSectorEmissions::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElementCheckDefault( mType, "type", aOut, aTabs, string( "none" ) );

    if( mAggregateEmissions.isInited() ){
        XMLWriteElement( mAggregateEmissions, "value", aOut, aTabs );
    }

    XMLWriteElementCheckDefault( mApplicableYear, "year", aOut, aTabs );
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void TotalSectorEmissions::toDebugXML( const int aPeriod, ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );
    XMLWriteElement( mType, "type", aOut, aTabs );
    XMLWriteElement( mAggregateEmissions, "value", aOut, aTabs );
    XMLWriteElement( mApplicableYear, "year", aOut, aTabs );
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
const string& TotalSectorEmissions::getXMLNameStatic() {
    const static string XML_NAME = "TotalSectorEmissions";
	return XML_NAME;
}

/*! \brief Returns the name.
* \author James Blackwood
* \return The name.
*/
const string& TotalSectorEmissions::getName() const {
    return mName;
}

/*! \brief Returns prefix to use for storing aggregate emissions factor in info
*          object.
* \author Steve Smith
* \return The prefix string
*/
const string& TotalSectorEmissions::aggrEmissionsPrefix() {
    const static string AGGR_EM_PREFIX = "AggrEmissionFactor"; 
    return AGGR_EM_PREFIX;
}

/*! \brief Returns prefix to use for storing aggregate emissions factor 
*          applicable year in info object.
* \author Steve Smith
* \return The prefix string
*/
const string& TotalSectorEmissions::aggrEmissionsYearPrefix() {
    const static string AGGR_EM_PREFIX_YEAR = "AggrEmissionYear"; 
    return AGGR_EM_PREFIX_YEAR;
}

/*! \brief Sets aggregate emissions factor into the region info object.
* \details Calculates the total calibrated output for sectors of the specified
*          type in the specified year from the list of sectors. This is used to
*          calculate an emissions coefficient for the gas given the total
*          emissions for the set of sectors of the specified type.
* \author James Blackwood
* \param aRegionName Name of the containing region.
* \param aSectors Vector of supply sectors used to determine the total
*        calibrated output.
* \param aRegionInfo Region info object which will be updated to contain the
*        emissions coefficients.
* \param aPeriod Model period.
*/
void TotalSectorEmissions::setAggregateEmissionFactor( const string& aRegionName,
                                                       const std::vector<Sector*>& aSectors,
                                                       IInfo* aRegionInfo,
                                                       const int aPeriod ) const
{
    // Write applicable year to info object
    aRegionInfo->setDouble( aggrEmissionsYearPrefix() + mName, mApplicableYear );
    
    // Check if the object is applicable in the current period.
    const Modeltime* modeltime = scenario->getModeltime();
    if ( mApplicableYear <= 0 || modeltime->getyr_to_per( mApplicableYear ) != aPeriod ) {
        return;
    }

    if( !mAggregateEmissions.isInited() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "No aggregate emissions were read-in for " << mName << "." << endl;
        return;
    }

    typedef vector<Sector*>::const_iterator SectorIterator;
    // Visit all appropriate sectors and sum their emissions using a cal
    // quantity tabulator.
    CalQuantityTabulator tabulator( aRegionName );
    tabulator.setApplicableSectorType( mType );
    for( SectorIterator currSector = aSectors.begin(); currSector != aSectors.end(); ++currSector ){
        (*currSector)->accept( &tabulator, aPeriod );
    }

    CalQuantityTabulator::CalInfoMap calInfo = tabulator.getSupplyInfo();
    // Loop through the calibrated info map and sum output and fixed output
    // for all sectors. Check for any uncalibrated sectors.
    bool isAllCalibrated = true;
    double totalSectorOutput = 0;
    for( CalQuantityTabulator::CalInfoMap::const_iterator i = calInfo.begin(); i != calInfo.end(); ++i ){
        CalQuantityTabulator::CalInfo currInfo = i->second;
        isAllCalibrated &= currInfo.mAllFixed;
        totalSectorOutput += currInfo.mFixedQuantity + currInfo.mCalQuantity;
    }

    double emissionsFactor;
    if( !isAllCalibrated || totalSectorOutput < util::getSmallNumber() ){
        ILogger& mainLog = ILogger::getLogger( "main_log" );
        mainLog.setLevel( ILogger::WARNING );
        mainLog << "Cannot set aggregate emissions factor for GHG " << mName
                << " because sectors are not fully calibrated or have zero output." << endl;
        // Set an emissions factor of zero after printing the warning for each
        // use of the emissions factor.
        emissionsFactor = 0;
    }
    else {
        emissionsFactor = mAggregateEmissions / totalSectorOutput;
    }

    // Store aggregate emissions factor in region info object
    aRegionInfo->setDouble( aggrEmissionsPrefix() + mName, emissionsFactor );
}
