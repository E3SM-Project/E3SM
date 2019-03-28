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
 * \file depleting_fixed_resource.cpp
 * \ingroup Objects
 * \brief DepletingFixedResource class source file.
 * \author Pralit Patel
 */

#include "util/base/include/definitions.h"

#include <string>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "resources/include/depleting_fixed_resource.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "marketplace/include/imarket_type.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Get the XML name of the class.
 * \return The XML name of the class.
 */
const string& DepletingFixedResource::getXMLNameStatic(){
    static const string XML_NAME = "depleting-fixed-resource";
    return XML_NAME;
}

//! Constructor.
DepletingFixedResource::DepletingFixedResource()
: mDepletionRate( 0 ), mInitialPrice( 1 ) {
}

//! Destructor.
DepletingFixedResource::~DepletingFixedResource() {
}

void DepletingFixedResource::XMLParse( const DOMNode* node ){

    // make sure we were passed a valid node.
    assert( node );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    const DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        const DOMNode* curr = nodeList->item( i );
        const string nodeName =
            XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "output-unit" ){
            mOutputUnit = XMLHelper<string>::getValue( curr );
        }
        else if( nodeName == "price-unit" ){
            mPriceUnit = XMLHelper<string>::getValue( curr );
        }
        else if( nodeName == "market" ){
            mMarket = XMLHelper<string>::getValue( curr );
        }
        else if( nodeName == "fixed-resource" ){
            mFixedResource = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "depletion-rate" ){
            mDepletionRate = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "price" ){
            mInitialPrice = XMLHelper<double>::getValue( curr );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName
                    << " found while parsing " << getXMLNameStatic() << "."
                    << endl;
        }
    }
}

const string& DepletingFixedResource::getXMLName() const {
    return getXMLNameStatic();
}

void DepletingFixedResource::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );

    // write the xml for the class members.
    XMLWriteElement( mOutputUnit, "output-unit", aOut, aTabs );
    XMLWriteElement( mPriceUnit, "price-unit", aOut, aTabs );
    XMLWriteElement( mMarket, "market", aOut, aTabs );
    XMLWriteElement( mFixedResource, "fixed-resource", aOut, aTabs );
    XMLWriteElement( mDepletionRate, "depletion-rate", aOut, aTabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void DepletingFixedResource::toDebugXML( const int aPeriod,
                                    ostream& aOut,
                                    Tabs* aTabs ) const
{
    // this object does not currently maintain state for past periods
    toInputXML( aOut, aTabs );
}

void DepletingFixedResource::completeInit( const string& aRegionName,
                                      const IInfo* aRegionInfo )
{
    // default unit to EJ
    if ( mOutputUnit.empty() ) {
        mOutputUnit = "EJ"; 
    }
    // default unit to $/GJ
    if ( mPriceUnit.empty() ) {
        mPriceUnit = "1975$/GJ"; 
    }
    // Setup markets for this resource.
    setMarket( aRegionName );
}

void DepletingFixedResource::initCalc( const string& aRegionName,
                                  const int aPeriod )
{
    Marketplace* marketplace = scenario->getMarketplace();
    IInfo* marketInfo = marketplace->getMarketInfo( mName,
                                                          aRegionName,
                                                          aPeriod,
                                                          true );
    assert( marketInfo );
    marketInfo->setDouble( "depleted-resource", 0 );
    if( aPeriod > 0 ) {
        // get the depletion from the previous period and subtract it
        // from the current quantity
        marketInfo = marketplace->getMarketInfo( mName,
                                                 aRegionName,
                                                 aPeriod - 1,
                                                 true );

        double depletedResourceAnnual = marketInfo->getDouble( "depleted-resource", true );

        double depletedResourcePrevAnnual;
        if( aPeriod > 1 ) {
            marketInfo = marketplace->getMarketInfo( mName,
                                                     aRegionName,
                                                     aPeriod - 2,
                                                     true );
            depletedResourcePrevAnnual = marketInfo->getDouble( "depleted-resource", true );
        }
        else {
            // TODO: read something maybe?
            depletedResourcePrevAnnual = depletedResourceAnnual * 0.95;
        }

        // Cumulative production for the current period
        // is equal to the cumulative production of the previous period plus the
        // trapezoidal area formed by the previous annual production and the
        // current annual production.
        // Cumulative(t) = Cumulative(t-1) + 1/2 * (Annual(t) - Annual(t-1))* timestep.
        mFixedResource -= ( .5 * ( depletedResourceAnnual  - depletedResourcePrevAnnual ) + depletedResourcePrevAnnual )
            * scenario->getModeltime()->gettimestep( aPeriod - 1 );

        // check for overuse of the resource which may happen due to no formal
        // constraint in the extraction of this resource
        if( mFixedResource < 0 ) {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << aRegionName << " " << mName << " was depleted beyond zero in period "
                << aPeriod - 1 << " by " << mFixedResource << " " << mOutputUnit << "." << endl;

            mFixedResource = 0;
        }
    }
}

void DepletingFixedResource::postCalc( const string& aRegionName,
                                  const int aPeriod )
{
}


const string& DepletingFixedResource::getName() const {
    return mName;
}

void DepletingFixedResource::calcSupply( const string& aRegionName,
                                    const GDP* aGDP,
                                    const int aPeriod )
{
    Marketplace* marketplace = scenario->getMarketplace();

    // the supply is just the fixed amount that is left.
    marketplace->addToSupply( mName, aRegionName, mFixedResource, aPeriod );
}

double DepletingFixedResource::getAnnualProd( const string& aRegionName,
                                         const int aPeriod ) const
{
    // Return the market supply.
    return scenario->getMarketplace()->getSupply( mName, aRegionName, aPeriod );
}

void DepletingFixedResource::csvOutputFile( const string& aRegionName )
{
    // function protocol
    void fileoutput3( string var1name,string var2name,string var3name,
        string var4name,string var5name,string uname,vector<double> dout);

    const int maxper = scenario->getModeltime()->getmaxper();
    vector<double> temp( maxper );
    for( int i = 0; i < maxper; ++i ){
        temp[ i ] = getAnnualProd( aRegionName, i );
    }
    fileoutput3( aRegionName , mName," "," ","production", mOutputUnit, temp );
}

void DepletingFixedResource::dbOutput( const string& aRegionName ){
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    vector<double> temp(maxper);
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // Subsectors do not exist for depleting fixed Resource.
    for (int m=0;m<maxper;m++) {
        temp[m] += getAnnualProd(aRegionName, m);
    }
    dboutput4( aRegionName, "Resource", "annual-production", mName, mOutputUnit, temp );
}

void DepletingFixedResource::setCalibratedSupplyInfo( const int aPeriod,
                                                 const string& aRegionName )
{
    // TODO: is this fine?
    const double MKT_NOT_ALL_FIXED = -1;
    Marketplace* marketplace = scenario->getMarketplace();

    IInfo* resourceInfo = marketplace->getMarketInfo( mName, aRegionName,
                                                      aPeriod, true );
    assert( resourceInfo );

    resourceInfo->setDouble( "calSupply", MKT_NOT_ALL_FIXED );
}

/*
* \brief Create the resource market.
* \details The resource creates a single solved market for the
*          resource.
* \param aRegionName Region name.
*/
void DepletingFixedResource::setMarket( const string& aRegionName ) {
    // Setup the market for the resource. This market will not be solved. Note
    // that in a standard Resource setMarketToSolve would be called here.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->createMarket( aRegionName, mMarket, mName,
                               IMarketType::NORMAL );

    // Set price and output units for period 0 market info
    IInfo* marketInfo = marketplace->getMarketInfo( mName, aRegionName, 0, true );
    marketInfo->setString( "price-unit", mPriceUnit );
    marketInfo->setString( "output-unit", mOutputUnit );
    // Set the depletion rate for each period.
    const Modeltime* modeltime = scenario->getModeltime();
    bool calibrationActive = Configuration::getInstance()->getBool( "CalibrationActive" );
    for( int i = 0; i < modeltime->getmaxper(); ++i ){
        if( !calibrationActive ) {
            marketplace->setMarketToSolve( mName, aRegionName, i );
        }
        marketplace->setPrice( mName, aRegionName, mInitialPrice, i );

        marketInfo = marketplace->getMarketInfo( mName, aRegionName, i, true );
        marketInfo->setDouble( "depletion-rate", mDepletionRate );
    }
}

void DepletingFixedResource::accept( IVisitor* aVisitor,
                                const int aPeriod ) const
{
    aVisitor->startVisitResource( this, aPeriod );
    aVisitor->endVisitResource( this, aPeriod );
}
