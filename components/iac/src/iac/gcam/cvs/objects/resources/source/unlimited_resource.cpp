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
 * \file unlimited_resource.cpp
 * \ingroup Objects
 * \brief UnlimitedResource class source file.
 * \author Josh Lurz
 */

#include "util/base/include/definitions.h"

#include <string>
#include <vector>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "resources/include/unlimited_resource.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "marketplace/include/imarket_type.h"
#include "containers/include/iinfo.h"
#include "util/base/include/ivisitor.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

/*!
 * \brief Get the XML name of the class.
 * \return The XML name of the class.
 */
const string& UnlimitedResource::getXMLNameStatic(){
    static const string XML_NAME = "unlimited-resource";
    return XML_NAME;
}

//! Constructor.
UnlimitedResource::UnlimitedResource()
: mFixedPrices( scenario->getModeltime()->getmaxper() ){
}

//! Destructor.
UnlimitedResource::~UnlimitedResource() {
}

void UnlimitedResource::XMLParse( const DOMNode* node ){

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
        else if( nodeName == "price" ){
            XMLHelper<double>::insertValueIntoVector( curr, mFixedPrices,
                                                      scenario->getModeltime() );
        }
        else if( nodeName == "capacity-factor" ){
            mCapacityFactor = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == "variance" ){
            mVariance = XMLHelper<double>::getValue( curr );
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

const string& UnlimitedResource::getXMLName() const {
    return getXMLNameStatic();
}

void UnlimitedResource::toInputXML( ostream& aOut, Tabs* aTabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );

    // write the xml for the class members.
    XMLWriteElement( mOutputUnit, "output-unit", aOut, aTabs );
    XMLWriteElement( mPriceUnit, "price-unit", aOut, aTabs );
    XMLWriteElement( mMarket, "market", aOut, aTabs );
    
    XMLWriteElementCheckDefault( mCapacityFactor, "capacity-factor", aOut,
                                 aTabs, Value( 0 ) );

    XMLWriteElementCheckDefault( mVariance, "variance", aOut,
                                 aTabs, Value( 0 ) );
    
    const Modeltime* modeltime = scenario->getModeltime();
    XMLWriteVector( mFixedPrices, "price", aOut, aTabs, modeltime, 0.0 );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void UnlimitedResource::toDebugXML( const int aPeriod,
                                    ostream& aOut,
                                    Tabs* aTabs ) const
{
    XMLWriteOpeningTag( getXMLNameStatic(), aOut, aTabs, mName );

    // Write the xml for the class members.
    // Write out the market string.
    XMLWriteElement( mOutputUnit, "output-unit", aOut, aTabs );
    XMLWriteElement( mPriceUnit, "price-unit", aOut, aTabs );
    XMLWriteElement( mMarket, "market", aOut, aTabs );
    XMLWriteElement( mCapacityFactor, "capacity-factor", aOut, aTabs );
    XMLWriteElement( mVariance, "variance", aOut, aTabs );

    // Write out resource prices for debugging period.
    XMLWriteElement( mFixedPrices[ aPeriod ], "price", aOut, aTabs );

    // finished writing xml for the class members.

    XMLWriteClosingTag( getXMLNameStatic(), aOut, aTabs );
}

void UnlimitedResource::completeInit( const string& aRegionName,
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
    // If not read-in, set variance to 0 and isInited to true.
    // The set() call sets isInited to true.
    if( !mVariance.isInited() ){
        mVariance.set( 0 );
    }
    // Setup markets for this resource.
    setMarket( aRegionName );
}

void UnlimitedResource::initCalc( const string& aRegionName,
                                  const int aPeriod )
{
    // Set the capacity factor and variance.
    IInfo* marketInfo = scenario->getMarketplace()->getMarketInfo( mName,
                                                                   aRegionName,
                                                                   aPeriod,
                                                                   true );
    assert( marketInfo );

    if( mCapacityFactor.isInited() ){
        marketInfo->setDouble( "resourceCapacityFactor", mCapacityFactor );
    }
    if( mVariance.isInited() ){
        marketInfo->setDouble( "resourceVariance", mVariance );
    }
}

void UnlimitedResource::postCalc( const string& aRegionName,
                                  const int aPeriod )
{
    // Reset initial resource prices to solved prices
    mFixedPrices[ aPeriod ] = scenario->getMarketplace()->getPrice( mName, aRegionName,
                              aPeriod, true );
}


const string& UnlimitedResource::getName() const {
    return mName;
}

void UnlimitedResource::calcSupply( const string& aRegionName,
                                    const GDP* aGDP,
                                    const int aPeriod )
{
    Marketplace* marketplace = scenario->getMarketplace();

    // Get the current demand and add the difference between current supply and
    // demand to the market.
    double currDemand = marketplace->getDemand( mName, aRegionName, aPeriod );
    double currSupply = marketplace->getSupply( mName, aRegionName, aPeriod );
    assert( currDemand >= currSupply );
    marketplace->addToSupply( mName, aRegionName, currDemand - currSupply,
                              aPeriod );
}

double UnlimitedResource::getAnnualProd( const string& aRegionName,
                                         const int aPeriod ) const
{
    // Return the market supply.
    return scenario->getMarketplace()->getSupply( mName, aRegionName, aPeriod );
}

void UnlimitedResource::csvOutputFile( const string& aRegionName )
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

void UnlimitedResource::dbOutput( const string& aRegionName ){
    const Modeltime* modeltime = scenario->getModeltime();
    const int maxper = modeltime->getmaxper();
    vector<double> temp(maxper);
    // function protocol
    void dboutput4(string var1name,string var2name,string var3name,string var4name,
        string uname,vector<double> dout);

    // Subsectors do not exist for Unlimited Resource.
    for (int m=0;m<maxper;m++) {
        temp[m] += getAnnualProd(aRegionName, m);
    }
    dboutput4( aRegionName, "Resource", "annual-production", mName, mOutputUnit, temp );
}

void UnlimitedResource::setCalibratedSupplyInfo( const int aPeriod,
                                                 const string& aRegionName )
{
    const double MKT_NOT_ALL_FIXED = -1;
    Marketplace* marketplace = scenario->getMarketplace();

    IInfo* resourceInfo = marketplace->getMarketInfo( mName, aRegionName,
                                                      aPeriod, true );
    assert( resourceInfo );

    resourceInfo->setDouble( "calSupply", MKT_NOT_ALL_FIXED );
}

/*
* \brief Create the resource market.
* \details The unlimited resource creates a single unsolved market for the
*          resource. The object will ensure that supply is always equal to
*          demand.
* \param aRegionName Region name.
*/
void UnlimitedResource::setMarket( const string& aRegionName ) {
    // Setup the market for the resource. This market will not be solved. Note
    // that in a standard Resource setMarketToSolve would be called here.
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->createMarket( aRegionName, mMarket, mName, IMarketType::NORMAL );

    // Set price and output units for period 0 market info
    IInfo* marketInfo = marketplace->getMarketInfo( mName, aRegionName, 0, true );
    marketInfo->setString( "price-unit", mPriceUnit );
    marketInfo->setString( "output-unit", mOutputUnit );
    // Need to set resource variance here because initCalc of technology is called
    // before that of resource. shk 2/27/07
    marketInfo->setDouble( "resourceVariance", mVariance );
    // Set the read-in fixed prices for each period.
    const Modeltime* modeltime = scenario->getModeltime();
    for( int i = 0; i < modeltime->getmaxper(); ++i ){
        if( mFixedPrices[ i ] != 0 ){
            marketplace->setPrice( mName, aRegionName, mFixedPrices[ i ], i, true );
        }
    }
}

void UnlimitedResource::accept( IVisitor* aVisitor,
                                const int aPeriod ) const
{
    aVisitor->startVisitResource( this, aPeriod );
    aVisitor->endVisitResource( this, aPeriod );
}
