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
* \file factor_supply.cpp
* \ingroup Objects
* \brief The FactorSupply class source file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "sectors/include/factor_supply.h"
#include "sectors/include/more_sector_info.h"
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "util/base/include/model_time.h"
#include "marketplace/include/marketplace.h"
#include "marketplace/include/imarket_type.h"
#include "util/base/include/ivisitor.h"
#include "util/logger/include/ilogger.h"
#include "containers/include/iinfo.h"
#include "functions/include/function_utils.h"
#include "util/base/include/configuration.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;

//! Default Constructor
FactorSupply::FactorSupply() {
    // resize vectors
    mBasePrice = 0;
    mBaseSupply = 0;
    // instantiate moreSectorInfo object here,
    // override with data read in
    moreSectorInfo.reset( new MoreSectorInfo() );
}

//! Get the name of the factor supply
const string& FactorSupply::getName() const {
    return name;
}
    
//! Parses SOME xml object
bool FactorSupply::XMLParse( const xercesc::DOMNode* node) {
    /*! \pre make sure we were passed a valid node. */
    assert( node );

    // get the name attribute.
    name = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "supply" ) {
            mBaseSupply = XMLHelper<double>::getValue( curr );
        }
        else if (nodeName == "price" ) {
            mBasePrice = XMLHelper<double>::getValue( curr );
        }
        else if( nodeName == MoreSectorInfo::getXMLNameStatic() ) {
            moreSectorInfo->XMLParse( curr );
        }
        else {
            return false;
            // Should print an error here.
        }
    }
    return true;
}

//! Write to XML
void FactorSupply::toInputXML( ostream &out, Tabs* tabs ) const {

    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLName(), out, tabs, name );
    
    // Write the data.
    XMLWriteElement( mBaseSupply, "supply", out, tabs );
    XMLWriteElement( mBasePrice, "price", out, tabs );

    // write the closing tag.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Debug info written to XML
void FactorSupply::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    // write the beginning tag.
    XMLWriteOpeningTag ( getXMLName(), out, tabs, name );

    XMLWriteElement( mBaseSupply, "supply", out, tabs );
    XMLWriteElement( mBasePrice, "price", out, tabs );

    // write the closing tag.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

double FactorSupply::getSupply( const string& aRegionName, const int period ) const {
    Marketplace* marketplace = scenario->getMarketplace();
    return marketplace->getSupply( name, aRegionName, period );
}

const string& FactorSupply::getXMLName() const {
    return getXMLNameStatic();
}

const string& FactorSupply::getXMLNameStatic() {
    const static string XML_NAME = "factorSupply";
    return XML_NAME;
}

void FactorSupply::completeInit( const string& aRegionName ) {
    setMarket( aRegionName );
}
/*! \brief Perform any initializations needed for each period.
*
* Any initializations or calculations that only need to be done once per period (instead of every iteration) should be placed in this function.
*
* \author Sonny Kim
* \param aRegionName region name
* \param aPeriod Model period
*/
void FactorSupply::initCalc( const string& aRegionName, const int aPeriod ) {
    Marketplace* marketplace = scenario->getMarketplace();
    if( aPeriod == 0 ) {
        // if the user does not read in a base price calculate it from the
        // base physical quantity read in and the total demand which will be in
        // values
        // note that the price for Capital (the intrest rate) can not be backed out
        // here, it must be-pre solved for from the input dataset and read in
        if( mBasePrice == 0 && name != "Capital" ) {
            assert( mBaseSupply != 0 );
            double baseDemand = marketplace->getDemand( name, aRegionName, aPeriod );
            mBasePrice = marketplace->getDemand( name, aRegionName, aPeriod ) / mBaseSupply;
        }
        // if we still do not have a price set it to 1
        if( mBasePrice == 0 ){
            mBasePrice = 1;
        }
        marketplace->setPrice( name, aRegionName, mBasePrice, aPeriod );
    }
    // set market prices before calling calcPricePaid
    calcPricePaid( aRegionName, aPeriod );
}

void FactorSupply::setMarket( const string& aRegionName ) {
    Marketplace* marketplace = scenario->getMarketplace();
    // Change modeltime to get start period
    const int START_PERIOD = scenario->getModeltime()->getBasePeriod();

    // Set the market name to the region name if it is not read in.
    // Currently this is never read in, not sure if this should be fixed.
    if( marketName.empty() ){
        marketName = aRegionName;
    }

    // name is factor supply name (name of primary factors)
    // market is the name of the regional market from the input file (i.e., global, region, regional group, etc.)
    if( marketplace->createMarket( aRegionName, marketName, name, IMarketType::NORMAL ) ) {
        if( !Configuration::getInstance()->getBool( "CalibrationActive" ) ){
            for( int per = START_PERIOD; per < scenario->getModeltime()->getmaxper(); ++per ){
                marketplace->setMarketToSolve( name, aRegionName, per );
            }
        }
    }
    // TODO: maybe a better place to put this
    if( name == "Capital" ) {
        // initialize the price of the capital good as 1
        FunctionUtils::setCapitalGoodPrice( aRegionName, 0, 1 );
    }
}

/*! \brief For calculating the price paid for primary factors of input
 * \author Sonny Kim
 * Price paid for factor inputs should be calculated in production and demand
 * sectors.  These price paid values are place holders for the final demand
 * sectors as they do not have demands for factor inputs.
 * \param period 
*/
void FactorSupply::calcPricePaid( const string& aRegionName, const int period ) {
    Marketplace* marketplace = scenario->getMarketplace();
    // set price received in market info
    double pricePaid = ( marketplace->getPrice(name, aRegionName, period) + 
        ( moreSectorInfo->getValue(MoreSectorInfo::TRANSPORTATION_COST)
        * moreSectorInfo->getValue(MoreSectorInfo::TRAN_COST_MULT) )
        * moreSectorInfo->getValue(MoreSectorInfo::PROPORTIONAL_TAX_RATE)
        + moreSectorInfo->getValue(MoreSectorInfo::ADDITIVE_TAX) ) // add carbon taxes
        * 1 ;
        //* moreSectorInfo->getValue(MoreSectorInfo::PRICE_ADJUST_MULT);
        
    // set price paid in market info.
    FunctionUtils::setPricePaid( aRegionName, name, period, pricePaid );
}

/*! \brief For outputing SGM data to a flat csv File
 * 
 * \author Pralit Patel, Sonny Kim
 * \param period The period which we are outputing for
 */
void FactorSupply::csvSGMOutputFile( ostream& aFile, const int period ) const {
    // Write factor supply output
    // aFile << "Factor Supply Results" << endl << endl;
}

void FactorSupply::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitFactorSupply( this, aPeriod );
    aVisitor->endVisitFactorSupply( this, aPeriod );
}
