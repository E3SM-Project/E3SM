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
* \file policy_ghg.cpp
* \ingroup Objects
* \brief GHGPolicy class source file.
* \author Sonny Kim
* \date $Date: 2007/01/11 00:12:30 $
* \version $Revision: 1.12.2.4 $
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <iostream>
#include <string>

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include "util/base/include/xml_helper.h"
#include "containers/include/scenario.h"
#include "containers/include/iinfo.h"
#include "util/base/include/model_time.h"
#include "policy/include/policy_ghg.h"
#include "marketplace/include/marketplace.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
// static initialize.
const string GHGPolicy::XML_NAME = "ghgpolicy";

/*! \brief Default constructor. */
GHGPolicy::GHGPolicy():
    isFixedTax( false ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1 ),
    mFixedTax( scenario->getModeltime()->getmaxper(), -1 ),
    mProportionalTaxRate( scenario->getModeltime()->getmaxper(), -1.0 )
{
}

/*!
* \brief Constructor which initializes a GHG policy without setting a tax or
*        constraint.
*/
GHGPolicy::GHGPolicy( const string aName, const string aMarket ):
    mName( aName ), 
    mMarket( aMarket ),
    isFixedTax( false ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1 ),
    mFixedTax( scenario->getModeltime()->getmaxper(), -1 ),
    mProportionalTaxRate( scenario->getModeltime()->getmaxper(), -1.0 )
{
}

/*! \brief Constructor used when explicitly constructing a fixed tax.
*/
GHGPolicy::GHGPolicy( const string aName, const string aMarket,
                      const vector<double>& aTaxes ):
    mName( aName ), 
    mMarket( aMarket ),
    mFixedTax( aTaxes ),
    isFixedTax( true ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1 ),
    mProportionalTaxRate( scenario->getModeltime()->getmaxper(), -1.0 ){
    // Ensure that the taxes vector passed in is the right size.
    assert( aTaxes.size() == mConstraint.size() );
}

/*! \brief Create a copy of the GHG policy.
* \return An exact copy of the policy.
*/
GHGPolicy* GHGPolicy::clone() const {
    return new GHGPolicy( *this );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Josh Lurz, James Blackwood
* \return The constant XML_NAME.
*/
const string& GHGPolicy::getXMLName() const {
    return XML_NAME;
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
const string& GHGPolicy::getXMLNameStatic() {
    return XML_NAME;
}

//! Get the ghg policy name. 
const string& GHGPolicy::getName() const {
    return mName;
}

//! Initializes data members from XML.
void GHGPolicy::XMLParse( const DOMNode* node ){

    /*! \pre assume we are passed a valid node.*/
    assert( node );

    // get the name attribute.
    mName = XMLHelper<string>::getAttr( node, "name" );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();
    const Modeltime* modeltime = scenario->getModeltime();
    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if( nodeName == "market" ){
            mMarket = XMLHelper<string>::getValue( curr ); // should be only one market
        }
        else if( nodeName == "isFixedTax" ) {
            isFixedTax = XMLHelper<bool>::getValue( curr );
        }
        else if( nodeName == "constraint" ){
            XMLHelper<double>::insertValueIntoVector( curr, mConstraint, modeltime );
        }
        else if( nodeName == "fixedTax" ){
            XMLHelper<double>::insertValueIntoVector( curr, mFixedTax, modeltime );
        }
        else if( nodeName == "proportional-tax-rate" ){
            XMLHelper<double>::insertValueIntoVector( curr, mProportionalTaxRate, modeltime );
            // Check to see if proportional tax rate is within valid range (between 0 and 1).
            // This could be used in the future for subsidies (negative rates) or multipliers
            // (greater than 1 ) on taxes.
            double tempProportionalTaxRate = XMLHelper<double>::getValue( curr );
            if ( tempProportionalTaxRate < 0 || tempProportionalTaxRate > 1 ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Proportional tax rate on ghg policy out of range (rate < 0 or rate > 1)." << endl;
            }
        }
        else {
            cout << "Unrecognized text string: " << nodeName << " found while parsing ghgmarket." << endl;
        }
    }
}

//! Writes data members to data stream in XML format.
void GHGPolicy::toInputXML( ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );
    XMLWriteElement( mMarket, "market", out, tabs );
    XMLWriteElement( isFixedTax, "isFixedTax", out, tabs );
    
    const Modeltime* modeltime = scenario->getModeltime();    
    XMLWriteVector( mConstraint, "constraint", out, tabs, modeltime, -1.0 );
    XMLWriteVector( mFixedTax, "fixedTax", out, tabs, modeltime, 0.0 );
    for( int per = 0; per < modeltime->getmaxper(); ++per ){
        XMLWriteElementCheckDefault( mProportionalTaxRate[ per ],
            "proportional-tax-rate", out, tabs, -1.0 );
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Writes data members to data stream in XML format.
void GHGPolicy::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );

    // write out the market string.
    XMLWriteElement( mMarket, "market", out, tabs );

    // Write whether we are a fixed tax policy.
    XMLWriteElement( isFixedTax, "isFixedTax", out, tabs );

    // Write the mConstraint for the current year
    XMLWriteElement( mConstraint[ period ], "constraint", out, tabs );
    
    // Write out the fixed tax for the current year.
    XMLWriteElement( mFixedTax[ period ], "fixedTax", out, tabs );
    
    // Write out the proportional tax rate for the current year.
    XMLWriteElement( mProportionalTaxRate[ period ], "proportional-tax-rate", out, tabs );

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Complete the initialization of the GHG policy.
* \details This function initializes a ghg market for the policy.
* GHG markets are created for both mConstraint and fixed tax policies.
* In the fixed tax policy, market prices are set to the fixed taxes, but
* the markets are not solved.  Also for the fixed tax policy, if the market name
* is the same for all regions, the fixed tax vector of the last region overrides
* the market prices.
* \author Sonny Kim and Josh Lurz
* \param regionName The name of the region the policy controls. 
*/
void GHGPolicy::completeInit( const string& aRegionName ) {
    const Modeltime* modeltime = scenario->getModeltime();
    Marketplace* marketplace = scenario->getMarketplace();
    marketplace->createMarket( aRegionName, mMarket, mName, IMarketType::TAX );

    // Set price and output units for period 0 market info
    IInfo* marketInfo = marketplace->getMarketInfo( mName, aRegionName, 0, true );
    //TODO: read-in as data the units of tax and emissions
    marketInfo->setString( "price-unit", "1990$/tC" );
    marketInfo->setString( "output-unit", "MTC" );

    // Set the proportional tax rate into the policy market info object for
    // retrieval by the emissions object.
    for( int per = 0; per < modeltime->getmaxper(); ++per ){
        IInfo* currMarketInfo = marketplace->getMarketInfo( mName, aRegionName, per, true );
        // Note: the key includes the region name.
        if( mProportionalTaxRate[ per ] != -1 ){
            currMarketInfo->setDouble( "proportional-tax-rate" + aRegionName, mProportionalTaxRate[ per ] );
        }
    }
    // check for missing periods in which case interpolate
    for( int i = 1; i < modeltime->getmaxper(); ++i ) {
        if( mFixedTax[ i ] == -1 && mFixedTax[ i - 1 ] != -1 ) {
            int j;
            for( j = i + 1; j < modeltime->getmaxper() && mFixedTax[ j ] == -1; ++j ) {
            }
            if( j < modeltime->getmaxper() ) {
                mFixedTax[ i ] = util::linearInterpolateY( modeltime->getper_to_yr( i ),
                                                           modeltime->getper_to_yr( i - 1 ),
                                                           modeltime->getper_to_yr( j ),
                                                           mFixedTax[ i - 1 ],
                                                           mFixedTax[ j ] );
            }
        }
        if( mConstraint[ i ] == -1 && mConstraint[ i - 1 ] != -1 ) {
            int j;
            for( j = i + 1; j < modeltime->getmaxper() && mConstraint[ j ] == -1; ++j ) {
            }
            if( j < modeltime->getmaxper() ) {
                mConstraint[ i ] = util::linearInterpolateY( modeltime->getper_to_yr( i ),
                                                             modeltime->getper_to_yr( i - 1 ),
                                                             modeltime->getper_to_yr( j ),
                                                             mConstraint[ i - 1 ],
                                                             mConstraint[ j ] );
            }
        }
    }

    // Loop through each period
    // If it is a fixed tax, set the tax level and set the market not to solve
    // If it is a constraint, add the constraint to the market and set the 
    // market to solve.
    for( unsigned int i = 0; i < modeltime->getmaxper(); ++i ){
        if( mFixedTax[ i ] != -1 ){
            marketplace->unsetMarketToSolve( mName, aRegionName, i );
            marketplace->setPrice( mName, aRegionName, mFixedTax[ i ], i );
        }
        else if( mConstraint[ i ] != -1 ){
            marketplace->setMarketToSolve( mName, aRegionName, i );
            // Adding the difference between the constraint for this period
            // and the current supply because addToSupply adds to the current
            // supply.  Passing false to suppress a warning the first time through.
            marketplace->addToSupply( mName, aRegionName, mConstraint[ i ] - 
                marketplace->getSupply( mName, aRegionName, i ), i, false );
        }
    }
}

/*! \brief Determine if the tax is applicable for a given region.
* \param aRegion Region name.
* \return Whether the tax is applicable.
* \todo This is not entirely correct for multiple regions within a market.
*/
bool GHGPolicy::isApplicable( const string& aRegion ) const {
    return mMarket == "global" || mMarket == aRegion;
}

/*!
* \brief Set the mConstraint to the vector passed in.
* \param aConstraint new mConstraint vector
*/
void GHGPolicy::setConstraint( const vector<double>& aConstraint ){
    isFixedTax = false;
    mConstraint = aConstraint;
}
