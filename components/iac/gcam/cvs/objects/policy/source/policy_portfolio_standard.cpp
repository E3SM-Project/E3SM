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
* \file policy_portfolio_standard.cpp
* \ingroup Objects
* \brief PolicyPortfolioStandard class source file.
* \author Sonny Kim
* \date $Date: 2007/11/21 00:12:30 $
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
#include "policy/include/policy_portfolio_standard.h"
#include "marketplace/include/marketplace.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

extern Scenario* scenario;
// static initialize.
const string PolicyPortfolioStandard::XML_NAME = "policy-portfolio-standard";

/*! \brief Default constructor. */
PolicyPortfolioStandard::PolicyPortfolioStandard():
    isFixedTax( false ),
    mIsShareBased( false ),
    mFixedTax( scenario->getModeltime()->getmaxper(), -1 ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1.0 ),
    mShareOfSectorOutput( scenario->getModeltime()->getmaxper(), -1.0 )
{
}

/*!
* \brief Constructor which initializes a portfolio standard policy without setting 
*  an output share or constraint.
*/
PolicyPortfolioStandard::PolicyPortfolioStandard( const string aName, const string aMarket ):
    mName( aName ), 
    mMarket( aMarket ),
    isFixedTax( false ),
    mIsShareBased( false ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1.0 ),
    mFixedTax( scenario->getModeltime()->getmaxper(), -1 ),
    mShareOfSectorOutput( scenario->getModeltime()->getmaxper(), -1.0 )
{
}

/*! \brief Constructor used when constructing an output share .
*/
PolicyPortfolioStandard::PolicyPortfolioStandard( const string aName, const string aMarket,
                      const vector<double>& aShareOfTotal ):
    mName( aName ), 
    mMarket( aMarket ),
    mShareOfSectorOutput( aShareOfTotal ),
    mIsShareBased( true ),
    isFixedTax( false ),
    mConstraint( scenario->getModeltime()->getmaxper(), -1.0 )
{
    // Ensure that the share vector passed in is the right size.
    assert( aShareOfTotal.size() == mConstraint.size() );
}

/*! \brief Create a copy of the portfolio standard policy.
* \return An exact copy of the policy.
*/
PolicyPortfolioStandard* PolicyPortfolioStandard::clone() const {
    return new PolicyPortfolioStandard( *this );
}

/*! \brief Get the XML node name for output to XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* This function may be virtual to be overridden by derived class pointers.
* \author Sonny Kim
* \return The constant XML_NAME.
*/
const string& PolicyPortfolioStandard::getXMLName() const {
    return XML_NAME;
}

/*! \brief Get the XML node name in static form for comparison when parsing XML.
*
* This public function accesses the private constant string, XML_NAME.
* This way the tag is always consistent for both read-in and output and can be easily changed.
* The "==" operator that is used when parsing, required this second function to return static.
* \note A function cannot be static and virtual.
* \author Sonny Kim
* \return The constant XML_NAME as a static.
*/
const string& PolicyPortfolioStandard::getXMLNameStatic() {
    return XML_NAME;
}

//! Get the ghg policy name. 
const string& PolicyPortfolioStandard::getName() const {
    return mName;
}

//! Initializes data members from XML.
void PolicyPortfolioStandard::XMLParse( const DOMNode* node ){

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
        else if( nodeName == "policyType" ){
            mPolicyType = XMLHelper<string>::getValue( curr ); // should be only one market
        }
        else if( nodeName == "isShareBased" ) {
            mIsShareBased = XMLHelper<bool>::getValue( curr );
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
        else if( nodeName == "share-of-sector-output" ){
            XMLHelper<double>::insertValueIntoVector( curr, mShareOfSectorOutput, modeltime );
            // Check to see if the output share is within valid range (between 0 and 1).
            double tempShare = XMLHelper<double>::getValue( curr );
            if ( tempShare <= 0 || tempShare > 1 ){
                ILogger& mainLog = ILogger::getLogger( "main_log" );
                mainLog.setLevel( ILogger::ERROR );
                mainLog << "Output share for portfolio standard policy out of range "
                    <<"(rate <= 0 or rate > 1)." << endl;
            }
        }
        else {
            cout << "Unrecognized text string: " << nodeName 
                 << " found while parsing portfolio standard policy." << endl;
        }
    }
}

//! Writes data members to data stream in XML format.
void PolicyPortfolioStandard::toInputXML( ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );
    XMLWriteElement( mMarket, "market", out, tabs );
    XMLWriteElement( mPolicyType, "policyType", out, tabs );
    XMLWriteElement( mIsShareBased, "isShareBased", out, tabs );
    XMLWriteElement( isFixedTax, "isFixedTax", out, tabs );
    
    const Modeltime* modeltime = scenario->getModeltime();    
    for( int per = 0; per < modeltime->getmaxper(); ++per ){
        XMLWriteElementCheckDefault( mConstraint[ per ],
            "constraint", out, tabs, -1.0 );
    }
    XMLWriteVector( mFixedTax, "fixedTax", out, tabs, modeltime, 0.0 );
    for( int per = 0; per < modeltime->getmaxper(); ++per ){
        XMLWriteElementCheckDefault( mShareOfSectorOutput[ per ],
            "share-of-sector-output", out, tabs, -1.0 );
    }

    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

//! Writes data members to data stream in XML format.
void PolicyPortfolioStandard::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {

    XMLWriteOpeningTag( getXMLName(), out, tabs, mName );

    // write out the market string.
    XMLWriteElement( mMarket, "market", out, tabs );
    
    // write whether we have a tax or subsidy
    XMLWriteElement( mPolicyType, "policyType", out, tabs );
    
    // Write whether we are a fixed tax policy.
    XMLWriteElement( mIsShareBased, "isShareBased", out, tabs );
    
    // Write whether we are a fixed tax policy.
    XMLWriteElement( isFixedTax, "isFixedTax", out, tabs );
    
    // Write the constraint for the current year
    XMLWriteElement( mConstraint[ period ], "constraint", out, tabs );
    
    // Write out the fixed tax for the current year.
    XMLWriteElement( mFixedTax[ period ], "fixedTax", out, tabs );
    
    // Write out the share for the current year.
    XMLWriteElement( mShareOfSectorOutput[ period ], "share-of-sector-output", out, tabs );
    
    // finished writing xml for the class members.
    XMLWriteClosingTag( getXMLName(), out, tabs );
}

/*! \brief Complete the initialization of the portfolio standard policy.
* \details This function initializes a market for the policy.
* Policy markets are created for both constraint and share policies.
* Both markets are solved.
* \author Sonny Kim
* \param regionName The name of the region the policy controls. 
*/
void PolicyPortfolioStandard::completeInit( const string& aRegionName ) {
    const Modeltime* modeltime = scenario->getModeltime();
    Marketplace* marketplace = scenario->getMarketplace();
    // Create the policy market, a solved market of GHG type which
    // sets the supply side as the constraint and the demand side
    // as the calculated value.

	if ( mPolicyType == "tax") {
        marketplace->createMarket( aRegionName, mMarket, mName, IMarketType::TAX );
    }
	else if ( mPolicyType == "RES") {
        marketplace->createMarket( aRegionName, mMarket, mName, IMarketType::RES );	
	} 
    else {
        marketplace->createMarket( aRegionName, mMarket, mName, IMarketType::SUBSIDY );
    }

    // Set price and output units for period 0 market info.
    IInfo* marketInfo = marketplace->getMarketInfo( mName, aRegionName, 0, true );
    marketInfo->setString( "price-unit", "1975$/GJ" );
    marketInfo->setString( "output-unit", "EJ_or_Share" );

        // Put the taxes in the market as the market prices if it is a fixed tax policy.
    if( isFixedTax ){
        // Set any taxes that are not unset.
        for( unsigned int i = 0; i < mFixedTax.size(); ++i ){
            // Make sure that the market is not solved. It could have been set
            // to solve by an earlier run.
            marketplace->unsetMarketToSolve( mName, aRegionName, i );
            if( mFixedTax[ i ] != -1 ){
                marketplace->setPrice( mName, aRegionName, mFixedTax[ i ], i );
            }
        }
    }       
    // Otherwise solve the market, given the read-in constraint.
    else {
        // Initialize temporary vector of contraints to constraint
        vector<double> tempConstraint( mConstraint );

        // Override tempConstraint with shares if shared based.
        // Note: the share is based on the total output of the sector that the
        // technology is in.
        if( mIsShareBased ){
            tempConstraint = mShareOfSectorOutput;
            marketInfo->setBoolean( "isShareBased", true );
        }
        // Set either of the constraints, quantity or share, into the
        // DEMAND side of the market for a subsidy, so that increasing subsidy 
        // (market price) increases supply, and into the SUPPLY side of the 
        // market for a tax, so increasing tax decreases demand.
        // USING SUPPLY SIDE AS THE DEMAND FOR SUBSIDY AND DEMAND SIDE AS THE
        // CONSTRAINT. TAXES USE DEMAND SIDE FOR DEMAND AND SUPPLY SIDE FOR 
        // CONSTRAINT. DEFAULT IS SUBSIDY
        for( int per = 1; per < modeltime->getmaxper(); ++per ){
            // Subtracting the current demand for this period to set the constraint
            // because addToDemand adds to any existing demand in the market.
            // Passing false to suppress a warning the first time through.
            if( tempConstraint[ per ] != -1 ){
                if ( mPolicyType == "tax" ){
                    marketplace->setMarketToSolve( mName, aRegionName, per );
					marketplace->addToSupply( mName, aRegionName, tempConstraint[ per ] - 
                        marketplace->getSupply( mName, aRegionName, per ), per, false );
                }
				else if ( mPolicyType == "RES" ){  // maw doesn't understand this
                    marketplace->setMarketToSolve( mName, aRegionName, per );
				//	maw doesn't understand this.  But it doesn;t work otherwise
					marketplace->addToSupply( mName, aRegionName, tempConstraint[ per ] - 
                        marketplace->getSupply( mName, aRegionName, per ), per, false );
                }
                else {
                    marketplace->setMarketToSolve( mName, aRegionName, per );
                    marketplace->addToDemand( mName, aRegionName, tempConstraint[ per ] - 
                        marketplace->getDemand( mName, aRegionName, per ), per, false );
                }
            }
        }
    }
}

/*!
* \brief Set the constraint to the vector passed in.
* \param aConstraint new constraint vector
*/
void PolicyPortfolioStandard::setShareConstraint( const vector<double>& aConstraint ){
    mIsShareBased = true;
    mShareOfSectorOutput = aConstraint;
}

/*!
* \brief Set the constraint to the vector passed in.
* \param aConstraint new constraint vector
*/
void PolicyPortfolioStandard::setQuantityConstraint( const vector<double>& aConstraint ){
    mIsShareBased = false;
    isFixedTax = false;
    mConstraint = aConstraint;
}
