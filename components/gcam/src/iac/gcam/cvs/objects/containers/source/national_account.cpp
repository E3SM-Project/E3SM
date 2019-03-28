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
* \file national_account.cpp
* \ingroup Objects
* \brief The National ccount class source file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include "util/base/include/definitions.h"
#include <cassert>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>

#include "containers/include/national_account.h"
#include "util/base/include/ivisitor.h"
#include "util/base/include/xml_helper.h"
#include "util/base/include/util.h"
#include "util/logger/include/ilogger.h"

using namespace std;
using namespace xercesc;

//! Default Constructor
NationalAccount::NationalAccount():
// Size the accounts to one after the last valid value in the vector,
// represented by END.
mAccounts( END )
{
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
const std::string& NationalAccount::getXMLNameStatic() {
    const static string XML_NAME = "nationalAccount";
    return XML_NAME;
}


void NationalAccount::XMLParse( const DOMNode* node ) {
    /*! \pre make sure we were passed a valid node. */
    assert( node );

    // get all child nodes.
    DOMNodeList* nodeList = node->getChildNodes();

    // loop through the child nodes.
    for( unsigned int i = 0; i < nodeList->getLength(); i++ ){
        DOMNode* curr = nodeList->item( i );
        const string nodeName = XMLHelper<string>::safeTranscode( curr->getNodeName() );

        if( nodeName == "#text" ) {
            continue;
        }
        else if  ( nodeName == enumToXMLName( RETAINED_EARNINGS ) ) {
            addToAccount( RETAINED_EARNINGS, XMLHelper<double>::getValue( curr ) );
        }
        else if  ( nodeName == enumToXMLName( DIVIDENDS ) ) {
            addToAccount( DIVIDENDS, XMLHelper<double>::getValue( curr ) );
        }
        else if  ( nodeName == enumToXMLName( TRANSFERS ) ) {
            addToAccount( TRANSFERS, XMLHelper<double>::getValue( curr ) );
        }
        else if  ( nodeName == enumToXMLName( CORPORATE_INCOME_TAX_RATE ) ) {
            setAccount( CORPORATE_INCOME_TAX_RATE, XMLHelper<double>::getValue( curr ) );
        }
        else if  ( nodeName == enumToXMLName( EXCHANGE_RATE ) ) {
            addToAccount( EXCHANGE_RATE, XMLHelper<double>::getValue( curr ) );
        }
        else if  ( nodeName == enumToXMLName( CORPORATE_PROFITS ) ) {
            addToAccount( CORPORATE_PROFITS, XMLHelper<double>::getValue( curr ) );
        }
        else {
            ILogger& mainLog = ILogger::getLogger( "main_log" );
            mainLog.setLevel( ILogger::WARNING );
            mainLog << "Unrecognized text string: " << nodeName << " found while parsing " << getXMLNameStatic() << "." << endl;
        }
    }
}

//! Output debug info to XML
void NationalAccount::toDebugXML( const int period, ostream& out, Tabs* tabs ) const {
    XMLWriteOpeningTag( getXMLNameStatic(), out, tabs );

        for( int i = 0; i < NationalAccount::END; ++i ) {
            XMLWriteElement( getAccountValue( static_cast< NationalAccount::AccountType >( i ) ),
            enumToXMLName( static_cast< NationalAccount::AccountType >( i ) ),
            out, tabs );
    }
    
    XMLWriteClosingTag( getXMLNameStatic(), out, tabs );
}

/*!
 * \brief Add to the value for the national account specified by the account type key.
 * \param aType The account which will be added to.
 * \param aValue The value which will be added to account aType.
 */
void NationalAccount::addToAccount( const AccountType aType, const double aValue ){
    mAccounts[ aType ]+= aValue;
}

/*!
 * \brief Set the value for the national account specified by the account type key.
 * \param aType The account for which value will be set.
 * \param aValue The value which wil be set to account aType.
 */
void NationalAccount::setAccount( const AccountType aType, const double aValue ){
    mAccounts[ aType ] = aValue;
}

/*!
 * \brief Get the value for the national account specified by the account type key.
 * \param aType The account who's value will be returned.
 * \return The current value of the account aType.
 */
double NationalAccount::getAccountValue( const AccountType aType ) const {
    assert( aType < END );
    return mAccounts[ aType ];
}

/*!
 * \brief Reset the national account values
 * \note The accounts TRANSFERS and INVESTMENT_TAX_CREDIT do not get reset.
 */
void NationalAccount::reset() {
    // save transfer amount calculated from government's initCalc
    // this is a giant hack. store the values that should not be reset.
    double tempTransfers = mAccounts[ TRANSFERS ];
    assert( tempTransfers > 0 );
    double tempInvTaxCreditRate = mAccounts[ INVESTMENT_TAX_CREDIT ];
    
    // Clear the values.
    mAccounts.clear();
    mAccounts.resize( END );

    // Restore the two values that were saved.
    mAccounts[ TRANSFERS ] = tempTransfers;
    mAccounts[ INVESTMENT_TAX_CREDIT ] = tempInvTaxCreditRate;
}

/*!
 * \brief For outputting SGM data to a flat csv File
 * \param aFile The file in which the results will be written.
 * \param period The period which we are outputting for
 * \author Pralit Patel
 */
void NationalAccount::csvSGMOutputFile( ostream& aFile, const int period ) const {
    // write headers
    aFile << "GNP Accounts" << endl;
    aFile << enumToName( GNP_NOMINAL ) << ',';
    aFile << enumToName( CONSUMPTION_NOMINAL ) << ',';
    aFile << enumToName( GOVERNMENT_NOMINAL ) << ',';
    aFile << enumToName( INVESTMENT_NOMINAL ) << ',';
    aFile << enumToName( NET_EXPORT_NOMINAL ) << ',';
    aFile << endl;
    // write data
    aFile << getAccountValue( GNP_NOMINAL ) << ',';
    aFile << getAccountValue( CONSUMPTION_NOMINAL ) << ',';
    aFile << getAccountValue( GOVERNMENT_NOMINAL ) << ',';
    aFile << getAccountValue( INVESTMENT_NOMINAL ) << ',';
    aFile << getAccountValue( NET_EXPORT_NOMINAL ) << ',';
    aFile << endl << endl;
    aFile << enumToName( GNP_REAL ) << ',';
    aFile << enumToName( CONSUMPTION_REAL ) << ',';
    aFile << enumToName( GOVERNMENT_REAL ) << ',';
    aFile << enumToName( INVESTMENT_REAL ) << ',';
    aFile << enumToName( NET_EXPORT_REAL ) << ',';
    aFile << endl;
    // write data
    aFile << getAccountValue( GNP_REAL ) << ',';
    aFile << getAccountValue( CONSUMPTION_REAL ) << ',';
    aFile << getAccountValue( GOVERNMENT_REAL ) << ',';
    aFile << getAccountValue( INVESTMENT_REAL ) << ',';
    aFile << getAccountValue( NET_EXPORT_REAL ) << ',';
    aFile << endl << endl;
    // write tax accounts
    aFile << "National Tax Accounts" << endl;
    // write headers
    aFile << enumToName( CORPORATE_INCOME_TAXES ) << ',';
    aFile << enumToName( INDIRECT_BUSINESS_TAX ) << ',';
    aFile << enumToName( SOCIAL_SECURITY_TAX ) << ',';
    aFile << enumToName( PERSONAL_INCOME_TAXES ) << ',';
    aFile << enumToName( INVESTMENT_TAX_CREDIT ) << ',';
    aFile << "Total Taxes";
    aFile << endl;
    // write data
    aFile << getAccountValue( CORPORATE_INCOME_TAXES ) << ',';
    aFile << getAccountValue( INDIRECT_BUSINESS_TAX ) << ',';
    aFile << getAccountValue( SOCIAL_SECURITY_TAX ) << ',';
    aFile << getAccountValue( PERSONAL_INCOME_TAXES ) << ',';
    aFile << getAccountValue( INVESTMENT_TAX_CREDIT ) << ',';
    aFile << getAccountValue( CORPORATE_INCOME_TAXES )
            + getAccountValue( INDIRECT_BUSINESS_TAX )
            + getAccountValue( SOCIAL_SECURITY_TAX )
            + getAccountValue( PERSONAL_INCOME_TAXES )
            + getAccountValue( INVESTMENT_TAX_CREDIT );
    aFile << endl << endl;
    // write import export breakdown
    aFile << "National Imports and Exports" << endl;
    // write headers
    aFile << enumToName( EXPORT_NOMINAL ) << ',';
    aFile << enumToName( EXPORT_REAL ) << ',';
    aFile << enumToName( IMPORT_NOMINAL ) << ',';
    aFile << enumToName( IMPORT_REAL );
    aFile << endl;
    // write data
    aFile << getAccountValue( EXPORT_NOMINAL ) << ',';
    aFile << getAccountValue( EXPORT_REAL ) << ',';
    aFile << getAccountValue( IMPORT_NOMINAL ) << ',';
    aFile << getAccountValue( IMPORT_REAL );
    aFile << endl << endl;

}

/*!
 * \brief Convert between the NationalAccount enum type to the String representation
 * \param aType The enum NationalAccount type
 * \return The string representation of the type
 * \author Pralit Patel
 */
const string& NationalAccount::enumToName( const AccountType aType ) const {
    /*! \pre aType is a valid account type. */
    assert( aType < END );
    // Create a static array of values. This will only on the first entrance to
    // the function since this is a const static.
    const static string names[] = {
            "Retained Earnings",
            "Subsidy",
            "Corporate Profits",
            "Corporate Retained Earnings",
            "Corporate Income Taxes",
            "Corporate Income Tax Rate",
            "Personal Income Taxes",
            "Investment Tax Credit",
            "Dividends",
            "Labor Wages",
            "Land Rents",
            "Transfers",
            "Social Security Tax",
            "Indirect Buisness Tax",
            "GNP Nominal",
            "GNP VA",
            "GNP Real",
            "Consumption Nominal",
            "Government Nominal",
            "Investment Nominal",
            "Net Export Nominal",
            "Consumption Real",
            "Government Real",
            "Investment Real",
            "Net Export Real",
            "Exchange Rate",
            "Annual Investment",
            "Carbon Tax",
            "Export Nominal",
            "Export Real",
            "Import Nominal",
            "Import Real"
    };
    // Return the string in the static array at the index.
    return names[ aType ];
}

/*!
 * \brief Convert between the NationalAccount enum type to the XML String representation
 * \param aType The enum NationalAccount type
 * \return The XML string representation of the type
 * \author Pralit Patel
 */
const string& NationalAccount::enumToXMLName( const AccountType aType ) const {
    /*! \pre aType is a valid account type. */
    assert( aType < END );
    // Create a static array of values. This will only on the first entrance to
    // the function since this is a const static.

    const static string names[] = {
            "retainedEarning",
            "subsidy",
            "corpProfits",
            "corpRetainedEarnings",
            "corporateIncomeTax",
            "corporateIncomeTaxRate",
            "personalIncomeTax",
            "investmentTaxCredit",
            "dividends",
            "laborWages",
            "landRents",
            "transfers",
            "socialSecurityTax",
            "indirectBusinessTax",
            "GNP-nominal",
            "GNPVA",
            "GNP-real",
            "consumption-nominal",
            "government-nominal",
            "investment-nominal",
            "netExport-nominal",
            "consumption-real",
            "government-real",
            "investment-real",
            "netExport-real",
            "exchangeRate",
            "annualInvestment",
            "carbon-tax",
            "export-nominal",
            "export-real",
            "import-nominal",
            "import-real"
    };
    // Return the string in the static array at the index.
    return names[ aType ];
}

// for reporting National Account information
void NationalAccount::accept( IVisitor* aVisitor, const int aPeriod ) const {
    aVisitor->startVisitNationalAccount( this, aPeriod );
    aVisitor->endVisitNationalAccount( this, aPeriod );
}
