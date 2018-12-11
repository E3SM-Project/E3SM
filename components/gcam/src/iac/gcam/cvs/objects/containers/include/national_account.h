#ifndef _NATIONAL_ACCOUNT_H_
#define _NATIONAL_ACCOUNT_H_
#if defined(_MSC_VER)
#pragma once
#endif

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
* \file national_account.h
* \ingroup Objects
* \brief The NationalAccount class header file.
* \author Pralit Patel
* \author Sonny Kim
*/

#include <vector>
#include <string>
#include <xercesc/dom/DOMNode.hpp>
#include "util/base/include/ivisitable.h"

class Tabs;

/*! 
 * \ingroup Objects
 * \brief A container of national accounts information for the Region.
 * \details This object defines an accoutning structure for a Region in a given
 *          period.  The accounts are used for both reporting as well as model
 *          operation.  The accounts critical for model operation would be:
 *           -ANNUAL_INVESTMENT: To account for the demand of Capital investment.
 *           -LABOR_WAGES: Value paid by producers to housholds for labor.
 *           -LAND_RENTS: Value paid by producers to housholds for access to land and
 *                        resources.
 *           -CARBON_TAX: Account for taxes paid for emissions by producers and consumers so
 *                        they may be added into a pool of taxes that will go the the Government
 *                        consumer.
 *           -INDIRECT_BUSINESS_TAX: Account for taxes paid by producers for both consumption (not
 *                                   including labor) and production.  This account will then be
 *                                   added to the pool of taxes given to he Government consumer.
 *           -SOCIAL_SECURITY_TAX:  Taxes paid for labor inputs to be added to the pool of taxes
 *                                  given to the Government consumer.
 *           -TRANSFERS: The value given back to the housholds by the government consumer.
 *           -CORPORATE_INCOME_TAXES: Taxes paid by producers on any profits they may have.
 *           -CORPORATE_INCOME_TAX_RATE:  Rate at which producers are taxed on profits.  This paramater
 *                                        is read in.
 *           -CORPORATE_PROFITS: Read in value used to back out the split between retained earnings
 *                               and profits for a producer.
 *           -RETAINED_EARNINGS: Read in value used to back out the split between retained earnings
 *                               and profits for a producer.
 * 
 *          <b>XML specification for NationalAccount</b>
 *          - XML name: \c nationalAccount
 *          - Contained by: RegionCGE
 *          - Parsing inherited from class: None.
 *          - Attributes: year
 *          - Elements:
 *              - \c retainedEarning NationalAccount::mAccounts[ RETAINED_EARNINGS ]
 *              - \c dividends NationalAccount::mAccounts[ DIVIDENDS ]
 *              - \c transfers NationalAccount::mAccounts[ TRANSFERS ]
 *              - \c corporateIncomeTaxRate NationalAccount::mAccounts[ CORPORATE_INCOME_TAX_RATE ]
 *              - \c exchangeRate NationalAccount::mAccounts[ EXCHANGE_RATE ]
 *              - \c corpProfits NationalAccount::mAccounts[ CORPORATE_PROFITS ]
 *
 * \author Pralit Patel, Sonny Kim
 */
class NationalAccount: public IVisitable
{
    friend class IVisitor;
    friend class XMLDBOutputter;
public:
    //! All the accounts contained in a NationalAccount
    enum AccountType {
        RETAINED_EARNINGS,
        SUBSIDY,
        CORPORATE_PROFITS,
        CORPORATE_RETAINED_EARNINGS,
        CORPORATE_INCOME_TAXES,
        CORPORATE_INCOME_TAX_RATE,
        PERSONAL_INCOME_TAXES,
        INVESTMENT_TAX_CREDIT,
        DIVIDENDS,
        LABOR_WAGES,
        LAND_RENTS,
        TRANSFERS,
        SOCIAL_SECURITY_TAX,
        INDIRECT_BUSINESS_TAX,
        GNP_NOMINAL,
        GNP_VA,
        GNP_REAL,
        CONSUMPTION_NOMINAL,
        GOVERNMENT_NOMINAL,
        INVESTMENT_NOMINAL,
        NET_EXPORT_NOMINAL,
        CONSUMPTION_REAL,
        GOVERNMENT_REAL,
        INVESTMENT_REAL,
        NET_EXPORT_REAL,
        EXCHANGE_RATE,
        ANNUAL_INVESTMENT,
        CARBON_TAX,
        EXPORT_NOMINAL,
        EXPORT_REAL,
        IMPORT_NOMINAL,
        IMPORT_REAL,
        // Insert new values before END marker.
        END
    };
    NationalAccount();
    static const std::string& getXMLNameStatic();
    void XMLParse( const xercesc::DOMNode* node );
    void toDebugXML( const int period, std::ostream& out, Tabs* tabs ) const;
    void reset();
    void addToAccount( const AccountType aType, const double aValue );
    void setAccount( const AccountType aType, const double aValue );
    double getAccountValue( const AccountType aType ) const;
    void csvSGMOutputFile( std::ostream& aFile, const int period ) const;
    void accept( IVisitor* aVisitor, const int aPeriod ) const;
private:
    const std::string& enumToName( const AccountType aType ) const;
    const std::string& enumToXMLName( const AccountType aType ) const;
    //! Vector to hold national account values
    std::vector<double> mAccounts;
};

#endif // _NATIONAL_ACCOUNT_H_
